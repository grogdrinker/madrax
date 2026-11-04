"""mutate.rotator and StructureOptimizer: every torsion parameter must move exactly its own dihedral.

The oracle is geometric and independent of MadraX: dihedrals are measured with numpy.
"""
import math

import numpy as np
import pytest
import torch

from conftest import dataStructures
from madrax.mutate import rotator, StructureOptimizer
from madrax.mutate import mutatingEngine
from madrax.sources import hashings

AD = hashings.atom_description_hash


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def wrap(d):
    return (d + 180.0) % 360.0 - 180.0


class Struct:
    """Flat per-atom view of one protein: lookup by (chain, resnum, atom name)."""

    def __init__(self, coords, names):
        self.c, self.n = coords, names
        self.info = dataStructures.create_info_tensors(names)
        ad = self.info[1]
        self.ad = ad
        self.names = [hashings.atom_hash_inverse[int(a)] for a in ad[:, AD["at_name"]]]
        self.names = ["N" if a == "tN" else "O" if a == "OXT" else a for a in self.names]
        self.res = [hashings.resi_hash_inverse[int(r)] for r in ad[:, AD["resname"]]]
        self.key = {(int(ad[i, AD["chain"]]), int(ad[i, AD["resnum"]]), self.names[i]): i for i in range(len(self.names))}
        self.residues = sorted({(int(ad[i, AD["chain"]]), int(ad[i, AD["resnum"]]), self.res[i]) for i in range(len(self.names))})
        self.shape = (1, int(ad[:, AD["chain"]].max()) + 1, int(ad[:, AD["resnum"]].max()) + 1, info_alt(self.info), 8)

    def flat(self, coords):
        return coords[self.ad[:, 0].long(), self.info[2]].detach().numpy()

    def dihedral(self, x, ch, atoms):
        """atoms: list of (resnum, name)."""
        if not all((ch, r, a) in self.key for r, a in atoms):
            return None
        return dihedral(*[x[self.key[(ch, r, a)]] for r, a in atoms])


def info_alt(info):
    return info[5].shape[-1]


def backbone_dihedrals(s, x):
    out = {}
    for ch, r, _ in s.residues:
        out[(ch, r, "phi")] = s.dihedral(x, ch, [(r - 1, "C"), (r, "N"), (r, "CA"), (r, "C")])
        out[(ch, r, "psi")] = s.dihedral(x, ch, [(r, "N"), (r, "CA"), (r, "C"), (r + 1, "N")])
        out[(ch, r, "omega")] = s.dihedral(x, ch, [(r, "CA"), (r, "C"), (r + 1, "N"), (r + 1, "CA")])
    return {k: v for k, v in out.items() if v is not None}


def bond_lengths(s, x):
    """All intra-residue distances < 1.9 A plus the C(i)-N(i+1) peptide bonds."""
    idx = np.arange(len(s.names))
    d = {}
    by_res = {}
    for (ch, r, a), i in s.key.items():
        by_res.setdefault((ch, r), []).append(i)
    for (ch, r), ii in by_res.items():
        for a in ii:
            for b in ii:
                if a < b and np.linalg.norm(x[a] - x[b]) < 1.9:
                    d[(a, b)] = np.linalg.norm(x[a] - x[b])
        if (ch, r, "C") in s.key and (ch, r + 1, "N") in s.key:
            a, b = s.key[(ch, r, "C")], s.key[(ch, r + 1, "N")]
            d[(a, b)] = np.linalg.norm(x[a] - x[b])
    return d


# ---------------------------------------------------------------- backbone
@pytest.mark.xfail_upstream("phi and psi read the omega slot; psi moves residue i's N and side chain")
@pytest.mark.parametrize("slot,which", [(0, "phi"), (1, "psi"), (2, "omega")])
def test_backbone_torsion_moves_only_its_dihedral(parse, slot, which):
    c, n = parse("small")
    s = Struct(c, n)
    target = (0, 190)
    rot = torch.zeros(s.shape)
    rot[0, target[0], target[1], 0, slot] = math.radians(30)
    moving = s.ad[:, AD["resnum"]] == target[1]
    x0 = s.flat(c)
    coords0 = torch.tensor(x0)
    out = rotator.RotateStruct().rotateBackBone(coords0.clone(), s.ad, rot, moving).numpy()
    before, after = backbone_dihedrals(s, x0), backbone_dihedrals(s, out)
    # the rotator's omega slot at residue i turns the C(i-1)-N(i) bond, which the scorer calls omega(i-1)
    expected = (target[0], target[1] - 1 if which == "omega" else target[1], which)
    for k in before:
        delta = wrap(after[k] - before[k])
        if k == expected:
            assert abs(abs(delta) - 30) < 0.01, "%s changed by %.3f deg, expected 30" % (k, delta)
        else:
            assert abs(delta) < 0.01, "%s changed by %.3f deg, expected 0" % (k, delta)
    b0, b1 = bond_lengths(s, x0), bond_lengths(s, out)
    worst = max(abs(b1[k] - b0[k]) for k in b0)
    assert worst < 1e-3, "a bond length changed by %.4f A" % worst
    upstream = (s.ad[:, AD["resnum"]] < target[1]).numpy()
    assert np.allclose(out[upstream], x0[upstream], atol=1e-5), "atoms before the rotated residue moved"


@pytest.mark.xfail_upstream("raise('...') raises TypeError")
def test_backbone_rotation_rejects_mutants(parse):
    c, n = parse("small")
    cm, nm = mutatingEngine.mutate(c, n, [[["178_A_GLU"]]])
    s = Struct(cm, nm)
    with pytest.raises(ValueError):
        rotator.RotateStruct().rotateBackBone(torch.tensor(s.flat(cm)), s.ad, torch.zeros(s.shape)[..., :3],
                                              s.ad[:, AD["resnum"]] == 190)


# ---------------------------------------------------------------- side chains
CHIS = ["chi1", "chi2", "chi3", "chi4", "chi5"]
NOT_ROTATABLE = {("PRO", "chi1"), ("PRO", "chi2")}   # ring closure: excluded by design
KNOWN_GAP = {("ILE", "chi2")}                        # rotator.py has no ILE chi2 axis


def _chi_cases(s, chi):
    for ch, r, res in s.residues:
        if res in hashings.chi_atoms[chi]:
            yield ch, r, res, [(r, a) for a in hashings.chi_atoms[chi][res]]


@pytest.fixture(scope="module")
def zoo(parse):
    """The 2-chain structure: contains all 20 residue types."""
    return Struct(*parse("large"))


@pytest.mark.slow
@pytest.mark.parametrize("k,chi", list(enumerate(CHIS)))
def test_chi_rotation_moves_exactly_that_chi(parse, zoo, k, chi):
    s = zoo
    rot = torch.zeros(s.shape)
    rot[..., 3 + k] = math.radians(30)
    with torch.no_grad():
        out = s.flat(rotator.RotateStruct()(s.c, s.info, rot, backbone_rotation=False))
    x0 = s.flat(s.c)
    wrong = []
    for ch, r, res, atoms in _chi_cases(s, chi):
        if (res, chi) in NOT_ROTATABLE | KNOWN_GAP:
            continue
        d0, d1 = s.dihedral(x0, ch, atoms), s.dihedral(out, ch, atoms)
        if d0 is None:
            continue
        if abs(abs(wrap(d1 - d0)) - 30) > 0.05:
            wrong.append((res, ch, r, round(wrap(d1 - d0), 3)))
        for other in CHIS:  # every other chi of the residue is rigid
            if other == chi or res not in hashings.chi_atoms[other]:
                continue
            o = [(r, a) for a in hashings.chi_atoms[other][res]]
            e0, e1 = s.dihedral(x0, ch, o), s.dihedral(out, ch, o)
            if e0 is not None and abs(wrap(e1 - e0)) > 0.05:
                wrong.append((res, ch, r, other, round(wrap(e1 - e0), 3)))
    assert not wrong, "%d residues wrong, e.g. %s" % (len(wrong), wrong[:5])
    assert np.allclose(np.array([v for v in backbone_dihedrals(s, out).values()]),
                       np.array([v for v in backbone_dihedrals(s, x0).values()]), atol=0.01)
    b0, b1 = bond_lengths(s, x0), bond_lengths(s, out)
    assert max(abs(b1[k] - b0[k]) for k in b0) < 1e-3


# ---------------------------------------------------------------- optimizer
@pytest.mark.xfail_upstream("ReduceLROnPlateau no longer accepts verbose=")
def test_optimize_runs(parse, ff):
    c, n = parse("small")
    info = dataStructures.create_info_tensors(n)
    with torch.no_grad():
        e0 = ff(c, info).sum()
    yp, coords = StructureOptimizer.optimize(ff, c, info, epochs=2, learning_rate=0.01)
    assert torch.isfinite(yp).all()
    assert coords.shape == c.shape
    moved = (coords - c)[c[..., 0] != -999].abs().max()
    assert moved > 0
    assert yp.sum() < e0 + 1.0
