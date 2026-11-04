"""utils.parsePDB: ring centroids, altlocs, insertion codes, CA-less residues."""
import warnings

import pytest
import torch

from conftest import utils, dataStructures, assert_energy_close


def _rc_atoms(names):
    return [a for a in names if a.split("_")[2] == "RC"]


@pytest.mark.xfail_upstream("ring centroids keyed by residue number only: 2 of 10 built, at the wrong ring")
def test_ring_centroid_per_chain(parse):
    coords, names = parse("multichain")
    names = names[0]
    rings = {(a.split("_")[3], a.split("_")[1]) for a in names
             if a.split("_")[0] in ("PHE", "TYR") and a.split("_")[2] == "CA"}
    rc = _rc_atoms(names)
    assert len(rc) == len(rings) == 10
    idx = {a: i for i, a in enumerate(names)}
    for a in rc:
        r, num, _, ch, mut, alt = a.split("_")
        cd1 = coords[0, idx["_".join([r, num, "CD1", ch, mut, alt])]]
        ce2 = coords[0, idx["_".join([r, num, "CE2", ch, mut, alt])]]
        assert torch.allclose(coords[0, idx[a]], (cd1 + ce2) / 2, atol=1e-5), a


@pytest.mark.xfail_upstream("altloc A and B both parsed as real atoms")
def test_altloc_tie_keeps_first_label(parse, energy):
    c0, n0 = parse("small")
    c1, n1 = parse("altloc")          # A and B both at occupancy 1.00
    assert n1 == n0
    assert torch.equal(c1, c0)
    assert_energy_close(energy(c1, n1), energy(c0, n0))


def _residue(coords, names, num):
    return {a.split("_")[2]: coords[0, i] for i, a in enumerate(names[0]) if a.split("_")[1] == str(num)}


@pytest.mark.xfail_upstream("altloc A and B both parsed as real atoms")
def test_altloc_highest_occupancy_wins(parse):
    c0, n0 = parse("small")
    c1, n1 = parse("altloc_b_major")  # B carries 0.6, A 0.4
    assert len(n1[0]) == len(n0[0])
    r0, r1 = _residue(c0, n0, 170), _residue(c1, n1, 170)
    for at in r1:
        shift = 0.4 if at not in ("N", "CA", "C", "O") else 0.0   # the fixture moved B by +0.4 A in x
        assert abs(float(r1[at][0] - r0[at][0]) - shift) < 1e-3, at


@pytest.mark.xfail_upstream("altloc A and B both parsed as real atoms")
def test_altloc_never_mixes_conformers(parse):
    c0, n0 = parse("small")
    c1, n1 = parse("altloc_mixed")    # A wins in total, one atom alone prefers B
    r0, r1 = _residue(c0, n0, 170), _residue(c1, n1, 170)
    assert r0.keys() == r1.keys()
    for at in r1:
        assert torch.allclose(r1[at], r0[at], atol=1e-3), "%s taken from conformer B" % at


@pytest.mark.xfail_upstream("both residue types are kept")
def test_altloc_microheterogeneity(tmp_path):
    """Two residue types at one position (SER in A, THR in B): one residue type survives, never both."""
    from conftest import read_lines, write, resnum, atname
    lines = []
    for l in read_lines("1A8O"):
        if resnum(l) == 178 and atname(l) not in ("N", "CA", "C", "O"):   # SER178 side chain
            lines.append(l[:16] + "A" + l[17:54] + "%6.2f" % 0.3 + l[60:])
            lines.append(l[:16] + "B" + "THR" + l[20:54] + "%6.2f" % 0.7 + l[60:])
        else:
            lines.append(l)
    path = write(tmp_path / "micro.pdb", lines)
    try:
        c, n = utils.parsePDB(path)[:2]
    except ValueError as e:              # acceptable: refuse a residue with two types
        assert "more than one residue type" in str(e)
        return
    types = {a.split("_")[0] for a in n[0] if a.split("_")[1] == "178"}
    assert len(types) == 1


@pytest.mark.xfail_upstream("insertion codes are silently merged")
def test_insertion_code_raises(pdbdir):
    with pytest.raises(ValueError, match="insertion codes"):
        utils.parsePDB(pdbdir["icode"])


@pytest.mark.xfail_upstream("a CA-less last residue crashes AngleScorer; one in the middle is silently scored")
@pytest.mark.parametrize("name", ["truncated_last", "missing_ca_middle"])
def test_ca_less_residue_raises(pdbdir, name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c, n = utils.parsePDB(pdbdir[name])[:2]
    with pytest.raises(ValueError, match="without a CA"):
        dataStructures.create_info_tensors(n)


@pytest.mark.slow
@pytest.mark.xfail_upstream("crashes with IndexError in AngleScorer")
def test_full_2xhe_raises_on_ca_less_residue(tmp_path):
    """2XHE as deposited: its last residue in chain A is a lone N."""
    from conftest import read_lines, write
    lines = read_lines("2XHE_trimmed")
    top = max(range(len(lines)), key=lambda i: int(lines[i][22:26]))
    last = lines[top]
    extra = last[:12] + " N  " + last[16:17] + "ALA" + last[20:22] + "%4d" % (int(last[22:26]) + 1) + last[26:]
    path = write(tmp_path / "with_lone_n.pdb", lines[:top + 1] + [extra] + lines[top + 1:])
    c, n = utils.parsePDB(path)[:2]
    with pytest.raises(ValueError, match="A:617"):
        dataStructures.create_info_tensors(n)
