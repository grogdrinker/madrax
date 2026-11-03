"""Shared fixtures for the MadraX test suite.

Run against the patched tree (default):      pytest
Run the bug tests against upstream master:    MADRAX_PATH=/path/to/master pytest --target=master
In master mode every test marked ``xfail_master`` becomes a *strict* xfail, so the
run is green only if each of those tests really fails on master, i.e. it detects its bug.
"""
import gzip
import os
import sys
import warnings

import pytest
import torch

if os.environ.get("MADRAX_PATH"):
    sys.path.insert(0, os.environ["MADRAX_PATH"])
else:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

warnings.filterwarnings("ignore", category=FutureWarning)  # torch.load in ForceField.__init__

from madrax.ForceField import ForceField  # noqa: E402
from madrax import utils, dataStructures  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TERMS = ["disulfide", "hbond", "electro", "clash", "solv_polar", "solv_hydrophobic",
         "vdw", "entropy_mc", "entropy_sc", "peptide_bond", "rotamer"]
ENERGY_ATOL = 1e-5  # master itself is non-deterministic at ~5e-7 (scatter order)


# ---------------------------------------------------------------- options / markers
def pytest_addoption(parser):
    parser.addoption("--target", choices=["patch", "master"], default="patch",
                     help="'master' turns xfail_master tests into strict xfails")


def pytest_configure(config):
    for m in ("xfail_master(reason): bug test - must fail on upstream master",
              "known_gap(reason): documented gap not fixed by this PR (strict xfail everywhere)",
              "slow: uses the 6k-atom structure or large scans",
              "gpu: needs CUDA",
              "benchmark: timing / memory, never part of the PR gate"):
        config.addinivalue_line("markers", m)


def pytest_collection_modifyitems(config, items):
    master = config.getoption("--target") == "master"
    for item in items:
        m = item.get_closest_marker("xfail_master")
        if master and m:
            item.add_marker(pytest.mark.xfail(strict=True, reason="bug present on master: " + (m.args[0] if m.args else "")))
        g = item.get_closest_marker("known_gap")
        if g:
            item.add_marker(pytest.mark.xfail(strict=True, reason="known gap: " + (g.args[0] if g.args else "")))
        if item.get_closest_marker("gpu") and not torch.cuda.is_available():
            item.add_marker(pytest.mark.skip(reason="no CUDA device"))


@pytest.fixture(scope="session")
def target(pytestconfig):
    return pytestconfig.getoption("--target")


# ---------------------------------------------------------------- PDB text helpers
def read_lines(name):
    with gzip.open(os.path.join(DATA, name + ".pdb.gz"), "rt") as fh:
        return [l for l in fh if l.startswith("ATOM")]


def resnum(line):
    return int(line[22:26])


def atname(line):
    return line[12:16].strip()


def set_resnum(line, n):
    return line[:22] + "%4d" % n + line[26:]


def write(path, lines):
    with open(path, "w") as fh:
        fh.writelines(l if l.endswith("\n") else l + "\n" for l in lines)
        fh.write("END\n")
    return str(path)


# ---------------------------------------------------------------- structure files
@pytest.fixture(scope="session")
def pdbdir(tmp_path_factory):
    """Every fixture structure written as a plain .pdb file (the parser reads paths)."""
    d = tmp_path_factory.mktemp("pdbs")
    small = read_lines("1A8O")          # 1 chain, residues 152-220, breaks at 185 and 214-215
    multi = read_lines("2BEG")          # 5 identical chains, PHE19/PHE20 in each
    large = read_lines("2XHE_trimmed")  # 2 chains, residues 0-616 / 2-261
    out = {
        "small": write(d / "small.pdb", small),
        "multichain": write(d / "multichain.pdb", multi),
        "large": write(d / "large.pdb", large),
    }
    for k in (-151, -160, 1000):
        out["renumbered%+d" % k] = write(d / ("renumbered%+d.pdb" % k), [set_resnum(l, resnum(l) + k) for l in small])
    last = max(resnum(l) for l in small)
    out["truncated_last"] = write(d / "truncated_last.pdb",
                                  [l for l in small if not (resnum(l) == last and atname(l) != "N")])
    out["missing_ca_middle"] = write(d / "missing_ca_middle.pdb",
                                     [l for l in small if not (resnum(l) == 200 and atname(l) == "CA")])
    cut = (154, 162, 167, 170, 182, 199, 203)
    out["truncated_sidechains"] = write(d / "truncated_sidechains.pdb",
                                        [l for l in small if not (resnum(l) in cut and atname(l) not in ("N", "CA", "C", "O", "CB"))])
    # ASP163 with OD1 missing but OD2 present: OD2's second partner (OD1) is absent
    out["missing_partner"] = write(d / "missing_partner.pdb",
                                   [l for l in small if not (resnum(l) == 163 and atname(l) == "OD1")])
    alt = []
    for l in small:
        if resnum(l) == 170 and atname(l) not in ("N", "CA", "C", "O"):
            alt.append(l[:16] + "A" + l[17:])
            alt.append(l[:16] + "B" + l[17:30] + "%8.3f" % (float(l[30:38]) + 0.4) + l[38:])
        else:
            alt.append(l)
    out["altloc"] = write(d / "altloc.pdb", alt)

    def occ(line, o):
        return line[:54] + "%6.2f" % o + line[60:]
    # residue 170 side chain: conformer B has the higher total occupancy (0.6 vs 0.4)
    out["altloc_b_major"] = write(d / "altloc_b_major.pdb",
                                  [occ(l, 0.4 if l[16] == "A" else 0.6) if l[16] in "AB" else l for l in alt])
    # A wins on total occupancy, but one atom (the first side-chain atom) has B > A:
    # per-atom picking would mix conformers, per-residue picking must not
    mixed, first = [], None
    for l in alt:
        if l[16] in "AB":
            first = first or atname(l)
            o = (0.3 if l[16] == "A" else 0.7) if atname(l) == first else (0.6 if l[16] == "A" else 0.4)
            l = occ(l, o)
        mixed.append(l)
    out["altloc_mixed"] = write(d / "altloc_mixed.pdb", mixed)
    # residue 190 split into 190 + 190A (insertion code): side chain moved to an inserted residue
    out["icode"] = write(d / "icode.pdb",
                         [set_resnum(l, 190)[:26] + "A" + l[27:] if resnum(l) == 191 else l for l in small])
    hyd = []
    for l in small:
        hyd.append(l)
        if atname(l) == "CB":
            for k, nm in enumerate(("1HB ", "2HB ")):
                x = float(l[30:38]) + 0.6 * (k + 1)
                hyd.append(l[:12] + nm + l[16:30] + "%8.3f" % x + l[38:76] + " H  \n")
    out["h_named"] = write(d / "h_named.pdb", hyd)
    return out


@pytest.fixture(scope="session")
def parse(pdbdir):
    """parse('small') -> (coords, atom_names). Works with both 2- and 3-value parsePDB returns."""
    cache = {}

    def _parse(name):
        if name not in cache:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out = utils.parsePDB(pdbdir[name])
            cache[name] = (out[0], out[1])
        return cache[name]
    return _parse


# ---------------------------------------------------------------- force field
@pytest.fixture(scope="session")
def ff():
    torch.manual_seed(0)
    return ForceField(device="cpu")


@pytest.fixture(scope="session")
def energy(ff):
    """energy(coords, names) -> full (B, chains, res, alt, 11) tensor, no grad."""
    def _energy(coords, names):
        info = dataStructures.create_info_tensors(names)
        with torch.no_grad():
            return ff(coords, info)
    return _energy


def pad_batch(*parsed):
    coords = torch.nn.utils.rnn.pad_sequence([c[0][c[0][:, 0] != -999] for c, _ in parsed],
                                             batch_first=True, padding_value=-999.0)
    names = [n[0] for _, n in parsed]
    return coords, names


def per_residue(e, b=0, alt=0):
    """(chains, res, 11) energies of one protein / alternative."""
    return e[b, :, :, alt, :]


def assert_energy_close(a, b, atol=ENERGY_ATOL, what="energy"):
    assert a.shape == b.shape, "%s: shape %s vs %s" % (what, tuple(a.shape), tuple(b.shape))
    d = (a - b).abs()
    if d.max() > atol:
        idx = tuple(int(i) for i in (d == d.max()).nonzero()[0])
        raise AssertionError("%s differs: max |diff| %.5g at index %s (term %s); totals %.4f vs %.4f"
                             % (what, float(d.max()), idx, TERMS[idx[-1]], float(a.sum()), float(b.sum())))
