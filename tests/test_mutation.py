"""mutate.mutatingEngine + the alternative machinery in the force field."""
import pytest
import torch

from conftest import dataStructures, assert_energy_close
from madrax.mutate import mutatingEngine

NON_LYS = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
           "LEU", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
POSITIONS = [155, 163, 170, 178, 190, 201]


def _mut(parse, muts, name="small"):
    c, n = parse(name)
    return mutatingEngine.mutate(c, n, [muts])


def _standalone(cm, nm, alt):
    """The mutant of alternative `alt` rebuilt as a plain wild-type structure."""
    pos = {(a.split("_")[3], a.split("_")[1]) for a in nm[0] if a.split("_")[4] == str(alt)}
    keep, names = [], []
    for i, a in enumerate(nm[0]):
        f = a.split("_")
        if ((f[3], f[1]) in pos and f[4] == str(alt)) or ((f[3], f[1]) not in pos and f[4] == "0"):
            f[4] = "0"
            keep.append(i)
            names.append("_".join(f))
    return cm[0, keep].unsqueeze(0), [names]


# ---------------------------------------------------------------- wild type untouched by mutants
@pytest.mark.xfail_upstream("WT atoms take partners from the mutant")
@pytest.mark.parametrize("mut", ["155_A_GLN", "155_A_HIS", "163_A_ASN", "155_A_GLU"])
def test_wt_unchanged_by_mutant(parse, energy, mut):
    wt = energy(*parse("small"))[..., 0, :]
    e = energy(*_mut(parse, [[mut]]))
    assert_energy_close(e[..., 0, :], wt, what="WT row with %s present" % mut)


@pytest.mark.slow
@pytest.mark.xfail_upstream("WT atoms take partners from the mutant")
def test_wt_unchanged_full_scan(parse, energy):
    wt = energy(*parse("small"))[..., 0, :]
    for p in POSITIONS:
        for a in NON_LYS + ["LYS"]:
            e = energy(*_mut(parse, [["%d_A_%s" % (p, a)]]))
            assert_energy_close(e[..., 0, :], wt, what="WT row with %d->%s present" % (p, a))


@pytest.mark.xfail_upstream("a LYS mutant's value overwrites the WT peptide-bond term (compiled BondLenConstrain)")
@pytest.mark.parametrize("pos", [155, 163, 178, 190, 201])
def test_wt_unchanged_by_lys_mutant(parse, energy, pos):
    wt = energy(*parse("small"))[..., 0, :]
    e = energy(*_mut(parse, [["%d_A_LYS" % pos]]))
    assert_energy_close(e[..., 0, :], wt)


# ---------------------------------------------------------------- each mutant is right on its own
@pytest.mark.parametrize("mut", ["178_A_GLU", "163_A_LYS", "201_A_TRP", "155_A_PRO"])
def test_mutant_matches_standalone_build(parse, energy, mut):
    """Guard: alternative k == the same mutant built as a plain structure (true upstream too)."""
    cm, nm = _mut(parse, [[mut]])
    e = energy(cm, nm)[..., 1, :]
    es = energy(*_standalone(cm, nm, 1))[..., 0, :]
    assert_energy_close(e, es, atol=1e-4)


def _atoms(c, n, alt):
    return {a.split("_")[0] + ":" + a.split("_")[2]: c[0, i] for i, a in enumerate(n[0]) if a.split("_")[4] == str(alt)}


@pytest.mark.xfail_upstream("a mutant reuses atoms built for an earlier mutant at the same position")
def test_mutant_coordinates_order_independent(parse):
    order1 = [["178_A_LYS"], ["178_A_GLU"], ["178_A_TRP"]]
    c1, n1 = _mut(parse, order1)
    c2, n2 = _mut(parse, order1[::-1])
    for a1, a2 in [(1, 3), (2, 2), (3, 1)]:
        x, y = _atoms(c1, n1, a1), _atoms(c2, n2, a2)
        assert x.keys() == y.keys()
        for k in x:
            assert torch.allclose(x[k], y[k], atol=1e-5), "%s moved by %.3f A" % (k, (x[k] - y[k]).norm())


@pytest.mark.xfail_upstream("a mutant reuses atoms built for an earlier mutant at the same position")
def test_mutant_energy_order_independent(parse, energy):
    order1 = [["178_A_GLN"], ["178_A_GLU"], ["178_A_TRP"]]
    e1 = energy(*_mut(parse, order1))
    e2 = energy(*_mut(parse, order1[::-1]))
    for a1, a2 in [(1, 3), (2, 2), (3, 1)]:
        assert_energy_close(e1[..., a1, :], e2[..., a2, :], atol=1e-4, what="alt %d vs %d" % (a1, a2))


@pytest.mark.xfail_upstream("a LYS alternative shifts other alternatives' peptide-bond term (compiled BondLenConstrain)")
def test_mutant_energy_order_independent_with_lys(parse, energy):
    order1 = [["178_A_LYS"], ["178_A_GLU"], ["178_A_TRP"]]
    e1 = energy(*_mut(parse, order1))
    e2 = energy(*_mut(parse, order1[::-1]))
    for a1, a2 in [(1, 3), (2, 2), (3, 1)]:
        assert_energy_close(e1[..., a1, :], e2[..., a2, :], atol=1e-4)


@pytest.mark.xfail_upstream("second mutate() call reuses alternative index 1")
def test_mutate_twice_gives_new_alternatives(parse):
    c, n = _mut(parse, [["178_A_GLU"]])
    c2, n2 = mutatingEngine.mutate(c, n, [[["163_A_TYR"]]])   # TYR: kept, built and ring-centroid atoms
    by_alt = {}
    for a in n2[0]:
        r, num, at, ch, mut, _ = a.split("_")
        by_alt.setdefault(int(mut), set()).add((r, int(num), at))
    assert sorted(by_alt) == [0, 1, 2]
    # every atom of each new residue carries that residue's alternative index
    assert by_alt[1] == {("GLU", 178, at) for at in mutatingEngine.atomsOfEachresi["GLU"]}
    assert by_alt[2] == {("TYR", 163, at) for at in mutatingEngine.atomsOfEachresi["TYR"]}
    info = dataStructures.create_info_tensors(n2)
    assert info[5].shape[-1] == 3


@pytest.mark.xfail_upstream("raise('...') raises TypeError")
def test_mutation_list_length_error(parse):
    c, n = parse("small")
    with pytest.raises(ValueError):
        mutatingEngine.mutate(c, n, [[["178_A_GLU"]], [["163_A_ASN"]]])


# ---------------------------------------------------------------- double mutants (one alternative, two residues)
DOUBLES = [["163_A_LYS", "178_A_GLU"], ["155_A_GLN", "190_A_TRP"], ["170_A_ALA", "201_A_PHE"]]


@pytest.mark.parametrize("double", DOUBLES, ids=["+".join(d) for d in DOUBLES])
def test_double_mutant_matches_standalone_build(parse, energy, double):
    """Guard: a double-mutant alternative equals the same double mutant built as a plain structure."""
    cm, nm = _mut(parse, [double])
    assert {a.split("_")[1] for a in nm[0] if a.split("_")[4] == "1"} == {m.split("_")[0] for m in double}
    e = energy(cm, nm)[..., 1, :]
    es = energy(*_standalone(cm, nm, 1))[..., 0, :]
    assert_energy_close(e, es, atol=1e-4)


@pytest.mark.xfail_upstream("WT row changes when (LYS-containing or partner-sharing) double mutants are present")
def test_wt_unchanged_by_double_mutants(parse, energy):
    wt = energy(*parse("small"))[..., 0, :]
    e = energy(*_mut(parse, DOUBLES))
    assert_energy_close(e[..., 0, :], wt, what="WT row with %d double mutants present" % len(DOUBLES))
    for k, d in enumerate(DOUBLES, start=1):     # and each double is unaffected by the others in the call
        alone = energy(*_mut(parse, [d]))[..., 1, :]
        assert_energy_close(e[..., k, :], alone, atol=1e-4, what="double %s next to others" % "+".join(d))
