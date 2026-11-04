"""dataStructures.create_info_tensors: termini, numbering, partners, angle indices."""
import pytest
import torch

from conftest import dataStructures, pad_batch
from madrax.mutate import mutatingEngine
from madrax.sources import hashings
from madrax.sources.globalVariables import PADDING_INDEX
from reference.dense import reference_angle_indices

AD = hashings.atom_description_hash
TN = {hashings.atom_hash[r]["tN"] for r in hashings.resi_hash}
OXT = {hashings.atom_hash[r]["OXT"] for r in hashings.resi_hash}


def _scan(parse, positions=(163, 178, 190), aas=("LYS", "GLU", "TRP")):
    c, n = parse("small")
    return mutatingEngine.mutate(c, n, [[[f"{p}_A_{a}"] for p in positions for a in aas]])


@pytest.mark.xfail_upstream("terminus bookkeeping is shared across proteins in a batch")
@pytest.mark.parametrize("order", [("small", "renumbered+1000"), ("renumbered+1000", "small"),
                                   ("multichain", "small")])
def test_termini_tagged_per_protein(parse, order):
    _, names = pad_batch(*[parse(o) for o in order])
    ad = dataStructures.create_info_tensors(names)[1]
    for b in range(len(order)):
        alone = dataStructures.create_info_tensors([parse(order[b])[1][0]])[1]
        mine = ad[ad[:, AD["batch"]] == b]
        # atom typing (incl. tN / OXT tags) must not depend on the rest of the batch
        assert torch.equal(mine[:, AD["at_name"]], alone[:, AD["at_name"]]), order[b]
        for ch in mine[:, AD["chain"]].unique():
            m = mine[:, AD["chain"]] == ch
            res, at = mine[m, AD["resnum"]], mine[m, AD["at_name"]].tolist()
            tn = [i for i, a in enumerate(at) if a in TN]
            ox = [i for i, a in enumerate(at) if a in OXT]   # an explicit OXT in the file adds a second one
            assert len(tn) == 1 and res[tn[0]] == res.min(), (order[b], int(ch))
            assert ox and all(res[i] == res.max() for i in ox), (order[b], int(ch))


@pytest.mark.xfail_upstream("negative residue numbers wrap around in dense residue tensors")
def test_negative_residue_numbers_raise(parse):
    _, n = parse("renumbered-160")
    with pytest.raises(ValueError, match="negative residue"):
        dataStructures.create_info_tensors(n)


@pytest.mark.xfail_upstream("partner lookup returns the last-written alternative's atom")
def test_partners_stay_within_alternative(parse):
    _, nm = _scan(parse)
    info = dataStructures.create_info_tensors(nm)
    ad, partners = info[1], info[3]
    alt = ad[:, AD["alternative"]]
    # partner indices are positions in the input atom list; map them to atom_description rows
    row_of = {int(j): i for i, j in enumerate(info[2].tolist())}
    for k in range(2):
        p = partners[:, k]
        has = p != PADDING_INDEX
        src = has.nonzero().flatten()
        dst = torch.tensor([row_of[int(j)] for j in p[has]])
        ok = (alt[dst] == 0) | (alt[dst] == alt[src])
        assert ok.all(), "%d atoms take partner %d from another alternative" % (int((~ok).sum()), k + 1)
        # and WT atoms only ever use WT partners
        assert (alt[dst][alt[src] == 0] == 0).all()


@pytest.mark.parametrize("case", ["small", "multichain", "batch", "scan"])
def test_angle_indices_match_original_loop(parse, case):
    """Differential: vectorised angle-index build == the original per-atom loop."""
    if case == "batch":
        _, names = pad_batch(parse("small"), parse("multichain"))
    elif case == "scan":
        _, names = _scan(parse)
    else:
        names = parse(case)[1]
    info = dataStructures.create_info_tensors(names)
    ref = reference_angle_indices(names, info[5], info[1])
    assert info[4].shape == ref.shape
    assert torch.equal(info[4].long(), ref.long())


@pytest.mark.slow
def test_angle_indices_match_original_loop_large(parse):
    names = parse("large")[1]
    info = dataStructures.create_info_tensors(names)
    assert torch.equal(info[4].long(), reference_angle_indices(names, info[5], info[1]).long())
