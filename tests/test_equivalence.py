"""Speedups change nothing: fast paths vs master's code (kept in tests/reference) and vs golden output from master."""
import os

import numpy as np
import pytest
import torch

from conftest import DATA, dataStructures, assert_energy_close
from madrax.mutate import mutatingEngine, rotator
from reference.dense import dense_pairwise_representation

GOLDEN = os.path.join(DATA, "golden")
# one mutant per position: master's mutant-order bug would otherwise corrupt the reference itself
SCAN = [["163_A_LYS"], ["178_A_TRP"], ["190_A_GLU"], ["201_A_ARG"], ["170_A_PHE"], ["155_A_ILE"]]


def _pairs(ff, coords, names, dense=False):
    info = dataStructures.create_info_tensors(names)
    x = coords[info[1][:, 0].long(), info[2]]
    with torch.no_grad():
        if dense:
            return dense_pairwise_representation(ff, x, info[0], info[1])
        return ff.getPairwiseRepresentation(x, info[0], info[1])


def _perturb(coords, names, seed, drop=0.05, jitter=0.5):
    g = torch.Generator().manual_seed(seed)
    is_ca = torch.tensor([a.split("_")[2] == "CA" for a in names[0]])   # every residue must keep its CA
    keep = ((torch.rand(len(names[0]), generator=g) > drop) | is_ca).nonzero().flatten()
    c = coords[:, keep] + torch.randn(1, len(keep), 3, generator=g) * jitter
    return c, [[names[0][i] for i in keep]]


# ---------------------------------------------------------------- neighbour search
@pytest.mark.parametrize("name", ["small", "multichain"])
def test_pairs_identical_to_dense(parse, ff, name):
    c, n = parse(name)
    fast, dense = _pairs(ff, c, n), _pairs(ff, c, n, dense=True)
    assert fast.shape == dense.shape and torch.equal(fast, dense)


@pytest.mark.parametrize("seed", range(20))
def test_pairs_identical_to_dense_random(parse, ff, seed):
    """Jittered (0.5 A) structures with 5% of atoms dropped: atoms sit right at the 7 A cutoff."""
    c, n = _perturb(*parse("small" if seed % 2 else "multichain"), seed)
    assert torch.equal(_pairs(ff, c, n), _pairs(ff, c, n, dense=True))


def test_pairs_identical_with_mutants(parse, ff):
    c, n = parse("small")
    cm, nm = mutatingEngine.mutate(c, n, [[["178_A_GLU"], ["178_A_LYS"], ["163_A_ASN"]]])
    assert torch.equal(_pairs(ff, cm, nm), _pairs(ff, cm, nm, dense=True))


@pytest.mark.slow
def test_pairs_identical_to_dense_large(parse, ff):
    c, n = parse("large")
    assert torch.equal(_pairs(ff, c, n), _pairs(ff, c, n, dense=True))


# ---------------------------------------------------------------- golden output from master
def _golden(name):
    path = os.path.join(GOLDEN, name + ".npz")
    if not os.path.exists(path):
        pytest.skip("golden file missing: run tools/make_golden.py against master")
    return np.load(path)


@pytest.mark.parametrize("name", ["small", "multichain", pytest.param("large", marks=pytest.mark.slow)])
def test_energy_matches_master_golden(parse, energy, name):
    g = _golden(name)
    assert_energy_close(energy(*parse(name)), torch.tensor(g["energy"]), what="%s vs master golden" % name)


def test_gradient_matches_master_golden(parse, ff):
    g = _golden("small")
    c, n = parse("small")
    c = c.clone().requires_grad_(True)
    ff(c, dataStructures.create_info_tensors(n)).sum().backward()
    ref = torch.tensor(g["grad"])
    assert torch.allclose(c.grad, ref, atol=1e-3, rtol=1e-4), "max |diff| %.3g" % (c.grad - ref).abs().max()


def test_sidechain_rotation_matches_master_golden(parse):
    """Rotator output on a mutant scan (exercises the topology cache) equals master's."""
    g = _golden("rotator_scan")
    c, n = parse("small")
    cm, nm = mutatingEngine.mutate(c, n, [SCAN])
    info = dataStructures.create_info_tensors(nm)
    rot = torch.tensor(g["rotation"])
    r = rotator.RotateStruct()
    with torch.no_grad():
        first = r(cm, info, rot, backbone_rotation=False)
        second = r(cm, info, rot, backbone_rotation=False)   # served from the cache
    assert torch.equal(first, second)
    # ILE CD1 now rotates with chi2 (fixed here, frozen on master); test_rotator checks it geometrically
    ile_cd1 = torch.tensor([a.split("_")[0] == "ILE" and a.split("_")[2] == "CD1" for a in nm[0]])
    assert ile_cd1.any()
    assert torch.allclose(first[0, ~ile_cd1], torch.tensor(g["coords"])[0, ~ile_cd1], atol=1e-5)


# ---------------------------------------------------------------- GPU
@pytest.mark.gpu
def test_cuda_matches_cpu(parse, energy):
    from madrax.ForceField import ForceField
    c, n = parse("small")
    cpu = energy(c, n)
    try:
        ffg = ForceField(device="cuda")   # moves module-global property tensors to the GPU
        with torch.no_grad():
            gpu = ffg(c.cuda(), dataStructures.create_info_tensors(n, device="cuda")).cpu()
    finally:
        ForceField(device="cpu")           # move them back for the rest of the session
    assert_energy_close(gpu, cpu, atol=1e-3)
