"""Physics the force field must respect whatever the implementation: batching, numbering, rigid motion."""
import pytest
import torch

from conftest import dataStructures, pad_batch, assert_energy_close


def _alone(energy, parse, name):
    return energy(*parse(name))


@pytest.mark.xfail_upstream("protein energies depend on the other proteins in the batch")
@pytest.mark.parametrize("order", [("small", "renumbered+1000"), ("renumbered+1000", "small"),
                                   ("multichain", "small"), ("small", "multichain", "renumbered-151")])
def test_batch_invariance(parse, energy, order):
    eb = energy(*pad_batch(*[parse(o) for o in order]))
    for b, name in enumerate(order):
        ea = _alone(energy, parse, name)[0]
        ch, res = ea.shape[0], ea.shape[1]
        assert_energy_close(eb[b, :ch, :res], ea, what="%s in batch %s" % (name, order))
        if eb.shape[1] > ch:
            assert eb[b, ch:].abs().max() == 0
        if eb.shape[2] > res:
            assert eb[b, :, res:].abs().max() == 0


@pytest.mark.parametrize("shift", [-151, 1000])
def test_residue_number_shift_invariance(parse, energy, shift):
    """Guard: renumbering by a constant (staying >= 0) must not change any residue's energy."""
    e0 = _alone(energy, parse, "small")[0]
    e1 = _alone(energy, parse, "renumbered%+d" % shift)[0]
    lo, hi = 152, 221
    assert_energy_close(e1[:, lo + shift:hi + shift], e0[:, lo:hi])
    assert e1[:, :lo + shift].abs().max() == 0


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_rigid_motion_invariance(parse, energy, seed):
    """Guard: rotating and translating the whole protein leaves every energy unchanged (float32 noise ~4e-4)."""
    c, n = parse("small")
    g = torch.Generator().manual_seed(seed)
    q = torch.randn(4, generator=g)
    a, b, cc, d = (q / q.norm()).tolist()
    R = torch.tensor([[a*a+b*b-cc*cc-d*d, 2*(b*cc-a*d), 2*(b*d+a*cc)],
                      [2*(b*cc+a*d), a*a-b*b+cc*cc-d*d, 2*(cc*d-a*b)],
                      [2*(b*d-a*cc), 2*(cc*d+a*b), a*a-b*b-cc*cc+d*d]])
    t = torch.randn(3, generator=g) * 20
    c2 = c.clone()
    m = c2[..., 0] != -999
    c2[m] = c2[m] @ R.T + t
    assert_energy_close(energy(c2, n), energy(c, n), atol=1e-3)


def test_gradients_finite_and_nonzero(parse, ff):
    c, n = parse("small")
    c = c.clone().requires_grad_(True)
    ff(c, dataStructures.create_info_tensors(n)).sum().backward()
    g = c.grad[c[..., 0] != -999]
    assert torch.isfinite(g).all()
    assert (g.abs().sum(-1) > 0).float().mean() > 0.9
