"""fakeAtomsGeneration: hydrogens / lone pairs must never be oriented from a missing (padding) partner."""
import pytest
import torch

from conftest import dataStructures
from madrax.sources import fakeAtomsGeneration, hashings
from madrax.sources.globalVariables import PADDING_INDEX


def _fake_atoms(coords, names):
    """Same steps as ForceField.forward up to generateFakeAtomTensor."""
    info = dataStructures.create_info_tensors(names)
    atom_number, ad, cidx, pidx = info[:4]
    parts = []
    for k in range(2):
        out = torch.full((ad.shape[0], 3), float(PADDING_INDEX))
        has = pidx[:, k] != PADDING_INDEX
        out[has] = coords[ad[:, 0].long()[has], pidx[:, k][has]]
        parts.append(out.unsqueeze(1))
    x = coords[ad[:, 0].long(), cidx]
    fake = fakeAtomsGeneration.generateFakeAtomTensor(x, torch.cat(parts, 1), ad, hashings.fake_atom_Properties)
    return fake, pidx, ad


@pytest.mark.xfail_master("second partner never checked: fake atoms built from the (-999,-999,-999) padding point")
@pytest.mark.parametrize("name", ["missing_partner", "truncated_sidechains"])
def test_no_fake_atoms_from_missing_partner(parse, name):
    fake, pidx, ad = _fake_atoms(*parse(name))
    missing = (pidx == PADDING_INDEX).any(-1)
    built = (fake[..., 0] != PADDING_INDEX).any(-1)
    bad = (missing & built).nonzero().flatten()
    names = [(hashings.atom_hash_inverse[int(ad[i, 4])], int(ad[i, 2])) for i in bad[:5]]
    assert len(bad) == 0, "%d atoms with a missing partner still get fake atoms, e.g. %s" % (len(bad), names)


def test_complete_residues_keep_their_fake_atoms(parse):
    """Guard: with every partner present the fix must not drop anything."""
    fake, pidx, ad = _fake_atoms(*parse("small"))
    complete = (pidx != PADDING_INDEX).all(-1)
    expected = (hashings.fake_atom_Properties[ad[:, 4].long()][..., 0] != PADDING_INDEX)
    built = fake[..., 0] != PADDING_INDEX
    assert torch.equal(built[complete], expected[complete])
    assert torch.isfinite(fake[built]).all()
