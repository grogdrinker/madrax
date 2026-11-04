"""utils.writepdb: paths and round trip."""
import os
import warnings

import pytest
import torch

from conftest import utils
from madrax.mutate import mutatingEngine


def _by_atom(coords, names):
    out = {}
    for i, a in enumerate(names[0]):
        r, num, at, ch = a.split("_")[:4]
        if at != "RC":
            out[(ch, int(num), r, at)] = coords[0, i]
    return out


@pytest.mark.parametrize("trailing_slash", [True, pytest.param(False, marks=pytest.mark.xfail_master(
    "folder and file name concatenated without a separator"))])
def test_writepdb_round_trip(parse, tmp_path, trailing_slash):
    c, n = parse("multichain")
    folder = str(tmp_path / "out") + ("/" if trailing_slash else "")
    utils.writepdb(c, n, pdb_names=["multi"], output_folder=folder)
    path = os.path.join(str(tmp_path / "out"), "multi.pdb")
    assert os.path.exists(path), os.listdir(str(tmp_path))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c2, n2 = utils.parsePDB(path)[:2]
    a, b = _by_atom(c, n), _by_atom(c2, n2)
    assert a.keys() == b.keys()
    worst = max(float((a[k] - b[k]).abs().max()) for k in a)
    assert worst <= 6e-4, "coordinates drift by %.4f A after a round trip" % worst


def test_writepdb_mutants_named(parse, tmp_path):
    c, n = parse("small")
    cm, nm = mutatingEngine.mutate(c, n, [[["178_A_GLU"]]])
    utils.writepdb(cm, nm, pdb_names=["s"], output_folder=str(tmp_path) + "/")
    assert sorted(os.listdir(str(tmp_path))) == ["s.pdb", "s_SERA178GLU.pdb"]


@pytest.mark.xfail_master("os.system('mkdir -p ' + folder) splits a folder name containing a space")
def test_writepdb_folder_with_space(parse, tmp_path):
    c, n = parse("small")
    folder = str(tmp_path / "my results") + "/"
    utils.writepdb(c, n, pdb_names=["s"], output_folder=folder)
    assert os.listdir(folder) == ["s.pdb"]
    assert not os.path.exists(str(tmp_path / "my"))
