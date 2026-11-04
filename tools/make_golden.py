"""Regenerate tests/data/golden/*.npz from an UPSTREAM checkout. Never run as part of the test suite.

    MADRAX_PATH=/path/to/upstream/madrax python tools/make_golden.py

Only inputs that none of the fixed bugs touch are frozen: single proteins with residue
numbers >= 0, no altlocs, no CA-less residues, at most one mutant per position.
"""
import os
import subprocess
import sys
import warnings

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(os.path.dirname(HERE), "tests")
upstream = os.environ.get("MADRAX_PATH")
if not upstream:
    sys.exit("set MADRAX_PATH to an upstream (unpatched) madrax checkout")
sys.path.insert(0, upstream)
sys.path.insert(1, TESTS)
warnings.filterwarnings("ignore")

from madrax.ForceField import ForceField  # noqa: E402
from madrax import utils, dataStructures  # noqa: E402
from conftest import read_lines, write  # noqa: E402

out_dir = os.path.join(TESTS, "data", "golden")
os.makedirs(out_dir, exist_ok=True)
tmp = os.path.join(out_dir, "_tmp")
os.makedirs(tmp, exist_ok=True)
try:
    commit = subprocess.check_output(["git", "-C", upstream, "rev-parse", "HEAD"], text=True).strip()
except Exception:
    commit = "unknown"
meta = dict(upstream_commit=commit, torch=torch.__version__, python=sys.version.split()[0])

torch.manual_seed(0)
ff = ForceField(device="cpu")
for name, src in (("small", "1A8O"), ("multichain", "2BEG"), ("large", "2XHE_trimmed")):
    path = write(os.path.join(tmp, name + ".pdb"), read_lines(src))
    c, n = utils.parsePDB(path)[:2]
    info = dataStructures.create_info_tensors(n)
    extra = {}
    if name == "small":
        cg = c.clone().requires_grad_(True)
        ff(cg, info).sum().backward()
        extra["grad"] = cg.grad.numpy()
    with torch.no_grad():
        e = ff(c, info)
    np.savez_compressed(os.path.join(out_dir, name + ".npz"), energy=e.numpy(), **extra, **meta)
    print(name, tuple(e.shape), "total %.4f" % float(e.sum()))

for f in os.listdir(tmp):
    os.remove(os.path.join(tmp, f))
os.rmdir(tmp)
print("written to", out_dir, "from upstream", commit)
