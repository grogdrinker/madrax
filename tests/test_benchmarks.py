"""Timing and memory. Never part of the PR gate:  pytest -m benchmark

Each measurement runs in a fresh subprocess so peak RSS is honest. Results go to
bench_results.json; if tests/data/bench_baseline.json exists (same runner), a run more than
2x slower than its baseline fails.
"""
import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
pytestmark = pytest.mark.benchmark

SCRIPT = r'''
import os, sys, time, resource, json, warnings
sys.path.insert(0, os.environ["MADRAX_ROOT"]); sys.path.insert(1, os.environ["TESTS"])
warnings.filterwarnings("ignore")
import torch
torch.set_num_threads(int(os.environ.get("BENCH_THREADS", "2")))
from madrax.ForceField import ForceField
from madrax import utils, dataStructures
from conftest import read_lines, write
what, copies = sys.argv[1], int(sys.argv[2])
lines = read_lines("2XHE_trimmed")
if copies > 1:   # tile copies side by side as extra chains, 200 A apart
    out = []
    for k in range(copies):
        for l in lines:
            ch = chr(ord("A") + 2 * k + (0 if l[21] == "A" else 1))
            out.append(l[:21] + ch + l[22:30] + "%8.3f" % (float(l[30:38]) + 200 * k) + l[38:])
    lines = out
path = write(os.path.join(os.environ["TMP"], "bench_%d.pdb" % copies), lines)
c, n = utils.parsePDB(path)[:2]
info = dataStructures.create_info_tensors(n)
ff = ForceField()
base = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
t = time.time()
with torch.no_grad():
    if what == "pairs":
        x = c[info[1][:, 0].long(), info[2]]
        ff.getPairwiseRepresentation(x, info[0], info[1])
    else:
        ff(c, info)
dt = time.time() - t
print(json.dumps({"atoms": len(n[0]), "seconds": dt,
                  "peak_rss_gb": (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss - base) / 1e6}))
'''


def _run(what, copies, tmp_path):
    env = dict(os.environ, MADRAX_ROOT=os.environ.get("MADRAX_PATH", os.path.dirname(HERE)),
               TESTS=HERE, TMP=str(tmp_path))
    out = subprocess.run([sys.executable, "-c", SCRIPT, what, str(copies)], env=env,
                         capture_output=True, text=True, timeout=3600)
    assert out.returncode == 0, out.stderr[-2000:]
    return json.loads(out.stdout.strip().splitlines()[-1])


def _record(key, result):
    path = "bench_results.json"
    data = json.load(open(path)) if os.path.exists(path) else {}
    data[key] = result
    json.dump(data, open(path, "w"), indent=1)
    base = os.path.join(HERE, "data", "bench_baseline.json")
    if os.path.exists(base):
        ref = json.load(open(base)).get(key)
        if ref:
            assert result["seconds"] < 2 * ref["seconds"], "%s: %.2fs vs baseline %.2fs" % (key, result["seconds"], ref["seconds"])


@pytest.mark.parametrize("copies", [1, 2, 4])
def test_pair_search_scaling(tmp_path, copies):
    r = _run("pairs", copies, tmp_path)
    print("pairs x%d: %d atoms %.2fs %.2f GB" % (copies, r["atoms"], r["seconds"], r["peak_rss_gb"]))
    _record("pairs_x%d" % copies, r)


def test_forward_large(tmp_path):
    r = _run("forward", 1, tmp_path)
    print("forward: %d atoms %.2fs %.2f GB" % (r["atoms"], r["seconds"], r["peak_rss_gb"]))
    _record("forward_x1", r)
