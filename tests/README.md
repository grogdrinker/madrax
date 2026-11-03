# MadraX tests

```
pytest                                   # PR gate: everything except gpu / benchmark (~2 min on 2 CPU cores)
pytest -m "not slow"                     # quick run (~1 min)
MADRAX_PATH=/path/to/upstream pytest --target=master   # bug tests must fail upstream (strict xfail)
```

Test layers

| layer | files | oracle |
|---|---|---|
| invariance | `test_invariance.py`, `test_mutation.py`, `test_rotator.py` | physics: batching, renumbering, rigid motion, mutant order; dihedrals measured with numpy |
| bug regressions | `@pytest.mark.xfail_master` | each fails on upstream, passes here |

Markers: `xfail_master` (bug test), `known_gap` (documented, not fixed by this PR; strict xfail everywhere),
`slow`, `gpu`, `benchmark`.

Fixtures are built in `conftest.py` from three vendored PDB entries (1A8O, 2BEG, 2XHE without its CA-less last
residue); every edge case (renumbered, truncated, altloc, insertion code, hydrogens, missing partner) is generated
from them.
