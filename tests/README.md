# MadraX tests

```
pytest                                                  # default run (about a minute on a laptop CPU)
pytest -m "not slow"                                    # quick run
MADRAX_PATH=/path/to/upstream pytest --target=upstream  # bug tests must fail on upstream (strict xfail)
```

| layer | files | oracle |
|---|---|---|
| invariance | `test_invariance.py`, `test_mutation.py`, `test_rotator.py` | batching, renumbering, rigid motion, mutant order; dihedrals measured with numpy |
| bug regressions | `@pytest.mark.xfail_upstream` | each fails on upstream, passes here |

Markers: `xfail_upstream` (bug test), `slow`, `gpu`.

Fixtures are built in `conftest.py` from three vendored PDB entries (1A8O, 2BEG, 2XHE without its CA-less last
residue); every edge case (renumbered, truncated, altloc, insertion code) is generated from them.
