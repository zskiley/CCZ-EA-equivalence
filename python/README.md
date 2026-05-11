# Python/Sage usage

## Build

From repository root:

```bash
sage -python python/build.py
```

## Example

```python
import sys
sys.path.append(r"/path/to/repo/python")

import ccz

f = [0, 1, 2, 3, 4, 5, 6, 7]
g = [0, 1, 2, 3, 4, 5, 6, 7]

G, complete = ccz.ea_auto((f, 3, 3))
eq = ccz.ea_equivalence((f, 3, 3), (g, 3, 3), right_auto=G)

print(G.order(), complete, eq is not None)
```

## Return value of `ccz_auto(...)` / `ea_auto(...)`

Automorphism functions require Sage and return:

```python
(G, complete)
```

`G` is a Sage `MatrixGroup` over `GF(2)` representing affine maps as
homogeneous matrices, and `complete` reports whether the auto search finished
before timeout. `G` or the full `(G, complete)` tuple can be passed back into
`ccz.ccz_equivalence(..., right_auto=...)` or
`ccz.ea_equivalence(..., right_auto=...)`.

Equivalence functions require Sage and return either `None` or a Sage matrix
over `GF(2)` representing an affine equivalence. By default, equivalence runs
left auto, right auto, and an unseeded equivalence search in parallel; as auto
groups finish, additional seeded equivalence searches are started.
