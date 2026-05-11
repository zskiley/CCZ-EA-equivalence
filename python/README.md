# Python/Sage usage

## Build

From repository root:

```bash
python python/build.py
```

For Sage compatibility:

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
eq = ccz.ea_equivalence(f, g, auto_group=G)

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
`ccz.ccz_equivalence(..., auto_group=...)` or
`ccz.ea_equivalence(..., auto_group=...)`.

By default, equivalence also runs a parallel auto-seeding heuristic: it
computes auto seeds for both inputs, prefers the larger discovered seed on the
right-hand side, and if both auto groups finish with different orders it
returns `None` immediately.
