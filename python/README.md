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
before timeout. `G` can be passed back into
`ccz.ccz_equivalence(..., right_auto=...)` or
`ccz.ea_equivalence(..., right_auto=...)`.

Equivalence functions require Sage and return either `None` or a Sage matrix
over `GF(2)` representing an affine equivalence. With `parallel=True` and no
supplied automorphism group, equivalence races left auto, right auto, and an
equivalence search in parallel. If an auto search finishes first, the other
running tasks are stopped and one seeded equivalence search is started with that
auto group.

If `left_auto` or `right_auto` is a supplied group, the wrapper skips the
parallel race and runs one seeded equivalence search. If both supplied groups
are available, the larger-order group is used as the seed.

`time_limit_auto_search` is measured in seconds and only limits automorphism
searches. It does not kill the equivalence search itself. If a supplied auto
group is passed with `left_auto=` or `right_auto=`, no auto search is run, so
the timeout has no effect for that call. The old names `time_limit` and
`time_limit_seconds` are still accepted as aliases.
