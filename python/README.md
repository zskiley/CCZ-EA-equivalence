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

auto = ccz.ea_auto(f)
eq = ccz.ea_equivalence(f, g, auto_group=auto)

print(auto["order"], eq is not None)
```

## Return value of `ccz_auto(...)` / `ea_auto(...)`

Returns a dictionary with:

- `order`
- `found_entire_group`
- `generators`
- `graph_generators`

`generators` is a list of ambient affine generators. Each entry has
`translation` and `linear_cols`; the columns encode the linear part by
`linear_cols[i] = A(1 << i)`.

`graph_generators` is the same discovered group action represented as graph
point maps using the encoding `p = x | (y << n)`.

These same generators can be passed back into
`ccz.ccz_equivalence(..., auto_group=...)` or
`ccz.ea_equivalence(..., auto_group=...)`. The full auto-result dict carries
both `generators` and `graph_generators`; equivalence will use whichever seed
data is available.

By default, equivalence also runs a parallel auto-seeding heuristic: it
computes auto seeds for both inputs, prefers the larger discovered seed on the
right-hand side, and if both auto groups finish with different orders it
returns `None` immediately.
