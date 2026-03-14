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

`generators` is a list of ambient affine generators:

```python
{
    "translation": int,
    "linear_cols": list[int],
}
```

These same generators can be passed back into
`ccz.ccz_equivalence(..., auto_group=...)` or
`ccz.ea_equivalence(..., auto_group=...)`.
