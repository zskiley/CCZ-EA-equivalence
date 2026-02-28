# Python/Sage usage

## Install from GitHub

```bash
pip install "git+https://github.com/<you>/<repo>.git"
```

## Manual build (developer mode)

From repository root:

```bash
python python/build.py
```

If you want Sage compatibility, build with Sage's Python:

```bash
sage -python python/build.py
```

This produces `python/ccz_bindings...` (`.pyd` on Windows, `.so` on Linux/macOS).

## Use from Python or Sage

```python
import sys
sys.path.append(r"/path/to/repo/python")

import ccz

n = 9
f = lambda x: x**3
g = lambda x: x**6

autos = ccz.ccz_auto(f, n, time_limit_seconds=2.0)
eq = ccz.ccz_equivalence(f, g, n, time_limit_seconds=2.0)
print(len(autos), eq is not None)
```

### API

- `ccz.ccz_auto(values_or_fn, n_bits, m_bits=None, time_limit_seconds=None)`
- `ccz.ea_auto(values_or_fn, n_bits, m_bits=None, time_limit_seconds=None)`
- `ccz.ccz_equivalence(f_values_or_fn, g_values_or_fn, n_bits, m_bits=None, time_limit_seconds=None)`
- `ccz.ea_equivalence(f_values_or_fn, g_values_or_fn, n_bits, m_bits=None, time_limit_seconds=None)`

`values_or_fn` can be:
- a truth-table iterable of length `2^n_bits`, or
- a callable `f(x)` returning integer outputs.
