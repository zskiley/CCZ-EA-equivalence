# ccz

CCZ/EA automorphism and equivalence algorithms with Python bindings.

## Install directly from GitHub

```bash
pip install "git+https://github.com/zskiley/CCZ-EA-equivalence.git"
```

If you install from a specific branch/tag:

```bash
pip install "git+https://github.com/zskiley/CCZ-EA-equivalence.git@main"
```

## Quick Python usage

Supported Python inputs:
- truth tables (`list[int]` / sequence of integers)
- `galois.Poly`
- Sage polynomials

Not supported:
- plain Python callables/lambdas like `lambda x: x**3`

### 1) galois (finite-field polynomial)

```python
import ccz
import galois

n = 9
F = galois.GF(2**n)
x = galois.Poly.Identity(F)

f = x**3
g = x**5

auto = ccz.ccz_auto(f, field=F)
eq = ccz.ccz_equivalence(f, g, field=F)

print(auto["order"], eq is not None)
```

### 2) Truth table

```python
import ccz

# n = 3 truth tables (length = 2^3 = 8)
f_tt = [0, 1, 2, 3, 4, 5, 6, 7]
g_tt = [0, 1, 2, 3, 4, 5, 6, 7]

auto = ccz.ccz_auto(f_tt)
eq = ccz.ccz_equivalence(f_tt, g_tt)

print(auto["order"], eq is not None)
```

### 3) SageMath polynomial

```python
import ccz
from sage.all import GF, PolynomialRing

n = 9
K = GF(2**n, name="a")
R = PolynomialRing(K, "x")
x = R.gen()

f = x**3
g = x**5

auto = ccz.ccz_auto(f)
eq = ccz.ccz_equivalence(f, g)

print(auto["order"], eq is not None)
```

## Available API

- `ccz.ccz_auto(values_or_fn, n_bits=None, m_bits=None, time_limit_seconds=None, field=None, min_active_hyperplanes=None)`
- `ccz.ea_auto(values_or_fn, n_bits=None, m_bits=None, time_limit_seconds=None, field=None, min_active_hyperplanes=None)`
- `ccz.ccz_equivalence(f_values_or_fn, g_values_or_fn, n_bits=None, m_bits=None, time_limit_seconds=None, field=None, min_active_hyperplanes=None, auto_group=None)`
- `ccz.ea_equivalence(f_values_or_fn, g_values_or_fn, n_bits=None, m_bits=None, time_limit_seconds=None, field=None, min_active_hyperplanes=None, auto_group=None)`

## Return formats

### `ccz_auto(...)` / `ea_auto(...)`

Returns a dictionary:

- `order: int`
  : current group order found by the algorithm.
- `found_entire_group: bool`
  : `True` if search finished before timeout, `False` if timeout was hit.
- `generators: list[dict[int, int]]`
  : generator maps on graph-point integers.

Each generator is a dictionary `p -> q`, where `p` and `q` are graph points
encoded as integers in `F_2^{n+m}`:

`p = x | (y << n)`

with low `n` bits as input part `x`, and upper bits as output part `y`.

### `ccz_equivalence(...)` / `ea_equivalence(...)`

Returns either:

- `None` (no equivalence found), or
- `dict[int, int]` mapping graph points from `F` to graph points in `G`
  using the same encoding `p = x | (y << n)`.

### Optional `auto_group=` input format for equivalence

You may pass:

- a `list[dict[int, int]]` of generator maps, or
- the full auto result dict from `ccz_auto`/`ea_auto` (it will use the
  `"generators"` field).
