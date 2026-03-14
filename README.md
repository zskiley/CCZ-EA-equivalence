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

### 1) galois (finite-field polynomial)

```python
import ccz
import galois

n = 9
F = galois.GF(2**n)
x = galois.Poly.Identity(F)

f = x**3
g = x**5

auto = ccz.ccz_auto(f)
eq = ccz.ccz_equivalence(f, g)

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

### `ccz.ccz_auto(...)`

```python
ccz.ccz_auto(
    values_or_fn,
    time_limit_seconds=None,
    min_active_hyperplanes=None,
)
```

- `values_or_fn`: truth table, `galois.Poly`, or Sage polynomial.
- Input/output dimensions are inferred automatically:
  truth-table length gives `n`, and table values (or polynomial field degree)
  determine `m`.
- `time_limit_seconds`: optional auto-search timeout.
- `min_active_hyperplanes`: optional refinement budget override.

### `ccz.ea_auto(...)`

```python
ccz.ea_auto(
    values_or_fn,
    time_limit_seconds=None,
    min_active_hyperplanes=None,
)
```

Parameters are the same as `ccz.ccz_auto(...)`.

### `ccz.ccz_equivalence(...)`

```python
ccz.ccz_equivalence(
    f_values_or_fn,
    g_values_or_fn,
    time_limit_seconds=None,
    min_active_hyperplanes=None,
    auto_group=None,
)
```

- `f_values_or_fn`, `g_values_or_fn`: truth tables, `galois.Poly`, or Sage
  polynomials.
- dimensions are inferred automatically (as in `ccz.ccz_auto(...)`).
- `time_limit_seconds`, `min_active_hyperplanes`:
  same meaning as in `ccz.ccz_auto(...)`.
- `auto_group`: optional precomputed automorphism group seed.
  Accepts either:
  - a `list[dict]` of ambient affine generators, or
  - the full auto-result dict from `ccz_auto`/`ea_auto`
    (uses its `"generators"` field).

### `ccz.ea_equivalence(...)`

```python
ccz.ea_equivalence(
    f_values_or_fn,
    g_values_or_fn,
    time_limit_seconds=None,
    min_active_hyperplanes=None,
    auto_group=None,
)
```

Parameters are the same as `ccz.ccz_equivalence(...)`.

## Return formats

### `ccz_auto(...)` / `ea_auto(...)`

Returns a dictionary:

- `order: int`
  : full ambient affine automorphism-group order found by the algorithm.
- `found_entire_group: bool`
  : `True` if search finished before timeout, `False` if timeout was hit.
- `generators: list[dict]`
  : ambient affine generators.

Each generator has the form:

```python
{
    "translation": int,
    "linear_cols": list[int],
}
```

This represents the affine map

`z |-> L(z) + t`

on the ambient space `F_2^{n+m}`, where:

- `translation` is `t`
- `linear_cols[i]` is the image of basis vector `e_i` under `L`

### `ccz_equivalence(...)` / `ea_equivalence(...)`

Returns either:

- `None` (no equivalence found), or
- `dict[int, int]` mapping graph points from `F` to graph points in `G`
  using the same encoding `p = x | (y << n)`.

### Optional `auto_group=` input format for equivalence

You may pass:

- a `list[dict]` of ambient affine generators, or
- the full auto result dict from `ccz_auto`/`ea_auto` (it will use the
  `"generators"` field).

For backward compatibility, a list of graph-point generator maps is also
accepted, but the default auto output now uses ambient affine generators.
