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

## Quick usage

Supported Python inputs:
- truth tables (`list[int]` / sequence of integers)
- `galois.Poly`
- Sage polynomials

Automorphism functions return Sage matrix groups, so run them in Sage or with
`sage -python`.

### 1) Automorphism group in Sage

```python
import ccz
from sage.all import GF, PolynomialRing

n = 9
K = GF(2**n, name="a")
R = PolynomialRing(K, "x")
x = R.gen()

G, complete = ccz.ccz_auto(x**3)

print(G.order(), complete)
```

For rectangular functions, pass the explicit function tuple `(f, n, m)`:

```python
import ccz

f = [0, 1, 2, 3, 0, 1, 2, 3]  # F_2^3 -> F_2^2

G, complete = ccz.ea_auto((f, 3, 2))

print(G.order(), complete)
```

### 2) Equivalence

```python
import ccz

f_tt = [0, 1, 2, 3, 4, 5, 6, 7]
g_tt = [0, 1, 2, 3, 4, 5, 6, 7]

eq = ccz.ccz_equivalence(f_tt, g_tt)

print(eq is not None)
```

## Recommendations
If the goal is to prove that two functions are equivalent, setting
`time_limit_seconds` to a low value can be useful. If the goal is to prove
inequivalence, the automorphism group is often essential.

## Available API

### `ccz.ccz_auto(...)`

```python
ccz.ccz_auto(
    function_input,
    time_limit=None,
    min_active_hyperplanes=None,
)
```

- `function_input`: either `(f, n, m)` or a bare truth table / polynomial.
  Bare inputs infer `n` and use `m = n`.
- `time_limit`: optional auto-search timeout in seconds.
- `min_active_hyperplanes`: optional refinement budget override.
- returns `(G, complete)`, where `G` is a Sage `MatrixGroup` over `GF(2)`
  acting on homogeneous vectors of length `n + m + 1`, and `complete` reports
  whether the search finished before timeout.

### `ccz.ea_auto(...)`

```python
ccz.ea_auto(
    function_input,
    time_limit=None,
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
    parallel_auto_seed=True,
)
```

- `f_values_or_fn`, `g_values_or_fn`: truth tables, `galois.Poly`, or Sage
  polynomials.
- dimensions are inferred automatically from truth-table length or polynomial
  field degree.
- `time_limit_seconds`, `min_active_hyperplanes`:
  same meaning as in `ccz.ccz_auto(...)`.
- `auto_group`: optional precomputed automorphism group seed.
  Accepts either:
  - a Sage `MatrixGroup` returned by `ccz_auto`/`ea_auto`,
  - the `(G, complete)` tuple returned by `ccz_auto`/`ea_auto`,
  - a `list[dict[int, int]]` of graph-point generators,
  - a `list[{"translation": int, "linear_cols": list[int]}]` of ambient
    affine generators, or
  - a raw auto-result dict from the native binding.
- `parallel_auto_seed`: enabled by default for equivalence.
  When `auto_group` is not supplied, the wrapper computes automorphism seeds for
  both inputs in parallel, prefers the larger discovered seed on the right-hand
  side, and if both auto searches finish with different group orders it
  immediately concludes the functions are not equivalent.

### `ccz.ea_equivalence(...)`

```python
ccz.ea_equivalence(
    f_values_or_fn,
    g_values_or_fn,
    time_limit_seconds=None,
    min_active_hyperplanes=None,
    auto_group=None,
    parallel_auto_seed=True,
)
```

Parameters are the same as `ccz.ccz_equivalence(...)`.

## Return formats

### `ccz_auto(...)` / `ea_auto(...)`

Returns `(G, complete)`.

- `G`: Sage `MatrixGroup` over `GF(2)`.
- `complete`: `True` if search finished before timeout, `False` if timeout was
  hit.

Affine maps `z -> Lz + t` are represented as homogeneous matrices:

```text
[ L  t ]
[ 0  1 ]
```

### `ccz_equivalence(...)` / `ea_equivalence(...)`

Returns either:

- `None` (no equivalence found), or
- `dict[int, int]` mapping graph points from `F` to graph points in `G`
  using the same encoding `p = x | (y << n)`.

### Optional `auto_group=` input format for equivalence

You may pass:

- a `list[dict[int, int]]` of graph-point generators,
- a `list[{"translation": int, "linear_cols": list[int]}]` of ambient affine
  generators, or
- a Sage `MatrixGroup`, or the `(G, complete)` tuple returned by
  `ccz_auto`/`ea_auto`.

If `auto_group` is provided, it is used directly and the default parallel
auto-seeding heuristic is skipped.

### Default `parallel_auto_seed=True` heuristic for equivalence

When enabled, the wrapper:

- computes auto seeds for both inputs in separate Python subprocesses,
- compares the discovered group orders,
- swaps the equivalence direction when the left input yields the larger seed so
  that larger seed is used on the right-hand side,
- and returns `None` immediately if both complete auto searches finish and the
  resulting group orders differ.

If the wrapper swaps the search direction, it inverts the returned point map
before returning it to the caller.
