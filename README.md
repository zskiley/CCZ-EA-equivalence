# CCZ/EA Equivalence

Algorithms for CCZ and EA automorphism groups and equivalence of vectorial
Boolean functions.

The public Python API is Sage-first: automorphism groups are returned as Sage
matrix groups, and equivalences are returned as Sage matrices over `GF(2)`.

## Build

From a source checkout:

```bash
sage -python python/build.py
```

The public `ccz_auto`, `ea_auto`, `ccz_equivalence`, and `ea_equivalence`
functions require Sage at runtime because they return Sage objects.

## Function Inputs

The preferred input format is an explicit tuple:

```python
F = (f, n, m)
```

where `f` is the function, `n` is the input dimension, and `m` is the output
dimension. This represents a function

```text
F : F_2^n -> F_2^m
```

The function `f` may be:

- a truth table, as a sequence of `2^n` non-negative integers,
- a Sage polynomial/callable over a binary field,
- a `galois.Poly`.

Bare function inputs are also accepted for convenience:

```python
ccz.ccz_auto(f)
```

For bare inputs, the wrapper infers `n` and assumes `m = n`. Use `(f, n, m)`
for rectangular functions.

Truth-table values are encoded as integers. The graph point for `x` is encoded
internally as:

```text
p = x | (f[x] << n)
```

## Automorphism Groups

```python
G, complete = ccz.ccz_auto(F, time_limit=None, min_active_hyperplanes=None)
G, complete = ccz.ea_auto(F, time_limit=None, min_active_hyperplanes=None)
```

`G` is a Sage `MatrixGroup` over `GF(2)`. It acts on homogeneous ambient
vectors of length `n + m + 1`.

Each affine generator

```text
z -> L z + t
```

is represented as the homogeneous matrix

```text
[ L  t ]
[ 0  1 ]
```

`complete` is `True` if the automorphism search finished before the timeout.
If `complete` is `False`, `G` is only the subgroup found before timeout.

Example:

```python
import ccz
from sage.all import GF, PolynomialRing

n = 8
K = GF(2**n, name="a")
R = PolynomialRing(K, "x")
x = R.gen()

G, complete = ccz.ccz_auto((x**3, n, n), time_limit=20)

print(G.order())
print(complete)
print(G.gens())
```

Rectangular truth-table example:

```python
import ccz

f = [0, 1, 2, 3, 0, 1, 2, 3]  # F_2^3 -> F_2^2

G, complete = ccz.ea_auto((f, 3, 2))
```

## Equivalence

```python
T = ccz.ccz_equivalence(
    F,
    H,
    parallel=True,
    left_auto=True,
    right_auto=True,
    time_limit=None,
    min_active_hyperplanes=None,
    verbose=True,
)

T = ccz.ea_equivalence(
    F,
    H,
    parallel=True,
    left_auto=True,
    right_auto=True,
    time_limit=None,
    min_active_hyperplanes=None,
    verbose=True,
)
```

The return value is either:

- `None`, if no equivalence is found, or
- a Sage matrix over `GF(2)` representing one affine equivalence as a
  homogeneous matrix.

Example:

```python
import ccz

f = [0, 1, 2, 3, 4, 5, 6, 7]
g = [0, 1, 2, 3, 4, 5, 6, 7]

T = ccz.ccz_equivalence((f, 3, 3), (g, 3, 3))

if T is not None:
    print("Equivalent")
    print(T)
```

## Parallel Equivalence Behavior

By default, equivalence starts three searches in parallel:

- an unseeded equivalence search,
- the left automorphism group search,
- the right automorphism group search.

As automorphism groups finish, the wrapper starts additional seeded equivalence
searches. If the complete left and right automorphism groups finish with
different orders, the wrapper returns `None`.

The auto-search switches can be booleans or precomputed groups:

```python
left = True        # compute the left auto group
right = False      # skip the right auto group
right = G          # use a supplied Sage MatrixGroup
right = (G, True)  # use the tuple returned by ccz_auto / ea_auto
```

Examples:

```python
# Run only the equivalence search.
T = ccz.ccz_equivalence(F, H, left_auto=False, right_auto=False)

# Precompute and reuse an automorphism group.
AutH, complete = ccz.ccz_auto(H)
T = ccz.ccz_equivalence(F, H, right_auto=AutH)

# Disable subprocess orchestration.
T = ccz.ccz_equivalence(F, H, parallel=False)
```

## Parameters

### `time_limit`

`time_limit` is measured in seconds and currently applies to automorphism
searches. `None` means no auto-search timeout.

The native equivalence DFS is not yet time-limited.

### `min_active_hyperplanes`

Advanced refinement parameter. Leave it as `None` unless benchmarking or
debugging the search.

### `verbose`

When `True`, the equivalence wrapper prints progress for started/finished
parallel tasks. Set `verbose=False` for quiet library use.

## Public API

```python
ccz.ccz_auto(F, time_limit=None, min_active_hyperplanes=None)
ccz.ea_auto(F, time_limit=None, min_active_hyperplanes=None)

ccz.ccz_equivalence(
    F,
    H,
    parallel=True,
    left_auto=True,
    right_auto=True,
    time_limit=None,
    min_active_hyperplanes=None,
    verbose=True,
)

ccz.ea_equivalence(
    F,
    H,
    parallel=True,
    left_auto=True,
    right_auto=True,
    time_limit=None,
    min_active_hyperplanes=None,
    verbose=True,
)
```
