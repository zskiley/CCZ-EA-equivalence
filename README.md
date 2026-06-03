# CCZ/EA Equivalence

Algorithms for CCZ and EA automorphism groups and equivalence of vectorial
Boolean functions. From the paper "Efficiently deciding and recovering CCZ and EA equivalence for arbitrary vectorial Boolean functions using the partition refinement framework".

The public Python API is Sage-first: automorphism groups are returned as Sage
matrix groups, and equivalences are returned as Sage matrices over `GF(2)`.

## Build

From a source checkout:

```bash
sage -python python/build.py
```

## Example Usage

```python
import ccz

f = [0, 1, 2, 3, 4, 5, 6, 7]
g = [0, 1, 2, 3, 4, 5, 6, 7]

T = ccz.ccz_equivalence(f, g)
print(T is not None)
```

```python
import ccz

f = [0, 1, 2, 3, 4, 5, 6, 7]
g = [0, 1, 2, 3, 4, 5, 6, 7]

AutG, complete = ccz.ccz_auto(g)

print(AutG.order())
print(complete)

T = ccz.ccz_equivalence(f, g, right_auto=AutG)
print(T is not None)
```
Note that there is some overhead when passing in the group like this, construction the sage MatrixGroup for example. It is usually more efficent to just let the equivalence algorithm compute it, unless one of the functions stays constant between tests. It is of course also possible to use the c++ code directly without sage. 

## Function Inputs

For truth tables and rectangular functions, use an explicit tuple:

```python
F = (f, n, m)
```

where `f` is the function, `n` is the input dimension, and `m` is the output
dimension. This represents a function

```text
F : F_2^n -> F_2^m
```

Here `f` is a truth table, as a sequence of `2^n` non-negative integers.

Square truth tables can also be passed directly:

```python
ccz.ccz_auto(f)
```

For bare truth tables, the wrapper infers `n` and assumes `m = n`. Use
`(f, n, m)` for rectangular truth tables.

Sage polynomials are passed directly:

```python
from sage.all import GF, PolynomialRing

K = GF(2**8, name="a")
R = PolynomialRing(K, "x")
x = R.gen()

G, complete = ccz.ccz_auto(x**3)
```

Truth-table values are encoded as integers. The graph point for `x` is encoded
internally as:

```text
p = x | (f[x] << n)
```

## Automorphism Groups

```python
G, complete = ccz.ccz_auto(
    F, time_limit_auto_search=None, min_active_hyperplanes=None
)
G, complete = ccz.ea_auto(
    F, time_limit_auto_search=None, min_active_hyperplanes=None
)
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
    left_auto=False,
    right_auto=True,
    time_limit_auto_search=None,
    min_active_hyperplanes=None,
    verbose=True,
)

T = ccz.ea_equivalence(
    F,
    H,
    parallel=True,
    left_auto=False,
    right_auto=True,
    time_limit_auto_search=None,
    min_active_hyperplanes=None,
    verbose=True,
)
```

To reuse a precomputed automorphism group, pass the Sage group directly:

```python
AutH, complete = ccz.ccz_auto(H)
T = ccz.ccz_equivalence(F, H, right_auto=AutH)
```

The return value is either:

- `None`, if no equivalence is found, or
- a Sage matrix over `GF(2)` representing one affine equivalence as a
  homogeneous matrix.

## Parallel Equivalence Behavior

With `parallel=True` and no supplied automorphism group, equivalence starts
three searches in parallel:

- an equivalence search,
- the left automorphism group search,
- the right automorphism group search.

The first successful task controls the next step. If the equivalence search
finishes first, the wrapper stops the auto searches and returns that result. If
one of the auto searches finishes first, the wrapper stops the other running
tasks and starts one new equivalence search seeded by the finished auto group.

If `left_auto` or `right_auto` is a supplied group, the wrapper skips the
parallel race and runs one seeded equivalence search. If both supplied groups
are available, the larger-order group is used as the seed.

With `parallel=False`, the auto-search switches control the sequential path:

```python
left_auto = False  # do not compute/use the left auto group
right_auto = True  # compute the right auto group
right_auto = G     # use a supplied Sage MatrixGroup
```

Examples:

```python
# Run only the equivalence search.
T = ccz.ccz_equivalence(F, H, parallel=False, right_auto=False)

# Disable subprocess orchestration.
T = ccz.ccz_equivalence(F, H, parallel=False)
```

## Parameters

### `time_limit_auto_search`

`time_limit_auto_search` is measured in seconds and applies only to
automorphism searches. `None` means no auto-search timeout.

The native equivalence DFS is not yet time-limited.

In `ccz_auto(...)` and `ea_auto(...)`, it limits that one auto search.

In `ccz_equivalence(...)` and `ea_equivalence(...)`, it limits only auto
searches started by the wrapper:

- with `parallel=True` and no supplied group, it is used by the left and right
  auto worker processes;
- with `parallel=False`, it is used only if `left_auto=True` or
  `right_auto=True`;
- with a supplied `left_auto` or `right_auto` group, no auto search is run, so
  this parameter has no effect.

The old names `time_limit` and `time_limit_seconds` are still accepted as
compatibility aliases.

### `min_active_hyperplanes`

Advanced refinement parameter. Leave it as `None` unless benchmarking or
debugging the search.

### `verbose`

When `True`, the equivalence wrapper prints progress for started/finished
parallel tasks. Set `verbose=False` for quiet library use.

## Public API

```python
ccz.ccz_auto(F, time_limit_auto_search=None, min_active_hyperplanes=None)
ccz.ea_auto(F, time_limit_auto_search=None, min_active_hyperplanes=None)

ccz.ccz_equivalence(
    F,
    H,
    parallel=True,
    left_auto=False,
    right_auto=True,
    time_limit_auto_search=None,
    min_active_hyperplanes=None,
    verbose=True,
)

ccz.ea_equivalence(
    F,
    H,
    parallel=True,
    left_auto=False,
    right_auto=True,
    time_limit_auto_search=None,
    min_active_hyperplanes=None,
    verbose=True,
)
```
