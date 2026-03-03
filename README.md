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
n = 9
mask = (1 << n) - 1

# Example table (replace with your own finite-field table if needed)
f_tt = [(x ^ (x << 1)) & mask for x in range(1 << n)]
g_tt = [(x ^ (x << 2)) & mask for x in range(1 << n)]

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
