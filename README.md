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

### Automorphism + CCZ equivalence

```python
import ccz

n = 9
f = lambda x: x**3
g = lambda x: x**5

auto = ccz.ccz_auto(f, n_bits=n)
eq = ccz.ccz_equivalence(f, g, n_bits=n)

print(auto["order"], eq is not None)
```

### Equivalence test only

```python
import ccz

n = 9
f = lambda x: x**3
g = lambda x: x**6

eq = ccz.ccz_equivalence(f, g, n_bits=n)

if eq is None:
    print("Not CCZ equivalent")
else:
    print("CCZ equivalent; map size:", len(eq))
```

## Available API

- `ccz.ccz_auto(values_or_fn, n_bits=None, m_bits=None, time_limit_seconds=None, field=None, min_active_hyperplanes=None)`
- `ccz.ea_auto(values_or_fn, n_bits=None, m_bits=None, time_limit_seconds=None, field=None, min_active_hyperplanes=None)`
- `ccz.ccz_equivalence(f_values_or_fn, g_values_or_fn, n_bits=None, m_bits=None, time_limit_seconds=None, field=None, min_active_hyperplanes=None, auto_group=None)`
- `ccz.ea_equivalence(f_values_or_fn, g_values_or_fn, n_bits=None, m_bits=None, time_limit_seconds=None, field=None, min_active_hyperplanes=None, auto_group=None)`
