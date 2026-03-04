#include "field_basics.h"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace {
int PolyDegree(uint32_t p) {
  if (p == 0u) return -1;
  for (int b = 31; b >= 0; --b) {
    if (((p >> b) & 1u) != 0u) return b;
  }
  return -1;
}

uint32_t PolyMod(uint32_t a, uint32_t b) {
  if (b == 0u) return a;
  const int deg_b = PolyDegree(b);
  int deg_a = PolyDegree(a);
  while (deg_a >= deg_b && deg_a >= 0) {
    a ^= (b << (deg_a - deg_b));
    deg_a = PolyDegree(a);
  }
  return a;
}

bool IsIrreducibleOverGF2(uint32_t poly, int n) {
  if (n <= 0 || n >= 31) return false;
  if (((poly >> n) & 1u) == 0u) return false;
  if ((poly & 1u) == 0u) return false;
  if (PolyDegree(poly) != n) return false;

  for (int d = 1; d <= n / 2; ++d) {
    const uint32_t begin = (1u << d);
    const uint32_t end = (1u << (d + 1));
    for (uint32_t g = begin; g < end; ++g) {
      if ((g & 1u) == 0u) continue;
      if (PolyMod(poly, g) == 0u) return false;
    }
  }
  return true;
}

struct PrimePowerFactor {
  uint32_t prime = 0u;
  uint32_t exponent = 0u;
};

std::vector<PrimePowerFactor> FactorizePrimePowers(uint32_t value) {
  std::vector<PrimePowerFactor> factors;
  if (value <= 1u) return factors;

  uint32_t x = value;
  for (uint32_t p = 2u; p * p <= x; ++p) {
    if (x % p != 0u) continue;
    uint32_t e = 0u;
    while (x % p == 0u) {
      x /= p;
      ++e;
    }
    factors.push_back(PrimePowerFactor{p, e});
  }
  if (x > 1u) factors.push_back(PrimePowerFactor{x, 1u});
  return factors;
}
}  // namespace

uint32_t GF2n::Mul(uint32_t a, uint32_t b) const {
  uint32_t res = 0;
  uint32_t x = a;
  uint32_t y = b;
  for (int i = 0; i < n; ++i) {
    if ((y & 1u) != 0u) res ^= x;
    y >>= 1u;
    const bool carry = ((x >> (n - 1)) & 1u) != 0u;
    x <<= 1u;
    if (carry) x ^= mod_poly;
    x &= mask;
  }
  return res & mask;
}

uint32_t GF2n::Pow(uint32_t a, int e) const {
  uint32_t result = 1u;
  uint32_t base = a;
  int exp = e;
  while (exp > 0) {
    if (exp & 1) result = Mul(result, base);
    base = Mul(base, base);
    exp >>= 1;
  }
  return result;
}

uint32_t PowerMap(const GF2n& f, uint32_t x, int exponent) {
  if (x == 0u) return 0u;
  const int mod = (1 << f.n) - 1;
  int e = exponent % mod;
  if (e < 0) e += mod;
  return f.Pow(x, e);
}

std::vector<uint32_t> BuildPowerTruthTable(const GF2n& f, int exponent) {
  const uint32_t q = static_cast<uint32_t>(1u << f.n);
  std::vector<uint32_t> truth_table(q, 0u);
  for (uint32_t x = 0u; x < q; ++x) {
    truth_table[x] = PowerMap(f, x, exponent);
  }
  return truth_table;
}

uint32_t FindPrimitiveElement(const GF2n& f) {
  if (f.n <= 0 || f.n >= 31) return 0u;
  const uint32_t order = (static_cast<uint32_t>(1u) << f.n) - 1u;
  if (order <= 1u) return 1u;

  const std::vector<PrimePowerFactor> factors = FactorizePrimePowers(order);
  if (factors.empty()) return 0u;

  for (uint32_t g = 2u; g <= order; ++g) {
    bool primitive = true;
    for (const PrimePowerFactor& factor : factors) {
      const uint32_t test_exp = order / factor.prime;
      if (f.Pow(g, static_cast<int>(test_exp)) == 1u) {
        primitive = false;
        break;
      }
    }
    if (primitive) return g;
  }
  return 0u;
}

uint32_t DefaultModPoly(int n) {
  if (n <= 0 || n >= 31) {
    throw std::invalid_argument("DefaultModPoly: n must be in [1, 30]");
  }

  const uint32_t lower = (1u << n) | 1u;
  const uint32_t upper = (1u << (n + 1));
  for (uint32_t poly = lower; poly < upper; poly += 2u) {
    if (IsIrreducibleOverGF2(poly, n)) return poly;
  }

  throw std::runtime_error("DefaultModPoly: no irreducible polynomial found");
}
