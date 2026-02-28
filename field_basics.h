#ifndef FIELD_BASICS_H
#define FIELD_BASICS_H

#include <cstdint>
#include <vector>

struct GF2n {
  int n;
  uint32_t mod_poly;
  uint32_t mask;

  uint32_t Mul(uint32_t a, uint32_t b) const;
  uint32_t Pow(uint32_t a, int e) const;
};

uint32_t PowerMap(const GF2n& f, uint32_t x, int exponent);
std::vector<uint32_t> BuildPowerTruthTable(const GF2n& f, int exponent);

// Returns a primitive element of GF(2^n)^*.
// Returns 1 when n=1 (trivial multiplicative group), and 0 on failure.
uint32_t FindPrimitiveElement(const GF2n& f);

// Returns a deterministic irreducible polynomial of degree n over GF(2).
// Throws std::invalid_argument when n is out of supported range.
uint32_t DefaultModPoly(int n);

#endif
