#include "weighted_fft.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace {
bool IsValidDimension(int dimension_bits) {
  return dimension_bits >= 0 && dimension_bits < 31;
}

uint64_t AbsToU64(int64_t v) {
  if (v >= 0) return static_cast<uint64_t>(v);
  return static_cast<uint64_t>(-(v + 1)) + 1ull;
}

uint32_t Parity32(uint32_t v) {
  // Return popcount(v) mod 2 without relying on C++20 <bit>::popcount.
  v ^= (v >> 16);
  v ^= (v >> 8);
  v ^= (v >> 4);
  v &= 0xFu;
  return static_cast<uint32_t>((0x6996u >> v) & 1u);
}

bool FWHTInPlaceInt64(std::vector<int64_t>* values) {
  if (values == nullptr) return false;
  const std::size_t n = values->size();
  if (!IsPowerOfTwo(n)) return false;

  int64_t* data = values->data();
  for (std::size_t len = 1; len < n; len <<= 1) {
    const std::size_t step = (len << 1);
    for (std::size_t i = 0; i < n; i += step) {
      int64_t* a = data + i;
      int64_t* b = a + len;
      for (std::size_t j = 0; j < len; ++j) {
        const int64_t u = a[j];
        const int64_t v = b[j];
        a[j] = u + v;
        b[j] = u - v;
      }
    }
  }

  return true;
}

bool FWHTInPlaceInt32(std::vector<int32_t>* values) {
  if (values == nullptr) return false;
  const std::size_t n = values->size();
  if (!IsPowerOfTwo(n)) return false;

  int32_t* data = values->data();
  for (std::size_t len = 1; len < n; len <<= 1) {
    const std::size_t step = (len << 1);
    for (std::size_t i = 0; i < n; i += step) {
      int32_t* a = data + i;
      int32_t* b = a + len;
      for (std::size_t j = 0; j < len; ++j) {
        const int32_t u = a[j];
        const int32_t v = b[j];
        a[j] = u + v;
        b[j] = u - v;
      }
    }
  }
  return true;
}
}  // namespace

bool IsPowerOfTwo(std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

bool FWHTInPlace(std::vector<int64_t>* values) {
  if (values == nullptr) return false;
  const std::size_t n = values->size();
  if (!IsPowerOfTwo(n)) return false;
  return FWHTInPlaceInt64(values);
}

bool BuildWeightedVector(int dimension_bits,
                         const std::vector<WeightedPoint>& weighted_points,
                         std::vector<int64_t>* values) {
  if (values == nullptr) return false;
  if (!IsValidDimension(dimension_bits)) return false;
  const std::size_t n = static_cast<std::size_t>(1u << dimension_bits);
  values->assign(n, 0);

  int64_t* out = values->data();
  for (const auto& entry : weighted_points) {
    const uint32_t x = entry.first;
    const int64_t w = entry.second;
    if (x >= n) return false;
    out[x] += w;
  }

  return true;
}

void WeightedFWHT(int dimension_bits,
                  const std::vector<WeightedPoint>& weighted_points,
                  std::vector<int64_t>* spectrum) {
  (void)dimension_bits;

  uint64_t sum_abs = 0;
  for (const auto& entry : weighted_points) {
    sum_abs += AbsToU64(entry.second);
  }

  if (sum_abs <= static_cast<uint64_t>(std::numeric_limits<int32_t>::max())) {
    const std::size_t n = static_cast<std::size_t>(1u << dimension_bits);
    thread_local std::vector<int32_t> buf32;
    if (buf32.size() != n) {
      buf32.assign(n, 0);
    } else {
      std::fill(buf32.begin(), buf32.end(), 0);
    }

    for (const auto& entry : weighted_points) {
      const uint32_t x = entry.first;
      buf32[x] += static_cast<int32_t>(entry.second);
    }

    FWHTInPlaceInt32(&buf32);

    if (spectrum->size() != n) spectrum->resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      (*spectrum)[i] = static_cast<int64_t>(buf32[i]);
    }
    return;
  }

  BuildWeightedVector(dimension_bits, weighted_points, spectrum);
  FWHTInPlace(spectrum);
}

void WeightedWalshAt(int dimension_bits,
                     const std::vector<WeightedPoint>& weighted_points,
                     const std::vector<uint32_t>& queries,
                     std::vector<int64_t>* values_at_queries) {
  (void)dimension_bits;
  values_at_queries->assign(queries.size(), 0);
  int64_t* out = values_at_queries->data();

  for (std::size_t qi = 0; qi < queries.size(); ++qi) {
    const uint32_t u = queries[qi];
    int64_t acc = 0;
    for (const auto& entry : weighted_points) {
      const uint32_t x = entry.first;
      const int64_t w = entry.second;
      const uint32_t parity = Parity32(u & x);
      acc += (parity == 0u) ? w : -w;
    }
    out[qi] = acc;
  }
}
