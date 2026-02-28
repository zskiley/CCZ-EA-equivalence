#ifndef WEIGHTED_FFT_H
#define WEIGHTED_FFT_H

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

using WeightedPoint = std::pair<uint32_t, int64_t>;

// Returns true iff n is a power of two and n > 0.
bool IsPowerOfTwo(std::size_t n);

// In-place unnormalized Walsh-Hadamard transform over F_2^m.
// Requires values != nullptr and values->size() = 2^m for some m >= 0.
bool FWHTInPlace(std::vector<int64_t>* values);

// Builds a dense weighted function f : F_2^m -> Z from sparse points.
// weighted_points stores (x, weight). Repeated x values are accumulated.
// Returns false on invalid dimension or out-of-range x.
bool BuildWeightedVector(int dimension_bits,
                         const std::vector<WeightedPoint>& weighted_points,
                         std::vector<int64_t>* values);

// Computes the unnormalized Walsh-Hadamard spectrum of a sparse weighted input.
// Equivalent to:
// 1) BuildWeightedVector(...)
// 2) FWHTInPlace(...)
void WeightedFWHT(int dimension_bits,
                  const std::vector<WeightedPoint>& weighted_points,
                  std::vector<int64_t>* spectrum);

// Computes selected Walsh-Hadamard spectrum values without computing the full
// 2^m-sized transform.
// For each query u, returns:
//   W(u) = sum_{(x,w)} w * (-1)^{<u,x>}
// where <u,x> is the F_2 dot-product (parity of popcount(u & x)).
// Runs in O(|queries| * |weighted_points|) time.
void WeightedWalshAt(int dimension_bits,
                     const std::vector<WeightedPoint>& weighted_points,
                     const std::vector<uint32_t>& queries,
                     std::vector<int64_t>* values_at_queries);

#endif
