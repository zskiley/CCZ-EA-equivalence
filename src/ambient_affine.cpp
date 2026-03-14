#include "ambient_affine.h"

#include "graph.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

namespace {

int MostSignificantBit(uint32_t x) {
  for (int b = 31; b >= 0; --b) {
    if (((x >> b) & 1u) != 0u) return b;
  }
  return -1;
}

void ReduceWithBasis(const std::vector<uint32_t>& basis, uint32_t* value) {
  if (value == nullptr) return;
  for (uint32_t pivot : basis) {
    const int pivot_bit = MostSignificantBit(pivot);
    if (pivot_bit >= 0 && (((*value >> pivot_bit) & 1u) != 0u)) {
      *value ^= pivot;
    }
  }
}

bool AddToEchelonBasis(uint32_t value, std::vector<uint32_t>* basis) {
  if (basis == nullptr) return false;
  ReduceWithBasis(*basis, &value);
  if (value == 0u) return false;
  const int value_bit = MostSignificantBit(value);
  std::size_t insert_at = 0;
  while (insert_at < basis->size() &&
         MostSignificantBit((*basis)[insert_at]) > value_bit) {
    ++insert_at;
  }
  basis->insert(basis->begin() + static_cast<std::ptrdiff_t>(insert_at), value);
  for (std::size_t i = 0; i < basis->size(); ++i) {
    if (i == insert_at) continue;
    if ((((*basis)[i] >> value_bit) & 1u) != 0u) (*basis)[i] ^= value;
  }
  return true;
}

bool IsInSpan(const std::vector<uint32_t>& basis, uint32_t value) {
  ReduceWithBasis(basis, &value);
  return value == 0u;
}

uint32_t ApplyLinearColumns(const std::vector<uint32_t>& cols, uint32_t x) {
  uint32_t y = 0u;
  std::size_t bit = 0;
  while (x != 0u) {
    if ((x & 1u) != 0u) y ^= cols[bit];
    x >>= 1u;
    ++bit;
  }
  return y;
}

bool InvertLinearColumns(const std::vector<uint32_t>& cols, int dimension_bits,
                         std::vector<uint32_t>* inv_cols) {
  if (inv_cols == nullptr) return false;
  if (cols.size() != static_cast<std::size_t>(dimension_bits)) return false;

  std::vector<uint32_t> rows(static_cast<std::size_t>(dimension_bits), 0u);
  std::vector<uint32_t> inv_rows(static_cast<std::size_t>(dimension_bits), 0u);
  for (int row = 0; row < dimension_bits; ++row) {
    uint32_t row_bits = 0u;
    for (int col = 0; col < dimension_bits; ++col) {
      if (((cols[static_cast<std::size_t>(col)] >> row) & 1u) != 0u) {
        row_bits |= (1u << col);
      }
    }
    rows[static_cast<std::size_t>(row)] = row_bits;
    inv_rows[static_cast<std::size_t>(row)] = (1u << row);
  }

  for (int col = 0; col < dimension_bits; ++col) {
    int pivot_row = col;
    while (pivot_row < dimension_bits &&
           ((rows[static_cast<std::size_t>(pivot_row)] >> col) & 1u) == 0u) {
      ++pivot_row;
    }
    if (pivot_row == dimension_bits) return false;
    if (pivot_row != col) {
      std::swap(rows[static_cast<std::size_t>(pivot_row)],
                rows[static_cast<std::size_t>(col)]);
      std::swap(inv_rows[static_cast<std::size_t>(pivot_row)],
                inv_rows[static_cast<std::size_t>(col)]);
    }
    for (int row = 0; row < dimension_bits; ++row) {
      if (row == col) continue;
      if (((rows[static_cast<std::size_t>(row)] >> col) & 1u) == 0u) continue;
      rows[static_cast<std::size_t>(row)] ^= rows[static_cast<std::size_t>(col)];
      inv_rows[static_cast<std::size_t>(row)] ^=
          inv_rows[static_cast<std::size_t>(col)];
    }
  }

  inv_cols->assign(static_cast<std::size_t>(dimension_bits), 0u);
  for (int col = 0; col < dimension_bits; ++col) {
    uint32_t image = 0u;
    for (int row = 0; row < dimension_bits; ++row) {
      if (((inv_rows[static_cast<std::size_t>(row)] >> col) & 1u) != 0u) {
        image |= (1u << row);
      }
    }
    (*inv_cols)[static_cast<std::size_t>(col)] = image;
  }
  return true;
}

bool BuildLinearMapFromBasisPairs(const std::vector<uint32_t>& domain_basis,
                                  const std::vector<uint32_t>& image_basis,
                                  int dimension_bits,
                                  std::vector<uint32_t>* cols) {
  if (cols == nullptr) return false;
  if (domain_basis.size() != static_cast<std::size_t>(dimension_bits)) {
    return false;
  }
  if (image_basis.size() != static_cast<std::size_t>(dimension_bits)) {
    return false;
  }
  std::vector<uint32_t> inv_domain_basis;
  if (!InvertLinearColumns(domain_basis, dimension_bits, &inv_domain_basis)) {
    return false;
  }
  cols->assign(static_cast<std::size_t>(dimension_bits), 0u);
  for (int bit = 0; bit < dimension_bits; ++bit) {
    const uint32_t coeffs = ApplyLinearColumns(inv_domain_basis, (1u << bit));
    (*cols)[static_cast<std::size_t>(bit)] =
        ApplyLinearColumns(image_basis, coeffs);
  }
  return true;
}

void BuildAffineSpanBasis(const GraphData& F, uint32_t* anchor,
                          std::vector<uint32_t>* span_basis) {
  if (anchor == nullptr || span_basis == nullptr || F.points.empty()) return;
  *anchor = F.points.front();
  span_basis->clear();
  for (uint32_t point : F.points) {
    (void)AddToEchelonBasis(point ^ *anchor, span_basis);
  }
}

uint64_t SaturatingMulU64(uint64_t a, uint64_t b) {
  if (a == 0u || b == 0u) return 0u;
  const uint64_t limit = std::numeric_limits<uint64_t>::max();
  if (a > limit / b) return limit;
  return a * b;
}

std::vector<uint32_t> ExtendToBasis(const std::vector<uint32_t>& basis,
                                    int dimension_bits) {
  std::vector<uint32_t> full_basis = basis;
  std::vector<uint32_t> span = basis;
  for (int bit = 0; bit < dimension_bits; ++bit) {
    const uint32_t e = (1u << bit);
    if (IsInSpan(span, e)) continue;
    full_basis.push_back(e);
    (void)AddToEchelonBasis(e, &span);
  }
  return full_basis;
}

}  // namespace

uint64_t KernelAffineOrder(const GraphData& F) {
  if (F.d_bits <= 0 || F.d_bits >= 31 || F.points.empty()) return 1u;
  uint32_t anchor = 0u;
  std::vector<uint32_t> span_basis;
  BuildAffineSpanBasis(F, &anchor, &span_basis);
  const int r = static_cast<int>(span_basis.size());
  const int k = F.d_bits - r;
  if (k <= 0) return 1u;

  uint64_t order = 1u;
  const uint64_t shear_factor =
      (r >= 64 || k >= 64) ? std::numeric_limits<uint64_t>::max()
                           : (1ull << static_cast<unsigned>(r));
  for (int j = 0; j < k; ++j) {
    order = SaturatingMulU64(order, shear_factor);
  }

  const uint64_t q = (1ull << static_cast<unsigned>(k));
  for (int i = 0; i < k; ++i) {
    order = SaturatingMulU64(order, q - (1ull << static_cast<unsigned>(i)));
  }
  return order;
}

std::vector<AffineMapData> BuildKernelAffineGenerators(const GraphData& F) {
  std::vector<AffineMapData> generators;
  if (F.d_bits <= 0 || F.d_bits >= 31 || F.points.empty()) return generators;

  uint32_t anchor = 0u;
  std::vector<uint32_t> span_basis;
  BuildAffineSpanBasis(F, &anchor, &span_basis);
  const std::vector<uint32_t> full_basis = ExtendToBasis(span_basis, F.d_bits);
  const std::size_t r = span_basis.size();
  const std::size_t d = full_basis.size();
  if (d <= r) return generators;

  std::vector<uint32_t> complement_basis(
      full_basis.begin() + static_cast<std::ptrdiff_t>(r), full_basis.end());

  for (std::size_t j = 0; j < complement_basis.size(); ++j) {
    for (std::size_t i = 0; i < span_basis.size(); ++i) {
      std::vector<uint32_t> image_basis = full_basis;
      image_basis[r + j] ^= span_basis[i];
      AffineMapData map;
      map.dimension_bits = F.d_bits;
      if (!BuildLinearMapFromBasisPairs(full_basis, image_basis, F.d_bits,
                                        &map.linear_cols)) {
        continue;
      }
      map.translation = anchor ^ ApplyLinearColumns(map.linear_cols, anchor);
      if (!map.IsIdentity()) generators.push_back(std::move(map));
    }
  }

  for (std::size_t j = 0; j < complement_basis.size(); ++j) {
    for (std::size_t i = 0; i < complement_basis.size(); ++i) {
      if (i == j) continue;
      std::vector<uint32_t> image_basis = full_basis;
      image_basis[r + j] ^= complement_basis[i];
      AffineMapData map;
      map.dimension_bits = F.d_bits;
      if (!BuildLinearMapFromBasisPairs(full_basis, image_basis, F.d_bits,
                                        &map.linear_cols)) {
        continue;
      }
      map.translation = anchor ^ ApplyLinearColumns(map.linear_cols, anchor);
      if (!map.IsIdentity()) generators.push_back(std::move(map));
    }
  }

  return generators;
}

bool LiftGraphPermutationToAffineMap(const GraphData& F,
                                     const groups::Permutation& perm,
                                     AffineMapData* out) {
  if (out == nullptr) return false;
  if (perm.Degree() != F.points.size()) return false;

  PartialAffineMap A(F.d_bits);
  for (std::size_t i = 0; i < F.points.size(); ++i) {
    if (!A.Update(F.points[i], F.points[perm.Apply(static_cast<uint32_t>(i))])) {
      return false;
    }
  }
  return A.ExtractRepresentativeAffineMap(out);
}

std::vector<AffineMapData> BuildAmbientAutoGenerators(
    const GraphData& F, const std::vector<groups::Permutation>& graph_generators) {
  std::vector<AffineMapData> generators = BuildKernelAffineGenerators(F);
  for (const groups::Permutation& perm : graph_generators) {
    if (perm.IsIdentity()) continue;
    AffineMapData map;
    if (!LiftGraphPermutationToAffineMap(F, perm, &map)) continue;
    if (!map.IsIdentity()) generators.push_back(std::move(map));
  }
  return generators;
}
