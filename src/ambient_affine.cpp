#include "ambient_affine.h"

#include "gf2_linear.h"
#include "graph.h"

#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace {

struct EchelonMaskRow {
  uint32_t value = 0u;
  uint64_t coeff_mask = 0u;
};

struct EaKernelDecomposition {
  std::vector<uint32_t> a;
  std::vector<uint32_t> b;
  std::vector<uint32_t> c;
  std::vector<uint32_t> d;
  std::vector<uint32_t> full_basis;
  bool valid = false;
};

void BuildAffineSpanBasis(const GraphData& F, uint32_t* anchor,
                          std::vector<uint32_t>* span_basis) {
  if (anchor == nullptr || span_basis == nullptr || F.points.empty()) return;
  *anchor = F.points.front();
  span_basis->clear();
  for (uint32_t point : F.points) {
    (void)gf2_linear::AddToEchelonBasis(point ^ *anchor, span_basis);
  }
}

uint64_t SaturatingMulU64(uint64_t a, uint64_t b) {
  if (a == 0u || b == 0u) return 0u;
  const uint64_t limit = std::numeric_limits<uint64_t>::max();
  if (a > limit / b) return limit;
  return a * b;
}

bool HasEaBlockForm(const AffineMapData& map, int n_bits) {
  if (n_bits <= 0 || n_bits >= map.dimension_bits) return false;
  const uint32_t x_mask = (1u << n_bits) - 1u;
  for (int bit = n_bits; bit < map.dimension_bits; ++bit) {
    if ((map.linear_cols[static_cast<std::size_t>(bit)] & x_mask) != 0u) {
      return false;
    }
  }
  return true;
}

uint64_t GLOrder(int dimension) {
  if (dimension <= 0) return 1u;
  if (dimension >= 64) return std::numeric_limits<uint64_t>::max();
  const uint64_t q = (1ull << static_cast<unsigned>(dimension));
  uint64_t order = 1u;
  for (int i = 0; i < dimension; ++i) {
    order = SaturatingMulU64(order, q - (1ull << static_cast<unsigned>(i)));
  }
  return order;
}

uint64_t PowerOfTwoOrder(int exponent) {
  if (exponent <= 0) return 1u;
  if (exponent >= 64) return std::numeric_limits<uint64_t>::max();
  return (1ull << static_cast<unsigned>(exponent));
}

bool AddMaskRow(uint32_t value, uint64_t coeff_mask,
                std::vector<EchelonMaskRow>* basis) {
  if (basis == nullptr) return false;
  for (const EchelonMaskRow& row : *basis) {
    const int pivot_bit = gf2_linear::MostSignificantBit(row.value);
    if (pivot_bit >= 0 && (((value >> pivot_bit) & 1u) != 0u)) {
      value ^= row.value;
      coeff_mask ^= row.coeff_mask;
    }
  }
  if (value == 0u) return false;

  const int value_bit = gf2_linear::MostSignificantBit(value);
  std::size_t insert_at = 0;
  while (insert_at < basis->size() &&
         gf2_linear::MostSignificantBit((*basis)[insert_at].value) >
             value_bit) {
    ++insert_at;
  }
  basis->insert(basis->begin() + static_cast<std::ptrdiff_t>(insert_at),
                EchelonMaskRow{value, coeff_mask});
  for (std::size_t i = 0; i < basis->size(); ++i) {
    if (i == insert_at) continue;
    if ((((*basis)[i].value >> value_bit) & 1u) == 0u) continue;
    (*basis)[i].value ^= value;
    (*basis)[i].coeff_mask ^= coeff_mask;
  }
  return true;
}

std::vector<uint32_t> IntersectionBasis(const std::vector<uint32_t>& left,
                                        const std::vector<uint32_t>& right) {
  std::vector<EchelonMaskRow> combined_basis;
  std::vector<uint32_t> intersection;

  const auto add_combined = [&](uint32_t value, uint64_t coeff_mask) {
    uint32_t reduced_value = value;
    uint64_t reduced_coeff_mask = coeff_mask;
    for (const EchelonMaskRow& row : combined_basis) {
      const int pivot_bit = gf2_linear::MostSignificantBit(row.value);
      if (pivot_bit >= 0 && (((reduced_value >> pivot_bit) & 1u) != 0u)) {
        reduced_value ^= row.value;
        reduced_coeff_mask ^= row.coeff_mask;
      }
    }
    if (reduced_value != 0u) {
      (void)AddMaskRow(value, coeff_mask, &combined_basis);
      return;
    }

    uint32_t intersection_vector = 0u;
    for (std::size_t i = 0; i < left.size(); ++i) {
      if (((reduced_coeff_mask >> i) & 1ull) != 0u) {
        intersection_vector ^= left[i];
      }
    }
    if (intersection_vector != 0u) {
      (void)gf2_linear::AddToEchelonBasis(intersection_vector, &intersection);
    }
  };

  for (std::size_t i = 0; i < left.size(); ++i) {
    add_combined(left[i], (1ull << static_cast<unsigned>(i)));
  }
  for (std::size_t i = 0; i < right.size(); ++i) {
    const std::size_t bit = left.size() + i;
    if (bit >= 64) break;
    add_combined(right[i], (1ull << static_cast<unsigned>(bit)));
  }

  return intersection;
}

void ExtendComponent(const std::vector<uint32_t>& candidates,
                     std::vector<uint32_t>* span,
                     std::vector<uint32_t>* component) {
  if (span == nullptr || component == nullptr) return;
  for (uint32_t candidate : candidates) {
    std::vector<uint32_t> next_span = *span;
    if (!gf2_linear::AddToEchelonBasis(candidate, &next_span)) continue;
    *span = std::move(next_span);
    component->push_back(candidate);
  }
}

std::vector<uint32_t> OutputSubspaceBasis(const GraphData& F) {
  std::vector<uint32_t> basis;
  if (F.n_bits <= 0 || F.m_bits <= 0 || F.n_bits + F.m_bits != F.d_bits) {
    return basis;
  }
  basis.reserve(static_cast<std::size_t>(F.m_bits));
  for (int bit = F.n_bits; bit < F.d_bits; ++bit) {
    basis.push_back(1u << bit);
  }
  return basis;
}

EaKernelDecomposition BuildEaKernelDecomposition(
    const GraphData& F, const std::vector<uint32_t>& span_basis) {
  EaKernelDecomposition decomposition;
  const std::vector<uint32_t> y_basis = OutputSubspaceBasis(F);
  if (y_basis.empty()) return decomposition;

  decomposition.a = IntersectionBasis(span_basis, y_basis);
  std::vector<uint32_t> span = decomposition.a;

  ExtendComponent(span_basis, &span, &decomposition.b);
  ExtendComponent(y_basis, &span, &decomposition.c);

  std::vector<uint32_t> standard_basis;
  standard_basis.reserve(static_cast<std::size_t>(F.d_bits));
  for (int bit = 0; bit < F.d_bits; ++bit) {
    standard_basis.push_back(1u << bit);
  }
  ExtendComponent(standard_basis, &span, &decomposition.d);

  decomposition.full_basis.reserve(static_cast<std::size_t>(F.d_bits));
  decomposition.full_basis.insert(decomposition.full_basis.end(),
                                  decomposition.a.begin(),
                                  decomposition.a.end());
  decomposition.full_basis.insert(decomposition.full_basis.end(),
                                  decomposition.b.begin(),
                                  decomposition.b.end());
  decomposition.full_basis.insert(decomposition.full_basis.end(),
                                  decomposition.c.begin(),
                                  decomposition.c.end());
  decomposition.full_basis.insert(decomposition.full_basis.end(),
                                  decomposition.d.begin(),
                                  decomposition.d.end());
  decomposition.valid =
      decomposition.full_basis.size() == static_cast<std::size_t>(F.d_bits);
  return decomposition;
}

void AddKernelGenerator(std::vector<AffineMapData>* generators,
                        AffineMapData map, bool require_ea, int n_bits) {
  if (generators == nullptr || map.IsIdentity()) return;
  if (require_ea && !HasEaBlockForm(map, n_bits)) return;
  for (const AffineMapData& existing : *generators) {
    if (existing.translation == map.translation &&
        existing.linear_cols == map.linear_cols) {
      return;
    }
  }
  generators->push_back(std::move(map));
}

void AddBasisTransvection(const GraphData& F, uint32_t anchor,
                          const std::vector<uint32_t>& full_basis,
                          std::size_t target_col, std::size_t add_row,
                          bool require_ea,
                          std::vector<AffineMapData>* generators) {
  if (target_col >= full_basis.size() || add_row >= full_basis.size()) return;
  if (target_col == add_row) return;

  std::vector<uint32_t> image_basis = full_basis;
  image_basis[target_col] ^= full_basis[add_row];

  AffineMapData map;
  map.dimension_bits = F.d_bits;
  if (!gf2_linear::BuildLinearMapFromBasisPairs(
          full_basis, image_basis, F.d_bits, &map.linear_cols)) {
    return;
  }
  map.translation =
      anchor ^ gf2_linear::ApplyLinearColumns(map.linear_cols, anchor);
  AddKernelGenerator(generators, std::move(map), require_ea, F.n_bits);
}

std::vector<AffineMapData> BuildCczKernelAffineGenerators(
    const GraphData& F, uint32_t anchor,
    const std::vector<uint32_t>& span_basis) {
  std::vector<AffineMapData> generators;
  const std::vector<uint32_t> full_basis =
      gf2_linear::ExtendToBasis(span_basis, F.d_bits);
  const std::size_t r = span_basis.size();
  const std::size_t d = full_basis.size();
  if (d <= r) return generators;

  for (std::size_t j = r; j < d; ++j) {
    for (std::size_t i = 0; i < r; ++i) {
      AddBasisTransvection(F, anchor, full_basis, j, i,
                           /*require_ea=*/false, &generators);
    }
  }

  for (std::size_t j = r; j < d; ++j) {
    for (std::size_t i = r; i < d; ++i) {
      AddBasisTransvection(F, anchor, full_basis, j, i,
                           /*require_ea=*/false, &generators);
    }
  }

  return generators;
}

std::vector<AffineMapData> BuildEaKernelAffineGenerators(
    const GraphData& F, uint32_t anchor,
    const std::vector<uint32_t>& span_basis) {
  std::vector<AffineMapData> generators;
  const EaKernelDecomposition decomposition =
      BuildEaKernelDecomposition(F, span_basis);
  if (!decomposition.valid) return generators;

  const std::vector<uint32_t>& full_basis = decomposition.full_basis;
  const std::size_t a = decomposition.a.size();
  const std::size_t b = decomposition.b.size();
  const std::size_t c = decomposition.c.size();
  const std::size_t d0 = decomposition.d.size();
  const std::size_t c_begin = a + b;
  const std::size_t d_begin = c_begin + c;

  for (std::size_t j = 0; j < c; ++j) {
    const std::size_t col = c_begin + j;
    for (std::size_t i = 0; i < a; ++i) {
      AddBasisTransvection(F, anchor, full_basis, col, i,
                           /*require_ea=*/true, &generators);
    }
    for (std::size_t i = 0; i < c; ++i) {
      AddBasisTransvection(F, anchor, full_basis, col, c_begin + i,
                           /*require_ea=*/true, &generators);
    }
  }

  for (std::size_t j = 0; j < d0; ++j) {
    const std::size_t col = d_begin + j;
    for (std::size_t i = 0; i < d_begin; ++i) {
      AddBasisTransvection(F, anchor, full_basis, col, i,
                           /*require_ea=*/true, &generators);
    }
    for (std::size_t i = 0; i < d0; ++i) {
      AddBasisTransvection(F, anchor, full_basis, col, d_begin + i,
                           /*require_ea=*/true, &generators);
    }
  }

  return generators;
}

}  // namespace

uint64_t KernelAffineOrder(const GraphData& F, bool require_ea) {
  if (F.d_bits <= 0 || F.d_bits >= 31 || F.points.empty()) return 1u;
  uint32_t anchor = 0u;
  std::vector<uint32_t> span_basis;
  BuildAffineSpanBasis(F, &anchor, &span_basis);
  const int r = static_cast<int>(span_basis.size());
  const int k = F.d_bits - r;
  if (k <= 0) return 1u;

  if (require_ea) {
    const EaKernelDecomposition decomposition =
        BuildEaKernelDecomposition(F, span_basis);
    if (!decomposition.valid) return 1u;

    const int a = static_cast<int>(decomposition.a.size());
    const int b = static_cast<int>(decomposition.b.size());
    const int c = static_cast<int>(decomposition.c.size());
    const int d0 = static_cast<int>(decomposition.d.size());

    uint64_t order = PowerOfTwoOrder((a * c) + ((a + b + c) * d0));
    order = SaturatingMulU64(order, GLOrder(c));
    order = SaturatingMulU64(order, GLOrder(d0));
    return order;
  }

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

std::vector<AffineMapData> BuildKernelAffineGenerators(const GraphData& F,
                                                       bool require_ea) {
  std::vector<AffineMapData> generators;
  if (F.d_bits <= 0 || F.d_bits >= 31 || F.points.empty()) return generators;

  uint32_t anchor = 0u;
  std::vector<uint32_t> span_basis;
  BuildAffineSpanBasis(F, &anchor, &span_basis);
  if (require_ea) {
    return BuildEaKernelAffineGenerators(F, anchor, span_basis);
  }
  return BuildCczKernelAffineGenerators(F, anchor, span_basis);
}
