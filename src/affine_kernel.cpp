#include "affine_kernel.h"

#include "gf2_linear.h"
#include "graph.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace {

struct KernelContext {
  uint32_t anchor = 0u;
  int ambient_dim = 0;
  std::vector<uint32_t> graph_span_basis;
  bool valid = false;

  int GraphSpanDim() const {
    return static_cast<int>(graph_span_basis.size());
  }

  int ComplementDim() const { return ambient_dim - GraphSpanDim(); }
};

struct EaKernelDecomposition {
  std::vector<uint32_t> span_output;
  std::vector<uint32_t> span_only;
  std::vector<uint32_t> output_only;
  std::vector<uint32_t> outside;
  std::vector<uint32_t> basis;
  bool valid = false;
};

KernelContext BuildKernelContext(const GraphData& F) {
  KernelContext context;
  context.ambient_dim = F.d_bits;
  if (F.d_bits <= 0 || F.d_bits >= 31 || F.points.empty()) return context;

  context.anchor = F.points.front();
  for (uint32_t point : F.points) {
    (void)gf2_linear::AddToEchelonBasis(point ^ context.anchor,
                                        &context.graph_span_basis);
  }
  context.valid = true;
  return context;
}

uint64_t SaturatingMulU64(uint64_t a, uint64_t b) {
  if (a == 0u || b == 0u) return 0u;
  const uint64_t limit = std::numeric_limits<uint64_t>::max();
  if (a > limit / b) return limit;
  return a * b;
}

uint64_t PowerOfTwoOrder(int exponent) {
  if (exponent <= 0) return 1u;
  if (exponent >= 64) return std::numeric_limits<uint64_t>::max();
  return (1ull << static_cast<unsigned>(exponent));
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

std::vector<uint32_t> StandardBasis(int dimension_bits) {
  std::vector<uint32_t> basis;
  basis.reserve(static_cast<std::size_t>(dimension_bits));
  for (int bit = 0; bit < dimension_bits; ++bit) {
    basis.push_back(1u << bit);
  }
  return basis;
}

void AppendIndependentVectors(const std::vector<uint32_t>& candidates,
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

void AppendAll(const std::vector<uint32_t>& source,
               std::vector<uint32_t>* target) {
  target->insert(target->end(), source.begin(), source.end());
}

EaKernelDecomposition BuildEaKernelDecomposition(
    const GraphData& F, const KernelContext& context) {
  EaKernelDecomposition parts;
  if (!context.valid) return parts;

  const std::vector<uint32_t> output_basis = OutputSubspaceBasis(F);
  if (output_basis.empty()) return parts;

  parts.span_output =
      gf2_linear::IntersectionBasis(context.graph_span_basis, output_basis);
  std::vector<uint32_t> span = parts.span_output;

  AppendIndependentVectors(context.graph_span_basis, &span, &parts.span_only);
  AppendIndependentVectors(output_basis, &span, &parts.output_only);
  AppendIndependentVectors(StandardBasis(context.ambient_dim), &span,
                           &parts.outside);

  parts.basis.reserve(static_cast<std::size_t>(context.ambient_dim));
  AppendAll(parts.span_output, &parts.basis);
  AppendAll(parts.span_only, &parts.basis);
  AppendAll(parts.output_only, &parts.basis);
  AppendAll(parts.outside, &parts.basis);
  parts.valid =
      parts.basis.size() == static_cast<std::size_t>(context.ambient_dim);
  return parts;
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

void AddBasisTransvection(const GraphData& F, const KernelContext& context,
                          const std::vector<uint32_t>& basis,
                          std::size_t target_col, std::size_t add_col,
                          bool require_ea,
                          std::vector<AffineMapData>* generators) {
  if (target_col >= basis.size() || add_col >= basis.size()) return;
  if (target_col == add_col) return;

  std::vector<uint32_t> image_basis = basis;
  image_basis[target_col] ^= basis[add_col];

  AffineMapData map;
  map.dimension_bits = F.d_bits;
  if (!gf2_linear::BuildLinearMapFromBasisPairs(
          basis, image_basis, F.d_bits, &map.linear_cols)) {
    return;
  }
  map.translation =
      context.anchor ^ gf2_linear::ApplyLinearColumns(map.linear_cols,
                                                      context.anchor);
  AddKernelGenerator(generators, std::move(map), require_ea, F.n_bits);
}

void AddTransvections(const GraphData& F, const KernelContext& context,
                      const std::vector<uint32_t>& basis,
                      std::size_t target_begin, std::size_t target_end,
                      std::size_t add_begin, std::size_t add_end,
                      bool require_ea,
                      std::vector<AffineMapData>* generators) {
  for (std::size_t target = target_begin; target < target_end; ++target) {
    for (std::size_t add = add_begin; add < add_end; ++add) {
      AddBasisTransvection(F, context, basis, target, add, require_ea,
                           generators);
    }
  }
}

uint64_t CczKernelOrder(const KernelContext& context) {
  const int span_dim = context.GraphSpanDim();
  const int complement_dim = context.ComplementDim();
  if (complement_dim <= 0) return 1u;

  uint64_t order = PowerOfTwoOrder(span_dim * complement_dim);
  order = SaturatingMulU64(order, GLOrder(complement_dim));
  return order;
}

uint64_t EaKernelOrder(const EaKernelDecomposition& parts) {
  if (!parts.valid) return 1u;

  const int span_output_dim = static_cast<int>(parts.span_output.size());
  const int span_only_dim = static_cast<int>(parts.span_only.size());
  const int output_only_dim = static_cast<int>(parts.output_only.size());
  const int outside_dim = static_cast<int>(parts.outside.size());

  const int unipotent_exponent =
      (span_output_dim * output_only_dim) +
      ((span_output_dim + span_only_dim + output_only_dim) * outside_dim);

  uint64_t order = PowerOfTwoOrder(unipotent_exponent);
  order = SaturatingMulU64(order, GLOrder(output_only_dim));
  order = SaturatingMulU64(order, GLOrder(outside_dim));
  return order;
}

std::vector<AffineMapData> BuildCczKernelGenerators(
    const GraphData& F, const KernelContext& context) {
  std::vector<AffineMapData> generators;
  const std::vector<uint32_t> basis =
      gf2_linear::ExtendToBasis(context.graph_span_basis, context.ambient_dim);

  const std::size_t span_dim = context.graph_span_basis.size();
  AddTransvections(F, context, basis, span_dim, basis.size(), 0, basis.size(),
                   /*require_ea=*/false, &generators);
  return generators;
}

std::vector<AffineMapData> BuildEaKernelGenerators(
    const GraphData& F, const KernelContext& context) {
  std::vector<AffineMapData> generators;
  const EaKernelDecomposition parts = BuildEaKernelDecomposition(F, context);
  if (!parts.valid) return generators;

  const std::size_t span_output_dim = parts.span_output.size();
  const std::size_t span_only_dim = parts.span_only.size();
  const std::size_t output_only_dim = parts.output_only.size();
  const std::size_t output_only_begin = span_output_dim + span_only_dim;
  const std::size_t outside_begin = output_only_begin + output_only_dim;

  AddTransvections(F, context, parts.basis, output_only_begin, outside_begin, 0,
                   span_output_dim, /*require_ea=*/true, &generators);
  AddTransvections(F, context, parts.basis, output_only_begin, outside_begin,
                   output_only_begin, outside_begin, /*require_ea=*/true,
                   &generators);
  AddTransvections(F, context, parts.basis, outside_begin, parts.basis.size(),
                   0, parts.basis.size(), /*require_ea=*/true, &generators);
  return generators;
}

}  // namespace

uint64_t AffineKernelOrder(const GraphData& F, bool require_ea) {
  const KernelContext context = BuildKernelContext(F);
  if (!context.valid || context.ComplementDim() <= 0) return 1u;

  if (!require_ea) return CczKernelOrder(context);
  return EaKernelOrder(BuildEaKernelDecomposition(F, context));
}

std::vector<AffineMapData> BuildAffineKernelGenerators(const GraphData& F,
                                                       bool require_ea) {
  const KernelContext context = BuildKernelContext(F);
  if (!context.valid || context.ComplementDim() <= 0) return {};

  if (!require_ea) return BuildCczKernelGenerators(F, context);
  return BuildEaKernelGenerators(F, context);
}
