#include "quadratic.h"

#include "graph.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace {

int HammingWeight32(uint32_t x) {
  int w = 0;
  while (x != 0u) {
    x &= (x - 1u);
    ++w;
  }
  return w;
}

struct QuadraticGraphContext {
  int n_bits = 0;
  uint32_t x_mask = 0u;
  std::vector<uint32_t> f_values;
  std::vector<uint32_t> point_index_by_x;
};

bool DecodeGraphAsFunction(const GraphData& F,
                           QuadraticGraphContext* context) {
  if (context == nullptr) return false;
  if (F.d_bits <= 0 || (F.d_bits % 2) != 0 || F.d_bits >= 31) return false;

  const int n = F.d_bits / 2;
  const uint32_t q = (1u << n);
  if (F.points.size() != static_cast<std::size_t>(q)) return false;

  const uint32_t unassigned = std::numeric_limits<uint32_t>::max();
  QuadraticGraphContext next;
  next.n_bits = n;
  next.x_mask = q - 1u;
  next.f_values.assign(q, unassigned);
  next.point_index_by_x.assign(q, unassigned);

  for (std::size_t i = 0; i < F.points.size(); ++i) {
    const uint32_t p = F.points[i];
    const uint32_t x = p & next.x_mask;
    const uint32_t y = (p >> n);
    if (y >= q) return false;
    uint32_t& slot = next.f_values[x];
    if (slot != unassigned && slot != y) return false;
    slot = y;
    next.point_index_by_x[x] = static_cast<uint32_t>(i);
  }

  for (uint32_t y : next.f_values) {
    if (y == unassigned) return false;
  }

  *context = std::move(next);
  return true;
}

bool IsQuadraticFunction(const std::vector<uint32_t>& f_values, int n_bits) {
  if (n_bits <= 0 || n_bits >= 31) return false;
  const uint32_t q = (1u << n_bits);
  if (f_values.size() != static_cast<std::size_t>(q)) return false;

  std::vector<uint8_t> anf(q, 0u);
  for (int out_bit = 0; out_bit < n_bits; ++out_bit) {
    for (uint32_t x = 0; x < q; ++x) {
      anf[x] = static_cast<uint8_t>((f_values[x] >> out_bit) & 1u);
    }

    for (int b = 0; b < n_bits; ++b) {
      const uint32_t bit = (1u << b);
      for (uint32_t mask = 0; mask < q; ++mask) {
        if ((mask & bit) != 0u) anf[mask] ^= anf[mask ^ bit];
      }
    }

    for (uint32_t mask = 0; mask < q; ++mask) {
      if (anf[mask] == 0u) continue;
      if (HammingWeight32(mask) > 2) return false;
    }
  }
  return true;
}

bool BuildQuadraticGraphContext(const GraphData& F,
                                QuadraticGraphContext* context) {
  if (context == nullptr) return false;

  QuadraticGraphContext next;
  if (!DecodeGraphAsFunction(F, &next)) return false;
  if (!IsQuadraticFunction(next.f_values, next.n_bits)) return false;

  *context = std::move(next);
  return true;
}

uint32_t AnchorPoint(const QuadraticGraphContext& context) {
  return context.f_values[0u] << context.n_bits;
}

bool BuildTranslationGenerators(
    const GraphData& F, const QuadraticGraphContext& context,
    std::vector<groups::Permutation>* generators) {
  if (generators == nullptr) return false;
  generators->clear();

  const std::size_t degree = F.points.size();
  if (degree > static_cast<std::size_t>(
                   std::numeric_limits<uint16_t>::max()) +
                   1u) {
    return false;
  }

  generators->reserve(static_cast<std::size_t>(context.n_bits));

  for (int bit = 0; bit < context.n_bits; ++bit) {
    const uint32_t a = (1u << bit);
    std::vector<uint16_t> images(degree, 0u);

    for (std::size_t i = 0; i < degree; ++i) {
      const uint32_t p = F.points[i];
      const uint32_t x = p & context.x_mask;
      const uint32_t xa = x ^ a;
      const uint32_t image_idx = context.point_index_by_x[xa];
      if (image_idx >= degree) return false;
      images[i] = static_cast<uint16_t>(image_idx);
    }
    generators->emplace_back(std::move(images));
  }

  return true;
}

bool BuildTranslationAffineGenerators(
    const GraphData& F, const QuadraticGraphContext& context,
    std::vector<AffineMapData>* generators) {
  if (generators == nullptr) return false;
  generators->clear();

  generators->reserve(static_cast<std::size_t>(context.n_bits));

  for (int bit = 0; bit < context.n_bits; ++bit) {
    const uint32_t a = (1u << bit);
    const uint32_t d0 = context.f_values[a] ^ context.f_values[0u];

    AffineMapData map;
    map.dimension_bits = F.d_bits;
    map.translation = a | (d0 << context.n_bits);
    map.linear_cols.assign(static_cast<std::size_t>(F.d_bits), 0u);

    for (int x_bit = 0; x_bit < context.n_bits; ++x_bit) {
      const uint32_t e = (1u << x_bit);
      const uint32_t de = context.f_values[e ^ a] ^ context.f_values[e];
      const uint32_t linear_output = de ^ d0;
      map.linear_cols[static_cast<std::size_t>(x_bit)] =
          e | (linear_output << context.n_bits);
    }
    for (int y_bit = 0; y_bit < context.n_bits; ++y_bit) {
      map.linear_cols[static_cast<std::size_t>(context.n_bits + y_bit)] =
          (1u << (context.n_bits + y_bit));
    }

    if (!map.IsIdentity()) generators->push_back(std::move(map));
  }

  return true;
}

}  // namespace

bool TryQuadraticAnchorPoint(const GraphData& F, uint32_t* anchor_point) {
  if (anchor_point == nullptr) return false;
  QuadraticGraphContext context;
  if (!BuildQuadraticGraphContext(F, &context)) return false;

  *anchor_point = AnchorPoint(context);
  return true;
}

bool BuildQuadraticTranslationGenerators(
    const GraphData& F, std::vector<groups::Permutation>* generators) {
  if (generators == nullptr) return false;
  generators->clear();

  QuadraticGraphContext context;
  if (!BuildQuadraticGraphContext(F, &context)) return false;

  return BuildTranslationGenerators(F, context, generators);
}

bool BuildQuadraticTranslationAffineGenerators(
    const GraphData& F, std::vector<AffineMapData>* generators) {
  if (generators == nullptr) return false;
  generators->clear();

  QuadraticGraphContext context;
  if (!BuildQuadraticGraphContext(F, &context)) return false;

  return BuildTranslationAffineGenerators(F, context, generators);
}

bool BuildQuadraticTranslationData(
    const GraphData& F, uint32_t* anchor_point,
    std::vector<groups::Permutation>* translation_generators,
    std::vector<AffineMapData>* translation_affine_generators) {
  if (anchor_point == nullptr || translation_generators == nullptr ||
      translation_affine_generators == nullptr) {
    return false;
  }
  translation_generators->clear();
  translation_affine_generators->clear();

  QuadraticGraphContext context;
  if (!BuildQuadraticGraphContext(F, &context)) return false;

  std::vector<groups::Permutation> new_translation_generators;
  std::vector<AffineMapData> new_translation_affine_generators;
  if (!BuildTranslationGenerators(F, context, &new_translation_generators)) {
    return false;
  }
  if (!BuildTranslationAffineGenerators(
          F, context, &new_translation_affine_generators)) {
    return false;
  }

  *anchor_point = AnchorPoint(context);
  *translation_generators = std::move(new_translation_generators);
  *translation_affine_generators =
      std::move(new_translation_affine_generators);
  return true;
}
