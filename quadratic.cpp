#include "quadratic.h"

#include "graph.h"
#include "groups/graph_point_perm.h"

#include <cstddef>
#include <cstdint>
#include <limits>
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

bool DecodeGraphAsFunction(const GraphData& F, int* n_bits,
                           std::vector<uint32_t>* f_values) {
  if (n_bits == nullptr || f_values == nullptr) return false;
  if (F.d_bits <= 0 || (F.d_bits % 2) != 0 || F.d_bits >= 31) return false;

  const int n = F.d_bits / 2;
  const uint32_t q = (1u << n);
  if (F.points.size() != static_cast<std::size_t>(q)) return false;

  const uint32_t mask = q - 1u;
  const uint32_t unassigned = std::numeric_limits<uint32_t>::max();
  f_values->assign(q, unassigned);

  for (uint32_t p : F.points) {
    const uint32_t x = p & mask;
    const uint32_t y = (p >> n);
    if (x >= q || y >= q) return false;
    uint32_t& slot = (*f_values)[x];
    if (slot != unassigned && slot != y) return false;
    slot = y;
  }

  for (uint32_t y : *f_values) {
    if (y == unassigned) return false;
  }

  *n_bits = n;
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

}  // namespace

bool TryQuadraticAnchorPoint(const GraphData& F, uint32_t* anchor_point) {
  if (anchor_point == nullptr) return false;
  int n_bits = 0;
  std::vector<uint32_t> f_values;
  if (!DecodeGraphAsFunction(F, &n_bits, &f_values)) return false;
  if (!IsQuadraticFunction(f_values, n_bits)) return false;

  const uint32_t f0 = f_values[0u];
  *anchor_point = (f0 << n_bits);
  return true;
}

bool BuildQuadraticTranslationGenerators(
    const GraphData& F, std::vector<groups::Permutation>* generators) {
  if (generators == nullptr) return false;
  generators->clear();

  int n_bits = 0;
  std::vector<uint32_t> f_values;
  if (!DecodeGraphAsFunction(F, &n_bits, &f_values)) return false;
  if (!IsQuadraticFunction(f_values, n_bits)) return false;

  groups::GraphPointIndex point_index(F);

  const std::size_t degree = F.points.size();
  if (degree > static_cast<std::size_t>(
                   std::numeric_limits<uint16_t>::max()) +
                   1u) {
    return false;
  }

  const uint32_t q = (1u << n_bits);
  const uint32_t x_mask = q - 1u;
  generators->reserve(static_cast<std::size_t>(n_bits));

  for (int bit = 0; bit < n_bits; ++bit) {
    const uint32_t a = (1u << bit);
    std::vector<uint16_t> images(degree, 0u);

    for (std::size_t i = 0; i < degree; ++i) {
      const uint32_t p = F.points[i];
      const uint32_t x = p & x_mask;
      const uint32_t xa = x ^ a;
      const uint32_t image_point = xa | (f_values[xa] << n_bits);
      const auto image_idx = point_index.IndexOf(image_point);
      if (!image_idx.has_value() || *image_idx >= degree) return false;
      images[i] = static_cast<uint16_t>(*image_idx);
    }
    generators->emplace_back(std::move(images));
  }

  return true;
}
