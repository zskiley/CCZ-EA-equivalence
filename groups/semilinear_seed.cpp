#include "semilinear_seed.h"

#include "../field_basics.h"
#include "../graph.h"
#include "graph_point_perm.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace groups {
namespace {

void BuildYOfXTable(const GraphData& graph, int n,
                    std::vector<uint32_t>& out_y_of_x) {
  const std::size_t size = static_cast<std::size_t>(1u) << n;
  const uint32_t mask = (static_cast<uint32_t>(1u) << n) - 1u;

  out_y_of_x.assign(size, 0u);
  for (uint32_t pt : graph.points) {
    const uint32_t x = pt & mask;
    const uint32_t y = pt >> n;
    out_y_of_x[x] = y;
  }
}

// For fixed nonzero a, find b (if it exists) such that
//   f(a*x) = b*f(x) for all x,
// where y_of_x[x] = f(x). We derive b from one nonzero sample f(x0) when
// available, then verify the identity for all x.
bool ComputeScalingBAndVerify(const GF2n& field,
                              const std::vector<uint32_t>& y_of_x, uint32_t a,
                              uint32_t* out_b) {
  const std::size_t size = y_of_x.size();

  uint32_t x0 = 0u;
  while (static_cast<std::size_t>(x0) < size && y_of_x[x0] == 0u) ++x0;

  uint32_t b = 1u;
  if (static_cast<std::size_t>(x0) < size) {
    const uint32_t y0 = y_of_x[x0];
    const uint32_t ax0 = field.Mul(a, x0);
    const uint32_t inv_y0 =
        field.Pow(y0, static_cast<int>(static_cast<uint32_t>(size) - 2u));
    b = field.Mul(y_of_x[ax0], inv_y0);
  }

  for (uint32_t x = 0u; static_cast<std::size_t>(x) < size; ++x) {
    const uint32_t ax = field.Mul(a, x);
    const uint32_t lhs = y_of_x[ax];
    const uint32_t rhs = field.Mul(b, y_of_x[x]);
    if (lhs != rhs) return false;
  }

  if (out_b != nullptr) *out_b = b;
  return true;
}

std::optional<Permutation> BuildDiagonalPermutation(const GraphData& graph,
                                                    const GF2n& field, int n,
                                                    uint32_t a, uint32_t b) {
  if (graph.is_graph.empty()) return std::nullopt;
  GraphPointIndex index(graph);

  const std::size_t degree = index.Size();
  if (degree != graph.points.size()) return std::nullopt;

  const uint32_t mask = (static_cast<uint32_t>(1u) << n) - 1u;
  std::vector<uint16_t> images(degree, 0u);
  std::vector<uint8_t> used_image(degree, 0u);

  for (std::size_t i = 0; i < degree; ++i) {
    const uint32_t pt = graph.points[i];
    const uint32_t x = pt & mask;
    const uint32_t y = pt >> n;
    const uint32_t x2 = field.Mul(a, x);
    const uint32_t y2 = field.Mul(b, y);
    const uint32_t img_pt = x2 | (y2 << n);
    if (img_pt >= graph.is_graph.size() || graph.is_graph[img_pt] == 0u) {
      return std::nullopt;
    }
    const auto img_idx = index.IndexOf(img_pt);
    if (!img_idx.has_value()) return std::nullopt;
    if (*img_idx >= degree) return std::nullopt;
    if (used_image[*img_idx] != 0u) return std::nullopt;
    used_image[*img_idx] = 1u;
    images[i] = static_cast<uint16_t>(*img_idx);
  }

  return Permutation(std::move(images));
}

}  // namespace

std::optional<Permutation> TryFrobeniusGenerator(const GraphData& graph,
                                                 const GF2n& field) {
  const int n = graph.n_bits;

  GraphPointIndex index(graph);

  const std::size_t degree = index.Size();
  if (degree != graph.points.size()) return std::nullopt;

  const uint32_t mask = (static_cast<uint32_t>(1u) << n) - 1u;
  std::vector<uint16_t> images(degree, 0u);
  std::vector<uint8_t> used_image(degree, 0u);

  for (std::size_t i = 0; i < degree; ++i) {
    const uint32_t pt = graph.points[i];
    const uint32_t x = pt & mask;
    const uint32_t y = pt >> n;
    const uint32_t x2 = field.Mul(x, x);
    const uint32_t y2 = field.Mul(y, y);
    const uint32_t img_pt = x2 | (y2 << n);
    if (img_pt >= graph.is_graph.size() || graph.is_graph[img_pt] == 0u) {
      return std::nullopt;
    }
    const auto img_idx = index.IndexOf(img_pt);
    if (!img_idx.has_value()) return std::nullopt;
    if (*img_idx >= degree) return std::nullopt;
    if (used_image[*img_idx] != 0u) return std::nullopt;
    used_image[*img_idx] = 1u;
    images[i] = static_cast<uint16_t>(*img_idx);
  }

  return Permutation(std::move(images));
}

std::vector<Permutation> FindSemilinearSeedGenerators(
    const GraphData& graph, const GF2n& field,
    std::size_t max_scaling_generators) {
  std::vector<Permutation> gens;

  const auto frob = TryFrobeniusGenerator(graph, field);
  if (frob.has_value() && !frob->IsIdentity()) gens.push_back(*frob);

  const int n = graph.n_bits;
  if (n <= 0 || n >= 31) return gens;
  if (graph.m_bits != n) return gens;
  if (field.n != n) return gens;
  if (max_scaling_generators == 0) return gens;

  std::vector<uint32_t> y_of_x;
  BuildYOfXTable(graph, n, y_of_x);

  const uint32_t alpha = ::FindPrimitiveElement(field);
  if (alpha != 0u && alpha != 1u) {
    uint32_t b = 0u;
    if (ComputeScalingBAndVerify(field, y_of_x, alpha, &b)) {
      const auto perm = BuildDiagonalPermutation(graph, field, n, alpha, b);
      if (perm.has_value() && !perm->IsIdentity()) gens.push_back(*perm);
    }
  }

  return gens;
}

}  // namespace groups
