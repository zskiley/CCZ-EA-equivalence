#include "graph_point_perm.h"

#include "../graph.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace groups {

GraphPointIndex::GraphPointIndex(const GraphData& graph) {
  if (graph.points.empty()) {
    size_ = 0;
    return;
  }

  size_ = graph.points.size();
  if (graph.is_graph.empty()) {
    return;
  }

  point_to_index_.assign(graph.is_graph.size(), -1);
  for (std::size_t i = 0; i < graph.points.size(); ++i) {
    const uint32_t p = graph.points[i];
    if (p >= point_to_index_.size()) {
      point_to_index_.clear();
      return;
    }
    if (point_to_index_[p] >= 0) {
      point_to_index_.clear();
      return;
    }
    point_to_index_[p] = static_cast<int32_t>(i);
  }
}

std::size_t GraphPointIndex::Size() const {
  return size_;
}

std::optional<uint32_t> GraphPointIndex::IndexOf(uint32_t point) const {
  if (point >= point_to_index_.size()) return std::nullopt;
  const int32_t idx = point_to_index_[point];
  if (idx < 0) return std::nullopt;
  return static_cast<uint32_t>(idx);
}

bool GraphPointMapToPermutation(const GraphPointMap& map,
                                const GraphPointIndex& index,
                                Permutation* out_permutation) {
  if (out_permutation == nullptr) return false;

  const std::size_t n = index.Size();
  std::vector<int32_t> image_by_domain(n, -1);

  for (const auto& kv : map) {
    const auto x_index = index.IndexOf(kv.first);
    const auto y_index = index.IndexOf(kv.second);
    if (!x_index.has_value() || !y_index.has_value()) return false;
    if (*x_index >= n || *y_index >= n) return false;

    int32_t& slot = image_by_domain[*x_index];
    if (slot >= 0 && slot != static_cast<int32_t>(*y_index)) return false;
    slot = static_cast<int32_t>(*y_index);
  }

  std::vector<uint8_t> used_image(n, 0u);
  std::vector<uint16_t> images(n, 0u);
  for (std::size_t i = 0; i < n; ++i) {
    const int32_t y = image_by_domain[i];
    if (y < 0) return false;
    if (y >= static_cast<int32_t>(n)) return false;
    if (used_image[y] != 0u) return false;
    used_image[y] = 1u;
    images[i] = static_cast<uint16_t>(y);
  }

  *out_permutation = Permutation(std::move(images));
  return true;
}

}  // namespace groups
