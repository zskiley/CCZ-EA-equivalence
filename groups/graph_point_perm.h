#ifndef GROUPS_GRAPH_POINT_PERM_H
#define GROUPS_GRAPH_POINT_PERM_H

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "perm.h"

struct GraphData;

namespace groups {

using GraphPointMap = std::vector<std::pair<uint32_t, uint32_t>>;

class GraphPointIndex {
 public:
  explicit GraphPointIndex(const GraphData& graph);

  std::size_t Size() const;
  std::optional<uint32_t> IndexOf(uint32_t point) const;

 private:
  std::size_t size_ = 0;
  std::vector<int32_t> point_to_index_;
};

bool GraphPointMapToPermutation(const GraphPointMap& map,
                                const GraphPointIndex& index,
                                Permutation* out_permutation);

}  // namespace groups

#endif
