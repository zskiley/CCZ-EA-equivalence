#ifndef DFS_EQUIVALENCE_H
#define DFS_EQUIVALENCE_H

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "groups/perm.h"
#include "graph.h"
#include "hyperplane.h"
#include "ordered_partition.h"
#include "partial_map.h"

using EquivalencePointMap = std::vector<std::pair<uint32_t, uint32_t>>;

struct DfsEquivalenceGroupState {
  std::vector<groups::Permutation> H_generators;
  std::vector<uint32_t> fixed_right_indices;
};

void ResetEquivalenceSearch();
void SetUseEaEquivalenceValidation(bool enabled);
const std::vector<EquivalencePointMap>& GetFoundEquivalences();

void CCZ_DFS_equivalence(const GraphData& F_left, const GraphData& F_right,
                         const std::vector<Hyperplane>& planes_left,
                         const std::vector<Hyperplane>& planes_right,
                         PartialAffineMap A,
                         OrderedPartition points_left,
                         OrderedPartition points_right,
                         OrderedPartition hyperplanes_left,
                         OrderedPartition hyperplanes_right,
                         std::size_t min_active_hyperplanes,
                         DfsEquivalenceGroupState group_state);

#endif
