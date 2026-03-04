#ifndef DFS_AUTO_H
#define DFS_AUTO_H

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "groups/perm.h"
#include "graph.h"
#include "hyperplane.h"
#include "ordered_partition.h"
#include "partial_map.h"

using GraphPointMap = std::vector<std::pair<uint32_t, uint32_t>>;

struct DfsGroupState {
  std::vector<groups::Permutation> H_generators;
  std::vector<uint32_t> fixed_right_indices;
  uint64_t H_epoch = 0;
};

void ResetAutoSearch();
void SetAutoSearchTimeLimitSeconds(double seconds);
void SetUseEaValidation(bool enabled);
bool FoundEntireAutoGroup();
uint64_t GetTotalAutoGroup();
const std::vector<groups::Permutation>& GetAutoGroupGenerators();
const std::vector<GraphPointMap>& GetFoundAutos();
void InitializeGroupSearch(const GraphData& F);
void AddInitialGroupGenerators(std::vector<groups::Permutation> generators);
DfsGroupState MakeRootGroupState();

void CCZ_DFS_auto(const GraphData& F, const std::vector<Hyperplane>& planes,
                  PartialAffineMap A,
                  OrderedPartition points_left,
                  OrderedPartition points_right,
                  OrderedPartition hyperplanes_left,
                  OrderedPartition hyperplanes_right,
                  std::size_t min_active_hyperplanes,
                  DfsGroupState group_state);

#endif
