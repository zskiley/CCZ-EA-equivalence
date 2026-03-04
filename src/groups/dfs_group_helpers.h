#ifndef GROUPS_DFS_GROUP_HELPERS_H
#define GROUPS_DFS_GROUP_HELPERS_H

#include <cstdint>
#include <optional>
#include <vector>

#include "graph_point_perm.h"
#include "perm.h"
#include "schreier_sims.h"

struct DfsGroupState;

namespace groups {

void AddInitialGroupGeneratorsImpl(
    const std::optional<GraphPointIndex>& graph_point_index,
    std::optional<SchreierSims>* global_group,
    std::vector<Permutation>* global_generators, uint64_t* group_epoch,
    std::vector<Permutation> generators);

void FillRootGroupStateImpl(const std::vector<Permutation>& global_generators,
                            uint64_t group_epoch, DfsGroupState* out_state);

}  // namespace groups

#endif
