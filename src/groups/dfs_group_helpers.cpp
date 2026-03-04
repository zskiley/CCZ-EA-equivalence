#include "dfs_group_helpers.h"

#include "../dfs_auto.h"

#include <cstddef>
#include <utility>
#include <vector>

namespace groups {

void AddInitialGroupGeneratorsImpl(
    const std::optional<GraphPointIndex>& graph_point_index,
    std::optional<SchreierSims>* global_group,
    std::vector<Permutation>* global_generators, uint64_t* group_epoch,
    std::vector<Permutation> generators) {
  if (!graph_point_index.has_value() || global_group == nullptr ||
      global_generators == nullptr || group_epoch == nullptr ||
      !global_group->has_value()) {
    return;
  }

  const std::size_t degree = graph_point_index->Size();
  if (!IsValidGenerators(generators, degree)) return;

  bool changed = false;
  for (const Permutation& g : generators) {
    if (g.IsIdentity()) continue;
    if ((*global_group)->Contains(g)) continue;
    global_generators->push_back(g);
    changed = true;
  }
  if (!changed) return;

  *global_generators = DeduplicateGenerators(std::move(*global_generators));
  *global_generators =
      ReduceGeneratorsByMembership(std::move(*global_generators), degree);
  (*global_group)->SetGenerators(*global_generators);
  (*global_group)->Build();
  ++(*group_epoch);
}

void FillRootGroupStateImpl(const std::vector<Permutation>& global_generators,
                            uint64_t group_epoch, DfsGroupState* out_state) {
  if (out_state == nullptr) return;
  out_state->H_generators = global_generators;
  out_state->H_epoch = group_epoch;
}

}  // namespace groups
