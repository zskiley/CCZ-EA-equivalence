#ifndef CCZ_EQUIVALENCE_H
#define CCZ_EQUIVALENCE_H

#include <cstddef>
#include <optional>
#include <vector>

#include "dfs_equivalence.h"

std::optional<EquivalencePointMap> RunCCZEquivalence(
    const GraphData& F_left, const GraphData& F_right,
    const std::vector<groups::Permutation>& optimal_auto_group_generators,
    std::size_t min_active_hyperplanes = 0);

#endif
