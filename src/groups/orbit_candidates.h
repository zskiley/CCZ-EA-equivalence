#ifndef GROUPS_ORBIT_CANDIDATES_H
#define GROUPS_ORBIT_CANDIDATES_H

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "../graph.h"
#include "../partial_map.h"
#include "perm.h"

namespace groups {

std::vector<uint32_t> OrbitPrunedCandidatesExcludingTried(
    const GraphData& F, const std::vector<uint32_t>& right_cell,
    const PartialAffineMap& A,
    const std::vector<Permutation>& stabilizer_generators,
    const std::unordered_set<uint32_t>& tried_y_indices);

}  // namespace groups

#endif
