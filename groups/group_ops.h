#ifndef GROUPS_GROUP_OPS_H
#define GROUPS_GROUP_OPS_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "orbit_traversal.h"
#include "perm.h"

namespace groups {

bool IsValidGenerators(const std::vector<Permutation>& generators,
                       std::size_t degree);

std::vector<Permutation> BuildSchreierGenerators(
    const std::vector<Permutation>& generators, const OrbitTraversal& orbit,
    std::size_t degree);

// Removes identity and duplicate permutations from `generators`.
// Consumes the input vector so implementations can move (not copy) large
// permutations when possible.
std::vector<Permutation> DeduplicateGenerators(std::vector<Permutation> generators);

std::vector<Permutation> FilterGeneratorsFixingPoints(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& fixed_points);

}  // namespace groups

#endif
