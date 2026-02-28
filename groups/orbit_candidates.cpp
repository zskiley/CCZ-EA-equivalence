#include "orbit_candidates.h"

#include "orbit_traversal.h"

#include <cstddef>
#include <vector>

namespace groups {

std::vector<uint32_t> OrbitPrunedCandidatesExcludingTried(
    const GraphData& F, const std::vector<uint32_t>& right_cell,
    const PartialAffineMap& A,
    const std::vector<Permutation>& stabilizer_generators,
    const std::unordered_set<uint32_t>& tried_y_indices) {
  std::vector<uint32_t> candidates;
  candidates.reserve(right_cell.size());
  const std::size_t degree = F.points.size();

  for (uint32_t y_index : right_cell) {
    if (y_index >= F.points.size()) continue;
    const uint32_t y = F.points[y_index];
    if (A.HasImage(y)) continue;
    candidates.push_back(y_index);
  }

  if (candidates.empty()) return candidates;
  return OrbitRepresentativesInSubsetExcludingTouched(
      stabilizer_generators, candidates, tried_y_indices, degree);
}

}  // namespace groups
