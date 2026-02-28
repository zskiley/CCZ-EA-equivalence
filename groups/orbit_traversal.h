#ifndef GROUPS_ORBIT_TRAVERSAL_H
#define GROUPS_ORBIT_TRAVERSAL_H

#include <cstddef>
#include <cstdint>
#include <optional>
#include <unordered_set>
#include <vector>

#include "perm.h"

namespace groups {

struct OrbitTraversal {
  uint32_t base_point = 0u;
  std::vector<uint32_t> orbit_points;
  std::vector<int32_t> point_to_orbit_index;
  std::vector<int32_t> parent_orbit_index;
  std::vector<int32_t> parent_generator_index;

  bool Contains(uint32_t point) const;
  std::size_t OrbitSize() const;
};

OrbitTraversal BuildOrbitTraversal(const std::vector<Permutation>& generators,
                                   uint32_t base_point, std::size_t degree);

bool RepresentativeToPoint(const OrbitTraversal& orbit, uint32_t point,
                           const std::vector<Permutation>& generators,
                           Permutation* representative);

// Picks one representative from `subset` using H-orbits: maximize orbit size
// inside `subset`, tie-breaking by smallest representative index.
std::optional<uint32_t> ChooseOrbitRepresentativeByMaxSubsetOrbit(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& subset, std::size_t degree);

// Returns one representative per H-orbit that intersects `subset`.
std::vector<uint32_t> OrbitRepresentativesInSubset(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& subset, std::size_t degree);

// Like OrbitRepresentativesInSubset, but drops any orbit that contains at
// least one point in `excluded_points`. Representative is the smallest index
// in the subset-intersection of the orbit.
std::vector<uint32_t> OrbitRepresentativesInSubsetExcludingTouched(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& subset,
    const std::unordered_set<uint32_t>& excluded_points, std::size_t degree);

}  // namespace groups

#endif
