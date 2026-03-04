#include "orbit_traversal.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace groups {
namespace {

bool IsValidGeneratorsForOrbit(const std::vector<Permutation>& generators,
                               std::size_t degree) {
  for (const Permutation& g : generators) {
    if (g.Degree() != degree) return false;
  }
  return true;
}

}  // namespace

bool OrbitTraversal::Contains(uint32_t point) const {
  if (point >= point_to_orbit_index.size()) return false;
  return point_to_orbit_index[point] >= 0;
}

std::size_t OrbitTraversal::OrbitSize() const {
  return orbit_points.size();
}

OrbitTraversal BuildOrbitTraversal(const std::vector<Permutation>& generators,
                                   uint32_t base_point, std::size_t degree) {
  if (degree == 0) {
    throw std::invalid_argument("BuildOrbitTraversal: degree must be > 0");
  }
  if (base_point >= degree) {
    throw std::invalid_argument("BuildOrbitTraversal: base point out of range");
  }
  if (!IsValidGeneratorsForOrbit(generators, degree)) {
    throw std::invalid_argument(
        "BuildOrbitTraversal: generator degree mismatch");
  }

  OrbitTraversal orbit;
  orbit.base_point = base_point;
  orbit.point_to_orbit_index.assign(degree, -1);
  orbit.orbit_points.reserve(degree);
  orbit.parent_orbit_index.reserve(degree);
  orbit.parent_generator_index.reserve(degree);

  orbit.point_to_orbit_index[base_point] = 0;
  orbit.orbit_points.push_back(base_point);
  orbit.parent_orbit_index.push_back(-1);
  orbit.parent_generator_index.push_back(-1);

  std::size_t bfs_index = 0;
  while (bfs_index < orbit.orbit_points.size()) {
    const uint32_t alpha = orbit.orbit_points[bfs_index];
    for (std::size_t g_idx = 0; g_idx < generators.size(); ++g_idx) {
      const uint32_t beta = generators[g_idx].Apply(alpha);
      if (orbit.point_to_orbit_index[beta] >= 0) continue;
      orbit.point_to_orbit_index[beta] =
          static_cast<int32_t>(orbit.orbit_points.size());
      orbit.orbit_points.push_back(beta);
      orbit.parent_orbit_index.push_back(static_cast<int32_t>(bfs_index));
      orbit.parent_generator_index.push_back(static_cast<int32_t>(g_idx));
    }
    ++bfs_index;
  }

  return orbit;
}

bool RepresentativeToPoint(const OrbitTraversal& orbit, uint32_t point,
                           const std::vector<Permutation>& generators,
                           Permutation* representative) {
  if (representative == nullptr) return false;
  const std::size_t degree = orbit.point_to_orbit_index.size();
  if (point >= degree) return false;
  if (!IsValidGeneratorsForOrbit(generators, degree)) return false;

  const int32_t target_idx = orbit.point_to_orbit_index[point];
  if (target_idx < 0) return false;

  int32_t idx = target_idx;
  std::vector<int32_t> word;
  while (idx > 0) {
    if (idx >= static_cast<int32_t>(orbit.parent_orbit_index.size())) return false;
    const int32_t gen_idx = orbit.parent_generator_index[idx];
    const int32_t parent_idx = orbit.parent_orbit_index[idx];
    if (gen_idx < 0 || gen_idx >= static_cast<int32_t>(generators.size())) {
      return false;
    }
    if (parent_idx < 0 || parent_idx >= idx) return false;
    word.push_back(gen_idx);
    idx = parent_idx;
  }

  Permutation rep = Permutation::Identity(degree);
  for (auto it = word.rbegin(); it != word.rend(); ++it) {
    rep = Compose(generators[*it], rep);
  }
  *representative = std::move(rep);
  return true;
}

std::optional<uint32_t> ChooseOrbitRepresentativeByMaxSubsetOrbit(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& subset, std::size_t degree) {
  if (subset.empty()) return std::nullopt;
  if (generators.empty() || !IsValidGeneratorsForOrbit(generators, degree)) {
    return subset.front();
  }

  std::vector<uint8_t> in_subset(degree, 0u);
  for (uint32_t idx : subset) {
    if (idx < degree) in_subset[idx] = 1u;
  }

  std::vector<uint8_t> visited(degree, 0u);
  std::size_t best_orbit_size = 0;
  uint32_t best_rep_idx = subset.front();

  for (uint32_t idx : subset) {
    if (idx >= degree) continue;
    if (visited[idx] != 0u) continue;

    OrbitTraversal orbit;
    try {
      orbit = BuildOrbitTraversal(generators, idx, degree);
    } catch (...) {
      return subset.front();
    }

    std::size_t orbit_size_in_subset = 0;
    uint32_t orbit_rep_idx = idx;
    for (uint32_t z : orbit.orbit_points) {
      if (z >= degree) continue;
      if (in_subset[z] == 0u) continue;
      if (visited[z] == 0u) ++orbit_size_in_subset;
      visited[z] = 1u;
      if (z < orbit_rep_idx) orbit_rep_idx = z;
    }
    if (orbit_size_in_subset == 0) orbit_size_in_subset = 1;

    if (orbit_size_in_subset > best_orbit_size ||
        (orbit_size_in_subset == best_orbit_size &&
         orbit_rep_idx < best_rep_idx)) {
      best_orbit_size = orbit_size_in_subset;
      best_rep_idx = orbit_rep_idx;
    }
  }

  if (best_rep_idx >= degree) return subset.front();
  return best_rep_idx;
}

std::vector<uint32_t> OrbitRepresentativesInSubset(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& subset, std::size_t degree) {
  if (subset.empty()) return {};
  if (generators.empty() || !IsValidGeneratorsForOrbit(generators, degree)) {
    return subset;
  }

  std::vector<uint8_t> in_subset(degree, 0u);
  for (uint32_t idx : subset) {
    if (idx < degree) in_subset[idx] = 1u;
  }
  std::vector<uint8_t> visited(degree, 0u);

  std::vector<uint32_t> representatives;
  representatives.reserve(subset.size());
  for (uint32_t idx : subset) {
    if (idx >= degree) continue;
    if (visited[idx] != 0u) continue;

    representatives.push_back(idx);
    OrbitTraversal orbit;
    try {
      orbit = BuildOrbitTraversal(generators, idx, degree);
    } catch (...) {
      return subset;
    }
    for (uint32_t z : orbit.orbit_points) {
      if (z >= degree) continue;
      if (in_subset[z] == 0u) continue;
      visited[z] = 1u;
    }
    visited[idx] = 1u;
  }

  return representatives;
}

std::vector<uint32_t> OrbitRepresentativesInSubsetExcludingTouched(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& subset,
    const std::unordered_set<uint32_t>& excluded_points, std::size_t degree) {
  if (subset.empty()) return {};
  if (generators.empty() || !IsValidGeneratorsForOrbit(generators, degree)) {
    std::vector<uint32_t> out;
    out.reserve(subset.size());
    for (uint32_t idx : subset) {
      if (excluded_points.find(idx) == excluded_points.end()) out.push_back(idx);
    }
    return out;
  }

  std::vector<uint8_t> in_subset(degree, 0u);
  for (uint32_t idx : subset) {
    if (idx < degree) in_subset[idx] = 1u;
  }

  std::vector<uint8_t> in_excluded(degree, 0u);
  for (uint32_t idx : excluded_points) {
    if (idx < degree) in_excluded[idx] = 1u;
  }

  std::vector<uint8_t> visited(degree, 0u);
  std::vector<uint32_t> representatives;
  representatives.reserve(subset.size());

  for (uint32_t idx : subset) {
    if (idx >= degree) continue;
    if (visited[idx] != 0u) continue;

    OrbitTraversal orbit;
    try {
      orbit = BuildOrbitTraversal(generators, idx, degree);
    } catch (...) {
      std::vector<uint32_t> out;
      out.reserve(subset.size());
      for (uint32_t y_index : subset) {
        if (excluded_points.find(y_index) == excluded_points.end()) {
          out.push_back(y_index);
        }
      }
      return out;
    }

    bool orbit_has_excluded = false;
    uint32_t rep_idx = idx;
    bool rep_set = false;
    for (uint32_t z : orbit.orbit_points) {
      if (z >= degree) continue;
      if (in_subset[z] == 0u) continue;
      visited[z] = 1u;
      if (in_excluded[z] != 0u) orbit_has_excluded = true;
      if (!rep_set || z < rep_idx) {
        rep_idx = z;
        rep_set = true;
      }
    }
    if (!rep_set) continue;
    if (orbit_has_excluded) continue;
    representatives.push_back(rep_idx);
  }

  return representatives;
}

}  // namespace groups
