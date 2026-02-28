#include "group_ops.h"

#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace groups {
namespace {

bool FixesAll(const Permutation& g, const std::vector<uint32_t>& fixed_points) {
  for (uint32_t x : fixed_points) {
    if (x >= g.Degree()) return false;
    if (g.Apply(x) != x) return false;
  }
  return true;
}

}  // namespace

bool IsValidGenerators(const std::vector<Permutation>& generators,
                       std::size_t degree) {
  for (const Permutation& g : generators) {
    if (g.Degree() != degree) return false;
  }
  return true;
}

std::vector<Permutation> BuildSchreierGenerators(
    const std::vector<Permutation>& generators, const OrbitTraversal& orbit,
    std::size_t degree) {
  if (!IsValidGenerators(generators, degree)) {
    throw std::invalid_argument(
        "BuildSchreierGenerators: generator degree mismatch");
  }

  std::vector<Permutation> reps;
  reps.reserve(orbit.orbit_points.size());
  for (uint32_t point : orbit.orbit_points) {
    Permutation rep;
    if (!RepresentativeToPoint(orbit, point, generators, &rep)) {
      throw std::invalid_argument(
          "BuildSchreierGenerators: failed to recover transversal element");
    }
    reps.push_back(std::move(rep));
  }

  std::vector<Permutation> reps_inv;
  reps_inv.reserve(reps.size());
  for (const auto& rep : reps) {
    reps_inv.push_back(rep.Inverse());
  }

  std::vector<Permutation> schreier;
  schreier.reserve(orbit.orbit_points.size() * generators.size());
  for (std::size_t orbit_idx = 0; orbit_idx < orbit.orbit_points.size();
       ++orbit_idx) {
    const uint32_t alpha = orbit.orbit_points[orbit_idx];
    for (const Permutation& s : generators) {
      const uint32_t beta = s.Apply(alpha);
      const int32_t beta_idx = orbit.point_to_orbit_index[beta];
      if (beta_idx < 0) continue;
      Permutation h = Compose3(reps_inv[beta_idx], s, reps[orbit_idx]);
      if (h.IsIdentity()) continue;
      schreier.push_back(std::move(h));
    }
  }

  return DeduplicateGenerators(std::move(schreier));
}

std::vector<Permutation> DeduplicateGenerators(
    std::vector<Permutation> generators) {
  std::vector<Permutation> out;
  out.reserve(generators.size());

  // Avoid copying full permutations into an unordered_set key. We deduplicate
  // by hashing into buckets of indices into `out`, and only do a full equality
  // compare on hash collisions.
  std::unordered_map<std::size_t, std::vector<std::size_t>> buckets;
  buckets.reserve(generators.size() * 2u + 1u);
  PermutationHash hasher;

  for (Permutation& g : generators) {
    if (g.IsIdentity()) continue;
    const std::size_t h = hasher(g);
    auto& idxs = buckets[h];
    bool dup = false;
    for (std::size_t idx : idxs) {
      if (idx < out.size() && out[idx] == g) {
        dup = true;
        break;
      }
    }
    if (dup) continue;
    idxs.push_back(out.size());
    out.push_back(std::move(g));
  }
  return out;
}

std::vector<Permutation> FilterGeneratorsFixingPoints(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& fixed_points) {
  std::vector<Permutation> filtered;
  filtered.reserve(generators.size());
  for (const Permutation& g : generators) {
    if (FixesAll(g, fixed_points)) {
      filtered.push_back(g);
    }
  }
  return DeduplicateGenerators(std::move(filtered));
}

}  // namespace groups
