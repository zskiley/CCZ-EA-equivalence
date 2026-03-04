#include "schreier_sims.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace groups {
namespace {

uint64_t SaturatingMulU64(uint64_t a, uint64_t b) {
  if (a == 0u || b == 0u) return 0u;
  const uint64_t limit = std::numeric_limits<uint64_t>::max();
  if (a > limit / b) return limit;
  return a * b;
}

uint32_t FindNextBasePoint(const std::vector<Permutation>& generators,
                           const std::vector<uint8_t>& fixed_mask) {
  if (generators.empty()) return static_cast<uint32_t>(fixed_mask.size());
  const std::size_t degree = generators.front().Degree();

  for (uint32_t x = 0; x < degree; ++x) {
    if (fixed_mask[x] != 0u) continue;
    bool moved = false;
    for (const Permutation& g : generators) {
      if (g.Apply(x) != x) {
        moved = true;
        break;
      }
    }
    if (moved) return x;
  }
  return static_cast<uint32_t>(degree);
}

}  // namespace

SchreierSims::SchreierSims(std::size_t degree,
                           std::vector<Permutation> generators)
    : SchreierSims(degree) {
  SetGenerators(std::move(generators));
}

void SchreierSims::SetGenerators(std::vector<Permutation> generators) {
  generators_ = DeduplicateGenerators(std::move(generators));
  base_.clear();
  levels_.clear();
}

void SchreierSims::AddGenerator(const Permutation& generator) {
  generators_.push_back(generator);
  generators_ = DeduplicateGenerators(std::move(generators_));
  base_.clear();
  levels_.clear();
}

void SchreierSims::Build() {
  generators_ = DeduplicateGenerators(std::move(generators_));
  base_.clear();
  levels_.clear();

  std::vector<Permutation> current = generators_;
  std::vector<uint8_t> fixed_mask(degree_, 0u);

  while (!current.empty()) {
    const uint32_t beta = FindNextBasePoint(current, fixed_mask);
    if (beta >= degree_) break;

    OrbitTraversal orbit = BuildOrbitTraversal(current, beta, degree_);
    StabilizerLevel level;
    level.base_point = beta;
    level.orbit = std::move(orbit);
    level.generators = current;
    levels_.push_back(std::move(level));
    base_.push_back(beta);

    fixed_mask[beta] = 1u;

    std::vector<Permutation> schreier =
        BuildSchreierGenerators(current, levels_.back().orbit, degree_);
    // Schreier generators are already in the stabilizer of `beta`. Since every
    // element of `current` fixes all previously chosen base points, the new
    // stabilizer generators do as well. Filtering by earlier base points here is
    // redundant and expensive.
    current = ReduceGeneratorsByMembership(std::move(schreier), degree_);
  }
}

bool SchreierSims::Contains(const Permutation& element) const {
  if (element.Degree() != degree_) return false;
  if (levels_.empty()) return element.IsIdentity();

  Permutation residue = element;
  for (const StabilizerLevel& level : levels_) {
    const uint32_t image = residue.Apply(level.base_point);
    if (!level.orbit.Contains(image)) return false;

    Permutation rep;
    if (!RepresentativeToPoint(level.orbit, image, level.generators, &rep)) {
      return false;
    }
    residue = Compose(rep.Inverse(), residue);
  }
  return residue.IsIdentity();
}

uint64_t SchreierSims::Order() const {
  if (levels_.empty()) return 1u;
  uint64_t order = 1u;
  for (const StabilizerLevel& level : levels_) {
    order = SaturatingMulU64(order, static_cast<uint64_t>(level.orbit.OrbitSize()));
  }
  return order;
}

std::size_t SchreierSims::Degree() const {
  return degree_;
}

const std::vector<uint32_t>& SchreierSims::Base() const {
  return base_;
}

const std::vector<StabilizerLevel>& SchreierSims::Levels() const {
  return levels_;
}

std::vector<Permutation> ReduceGeneratorsByMembership(
    std::vector<Permutation> generators, std::size_t degree) {
  generators = DeduplicateGenerators(std::move(generators));
  if (generators.empty()) return generators;

  std::vector<Permutation> reduced;
  reduced.reserve(std::min<std::size_t>(generators.size(), 32u));

  SchreierSims ss(degree);
  ss.SetGenerators({});
  ss.Build();

  for (Permutation& g : generators) {
    if (g.IsIdentity()) continue;
    if (ss.Contains(g)) continue;
    reduced.push_back(std::move(g));
    ss.AddGenerator(reduced.back());
    ss.Build();
  }

  return reduced;
}

std::vector<Permutation> BuildPointwiseStabilizerGenerators(
    const std::vector<Permutation>& generators,
    const std::vector<uint32_t>& fixed_points, std::size_t degree) {
  if (!IsValidGenerators(generators, degree)) return {};

  std::vector<Permutation> gens = DeduplicateGenerators(generators);
  if (gens.empty() || fixed_points.empty()) return gens;

  // The pointwise stabilizer is order-independent; choosing the next fixed
  // point by smallest current orbit tends to keep Schreier sets smaller.
  std::vector<uint32_t> remaining = fixed_points;
  while (!remaining.empty() && !gens.empty()) {
    std::size_t best_pos = remaining.size();
    std::size_t best_orbit_size = static_cast<std::size_t>(-1);

    for (std::size_t i = 0; i < remaining.size(); ++i) {
      const uint32_t fixed = remaining[i];
      if (fixed >= degree) {
        gens.clear();
        return gens;
      }

      bool moved = false;
      for (const Permutation& g : gens) {
        if (g.Images()[static_cast<std::size_t>(fixed)] != fixed) {
          moved = true;
          break;
        }
      }
      if (!moved) {
        best_pos = i;
        best_orbit_size = 1;
        break;
      }

      std::size_t orbit_size = 0;
      try {
        orbit_size = BuildOrbitTraversal(gens, fixed, degree).OrbitSize();
      } catch (...) {
        gens.clear();
        return gens;
      }

      if (orbit_size < best_orbit_size ||
          (orbit_size == best_orbit_size &&
           (best_pos == remaining.size() || fixed < remaining[best_pos]))) {
        best_pos = i;
        best_orbit_size = orbit_size;
        if (best_orbit_size <= 2) break;
      }
    }

    if (best_pos >= remaining.size()) break;
    const uint32_t chosen = remaining[best_pos];
    remaining[best_pos] = remaining.back();
    remaining.pop_back();

    if (best_orbit_size <= 1) continue;

    try {
      OrbitTraversal orbit = BuildOrbitTraversal(gens, chosen, degree);
      gens = BuildSchreierGenerators(gens, orbit, degree);
      gens = ReduceGeneratorsByMembership(std::move(gens), degree);
    } catch (...) {
      gens.clear();
      return gens;
    }
  }

  return gens;
}

void RefreshPointwiseStabilizerFromGlobal(
    const std::vector<Permutation>& global_generators,
    const std::vector<uint32_t>& fixed_points, std::size_t degree,
    uint64_t global_epoch, std::vector<Permutation>* cached_generators,
    uint64_t* cached_epoch) {
  if (cached_generators == nullptr || cached_epoch == nullptr) return;
  if (*cached_epoch == global_epoch) return;
  if (!IsValidGenerators(global_generators, degree)) return;

  *cached_generators = BuildPointwiseStabilizerGenerators(
      global_generators, fixed_points, degree);
  *cached_epoch = global_epoch;
}

bool BuildChildPointStabilizerGenerators(
    const std::vector<Permutation>& parent_stabilizer_generators,
    uint32_t new_fixed_point, std::size_t degree,
    std::vector<Permutation>* child_stabilizer_generators) {
  if (child_stabilizer_generators == nullptr) return false;
  if (new_fixed_point >= degree) return false;
  if (!IsValidGenerators(parent_stabilizer_generators, degree)) return false;

  if (parent_stabilizer_generators.empty()) {
    child_stabilizer_generators->clear();
    return true;
  }

  try {
    OrbitTraversal orbit = BuildOrbitTraversal(parent_stabilizer_generators,
                                               new_fixed_point, degree);
    *child_stabilizer_generators = BuildSchreierGenerators(
        parent_stabilizer_generators, orbit, degree);
    *child_stabilizer_generators = ReduceGeneratorsByMembership(
        std::move(*child_stabilizer_generators), degree);
    return true;
  } catch (...) {
    // Keep DFS robust if orbit/Schreier construction fails unexpectedly.
    child_stabilizer_generators->clear();
    return true;
  }
}

}  // namespace groups
