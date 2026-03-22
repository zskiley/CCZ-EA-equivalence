#include "algorithms.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "groups/graph_point_perm.h"
#include "groups/group_ops.h"

namespace algorithms {

namespace {

double DefaultAutoSeedTimeLimitSeconds(const GraphData& F) {
  const int n = F.n_bits;

  if (n <= 7) return 10.0;
  if (n == 8) return 20.0;
  if (n == 9) return 150.0;
  if (n == 10) return 150.0;
  if (n == 11) return 300.0;
  if (n == 12) return 300.0;
  if (n==13) return 1000;
  if (n==14) return 1000;
  if (n==15) return 3600;
  if (n==16) return 3600;
  return 90.0;
}

std::vector<groups::Permutation> ExtractAutoGroupGenerators(
    std::size_t degree) {
  std::vector<groups::Permutation> seed_generators = GetAutoGroupGenerators();
  seed_generators.erase(
      std::remove_if(seed_generators.begin(), seed_generators.end(),
                     [&](const groups::Permutation& g_perm) {
                       return g_perm.Degree() != degree || g_perm.IsIdentity();
                     }),
      seed_generators.end());
  return groups::DeduplicateGenerators(std::move(seed_generators));
}

}  // namespace

std::vector<GraphPointMap> ccz_auto(const GraphData& F, double timelimit_seconds,
                                    std::size_t min_active_hyperplanes) {
  return RunCCZAuto(F, min_active_hyperplanes, timelimit_seconds);
}

std::vector<GraphPointMap> ea_auto(const GraphData& F, double timelimit_seconds,
                                   std::size_t min_active_hyperplanes) {
  return RunEAAuto(F, min_active_hyperplanes, timelimit_seconds);
}

std::optional<EquivalencePointMap> ccz_equivalence(const GraphData& F,
                                                   const GraphData& G,
                                                   std::optional<double> timelimit_seconds,
                                                   std::size_t min_active_hyperplanes) {
  const double auto_timelimit = timelimit_seconds.has_value()
                                    ? *timelimit_seconds
                                    : DefaultAutoSeedTimeLimitSeconds(F);
  // Orbit pruning during equivalence search acts on the right-hand graph,
  // so seed with a subgroup of Aut(G), not Aut(F).
  (void)ccz_auto(G, auto_timelimit, min_active_hyperplanes);
  if (FoundEntireAutoGroup()) {
    std::cout << "(Found entire auto group)\n";
  } else {
    std::cout << "(potentially incoplete auto group)\n";
  }
  std::cout << "Auto group size before equivalence search: "
            << GetTotalAutoGroup() << "\n";

  std::vector<groups::Permutation> seed_generators =
      ExtractAutoGroupGenerators(F.points.size());
  return RunCCZEquivalence(F, G, seed_generators, min_active_hyperplanes);
}

std::optional<EquivalencePointMap> ea_equivalence(
    const GraphData& F, const GraphData& G,
    std::optional<double> timelimit_seconds,
    std::size_t min_active_hyperplanes) {
  const double auto_timelimit = timelimit_seconds.has_value()
                                    ? *timelimit_seconds
                                    : DefaultAutoSeedTimeLimitSeconds(F);
  // Orbit pruning during equivalence search acts on the right-hand graph,
  // so seed with a subgroup of Aut(G), not Aut(F).
  (void)ea_auto(G, auto_timelimit, min_active_hyperplanes);
  if (FoundEntireAutoGroup()) {
    std::cout << "(Found entire auto group)\n";
  } else {
    std::cout << "(potentially incoplete auto group)\n";
  }
  std::cout << "EA auto group size before equivalence search: "
            << GetTotalAutoGroup() << "\n";

  std::vector<groups::Permutation> seed_generators =
      ExtractAutoGroupGenerators(F.points.size());
  return RunEAEquivalence(F, G, seed_generators, min_active_hyperplanes);
}

}  // namespace algorithms
