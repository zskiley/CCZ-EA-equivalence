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

void PrintAutoSeedStatus(const char* label_prefix,
                         bool found_entire_group,
                         uint64_t total_auto_group,
                         const std::vector<groups::Permutation>& seed_generators,
                         double auto_timelimit) {
  if (found_entire_group) {
    std::cout << "(Found entire auto group)\n";
  } else {
    std::cout << "(potentially incoplete auto group)\n";
  }
  std::cout << label_prefix << total_auto_group << "\n";
  if (!found_entire_group && seed_generators.empty()) {
    std::cout << "Consider increasing the time limit for the automorphism "
                 "search, time_limit_seconds = "
              << auto_timelimit << "\n";
  }
}

void PrintEquivalenceSearchStart(const char* label_prefix,
                                 uint64_t total_auto_group) {
  std::cout << label_prefix << total_auto_group << "\n";
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
                                    : 0.0;
  // Orbit pruning during equivalence search acts on the right-hand graph,
  // so seed with a subgroup of Aut(G), not Aut(F).
  (void)ccz_auto(G, auto_timelimit, min_active_hyperplanes);

  std::vector<groups::Permutation> seed_generators =
      ExtractAutoGroupGenerators(F.points.size());
  PrintAutoSeedStatus("Auto group size before equivalence search: ",
                      FoundEntireAutoGroup(), GetTotalAutoGroup(),
                      seed_generators, auto_timelimit);
  PrintEquivalenceSearchStart("Starting equivalence search with auto group size: ",
                              GetTotalAutoGroup());
  return RunCCZEquivalence(F, G, seed_generators, min_active_hyperplanes);
}

std::optional<EquivalencePointMap> ea_equivalence(
    const GraphData& F, const GraphData& G,
    std::optional<double> timelimit_seconds,
    std::size_t min_active_hyperplanes) {
  const double auto_timelimit = timelimit_seconds.has_value()
                                    ? *timelimit_seconds
                                    : 0.0;
  // Orbit pruning during equivalence search acts on the right-hand graph,
  // so seed with a subgroup of Aut(G), not Aut(F).
  (void)ea_auto(G, auto_timelimit, min_active_hyperplanes);

  std::vector<groups::Permutation> seed_generators =
      ExtractAutoGroupGenerators(F.points.size());
  PrintAutoSeedStatus("EA auto group size before equivalence search: ",
                      FoundEntireAutoGroup(), GetTotalAutoGroup(),
                      seed_generators, auto_timelimit);
  PrintEquivalenceSearchStart(
      "Starting EA equivalence search with auto group size: ",
      GetTotalAutoGroup());
  return RunEAEquivalence(F, G, seed_generators, min_active_hyperplanes);
}

}  // namespace algorithms
