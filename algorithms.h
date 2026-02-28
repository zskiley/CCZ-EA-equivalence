#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <cstddef>
#include <optional>
#include <vector>

#include "ccz_equivalence.h"
#include "ccz_auto.h"
#include "ea_auto.h"
#include "ea_equivalence.h"

namespace algorithms {

std::vector<GraphPointMap> ccz_auto(const GraphData& F,
                                    double timelimit_seconds = 0.0,
                                    std::size_t min_active_hyperplanes = 0);

std::vector<GraphPointMap> ea_auto(const GraphData& F,
                                   double timelimit_seconds = 0.0,
                                   std::size_t min_active_hyperplanes = 0);

std::optional<EquivalencePointMap> ccz_equivalence(const GraphData& F,
                                                   const GraphData& G,
                                                   std::optional<double> timelimit_seconds = std::nullopt,
                                                   std::size_t min_active_hyperplanes = 0);

std::optional<EquivalencePointMap> ea_equivalence(const GraphData& F,
                                                   const GraphData& G,
                                                   std::optional<double> timelimit_seconds = std::nullopt,
                                                   std::size_t min_active_hyperplanes = 0);

}  // namespace algorithms

#endif
