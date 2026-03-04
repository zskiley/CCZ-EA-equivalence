#ifndef GROUPS_SEMILINEAR_SEED_H
#define GROUPS_SEMILINEAR_SEED_H

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "perm.h"

struct GraphData;
struct GF2n;

namespace groups {

// Tries to verify that the Frobenius map (x,y) -> (x^2, y^2) (in GF(2^n))
// preserves `graph` and returns the induced permutation on graph-point indices.
std::optional<Permutation> TryFrobeniusGenerator(const GraphData& graph,
                                                const GF2n& field);

// Finds a small verified set of "semilinear-style" seed generators:
// - Frobenius (if it preserves the graph)
// - one primitive diagonal scaling (a,b), where a is a primitive field element
//   and b is solved from f(a*x)=b*f(x), if it preserves the graph.
//
std::vector<Permutation> FindSemilinearSeedGenerators(
    const GraphData& graph, const GF2n& field,
    std::size_t max_scaling_generators);

}  // namespace groups

#endif  // GROUPS_SEMILINEAR_SEED_H
