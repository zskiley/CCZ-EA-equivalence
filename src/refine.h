#ifndef REFINE_H
#define REFINE_H

#include <cstddef>
#include <vector>

#include "ordered_partition.h"

struct GraphData;
struct Hyperplane;

// Runs FWHT-based partition refinement to a fixpoint for paired
// left/right point and hyperplane partitions.
// For automorphisms use identical left/right inputs.
// Returns false on invalid input or shape mismatch during refinement.
bool RefineToFixpoint(const GraphData& F_left, const GraphData& F_right,
                      const std::vector<Hyperplane>& planes_left,
                      const std::vector<Hyperplane>& planes_right,
                      OrderedPartition* points_left,
                      OrderedPartition* points_right,
                      OrderedPartition* hyperplanes_left,
                      OrderedPartition* hyperplanes_right,
                      std::size_t min_active_hyperplanes);

#endif
