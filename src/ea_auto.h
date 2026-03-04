#ifndef EA_AUTO_H
#define EA_AUTO_H

#include <cstddef>
#include <vector>

#include "dfs_auto.h"

std::vector<GraphPointMap> RunEAAuto(
    const GraphData& F, std::size_t min_active_hyperplanes = 0,
    double time_limit_seconds = 0.0);

#endif
