#ifndef CCZ_AUTO_H
#define CCZ_AUTO_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "dfs_auto.h"

std::vector<GraphPointMap> RunCCZAuto(
    const GraphData& F, std::size_t min_active_hyperplanes = 0,
    double time_limit_seconds = 0.0);

#endif
