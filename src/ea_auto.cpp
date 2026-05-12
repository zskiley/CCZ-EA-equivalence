#include "ea_auto.h"

#include "ccz_auto.h"

std::vector<GraphPointMap> RunEAAuto(const GraphData& F,
                                     std::size_t min_active_hyperplanes,
                                     double time_limit_seconds) {
  return RunAuto(F, /*ea_mode=*/true, min_active_hyperplanes,
                 time_limit_seconds);
}
