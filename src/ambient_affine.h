#ifndef AMBIENT_AFFINE_H
#define AMBIENT_AFFINE_H

#include <cstdint>
#include <vector>

#include "partial_map.h"

struct GraphData;

uint64_t KernelAffineOrder(const GraphData& F, bool require_ea = false);
std::vector<AffineMapData> BuildKernelAffineGenerators(const GraphData& F,
                                                       bool require_ea = false);

#endif
