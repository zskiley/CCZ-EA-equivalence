#ifndef AFFINE_KERNEL_H
#define AFFINE_KERNEL_H

#include <cstdint>
#include <vector>

#include "partial_map.h"

struct GraphData;

uint64_t AffineKernelOrder(const GraphData& F, bool require_ea = false);
std::vector<AffineMapData> BuildAffineKernelGenerators(
    const GraphData& F, bool require_ea = false);

#endif
