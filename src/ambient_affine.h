#ifndef AMBIENT_AFFINE_H
#define AMBIENT_AFFINE_H

#include <cstdint>
#include <vector>

#include "groups/perm.h"
#include "partial_map.h"

struct GraphData;

uint64_t KernelAffineOrder(const GraphData& F);
std::vector<AffineMapData> BuildKernelAffineGenerators(const GraphData& F);
bool LiftGraphPermutationToAffineMap(const GraphData& F,
                                     const groups::Permutation& perm,
                                     AffineMapData* out);
std::vector<AffineMapData> BuildAmbientAutoGenerators(
    const GraphData& F, const std::vector<groups::Permutation>& graph_generators);

#endif
