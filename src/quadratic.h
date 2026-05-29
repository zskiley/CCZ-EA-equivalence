#ifndef QUADRATIC_H
#define QUADRATIC_H

#include <cstdint>
#include <vector>

#include "groups/perm.h"
#include "partial_map.h"

struct GraphData;

bool TryQuadraticAnchorPoint(const GraphData& F, uint32_t* anchor_point);

bool BuildQuadraticTranslationGenerators(
    const GraphData& F, std::vector<groups::Permutation>* generators);

bool BuildQuadraticTranslationAffineGenerators(
    const GraphData& F, std::vector<AffineMapData>* generators);

bool BuildQuadraticTranslationData(
    const GraphData& F, uint32_t* anchor_point,
    std::vector<groups::Permutation>* translation_generators,
    std::vector<AffineMapData>* translation_affine_generators);

#endif
