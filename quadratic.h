#ifndef QUADRATIC_H
#define QUADRATIC_H

#include <cstdint>
#include <vector>

#include "groups/perm.h"

struct GraphData;

bool TryQuadraticAnchorPoint(const GraphData& F, uint32_t* anchor_point);

bool BuildQuadraticTranslationGenerators(
    const GraphData& F, std::vector<groups::Permutation>* generators);

#endif
