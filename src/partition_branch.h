#ifndef PARTITION_BRANCH_H
#define PARTITION_BRANCH_H

#include <cstddef>
#include <optional>

#include "graph.h"
#include "ordered_partition.h"
#include "partial_map.h"

namespace partition_branch {

bool ApplyMatchedSingletonCells(const GraphData& F,
                                const OrderedPartition& points_left,
                                const OrderedPartition& points_right,
                                PartialAffineMap* A);

std::optional<std::size_t> NextUndeterminedBranchCell(
    const GraphData& F, const OrderedPartition& points_left,
    const OrderedPartition& points_right, const PartialAffineMap& A);

}  // namespace partition_branch

#endif
