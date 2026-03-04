#include "partition_branch.h"

#include <cstdint>
#include <optional>

namespace partition_branch {

bool ApplyMatchedSingletonCells(const GraphData& F,
                                const OrderedPartition& points_left,
                                const OrderedPartition& points_right,
                                PartialAffineMap* A) {
  if (A == nullptr) return false;

  const std::size_t n_cells = points_left.NumCells();
  for (std::size_t i = 0; i < n_cells; ++i) {
    const auto& left_cell = points_left.Cells()[i];
    const auto& right_cell = points_right.Cells()[i];
    if (left_cell.size() != 1) continue;

    const uint32_t x_index = left_cell[0];
    const uint32_t y_index = right_cell[0];
    const uint32_t x = F.points[x_index];
    const uint32_t y = F.points[y_index];

    const auto current = A->GetImage(x);
    if (current.has_value()) {
      if (*current != y) return false;
      continue;
    }

    if (A->HasImage(y)) return false;
    if (!A->Update(x, y)) return false;
  }
  return true;
}

std::optional<std::size_t> NextUndeterminedBranchCell(
    const GraphData& F, const OrderedPartition& points_left,
    const OrderedPartition& points_right, const PartialAffineMap& A) {
  (void)points_right;

  std::optional<std::size_t> best_idx;
  std::size_t best_size = 0;
  const std::size_t n_cells = points_left.NumCells();

  for (std::size_t i = 0; i < n_cells; ++i) {
    const auto& cell = points_left.Cells()[i];
    const std::size_t sz = cell.size();
    if (sz <= 1) continue;

    bool has_undetermined = false;
    for (uint32_t x_index : cell) {
      const uint32_t x = F.points[x_index];
      if (!A.GetImage(x).has_value()) {
        has_undetermined = true;
        break;
      }
    }
    if (!has_undetermined) continue;

    if (!best_idx.has_value() || sz < best_size) {
      best_idx = i;
      best_size = sz;
    }
  }

  return best_idx;
}

}  // namespace partition_branch
