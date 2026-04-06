#include "hyperplane.h"

#include "dfs_helpers.h"
#include "graph.h"
#include "ordered_partition.h"
#include "weighted_fft.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

bool BuildHyperplaneSubsetByWalsh(const GraphData& F, std::size_t threshold,
                                  std::vector<Hyperplane>* planes,
                                  OrderedPartition* hyperplanes_partition) {
  if (planes == nullptr || hyperplanes_partition == nullptr) return false;
  planes->clear();
  *hyperplanes_partition = OrderedPartition();

  const int dim = F.d_bits;
  if (dim <= 0 || dim >= 31) return false;
  const std::size_t normals = static_cast<std::size_t>(1u) << dim;
  const std::size_t total_points = F.points.size();
  if (total_points == 0) return false;

  std::vector<WeightedPoint> sparse;
  sparse.reserve(total_points);
  for (uint32_t p : F.points) sparse.push_back({p, 1});
  std::vector<int64_t> walsh;
  WeightedFWHT(dim, sparse, &walsh);
  if (walsh.size() != normals) return false;

  if (threshold == 0) threshold = 1;
  const int total_i = static_cast<int>(total_points);
  auto counts_from_walsh = [total_i](int64_t w64, int* c0, int* c1) -> bool {
    if (c0 == nullptr || c1 == nullptr) return false;
    if (w64 < std::numeric_limits<int>::min() ||
        w64 > std::numeric_limits<int>::max()) {
      return false;
    }
    const int w = static_cast<int>(w64);
    const int v0 = total_i + w;
    const int v1 = total_i - w;
    if ((v0 & 1) != 0 || (v1 & 1) != 0) return false;
    const int local_c0 = v0 / 2;
    const int local_c1 = v1 / 2;
    if (local_c0 < 0 || local_c0 > total_i) return false;
    if (local_c1 < 0 || local_c1 > total_i) return false;
    *c0 = local_c0;
    *c1 = local_c1;
    return true;
  };

  std::vector<uint32_t> hist(total_points + 1u, 0u);
  for (std::size_t normal = 1; normal < normals; ++normal) {
    int c0 = 0;
    int c1 = 0;
    if (!counts_from_walsh(walsh[normal], &c0, &c1)) return false;
    ++hist[static_cast<std::size_t>(c0)];
    ++hist[static_cast<std::size_t>(c1)];
  }

  std::vector<uint8_t> keep;
  std::size_t keep_total = 0;
  if (!SelectSmallestColorClasses(hist, threshold, &keep, &keep_total)) {
    return false;
  }

  std::vector<OrderedPartition::Cell> cells;
  std::vector<int32_t> color_to_cell(keep.size(), -1);
  cells.reserve(keep.size());
  for (std::size_t color = 0; color < keep.size(); ++color) {
    if (keep[color] == 0u) continue;
    color_to_cell[color] = static_cast<int32_t>(cells.size());
    cells.emplace_back();
  }
  if (cells.empty()) return false;

  planes->reserve(keep_total);
  for (std::size_t normal = 1; normal < normals; ++normal) {
    int c0 = 0;
    int c1 = 0;
    if (!counts_from_walsh(walsh[normal], &c0, &c1)) return false;

    if (keep[static_cast<std::size_t>(c0)] != 0u) {
      const uint32_t plane_index = static_cast<uint32_t>(planes->size());
      planes->push_back(Hyperplane{static_cast<uint32_t>(normal), 0});
      const int32_t cell_id = color_to_cell[static_cast<std::size_t>(c0)];
      if (cell_id < 0) return false;
      cells[static_cast<std::size_t>(cell_id)].push_back(plane_index);
    }
    if (keep[static_cast<std::size_t>(c1)] != 0u) {
      const uint32_t plane_index = static_cast<uint32_t>(planes->size());
      planes->push_back(Hyperplane{static_cast<uint32_t>(normal), 1});
      const int32_t cell_id = color_to_cell[static_cast<std::size_t>(c1)];
      if (cell_id < 0) return false;
      cells[static_cast<std::size_t>(cell_id)].push_back(plane_index);
    }
  }

  std::vector<OrderedPartition::Cell> nonempty;
  nonempty.reserve(cells.size());
  for (auto& cell : cells) {
    if (!cell.empty()) nonempty.push_back(std::move(cell));
  }

  *hyperplanes_partition = OrderedPartition(std::move(nonempty), false);
  return !planes->empty() && hyperplanes_partition->NumCells() > 0;
}
