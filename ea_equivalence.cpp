#include "ea_equivalence.h"

#include "groups/group_ops.h"
#include "hyperplane.h"
#include "ordered_partition.h"
#include "partial_map.h"

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <vector>

std::optional<EquivalencePointMap> RunEAEquivalence(
    const GraphData& F_left, const GraphData& F_right,
    const std::vector<groups::Permutation>& optimal_auto_group_generators,
    std::size_t min_active_hyperplanes) {
  if (F_left.d_bits != F_right.d_bits) return std::nullopt;
  if (F_left.n_bits != F_right.n_bits) return std::nullopt;
  if (F_left.m_bits != F_right.m_bits) return std::nullopt;
  if (F_left.points.size() != F_right.points.size()) return std::nullopt;
  if (F_left.points.empty()) return std::nullopt;

  const std::size_t point_count = F_left.points.size();
  if (min_active_hyperplanes == 0) {
    // Equivalence instances that are NOT equivalent can require deep search
    // if refinement is too weak. For small n, we can afford to keep a larger
    // hyperplane subset to strengthen refinement and reject earlier.
    min_active_hyperplanes = static_cast<std::size_t>(2u) * point_count;
    if (F_left.n_bits > 0 && F_left.n_bits <= 7) {
      min_active_hyperplanes = static_cast<std::size_t>(16u) * point_count;
    }
  }

  std::vector<uint32_t> point_indices(point_count, 0u);
  std::iota(point_indices.begin(), point_indices.end(), 0u);
  OrderedPartition points_left(point_indices);
  OrderedPartition points_right(point_indices);

  std::vector<Hyperplane> planes_left;
  OrderedPartition hyperplanes_left;
  BuildEAHyperplaneSubsetByWalsh(F_left, min_active_hyperplanes, &planes_left,
                                 &hyperplanes_left);

  std::vector<Hyperplane> planes_right;
  OrderedPartition hyperplanes_right;
  BuildEAHyperplaneSubsetByWalsh(F_right, min_active_hyperplanes, &planes_right,
                                 &hyperplanes_right);

  if (planes_left.empty() || planes_right.empty()) return std::nullopt;
  if (!hyperplanes_left.HasSameShape(hyperplanes_right)) return std::nullopt;

  PartialAffineMap A0(F_left.d_bits);
  ResetEquivalenceSearch();
  SetUseEaEquivalenceValidation(true);

  DfsEquivalenceGroupState root_group_state;
  root_group_state.H_generators.reserve(optimal_auto_group_generators.size());
  for (const auto& g : optimal_auto_group_generators) {
    if (g.Degree() != point_count) continue;
    if (g.IsIdentity()) continue;
    root_group_state.H_generators.push_back(g);
  }
  root_group_state.H_generators =
      groups::DeduplicateGenerators(std::move(root_group_state.H_generators));

  CCZ_DFS_equivalence(F_left, F_right, planes_left, planes_right, std::move(A0),
                      std::move(points_left), std::move(points_right),
                      std::move(hyperplanes_left), std::move(hyperplanes_right),
                      min_active_hyperplanes, std::move(root_group_state));

  const auto& found = GetFoundEquivalences();
  if (found.empty()) return std::nullopt;
  return found.front();
}
