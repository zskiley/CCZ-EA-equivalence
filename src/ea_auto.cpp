#include "ea_auto.h"

#include "dfs_auto.h"
#include "field_basics.h"
#include "groups/semilinear_seed.h"
#include "hyperplane.h"
#include "ordered_partition.h"
#include "quadratic.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

std::vector<GraphPointMap> RunEAAuto(const GraphData& F,
                                     std::size_t min_active_hyperplanes,
                                     double time_limit_seconds) {
  ResetAutoSearch();
  SetAutoSearchTimeLimitSeconds(time_limit_seconds);
  SetUseEaValidation(true);

  const std::size_t point_count = F.points.size();
  if (min_active_hyperplanes == 0) {
    min_active_hyperplanes = static_cast<std::size_t>(2u) * point_count;
  }

  std::vector<uint32_t> point_indices(point_count, 0u);
  std::iota(point_indices.begin(), point_indices.end(), 0u);
  OrderedPartition points_left(point_indices);
  OrderedPartition points_right(point_indices);

  std::vector<Hyperplane> planes;
  OrderedPartition hyperplanes_left;
  BuildEAHyperplaneSubsetByWalsh(F, min_active_hyperplanes, &planes,
                                 &hyperplanes_left);
  OrderedPartition hyperplanes_right = hyperplanes_left;

  PartialAffineMap A0(F.d_bits);

  InitializeGroupSearch(F);
  if (F.n_bits > 0 && F.n_bits < 31 && F.m_bits == F.n_bits) {
    const uint32_t mask = (static_cast<uint32_t>(1u) << F.n_bits) - 1u;
    const GF2n field{F.n_bits, DefaultModPoly(F.n_bits), mask};
    auto semilinear = groups::FindSemilinearSeedGenerators(
        F, field, /*max_scaling_generators=*/1);
    if (!semilinear.empty()) AddInitialGroupGenerators(std::move(semilinear));
  }

  DfsGroupState root_group_state = MakeRootGroupState();
  uint32_t quadratic_anchor_point = 0u;
  if (TryQuadraticAnchorPoint(F, &quadratic_anchor_point)) {
    uint32_t quadratic_anchor_index = 0u;
    (void)A0.Update(quadratic_anchor_point, quadratic_anchor_point);
    for (uint32_t i = 0; i < F.points.size(); ++i) {
      if (F.points[i] == quadratic_anchor_point) {
        quadratic_anchor_index = i;
        break;
      }
    }
    (void)points_left.Individualize(quadratic_anchor_index);
    (void)points_right.Individualize(quadratic_anchor_index);
    std::vector<groups::Permutation> translation_generators;
    BuildQuadraticTranslationGenerators(F, &translation_generators);
    AddInitialGroupGenerators(std::move(translation_generators));
    root_group_state = MakeRootGroupState();
    root_group_state.fixed_right_indices.push_back(quadratic_anchor_index);
    root_group_state.H_epoch = std::numeric_limits<uint64_t>::max();
  }

  CCZ_DFS_auto(F, planes, std::move(A0), std::move(points_left),
               std::move(points_right), std::move(hyperplanes_left),
               std::move(hyperplanes_right), min_active_hyperplanes,
               std::move(root_group_state));
  return GetFoundAutos();
}
