#include "dfs_equivalence.h"

#include "groups/orbit_candidates.h"
#include "groups/schreier_sims.h"
#include "refine.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <unordered_set>
#include <vector>

namespace {

std::vector<EquivalencePointMap> g_found_equivalences;
bool g_use_ea_validation = false;

bool IsValidCCZBetween(const PartialAffineMap& A, const GraphData& F_left,
                       const GraphData& F_right) {
  if (F_left.d_bits != F_right.d_bits) return false;
  if (F_left.d_bits <= 0 || F_left.d_bits >= 31) return false;
  if (F_right.is_graph.size() != (static_cast<std::size_t>(1u) << F_right.d_bits)) {
    return false;
  }
  if (F_left.is_graph.size() != (static_cast<std::size_t>(1u) << F_left.d_bits)) {
    return false;
  }

  for (const auto& kv : A.Mappings()) {
    const uint32_t x = kv.first;
    const uint32_t y = kv.second;
    if (x >= F_left.is_graph.size() || !F_left.is_graph[x]) return false;
    if (y >= F_right.is_graph.size() || !F_right.is_graph[y]) return false;
  }

  for (uint32_t x : F_left.points) {
    const auto y = A.GetImage(x);
    if (!y.has_value()) continue;
    if (*y >= F_right.is_graph.size()) return false;
    if (!F_right.is_graph[*y]) return false;
  }

  return true;
}

bool IsValidEABetween(const PartialAffineMap& A, const GraphData& F_left,
                      const GraphData& F_right) {
  if (!IsValidCCZBetween(A, F_left, F_right)) return false;
  if (F_left.n_bits != F_right.n_bits) return false;
  if (F_left.m_bits != F_right.m_bits) return false;

  const int n = F_left.n_bits;
  if (n <= 0 || n >= 31) return false;
  const uint32_t mask_x = (1u << n) - 1u;

  std::vector<int32_t> left_to_right(static_cast<std::size_t>(1u << n), -1);
  std::vector<int32_t> right_owner(static_cast<std::size_t>(1u << n), -1);

  for (const auto& kv : A.Mappings()) {
    const uint32_t x_left = kv.first & mask_x;
    const uint32_t y_left = kv.second & mask_x;

    const int32_t previous = left_to_right[x_left];
    if (previous >= 0 && static_cast<uint32_t>(previous) != y_left) return false;
    left_to_right[x_left] = static_cast<int32_t>(y_left);

    const int32_t owner = right_owner[y_left];
    if (owner >= 0 && static_cast<uint32_t>(owner) != x_left) return false;
    right_owner[y_left] = static_cast<int32_t>(x_left);
  }

  return true;
}

bool IsValidMapByMode(const PartialAffineMap& A, const GraphData& F_left,
                      const GraphData& F_right) {
  return g_use_ea_validation ? IsValidEABetween(A, F_left, F_right)
                             : IsValidCCZBetween(A, F_left, F_right);
}

bool TryFinalizeCurrentMap(const GraphData& F_left, const GraphData& F_right,
                           const PartialAffineMap& A) {
  if (!IsValidMapByMode(A, F_left, F_right)) return false;

  const int n_points = static_cast<int>(F_left.points.size());
  EquivalencePointMap map;
  map.reserve(n_points);
  for (int i = 0; i < n_points; ++i) {
    const uint32_t x = F_left.points[static_cast<std::size_t>(i)];
    const auto y = A.GetImage(x);
    if (!y.has_value()) return false;
    if (*y >= F_right.is_graph.size()) return false;
    if (!F_right.is_graph[*y]) return false;
    map.push_back({x, *y});
  }
  g_found_equivalences.push_back(std::move(map));
  return true;
}

bool ApplyMatchedSingletonCellsPair(const GraphData& F_left,
                                    const GraphData& F_right,
                                    const OrderedPartition& points_left,
                                    const OrderedPartition& points_right,
                                    PartialAffineMap* A) {
  if (A == nullptr) return false;
  if (!points_left.HasSameShape(points_right)) return false;

  const std::size_t n_cells = points_left.NumCells();
  for (std::size_t i = 0; i < n_cells; ++i) {
    const auto& left_cell = points_left.Cells()[i];
    const auto& right_cell = points_right.Cells()[i];
    if (left_cell.size() != 1) continue;

    const uint32_t x_index = left_cell[0];
    const uint32_t y_index = right_cell[0];
    if (x_index >= F_left.points.size() || y_index >= F_right.points.size()) {
      return false;
    }

    const uint32_t x = F_left.points[x_index];
    const uint32_t y = F_right.points[y_index];
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

std::optional<std::size_t> NextUndeterminedBranchCellPair(
    const GraphData& F_left, const OrderedPartition& points_left,
    const OrderedPartition& points_right, const PartialAffineMap& A) {
  if (!points_left.HasSameShape(points_right)) return std::nullopt;

  std::optional<std::size_t> best_index;
  std::size_t best_size = 0;

  for (std::size_t cell_index = 0; cell_index < points_left.NumCells();
       ++cell_index) {
    const auto& cell = points_left.Cells()[cell_index];
    if (cell.size() <= 1) continue;

    bool has_undetermined = false;
    for (uint32_t x_index : cell) {
      if (x_index >= F_left.points.size()) return std::nullopt;
      const uint32_t x = F_left.points[x_index];
      if (!A.GetImage(x).has_value()) {
        has_undetermined = true;
        break;
      }
    }
    if (!has_undetermined) continue;

    if (!best_index.has_value() || cell.size() < best_size) {
      best_index = cell_index;
      best_size = cell.size();
    }
  }

  return best_index;
}

}  // namespace

void ResetEquivalenceSearch() {
  g_found_equivalences.clear();
  g_use_ea_validation = false;
}

void SetUseEaEquivalenceValidation(bool enabled) {
  g_use_ea_validation = enabled;
}

const std::vector<EquivalencePointMap>& GetFoundEquivalences() {
  return g_found_equivalences;
}

void CCZ_DFS_equivalence(const GraphData& F_left, const GraphData& F_right,
                         const std::vector<Hyperplane>& planes_left,
                         const std::vector<Hyperplane>& planes_right,
                         PartialAffineMap A,
                         OrderedPartition points_left,
                         OrderedPartition points_right,
                         OrderedPartition hyperplanes_left,
                         OrderedPartition hyperplanes_right,
                         std::size_t min_active_hyperplanes,
                         DfsEquivalenceGroupState group_state) {
  if (!g_found_equivalences.empty()) return;
  if (!points_left.HasSameShape(points_right)) return;

  if (A.IsFullyDetermined()) {
    (void)TryFinalizeCurrentMap(F_left, F_right, A);
    return;
  }

  if (!RefineToFixpoint(F_left, F_right, planes_left, planes_right,
                        &points_left, &points_right, &hyperplanes_left,
                        &hyperplanes_right, min_active_hyperplanes)) {
    return;
  }

  if (!ApplyMatchedSingletonCellsPair(F_left, F_right, points_left, points_right,
                                      &A)) {
    return;
  }

  const auto branch_cell =
      NextUndeterminedBranchCellPair(F_left, points_left, points_right, A);
  if (!branch_cell.has_value()) {
    (void)TryFinalizeCurrentMap(F_left, F_right, A);
    return;
  }

  const std::size_t cell_index = *branch_cell;
  const auto& left_cell = points_left.Cells()[cell_index];
  const auto& right_cell = points_right.Cells()[cell_index];
  if (left_cell.empty() || right_cell.empty()) return;

  std::optional<uint32_t> x_index_opt;
  for (uint32_t x_index_candidate : left_cell) {
    if (x_index_candidate >= F_left.points.size()) continue;
    const uint32_t x_candidate = F_left.points[x_index_candidate];
    if (A.GetImage(x_candidate).has_value()) continue;
    x_index_opt = x_index_candidate;
    break;
  }
  if (!x_index_opt.has_value()) return;
  const uint32_t x_index = *x_index_opt;
  const uint32_t x = F_left.points[x_index];

  std::unordered_set<uint32_t> tried_y_indices;
  std::vector<uint32_t> pending_candidates =
      groups::OrbitPrunedCandidatesExcludingTried(
          F_right, right_cell, A, group_state.H_generators, tried_y_indices);

  for (uint32_t y_index : pending_candidates) {
    if (y_index >= F_right.points.size()) continue;
    if (!tried_y_indices.insert(y_index).second) continue;

    const uint32_t y = F_right.points[y_index];
    PartialAffineMap next_A = A;
    if (!next_A.Update(x, y)) continue;

    OrderedPartition next_points_left = points_left;
    OrderedPartition next_points_right = points_right;
    OrderedPartition next_hyperplanes_left = hyperplanes_left;
    OrderedPartition next_hyperplanes_right = hyperplanes_right;

    if (!next_points_left.Individualize(x_index)) continue;
    if (!next_points_right.Individualize(y_index)) continue;

    DfsEquivalenceGroupState child_group_state;
    child_group_state.fixed_right_indices = group_state.fixed_right_indices;
    child_group_state.fixed_right_indices.push_back(y_index);
    if (!groups::BuildChildPointStabilizerGenerators(
            group_state.H_generators, y_index, F_right.points.size(),
            &child_group_state.H_generators)) {
      child_group_state.H_generators.clear();
    }

    CCZ_DFS_equivalence(F_left, F_right, planes_left, planes_right,
                        std::move(next_A),
                        std::move(next_points_left), std::move(next_points_right),
                        std::move(next_hyperplanes_left),
                        std::move(next_hyperplanes_right),
                        min_active_hyperplanes, std::move(child_group_state));
    if (!g_found_equivalences.empty()) return;
  }
}
