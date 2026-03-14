#include "dfs_auto.h"

#include "ambient_affine.h"
#include "partition_branch.h"
#include "groups/dfs_group_helpers.h"
#include "groups/graph_point_perm.h"
#include "groups/orbit_candidates.h"
#include "groups/schreier_sims.h"
#include "refine.h"

#include <chrono>
#include <limits>
#include <optional>
#include <unordered_set>
#include <vector>

namespace {
std::vector<GraphPointMap> g_found_autos;
uint64_t g_group_epoch = 0;
std::optional<groups::GraphPointIndex> g_graph_point_index;
std::vector<groups::Permutation> g_global_generators;
std::optional<groups::SchreierSims> g_global_group;
bool g_use_ea_validation = false;
bool g_has_time_limit = false;
std::chrono::steady_clock::time_point g_deadline;
bool g_auto_search_timed_out = false;
uint64_t g_kernel_order = 1u;

uint64_t SaturatingMulU64(uint64_t a, uint64_t b) {
  if (a == 0u || b == 0u) return 0u;
  const uint64_t limit = std::numeric_limits<uint64_t>::max();
  if (a > limit / b) return limit;
  return a * b;
}

bool IsTimedOut() {
  if (!g_has_time_limit) return false;
  if (std::chrono::steady_clock::now() < g_deadline) return false;
  g_auto_search_timed_out = true;
  return true;
}

bool RegisterLeafPermutation(const GraphPointMap& map) {
  if (!g_graph_point_index.has_value() || !g_global_group.has_value()) return false;

  groups::Permutation perm;
  if (!groups::GraphPointMapToPermutation(map, *g_graph_point_index, &perm)) {
    return false;
  }

  if (g_global_group->Contains(perm)) return false;

  g_global_generators.push_back(std::move(perm));
  g_global_generators =
      groups::DeduplicateGenerators(std::move(g_global_generators));
  g_global_group->SetGenerators(g_global_generators);
  g_global_group->Build();
  ++g_group_epoch;
  return true;
}

bool IsValidMapByMode(const PartialAffineMap& A, const GraphData& F) {
  return g_use_ea_validation ? A.valid_ea(F) : A.valid_ccz(F);
}

bool TryFinalizeCurrentMap(const GraphData& F, const PartialAffineMap& A) {
  if (!IsValidMapByMode(A, F)) return false;
  const int n_pts = static_cast<int>(F.points.size());
  GraphPointMap map;
  map.reserve(n_pts);
  for (int i = 0; i < n_pts; ++i) {
    const auto y = A.GetImage(F.points[i]);
    if (!y.has_value()) return false;
    map.push_back({F.points[i], *y});
  }
  if (!RegisterLeafPermutation(map)) return false;
  g_found_autos.push_back(std::move(map));
  return true;
}
}

void ResetAutoSearch() {
  g_found_autos.clear();
  g_group_epoch = 0;
  g_graph_point_index.reset();
  g_global_generators.clear();
  g_global_group.reset();
  g_has_time_limit = false;
  g_auto_search_timed_out = false;
  g_kernel_order = 1u;
}

void SetAutoSearchTimeLimitSeconds(double seconds) {
  if (seconds > 0.0) {
    g_has_time_limit = true;
    g_deadline = std::chrono::steady_clock::now() +
                 std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                     std::chrono::duration<double>(seconds));
    return;
  }
  g_has_time_limit = false;
}

void SetUseEaValidation(bool enabled) { g_use_ea_validation = enabled; }

bool FoundEntireAutoGroup() { return !g_auto_search_timed_out; }

uint64_t GetTotalAutoGroup() {
  if (g_global_group.has_value()) {
    return SaturatingMulU64(g_global_group->Order(), g_kernel_order);
  }
  return 0;
}

uint64_t GetGraphActionAutoGroup() {
  if (g_global_group.has_value()) return g_global_group->Order();
  return 0;
}

const std::vector<groups::Permutation>& GetAutoGroupGenerators() {
  return g_global_generators;
}

const std::vector<GraphPointMap>& GetFoundAutos() { return g_found_autos; }

void InitializeGroupSearch(const GraphData& F) {
  groups::GraphPointIndex index(F);

  g_graph_point_index = std::move(index);
  g_global_generators.clear();
  g_global_group.emplace(F.points.size());
  g_global_group->SetGenerators({});
  g_global_group->Build();
  g_group_epoch = 0;
  g_kernel_order = KernelAffineOrder(F);
}

void AddInitialGroupGenerators(std::vector<groups::Permutation> generators) {
  groups::AddInitialGroupGeneratorsImpl(
      g_graph_point_index, &g_global_group, &g_global_generators,
      &g_group_epoch, std::move(generators));
}

DfsGroupState MakeRootGroupState() {
  DfsGroupState state;
  groups::FillRootGroupStateImpl(g_global_generators, g_group_epoch, &state);
  return state;
}

void SetComputedAutoGroup(std::vector<GraphPointMap> autos,
                          std::vector<groups::Permutation> generators,
                          bool timed_out) {
  g_found_autos = std::move(autos);
  g_global_generators = groups::DeduplicateGenerators(std::move(generators));
  if (g_global_group.has_value()) {
    g_global_group->SetGenerators(g_global_generators);
    g_global_group->Build();
  }
  g_group_epoch = 0;
  g_auto_search_timed_out = timed_out;
}

void CCZ_DFS_auto(const GraphData& F, const std::vector<Hyperplane>& planes,
                  PartialAffineMap A,
                  OrderedPartition points_left,
                  OrderedPartition points_right,
                  OrderedPartition hyperplanes_left,
                  OrderedPartition hyperplanes_right,
                  std::size_t min_active_hyperplanes,
                  DfsGroupState group_state) {
  if (IsTimedOut()) return;
  if (!g_graph_point_index.has_value()) return;
  if (A.IsFullyDetermined()) {
    if (!TryFinalizeCurrentMap(F, A)) return;
    return;
  }

  {
    if (!RefineToFixpoint(F, F, planes, planes, &points_left, &points_right,
                          &hyperplanes_left, &hyperplanes_right,
                          min_active_hyperplanes)) {
      return;
    }

    if (!partition_branch::ApplyMatchedSingletonCells(F, points_left,
                                                      points_right, &A)) {
      return;
    }
    if (!IsValidMapByMode(A, F)) return;
  }

  const auto branch_cell =
      partition_branch::NextUndeterminedBranchCell(F, points_left, points_right,
                                                   A);
  if (!branch_cell.has_value()) {
    if (!TryFinalizeCurrentMap(F, A)) return;
    return;
  }

  const std::size_t cell_index = *branch_cell;
  const auto& left_cell = points_left.Cells()[cell_index];
  const auto& right_cell = points_right.Cells()[cell_index];
  if (left_cell.empty() || right_cell.empty()) return;

  std::optional<uint32_t> x_index_opt;
  for (uint32_t x_index_candidate : left_cell) {
    const uint32_t x_candidate = F.points[x_index_candidate];
    if (A.GetImage(x_candidate).has_value()) continue;
    x_index_opt = x_index_candidate;
    break;
  }
  if (!x_index_opt.has_value()) return;
  const uint32_t x_index = *x_index_opt;
  const uint32_t x = F.points[x_index];

  std::unordered_set<uint32_t> tried_y_indices;
  std::vector<uint32_t> pending_candidates =
      groups::OrbitPrunedCandidatesExcludingTried(
          F, right_cell, A, group_state.H_generators, tried_y_indices);
  std::size_t pending_pos = 0;
  while (pending_pos < pending_candidates.size()) {
    if (IsTimedOut()) return;
    const uint32_t y_index = pending_candidates[pending_pos++];
    if (y_index >= F.points.size()) continue;
    if (!tried_y_indices.insert(y_index).second) continue;

    const uint32_t y = F.points[y_index];

    PartialAffineMap next_A = A;
    if (!next_A.Update(x, y)) continue;
    if (!IsValidMapByMode(next_A, F)) continue;

    OrderedPartition next_points_left = points_left;
    OrderedPartition next_points_right = points_right;
    OrderedPartition next_hyperplanes_left = hyperplanes_left;
    OrderedPartition next_hyperplanes_right = hyperplanes_right;

    if (!next_points_left.Individualize(x_index)) continue;
    if (!next_points_right.Individualize(y_index)) continue;

    DfsGroupState child_group_state;
    child_group_state.fixed_right_indices = group_state.fixed_right_indices;
    child_group_state.H_epoch = group_state.H_epoch;
    child_group_state.fixed_right_indices.push_back(y_index);
    if (!groups::BuildChildPointStabilizerGenerators(
            group_state.H_generators, y_index, g_graph_point_index->Size(),
            &child_group_state.H_generators)) {
      child_group_state.H_generators.clear();
    }
    const uint64_t group_epoch_before = g_group_epoch;
    CCZ_DFS_auto(F, planes, std::move(next_A), std::move(next_points_left),
                 std::move(next_points_right), std::move(next_hyperplanes_left),
                 std::move(next_hyperplanes_right), min_active_hyperplanes,
                 std::move(child_group_state));

    if (g_group_epoch != group_epoch_before) {
      groups::RefreshPointwiseStabilizerFromGlobal(
          g_global_generators, group_state.fixed_right_indices,
          g_graph_point_index->Size(), g_group_epoch,
          &group_state.H_generators, &group_state.H_epoch);
      pending_candidates = groups::OrbitPrunedCandidatesExcludingTried(
          F, right_cell, A, group_state.H_generators, tried_y_indices);
      pending_pos = 0;
    }
  }
}
