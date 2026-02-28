#include "refine.h"

#include "graph.h"
#include "hyperplane.h"
#include "weighted_fft.h"

#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace {
long double DirectWalshScale() {
  static const long double scale = []() {
    const char* env = std::getenv("CCZ_DIRECT_WALSH_SCALE");
    if (env == nullptr || *env == '\0') return 1.0L;
    char* end = nullptr;
    const long double parsed = std::strtold(env, &end);
    if (end == env || parsed < 0.0L) return 1.0L;
    return parsed;
  }();
  return scale;
}

bool ShouldUseDirect(std::size_t direct_est, std::size_t fwht_est) {
  const long double scale = DirectWalshScale();
  if (scale == 1.0L) return direct_est <= fwht_est;
  return static_cast<long double>(direct_est) <=
         static_cast<long double>(fwht_est) * scale;
}

// Converts signed FWHT outputs to a non-negative key used for label histograms.
uint64_t AbsSigned64ToU64(int64_t v) {
  if (v >= 0) return static_cast<uint64_t>(v);
  return static_cast<uint64_t>(-(v + 1)) + 1ull;
}

// Lightweight progress measure used to detect refinement fixpoints.
std::size_t TotalCells(const OrderedPartition& points_left,
                       const OrderedPartition& points_right,
                       const OrderedPartition& hyperplanes_left,
                       const OrderedPartition& hyperplanes_right) {
  return points_left.NumCells() + points_right.NumCells() +
         hyperplanes_left.NumCells() + hyperplanes_right.NumCells();
}

// Builds paired flat labels for any paired partition (points or hyperplanes)
// from score vectors. Also checks paired-cell histogram compatibility so both
// sides split consistently.
bool BuildPairedLabelsFromScoresDirect(
    const OrderedPartition& left, const OrderedPartition& right,
    const std::vector<int64_t>& left_scores, const std::vector<int64_t>& right_scores,
    std::vector<uint64_t>* left_labels, std::vector<uint64_t>* right_labels) {
  if (left_labels == nullptr || right_labels == nullptr) return false;
  left_labels->assign(left_scores.size(), OrderedPartition::kMissingLabel);
  right_labels->assign(right_scores.size(), OrderedPartition::kMissingLabel);

  const std::size_t n_cells = left.NumCells();

  std::unordered_map<uint64_t, uint32_t> sig_counts;
  std::vector<uint64_t> uniq;
  sig_counts.reserve(256);
  uniq.reserve(256);

  for (std::size_t cell_index = 0; cell_index < n_cells; ++cell_index) {
    const auto& lc = left.Cells()[cell_index];
    const auto& rc = right.Cells()[cell_index];
    if (lc.size() != rc.size()) return false;
    if (lc.size() <= 1) {
      // No refinement possible and RefineByFlatLabels does not require labels
      // for singleton cells.
      continue;
    }

    sig_counts.clear();
    uniq.clear();
    if (sig_counts.bucket_count() < lc.size() * 2u) {
      sig_counts.reserve(lc.size() * 2u);
    }

    for (std::size_t i = 0; i < lc.size(); ++i) {
      const uint32_t e = lc[i];
      if (e >= left_scores.size()) return false;
      const uint64_t sig = AbsSigned64ToU64(left_scores[e]);
      auto it = sig_counts.find(sig);
      if (it == sig_counts.end()) {
        sig_counts.emplace(sig, 1u);
        uniq.push_back(sig);
      } else {
        ++it->second;
      }
    }

    for (std::size_t i = 0; i < rc.size(); ++i) {
      const uint32_t e = rc[i];
      if (e >= right_scores.size()) return false;
      const uint64_t sig = AbsSigned64ToU64(right_scores[e]);
      auto it = sig_counts.find(sig);
      if (it == sig_counts.end()) return false;
      if (it->second <= 1u) {
        sig_counts.erase(it);
      } else {
        --it->second;
      }
    }
    if (!sig_counts.empty()) return false;

    std::sort(uniq.begin(), uniq.end());
    sig_counts.clear();
    if (sig_counts.bucket_count() < uniq.size() * 2u) {
      sig_counts.reserve(uniq.size() * 2u);
    }

    for (uint32_t rank = 0; rank < uniq.size(); ++rank) {
      sig_counts.emplace(uniq[rank], rank);
    }

    for (uint32_t e : lc) {
      const uint64_t sig = AbsSigned64ToU64(left_scores[e]);
      const auto it = sig_counts.find(sig);
      if (it == sig_counts.end()) return false;
      (*left_labels)[e] =
          (static_cast<uint64_t>(cell_index) << 32) ^ it->second;
    }
    for (uint32_t e : rc) {
      const uint64_t sig = AbsSigned64ToU64(right_scores[e]);
      const auto it = sig_counts.find(sig);
      if (it == sig_counts.end()) return false;
      (*right_labels)[e] =
          (static_cast<uint64_t>(cell_index) << 32) ^ it->second;
    }
  }

  return true;
}

std::size_t CountActiveElements(const OrderedPartition& partition) {
  std::size_t total = 0;
  for (std::size_t i = 0; i < partition.NumCells(); ++i) {
    total += partition.Cells()[i].size();
  }
  return total;
}

void BuildSignedNormalWeightsFromHyperplanes(
    const std::vector<Hyperplane>& planes,
    const OrderedPartition& hyperplanes, std::vector<WeightedPoint>* out) {
  out->clear();

  out->reserve(CountActiveElements(hyperplanes));
  for (std::size_t cell_id = 0; cell_id < hyperplanes.NumCells(); ++cell_id) {
    const int64_t w =
        static_cast<int64_t>(hyperplanes.CellWeightForCellIndex(cell_id));
    for (uint32_t plane_idx : hyperplanes.Cells()[cell_id]) {
      const Hyperplane& H = planes[plane_idx];
      const int64_t signed_w = (H.c == 0) ? w : -w;
      out->push_back({H.normal, signed_w});
    }
  }
}

void BuildPointWeightsFromPartition(const GraphData& F,
                                    const OrderedPartition& points,
                                    std::vector<WeightedPoint>* out) {
  out->clear();

  std::size_t reserve = 0;
  for (std::size_t i = 0; i < points.NumCells(); ++i) {
    reserve += points.Cells()[i].size();
  }
  out->reserve(reserve);
  for (std::size_t cell_id = 0; cell_id < points.NumCells(); ++cell_id) {
    const int64_t w =
        static_cast<int64_t>(points.CellWeightForCellIndex(cell_id));
    for (uint32_t point_idx : points.Cells()[cell_id]) {
      out->push_back({F.points[point_idx], w});
    }
  }
}

void BuildPointScoresFromWalsh(const GraphData& F,
                               const std::vector<int64_t>& walsh,
                               std::vector<int64_t>* scores) {
  scores->assign(F.points.size(), 0);
  for (std::size_t i = 0; i < F.points.size(); ++i) {
    const uint32_t p = F.points[i];
    (*scores)[i] = walsh[p];
  }
}

void BuildHyperplaneScoresFromWalsh(const std::vector<Hyperplane>& planes,
                                    const std::vector<int64_t>& walsh,
                                    std::vector<int64_t>* scores) {
  scores->assign(planes.size(), 0);
  for (std::size_t i = 0; i < planes.size(); ++i) {
    const Hyperplane& H = planes[i];
    const int64_t t = walsh[H.normal];
    (*scores)[i] = (H.c == 0) ? t : -t;
  }
}

// Refines point partitions from active hyperplane partitions via Walsh scores.
bool RefinePointPairFromHyperplanePair(
    const GraphData& F_left, const GraphData& F_right,
    const std::vector<Hyperplane>& planes_left,
    const std::vector<Hyperplane>& planes_right,
    const OrderedPartition& hyperplanes_left,
    const OrderedPartition& hyperplanes_right, OrderedPartition& points_left,
    OrderedPartition& points_right) {

  const int d = F_left.d_bits;
  if (F_right.d_bits != d) return false;
  if (d <= 0 || d >= 31) return false;
 
  std::vector<WeightedPoint> normals_left;
  std::vector<WeightedPoint> normals_right;
  BuildSignedNormalWeightsFromHyperplanes(planes_left, hyperplanes_left,
                                          &normals_left);
  BuildSignedNormalWeightsFromHyperplanes(planes_right, hyperplanes_right,
                                          &normals_right);

  std::vector<int64_t> scores_left;
  std::vector<int64_t> scores_right;
  static thread_local std::vector<int64_t> walsh_buf;

  {
    const std::size_t active_hyperplanes = CountActiveElements(hyperplanes_left);
    const std::size_t point_count = F_left.points.size();
    const std::size_t N = static_cast<std::size_t>(1u) << d;
    const std::size_t fwht_est = N * static_cast<std::size_t>(d);
    const std::size_t direct_est = point_count * active_hyperplanes;
    const bool use_direct = ShouldUseDirect(direct_est, fwht_est);

    if (use_direct) {
      std::vector<int64_t> walsh_left_at_points;
      std::vector<int64_t> walsh_right_at_points;
      WeightedWalshAt(d, normals_left, F_left.points, &walsh_left_at_points);
      WeightedWalshAt(d, normals_right, F_right.points, &walsh_right_at_points);
      scores_left = std::move(walsh_left_at_points);
      scores_right = std::move(walsh_right_at_points);
    } else {
      WeightedFWHT(d, normals_left, &walsh_buf);
      BuildPointScoresFromWalsh(F_left, walsh_buf, &scores_left);
      WeightedFWHT(d, normals_right, &walsh_buf);
      BuildPointScoresFromWalsh(F_right, walsh_buf, &scores_right);
    }
  }

  std::vector<uint64_t> labels_left;
  std::vector<uint64_t> labels_right;

  if (!BuildPairedLabelsFromScoresDirect(points_left, points_right, scores_left,
                                         scores_right, &labels_left,
                                         &labels_right)) {
    return false;
  }
  points_left.RefineByFlatLabels(labels_left);
  points_right.RefineByFlatLabels(labels_right);

  return true;
}

// Refines hyperplane partitions from active point partitions via Walsh scores.
bool RefineHyperplanePairFromPointPair(
    const GraphData& F_left, const GraphData& F_right,
    const std::vector<Hyperplane>& planes_left,
    const std::vector<Hyperplane>& planes_right,
    const OrderedPartition& points_left, const OrderedPartition& points_right,
    OrderedPartition& hyperplanes_left, OrderedPartition& hyperplanes_right) {
  const int d = F_left.d_bits;
  if (F_right.d_bits != d) return false;
  if (d <= 0 || d >= 31) return false;

  std::vector<WeightedPoint> weights_left;
  std::vector<WeightedPoint> weights_right;
  BuildPointWeightsFromPartition(F_left, points_left, &weights_left);
  BuildPointWeightsFromPartition(F_right, points_right, &weights_right);

  std::vector<int64_t> scores_left;
  std::vector<int64_t> scores_right;
  static thread_local std::vector<int64_t> walsh_buf;

  {
    const std::size_t active_hyperplanes = CountActiveElements(hyperplanes_left);
    const std::size_t point_count = F_left.points.size();
    const std::size_t N = static_cast<std::size_t>(1u) << d;
    const std::size_t fwht_est = N * static_cast<std::size_t>(d);
    const std::size_t direct_est = point_count * active_hyperplanes;
    const bool use_direct = ShouldUseDirect(direct_est, fwht_est);

    if (use_direct) {
      std::vector<uint32_t> plane_indices_left =
          hyperplanes_left.FlattenElementsInCellOrder();
      std::vector<uint32_t> plane_indices_right =
          hyperplanes_right.FlattenElementsInCellOrder();

      std::vector<uint32_t> queries_left;
      std::vector<uint32_t> queries_right;
      queries_left.reserve(plane_indices_left.size());
      queries_right.reserve(plane_indices_right.size());
      for (uint32_t idx : plane_indices_left) {
        if (idx >= planes_left.size()) return false;
        queries_left.push_back(planes_left[idx].normal);
      }
      for (uint32_t idx : plane_indices_right) {
        if (idx >= planes_right.size()) return false;
        queries_right.push_back(planes_right[idx].normal);
      }

      std::vector<int64_t> walsh_left_at_normals;
      std::vector<int64_t> walsh_right_at_normals;
      WeightedWalshAt(d, weights_left, queries_left, &walsh_left_at_normals);
      WeightedWalshAt(d, weights_right, queries_right, &walsh_right_at_normals);

      scores_left.assign(planes_left.size(), 0);
      scores_right.assign(planes_right.size(), 0);

      for (std::size_t j = 0; j < plane_indices_left.size(); ++j) {
        const uint32_t plane_idx = plane_indices_left[j];
        const Hyperplane& H = planes_left[plane_idx];
        const int64_t t = walsh_left_at_normals[j];
        scores_left[plane_idx] = (H.c == 0) ? t : -t;
      }
      for (std::size_t j = 0; j < plane_indices_right.size(); ++j) {
        const uint32_t plane_idx = plane_indices_right[j];
        const Hyperplane& H = planes_right[plane_idx];
        const int64_t t = walsh_right_at_normals[j];
        scores_right[plane_idx] = (H.c == 0) ? t : -t;
      }
    } else {
      WeightedFWHT(d, weights_left, &walsh_buf);
      BuildHyperplaneScoresFromWalsh(planes_left, walsh_buf, &scores_left);
      WeightedFWHT(d, weights_right, &walsh_buf);
      BuildHyperplaneScoresFromWalsh(planes_right, walsh_buf, &scores_right);
    }
  }

  std::vector<uint64_t> labels_left;
  std::vector<uint64_t> labels_right;

  if (!BuildPairedLabelsFromScoresDirect(hyperplanes_left, hyperplanes_right,
                                         scores_left, scores_right, &labels_left,
                                         &labels_right)) {
    return false;
  }
  hyperplanes_left.RefineByFlatLabels(labels_left);
  hyperplanes_right.RefineByFlatLabels(labels_right);

  return true;
}

// If total hyperplane elements exceed the budget, keep only the smallest cells
// up to that budget and compact both paired partitions to those kept cells.
// Otherwise leaves partitions unchanged.
void KeepSmallestCells(std::size_t min_active_elements, OrderedPartition* left,
                       OrderedPartition* right) {
  using CellInfo = std::pair<std::size_t, std::size_t>;  // (size, index)

  std::vector<CellInfo> cells;
  cells.reserve(left->NumCells());
  std::size_t total_active_elements = 0;
  for (std::size_t i = 0; i < left->NumCells(); ++i) {
    const std::size_t sz = left->Cells()[i].size();
    if (sz == 0) continue;
    cells.push_back({sz, i});
    total_active_elements += sz;
  }

  if (total_active_elements <= min_active_elements) return;

  std::sort(cells.begin(), cells.end());
  std::vector<std::size_t> kept_indices;
  kept_indices.reserve(cells.size());
  std::size_t kept_elements = 0;
  for (const auto& cell : cells) {
    if (kept_elements >= min_active_elements) break;
    kept_indices.push_back(cell.second);
    kept_elements += cell.first;
  }
  std::sort(kept_indices.begin(), kept_indices.end());

  std::vector<OrderedPartition::Cell> left_kept_cells;
  std::vector<OrderedPartition::Cell> right_kept_cells;
  left_kept_cells.reserve(kept_indices.size());
  right_kept_cells.reserve(kept_indices.size());
  for (std::size_t i : kept_indices) {
    left_kept_cells.push_back(left->Cells()[i]);
    right_kept_cells.push_back(right->Cells()[i]);
  }

  *left = OrderedPartition(std::move(left_kept_cells), false);
  *right = OrderedPartition(std::move(right_kept_cells), false);
}
}  // namespace

// Runs alternating point/hyperplane refinement until no further cell splits occur.
// Returns false on any shape mismatch or failed paired-label consistency check.
bool RefineToFixpoint(const GraphData& F_left, const GraphData& F_right,
                      const std::vector<Hyperplane>& planes_left,
                      const std::vector<Hyperplane>& planes_right,
                      OrderedPartition* points_left,
                      OrderedPartition* points_right,
                      OrderedPartition* hyperplanes_left,
                      OrderedPartition* hyperplanes_right,
                      std::size_t min_active_hyperplanes) {
  if (F_left.d_bits != F_right.d_bits) return false;
  if (F_left.points.size() != F_right.points.size()) return false;
  if (!points_left->HasSameShape(*points_right)) return false;
  if (!hyperplanes_left->HasSameShape(*hyperplanes_right)) return false;
  if (min_active_hyperplanes == 0) {
    min_active_hyperplanes = static_cast<std::size_t>(F_left.d_bits);
  }

  while (true) {
    const std::size_t before =
        TotalCells(*points_left, *points_right, *hyperplanes_left,
                   *hyperplanes_right);

    if (!RefinePointPairFromHyperplanePair(
            F_left, F_right, planes_left, planes_right, *hyperplanes_left,
            *hyperplanes_right, *points_left,
            *points_right)) {
      return false;
    }

    if (!RefineHyperplanePairFromPointPair(
            F_left, F_right, planes_left, planes_right, *points_left,
            *points_right, *hyperplanes_left,
            *hyperplanes_right)) {
      return false;
    }

    const std::size_t after =
        TotalCells(*points_left, *points_right, *hyperplanes_left,
                   *hyperplanes_right);
    if (after == before) break;
  }

  KeepSmallestCells(min_active_hyperplanes, hyperplanes_left,
                    hyperplanes_right);

  return true;
}

