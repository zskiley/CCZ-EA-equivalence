#include "ordered_partition.h"

#include <algorithm>
#include <cstdlib>

namespace {
constexpr uint32_t kCellWeightSeed = 0x9e3779b9u;
constexpr uint32_t kCellWeightMod = 8191u;

uint32_t Mix32(uint32_t x) {
  x ^= x >> 16;
  x *= 0x7feb352du;
  x ^= x >> 15;
  x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}

uint64_t CellWeight(std::size_t cell_index, std::size_t total_cells) {
  uint32_t x = static_cast<uint32_t>(cell_index);
  x ^= static_cast<uint32_t>(total_cells) * 0x85ebca6bu;
  x ^= kCellWeightSeed;
  const uint32_t mixed = Mix32(x);
  return static_cast<uint64_t>((mixed % kCellWeightMod) + 1u);
}
}  // namespace

OrderedPartition::OrderedPartition(const std::vector<Element>& elements) {
  if (!elements.empty()) {
    cells_.push_back(elements);
  }
  if (!indexing_enabled_) return;
  element_to_cell_.clear();
  if (!cells_.empty()) {
    element_to_cell_.reserve(cells_[0].size());
    for (Element e : cells_[0]) {
      element_to_cell_[e] = 0;
    }
  }
}

OrderedPartition::OrderedPartition(std::vector<Cell> cells) : cells_(std::move(cells)) {
  ValidateAndIndex();
}

OrderedPartition::OrderedPartition(std::vector<Cell> cells, bool indexing_enabled)
    : cells_(std::move(cells)), indexing_enabled_(indexing_enabled) {
  ValidateAndIndex();
}

const std::vector<OrderedPartition::Cell>& OrderedPartition::Cells() const {
  return cells_;
}

std::vector<OrderedPartition::Element>
OrderedPartition::FlattenElementsInCellOrder() const {
  std::size_t total = 0;
  for (const auto& cell : cells_) {
    total += cell.size();
  }

  std::vector<Element> out;
  out.reserve(total);
  for (const auto& cell : cells_) {
    out.insert(out.end(), cell.begin(), cell.end());
  }
  return out;
}

std::size_t OrderedPartition::NumCells() const { return cells_.size(); }

bool OrderedPartition::HasSameShape(const OrderedPartition& other) const {
  if (cells_.size() != other.cells_.size()) return false;
  for (std::size_t i = 0; i < cells_.size(); ++i) {
    if (cells_[i].size() != other.cells_[i].size()) return false;
  }
  return true;
}

uint64_t OrderedPartition::CellWeightForCellIndex(std::size_t cell_index) const {
  return CellWeight(cell_index, cells_.size());
}

bool OrderedPartition::Individualize(Element e) {
  if (!indexing_enabled_) return false;
  const auto it = element_to_cell_.find(e);
  if (it == element_to_cell_.end()) return false;
  const std::size_t cell_index = it->second;
  if (cells_[cell_index].size() == 1) return true;

  const Cell& original = cells_[cell_index];
  Cell singleton;
  Cell remainder;
  singleton.reserve(1);
  remainder.reserve(original.size() - 1);

  for (Element v : original) {
    if (v == e) {
      singleton.push_back(v);
    } else {
      remainder.push_back(v);
    }
  }
  if (singleton.empty() || remainder.empty()) return false;

  cells_[cell_index] = std::move(singleton);
  cells_.insert(cells_.begin() + static_cast<std::ptrdiff_t>(cell_index + 1),
                std::move(remainder));
  ValidateAndIndex();
  return true;
}

void OrderedPartition::RefineByFlatLabels(
    const std::vector<uint64_t>& labels, uint64_t missing_label) {
  std::vector<Cell> refined;
  refined.reserve(cells_.size());

  for (std::size_t c = 0; c < cells_.size(); ++c) {
    const auto& cell = cells_[c];
    if (cell.size() <= 1) {
      refined.push_back(cell);
      continue;
    }

    std::unordered_map<uint64_t, std::size_t> group_index;
    group_index.reserve(cell.size());
    std::vector<uint64_t> distinct_labels;
    distinct_labels.reserve(cell.size());
    std::vector<Cell> groups;
    groups.reserve(cell.size());

    for (Element e : cell) {
      if (e >= labels.size()) std::abort();
      const uint64_t label = labels[e];
      if (label == missing_label) std::abort();

      auto [it, inserted] = group_index.emplace(label, groups.size());
      if (inserted) {
        groups.emplace_back();
        distinct_labels.push_back(label);
      }
      groups[it->second].push_back(e);
    }

    std::sort(distinct_labels.begin(), distinct_labels.end());
    for (uint64_t label : distinct_labels) {
      const auto it = group_index.find(label);
      if (it == group_index.end()) std::abort();
      refined.push_back(std::move(groups[it->second]));
    }
  }

  cells_ = std::move(refined);
  ValidateAndIndex();
}

void OrderedPartition::ValidateAndIndex() {
  if (!indexing_enabled_) return;

  element_to_cell_.clear();
  for (std::size_t i = 0; i < cells_.size(); ++i) {
    for (Element e : cells_[i]) {
      element_to_cell_[e] = i;
    }
  }
}
