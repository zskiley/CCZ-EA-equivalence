#ifndef ORDERED_PARTITION_H
#define ORDERED_PARTITION_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

class OrderedPartition {
 public:
  using Element = uint32_t;
  using Cell = std::vector<Element>;

  // Creates an empty ordered partition.
  OrderedPartition() = default;

  // Creates a partition with one ordered cell containing all elements.
  explicit OrderedPartition(const std::vector<Element>& elements);

  // Creates a partition from explicit ordered cells.
  // Each element must appear exactly once and cells must be non-empty.
  explicit OrderedPartition(std::vector<Cell> cells);

  // Creates a partition from explicit ordered cells with optional indexing.
  // If indexing is disabled, element lookup (FindCell/Individualize) is not
  // supported but refinement/Cells() iteration still works.
  OrderedPartition(std::vector<Cell> cells, bool indexing_enabled);

  // Returns all cells in order.
  const std::vector<Cell>& Cells() const;

  // Flattens all partition elements in deterministic order:
  // by cell order, then by in-cell order.
  std::vector<Element> FlattenElementsInCellOrder() const;

  // Returns the number of cells.
  std::size_t NumCells() const;

  // Returns true iff both partitions have the same ordered cell-size pattern.
  bool HasSameShape(const OrderedPartition& other) const;

  // Returns a deterministic per-cell weight used by refinement code.
  uint64_t CellWeightForCellIndex(std::size_t cell_index) const;

  // Splits the cell containing e into [{e}, remainder].
  bool Individualize(Element e);

  // Fast refinement variant using a flat label array indexed by element id.
  // Label value `missing_label` means "no label provided".
  void RefineByFlatLabels(
      const std::vector<uint64_t>& labels,
      uint64_t missing_label = std::numeric_limits<uint64_t>::max());

  static constexpr uint64_t kMissingLabel =
      std::numeric_limits<uint64_t>::max();

 private:
  void ValidateAndIndex();

  std::vector<Cell> cells_;
  bool indexing_enabled_ = true;
  std::unordered_map<Element, std::size_t> element_to_cell_;
};

#endif
