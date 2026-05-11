#ifndef PARTIAL_MAP_H
#define PARTIAL_MAP_H

#include <cstdint>
#include <optional>
#include <cstddef>
#include <unordered_map>
#include <vector>

struct GraphData;

struct AffineMapData {
  int dimension_bits = 0;
  uint32_t translation = 0u;
  std::vector<uint32_t> linear_cols;

  uint32_t Apply(uint32_t x) const;
  bool IsIdentity() const;
};

class PartialAffineMap {
 public:
  // Creates an empty partial map on F_2^{dimension_bits}.
  explicit PartialAffineMap(int dimension_bits);

  // Update step (Algorithm 5 style): tries to add x -> y while preserving
  // affine-consistency constraints of the current partial map.
  // Returns false on conflict.
  bool Update(uint32_t x, uint32_t y);

  // Returns true if y is already used as an image by some mapped domain point.
  bool HasImage(uint32_t y) const;
  // Returns image of x if explicitly mapped or affine-determined,
  // otherwise std::nullopt.
  std::optional<uint32_t> GetImage(uint32_t x) const;
  // Exposes the current explicit assignments x -> y.
  const std::unordered_map<uint32_t, uint32_t>& Mappings() const;
  // Rank of the learned linear part from mapped differences.
  std::size_t LinearRank() const;
  // True iff current constraints determine a full-rank affine map.
  bool IsFullyDetermined() const;
  // Minimum number of additional point mappings still needed to potentially
  // determine a full-rank affine map.
  std::size_t MinAdditionalPointsForFullDetermination() const;
  // Checks whether the current partial map is still CCZ-valid on graph F.
  bool valid_ccz(const GraphData& F) const;
  // Checks whether the current partial map is still EA-valid on graph F.
  bool valid_ea(const GraphData& F) const;
  // Checks the EA block-form condition on the current learned linear span:
  // any known pure output-space direction must map back into the output space.
  bool PreservesOutputSubspace(int n_bits, int m_bits) const;
  // Returns one full affine extension consistent with the current partial map.
  bool ExtractRepresentativeAffineMap(AffineMapData* out, int n_bits = 0,
                                      int m_bits = 0) const;

 private:
  bool TryExpressDx(uint32_t dx, uint64_t* coeff_mask) const;
  bool IsInVSpan(uint32_t dy) const;
  uint32_t EvalDyFromCoeff(uint64_t coeff_mask) const;
  void InsertBasisPair(uint32_t dx, uint32_t dy);

  int dimension_bits_;
  std::unordered_map<uint32_t, uint32_t> mapping_;
  std::unordered_map<uint32_t, uint32_t> inverse_mapping_;
  std::vector<uint8_t> used_image_;

  bool has_anchor_;
  uint32_t anchor_x_;
  uint32_t anchor_y_;

  std::vector<uint32_t> basis_u_;
  std::vector<uint32_t> basis_v_;

  std::vector<uint8_t> u_has_pivot_;
  std::vector<uint32_t> u_pivot_row_;
  std::vector<uint64_t> u_pivot_expr_;

  std::vector<uint8_t> v_has_pivot_;
  std::vector<uint32_t> v_pivot_row_;
};

#endif
