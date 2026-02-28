#include "partial_map.h"

#include "graph.h"

#include <cstddef>

namespace {
int MostSignificantBit(uint32_t x) {
  for (int b = 31; b >= 0; --b) {
    if ((x >> b) & 1u) return b;
  }
  return -1;
}
}  // namespace

PartialAffineMap::PartialAffineMap(int dimension_bits)
    : dimension_bits_(dimension_bits),
      used_image_(1u << dimension_bits, 0),
      has_anchor_(false),
      anchor_x_(0u),
      anchor_y_(0u),
      u_has_pivot_(dimension_bits, 0),
      u_pivot_row_(dimension_bits, 0u),
      u_pivot_expr_(dimension_bits, 0ull),
      v_has_pivot_(dimension_bits, 0),
      v_pivot_row_(dimension_bits, 0u) {}

bool PartialAffineMap::TryExpressDx(uint32_t dx, uint64_t* coeff_mask) const {
  uint32_t r = dx;
  uint64_t m = 0ull;
  for (int b = dimension_bits_ - 1; b >= 0; --b) {
    if (((r >> b) & 1u) == 0u) continue;
    if (!u_has_pivot_[b]) return false;
    r ^= u_pivot_row_[b];
    m ^= u_pivot_expr_[b];
  }
  if (coeff_mask != nullptr) *coeff_mask = m;
  return true;
}

bool PartialAffineMap::IsInVSpan(uint32_t dy) const {
  uint32_t r = dy;
  for (int b = dimension_bits_ - 1; b >= 0; --b) {
    if (((r >> b) & 1u) == 0u) continue;
    if (!v_has_pivot_[b]) return false;
    r ^= v_pivot_row_[b];
  }
  return true;
}

uint32_t PartialAffineMap::EvalDyFromCoeff(uint64_t coeff_mask) const {
  uint32_t out = 0u;
  const size_t m = basis_v_.size();
  for (size_t i = 0; i < m; ++i) {
    if ((coeff_mask >> i) & 1ull) out ^= basis_v_[i];
  }
  return out;
}

void PartialAffineMap::InsertBasisPair(uint32_t dx, uint32_t dy) {
  const size_t idx = basis_u_.size();
  basis_u_.push_back(dx);
  basis_v_.push_back(dy);

  uint32_t u_row = dx;
  uint64_t u_expr = (1ull << idx);
  for (int b = dimension_bits_ - 1; b >= 0; --b) {
    if (((u_row >> b) & 1u) == 0u) continue;
    if (!u_has_pivot_[b]) continue;
    u_row ^= u_pivot_row_[b];
    u_expr ^= u_pivot_expr_[b];
  }

  const int u_pivot = MostSignificantBit(u_row);
  if (u_pivot >= 0) {
    u_has_pivot_[u_pivot] = 1;
    u_pivot_row_[u_pivot] = u_row;
    u_pivot_expr_[u_pivot] = u_expr;
    for (int b = 0; b < dimension_bits_; ++b) {
      if (b == u_pivot || !u_has_pivot_[b]) continue;
      if (((u_pivot_row_[b] >> u_pivot) & 1u) == 0u) continue;
      u_pivot_row_[b] ^= u_row;
      u_pivot_expr_[b] ^= u_expr;
    }
  }

  uint32_t v_row = dy;
  for (int b = dimension_bits_ - 1; b >= 0; --b) {
    if (((v_row >> b) & 1u) == 0u) continue;
    if (!v_has_pivot_[b]) continue;
    v_row ^= v_pivot_row_[b];
  }

  const int v_pivot = MostSignificantBit(v_row);
  if (v_pivot >= 0) {
    v_has_pivot_[v_pivot] = 1;
    v_pivot_row_[v_pivot] = v_row;
    for (int b = 0; b < dimension_bits_; ++b) {
      if (b == v_pivot || !v_has_pivot_[b]) continue;
      if (((v_pivot_row_[b] >> v_pivot) & 1u) == 0u) continue;
      v_pivot_row_[b] ^= v_row;
    }
  }
}

bool PartialAffineMap::Update(uint32_t x, uint32_t y) {
  const auto it = mapping_.find(x);
  if (it != mapping_.end()) return it->second == y;

  const auto inv = inverse_mapping_.find(y);
  if (inv != inverse_mapping_.end() && inv->second != x) return false;
  if (used_image_[y]) return false;

  if (!has_anchor_) {
    has_anchor_ = true;
    anchor_x_ = x;
    anchor_y_ = y;
    mapping_[x] = y;
    inverse_mapping_[y] = x;
    used_image_[y] = 1;
    return true;
  }

  const uint32_t dx = x ^ anchor_x_;
  const uint32_t dy = y ^ anchor_y_;

  uint64_t coeff_mask = 0ull;
  if (TryExpressDx(dx, &coeff_mask)) {
    if (EvalDyFromCoeff(coeff_mask) != dy) return false;
  } else {
    if (IsInVSpan(dy)) return false;
    if (basis_u_.size() >= 63) return false;
    InsertBasisPair(dx, dy);
  }

  mapping_[x] = y;
  inverse_mapping_[y] = x;
  used_image_[y] = 1;
  return true;
}


bool PartialAffineMap::HasImage(uint32_t y) const { return used_image_[y] != 0; }

std::optional<uint32_t> PartialAffineMap::GetImage(uint32_t x) const {
  const auto it = mapping_.find(x);
  if (it != mapping_.end()) return it->second;
  if (!has_anchor_) return std::nullopt;

  uint64_t coeff_mask = 0ull;
  if (!TryExpressDx(x ^ anchor_x_, &coeff_mask)) return std::nullopt;
  return anchor_y_ ^ EvalDyFromCoeff(coeff_mask);
}

const std::unordered_map<uint32_t, uint32_t>& PartialAffineMap::Mappings() const {
  return mapping_;
}

std::size_t PartialAffineMap::LinearRank() const { return basis_u_.size(); }

bool PartialAffineMap::IsFullyDetermined() const {
  return has_anchor_ &&
         basis_u_.size() >= static_cast<std::size_t>(dimension_bits_);
}

std::size_t PartialAffineMap::MinAdditionalPointsForFullDetermination() const {
  if (!has_anchor_) return static_cast<std::size_t>(dimension_bits_ + 1);
  const std::size_t rank = basis_u_.size();
  const std::size_t d = static_cast<std::size_t>(dimension_bits_);
  if (rank >= d) return 0;
  return d - rank;
}

bool PartialAffineMap::valid_ccz(const GraphData& F) const {
  if (F.d_bits <= 0 || F.d_bits >= 31) return false;
  if (F.is_graph.size() != (1u << F.d_bits)) return false;
  if (F.d_bits != dimension_bits_) return false;

  if (mapping_.empty()) return true;
  if (!has_anchor_) return false;

  for (const auto& kv : mapping_) {
    uint64_t coeff_mask = 0ull;
    if (!TryExpressDx(kv.first ^ anchor_x_, &coeff_mask)) return false;
    const uint32_t y = anchor_y_ ^ EvalDyFromCoeff(coeff_mask);
    if (y != kv.second) return false;
    if (y >= F.is_graph.size()) return false;
    if (!F.is_graph[y]) return false;
  }

  for (uint32_t x : F.points) {
    uint64_t coeff_mask = 0ull;
    if (!TryExpressDx(x ^ anchor_x_, &coeff_mask)) continue;
    const uint32_t y = anchor_y_ ^ EvalDyFromCoeff(coeff_mask);
    if (y >= F.is_graph.size()) return false;
    if (!F.is_graph[y]) return false;
  }

  return true;
}

bool PartialAffineMap::valid_ea(const GraphData& F) const {
  if (!valid_ccz(F)) return false;
  const int n = F.n_bits;
  const int m = F.m_bits;
  if (n <= 0 || m <= 0 || n >= 31 || m >= 31) return false;
  if (n + m != F.d_bits) return false;
  const uint32_t mask_x = (1u << n) - 1u;

  if (!mapping_.empty()) {
    std::vector<int32_t> x_left_image(static_cast<std::size_t>(1u << n), -1);
    std::vector<int32_t> x_image_owner(static_cast<std::size_t>(1u << n), -1);

    for (const auto& kv : mapping_) {
      const uint32_t x = kv.first & mask_x;
      const uint32_t xp = kv.second & mask_x;

      const int32_t left = x_left_image[static_cast<std::size_t>(x)];
      if (left >= 0 && static_cast<uint32_t>(left) != xp) return false;
      x_left_image[static_cast<std::size_t>(x)] = static_cast<int32_t>(xp);

      const int32_t owner = x_image_owner[static_cast<std::size_t>(xp)];
      if (owner >= 0 && static_cast<uint32_t>(owner) != x) return false;
      x_image_owner[static_cast<std::size_t>(xp)] = static_cast<int32_t>(x);
    }
  }

  for (int i = 0; i < m; ++i) {
    const uint32_t u = (1u << (n + i));
    uint64_t coeff_mask = 0ull;
    if (!TryExpressDx(u, &coeff_mask)) continue;
    const uint32_t v = EvalDyFromCoeff(coeff_mask);
    if ((v & mask_x) != 0u) return false;
  }

  return true;
}
