#include "gf2_linear.h"

#include <cstddef>
#include <utility>

namespace gf2_linear {

int MostSignificantBit(uint32_t x) {
  for (int b = 31; b >= 0; --b) {
    if (((x >> b) & 1u) != 0u) return b;
  }
  return -1;
}

void ReduceWithBasis(const std::vector<uint32_t>& basis, uint32_t* value) {
  if (value == nullptr) return;
  for (uint32_t pivot : basis) {
    const int pivot_bit = MostSignificantBit(pivot);
    if (pivot_bit >= 0 && (((*value >> pivot_bit) & 1u) != 0u)) {
      *value ^= pivot;
    }
  }
}

bool AddToEchelonBasis(uint32_t value, std::vector<uint32_t>* basis) {
  if (basis == nullptr) return false;
  ReduceWithBasis(*basis, &value);
  if (value == 0u) return false;
  const int value_bit = MostSignificantBit(value);
  std::size_t insert_at = 0;
  while (insert_at < basis->size() &&
         MostSignificantBit((*basis)[insert_at]) > value_bit) {
    ++insert_at;
  }
  basis->insert(basis->begin() + static_cast<std::ptrdiff_t>(insert_at), value);
  for (std::size_t i = 0; i < basis->size(); ++i) {
    if (i == insert_at) continue;
    if ((((*basis)[i] >> value_bit) & 1u) != 0u) (*basis)[i] ^= value;
  }
  return true;
}

bool IsInSpan(const std::vector<uint32_t>& basis, uint32_t value) {
  ReduceWithBasis(basis, &value);
  return value == 0u;
}

std::vector<uint32_t> ExtendToBasis(const std::vector<uint32_t>& basis,
                                    int dimension_bits) {
  std::vector<uint32_t> full_basis = basis;
  std::vector<uint32_t> span = basis;
  for (int bit = 0; bit < dimension_bits; ++bit) {
    const uint32_t e = (1u << bit);
    if (IsInSpan(span, e)) continue;
    full_basis.push_back(e);
    (void)AddToEchelonBasis(e, &span);
  }
  return full_basis;
}

uint32_t ApplyLinearColumns(const std::vector<uint32_t>& cols, uint32_t x) {
  uint32_t y = 0u;
  std::size_t bit = 0;
  while (x != 0u) {
    if ((x & 1u) != 0u) y ^= cols[bit];
    x >>= 1u;
    ++bit;
  }
  return y;
}

bool InvertLinearColumns(const std::vector<uint32_t>& cols, int dimension_bits,
                         std::vector<uint32_t>* inv_cols) {
  if (inv_cols == nullptr) return false;
  if (cols.size() != static_cast<std::size_t>(dimension_bits)) return false;

  std::vector<uint32_t> rows(static_cast<std::size_t>(dimension_bits), 0u);
  std::vector<uint32_t> inv_rows(static_cast<std::size_t>(dimension_bits), 0u);
  for (int row = 0; row < dimension_bits; ++row) {
    uint32_t row_bits = 0u;
    for (int col = 0; col < dimension_bits; ++col) {
      if (((cols[static_cast<std::size_t>(col)] >> row) & 1u) != 0u) {
        row_bits |= (1u << col);
      }
    }
    rows[static_cast<std::size_t>(row)] = row_bits;
    inv_rows[static_cast<std::size_t>(row)] = (1u << row);
  }

  for (int col = 0; col < dimension_bits; ++col) {
    int pivot_row = col;
    while (pivot_row < dimension_bits &&
           ((rows[static_cast<std::size_t>(pivot_row)] >> col) & 1u) == 0u) {
      ++pivot_row;
    }
    if (pivot_row == dimension_bits) return false;
    if (pivot_row != col) {
      std::swap(rows[static_cast<std::size_t>(pivot_row)],
                rows[static_cast<std::size_t>(col)]);
      std::swap(inv_rows[static_cast<std::size_t>(pivot_row)],
                inv_rows[static_cast<std::size_t>(col)]);
    }
    for (int row = 0; row < dimension_bits; ++row) {
      if (row == col) continue;
      if (((rows[static_cast<std::size_t>(row)] >> col) & 1u) == 0u) continue;
      rows[static_cast<std::size_t>(row)] ^= rows[static_cast<std::size_t>(col)];
      inv_rows[static_cast<std::size_t>(row)] ^=
          inv_rows[static_cast<std::size_t>(col)];
    }
  }

  inv_cols->assign(static_cast<std::size_t>(dimension_bits), 0u);
  for (int col = 0; col < dimension_bits; ++col) {
    uint32_t image = 0u;
    for (int row = 0; row < dimension_bits; ++row) {
      if (((inv_rows[static_cast<std::size_t>(row)] >> col) & 1u) != 0u) {
        image |= (1u << row);
      }
    }
    (*inv_cols)[static_cast<std::size_t>(col)] = image;
  }
  return true;
}

bool BuildLinearMapFromBasisPairs(const std::vector<uint32_t>& domain_basis,
                                  const std::vector<uint32_t>& image_basis,
                                  int dimension_bits,
                                  std::vector<uint32_t>* cols) {
  if (cols == nullptr) return false;
  if (domain_basis.size() != static_cast<std::size_t>(dimension_bits)) {
    return false;
  }
  if (image_basis.size() != static_cast<std::size_t>(dimension_bits)) {
    return false;
  }
  std::vector<uint32_t> inv_domain_basis;
  if (!InvertLinearColumns(domain_basis, dimension_bits, &inv_domain_basis)) {
    return false;
  }
  cols->assign(static_cast<std::size_t>(dimension_bits), 0u);
  for (int bit = 0; bit < dimension_bits; ++bit) {
    const uint32_t coeffs = ApplyLinearColumns(inv_domain_basis, (1u << bit));
    (*cols)[static_cast<std::size_t>(bit)] =
        ApplyLinearColumns(image_basis, coeffs);
  }
  return true;
}

}  // namespace gf2_linear
