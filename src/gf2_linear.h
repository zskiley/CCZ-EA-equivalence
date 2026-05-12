#ifndef GF2_LINEAR_H
#define GF2_LINEAR_H

#include <cstdint>
#include <vector>

namespace gf2_linear {

int MostSignificantBit(uint32_t x);
void ReduceWithBasis(const std::vector<uint32_t>& basis, uint32_t* value);
bool AddToEchelonBasis(uint32_t value, std::vector<uint32_t>* basis);
bool IsInSpan(const std::vector<uint32_t>& basis, uint32_t value);
std::vector<uint32_t> IntersectionBasis(const std::vector<uint32_t>& left,
                                        const std::vector<uint32_t>& right);
std::vector<uint32_t> ExtendToBasis(const std::vector<uint32_t>& basis,
                                    int dimension_bits);

uint32_t ApplyLinearColumns(const std::vector<uint32_t>& cols, uint32_t x);
bool InvertLinearColumns(const std::vector<uint32_t>& cols, int dimension_bits,
                         std::vector<uint32_t>* inv_cols);
bool BuildLinearMapFromBasisPairs(const std::vector<uint32_t>& domain_basis,
                                  const std::vector<uint32_t>& image_basis,
                                  int dimension_bits,
                                  std::vector<uint32_t>* cols);

}  // namespace gf2_linear

#endif
