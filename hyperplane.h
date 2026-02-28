#ifndef HYPERPLANE_H
#define HYPERPLANE_H

#include <cstddef>
#include <cstdint>
#include <vector>

struct GraphData;
class OrderedPartition;

// Affine hyperplane in F_2^d given by normal · x = c.
// The normal is a d-bit vector stored in the low bits of `normal`.
struct Hyperplane {
  uint32_t normal = 0;
  uint8_t c = 0;
};

bool BuildHyperplaneSubsetByWalsh(const GraphData& F, std::size_t threshold,
                                  std::vector<Hyperplane>* planes,
                                  OrderedPartition* hyperplanes_partition);

void BuildEAHyperplaneSubsetByWalsh(const GraphData& F, std::size_t threshold,
                                    std::vector<Hyperplane>* planes,
                                    OrderedPartition* hyperplanes_partition);

#endif
