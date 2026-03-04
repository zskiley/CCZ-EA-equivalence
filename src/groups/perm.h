#ifndef GROUPS_PERM_H
#define GROUPS_PERM_H

#include <cstddef>
#include <cstdint>
#include <vector>

namespace groups {

class Permutation {
 public:
  Permutation() = default;
  // Stores permutations on degrees up to 2^16 (65536 points) using uint16_t
  // images to reduce memory footprint. For our CCZ workloads (n <= 11),
  // degree = 2^n fits easily in this range.
  explicit Permutation(std::vector<uint16_t> images);

  static Permutation Identity(std::size_t degree);

  std::size_t Degree() const;
  uint32_t Apply(uint32_t x) const;
  const std::vector<uint16_t>& Images() const;

  bool IsIdentity() const;
  bool IsValid() const;
  Permutation Inverse() const;

  bool operator==(const Permutation& other) const;

 private:
  struct UncheckedTag {};
  // Constructs without validating bijectivity. Only use when validity is
  // guaranteed by construction (e.g., composition/inversion of valid perms).
  Permutation(std::vector<uint16_t> images, UncheckedTag);

  std::vector<uint16_t> images_;

  friend Permutation Compose(const Permutation& left, const Permutation& right);
  friend Permutation Compose3(const Permutation& left, const Permutation& mid,
                              const Permutation& right);
};

struct PermutationHash {
  std::size_t operator()(const Permutation& p) const;
};

Permutation Compose(const Permutation& left, const Permutation& right);
// Computes left ∘ mid ∘ right (apply right, then mid, then left).
Permutation Compose3(const Permutation& left, const Permutation& mid,
                     const Permutation& right);

}  // namespace groups

#endif
