#include "perm.h"

#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace groups {
namespace {

bool ValidateImages(const std::vector<uint16_t>& images) {
  const std::size_t n = images.size();
  std::vector<uint8_t> seen(n, 0u);
  for (uint16_t y : images) {
    if (y >= n) return false;
    if (seen[y] != 0u) return false;
    seen[y] = 1u;
  }
  return true;
}

std::size_t Mix64(std::size_t x) {
  x ^= (x >> 30);
  x *= static_cast<std::size_t>(0xbf58476d1ce4e5b9ull);
  x ^= (x >> 27);
  x *= static_cast<std::size_t>(0x94d049bb133111ebull);
  x ^= (x >> 31);
  return x;
}

}  // namespace

Permutation::Permutation(std::vector<uint16_t> images) : images_(std::move(images)) {
  if (images_.size() > static_cast<std::size_t>(
                          std::numeric_limits<uint16_t>::max()) +
                          1u) {
    throw std::invalid_argument(
        "Permutation: degree must be <= 65536 for uint16 storage");
  }
  if (!ValidateImages(images_)) {
    throw std::invalid_argument(
        "Permutation: images must be a bijection of [0..n-1]");
  }
}

Permutation::Permutation(std::vector<uint16_t> images, UncheckedTag)
    : images_(std::move(images)) {
}

Permutation Permutation::Identity(std::size_t degree) {
  if (degree >
      static_cast<std::size_t>(std::numeric_limits<uint16_t>::max()) + 1u) {
    throw std::invalid_argument(
        "Permutation::Identity: degree must be <= 65536 for uint16 storage");
  }
  std::vector<uint16_t> images(degree, 0u);
  for (std::size_t i = 0; i < degree; ++i) {
    images[i] = static_cast<uint16_t>(i);
  }
  return Permutation(std::move(images), UncheckedTag{});
}

std::size_t Permutation::Degree() const {
  return images_.size();
}

uint32_t Permutation::Apply(uint32_t x) const {
  if (x >= images_.size()) {
    throw std::out_of_range("Permutation::Apply: point out of range");
  }
  return static_cast<uint32_t>(images_[x]);
}

const std::vector<uint16_t>& Permutation::Images() const {
  return images_;
}

bool Permutation::IsIdentity() const {
  const std::size_t n = images_.size();
  if (n == 0u) return true;
  if (images_[0] != 0u) return false;
  if (n > 1u && images_[1] != 1u) return false;
  if (n > 2u && images_[2] != 2u) return false;
  for (std::size_t i = 3; i < n; ++i) {
    if (images_[i] != i) return false;
  }
  return true;
}

bool Permutation::IsValid() const {
  return ValidateImages(images_);
}

Permutation Permutation::Inverse() const {
  const std::size_t n = images_.size();
  std::vector<uint16_t> inverse(n, 0u);
  for (std::size_t i = 0; i < n; ++i) {
    inverse[images_[i]] = static_cast<uint16_t>(i);
  }
  return Permutation(std::move(inverse), UncheckedTag{});
}

bool Permutation::operator==(const Permutation& other) const {
  return images_ == other.images_;
}

std::size_t PermutationHash::operator()(const Permutation& p) const {
  // FNV-1a over 16-bit images, with a final avalanche mix.
  uint64_t h = 1469598103934665603ull;
  const auto& imgs = p.Images();
  const std::size_t n = imgs.size();
  const uint16_t* data = imgs.data();
  std::size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    const uint64_t word = static_cast<uint64_t>(data[i]) |
                          (static_cast<uint64_t>(data[i + 1]) << 16) |
                          (static_cast<uint64_t>(data[i + 2]) << 32) |
                          (static_cast<uint64_t>(data[i + 3]) << 48);
    h ^= word;
    h *= 1099511628211ull;
  }
  for (; i < n; ++i) {
    h ^= static_cast<uint64_t>(data[i]);
    h *= 1099511628211ull;
  }
  h ^= static_cast<uint64_t>(p.Degree());
  h *= 1099511628211ull;
  const std::size_t mixed = Mix64(static_cast<std::size_t>(h));
  return mixed;
}

Permutation Compose(const Permutation& left, const Permutation& right) {
  if (left.Degree() != right.Degree()) {
    throw std::invalid_argument("Compose: permutations must have same degree");
  }

  const std::size_t n = left.Degree();
  std::vector<uint16_t> images(n, 0u);
  for (std::size_t x = 0; x < n; ++x) {
    const uint16_t r = right.Images()[x];
    images[x] = left.Images()[static_cast<std::size_t>(r)];
  }
  return Permutation(std::move(images), Permutation::UncheckedTag{});
}

Permutation Compose3(const Permutation& left, const Permutation& mid,
                     const Permutation& right) {
  if (left.Degree() != mid.Degree() || left.Degree() != right.Degree()) {
    throw std::invalid_argument("Compose3: permutations must have same degree");
  }

  const std::size_t n = left.Degree();
  std::vector<uint16_t> images(n, 0u);
  const auto& l = left.Images();
  const auto& m = mid.Images();
  const auto& r = right.Images();
  for (std::size_t x = 0; x < n; ++x) {
    const uint16_t rx = r[x];
    const uint16_t mx = m[static_cast<std::size_t>(rx)];
    images[x] = l[static_cast<std::size_t>(mx)];
  }
  return Permutation(std::move(images), Permutation::UncheckedTag{});
}

}  // namespace groups
