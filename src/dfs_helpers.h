#ifndef DFS_HELPERS_H
#define DFS_HELPERS_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

inline bool SelectSmallestColorClasses(const std::vector<uint32_t>& hist,
                                       std::size_t threshold,
                                       std::vector<uint8_t>* keep,
                                       std::size_t* keep_total) {
  if (keep == nullptr || keep_total == nullptr) return false;
  keep->assign(hist.size(), 0u);
  *keep_total = 0u;
  if (hist.empty()) return false;

  if (threshold == 0) threshold = 1;

  std::vector<std::pair<uint32_t, uint32_t>> items;
  items.reserve(hist.size());
  for (std::size_t color = 0; color < hist.size(); ++color) {
    const uint32_t count = hist[color];
    if (count == 0u) continue;
    items.emplace_back(count, static_cast<uint32_t>(color));
  }
  if (items.empty()) return false;

  std::sort(items.begin(), items.end(),
            [](const auto& a, const auto& b) {
              if (a.first != b.first) return a.first < b.first;
              return a.second < b.second;
            });

  for (const auto& kv : items) {
    const uint32_t count = kv.first;
    const std::size_t color = static_cast<std::size_t>(kv.second);
    (*keep)[color] = 1u;
    *keep_total += count;
    if (*keep_total >= threshold) break;
  }

  return *keep_total > 0u;
}

#endif
