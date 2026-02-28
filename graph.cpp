#include "graph.h"

GraphData BuildGraph(const std::vector<uint32_t>& truth_table, int n_bits,
                     int m_bits) {
  const uint32_t size = static_cast<uint32_t>(1u << n_bits);
  const int d_bits = n_bits + m_bits;
  const uint32_t y_mask = static_cast<uint32_t>((1u << m_bits) - 1u);
  GraphData G;
  G.n_bits = n_bits;
  G.m_bits = m_bits;
  G.d_bits = d_bits;
  G.points.reserve(size);
  G.is_graph.assign(1u << d_bits, 0);
  for (uint32_t x = 0u; x < size; ++x) {
    const uint32_t y = truth_table[x] & y_mask;
    const uint32_t pt =
        static_cast<uint32_t>(x) | (static_cast<uint32_t>(y) << n_bits);
    G.points.push_back(pt);
    G.is_graph[pt] = 1;
  }
  return G;
}
