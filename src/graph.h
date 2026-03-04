#ifndef GRAPH_H
#define GRAPH_H

#include <cstdint>
#include <vector>

struct GraphData {
  std::vector<uint32_t> points;
  std::vector<uint8_t> is_graph;
  int n_bits = 0;
  int m_bits = 0;
  int d_bits = 0;
};

GraphData BuildGraph(const std::vector<uint32_t>& truth_table, int n_bits,
                     int m_bits);

#endif
