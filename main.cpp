#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "algorithms.h"
#include "field_basics.h"
#include "graph.h"

namespace {

template <typename Fn>
long long TimeMs(Fn&& fn) {
  const auto t0 = std::chrono::steady_clock::now();
  fn();
  const auto t1 = std::chrono::steady_clock::now();
  return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
}

GraphData BuildRandomGraph(const GF2n& field, std::uint32_t seed) {
  const std::uint32_t q = 1u << field.n;
  GraphData F;
  F.n_bits = field.n;
  F.m_bits = field.n;
  F.d_bits = 2 * field.n;
  F.points.reserve(q);
  F.is_graph.assign(1u << (2 * field.n), 0u);

  std::mt19937 rng(seed);
  std::uniform_int_distribution<std::uint32_t> dist(0u, q - 1u);

  for (std::uint32_t x = 0; x < q; ++x) {
    const std::uint32_t y = dist(rng);
    const std::uint32_t p = x | (y << field.n);
    F.points.push_back(p);
    F.is_graph[p] = 1u;
  }
  return F;
}

GraphData BuildRectangularProjectionGraph(int n_in, int m_out) {
  const std::uint32_t q = 1u << n_in;
  const std::uint32_t y_mask = (1u << m_out) - 1u;
  GraphData F;
  F.n_bits = n_in;
  F.m_bits = m_out;
  F.d_bits = n_in + m_out;
  F.points.reserve(q);
  F.is_graph.assign(1u << (n_in + m_out), 0u);

  for (std::uint32_t x = 0; x < q; ++x) {
    const std::uint32_t y = x & y_mask;  // f(x): F_2^n -> F_2^m
    const std::uint32_t p = x | (y << n_in);
    F.points.push_back(p);
    F.is_graph[p] = 1u;
  }
  return F;
}

struct FamilyCase {
  std::string label;
  int exponent = 0;
};

int ReduceExponentModOrder(int exponent, int n) {
  const int order = (1 << n) - 1;
  int e = exponent % order;
  if (e < 0) e += order;
  return (e == 0) ? order : e;
}

std::vector<FamilyCase> BuildSupportedPowerFamilyCases(int n) {
  if (n <= 0 || n >= 31) return {};
  const int order = (1 << n) - 1;
  std::vector<uint8_t> seen(static_cast<std::size_t>(order + 1), 0u);
  std::vector<FamilyCase> out;

  auto add_case = [&](const std::string& family_label, int exponent_raw) {
    const int exponent = ReduceExponentModOrder(exponent_raw, n);
    if (exponent <= 0 || exponent > order) return;
    if (seen[static_cast<std::size_t>(exponent)] != 0u) return;
    seen[static_cast<std::size_t>(exponent)] = 1u;
    out.push_back(
        FamilyCase{family_label + ",d=" + std::to_string(exponent), exponent});
  };

  for (int i = 1; i < n; ++i) {
    if (std::gcd(i, n) != 1) continue;
    const int gold_exp = (1 << i) + 1;
    add_case("Gold(n=" + std::to_string(n) + ",i=" + std::to_string(i) + ")",
             gold_exp);

    const int kasami_exp = (1 << (2 * i)) - (1 << i) + 1;
    add_case(
        "Kasami(n=" + std::to_string(n) + ",i=" + std::to_string(i) + ")",
        kasami_exp);
  }

  if ((n & 1) == 1) {
    const int t = (n - 1) / 2;
    const int welch_exp = (1 << t) + 3;
    add_case("Welch(n=" + std::to_string(n) + ")", welch_exp);

    int niho_exp = 0;
    if ((t & 1) == 0) {
      niho_exp = (1 << t) + (1 << (t / 2)) - 1;
    } else {
      niho_exp = (1 << t) + (1 << ((3 * t + 1) / 2)) - 1;
    }
    add_case("Niho(n=" + std::to_string(n) + ")", niho_exp);

    const int inverse_exp = (1 << n) - 2;
    add_case("Inverse(n=" + std::to_string(n) + ")", inverse_exp);
  }

  if ((n % 5) == 0) {
    const int i = n / 5;
    const int dobbertin_exp =
        (1 << (4 * i)) + (1 << (3 * i)) + (1 << (2 * i)) + (1 << i) - 1;
    add_case("Dobbertin(n=" + std::to_string(n) + ")", dobbertin_exp);
  }

  return out;
}

void RunAndPrint(const char* mode, const std::string& label, const GraphData& F,
                 bool ea_mode, double time_limit_seconds) {
  std::vector<GraphPointMap> autos;
  const auto run_ms = TimeMs([&]() {
    if (ea_mode) {
      autos =
          algorithms::ea_auto(F, /*timelimit_seconds=*/time_limit_seconds);
    } else {
      autos =
          algorithms::ccz_auto(F, /*timelimit_seconds=*/time_limit_seconds);
    }
  });
  std::cout << mode << " " << label << " automorphisms: " << GetTotalAutoGroup()
            << "\n";
  std::cout << mode << " " << label
            << " stabilizer size now: " << autos.size() << "\n";
  std::cout << mode << " " << label
            << (ea_mode ? " algorithms::ea_auto: "
                        : " algorithms::ccz_auto: ")
            << run_ms
            << " ms\n\n";
}

void RunEquivalenceAndPrint(const char* mode, const std::string& label_left,
                            const GraphData& F, const std::string& label_right,
                            const GraphData& G, bool ea_mode,
                            double time_limit_seconds) {
  std::optional<EquivalencePointMap> eq;
  const auto run_ms = TimeMs([&]() {
    if (ea_mode) {
      eq = algorithms::ea_equivalence(F, G, time_limit_seconds);
    } else {
      eq = algorithms::ccz_equivalence(F, G, time_limit_seconds);
    }
  });

  std::cout << mode << " equivalence " << label_left << " ~ " << label_right
            << ": " << (eq.has_value() ? "YES" : "NO");
  if (eq.has_value()) {
    std::cout << ", map size=" << eq->size();
  }
  std::cout << "\n";
  std::cout << mode << " " << label_left << " vs " << label_right
            << (ea_mode ? " algorithms::ea_equivalence: "
                        : " algorithms::ccz_equivalence: ")
            << run_ms << " ms\n\n";
}

void RunFamilyAutosForDimensions(const std::vector<int>& dimensions,
                                 double time_limit_seconds) {
  std::cout << "\nRunning supported family.txt power-map cases\n";
  for (int n : dimensions) {
    if (n <= 0 || n >= 31) continue;
    const std::uint32_t mod_poly = DefaultModPoly(n);
    GF2n field{n, mod_poly, (1u << n) - 1u};
    const std::vector<FamilyCase> cases = BuildSupportedPowerFamilyCases(n);

    std::cout << "\nn=" << n << " supported cases: " << cases.size() << "\n";
    std::cout << "\n";

    for (const FamilyCase& c : cases) {
      GraphData graph;
      const auto build_ms =
          TimeMs([&]() {
            graph = BuildGraph(BuildPowerTruthTable(field, c.exponent), n, n);
          });
      std::cout << "BuildGraph " << c.label << ": " << build_ms << " ms\n";
      RunAndPrint("CCZ", c.label, graph, false, time_limit_seconds);
      RunAndPrint("EA", c.label, graph, true, time_limit_seconds);
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  bool run_family_n12 = true;
  bool run_family_9_12 = true;
  double time_limit_seconds = 2.0;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--run_family_n12") {
      run_family_n12 = true;
      continue;
    }
    if (arg == "--run_family_9_12") {
      run_family_9_12 = true;
      continue;
    }
    if (arg == "--time_limit_seconds" && i + 1 < argc) {
      try {
        time_limit_seconds = std::stod(argv[++i]);
      } catch (...) {
        std::cerr << "invalid --time_limit_seconds value\n";
        return 2;
      }
      continue;
    }
    std::cerr
        << "usage: ccz_auto_run [--run_family_n12] [--run_family_9_12]"
        << " [--time_limit_seconds T]\n";
    return 2;
  }

  const int n = 11;
  const std::uint32_t mod_poly = DefaultModPoly(n);
  GF2n field{n, mod_poly, (1u << n) - 1u};

  GraphData graph_x3;
  const auto build_x3_ms =
      TimeMs([&]() { graph_x3 = BuildGraph(BuildPowerTruthTable(field, 3), n, n); });
  std::cout << "BuildGraph x^3: " << build_x3_ms << " ms\n";

  GraphData graph_inv;
  const auto build_inv_ms = TimeMs([&]() {
    graph_inv = BuildGraph(BuildPowerTruthTable(field, -1), n, n);
  });
  std::cout << "BuildGraph x^{-1}: " << build_inv_ms << " ms\n";

  const std::uint32_t seed = 123456789u;
  GraphData graph_rand;
  const auto build_rand_ms = TimeMs([&]() { graph_rand = BuildRandomGraph(field, seed); });
  std::cout << "BuildRandomGraph (seed=" << seed << "): " << build_rand_ms << " ms\n\n";

  GraphData graph_rect_3_2;
  const auto build_rect_ms =
      TimeMs([&]() { graph_rect_3_2 = BuildRectangularProjectionGraph(3, 2); });
  std::cout << "BuildRectangularGraph n=3,m=2 (projection): " << build_rect_ms
            << " ms\n\n";

  const int n_eq = 9;
  const std::uint32_t mod_poly_eq = DefaultModPoly(n_eq);
  GF2n field_eq{n_eq, mod_poly_eq, (1u << n_eq) - 1u};
  GraphData graph_eq_x3;
  const auto build_eq_x3_ms = TimeMs([&]() {
    graph_eq_x3 = BuildGraph(BuildPowerTruthTable(field_eq, 3), n_eq, n_eq);
  });
  GraphData graph_eq_x6;
  const auto build_eq_x6_ms = TimeMs([&]() {
    graph_eq_x6 = BuildGraph(BuildPowerTruthTable(field_eq, 6), n_eq, n_eq);
  });
  std::cout << "BuildGraph equivalence x^3 (n=9): " << build_eq_x3_ms << " ms\n";
  std::cout << "BuildGraph equivalence x^6 (n=9): " << build_eq_x6_ms << " ms\n\n";

  if (time_limit_seconds > 0.0) {
    std::cout << "time limit per auto run: " << time_limit_seconds << " s\n\n";
  }

  //RunAndPrint("CCZ", "x^3", graph_x3, false, time_limit_seconds);
  //RunAndPrint("EA", "x^3", graph_x3, true, time_limit_seconds);

  RunAndPrint("CCZ", "x^{-1}", graph_inv, false, time_limit_seconds);
  RunAndPrint("EA", "x^{-1}", graph_inv, true, time_limit_seconds);

  RunAndPrint("CCZ", "random", graph_rand, false, time_limit_seconds);
  RunAndPrint("EA", "random", graph_rand, true, time_limit_seconds);

  // Rectangular case n != m.
  RunAndPrint("CCZ", "rect_n3_m2_projection", graph_rect_3_2, false,
              time_limit_seconds);
  RunAndPrint("EA", "rect_n3_m2_projection", graph_rect_3_2, true,
              time_limit_seconds);

  RunEquivalenceAndPrint("CCZ", "x^3(n=9)", graph_eq_x3, "x^6(n=9)",
                         graph_eq_x6, false, time_limit_seconds);
  RunEquivalenceAndPrint("EA", "x^3(n=9)", graph_eq_x3, "x^6(n=9)", graph_eq_x6,
                         true, time_limit_seconds);

  if (run_family_9_12) {
    RunFamilyAutosForDimensions({9, 10, 11, 12}, time_limit_seconds);
  } else if (run_family_n12) {
    RunFamilyAutosForDimensions({12}, time_limit_seconds);
  }

  return 0;
}

