#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "schreier_sims.h"

namespace {

groups::Permutation MakeIdentity(std::size_t n) {
  return groups::Permutation::Identity(n);
}

groups::Permutation MakeCycle(std::size_t n) {
  if (n == 0u) throw std::invalid_argument("MakeCycle: n must be > 0");
  std::vector<uint16_t> images(n, 0u);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    images[i] = static_cast<uint16_t>(i + 1);
  }
  images[n - 1] = 0u;
  return groups::Permutation(std::move(images));
}

groups::Permutation MakeTransposition(std::size_t n, uint32_t a, uint32_t b) {
  if (n == 0u) throw std::invalid_argument("MakeTransposition: n must be > 0");
  if (a >= n || b >= n) {
    throw std::invalid_argument("MakeTransposition: points out of range");
  }
  std::vector<uint16_t> images(n, 0u);
  for (std::size_t i = 0; i < n; ++i) {
    images[i] = static_cast<uint16_t>(i);
  }
  const uint32_t tmp = images[a];
  images[a] = images[b];
  images[b] = tmp;
  return groups::Permutation(std::move(images));
}

std::vector<groups::Permutation> SymmetricGenerators(std::size_t n) {
  if (n < 2u) {
    return {MakeIdentity(n)};
  }
  return {MakeCycle(n), MakeTransposition(n, 0u, 1u)};
}

uint64_t FactorialU64(std::size_t n) {
  uint64_t out = 1u;
  for (std::size_t i = 2; i <= n; ++i) {
    const uint64_t prev = out;
    out *= static_cast<uint64_t>(i);
    if (out / static_cast<uint64_t>(i) != prev) {
      return 0u;
    }
  }
  return out;
}

groups::Permutation RandomWord(const std::vector<groups::Permutation>& gens,
                               const std::vector<groups::Permutation>& gen_inverses,
                               int word_len, std::mt19937& rng) {
  if (gens.empty()) throw std::invalid_argument("RandomWord: no generators");
  groups::Permutation out = groups::Permutation::Identity(gens[0].Degree());
  std::uniform_int_distribution<std::size_t> pick_gen(0u, gens.size() - 1u);
  std::uniform_int_distribution<int> pick_sign(0, 1);
  for (int step = 0; step < word_len; ++step) {
    const std::size_t idx = pick_gen(rng);
    if (pick_sign(rng) == 0) {
      out = groups::Compose(gens[idx], out);
    } else {
      out = groups::Compose(gen_inverses[idx], out);
    }
  }
  return out;
}

struct Config {
  int min_n = 8;
  int max_n = 12;
  int contains_trials = 2000;
  int word_len = 96;
  uint32_t seed = 123456789u;
};

bool ParseIntArg(const char* s, int* out) {
  if (s == nullptr || out == nullptr) return false;
  char* end = nullptr;
  const long v = std::strtol(s, &end, 10);
  if (end == s || *end != '\0') return false;
  if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max()) {
    return false;
  }
  *out = static_cast<int>(v);
  return true;
}

bool ParseU32Arg(const char* s, uint32_t* out) {
  if (s == nullptr || out == nullptr) return false;
  char* end = nullptr;
  const unsigned long v = std::strtoul(s, &end, 10);
  if (end == s || *end != '\0') return false;
  if (v > std::numeric_limits<uint32_t>::max()) return false;
  *out = static_cast<uint32_t>(v);
  return true;
}

Config ParseArgs(int argc, char** argv) {
  Config cfg;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--min_n" && i + 1 < argc) {
      if (!ParseIntArg(argv[++i], &cfg.min_n)) {
        throw std::invalid_argument("Invalid --min_n value");
      }
    } else if (arg == "--max_n" && i + 1 < argc) {
      if (!ParseIntArg(argv[++i], &cfg.max_n)) {
        throw std::invalid_argument("Invalid --max_n value");
      }
    } else if (arg == "--trials" && i + 1 < argc) {
      if (!ParseIntArg(argv[++i], &cfg.contains_trials)) {
        throw std::invalid_argument("Invalid --trials value");
      }
    } else if (arg == "--word_len" && i + 1 < argc) {
      if (!ParseIntArg(argv[++i], &cfg.word_len)) {
        throw std::invalid_argument("Invalid --word_len value");
      }
    } else if (arg == "--seed" && i + 1 < argc) {
      if (!ParseU32Arg(argv[++i], &cfg.seed)) {
        throw std::invalid_argument("Invalid --seed value");
      }
    } else {
      throw std::invalid_argument("Unknown argument: " + arg);
    }
  }

  if (cfg.min_n < 2 || cfg.max_n < cfg.min_n) {
    throw std::invalid_argument("Require 2 <= min_n <= max_n");
  }
  if (cfg.contains_trials <= 0 || cfg.word_len <= 0) {
    throw std::invalid_argument("Require positive --trials and --word_len");
  }
  return cfg;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Config cfg = ParseArgs(argc, argv);

    std::cout << "Benchmark Schreier-Sims on S_n\n";
    std::cout << "min_n=" << cfg.min_n << " max_n=" << cfg.max_n
              << " trials=" << cfg.contains_trials
              << " word_len=" << cfg.word_len
              << " seed=" << cfg.seed << "\n\n";

    std::cout << std::left << std::setw(5) << "n" << std::setw(16)
              << "expected|S_n|" << std::setw(16) << "computed|G|"
              << std::setw(12) << "build_ms" << std::setw(14)
              << "contains_ms" << std::setw(18) << "contains_us/op"
              << std::setw(10) << "ok_hits" << "\n";

    for (int n = cfg.min_n; n <= cfg.max_n; ++n) {
      std::mt19937 rng(cfg.seed + static_cast<uint32_t>(n));

      const auto gens = SymmetricGenerators(static_cast<std::size_t>(n));
      std::vector<groups::Permutation> gen_inv;
      gen_inv.reserve(gens.size());
      for (const auto& g : gens) {
        gen_inv.push_back(g.Inverse());
      }

      groups::SchreierSims ss(static_cast<std::size_t>(n), gens);
      const auto t_build_0 = std::chrono::steady_clock::now();
      ss.Build();
      const auto t_build_1 = std::chrono::steady_clock::now();
      const double build_ms =
          std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
              t_build_1 - t_build_0)
              .count();

      int ok_hits = 0;
      const auto t_contains_0 = std::chrono::steady_clock::now();
      for (int trial = 0; trial < cfg.contains_trials; ++trial) {
        const groups::Permutation x =
            RandomWord(gens, gen_inv, cfg.word_len, rng);
        if (ss.Contains(x)) ++ok_hits;
      }
      const auto t_contains_1 = std::chrono::steady_clock::now();
      const double contains_ms =
          std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
              t_contains_1 - t_contains_0)
              .count();
      const double us_per_op =
          (contains_ms * 1000.0) / static_cast<double>(cfg.contains_trials);

      const uint64_t expected_order = FactorialU64(static_cast<std::size_t>(n));
      const uint64_t computed_order = ss.Order();

      std::cout << std::left << std::setw(5) << n << std::setw(16)
                << expected_order << std::setw(16) << computed_order
                << std::setw(12) << std::fixed << std::setprecision(2) << build_ms
                << std::setw(14) << contains_ms << std::setw(18) << us_per_op
                << std::setw(10) << ok_hits << "\n";
    }
  } catch (const std::exception& ex) {
    std::cerr << "group_benchmark error: " << ex.what() << "\n";
    std::cerr << "usage: group_benchmark [--min_n N] [--max_n N] [--trials T] "
                 "[--word_len L] [--seed S]\n";
    return 1;
  }

  return 0;
}
