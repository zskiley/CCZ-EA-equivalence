#!/usr/bin/env python3
from __future__ import annotations

import argparse
import time

import ccz


def poly_degree(p: int) -> int:
    if p == 0:
        return -1
    return p.bit_length() - 1


def poly_mod(a: int, b: int) -> int:
    if b == 0:
        return a
    deg_b = poly_degree(b)
    deg_a = poly_degree(a)
    while deg_a >= deg_b and deg_a >= 0:
        a ^= b << (deg_a - deg_b)
        deg_a = poly_degree(a)
    return a


def is_irreducible_over_gf2(poly: int, n: int) -> bool:
    if n <= 0 or n >= 31:
        return False
    if ((poly >> n) & 1) == 0:
        return False
    if (poly & 1) == 0:
        return False
    if poly_degree(poly) != n:
        return False

    for d in range(1, n // 2 + 1):
        begin = 1 << d
        end = 1 << (d + 1)
        for g in range(begin, end, 2):
            if poly_mod(poly, g) == 0:
                return False
    return True


def default_mod_poly(n: int) -> int:
    if n <= 0 or n >= 31:
        raise ValueError("n must be in [1, 30]")
    lower = (1 << n) | 1
    upper = 1 << (n + 1)
    for poly in range(lower, upper, 2):
        if is_irreducible_over_gf2(poly, n):
            return poly
    raise RuntimeError("no irreducible polynomial found")


class GF2n:
    def __init__(self, n: int, mod_poly: int):
        self.n = n
        self.mod_poly = mod_poly
        self.mask = (1 << n) - 1

    def mul(self, a: int, b: int) -> int:
        res = 0
        x = a
        y = b
        for _ in range(self.n):
            if (y & 1) != 0:
                res ^= x
            y >>= 1
            carry = ((x >> (self.n - 1)) & 1) != 0
            x <<= 1
            if carry:
                x ^= self.mod_poly
            x &= self.mask
        return res & self.mask

    def pow(self, a: int, e: int) -> int:
        result = 1
        base = a
        exp = e
        while exp > 0:
            if (exp & 1) != 0:
                result = self.mul(result, base)
            base = self.mul(base, base)
            exp >>= 1
        return result


def power_truth_table(n: int, exponent: int) -> list[int]:
    field = GF2n(n, default_mod_poly(n))
    q = 1 << n
    mod = (1 << n) - 1
    e = exponent % mod
    if e < 0:
        e += mod
    table = [0] * q
    for x in range(q):
        if x == 0:
            table[x] = 0
        else:
            table[x] = field.pow(x, e)
    return table


def run_and_time(label: str, fn):
    t0 = time.perf_counter()
    out = fn()
    dt_ms = (time.perf_counter() - t0) * 1000.0
    print(f"{label}: {dt_ms:.1f} ms")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run CCZ/EA auto + equivalence via python bindings."
    )
    parser.add_argument("--n", type=int, default=9)
    parser.add_argument("--d1", type=int, default=3, help="Exponent for F(x)=x^d1")
    parser.add_argument("--d2", type=int, default=6, help="Exponent for G(x)=x^d2")
    parser.add_argument(
        "--mode",
        choices=[
            "all",
            "ccz_auto",
            "ea_auto",
            "ccz_equivalence",
            "ea_equivalence",
        ],
        default="all",
    )
    parser.add_argument(
        "--time-limit",
        type=float,
        default=None,
        help="Optional seconds for auto seeding in equivalence wrappers.",
    )
    args = parser.parse_args()

    print(f"building truth tables for n={args.n}, d1={args.d1}, d2={args.d2}")
    tt_f = run_and_time("Build F truth table", lambda: power_truth_table(args.n, args.d1))
    tt_g = run_and_time("Build G truth table", lambda: power_truth_table(args.n, args.d2))
    print("")

    if args.mode in ("all", "ccz_auto"):
        auto = run_and_time(
            "ccz_auto",
            lambda: ccz.ccz_auto(
                tt_f, time_limit_seconds=args.time_limit
            ),
        )
        print(f"ccz_auto order: {auto['order']}\n")

    if args.mode in ("all", "ea_auto"):
        auto = run_and_time(
            "ea_auto",
            lambda: ccz.ea_auto(
                tt_f, time_limit_seconds=args.time_limit
            ),
        )
        print(f"ea_auto order: {auto['order']}\n")

    if args.mode in ("all", "ccz_equivalence"):
        eq = run_and_time(
            "ccz_equivalence",
            lambda: ccz.ccz_equivalence(
                tt_f, tt_g, time_limit_seconds=args.time_limit
            ),
        )
        print(f"ccz_equivalence: {'YES' if eq is not None else 'NO'}")
        if eq is not None:
            print(f"map size: {len(eq)}")
        print("")

    if args.mode in ("all", "ea_equivalence"):
        eq = run_and_time(
            "ea_equivalence",
            lambda: ccz.ea_equivalence(
                tt_f, tt_g, time_limit_seconds=args.time_limit
            ),
        )
        print(f"ea_equivalence: {'YES' if eq is not None else 'NO'}")
        if eq is not None:
            print(f"map size: {len(eq)}")
        print("")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
