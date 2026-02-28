#!/usr/bin/env python3
"""Print only equivalence runtimes for available input representations."""

import time

import ccz


def _timed(label, fn):
    t0 = time.perf_counter()
    res = fn()
    ms = (time.perf_counter() - t0) * 1000.0
    print(f"{label}: {ms:.2f} ms")
    return res

def main() -> None:
    n_truth = 3
    tt_f = [x for x in range(1 << n_truth)]
    tt_g = [x for x in range(1 << n_truth)]
    _timed(
        "truth_table.ccz_equivalence",
        lambda: ccz.ccz_equivalence(tt_f, tt_g, n_truth),
    )
    _timed(
        "truth_table.ea_equivalence",
        lambda: ccz.ea_equivalence(tt_f, tt_g, n_truth),
    )

    n_callable = 3
    f_int = lambda x: x
    g_int = lambda x: x
    _timed(
        "callable.ccz_equivalence",
        lambda: ccz.ccz_equivalence(f_int, g_int, n_callable),
    )
    _timed(
        "callable.ea_equivalence",
        lambda: ccz.ea_equivalence(f_int, g_int, n_callable),
    )

    try:
        import galois  # type: ignore

        n_galois = 8
        F = galois.GF(2**n_galois)
        x = galois.Poly.Identity(F)
        f_ff = x**3
        g_ff = x**5
        print("galois.Poly x^3 vs x^5, n=", n_galois)
        auto_group_ccz = _timed(
            "galois.ccz_auto(g) (seed group, time-limited)",
            lambda: ccz.ccz_auto(g_ff),
        )
        print(
            "seed order=",
            auto_group_ccz["order"],
            "generators=",
            len(auto_group_ccz["generators"]),
            "complete=",
            auto_group_ccz["found_entire_group"],
        )
        auto_group_ea = _timed(
            "galois.ea_auto(g) (seed group, time-limited)",
            lambda: ccz.ea_auto(g_ff),
        )
        print(
            "ea seed order=",
            auto_group_ea["order"],
            "generators=",
            len(auto_group_ea["generators"]),
            "complete=",
            auto_group_ea["found_entire_group"],
        )
        res1 = _timed(
            "galois.ccz_equivalence (poly)",
            lambda: ccz.ccz_equivalence(
                f_ff, g_ff, time_limit_seconds=2.0, auto_group=auto_group_ccz
            ),
        )
        print(res1)
        res2 = _timed(
            "galois.ea_equivalence (poly)",
            lambda: ccz.ea_equivalence(
                f_ff, g_ff, time_limit_seconds=2.0, auto_group=auto_group_ea
            ),
        )
        print(res2)

        # Kasami exponent d = 2^(2i) - 2^i + 1 over GF(2^n).
        n_kasami = 9
        i_kasami = 2
        d_kasami = (1 << (2 * i_kasami)) - (1 << i_kasami) + 1
        F_kasami = galois.GF(2**n_kasami)
        x_kasami = galois.Poly.Identity(F_kasami)
        kasami = x_kasami**d_kasami
        print(
            f"galois.Poly Kasami, n={n_kasami}, i={i_kasami}, d={d_kasami}"
        )
        res3 = _timed(
            "galois.ccz_auto (Kasami poly, time-limited)",
            lambda: ccz.ccz_auto(kasami),
        )
        print(
            "kasami seed order=",
            res3["order"],
            "generators=",
            len(res3["generators"]),
            "complete=",
            res3["found_entire_group"],
        )
        _timed(
            "galois.ea_auto (Kasami poly, time-limited)",
            lambda: ccz.ea_auto(kasami, time_limit_seconds=2.0),
        )

    except Exception as exc:
        print("galois section failed:", exc)

    try:
        from sage.all import GF, PolynomialRing  # type: ignore

        n_sage = 3
        F_sage = GF(2**n_sage, name="a")
        R = PolynomialRing(F_sage, "x")
        x = R.gen()
        p = x**3
        q = x**6
        _timed(
            "sage.ccz_equivalence",
            lambda: ccz.ccz_equivalence(p, q, n_sage, field=F_sage),
        )
        _timed(
            "sage.ea_equivalence",
            lambda: ccz.ea_equivalence(p, q, n_sage, field=F_sage),
        )
    except Exception as exc:
        print("sage section failed:", exc)


if __name__ == "__main__":
    main()
