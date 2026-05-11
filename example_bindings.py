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
    tt_f = [x for x in range(8)]
    tt_g = [x for x in range(8)]
    _timed(
        "truth_table.ccz_equivalence",
        lambda: ccz.ccz_equivalence(tt_f, tt_g),
    )
    _timed(
        "truth_table.ea_equivalence",
        lambda: ccz.ea_equivalence(tt_f, tt_g),
    )

    try:
        import galois  # type: ignore

        n_galois = 8
        F = galois.GF(2**n_galois)
        x = galois.Poly.Identity(F)
        f_ff = x**3
        g_ff = x**5
        print("galois.Poly x^3 vs x^5, n=", n_galois)
        auto_group_ccz, ccz_complete = _timed(
            "galois.ccz_auto(g) (seed group, time-limited)",
            lambda: ccz.ccz_auto(g_ff),
        )
        print(
            "seed order=",
            auto_group_ccz.order(),
            "generators=",
            len(auto_group_ccz.gens()),
            "complete=",
            ccz_complete,
        )
        auto_group_ea, ea_complete = _timed(
            "galois.ea_auto(g) (seed group, time-limited)",
            lambda: ccz.ea_auto(g_ff),
        )
        print(
            "ea seed order=",
            auto_group_ea.order(),
            "generators=",
            len(auto_group_ea.gens()),
            "complete=",
            ea_complete,
        )
        res1 = _timed(
            "galois.ccz_equivalence (poly)",
            lambda: ccz.ccz_equivalence(
                f_ff, g_ff, time_limit=2.0, right_auto=auto_group_ccz
            ),
        )
        print(res1)
        res2 = _timed(
            "galois.ea_equivalence (poly)",
            lambda: ccz.ea_equivalence(
                f_ff, g_ff, time_limit=2.0, right_auto=auto_group_ea
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
        kasami_group, kasami_complete = _timed(
            "galois.ccz_auto (Kasami poly, time-limited)",
            lambda: ccz.ccz_auto(kasami),
        )
        print(
            "kasami seed order=",
            kasami_group.order(),
            "generators=",
            len(kasami_group.gens()),
            "complete=",
            kasami_complete,
        )
        _timed(
            "galois.ea_auto (Kasami poly, time-limited)",
            lambda: ccz.ea_auto(kasami, time_limit=2.0),
        )

    except Exception as exc:
        print("galois section failed:", exc)

    try:
        from sage.all import GF, PolynomialRing  # type: ignore

        F_sage = GF(2**3, name="a")
        R = PolynomialRing(F_sage, "x")
        x = R.gen()
        p = x**3
        q = x**6
        _timed(
            "sage.ccz_equivalence",
            lambda: ccz.ccz_equivalence(p, q),
        )
        _timed(
            "sage.ea_equivalence",
            lambda: ccz.ea_equivalence(p, q),
        )
    except Exception as exc:
        print("sage section failed:", exc)


if __name__ == "__main__":
    main()
