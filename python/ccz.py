"""High-level Python helpers for CCZ/EA algorithms.

Usage:
  1) Build extension once:
       python python/build.py
     (or from Sage's Python):
       sage -python python/build.py
  2) In Python/Sage:
       import sys
       sys.path.append(r".../repo/python")
       import ccz
"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import shutil
import sys
from typing import Any, Callable, Iterable, Optional

_HERE = Path(__file__).resolve().parent


def _prepare_windows_dll_paths() -> None:
    if os.name != "nt":
        return

    add_dir = getattr(os, "add_dll_directory", None)
    if add_dir is None:
        return

    seen: set[str] = set()

    def _try_add(path: Path) -> None:
        p = str(path)
        if not path.is_dir():
            return
        if p in seen:
            return
        try:
            add_dir(p)
            seen.add(p)
        except Exception:
            pass

    _try_add(_HERE)
    _try_add(Path(sys.base_prefix))
    _try_add(Path(sys.base_prefix) / "DLLs")
    _try_add(Path(sys.base_prefix) / "libs")

    gpp = shutil.which("g++")
    if gpp:
        _try_add(Path(gpp).resolve().parent)

    winget_root = Path.home() / "AppData" / "Local" / "Microsoft" / "WinGet" / "Packages"
    if winget_root.is_dir():
        for p in winget_root.glob("*\\mingw64\\bin"):
            _try_add(p)


_prepare_windows_dll_paths()


def _import_core_module():
    try:
        from . import ccz_bindings as _m
        return _m
    except Exception:
        import ccz_bindings as _m
        return _m


try:
    _core = _import_core_module()
except Exception:
    build_script = _HERE / "build.py"
    if build_script.is_file():
        try:
            subprocess.run(
                [sys.executable, str(build_script)],
                cwd=str(_HERE.parent),
                check=True,
            )
        except Exception:
            pass
    try:
        _core = _import_core_module()
    except Exception as exc:  # pragma: no cover
        raise ImportError(
            "Could not import ccz_bindings. Build it first with "
            "`python python/build.py` (or `sage -python python/build.py`)."
        ) from exc


IntFn = Callable[[int], int]


def _coerce_int_maybe(value: Any) -> Optional[int]:
    if callable(value):
        try:
            value = value()
        except Exception:
            return None
    try:
        return int(value)
    except Exception:
        return None


def _is_galois_poly(obj: Any) -> bool:
    cls = obj.__class__
    module = getattr(cls, "__module__", "")
    name = getattr(cls, "__name__", "")
    return module.startswith("galois") and name == "Poly" and hasattr(obj, "field")


def _infer_n_bits_from_field(field: Any) -> Optional[int]:
    degree = _coerce_int_maybe(getattr(field, "degree", None))
    if degree is not None and degree > 0:
        return degree

    order = _coerce_int_maybe(getattr(field, "order", None))
    if order is not None and order > 0 and (order & (order - 1)) == 0:
        n = order.bit_length() - 1
        if n > 0:
            return n

    cardinality = _coerce_int_maybe(getattr(field, "cardinality", None))
    if cardinality is not None and cardinality > 0 and (cardinality & (cardinality - 1)) == 0:
        n = cardinality.bit_length() - 1
        if n > 0:
            return n

    return None


def _field_from_callable(func: Any) -> Any:
    # galois.Poly exposes the field directly.
    ff = getattr(func, "field", None)
    if ff is not None:
        return ff

    try:
        parent = func.parent()
    except Exception:
        parent = None

    if parent is not None:
        try:
            base = parent.base_ring()
        except Exception:
            base = None

        if base is not None and hasattr(base, "fetch_int"):
            return base
        if hasattr(parent, "fetch_int"):
            return parent
    return None


def _int_to_field_element(x: int, field: Any) -> Any:
    if hasattr(field, "fetch_int"):
        return field.fetch_int(x)
    return field(x)


def _field_element_to_int(y: Any) -> int:
    if hasattr(y, "integer_representation"):
        return int(y.integer_representation())
    return int(y)


def _eval_callable_at_int(
    func: Any, x: int, n_bits: int, m_bits: int, field: Any
) -> int:
    y_mask = (1 << m_bits) - 1

    ff = field if field is not None else _field_from_callable(func)

    # If the caller explicitly provided a field, prefer field evaluation.
    # Otherwise, many simple callables like `lambda a: a**3` will happily
    # accept Python ints and compute integer arithmetic, which is not what we
    # want for finite-field/polynomial inputs.
    if ff is not None:
        if m_bits != n_bits:
            raise ValueError(
                "polynomial/field callable mode currently requires m_bits == n_bits"
            )
        y = func(_int_to_field_element(x, ff))
        return _field_element_to_int(y) & y_mask

    try:
        return int(func(x)) & y_mask
    except Exception:
        pass

    raise TypeError(
        "callable input could not be evaluated; pass `field=` for "
        "finite-field callables (or provide a callable/polynomial with "
        "built-in field metadata)"
    )


def _to_truth_table(
    values_or_fn: Iterable[int] | IntFn | Any,
    n_bits: int,
    m_bits: int,
    field: Any = None,
) -> list[int]:
    if n_bits <= 0:
        raise ValueError("n_bits must be > 0")
    if m_bits <= 0:
        raise ValueError("m_bits must be > 0")

    size = 1 << n_bits
    y_mask = (1 << m_bits) - 1

    if callable(values_or_fn):
        return [
            _eval_callable_at_int(values_or_fn, x, n_bits, m_bits, field)
            for x in range(size)
        ]

    values = [int(v) for v in values_or_fn]
    if len(values) != size:
        raise ValueError(f"truth table must have length 2^n_bits = {size}")
    return [v & y_mask for v in values]


def _infer_n_bits_from_truth_table_values(values_or_fn: Any) -> Optional[int]:
    if callable(values_or_fn):
        return None
    try:
        size = int(len(values_or_fn))
    except Exception:
        return None
    if size <= 0:
        return None
    if (size & (size - 1)) != 0:
        raise ValueError(
            "truth table length must be a power of two when `n_bits` is omitted"
        )
    return size.bit_length() - 1


def _resolve_bits_single(
    values_or_fn: Iterable[int] | IntFn | Any,
    n_bits: Optional[int],
    m_bits: Optional[int],
    field: Any = None,
) -> tuple[int, int, Any]:
    if _is_galois_poly(values_or_fn):
        if n_bits is not None:
            raise ValueError("`n_bits` must be omitted for galois polynomial input")
        ff = values_or_fn.field
        n = _infer_n_bits_from_field(ff)
        if n is None:
            raise ValueError("could not infer `n_bits` from galois polynomial field")
        m = n if m_bits is None else int(m_bits)
        if m != n:
            raise ValueError("galois polynomial mode requires m_bits == n_bits")
        return n, m, ff

    inferred_field = _field_from_callable(values_or_fn)
    if n_bits is None:
        candidate_field = inferred_field if inferred_field is not None else field
        if candidate_field is not None:
            n = _infer_n_bits_from_field(candidate_field)
            if n is not None:
                m = n if m_bits is None else int(m_bits)
                if m != n:
                    raise ValueError(
                        "polynomial/field callable mode currently requires m_bits == n_bits"
                    )
                return n, m, candidate_field
        inferred_n = _infer_n_bits_from_truth_table_values(values_or_fn)
        if inferred_n is not None:
            n = inferred_n
            m = n if m_bits is None else int(m_bits)
            return n, m, inferred_field

    if n_bits is None:
        raise ValueError(
            "`n_bits` is required unless input/field metadata is inferable, "
            "or the truth-table length is a power of two"
        )
    n = int(n_bits)
    m = n if m_bits is None else int(m_bits)
    return n, m, inferred_field


def _resolve_bits_pair(
    f_values_or_fn: Iterable[int] | IntFn | Any,
    g_values_or_fn: Iterable[int] | IntFn | Any,
    n_bits: Optional[int],
    m_bits: Optional[int],
    field: Any = None,
) -> tuple[int, int, Any, Any]:
    f_is_galois_poly = _is_galois_poly(f_values_or_fn)
    g_is_galois_poly = _is_galois_poly(g_values_or_fn)

    if f_is_galois_poly or g_is_galois_poly:
        if n_bits is not None:
            raise ValueError("`n_bits` must be omitted for galois polynomial input")

        inferred: list[int] = []
        f_field = f_values_or_fn.field if f_is_galois_poly else None
        g_field = g_values_or_fn.field if g_is_galois_poly else None

        if f_field is not None:
            f_n = _infer_n_bits_from_field(f_field)
            if f_n is None:
                raise ValueError("could not infer `n_bits` from first galois polynomial")
            inferred.append(f_n)
        if g_field is not None:
            g_n = _infer_n_bits_from_field(g_field)
            if g_n is None:
                raise ValueError("could not infer `n_bits` from second galois polynomial")
            inferred.append(g_n)

        if not inferred:
            raise ValueError("internal error: expected at least one galois polynomial")
        if len(set(inferred)) != 1:
            raise ValueError("galois polynomial inputs must have the same field degree")

        n = inferred[0]
        m = n if m_bits is None else int(m_bits)
        if m != n:
            raise ValueError("galois polynomial mode requires m_bits == n_bits")
        return n, m, f_field, g_field

    f_field = _field_from_callable(f_values_or_fn)
    g_field = _field_from_callable(g_values_or_fn)
    if n_bits is None:
        inferred: list[int] = []
        if f_field is not None:
            f_n = _infer_n_bits_from_field(f_field)
            if f_n is not None:
                inferred.append(f_n)
        if g_field is not None:
            g_n = _infer_n_bits_from_field(g_field)
            if g_n is not None:
                inferred.append(g_n)
        if not inferred and field is not None:
            field_n = _infer_n_bits_from_field(field)
            if field_n is not None:
                inferred.append(field_n)
        if inferred:
            if len(set(inferred)) != 1:
                raise ValueError("inputs must have the same inferred field degree")
            n = inferred[0]
            m = n if m_bits is None else int(m_bits)
            if m != n:
                raise ValueError(
                    "polynomial/field callable mode currently requires m_bits == n_bits"
                )
            return n, m, (f_field if f_field is not None else field), (
                g_field if g_field is not None else field
            )
        f_n = _infer_n_bits_from_truth_table_values(f_values_or_fn)
        g_n = _infer_n_bits_from_truth_table_values(g_values_or_fn)
        if f_n is not None and g_n is not None:
            if f_n != g_n:
                raise ValueError("truth tables must have the same inferred n_bits")
            n = f_n
            m = n if m_bits is None else int(m_bits)
            return n, m, f_field, g_field

    if n_bits is None:
        raise ValueError(
            "`n_bits` is required unless input/field metadata is inferable, "
            "or both truth-table lengths are powers of two"
        )
    n = int(n_bits)
    m = n if m_bits is None else int(m_bits)
    return n, m, f_field, g_field


def ccz_auto(
    values_or_fn: Iterable[int] | IntFn | Any,
    n_bits: Optional[int] = None,
    m_bits: Optional[int] = None,
    time_limit_seconds: Optional[float] = None,
    field: Any = None,
    min_active_hyperplanes: Optional[int] = None,
):
    n, m, inferred_field = _resolve_bits_single(values_or_fn, n_bits, m_bits, field)
    eval_field = inferred_field if inferred_field is not None else field
    tt = _to_truth_table(values_or_fn, n, m, field=eval_field)
    return _core.ccz_auto(tt, n, m, time_limit_seconds, min_active_hyperplanes)


def ea_auto(
    values_or_fn: Iterable[int] | IntFn | Any,
    n_bits: Optional[int] = None,
    m_bits: Optional[int] = None,
    time_limit_seconds: Optional[float] = None,
    field: Any = None,
    min_active_hyperplanes: Optional[int] = None,
):
    n, m, inferred_field = _resolve_bits_single(values_or_fn, n_bits, m_bits, field)
    eval_field = inferred_field if inferred_field is not None else field
    tt = _to_truth_table(values_or_fn, n, m, field=eval_field)
    return _core.ea_auto(tt, n, m, time_limit_seconds, min_active_hyperplanes)


def ccz_equivalence(
    f_values_or_fn: Iterable[int] | IntFn | Any,
    g_values_or_fn: Iterable[int] | IntFn | Any,
    n_bits: Optional[int] = None,
    m_bits: Optional[int] = None,
    time_limit_seconds: Optional[float] = None,
    field: Any = None,
    min_active_hyperplanes: Optional[int] = None,
    auto_group: Any = None,
):
    n, m, f_inferred_field, g_inferred_field = _resolve_bits_pair(
        f_values_or_fn, g_values_or_fn, n_bits, m_bits, field
    )
    f_eval_field = f_inferred_field if f_inferred_field is not None else field
    g_eval_field = g_inferred_field if g_inferred_field is not None else field
    f_tt = _to_truth_table(f_values_or_fn, n, m, field=f_eval_field)
    g_tt = _to_truth_table(g_values_or_fn, n, m, field=g_eval_field)
    return _core.ccz_equivalence(
        f_tt, g_tt, n, m, time_limit_seconds, min_active_hyperplanes, auto_group
    )


def ea_equivalence(
    f_values_or_fn: Iterable[int] | IntFn | Any,
    g_values_or_fn: Iterable[int] | IntFn | Any,
    n_bits: Optional[int] = None,
    m_bits: Optional[int] = None,
    time_limit_seconds: Optional[float] = None,
    field: Any = None,
    min_active_hyperplanes: Optional[int] = None,
    auto_group: Any = None,
):
    n, m, f_inferred_field, g_inferred_field = _resolve_bits_pair(
        f_values_or_fn, g_values_or_fn, n_bits, m_bits, field
    )
    f_eval_field = f_inferred_field if f_inferred_field is not None else field
    g_eval_field = g_inferred_field if g_inferred_field is not None else field
    f_tt = _to_truth_table(f_values_or_fn, n, m, field=f_eval_field)
    g_tt = _to_truth_table(g_values_or_fn, n, m, field=g_eval_field)
    return _core.ea_equivalence(
        f_tt, g_tt, n, m, time_limit_seconds, min_active_hyperplanes, auto_group
    )


__all__ = [
    "ccz_auto",
    "ea_auto",
    "ccz_equivalence",
    "ea_equivalence",
]
