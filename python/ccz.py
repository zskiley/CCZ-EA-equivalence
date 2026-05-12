"""Public Python API for CCZ/EA algorithms."""

from __future__ import annotations

from typing import Any, Optional

try:
    from ._ccz_core import _core_auto, _core_equivalence
    from ._ccz_inputs import _normalize_inputs
    from ._ccz_parallel import _run_parallel_equivalence
    from ._ccz_sage import (
        _affine_result_to_sage_matrix,
        _auto_group_for_core,
        _auto_result_order,
        _auto_result_to_sage_group,
        _sage_constructors,
    )
except ImportError:  # pragma: no cover - direct module import from python/
    from _ccz_core import _core_auto, _core_equivalence
    from _ccz_inputs import _normalize_inputs
    from _ccz_parallel import _run_parallel_equivalence
    from _ccz_sage import (
        _affine_result_to_sage_matrix,
        _auto_group_for_core,
        _auto_result_order,
        _auto_result_to_sage_group,
        _sage_constructors,
    )


def _print(verbose: bool, message: str) -> None:
    if verbose:
        print(message, flush=True)


def _supplied_group(auto_option: Any) -> bool:
    return not isinstance(auto_option, bool)


def _seed(auto_option: Any, n_bits: int, m_bits: int):
    if not _supplied_group(auto_option):
        return None
    return _auto_group_for_core(auto_option, n_bits, m_bits)


def _order(auto_option: Any) -> int:
    if isinstance(auto_option, dict):
        order = _auto_result_order(auto_option)
    elif hasattr(auto_option, "order"):
        order = int(auto_option.order())
    else:
        order = None
    if order is None:
        raise ValueError("could not determine supplied auto group order")
    return order


def _run_equivalence(
    mode: str,
    left_tt: list[int],
    right_tt: list[int],
    n_bits: int,
    m_bits: int,
    min_active_hyperplanes: Optional[int],
    seed: Any,
    sage_constructors: tuple[Any, Any, Any, Any],
    swapped: bool = False,
):
    f_tt, g_tt = (right_tt, left_tt) if swapped else (left_tt, right_tt)
    affine = _core_equivalence(
        mode, f_tt, g_tt, n_bits, m_bits, None, min_active_hyperplanes, seed
    )
    return _affine_result_to_sage_matrix(
        affine, n_bits, m_bits, sage_constructors, invert=swapped
    )


def _run_sequential_equivalence(
    mode: str,
    left_tt: list[int],
    right_tt: list[int],
    n_bits: int,
    m_bits: int,
    time_limit_auto_search: Optional[float],
    min_active_hyperplanes: Optional[int],
    left_auto: Any,
    right_auto: Any,
    verbose: bool,
    sage_constructors: tuple[Any, Any, Any, Any],
):
    left_seed = _seed(left_auto, n_bits, m_bits)
    right_seed = _seed(right_auto, n_bits, m_bits)
    seed = right_seed
    swapped = False
    if left_seed is not None and right_seed is not None:
        if _order(left_auto) > _order(right_auto):
            seed = left_seed
            swapped = True
    elif left_seed is not None:
        seed = left_seed
        swapped = True

    if seed is None:
        if right_auto is True:
            _print(verbose, "Computing right automorphism group")
            seed = _core_auto(
                mode,
                right_tt,
                n_bits,
                m_bits,
                time_limit_auto_search,
                min_active_hyperplanes,
            )
        elif left_auto is True:
            swapped = True
            _print(verbose, "Computing left automorphism group")
            seed = _core_auto(
                mode,
                left_tt,
                n_bits,
                m_bits,
                time_limit_auto_search,
                min_active_hyperplanes,
            )

    _print(
        verbose,
        "Running equivalence search with left auto seed"
        if swapped
        else "Running equivalence search",
    )
    return _run_equivalence(
        mode,
        left_tt,
        right_tt,
        n_bits,
        m_bits,
        min_active_hyperplanes,
        seed,
        sage_constructors,
        swapped,
    )


def _resolve_time_limit_auto_search(
    time_limit_auto_search: Optional[float],
    kwargs: dict[str, Any],
) -> Optional[float]:
    for legacy_name in ("time_limit", "time_limit_seconds"):
        if legacy_name not in kwargs:
            continue
        legacy_time_limit = kwargs.pop(legacy_name)
        if time_limit_auto_search is not None:
            raise TypeError(
                "pass only one of time_limit_auto_search, "
                "time_limit, and time_limit_seconds"
            )
        time_limit_auto_search = legacy_time_limit
    if kwargs:
        names = ", ".join(sorted(kwargs))
        raise TypeError(f"unexpected keyword argument(s): {names}")
    if time_limit_auto_search is None:
        return None
    return float(time_limit_auto_search)


def _resolve_equivalence_legacy_kwargs(
    kwargs: dict[str, Any],
    time_limit_auto_search: Optional[float],
    parallel: bool,
    right_auto: Any,
) -> tuple[Optional[float], bool, Any]:
    if "parallel_auto_seed" in kwargs:
        parallel = bool(kwargs.pop("parallel_auto_seed"))
    if "auto_group" in kwargs:
        right_auto = kwargs.pop("auto_group")
    time_limit_auto_search = _resolve_time_limit_auto_search(
        time_limit_auto_search, kwargs
    )
    return time_limit_auto_search, parallel, right_auto


def ccz_auto(
    function_input: Any,
    time_limit_auto_search: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    **kwargs: Any,
):
    time_limit_auto_search = _resolve_time_limit_auto_search(
        time_limit_auto_search, kwargs
    )
    sage_constructors = _sage_constructors()
    truth_table, n, m = _normalize_inputs(function_input)
    auto_result = _core_auto(
        "ccz", truth_table, n, m, time_limit_auto_search, min_active_hyperplanes
    )
    group = _auto_result_to_sage_group(auto_result, n, m, sage_constructors)
    return group, bool(auto_result.get("found_entire_group"))


def ea_auto(
    function_input: Any,
    time_limit_auto_search: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    **kwargs: Any,
):
    time_limit_auto_search = _resolve_time_limit_auto_search(
        time_limit_auto_search, kwargs
    )
    sage_constructors = _sage_constructors()
    truth_table, n, m = _normalize_inputs(function_input)
    auto_result = _core_auto(
        "ea", truth_table, n, m, time_limit_auto_search, min_active_hyperplanes
    )
    group = _auto_result_to_sage_group(auto_result, n, m, sage_constructors)
    return group, bool(auto_result.get("found_entire_group"))


def ccz_equivalence(
    left_input: Any,
    right_input: Any,
    *,
    parallel: bool = True,
    left_auto: Any = False,
    right_auto: Any = True,
    time_limit_auto_search: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    verbose: bool = True,
    **kwargs: Any,
):
    time_limit_auto_search, parallel, right_auto = _resolve_equivalence_legacy_kwargs(
        kwargs, time_limit_auto_search, parallel, right_auto
    )
    sage_constructors = _sage_constructors()
    left_tt, n, m = _normalize_inputs(left_input)
    right_tt, right_n, right_m = _normalize_inputs(right_input)
    if (right_n, right_m) != (n, m):
        raise ValueError("function inputs must have the same n and m dimensions")
    if parallel and not (
        _supplied_group(left_auto) or _supplied_group(right_auto)
    ):
        return _run_parallel_equivalence(
            "ccz",
            left_tt,
            right_tt,
            n,
            m,
            time_limit_auto_search,
            min_active_hyperplanes,
            verbose,
            sage_constructors,
        )
    return _run_sequential_equivalence(
        "ccz",
        left_tt,
        right_tt,
        n,
        m,
        time_limit_auto_search,
        min_active_hyperplanes,
        left_auto,
        right_auto,
        verbose,
        sage_constructors,
    )


def ea_equivalence(
    left_input: Any,
    right_input: Any,
    *,
    parallel: bool = True,
    left_auto: Any = False,
    right_auto: Any = True,
    time_limit_auto_search: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    verbose: bool = True,
    **kwargs: Any,
):
    time_limit_auto_search, parallel, right_auto = _resolve_equivalence_legacy_kwargs(
        kwargs, time_limit_auto_search, parallel, right_auto
    )
    sage_constructors = _sage_constructors()
    left_tt, n, m = _normalize_inputs(left_input)
    right_tt, right_n, right_m = _normalize_inputs(right_input)
    if (right_n, right_m) != (n, m):
        raise ValueError("function inputs must have the same n and m dimensions")
    if parallel and not (
        _supplied_group(left_auto) or _supplied_group(right_auto)
    ):
        return _run_parallel_equivalence(
            "ea",
            left_tt,
            right_tt,
            n,
            m,
            time_limit_auto_search,
            min_active_hyperplanes,
            verbose,
            sage_constructors,
        )
    return _run_sequential_equivalence(
        "ea",
        left_tt,
        right_tt,
        n,
        m,
        time_limit_auto_search,
        min_active_hyperplanes,
        left_auto,
        right_auto,
        verbose,
        sage_constructors,
    )


__all__ = [
    "ccz_auto",
    "ea_auto",
    "ccz_equivalence",
    "ea_equivalence",
]
