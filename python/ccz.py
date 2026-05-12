"""Public Python API for CCZ/EA algorithms."""

from __future__ import annotations

from typing import Any, Optional

try:
    from ._ccz_core import _core_auto
    from ._ccz_inputs import _normalize_inputs
    from ._ccz_parallel import _run_parallel_equivalence, _run_sequential_equivalence
    from ._ccz_sage import _auto_result_to_sage_group, _sage_constructors
except ImportError:  # pragma: no cover - direct module import from python/
    from _ccz_core import _core_auto
    from _ccz_inputs import _normalize_inputs
    from _ccz_parallel import _run_parallel_equivalence, _run_sequential_equivalence
    from _ccz_sage import _auto_result_to_sage_group, _sage_constructors


def _resolve_time_limit(
    time_limit: Optional[float],
    kwargs: dict[str, Any],
) -> Optional[float]:
    if "time_limit_seconds" in kwargs:
        legacy_time_limit = kwargs.pop("time_limit_seconds")
        if time_limit is not None:
            raise TypeError(
                "pass only one of time_limit and time_limit_seconds"
            )
        time_limit = legacy_time_limit
    if kwargs:
        names = ", ".join(sorted(kwargs))
        raise TypeError(f"unexpected keyword argument(s): {names}")
    if time_limit is None:
        return None
    return float(time_limit)


def _resolve_equivalence_legacy_kwargs(
    kwargs: dict[str, Any],
    time_limit: Optional[float],
    parallel: bool,
    right_auto: Any,
) -> tuple[Optional[float], bool, Any]:
    if "parallel_auto_seed" in kwargs:
        parallel = bool(kwargs.pop("parallel_auto_seed"))
    if "auto_group" in kwargs:
        right_auto = kwargs.pop("auto_group")
    time_limit = _resolve_time_limit(time_limit, kwargs)
    return time_limit, parallel, right_auto


def ccz_auto(
    function_input: Any,
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    **kwargs: Any,
):
    time_limit = _resolve_time_limit(time_limit, kwargs)
    sage_constructors = _sage_constructors()
    truth_table, n, m = _normalize_inputs(function_input)
    auto_result = _core_auto(
        "ccz", truth_table, n, m, time_limit, min_active_hyperplanes
    )
    group = _auto_result_to_sage_group(auto_result, n, m, sage_constructors)
    return group, bool(auto_result.get("found_entire_group"))


def ea_auto(
    function_input: Any,
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    **kwargs: Any,
):
    time_limit = _resolve_time_limit(time_limit, kwargs)
    sage_constructors = _sage_constructors()
    truth_table, n, m = _normalize_inputs(function_input)
    auto_result = _core_auto(
        "ea", truth_table, n, m, time_limit, min_active_hyperplanes
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
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    verbose: bool = True,
    **kwargs: Any,
):
    time_limit, parallel, right_auto = _resolve_equivalence_legacy_kwargs(
        kwargs, time_limit, parallel, right_auto
    )
    sage_constructors = _sage_constructors()
    left_tt, n, m = _normalize_inputs(left_input)
    right_tt, right_n, right_m = _normalize_inputs(right_input)
    if (right_n, right_m) != (n, m):
        raise ValueError("function inputs must have the same n and m dimensions")
    if parallel:
        return _run_parallel_equivalence(
            "ccz",
            left_tt,
            right_tt,
            n,
            m,
            time_limit,
            min_active_hyperplanes,
            left_auto,
            right_auto,
            verbose,
            sage_constructors,
        )
    return _run_sequential_equivalence(
        "ccz",
        left_tt,
        right_tt,
        n,
        m,
        time_limit,
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
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    verbose: bool = True,
    **kwargs: Any,
):
    time_limit, parallel, right_auto = _resolve_equivalence_legacy_kwargs(
        kwargs, time_limit, parallel, right_auto
    )
    sage_constructors = _sage_constructors()
    left_tt, n, m = _normalize_inputs(left_input)
    right_tt, right_n, right_m = _normalize_inputs(right_input)
    if (right_n, right_m) != (n, m):
        raise ValueError("function inputs must have the same n and m dimensions")
    if parallel:
        return _run_parallel_equivalence(
            "ea",
            left_tt,
            right_tt,
            n,
            m,
            time_limit,
            min_active_hyperplanes,
            left_auto,
            right_auto,
            verbose,
            sage_constructors,
        )
    return _run_sequential_equivalence(
        "ea",
        left_tt,
        right_tt,
        n,
        m,
        time_limit,
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
