"""Native binding loading and raw CCZ/EA wrappers."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Optional

try:
    from ._ccz_inputs import _normalize_function_input, _normalize_function_pair
    from ._ccz_sage import _auto_group_for_core
except ImportError:  # pragma: no cover - direct module import from python/
    from _ccz_inputs import _normalize_function_input, _normalize_function_pair
    from _ccz_sage import _auto_group_for_core


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
            "`sage -python python/build.py`."
        ) from exc


def _core_equivalence(
    mode: str,
    left_tt: list[int],
    right_tt: list[int],
    n_bits: int,
    m_bits: int,
    time_limit: Optional[float],
    min_active_hyperplanes: Optional[int],
    auto_group: Any,
):
    fn = _core.ccz_equivalence if mode == "ccz" else _core.ea_equivalence
    return fn(
        left_tt,
        right_tt,
        n_bits,
        m_bits,
        time_limit,
        min_active_hyperplanes,
        [] if auto_group is None else auto_group,
    )


def _core_auto(
    mode: str,
    truth_table: list[int],
    n_bits: int,
    m_bits: int,
    time_limit: Optional[float],
    min_active_hyperplanes: Optional[int],
) -> dict[str, Any]:
    fn = _core.ccz_auto if mode == "ccz" else _core.ea_auto
    return fn(truth_table, n_bits, m_bits, time_limit, min_active_hyperplanes)


def _raw_ccz_auto(
    function_input: Any,
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
) -> tuple[dict[str, Any], int, int]:
    tt, n, m = _normalize_function_input(function_input)
    result = _core.ccz_auto(tt, n, m, time_limit, min_active_hyperplanes)
    return result, n, m


def _raw_ea_auto(
    function_input: Any,
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
) -> tuple[dict[str, Any], int, int]:
    tt, n, m = _normalize_function_input(function_input)
    result = _core.ea_auto(tt, n, m, time_limit, min_active_hyperplanes)
    return result, n, m


def _raw_ccz_equivalence(
    left_input: Any,
    right_input: Any,
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    auto_group: Any = None,
) -> tuple[Optional[dict[str, Any]], int, int]:
    left_tt, right_tt, n, m = _normalize_function_pair(left_input, right_input)
    core_auto_group = (
        [] if auto_group is None else _auto_group_for_core(auto_group, n, m)
    )
    result = _core.ccz_equivalence(
        left_tt,
        right_tt,
        n,
        m,
        time_limit,
        min_active_hyperplanes,
        core_auto_group,
    )
    return result, n, m


def _raw_ea_equivalence(
    left_input: Any,
    right_input: Any,
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    auto_group: Any = None,
) -> tuple[Optional[dict[str, Any]], int, int]:
    left_tt, right_tt, n, m = _normalize_function_pair(left_input, right_input)
    core_auto_group = (
        [] if auto_group is None else _auto_group_for_core(auto_group, n, m)
    )
    result = _core.ea_equivalence(
        left_tt,
        right_tt,
        n,
        m,
        time_limit,
        min_active_hyperplanes,
        core_auto_group,
    )
    return result, n, m
