"""Subprocess orchestration for CCZ/EA equivalence searches."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import Any, Optional

try:
    from ._ccz_core import _core_auto, _core_equivalence
    from ._ccz_sage import (
        _affine_result_to_sage_matrix,
        _auto_group_for_core,
        _auto_result_order,
    )
except ImportError:  # pragma: no cover - direct module import from python/
    from _ccz_core import _core_auto, _core_equivalence
    from _ccz_sage import (
        _affine_result_to_sage_matrix,
        _auto_group_for_core,
        _auto_result_order,
    )

_HERE = Path(__file__).resolve().parent
_WORKER = _HERE / "auto_worker.py"


def _print(verbose: bool, message: str) -> None:
    if verbose:
        print(message, flush=True)


def _seed(auto_option: Any, n_bits: int, m_bits: int):
    if isinstance(auto_option, bool):
        return None
    return _auto_group_for_core(auto_option, n_bits, m_bits)


def _to_matrix(affine, n_bits, m_bits, sage_constructors, swapped=False):
    return _affine_result_to_sage_matrix(
        affine, n_bits, m_bits, sage_constructors, invert=swapped
    )


def _run_equivalence(
    mode, left_tt, right_tt, n_bits, m_bits, time_limit, min_active_hyperplanes,
    seed, sage_constructors, swapped=False,
):
    f_tt, g_tt = (right_tt, left_tt) if swapped else (left_tt, right_tt)
    affine = _core_equivalence(
        mode, f_tt, g_tt, n_bits, m_bits, time_limit, min_active_hyperplanes, seed
    )
    return _to_matrix(affine, n_bits, m_bits, sage_constructors, swapped)


def _spawn(name: str, kind: str, payload: dict[str, Any], swapped: bool = False):
    out_fd, out_name = tempfile.mkstemp(prefix="ccz-worker-", suffix=".stdout")
    err_fd, err_name = tempfile.mkstemp(prefix="ccz-worker-", suffix=".stderr")
    os.close(out_fd)
    os.close(err_fd)
    out = Path(out_name)
    err = Path(err_name)
    with out.open("w", encoding="utf-8") as stdout, err.open(
        "w", encoding="utf-8"
    ) as stderr:
        proc = subprocess.Popen(
            [sys.executable, str(_WORKER)],
            cwd=str(_HERE.parent),
            stdin=subprocess.PIPE,
            stdout=stdout,
            stderr=stderr,
            text=True,
        )
    if proc.stdin is not None:
        proc.stdin.write(json.dumps(payload))
        proc.stdin.close()
    return {"name": name, "kind": kind, "proc": proc, "out": out, "err": err, "swapped": swapped}


def _read(task):
    task["proc"].wait()
    try:
        response = json.loads(task["out"].read_text(encoding="utf-8"))
    except Exception as exc:
        return False, {"error_message": f"worker produced invalid JSON: {exc}"}
    return (
        task["proc"].returncode == 0 and isinstance(response, dict) and response.get("ok"),
        response.get("result") if isinstance(response, dict) and response.get("ok") else response,
    )


def _stop(tasks):
    for task in tasks:
        if task["proc"].poll() is None:
            task["proc"].terminate()
            try:
                task["proc"].wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                task["proc"].kill()
                task["proc"].wait()
        for key in ("out", "err"):
            try:
                task[key].unlink()
            except (FileNotFoundError, PermissionError):
                pass
    tasks.clear()


def _auto_payload(mode, truth_table, n_bits, m_bits, time_limit, min_active_hyperplanes):
    return {
        "mode": mode,
        "truth_table": truth_table,
        "n_bits": n_bits,
        "m_bits": m_bits,
        "time_limit_seconds": time_limit,
        "min_active_hyperplanes": min_active_hyperplanes,
    }


def _equiv_payload(mode, left_tt, right_tt, n_bits, m_bits, time_limit,
                   min_active_hyperplanes, seed):
    return {
        "mode": f"{mode}_equivalence",
        "truth_table_f": left_tt,
        "truth_table_g": right_tt,
        "n_bits": n_bits,
        "m_bits": m_bits,
        "time_limit_seconds": time_limit,
        "min_active_hyperplanes": min_active_hyperplanes,
        "auto_group": seed,
    }


def _run_sequential_equivalence(
    mode: str,
    left_tt: list[int],
    right_tt: list[int],
    n_bits: int,
    m_bits: int,
    time_limit: Optional[float],
    min_active_hyperplanes: Optional[int],
    left_auto: Any,
    right_auto: Any,
    verbose: bool,
    sage_constructors: tuple[Any, Any, Any, Any],
):
    seed = _seed(right_auto, n_bits, m_bits)
    swapped = False
    if seed is None:
        seed = _seed(left_auto, n_bits, m_bits)
        swapped = seed is not None
    if seed is None:
        if right_auto is True:
            _print(verbose, "Computing right automorphism group")
            seed = _core_auto(
                mode, right_tt, n_bits, m_bits, time_limit, min_active_hyperplanes
            )
        elif left_auto is True:
            swapped = True
            _print(verbose, "Computing left automorphism group")
            seed = _core_auto(
                mode, left_tt, n_bits, m_bits, time_limit, min_active_hyperplanes
            )

    _print(
        verbose,
        "Running equivalence search with left auto seed" if swapped else "Running equivalence search",
    )
    return _run_equivalence(
        mode, left_tt, right_tt, n_bits, m_bits, time_limit,
        min_active_hyperplanes, seed, sage_constructors, swapped,
    )


def _run_parallel_equivalence(
    mode: str,
    left_tt: list[int],
    right_tt: list[int],
    n_bits: int,
    m_bits: int,
    time_limit: Optional[float],
    min_active_hyperplanes: Optional[int],
    left_auto: Any,
    right_auto: Any,
    verbose: bool,
    sage_constructors: tuple[Any, Any, Any, Any],
):
    if _seed(left_auto, n_bits, m_bits) is not None or _seed(
        right_auto, n_bits, m_bits
    ) is not None:
        return _run_sequential_equivalence(
            mode, left_tt, right_tt, n_bits, m_bits, time_limit,
            min_active_hyperplanes, left_auto, right_auto, verbose, sage_constructors,
        )

    if not _WORKER.is_file():
        return _run_sequential_equivalence(
            mode, left_tt, right_tt, n_bits, m_bits, time_limit,
            min_active_hyperplanes, left_auto, right_auto, verbose, sage_constructors,
        )

    tasks = []

    def add_equiv(name, seed, swapped=False):
        f_tt, g_tt = (right_tt, left_tt) if swapped else (left_tt, right_tt)
        tasks.append(_spawn(
            name, "equiv",
            _equiv_payload(mode, f_tt, g_tt, n_bits, m_bits, time_limit,
                           min_active_hyperplanes, seed),
            swapped,
        ))
        _print(verbose, f"Started {name}")

    right_seed = _seed(right_auto, n_bits, m_bits)
    left_seed = _seed(left_auto, n_bits, m_bits)
    if right_seed is not None:
        add_equiv("right-seeded equivalence", right_seed)
    elif left_seed is not None:
        add_equiv("left-seeded equivalence", left_seed, True)
    else:
        add_equiv("equivalence", None)

    tasks.append(_spawn(
        "left auto", "auto_left",
        _auto_payload(mode, left_tt, n_bits, m_bits, time_limit,
                      min_active_hyperplanes),
    ))
    _print(verbose, "Started left auto")
    tasks.append(_spawn(
        "right auto", "auto_right",
        _auto_payload(mode, right_tt, n_bits, m_bits, time_limit,
                      min_active_hyperplanes),
    ))
    _print(verbose, "Started right auto")

    try:
        while tasks:
            done = [task for task in tasks if task["proc"].poll() is not None]
            if not done:
                time.sleep(0.05)
                continue
            done.sort(key=lambda task: task["kind"] != "equiv")
            task = done[0]
            tasks.remove(task)
            ok, result = _read(task)
            _stop([task])
            if not ok:
                _print(verbose, f"{task['name']} failed: {result.get('error_message', result)}")
                continue
            if task["kind"] == "equiv":
                _print(verbose, f"{task['name']} finished")
                _stop(tasks)
                return _to_matrix(
                    result, n_bits, m_bits, sage_constructors, task["swapped"]
                )

            swapped = task["kind"] == "auto_left"
            _print(verbose, f"{task['name']} finished with order {_auto_result_order(result)}")
            _stop(tasks)
            return _run_equivalence(
                mode, left_tt, right_tt, n_bits, m_bits, time_limit,
                min_active_hyperplanes, _auto_group_for_core(result, n_bits, m_bits),
                sage_constructors, swapped,
            )
    finally:
        _stop(tasks)
    return None
