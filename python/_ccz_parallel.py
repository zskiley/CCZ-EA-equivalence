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
        _auto_result_is_complete,
        _auto_result_order,
    )
except ImportError:  # pragma: no cover - direct module import from python/
    from _ccz_core import _core_auto, _core_equivalence
    from _ccz_sage import (
        _affine_result_to_sage_matrix,
        _auto_group_for_core,
        _auto_result_is_complete,
        _auto_result_order,
    )

_HERE = Path(__file__).resolve().parent


def _auto_worker_script_path() -> Path:
    return _HERE / "auto_worker.py"


def _start_worker_task(
    name: str,
    payload: dict[str, Any],
    *,
    task_type: str,
    swapped: bool = False,
) -> dict[str, Any]:
    stdout_fd, stdout_path_str = tempfile.mkstemp(
        prefix="ccz-worker-", suffix=".stdout"
    )
    stderr_fd, stderr_path_str = tempfile.mkstemp(
        prefix="ccz-worker-", suffix=".stderr"
    )
    os.close(stdout_fd)
    os.close(stderr_fd)
    stdout_path = Path(stdout_path_str)
    stderr_path = Path(stderr_path_str)

    with stdout_path.open("w", encoding="utf-8") as stdout_file, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_file:
        proc = subprocess.Popen(
            [sys.executable, str(_auto_worker_script_path())],
            cwd=str(_HERE.parent),
            stdin=subprocess.PIPE,
            stdout=stdout_file,
            stderr=stderr_file,
            text=True,
        )

    if proc.stdin is not None:
        proc.stdin.write(json.dumps(payload))
        proc.stdin.close()

    return {
        "name": name,
        "type": task_type,
        "swapped": swapped,
        "proc": proc,
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
    }


def _stop_worker_task(task: dict[str, Any]) -> None:
    proc = task["proc"]
    if proc.poll() is None:
        proc.terminate()
        try:
            proc.wait(timeout=2.0)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()


def _cleanup_worker_task(task: dict[str, Any]) -> None:
    for key in ("stdout_path", "stderr_path"):
        try:
            task[key].unlink()
        except (FileNotFoundError, PermissionError):
            pass


def _read_worker_task(task: dict[str, Any]) -> tuple[bool, Any]:
    proc = task["proc"]
    proc.wait()
    try:
        response = json.loads(task["stdout_path"].read_text(encoding="utf-8"))
    except Exception as exc:
        return False, {
            "error_type": type(exc).__name__,
            "error_message": "worker produced invalid JSON",
        }
    if proc.returncode != 0:
        return False, response
    if not isinstance(response, dict) or not response.get("ok"):
        return False, response
    return True, response.get("result")


def _auto_worker_payload(
    mode: str,
    truth_table: list[int],
    n_bits: int,
    m_bits: int,
    time_limit: Optional[float],
    min_active_hyperplanes: Optional[int],
) -> dict[str, Any]:
    return {
        "mode": mode,
        "truth_table": truth_table,
        "n_bits": n_bits,
        "m_bits": m_bits,
        "time_limit_seconds": time_limit,
        "min_active_hyperplanes": min_active_hyperplanes,
    }


def _equivalence_worker_payload(
    mode: str,
    left_tt: list[int],
    right_tt: list[int],
    n_bits: int,
    m_bits: int,
    time_limit: Optional[float],
    min_active_hyperplanes: Optional[int],
    auto_group: Any,
) -> dict[str, Any]:
    return {
        "mode": f"{mode}_equivalence",
        "truth_table_f": left_tt,
        "truth_table_g": right_tt,
        "n_bits": n_bits,
        "m_bits": m_bits,
        "time_limit_seconds": time_limit,
        "min_active_hyperplanes": min_active_hyperplanes,
        "auto_group": auto_group,
    }


def _supplied_auto_seed_and_info(auto_option: Any, n_bits: int, m_bits: int):
    if isinstance(auto_option, bool):
        return None, None

    candidate = auto_option
    complete = True
    if isinstance(candidate, tuple) and len(candidate) == 2:
        candidate, complete = candidate

    seed = _auto_group_for_core(candidate, n_bits, m_bits)
    order: Optional[int] = None
    if isinstance(candidate, dict):
        order = _auto_result_order(candidate)
        complete = _auto_result_is_complete(candidate)
    elif hasattr(candidate, "order"):
        try:
            order = int(candidate.order())
        except Exception:
            order = None
    return seed, {"order": order, "complete": bool(complete)}


def _print_verbose(verbose: bool, message: str) -> None:
    if verbose:
        print(message, flush=True)


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
    left_seed, left_info = _supplied_auto_seed_and_info(left_auto, n_bits, m_bits)
    right_seed, right_info = _supplied_auto_seed_and_info(right_auto, n_bits, m_bits)

    if left_auto is True:
        _print_verbose(verbose, "Computing left automorphism group")
        left_result = _core_auto(
            mode, left_tt, n_bits, m_bits, time_limit, min_active_hyperplanes
        )
        left_seed = left_result
        left_info = {
            "order": _auto_result_order(left_result),
            "complete": _auto_result_is_complete(left_result),
        }
    if right_auto is True:
        _print_verbose(verbose, "Computing right automorphism group")
        right_result = _core_auto(
            mode, right_tt, n_bits, m_bits, time_limit, min_active_hyperplanes
        )
        right_seed = right_result
        right_info = {
            "order": _auto_result_order(right_result),
            "complete": _auto_result_is_complete(right_result),
        }

    if (
        left_info is not None
        and right_info is not None
        and left_info.get("complete")
        and right_info.get("complete")
        and left_info.get("order") is not None
        and right_info.get("order") is not None
        and left_info["order"] != right_info["order"]
    ):
        _print_verbose(verbose, "Complete auto groups have different orders")
        return None

    use_left = False
    seed = right_seed
    if right_seed is None and left_seed is not None:
        use_left = True
        seed = left_seed
    elif right_seed is not None and left_seed is not None:
        left_order = left_info.get("order") if left_info else None
        right_order = right_info.get("order") if right_info else None
        if (
            left_order is not None
            and right_order is not None
            and left_order > right_order
        ):
            use_left = True
            seed = left_seed

    if use_left:
        _print_verbose(verbose, "Running equivalence search with left auto seed")
        affine = _core_equivalence(
            mode,
            right_tt,
            left_tt,
            n_bits,
            m_bits,
            time_limit,
            min_active_hyperplanes,
            seed,
        )
        return _affine_result_to_sage_matrix(
            affine, n_bits, m_bits, sage_constructors, invert=True
        )

    _print_verbose(verbose, "Running equivalence search")
    affine = _core_equivalence(
        mode,
        left_tt,
        right_tt,
        n_bits,
        m_bits,
        time_limit,
        min_active_hyperplanes,
        seed,
    )
    return _affine_result_to_sage_matrix(affine, n_bits, m_bits, sage_constructors)


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
    if not _auto_worker_script_path().is_file():
        return _run_sequential_equivalence(
            mode,
            left_tt,
            right_tt,
            n_bits,
            m_bits,
            time_limit,
            min_active_hyperplanes,
            left_auto,
            right_auto,
            verbose,
            sage_constructors,
        )

    tasks: list[dict[str, Any]] = []
    left_seed, _left_info = _supplied_auto_seed_and_info(left_auto, n_bits, m_bits)
    right_seed, _right_info = _supplied_auto_seed_and_info(right_auto, n_bits, m_bits)

    def launch_equivalence(name: str, seed: Any, swapped: bool = False) -> None:
        l_tt, r_tt = (right_tt, left_tt) if swapped else (left_tt, right_tt)
        tasks.append(
            _start_worker_task(
                name,
                _equivalence_worker_payload(
                    mode,
                    l_tt,
                    r_tt,
                    n_bits,
                    m_bits,
                    time_limit,
                    min_active_hyperplanes,
                    seed,
                ),
                task_type="equivalence",
                swapped=swapped,
            )
        )
        _print_verbose(verbose, f"Started {name}")

    def stop_remaining() -> None:
        for remaining in tasks:
            _stop_worker_task(remaining)
            _cleanup_worker_task(remaining)
        tasks.clear()

    def run_seeded_after_auto(
        name: str, seed: Any, *, swapped: bool = False
    ):
        seed = _auto_group_for_core(seed, n_bits, m_bits)
        _print_verbose(verbose, f"Running {name}")
        l_tt, r_tt = (right_tt, left_tt) if swapped else (left_tt, right_tt)
        affine = _core_equivalence(
            mode,
            l_tt,
            r_tt,
            n_bits,
            m_bits,
            time_limit,
            min_active_hyperplanes,
            seed,
        )
        return _affine_result_to_sage_matrix(
            affine,
            n_bits,
            m_bits,
            sage_constructors,
            invert=swapped,
        )

    launch_equivalence("equivalence", right_seed)
    if left_seed is not None:
        launch_equivalence("left-seeded equivalence", left_seed, swapped=True)

    if left_auto is True:
        tasks.append(
            _start_worker_task(
                "left auto",
                _auto_worker_payload(
                    mode,
                    left_tt,
                    n_bits,
                    m_bits,
                    time_limit,
                    min_active_hyperplanes,
                ),
                task_type="auto_left",
            )
        )
        _print_verbose(verbose, "Started left auto")
    if right_auto is True:
        tasks.append(
            _start_worker_task(
                "right auto",
                _auto_worker_payload(
                    mode,
                    right_tt,
                    n_bits,
                    m_bits,
                    time_limit,
                    min_active_hyperplanes,
                ),
                task_type="auto_right",
            )
        )
        _print_verbose(verbose, "Started right auto")

    try:
        while tasks:
            finished = [task for task in tasks if task["proc"].poll() is not None]
            if not finished:
                time.sleep(0.05)
                continue

            equivalence_finished = [
                task for task in finished if task["type"] == "equivalence"
            ]
            task = equivalence_finished[0] if equivalence_finished else finished[0]

            tasks.remove(task)
            ok, result = _read_worker_task(task)
            _cleanup_worker_task(task)
            if not ok:
                _print_verbose(
                    verbose,
                    f"{task['name']} failed: {result.get('error_message', result)}",
                )
                continue

            if task["type"] == "equivalence":
                _print_verbose(verbose, f"{task['name']} finished")
                matrix_obj = _affine_result_to_sage_matrix(
                    result,
                    n_bits,
                    m_bits,
                    sage_constructors,
                    invert=bool(task["swapped"]),
                )
                stop_remaining()
                return matrix_obj

            if task["type"] == "auto_left":
                left_order = _auto_result_order(result)
                _print_verbose(
                    verbose,
                    f"Left auto finished with order {left_order}",
                )
                stop_remaining()
                return run_seeded_after_auto(
                    "left-seeded equivalence", result, swapped=True
                )

            if task["type"] == "auto_right":
                right_order = _auto_result_order(result)
                _print_verbose(
                    verbose,
                    f"Right auto finished with order {right_order}",
                )
                stop_remaining()
                return run_seeded_after_auto("right-seeded equivalence", result)
    finally:
        stop_remaining()

    return None
