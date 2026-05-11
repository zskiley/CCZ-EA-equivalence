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

import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile
import time
from typing import Any, Iterable, Optional

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


def _is_supported_polynomial_callable(obj: Any) -> bool:
    if not callable(obj):
        return False
    if _is_galois_poly(obj):
        return True
    return _field_from_callable(obj) is not None


def _validate_supported_input(values_or_fn: Any, name: str) -> None:
    if callable(values_or_fn) and not _is_supported_polynomial_callable(values_or_fn):
        raise TypeError(
            f"{name} callable input is not supported. "
            "Pass a truth table, galois.Poly, or Sage polynomial."
        )


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
    if ff is None:
        raise TypeError(
            "callable input must be a field-backed polynomial "
            "(galois.Poly or Sage polynomial)"
        )
    if m_bits != n_bits:
        raise ValueError(
            "polynomial/field callable mode currently requires m_bits == n_bits"
        )
    y = func(_int_to_field_element(x, ff))
    return _field_element_to_int(y) & y_mask


def _to_truth_table(
    values_or_fn: Iterable[int] | Any,
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
        raise ValueError("truth table length must be a power of two")
    return size.bit_length() - 1


def _infer_m_bits_from_truth_table_values(values_or_fn: Any) -> Optional[int]:
    if callable(values_or_fn):
        return None
    max_value = 0
    seen = False
    try:
        for raw in values_or_fn:
            v = int(raw)
            if v < 0:
                raise ValueError("truth table values must be non-negative integers")
            if v > max_value:
                max_value = v
            seen = True
    except TypeError:
        return None
    if not seen:
        return None
    return max(1, max_value.bit_length())


def _resolve_bits_single(
    values_or_fn: Iterable[int] | Any,
    n_bits: Optional[int],
    m_bits: Optional[int],
    field: Any = None,
) -> tuple[int, int, Any]:
    _validate_supported_input(values_or_fn, "values_or_fn")

    if _is_galois_poly(values_or_fn):
        if n_bits is not None:
            raise ValueError("internal error: unexpected explicit n_bits")
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
            if m_bits is None:
                inferred_m = _infer_m_bits_from_truth_table_values(values_or_fn)
                m = max(n, inferred_m if inferred_m is not None else 1)
            else:
                m = int(m_bits)
            return n, m, inferred_field

    if n_bits is None:
        raise ValueError(
            "could not infer dimensions from input; use a truth table "
            "(length 2^n) or a field-backed polynomial (galois/Sage)"
        )
    n = int(n_bits)
    m = n if m_bits is None else int(m_bits)
    return n, m, inferred_field


def _resolve_bits_pair(
    f_values_or_fn: Iterable[int] | Any,
    g_values_or_fn: Iterable[int] | Any,
    n_bits: Optional[int],
    m_bits: Optional[int],
    field: Any = None,
) -> tuple[int, int, Any, Any]:
    _validate_supported_input(f_values_or_fn, "f_values_or_fn")
    _validate_supported_input(g_values_or_fn, "g_values_or_fn")

    f_is_galois_poly = _is_galois_poly(f_values_or_fn)
    g_is_galois_poly = _is_galois_poly(g_values_or_fn)

    if f_is_galois_poly or g_is_galois_poly:
        if n_bits is not None:
            raise ValueError("internal error: unexpected explicit n_bits")

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
            if m_bits is None:
                f_m = _infer_m_bits_from_truth_table_values(f_values_or_fn)
                g_m = _infer_m_bits_from_truth_table_values(g_values_or_fn)
                m = max(
                    n,
                    f_m if f_m is not None else 1,
                    g_m if g_m is not None else 1,
                )
            else:
                m = int(m_bits)
            return n, m, f_field, g_field

    if n_bits is None:
        raise ValueError(
            "could not infer dimensions from inputs; pass two truth tables "
            "with lengths 2^n, or two field-backed polynomials "
            "(galois/Sage) over the same field size"
        )
    n = int(n_bits)
    m = n if m_bits is None else int(m_bits)
    return n, m, f_field, g_field


def _normalize_function_input(function_input: Any) -> tuple[list[int], int, int]:
    if isinstance(function_input, tuple) and len(function_input) == 3:
        values_or_fn, n_bits, m_bits = function_input
        n = int(n_bits)
        m = int(m_bits)
        _validate_supported_input(values_or_fn, "function_input[0]")
        field = _field_from_callable(values_or_fn)
        if _is_galois_poly(values_or_fn):
            field = values_or_fn.field
        return _to_truth_table(values_or_fn, n, m, field=field), n, m

    n, _m, inferred_field = _resolve_bits_single(
        function_input, None, None, None
    )
    return _to_truth_table(function_input, n, n, field=inferred_field), n, n


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


def _equivalence_auto_seed_time_limit_seconds_arg(
    time_limit_seconds: Optional[float],
) -> Optional[float]:
    if time_limit_seconds is None:
        return None
    return float(time_limit_seconds)


def _auto_worker_script_path() -> Path:
    return _HERE / "auto_worker.py"


def _parallel_auto_worker_affinities() -> list[Optional[list[int]]]:
    sched_getaffinity = getattr(os, "sched_getaffinity", None)
    if sched_getaffinity is None:
        return [None, None]

    try:
        allowed_cores = sorted(int(core) for core in sched_getaffinity(0))
    except Exception:
        return [None, None]

    if len(allowed_cores) < 2:
        return [None, None]

    return [[allowed_cores[0]], [allowed_cores[1]]]


def _compute_pair_autos_in_parallel(
    mode: str,
    f_tt: list[int],
    g_tt: list[int],
    n_bits: int,
    m_bits: int,
    time_limit_seconds: Optional[float],
    min_active_hyperplanes: Optional[int],
) -> tuple[Optional[dict[str, Any]], Optional[dict[str, Any]]]:
    worker_script = _auto_worker_script_path()
    if not worker_script.is_file():
        return None, None

    auto_time_limit_seconds = _equivalence_auto_seed_time_limit_seconds_arg(
        time_limit_seconds
    )
    stop_on_first_success = auto_time_limit_seconds is None
    worker_affinities = _parallel_auto_worker_affinities()
    payloads = [
        {
            "mode": mode,
            "truth_table": f_tt,
            "n_bits": n_bits,
            "m_bits": m_bits,
            "time_limit_seconds": auto_time_limit_seconds,
            "min_active_hyperplanes": min_active_hyperplanes,
            "cpu_affinity": worker_affinities[0],
        },
        {
            "mode": mode,
            "truth_table": g_tt,
            "n_bits": n_bits,
            "m_bits": m_bits,
            "time_limit_seconds": auto_time_limit_seconds,
            "min_active_hyperplanes": min_active_hyperplanes,
            "cpu_affinity": worker_affinities[1],
        },
    ]
    processes: list[subprocess.Popen[str]] = []
    worker_output_paths: list[tuple[Path, Path]] = []

    for _ in payloads:
        stdout_fd, stdout_path_str = tempfile.mkstemp(
            prefix="ccz-auto-worker-", suffix=".stdout"
        )
        stderr_fd, stderr_path_str = tempfile.mkstemp(
            prefix="ccz-auto-worker-", suffix=".stderr"
        )
        os.close(stdout_fd)
        os.close(stderr_fd)
        stdout_path = Path(stdout_path_str)
        stderr_path = Path(stderr_path_str)
        worker_output_paths.append((stdout_path, stderr_path))
        with stdout_path.open("w", encoding="utf-8") as stdout_file, stderr_path.open(
            "w", encoding="utf-8"
        ) as stderr_file:
            proc = subprocess.Popen(
                [sys.executable, str(worker_script)],
                cwd=str(_HERE.parent),
                stdin=subprocess.PIPE,
                stdout=stdout_file,
                stderr=stderr_file,
                text=True,
            )
        processes.append(proc)

    for proc, payload in zip(processes, payloads):
        if proc.stdin is None:
            continue
        proc.stdin.write(json.dumps(payload))
        proc.stdin.close()

    def _read_worker_result(
        proc: subprocess.Popen[str], stdout_path: Path
    ) -> Optional[dict[str, Any]]:
        proc.wait()
        if proc.returncode != 0:
            return None
        try:
            response = json.loads(stdout_path.read_text(encoding="utf-8"))
        except Exception:
            return None
        if not isinstance(response, dict) or not response.get("ok"):
            return None
        result = response.get("result")
        if not isinstance(result, dict):
            return None
        result = _normalize_auto_result(result)
        if "generators" not in result or "order" not in result:
            return None
        return result

    def _stop_worker(proc: subprocess.Popen[str]) -> None:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                proc.kill()

    results: list[Optional[dict[str, Any]]] = [None, None]
    try:
        if stop_on_first_success:
            pending = set(range(len(processes)))
            while pending:
                finished_any = False
                for idx in list(pending):
                    proc = processes[idx]
                    if proc.poll() is None:
                        continue
                    pending.remove(idx)
                    finished_any = True
                    result = _read_worker_result(proc, worker_output_paths[idx][0])
                    results[idx] = result
                    if result is not None:
                        for other_idx in pending:
                            _stop_worker(processes[other_idx])
                        return results[0], results[1]
                if not finished_any:
                    time.sleep(0.05)
            return results[0], results[1]

        for idx, proc in enumerate(processes):
            results[idx] = _read_worker_result(proc, worker_output_paths[idx][0])

        return results[0], results[1]
    finally:
        for proc in processes:
            _stop_worker(proc)
        for stdout_path, stderr_path in worker_output_paths:
            for path in (stdout_path, stderr_path):
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass


def _auto_result_order(auto_result: Optional[dict[str, Any]]) -> Optional[int]:
    if auto_result is None:
        return None
    try:
        return int(auto_result["order"])
    except Exception:
        return None


def _auto_result_is_complete(auto_result: Optional[dict[str, Any]]) -> bool:
    if auto_result is None:
        return False
    return bool(auto_result.get("found_entire_group"))


def _normalize_auto_result(auto_result: dict[str, Any]) -> dict[str, Any]:
    def _normalize_generator_list(value: Any) -> Any:
        if not isinstance(value, list):
            return value
        normalized_generators: list[dict[Any, Any]] = []
        for generator in value:
            if not isinstance(generator, dict):
                continue
            if "translation" in generator and "linear_cols" in generator:
                try:
                    normalized_generators.append(
                        {
                            "translation": int(generator["translation"]),
                            "linear_cols": [
                                int(v) for v in generator["linear_cols"]
                            ],
                        }
                    )
                except Exception:
                    continue
                continue
            normalized_generator: dict[int, int] = {}
            for k, v in generator.items():
                try:
                    normalized_generator[int(k)] = int(v)
                except Exception:
                    continue
            normalized_generators.append(normalized_generator)
        return normalized_generators

    auto_result["generators"] = _normalize_generator_list(auto_result.get("generators"))
    auto_result["graph_generators"] = _normalize_generator_list(
        auto_result.get("graph_generators")
    )
    return auto_result


def _auto_result_has_usable_seed(auto_result: Optional[dict[str, Any]]) -> bool:
    if auto_result is None:
        return False
    generators = auto_result.get("generators")
    if isinstance(generators, list) and len(generators) > 0:
        return True
    graph_generators = auto_result.get("graph_generators")
    return isinstance(graph_generators, list) and len(graph_generators) > 0


def _invert_equivalence_point_map(point_map: Optional[dict[int, int]]):
    if point_map is None:
        return None
    return {int(v): int(k) for k, v in point_map.items()}


def _matrix_entry_bit(value: Any) -> int:
    return int(value) & 1


def _sage_matrix_to_affine_generator(matrix_like: Any, ambient_dimension: int):
    matrix_obj = matrix_like.matrix() if hasattr(matrix_like, "matrix") else matrix_like
    expected = ambient_dimension + 1
    if matrix_obj.nrows() != expected or matrix_obj.ncols() != expected:
        raise ValueError("auto_group Sage generator has wrong matrix dimension")

    for col in range(ambient_dimension):
        if _matrix_entry_bit(matrix_obj[ambient_dimension, col]) != 0:
            raise ValueError("auto_group Sage generator is not affine-homogeneous")
    if _matrix_entry_bit(matrix_obj[ambient_dimension, ambient_dimension]) != 1:
        raise ValueError("auto_group Sage generator is not affine-homogeneous")

    translation = 0
    for row in range(ambient_dimension):
        translation |= _matrix_entry_bit(matrix_obj[row, ambient_dimension]) << row

    linear_cols: list[int] = []
    for col in range(ambient_dimension):
        encoded_col = 0
        for row in range(ambient_dimension):
            encoded_col |= _matrix_entry_bit(matrix_obj[row, col]) << row
        linear_cols.append(encoded_col)

    return {"translation": translation, "linear_cols": linear_cols}


def _auto_group_for_core(auto_group: Any, n_bits: int, m_bits: int) -> Any:
    if auto_group is None:
        return None

    candidate = auto_group
    if isinstance(candidate, tuple) and len(candidate) == 2:
        candidate = candidate[0]

    if hasattr(candidate, "gens"):
        ambient_dimension = n_bits + m_bits
        return [
            _sage_matrix_to_affine_generator(generator, ambient_dimension)
            for generator in candidate.gens()
        ]

    return auto_group


def _print_python_auto_seed_status(
    auto_result: dict[str, Any],
    time_limit_seconds: Optional[float],
    start_label: str,
) -> None:
    found_entire_group = _auto_result_is_complete(auto_result)
    order = _auto_result_order(auto_result)
    if order is None:
        return
    if found_entire_group:
        print("(Found entire auto group)", flush=True)
    else:
        print("(potentially incoplete auto group)", flush=True)
    print(f"Auto group size before equivalence search: {order}", flush=True)
    if not found_entire_group and not _auto_result_has_usable_seed(auto_result):
        limit_text = (
            "until one auto search finishes"
            if time_limit_seconds is None
            else str(time_limit_seconds)
        )
        print(
            "Consider increasing the time limit for the automorphism "
            f"search, time_limit_seconds = {limit_text}",
            flush=True,
        )
    print(f"{start_label}{order}", flush=True)


def _sage_constructors():
    try:
        from sage.all import GF, MatrixGroup, identity_matrix, matrix
    except Exception as exc:  # pragma: no cover
        raise ImportError(
            "ccz_auto and ea_auto now return Sage MatrixGroup objects. "
            "Run them with Sage Python, for example `sage -python`, or from "
            "inside Sage."
        ) from exc
    return GF, MatrixGroup, identity_matrix, matrix


def _affine_generator_to_sage_matrix(
    generator: dict[str, Any],
    ambient_dimension: int,
    field: Any,
    matrix_ctor: Any,
):
    translation = int(generator["translation"])
    linear_cols = [int(col) for col in generator["linear_cols"]]
    if len(linear_cols) != ambient_dimension:
        raise ValueError("affine generator has wrong ambient dimension")

    rows: list[list[int]] = []
    for row in range(ambient_dimension):
        rows.append(
            [
                (linear_cols[col] >> row) & 1
                for col in range(ambient_dimension)
            ]
            + [(translation >> row) & 1]
        )
    rows.append([0] * ambient_dimension + [1])
    return matrix_ctor(field, rows)


def _auto_result_to_sage_group(
    auto_result: dict[str, Any],
    n_bits: int,
    m_bits: int,
    sage_constructors: tuple[Any, Any, Any, Any],
):
    GF, MatrixGroup, identity_matrix, matrix_ctor = sage_constructors
    field = GF(2)
    ambient_dimension = n_bits + m_bits
    homogeneous_dimension = ambient_dimension + 1
    generators = [
        _affine_generator_to_sage_matrix(
            generator, ambient_dimension, field, matrix_ctor
        )
        for generator in auto_result.get("generators", [])
    ]
    if not generators:
        generators = [identity_matrix(field, homogeneous_dimension)]
    return MatrixGroup(generators)


def ccz_auto(
    function_input: Any,
    time_limit: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    **kwargs: Any,
):
    time_limit = _resolve_time_limit(time_limit, kwargs)
    sage_constructors = _sage_constructors()
    auto_result, n, m = _raw_ccz_auto(
        function_input, time_limit, min_active_hyperplanes
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
    auto_result, n, m = _raw_ea_auto(
        function_input, time_limit, min_active_hyperplanes
    )
    group = _auto_result_to_sage_group(auto_result, n, m, sage_constructors)
    return group, bool(auto_result.get("found_entire_group"))


def ccz_equivalence(
    f_values_or_fn: Iterable[int] | Any,
    g_values_or_fn: Iterable[int] | Any,
    time_limit_seconds: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    auto_group: Any = None,
    parallel_auto_seed: bool = True,
):
    n, m, f_inferred_field, g_inferred_field = _resolve_bits_pair(
        f_values_or_fn, g_values_or_fn, None, None, None
    )
    f_tt = _to_truth_table(f_values_or_fn, n, m, field=f_inferred_field)
    g_tt = _to_truth_table(g_values_or_fn, n, m, field=g_inferred_field)

    if parallel_auto_seed and auto_group is None:
        f_auto, g_auto = _compute_pair_autos_in_parallel(
            "ccz",
            f_tt,
            g_tt,
            n,
            m,
            time_limit_seconds,
            min_active_hyperplanes,
        )
        f_order = _auto_result_order(f_auto)
        g_order = _auto_result_order(g_auto)

        if (
            _auto_result_is_complete(f_auto)
            and _auto_result_is_complete(g_auto)
            and f_order is not None
            and g_order is not None
            and f_order != g_order
        ):
            return None

        if f_auto is not None and (
            g_auto is None
            or (f_order is not None and g_order is not None and f_order > g_order)
        ):
            selected_auto_group = f_auto if _auto_result_has_usable_seed(f_auto) else None
            if selected_auto_group is not None:
                _print_python_auto_seed_status(
                    f_auto,
                    time_limit_seconds,
                    "Starting equivalence search with auto group size: ",
                )
            swapped = _core.ccz_equivalence(
                g_tt,
                f_tt,
                n,
                m,
                time_limit_seconds,
                min_active_hyperplanes,
                selected_auto_group,
            )
            return _invert_equivalence_point_map(swapped)

        if _auto_result_has_usable_seed(g_auto):
            auto_group = g_auto
            _print_python_auto_seed_status(
                g_auto,
                time_limit_seconds,
                "Starting equivalence search with auto group size: ",
            )

    core_auto_group = _auto_group_for_core(auto_group, n, m)
    return _core.ccz_equivalence(
        f_tt,
        g_tt,
        n,
        m,
        time_limit_seconds,
        min_active_hyperplanes,
        core_auto_group,
    )


def ea_equivalence(
    f_values_or_fn: Iterable[int] | Any,
    g_values_or_fn: Iterable[int] | Any,
    time_limit_seconds: Optional[float] = None,
    min_active_hyperplanes: Optional[int] = None,
    auto_group: Any = None,
    parallel_auto_seed: bool = True,
):
    n, m, f_inferred_field, g_inferred_field = _resolve_bits_pair(
        f_values_or_fn, g_values_or_fn, None, None, None
    )
    f_tt = _to_truth_table(f_values_or_fn, n, m, field=f_inferred_field)
    g_tt = _to_truth_table(g_values_or_fn, n, m, field=g_inferred_field)

    if parallel_auto_seed and auto_group is None:
        f_auto, g_auto = _compute_pair_autos_in_parallel(
            "ea",
            f_tt,
            g_tt,
            n,
            m,
            time_limit_seconds,
            min_active_hyperplanes,
        )
        f_order = _auto_result_order(f_auto)
        g_order = _auto_result_order(g_auto)

        if (
            _auto_result_is_complete(f_auto)
            and _auto_result_is_complete(g_auto)
            and f_order is not None
            and g_order is not None
            and f_order != g_order
        ):
            return None

        if f_auto is not None and (
            g_auto is None
            or (f_order is not None and g_order is not None and f_order > g_order)
        ):
            selected_auto_group = f_auto if _auto_result_has_usable_seed(f_auto) else None
            if selected_auto_group is not None:
                _print_python_auto_seed_status(
                    f_auto,
                    time_limit_seconds,
                    "Starting EA equivalence search with auto group size: ",
                )
            swapped = _core.ea_equivalence(
                g_tt,
                f_tt,
                n,
                m,
                time_limit_seconds,
                min_active_hyperplanes,
                selected_auto_group,
            )
            return _invert_equivalence_point_map(swapped)

        if _auto_result_has_usable_seed(g_auto):
            auto_group = g_auto
            _print_python_auto_seed_status(
                g_auto,
                time_limit_seconds,
                "Starting EA equivalence search with auto group size: ",
            )

    core_auto_group = _auto_group_for_core(auto_group, n, m)
    return _core.ea_equivalence(
        f_tt,
        g_tt,
        n,
        m,
        time_limit_seconds,
        min_active_hyperplanes,
        core_auto_group,
    )


__all__ = [
    "ccz_auto",
    "ea_auto",
    "ccz_equivalence",
    "ea_equivalence",
]
