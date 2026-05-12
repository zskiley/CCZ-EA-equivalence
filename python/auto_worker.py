#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys

try:
    from ._ccz_core import _core_auto, _core_equivalence
except ImportError:  # pragma: no cover - direct script execution
    from _ccz_core import _core_auto, _core_equivalence


def main() -> int:
    payload = json.loads(sys.stdin.read())
    mode = payload["mode"]
    truth_table = payload.get("truth_table")
    truth_table_f = payload.get("truth_table_f")
    truth_table_g = payload.get("truth_table_g")
    n_bits = int(payload["n_bits"])
    m_bits = int(payload["m_bits"])
    raw_time_limit_seconds = payload.get("time_limit_seconds")
    if raw_time_limit_seconds is None:
        time_limit_seconds = None
    else:
        time_limit_seconds = float(raw_time_limit_seconds)
    min_active_hyperplanes = payload["min_active_hyperplanes"]
    cpu_affinity = payload.get("cpu_affinity")
    auto_group = payload.get("auto_group")

    if isinstance(cpu_affinity, list) and hasattr(os, "sched_setaffinity"):
        try:
            os.sched_setaffinity(0, {int(core) for core in cpu_affinity})
        except Exception:
            pass

    try:
        if mode == "ccz":
            result = _core_auto(
                mode,
                truth_table,
                n_bits,
                m_bits,
                time_limit_seconds,
                min_active_hyperplanes,
            )
        elif mode == "ccz_equivalence":
            result = _core_equivalence(
                "ccz",
                truth_table_f,
                truth_table_g,
                n_bits,
                m_bits,
                time_limit_seconds,
                min_active_hyperplanes,
                [] if auto_group is None else auto_group,
            )
        elif mode == "ea":
            result = _core_auto(
                mode,
                truth_table,
                n_bits,
                m_bits,
                time_limit_seconds,
                min_active_hyperplanes,
            )
        elif mode == "ea_equivalence":
            result = _core_equivalence(
                "ea",
                truth_table_f,
                truth_table_g,
                n_bits,
                m_bits,
                time_limit_seconds,
                min_active_hyperplanes,
                [] if auto_group is None else auto_group,
            )
        else:
            raise ValueError(f"Unsupported mode: {mode}")
    except Exception as exc:
        sys.stdout.write(
            json.dumps(
                {
                    "ok": False,
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                }
            )
        )
        return 1

    sys.stdout.write(json.dumps({"ok": True, "result": result}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
