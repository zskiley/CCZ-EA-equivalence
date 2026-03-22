#!/usr/bin/env python3
from __future__ import annotations

import json
import sys

import ccz


def main() -> int:
    payload = json.loads(sys.stdin.read())
    mode = payload["mode"]
    truth_table = payload["truth_table"]
    n_bits = int(payload["n_bits"])
    m_bits = int(payload["m_bits"])
    time_limit_seconds = float(payload["time_limit_seconds"])
    min_active_hyperplanes = payload["min_active_hyperplanes"]

    try:
        if mode == "ccz":
            result = ccz._core.ccz_auto(  # pylint: disable=protected-access
                truth_table,
                n_bits,
                m_bits,
                time_limit_seconds,
                min_active_hyperplanes,
            )
        elif mode == "ea":
            result = ccz._core.ea_auto(  # pylint: disable=protected-access
                truth_table,
                n_bits,
                m_bits,
                time_limit_seconds,
                min_active_hyperplanes,
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
