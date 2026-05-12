"""Input normalization for CCZ/EA Python wrappers."""

from __future__ import annotations

from typing import Any


def _normalize_inputs(function_input: Any) -> tuple[list[int], int, int]:
    if isinstance(function_input, tuple) and len(function_input) == 3:
        values, n_bits, m_bits = function_input
        n = int(n_bits)
        m = int(m_bits)
        table = [int(v) for v in values]
    elif hasattr(function_input, "parent"):
        parent = function_input.parent()
        field = parent.base_ring()
        if not hasattr(field, "fetch_int"):
            field = parent
        if not hasattr(field, "fetch_int"):
            raise TypeError("polynomial input must be over a finite field")

        n = int(field.degree())
        m = n
        mask = (1 << n) - 1
        table = []
        for x in range(1 << n):
            field_x = field.fetch_int(x)
            y = function_input(field_x)
            table.append(int(y.integer_representation()) & mask)
    else:
        table = [int(v) for v in function_input]
        size = len(table)
        if size <= 0 or (size & (size - 1)) != 0:
            raise ValueError("truth table length must be a positive power of two")
        n = size.bit_length() - 1
        m = n

    if n <= 0 or m <= 0:
        raise ValueError("n_bits and m_bits must be positive")
    if len(table) != (1 << n):
        raise ValueError("truth table length must be 2^n_bits")
    if any(v < 0 or v >= (1 << m) for v in table):
        raise ValueError("truth table values must satisfy 0 <= value < 2^m_bits")
    return table, n, m
