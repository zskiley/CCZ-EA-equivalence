"""Input normalization for CCZ/EA Python wrappers."""

from __future__ import annotations

from typing import Any, Iterable, Optional


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


def _normalize_function_with_dimensions(
    function_input: Any, n_bits: int, m_bits: int
) -> tuple[list[int], int, int]:
    if isinstance(function_input, tuple) and len(function_input) == 3:
        return _normalize_function_input(function_input)

    _validate_supported_input(function_input, "function_input")
    field = _field_from_callable(function_input)
    if _is_galois_poly(function_input):
        field = function_input.field
    return _to_truth_table(function_input, n_bits, m_bits, field=field), n_bits, m_bits


def _normalize_function_pair(
    left_input: Any, right_input: Any
) -> tuple[list[int], list[int], int, int]:
    left_is_explicit = isinstance(left_input, tuple) and len(left_input) == 3
    right_is_explicit = isinstance(right_input, tuple) and len(right_input) == 3

    if left_is_explicit:
        left_tt, n, m = _normalize_function_input(left_input)
        right_tt, right_n, right_m = _normalize_function_with_dimensions(
            right_input, n, m
        )
    elif right_is_explicit:
        right_tt, n, m = _normalize_function_input(right_input)
        left_tt, right_n, right_m = _normalize_function_with_dimensions(
            left_input, n, m
        )
    else:
        left_tt, n, m = _normalize_function_input(left_input)
        right_tt, right_n, right_m = _normalize_function_input(right_input)

    if (right_n, right_m) != (n, m):
        raise ValueError("function inputs must have the same n and m dimensions")
    return left_tt, right_tt, n, m
