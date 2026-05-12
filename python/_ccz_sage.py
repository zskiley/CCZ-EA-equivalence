"""Sage object conversion for CCZ/EA Python wrappers."""

from __future__ import annotations

from typing import Any, Optional


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

    if isinstance(candidate, dict) and "generators" in candidate:
        return {"generators": candidate.get("generators") or []}

    if hasattr(candidate, "gens"):
        ambient_dimension = n_bits + m_bits
        return [
            _sage_matrix_to_affine_generator(generator, ambient_dimension)
            for generator in candidate.gens()
        ]

    return auto_group


def _sage_constructors():
    try:
        from sage.all import GF, MatrixGroup, identity_matrix, matrix
    except Exception as exc:  # pragma: no cover
        raise ImportError(
            "ccz_auto, ea_auto, ccz_equivalence, and ea_equivalence now "
            "return Sage matrix objects. "
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


def _affine_result_to_sage_matrix(
    affine_result: Optional[dict[str, Any]],
    n_bits: int,
    m_bits: int,
    sage_constructors: tuple[Any, Any, Any, Any],
    *,
    invert: bool = False,
):
    if affine_result is None:
        return None
    GF, _MatrixGroup, _identity_matrix, matrix_ctor = sage_constructors
    matrix_obj = _affine_generator_to_sage_matrix(
        affine_result, n_bits + m_bits, GF(2), matrix_ctor
    )
    if invert:
        matrix_obj = matrix_obj.inverse()
    return matrix_obj
