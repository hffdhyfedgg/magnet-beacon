from __future__ import annotations

import numpy as np
import pytest

from mag_nav import RotationMatrix, solve_orientation, solve_position
from mag_nav.navigation.position import REASON_DENOM, REASON_NONFINITE
from tests.primitives import nav_checks as chk


def test_roundtrip_reference_points():
    pts = np.array([[0.3, 0.4, 0.5], [1.0, 0.2, 0.8], [0.05, 1.2, 0.3]])
    chk.check_roundtrip(pts, 1e-7)


def test_roundtrip_random():
    chk.check_roundtrip_random(n=3000, seed=123)


def test_uvw_normalization(octant_points, dipole):
    from mag_nav import ab_from_dipole

    a, b = ab_from_dipole(octant_points, dipole.dipole_constant).as_tuple()
    sol = solve_position(a, b, dipole.dipole_constant)
    assert np.allclose(sol.uvw_norm, 1.0, rtol=1e-11)


def test_direction_independent_of_calibration():
    chk.check_direction_calibration_independence()


def test_range_calibration_scaling():
    chk.check_range_calibration_scaling()


def test_scalar_helpers_consistency():
    chk.check_scalar_helpers_consistency()


def test_noise_monotonicity():
    out = chk.check_noise_error_monotonic()
    assert out[0] == pytest.approx(0.0, abs=1e-12)


def test_invalid_detection():
    chk.check_invalid_detection()


def test_error_metrics_identity():
    chk.check_error_metrics_identity()


def test_negative_y_recovery():
    chk.check_negative_y_recovery()


def test_vectorization_matches_loop():
    chk.check_vectorization_matches_loop()


def test_batch_shapes_broadcast():
    from mag_nav import ab_from_dipole

    pts = np.abs(np.random.default_rng(3).uniform(0.1, 1.0, (4, 5, 3)))
    a, b = ab_from_dipole(pts, 1e-7).as_tuple()
    sol = solve_position(a, b, 1e-7)
    assert sol.points.shape == (4, 5, 3)
    assert sol.r.shape == (4, 5)
    assert sol.valid.shape == (4, 5)
    assert sol.valid.all()


def test_b0_must_be_positive():
    with pytest.raises(ValueError):
        solve_position(np.zeros((1, 3)), np.zeros((1, 3)), 0.0)
    with pytest.raises(ValueError):
        solve_position(np.zeros((1, 3)), np.zeros((1, 3)), -1e-7)


def test_orientation_not_implemented():
    with pytest.raises(NotImplementedError):
        solve_orientation(np.zeros((1, 3)), np.zeros((1, 3)), np.zeros((1, 3)))


def test_rotation_matrix_roundtrip():
    R = RotationMatrix.from_euler(0.21, -0.34, 1.11)
    R.validate()
    v = np.array([[0.3, 0.4, 0.5], [1.0, 0.0, 0.2]])
    assert np.allclose(R.apply_inverse(R.apply(v)), v, atol=1e-13)
    r, p, y = R.to_euler()
    R2 = RotationMatrix.from_euler(r, p, y)
    assert np.allclose(R.matrix, R2.matrix, atol=1e-12)


def test_rotation_matrix_composition():
    R1 = RotationMatrix.from_euler(0.3, 0.0, 0.0)
    R2 = RotationMatrix.from_euler(0.0, 0.4, 0.0)
    assert np.allclose((R1 @ R2).matrix, R1.matrix @ R2.matrix)
    assert np.allclose((R1 @ R1.inverse()).matrix, np.eye(3), atol=1e-14)


def test_align_vectors():
    src = np.array([1.0, 0.0, 0.0])
    dst = np.array([0.0, 0.6, 0.8])
    R = RotationMatrix.align_vectors(src, dst)
    assert np.allclose(R.apply(src), dst / np.linalg.norm(dst), atol=1e-12)


def test_position_solution_accessors():
    from mag_nav import ab_from_dipole

    pts = np.abs(np.random.default_rng(5).uniform(0.2, 1.0, (10, 3)))
    a, b = ab_from_dipole(pts, 1e-7).as_tuple()
    sol = solve_position(a, b, 1e-7)
    assert np.allclose(sol.x(), pts[:, 0])
    assert np.allclose(sol.y(), pts[:, 1])
    assert np.allclose(sol.z(), pts[:, 2])
    assert sol.n_valid == 10
    assert sol.validity_ratio == 1.0
    assert sol.valid_only().shape == (10, 3)
    assert sol.where_invalid().shape == (0, 3)