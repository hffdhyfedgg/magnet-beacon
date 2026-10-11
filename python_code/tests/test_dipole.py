from __future__ import annotations

import numpy as np
import pytest

from mag_nav import MU_0, RotatingDipole, dipole_field, norm
from tests.primitives import dipole_checks as chk


def test_inverse_cube_law_x_direction(dipole):
    chk.check_inverse_cube_law(dipole, np.array([1.0, 0.0, 0.0]), np.geomspace(0.1, 5.0, 12))


def test_inverse_cube_law_diagonal(dipole):
    chk.check_inverse_cube_law(dipole, np.array([0.6, 0.5, 0.62]), np.geomspace(0.1, 5.0, 12))


def test_analytic_on_axis_and_equator():
    chk.check_dipole_analytic_on_axis(b0=1e-7, r=0.5)


def test_divergence_free():
    chk.check_dipole_divergence_free()


def test_ab_quadrature_identity():
    chk.check_ab_quadrature_identity()


def test_ab_harmonic_purity():
    chk.check_ab_harmonic_purity()


def test_ab_unit_vector_invariants():
    chk.check_ab_unit_vectors()


def test_field_line_tangency():
    chk.check_field_line_is_tangent()


def test_moment_and_b0_roundtrip():
    src = RotatingDipole(origin=np.zeros(3), moment_amplitude=2.5, quadrant_policy="off")
    assert np.isclose(src.dipole_constant, MU_0 * 2.5 / (4 * np.pi))
    assert np.isclose(src.moment_amplitude_value, 2.5)


def test_rotation_phase_periodicity():
    src = RotatingDipole(origin=np.zeros(3), b0=1e-7, quadrant_policy="off")
    p = np.array([[0.3, 0.4, 0.5]])
    assert np.allclose(src.field_at_phase(p, 0.0), src.field_at_phase(p, 2 * np.pi))
    assert np.allclose(src.field_at_phase(p, 0.3), src.field_at_phase(p, 0.3 + 2 * np.pi))


def test_quadrant_policy_strict_raises():
    src = RotatingDipole(origin=np.zeros(3), b0=1e-7, quadrant_policy="strict")
    with pytest.raises(ValueError):
        # origin == 0 -> relative == absolute; z < 0 нарушает квадрант
        src.field_at_phase(np.array([[0.3, 0.3, -0.2]]), 0.0)


def test_origin_offset():
    src = RotatingDipole(origin=np.array([0.4, 0.0, 0.25]), b0=1e-7, quadrant_policy="off")
    p_abs = np.array([[0.4 + 0.3, 0.4, 0.25 + 0.5]])
    p_rel = np.array([[0.3, 0.4, 0.5]])
    assert np.allclose(
        src.field_at_phase(p_abs, 0.0), dipole_field(p_rel, np.array([1.0, 0.0, 0.0]), 1e-7)
    )


def test_time_series_shape(dipole):
    pts = np.abs(np.random.default_rng(0).uniform(0.2, 1.0, (7, 3)))
    ts = dipole.time_series(pts, np.linspace(0, 2 * np.pi, 33, endpoint=False))
    assert ts.shape == (33, 7, 3)


def test_superposition_linearity():
    """Поле диполя линейно по вектору момента."""
    b0 = 1e-7
    p = np.array([[0.3, 0.4, 0.5]])
    m1, m2 = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    lhs = dipole_field(p, m1 + m2, b0)
    rhs = dipole_field(p, m1, b0) + dipole_field(p, m2, b0)
    assert np.allclose(lhs, rhs, rtol=1e-12)