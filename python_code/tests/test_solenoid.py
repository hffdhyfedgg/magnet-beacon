from __future__ import annotations

import numpy as np
import pytest

from mag_nav import MU_0, Solenoid
from tests.primitives import dipole_checks as chk


def test_on_axis_matches_analytic():
    chk.check_solenoid_on_axis_vs_analytic()


def test_loop_far_field_is_dipole():
    chk.check_loop_far_field_dipole()


def test_convergence_to_dipole():
    dev = chk.check_solenoid_dipole_convergence()
    assert dev[-1] < 1e-3, dev


def test_mirror_symmetry():
    chk.check_solenoid_symmetry()


def test_quadrature_order_convergence():
    """Увеличение ``n_quad`` сходится к пределу (порядок квадратуры)."""
    base = Solenoid(
        origin=np.zeros(3), n_turns=100, current=2.0, radius=0.05,
        length=0.2, n_quad=256, quadrant_policy="off",
    )
    p = np.array([[0.35, 0.22, 0.41]])
    ref = base.field_at_phase(p, 0.0)
    for nq, tol in ((16, 1e-4), (32, 1e-7), (64, 1e-10), (96, 1e-12)):
        sol = Solenoid(
            origin=base.origin, n_turns=base.n_turns, current=base.current,
            radius=base.radius, length=base.length, n_quad=nq, quadrant_policy="off",
        )
        err = np.linalg.norm(sol.field_at_phase(p, 0.0) - ref) / np.linalg.norm(ref)
        assert err < tol, (nq, err)


def test_dipole_constant_formula():
    sol = Solenoid(n_turns=200, current=1.5, radius=0.03, length=0.1, quadrant_policy="off")
    m = 200 * 1.5 * np.pi * 0.03**2
    assert np.isclose(sol.magnetic_moment, m)
    assert np.isclose(sol.dipole_constant, MU_0 * m / (4 * np.pi))


def test_zero_length_limit_is_loop():
    """``L -> 0`` даёт поле одиночного витка."""
    from mag_nav.source.solenoid import loop_field_cylindrical

    a, I, N = 0.03, 2.0, 50
    sol = Solenoid(
        origin=np.zeros(3), n_turns=N, current=I, radius=a,
        length=1e-9, n_quad=8, quadrant_policy="off",
    )
    zeta = np.array([0.0, 0.05, 0.2, 0.7])
    got = sol.field_local(np.stack([np.zeros_like(zeta), np.zeros_like(zeta), zeta], -1))[:, 2]
    br, bz = loop_field_cylindrical(np.zeros_like(zeta), zeta, a, current=N * I)
    assert np.allclose(got, bz, rtol=1e-8, atol=1e-20)


def test_invalid_geometry_raises():
    with pytest.raises(ValueError):
        Solenoid(radius=-0.01)
    with pytest.raises(ValueError):
        Solenoid(length=0.0)


def test_field_is_finite_outside_winding():
    sol = Solenoid(n_turns=80, current=1.0, radius=0.05, length=0.15,
                   n_quad=64, quadrant_policy="off")
    pts = np.array([[0.3, 0.3, 0.3], [0.05, 0.4, 0.0], [1.0, 1.0, 1.0]])
    b = sol.field_at_phase(pts, 0.37)
    assert np.all(np.isfinite(b))


def test_phase_dependence_is_periodic():
    sol = Solenoid(n_turns=80, current=1.0, radius=0.05, length=0.15,
                   n_quad=64, quadrant_policy="off")
    p = np.array([[0.4, 0.3, 0.35]])
    assert np.allclose(sol.field_at_phase(p, 0.1), sol.field_at_phase(p, 0.1 + 2 * np.pi))


def test_not_quadrature_exact():
    sol = Solenoid(quadrant_policy="off")
    assert sol.is_quadrature_exact is False
    assert sol.as_dipole().is_quadrature_exact is True