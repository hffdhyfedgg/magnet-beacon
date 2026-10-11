from __future__ import annotations

import numpy as np
import pytest

from mag_nav import (
    RotatingDipole,
    ab_from_dipole,
    ab_from_source,
    ab_invariants,
    make_phase_grid,
    synthesize_series,
)
from tests.primitives import ab_checks as chk


def test_denominator_equals_two_kappa():
    chk.check_invariant_denominator_equals_2kappa()


def test_ab_from_source_matches_analytic():
    chk.check_ab_from_source_matches_analytic()


def test_first_harmonic_pure_tone():
    chk.check_first_harmonic_recovers_pure_tone()


def test_dipole_residual_zero():
    chk.check_harmonic_residual_dipole_is_zero()


def test_solenoid_residual_positive_near_field():
    chk.check_harmonic_residual_solenoid_positive_in_near_field()


def test_polarization_bounds():
    chk.check_polarization_ellipse_bounds()


def test_solenoid_far_field_ab():
    chk.check_solenoid_ab_close_to_dipole_far_field()


def test_shape_mismatch_raises():
    a = np.zeros((4, 3))
    b = np.zeros((5, 3))
    from mag_nav.signal.ab_params import ABParams

    with pytest.raises(ValueError):
        ABParams(a=a, b=b)


def test_phase_grid_properties():
    ph = make_phase_grid(32, 3, 0.4)
    assert ph.size == 96
    assert np.isclose(ph[0], 0.4)
    assert np.isclose(ph[1] - ph[0], 2 * np.pi / 32)


def test_synthesize_series_roundtrip():
    rng = np.random.default_rng(0)
    a, b = rng.standard_normal((6, 3)), rng.standard_normal((6, 3))
    ph = make_phase_grid(64, 2)
    s = synthesize_series(a, b, ph)
    assert s.shape == (128, 6, 3)
    # проекция обратно
    n = ph.size
    a_r = (2 / n) * np.sum(s * np.cos(ph).reshape(-1, 1, 1), axis=0)
    b_r = (2 / n) * np.sum(s * np.sin(ph).reshape(-1, 1, 1), axis=0)
    assert np.allclose(a_r, a, atol=1e-12)
    assert np.allclose(b_r, b, atol=1e-12)


def test_ab_scales_as_one_over_r3(dipole):
    p1 = np.array([[0.6, 0.5, 0.7]])
    p2 = p1 * 2.0
    a1, b1 = ab_from_dipole(p1, dipole.dipole_constant).as_tuple()
    a2, b2 = ab_from_dipole(p2, dipole.dipole_constant).as_tuple()
    assert np.allclose(np.linalg.norm(a2) / np.linalg.norm(a1), 1 / 8.0, rtol=1e-12)
    assert np.allclose(np.linalg.norm(b2) / np.linalg.norm(b1), 1 / 8.0, rtol=1e-12)