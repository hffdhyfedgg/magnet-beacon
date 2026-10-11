"""Проверки прямой задачи a/b-параметров."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from mag_nav.convention import TWO_PI, norm, safe_norm
from mag_nav.signal.ab_params import (
    ab_from_dipole,
    ab_from_source,
    ab_invariants,
    first_harmonic,
    harmonic_residual,
    invariants_from_geometry,
    phase_grid,
    polarization_ellipse,
)
from mag_nav.source.base import MagneticSource
from mag_nav.source.dipole import RotatingDipole, dipole_invariants
from mag_nav.source.solenoid import Solenoid


def check_invariant_denominator_equals_2kappa(rtol: float = 1e-12) -> None:
    """``D = asq + bsq - sqrt((asq-bsq)^2+4ab^2) == 2*kappa``."""
    rng = np.random.default_rng(3)
    b0 = 3.1e-7
    for _ in range(300):
        uvw = np.abs(rng.standard_normal(3))
        uvw /= np.linalg.norm(uvw)
        r = float(rng.uniform(0.1, 2.0))
        p = (r * uvw).reshape(1, 3)
        kappa, asq, bsq, ab = dipole_invariants(p, b0)
        inv = ab_invariants(*ab_from_dipole(p, b0).as_tuple())
        assert np.isclose(float(inv.denominator), 2.0 * float(kappa), rtol=rtol)


def check_ab_from_source_matches_analytic(rtol: float = 1e-12) -> None:
    """Универсальный путь ``ab_from_source`` == аналитика для диполя."""
    src = RotatingDipole(origin=np.zeros(3), b0=2e-7, quadrant_policy="off")
    rng = np.random.default_rng(11)
    pts = np.abs(rng.uniform(0.1, 1.5, size=(25, 3)))
    fast = ab_from_source(src, pts)
    slow = ab_from_source(src, pts, force_harmonic=True, n_harmonic=64)
    ref = ab_from_dipole(src.relative_points(pts), src.dipole_constant)
    assert fast.exact and not slow.exact
    assert np.allclose(fast.a, ref.a, rtol=rtol)
    assert np.allclose(fast.b, ref.b, rtol=rtol)
    assert np.allclose(slow.a, ref.a, rtol=1e-11)
    assert np.allclose(slow.b, ref.b, rtol=1e-11)


def check_first_harmonic_recovers_pure_tone(rtol: float = 1e-12) -> None:
    """``first_harmonic`` точно восстанавливает ``a cos + b sin``."""
    rng = np.random.default_rng(5)
    a = rng.standard_normal((4, 3))
    b = rng.standard_normal((4, 3))
    n = 64
    ph = phase_grid(n)
    samples = np.cos(ph)[:, None, None] * a[None] + np.sin(ph)[:, None, None] * b[None]
    a_r, b_r = first_harmonic(samples, ph)
    assert np.allclose(a_r, a, rtol=rtol, atol=1e-14)
    assert np.allclose(b_r, b, rtol=rtol, atol=1e-14)


def check_harmonic_residual_dipole_is_zero(atol: float = 1e-13) -> None:
    src = RotatingDipole(origin=np.zeros(3), b0=1e-7, quadrant_policy="off")
    pts = np.array([[0.3, 0.4, 0.5], [0.7, 0.2, 0.6]])
    rms, rel = harmonic_residual(src, pts, n_harmonic=64)
    scale = norm(src.field_at_phase(pts, 0.0), axis=-1)
    assert np.all(rms < atol * scale), (rms, scale)
    assert np.all(rel < atol)


def check_harmonic_residual_solenoid_positive_in_near_field() -> None:
    """Соленоид несинусоидален в ближней зоне → остаток заметно больше машинного."""
    sol = Solenoid(
        origin=np.zeros(3), n_turns=120, current=4.0, radius=0.06,
        length=0.18, n_quad=96, quadrant_policy="off",
    )
    near = np.array([[0.10, 0.12, 0.14]])
    far = np.array([[3.0, 3.0, 3.0]])
    _, rel_near = harmonic_residual(sol, near, n_harmonic=96)
    _, rel_far = harmonic_residual(sol, far, n_harmonic=96)
    assert float(rel_near.max()) > 1e-4, rel_near
    assert float(rel_far.max()) < float(rel_near.max()) / 10, (rel_near, rel_far)


def check_polarization_ellipse_bounds(rtol: float = 1e-10) -> None:
    """``s_min <= |B(phi)| <= s_max`` для всех фаз."""
    src = RotatingDipole(origin=np.zeros(3), b0=1e-7, quadrant_policy="off")
    rng = np.random.default_rng(2)
    pts = np.abs(rng.uniform(0.2, 1.2, size=(30, 3)))
    ab = ab_from_source(src, pts)
    s_min, s_max, _ = polarization_ellipse(ab.a, ab.b)
    for ph in np.linspace(0, TWO_PI, 257):
        mag = norm(src.field_at_phase(pts, ph), axis=-1)
        assert np.all(mag <= s_max * (1 + rtol) + 1e-300)
        assert np.all(mag >= s_min * (1 - rtol) - 1e-300)


def check_solenoid_ab_close_to_dipole_far_field(rtol: float = 5e-3) -> None:
    """a/b соленоида сходятся к дипольным в дальней зоне."""
    sol = Solenoid(
        origin=np.zeros(3), n_turns=100, current=2.0, radius=0.04,
        length=0.12, n_quad=128, quadrant_policy="off",
    )
    dip = sol.as_dipole()
    pts = np.array([[2.0, 1.5, 1.8], [3.0, 2.5, 2.0]])
    ab_s = ab_from_source(sol, pts, n_harmonic=96)
    ab_d = ab_from_source(dip, pts)
    da = norm(ab_s.a - ab_d.a, axis=-1) / safe_norm(ab_d.a, axis=-1)
    db = norm(ab_s.b - ab_d.b, axis=-1) / safe_norm(ab_d.b, axis=-1)
    assert np.all(da < rtol), da
    assert np.all(db < rtol), db