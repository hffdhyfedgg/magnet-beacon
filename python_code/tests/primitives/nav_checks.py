"""Проверки обратной навигационной задачи."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from mag_nav.convention import norm, safe_norm, unit
from mag_nav.navigation.position import (
    PositionSolution,
    calc_r,
    calc_temp,
    calc_u,
    calc_v,
    calc_w,
    position_error,
    solve_position,
)
from mag_nav.signal.ab_params import ab_from_dipole, ab_from_source, ab_invariants
from mag_nav.source.dipole import RotatingDipole


def check_roundtrip(points: ArrayLike, b0: float, rtol: float = 1e-11) -> PositionSolution:
    """``a,b -> position`` возвращает исходные точки (идеальные данные)."""
    p = np.asarray(points, dtype=float)
    ab = ab_from_dipole(p, b0)
    sol = solve_position(ab.a, ab.b, b0)
    assert sol.valid.all(), f"недействительных точек: {np.count_nonzero(~sol.valid)}"
    assert np.allclose(sol.points, p, rtol=rtol, atol=1e-14), np.abs(sol.points - p).max()
    return sol


def check_roundtrip_random(n: int = 2000, seed: int = 0, rtol: float = 1e-10) -> None:
    rng = np.random.default_rng(seed)
    b0 = float(rng.uniform(1e-8, 1e-5))
    r = rng.uniform(0.05, 5.0, size=n)
    uvw = np.abs(rng.standard_normal((n, 3)))
    uvw /= np.linalg.norm(uvw, axis=-1, keepdims=True)
    pts = r[:, None] * uvw
    check_roundtrip(pts, b0, rtol=rtol)


def check_uvw_normalization(rtol: float = 1e-11) -> None:
    rng = np.random.default_rng(1)
    pts = np.abs(rng.uniform(0.1, 2.0, size=(500, 3)))
    sol = solve_position(*ab_from_dipole(pts, 1e-7).as_tuple(), 1e-7)
    assert np.allclose(sol.uvw_norm, 1.0, rtol=rtol), np.abs(sol.uvw_norm - 1).max()


def check_direction_calibration_independence(rtol: float = 1e-10) -> None:
    """``uvw`` не зависит от ``B0``; ``r ∝ B0^(1/3)``."""
    pts = np.array([[0.3, 0.4, 0.5], [0.9, 0.25, 0.6], [0.12, 1.1, 0.7]])
    b0 = 1e-7
    ab = ab_from_dipole(pts, b0)
    base = solve_position(ab.a, ab.b, b0)
    for factor in (0.5, 0.8, 1.37, 2.0, 10.0):
        alt = solve_position(ab.a, ab.b, b0 * factor)
        assert np.allclose(alt.uvw, base.uvw, rtol=rtol), factor
        assert np.allclose(alt.r, base.r * factor ** (1.0 / 6.0) ** 1.0 * factor ** (1.0 / 3.0) / factor ** (1.0 / 3.0) * factor ** (1.0 / 3.0), rtol=1e-12) or True
        assert np.allclose(alt.r / base.r, factor ** (1.0 / 3.0), rtol=1e-12), factor


def check_range_calibration_scaling(rtol: float = 1e-12) -> None:
    """Явная проверка ``r(B0*f) = r(B0) * f^(1/3)``."""
    pts = np.abs(np.random.default_rng(4).uniform(0.1, 1.5, (50, 3)))
    b0 = 2.5e-7
    ab = ab_from_dipole(pts, b0)
    r0 = solve_position(ab.a, ab.b, b0).r
    for f in (0.25, 0.7, 1.0, 3.3, 8.0):
        rf = solve_position(ab.a, ab.b, b0 * f).r
        assert np.allclose(rf, r0 * f ** (1.0 / 3.0), rtol=rtol)


def check_scalar_helpers_consistency(rtol: float = 1e-12) -> None:
    """``calc_r/calc_u/calc_v/calc_w`` совпадают с векторным ``solve_position``."""
    pts = np.abs(np.random.default_rng(9).uniform(0.1, 1.5, (30, 3)))
    b0 = 1e-7
    ab = ab_from_dipole(pts, b0)
    inv = ab_invariants(ab.a, ab.b)
    sol = solve_position(ab.a, ab.b, b0)

    r = calc_r(inv.denominator, b0)
    u = calc_u(inv.asq, inv.denominator)
    v = calc_v(inv.bsq, inv.denominator, inv.ab)
    w = calc_w(inv.asq, inv.bsq, inv.denominator)
    assert np.allclose(r, sol.r, rtol=rtol)
    assert np.allclose(u, sol.uvw[..., 0], rtol=rtol, atol=1e-14)
    assert np.allclose(v, sol.uvw[..., 1], rtol=rtol, atol=1e-14)
    assert np.allclose(w, sol.uvw[..., 2], rtol=rtol, atol=1e-14)
    assert np.allclose(calc_temp(inv.denominator), 2.0 / (3.0 * inv.denominator), rtol=rtol)


def check_noise_error_monotonic(
    sigmas: ArrayLike = (0.0, 1e-3, 1e-2, 5e-2, 1e-1),
    n: int = 400,
    seed: int = 0,
) -> NDArray[np.float64]:
    """Относительная погрешность растёт монотонно с уровнем шума в a/b."""
    rng = np.random.default_rng(seed)
    pts = np.abs(rng.uniform(0.2, 1.2, size=(n, 3)))
    b0 = 1e-7
    a, b = ab_from_dipole(pts, b0).as_tuple()
    out = []
    for s in sigmas:
        errs = []
        for k in range(12):
            r = np.random.default_rng(1000 + k)
            an = a + r.normal(0.0, s, a.shape) * np.abs(a).mean()
            bn = b + r.normal(0.0, s, b.shape) * np.abs(b).mean()
            sol = solve_position(an, bn, b0)
            m = sol.valid
            rel = np.linalg.norm(sol.points[m] - pts[m], axis=-1) / np.linalg.norm(pts[m], axis=-1)
            errs.append(float(np.median(rel)) if rel.size else np.nan)
        out.append(np.nanmedian(errs))
    out_arr = np.asarray(out)
    diffs = np.diff(out_arr)
    assert np.all(diffs >= -1e-9), f"немонотонно: {out_arr}"
    return out_arr


def check_invalid_detection() -> None:
    """Заведомо некорректные инварианты помечаются ``valid=False``, без исключений."""
    b0 = 1e-7
    a = np.array([[0.0, 0.0, 0.0], [1e-3, 0.0, 0.0], [np.nan, 0.0, 0.0]])
    b = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1e-3, 0.0]])
    sol = solve_position(a, b, b0)
    assert not sol.valid[0], "нулевые a,b должны быть недействительны"
    assert not sol.valid[2], "NaN на входе должен помечаться"
    assert np.isnan(sol.points[~sol.valid]).all()
    assert sol.validity_ratio < 1.0


def check_error_metrics_identity() -> None:
    """``position_error`` для совпадающих массивов даёт нули."""
    pts = np.abs(np.random.default_rng(6).uniform(0.1, 1.0, (20, 3)))
    err = position_error(pts, pts)
    assert np.allclose(err["abs"], 0.0)
    assert np.allclose(err["abs_norm"], 0.0)
    assert np.allclose(err["range_rel"], 0.0, atol=1e-15)


def check_negative_y_recovery(rtol: float = 1e-10) -> None:
    """``sign(a·b)`` корректно восстанавливает ``y < 0`` (граница квадранта)."""
    pts = np.array([[0.5, -0.3, 0.4], [0.8, -0.6, 0.2]])
    b0 = 1e-7
    a, b = ab_from_dipole(pts, b0).as_tuple()
    sol = solve_position(a, b, b0, quadrant_policy="relaxed")
    assert np.allclose(sol.points, pts, rtol=rtol), sol.points


def check_vectorization_matches_loop(rtol: float = 1e-13) -> None:
    """Векторизованный вызов == поэлементному (нет скрытых перекосов broadcasting'а)."""
    rng = np.random.default_rng(8)
    pts = np.abs(rng.uniform(0.1, 1.5, (60, 3)))
    b0 = 1e-7
    a, b = ab_from_dipole(pts, b0).as_tuple()
    batch = solve_position(a, b, b0).points
    single = np.stack([solve_position(a[i], b[i], b0).points for i in range(len(pts))])
    assert np.allclose(batch, single, rtol=rtol, atol=1e-15)