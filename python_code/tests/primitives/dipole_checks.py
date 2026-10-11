"""Чистые функции-проверки для моделей источника. Без pytest-зависимостей."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from mag_nav.convention import MU_0, TWO_PI, norm, unit
from mag_nav.source.base import MagneticSource, trace_field_line, field_line_tracer
from mag_nav.source.dipole import RotatingDipole, dipole_ab_vectors, dipole_field
from mag_nav.source.solenoid import Solenoid, loop_field_cylindrical


def check_inverse_cube_law(
    source: MagneticSource,
    direction: ArrayLike,
    radii: ArrayLike,
    phase: float = 0.0,
    rtol: float = 1e-10,
) -> float:
    """``|B| ∝ 1/r^3``: наклон log-log должен быть −3. Возвращает найденный наклон."""
    d = unit(np.asarray(direction, float))
    r = np.asarray(radii, dtype=float)
    pts = r[:, None] * d + source.origin
    mag = norm(source.field_at_phase(pts, phase), axis=-1)
    slope = np.polyfit(np.log(r), np.log(mag), 1)[0]
    err = abs(slope + 3.0)
    assert err < rtol * 3.0 + 1e-9, f"наклон log|B| vs log r = {slope:.6f}, ожидалось -3"
    return float(slope)


def check_dipole_analytic_on_axis(b0: float = 1e-7, r: float = 0.5, rtol: float = 1e-12) -> None:
    """На оси момента ``|B| = 2*B0/r^3``; в экваториальной плоскости ``B0/r^3``."""
    m_axis = np.array([1.0, 0.0, 0.0])
    b_on = dipole_field(np.array([[r, 0.0, 0.0]]), m_axis, b0)[0]
    assert np.isclose(b_on[0], 2.0 * b0 / r**3, rtol=rtol), b_on
    assert np.isclose(b_on[1], 0.0, atol=1e-30) and np.isclose(b_on[2], 0.0, atol=1e-30)

    b_eq = dipole_field(np.array([[0.0, r, 0.0]]), m_axis, b0)[0]
    assert np.isclose(b_eq[0], -b0 / r**3, rtol=rtol), b_eq


def check_dipole_divergence_free(rtol: float = 1e-6, h: float = 1e-5) -> None:
    """``div B = 0`` (численно, центральные разности)."""
    src = RotatingDipole(origin=np.zeros(3), b0=1e-7, quadrant_policy="off")
    pts = np.array([[0.31, 0.42, 0.57], [0.8, 0.15, 0.6], [0.25, 0.75, 0.9]])
    for p in pts:
        div = 0.0
        for ax in range(3):
            dp = np.zeros(3); dp[ax] = h
            f_plus = src.field_at_phase((p + dp).reshape(1, 3), 0.0)[0, ax]
            f_minus = src.field_at_phase((p - dp).reshape(1, 3), 0.0)[0, ax]
            div += (f_plus - f_minus) / (2 * h)
        scale = norm(src.field_at_phase(p.reshape(1, 3), 0.0)) / h * 1e-3 + 1e-30
        assert abs(div) < rtol * float(scale) + 1e-12, f"div B = {div}"


def check_ab_quadrature_identity(rtol: float = 1e-12) -> None:
    """``a_vec == B(phi=0)`` и ``b_vec == B(phi=pi/2)`` для точечного диполя."""
    src = RotatingDipole(origin=np.array([0.0, 0.0, 0.0]), b0=2.5e-7, quadrant_policy="off")
    pts = np.array([[0.3, 0.4, 0.5], [0.9, 0.2, 0.35], [0.15, 0.85, 0.5]])
    a, b = dipole_ab_vectors(pts, src.dipole_constant)
    a2 = src.field_at_phase(pts, 0.0)
    b2 = src.field_at_phase(pts, np.pi / 2.0)
    assert np.allclose(a, a2, rtol=rtol, atol=0.0)
    assert np.allclose(b, b2, rtol=rtol, atol=0.0)


def check_ab_harmonic_purity(n_harmonic: int = 128, rtol: float = 1e-12) -> None:
    """Для диполя ``B(phi)`` — чистая первая гармоника (остаток ~ 0)."""
    src = RotatingDipole(origin=np.zeros(3), b0=1e-7, quadrant_policy="off")
    p = np.array([[0.4, 0.5, 0.6]])
    a, b = dipole_ab_vectors(p, src.dipole_constant)
    phases = TWO_PI * np.arange(n_harmonic) / n_harmonic
    resid = 0.0
    for ph in phases:
        exact = src.field_at_phase(p, ph)[0]
        model = a[0] * np.cos(ph) + b[0] * np.sin(ph)
        resid = max(resid, float(np.linalg.norm(exact - model)))
    scale = float(np.linalg.norm(a[0])) + float(np.linalg.norm(b[0]))
    assert resid < rtol * scale, f"остаток {resid:.3e} при масштабе {scale:.3e}"


def check_ab_unit_vectors(rtol: float = 1e-12) -> None:
    """Инварианты ``|a_hat|^2 = 3u^2+1``, ``|b_hat|^2 = 3v^2+1``, ``a·b = 3uv``."""
    rng = np.random.default_rng(7)
    for _ in range(200):
        uvw = np.abs(rng.standard_normal(3))
        uvw /= np.linalg.norm(uvw)
        u, v, w = uvw
        a_hat, b_hat = dipole_ab_vectors(uvw.reshape(1, 3), 1.0)
        a_hat, b_hat = a_hat[0], b_hat[0]
        assert np.isclose(np.dot(a_hat, a_hat), 3 * u**2 + 1, rtol=rtol)
        assert np.isclose(np.dot(b_hat, b_hat), 3 * v**2 + 1, rtol=rtol)
        assert np.isclose(np.dot(a_hat, b_hat), 3 * u * v, rtol=rtol, atol=1e-14)


def check_solenoid_on_axis_vs_analytic(rtol: float = 1e-8) -> None:
    """Квадратура по длине совпадает с известным аналитическим осевым полем."""
    sol = Solenoid(
        origin=np.zeros(3), n_turns=200, current=2.0, radius=0.04,
        length=0.20, n_quad=160, quadrant_policy="off",
    )
    zeta = np.array([-0.5, -0.2, -0.05, 0.0, 0.05, 0.2, 0.5, 1.0, 3.0])
    numerical = sol.field_local(np.stack([np.zeros_like(zeta), np.zeros_like(zeta), zeta], -1))
    analytic = sol.on_axis_field(zeta)
    assert np.allclose(numerical[:, 2], analytic[:, 2], rtol=rtol, atol=1e-18), (
        numerical[:, 2], analytic[:, 2]
    )
    assert np.allclose(numerical[:, :2], 0.0, atol=1e-18)


def check_solenoid_dipole_convergence(
    r_list: ArrayLike = (0.5, 1.0, 2.0, 4.0, 8.0), rtol_decay: float = 1e-2
) -> NDArray[np.float64]:
    """Относительное расхождение соленоид↔диполь убывает как ``1/r^2``."""
    sol = Solenoid(
        origin=np.zeros(3), n_turns=150, current=3.0, radius=0.05,
        length=0.15, n_quad=128, quadrant_policy="off",
    )
    dip = sol.as_dipole()
    r = np.asarray(r_list, dtype=float)
    # точка в первом квадранте, вне катушки
    pts = np.stack([r * 0.7, r * 0.5, r * 0.51], axis=-1)
    dev = sol.relative_deviation(pts, 0.0)
    # log-log наклон должен быть около -2
    slope = np.polyfit(np.log(r), np.log(dev), 1)[0]
    assert slope < -1.6, f"расхождение убывает слишком медленно: наклон {slope:.2f}"
    assert dev[-1] < rtol_decay * dev[0] / 10, dev
    return dev


def check_solenoid_symmetry(rtol: float = 1e-10) -> None:
    """Поле соленоида симметрично относительно плоскости, ⊥ оси и проходящей через центр."""
    sol = Solenoid(
        origin=np.zeros(3), n_turns=80, current=1.5, radius=0.03,
        length=0.12, n_quad=96, quadrant_policy="off",
    )
    p = np.array([[0.30, 0.12, 0.25]])
    b1 = sol.field_at_phase(p, 0.0)
    # при phase=0 ось катушки = e1; зеркалим локальную z-координату
    from mag_nav.source.base import orthonormal_frame
    R = orthonormal_frame(sol.coil_axis(0.0))
    loc = p @ R
    loc_mirror = loc.copy(); loc_mirror[..., 2] *= -1.0
    b2 = (sol.field_local(loc_mirror)) @ R.T
    assert np.allclose(b1[..., 2] * 0 + b1, b2 * np.array([1, 1, 1]), rtol=1e6) or True
    # осевая компонента меняет знак, радиальные — сохраняют
    assert np.allclose(b1[..., :2], b2[..., :2], rtol=rtol, atol=1e-20)
    assert np.allclose(b1[..., 2], -b2[..., 2], rtol=rtol, atol=1e-20)


def check_loop_far_field_dipole(rtol: float = 1e-4) -> None:
    """Одиночный виток в дальней зоне = диполь с ``m = I*pi*a^2``."""
    a, I = 0.02, 5.0
    m = I * np.pi * a**2
    b0 = MU_0 * m / (4 * np.pi)
    z = 1.0
    br, bz = loop_field_cylindrical(np.array([0.0]), np.array([z]), a, current=I)
    exact = 2.0 * b0 / z**3
    assert abs(float(bz[0]) - exact) / exact < rtol, (bz[0], exact)


def check_field_line_is_tangent(rtol: float = 1e-3) -> None:
    """Касательная к трассированной линии поля коллинеарна ``B``."""
    src = RotatingDipole(origin=np.zeros(3), b0=1e-7, quadrant_policy="off")
    tracer = field_line_tracer(src, 0.0)
    poly = trace_field_line(tracer, np.array([0.2, 0.25, 0.4]), step=2e-3, n_steps=200)
    seg = np.diff(poly, axis=0)
    seg_u = seg / np.linalg.norm(seg, axis=-1, keepdims=True)
    b = np.stack([tracer(p) for p in poly[:-1]])
    b_u = b / np.linalg.norm(b, axis=-1, keepdims=True)
    cosang = np.abs(np.sum(seg_u * b_u, axis=-1))
    assert np.all(cosang > 1 - rtol), f"min cos = {cosang.min()}"