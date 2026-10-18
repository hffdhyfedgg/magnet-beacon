from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mag_nav import GaussianNoise, SensorConfig, Solenoid, line
from mag_nav.signal.recovery import demodulate_lstsq, make_phase_grid, synthesize_series
from mag_nav.source.solenoid import gauss_legendre, loop_field_cylindrical
from mag_nav.sim import visualize as viz
from mag_nav.sim.monte_carlo_sim import MonteCarloSimulator


def _old_field_local(solenoid: Solenoid, points: np.ndarray) -> np.ndarray:
    """Reference implementation using the former sequential quadrature loop."""
    p = np.asarray(points, dtype=np.float64)
    px, py, pz = p[..., 0], p[..., 1], p[..., 2]
    rho = np.hypot(px, py)
    cos_phi = np.divide(px, rho, out=np.ones_like(px), where=rho > 0.0)
    sin_phi = np.divide(py, rho, out=np.zeros_like(py), where=rho > 0.0)
    nodes, weights = gauss_legendre(solenoid.n_quad)
    br_sum = np.zeros_like(rho)
    bz_sum = np.zeros_like(pz)
    for node, weight in zip(nodes, weights):
        br, bz = loop_field_cylindrical(
            rho, pz - 0.5 * solenoid.length * node, solenoid.radius, current=1.0
        )
        br_sum += weight * br
        bz_sum += weight * bz
    factor = 0.5 * solenoid.n_turns * solenoid.current
    return np.stack([br_sum * factor * cos_phi, br_sum * factor * sin_phi, bz_sum * factor], axis=-1)


def test_vectorized_solenoid_quadrature_matches_reference():
    sol = Solenoid(
        n_turns=80, current=0.75, radius=0.04, length=0.12,
        n_quad=48, quadrant_policy="off",
    )
    rng = np.random.default_rng(123)
    points = rng.uniform(0.08, 0.8, size=(3, 17, 3))
    assert np.array_equal(sol.field_local(points), _old_field_local(sol, points))


def test_synthesize_series_accepts_precomputed_trigonometry():
    rng = np.random.default_rng(3)
    a, b = rng.normal(size=(10, 3)), rng.normal(size=(10, 3))
    phases = make_phase_grid(24, 2, phase_offset=0.17)
    expected = synthesize_series(a, b, phases)
    got = synthesize_series(a, b, phases, cos_sin=(np.cos(phases), np.sin(phases)))
    assert np.array_equal(got, expected)


def _simulator(n_points: int = 24) -> MonteCarloSimulator:
    source = Solenoid(
        n_turns=50, current=0.5, radius=0.04, length=0.08,
        n_quad=24, quadrant_policy="off",
    )
    trajectory = line(np.array([1.0, 0.8, 0.6]), np.array([0.2, 0.16, 0.12]), n_points)
    noise = GaussianNoise(sigma=1e-9)
    sensor = SensorConfig(
        samples_per_period=16, n_periods=2, noise=noise,
        demod_method="lstsq", with_dc=True,
    )
    return MonteCarloSimulator(
        source=source, trajectory=trajectory, sensor=sensor, noise=noise,
        inject_at="samples", n_harmonic=16,
    )


def test_serial_and_process_monte_carlo_are_reproducible():
    sim = _simulator()
    serial = sim.run(3, base_seed=29, keep_results=False, n_workers=1)
    parallel = sim.run(3, base_seed=29, keep_results=False, n_workers=2)
    assert np.array_equal(serial.valid, parallel.valid)
    assert np.allclose(serial.errors, parallel.errors, rtol=1e-12, atol=1e-14, equal_nan=True)
    assert np.allclose(serial.l2_per_point, parallel.l2_per_point, equal_nan=True)


def test_monte_carlo_result_visualizations_without_full_results():
    sim = _simulator(n_points=20)
    mc = sim.run(3, base_seed=4, keep_results=False)
    figs = [
        viz.plot_monte_carlo(mc),
        viz.plot_error_histogram(mc, bins=12),
        viz.plot_monte_carlo_error_vs_range(mc, sim.trajectory, origin=sim.source.origin),
    ]
    try:
        assert [len(fig.axes) for fig in figs] == [2, 4, 1]
    finally:
        for fig in figs:
            plt.close(fig)
