"""Проверка, что визуализация не падает и не производит побочных эффектов показа."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: до импорта mag_nav.sim.visualize

import matplotlib.pyplot as plt
import numpy as np
import pytest

from mag_nav import Scenario, Solenoid, monte_carlo, random_first_octant
from mag_nav.sim import visualize as viz


@pytest.fixture(scope="module")
def result(dipole):
    traj = random_first_octant(120, 0.25, 1.2, seed=77)
    return Scenario(source=dipole, trajectory=traj, inject_at="none").run()


@pytest.fixture(scope="module")
def noisy_result(dipole):
    from mag_nav import GaussianNoise

    traj = random_first_octant(300, 0.25, 1.2, seed=78)
    sc = Scenario(source=dipole, trajectory=traj, noise=GaussianNoise(sigma=2e-2, relative=True),
                  inject_at="ab", seed=79)
    return sc.run()


def test_plot_field_magnitude(dipole):
    fig = viz.plot_field_magnitude(dipole, plane="xz", n=40)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_field_slice(dipole):
    fig = viz.plot_field_slice(dipole, plane="xz", n=40)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_field_lines(dipole):
    seeds = np.array([[0.0, 0.0, 0.35], [0.25, 0.0, 0.25]])
    fig = viz.plot_field_lines(dipole, seeds=seeds, plane="xz", step=8e-3, n_steps=120)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_trajectory_3d(result):
    fig = viz.plot_trajectory_3d(result)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_error_components(noisy_result):
    fig = viz.plot_error_components(noisy_result)
    assert len(fig.axes) == 3
    plt.close(fig)


def test_plot_error_histogram(noisy_result):
    fig = viz.plot_error_histogram(noisy_result, bins=30)
    assert len(fig.axes) == 4
    plt.close(fig)


def test_plot_error_vs_range(noisy_result):
    fig = viz.plot_error_vs_range(noisy_result)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_ab_invariants(noisy_result):
    fig = viz.plot_ab_invariants(noisy_result)
    assert len(fig.axes) == 3
    plt.close(fig)


def test_plot_polarization(dipole):
    fig = viz.plot_polarization(dipole, points=np.zeros((1, 3)), n=30)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_model_discrepancy(solenoid, dipole):
    fig = viz.plot_model_discrepancy(solenoid, dipole, n=(25, 25),
                                     u_range=(0.15, 1.0), v_range=(0.15, 1.0))
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_monte_carlo(dipole):
    traj = random_first_octant(40, 0.3, 1.0, seed=80)
    from mag_nav import GaussianNoise

    sc = Scenario(source=dipole, trajectory=traj,
                  noise=GaussianNoise(sigma=1e-2, relative=True), inject_at="ab", seed=81)
    mc = monte_carlo(sc, n_trials=3)
    fig = viz.plot_monte_carlo(mc)
    assert len(fig.axes) == 2
    plt.close(fig)


def test_animate_navigation(result):
    anim = viz.animate_navigation(result, interval_ms=10)
    assert anim is not None
    plt.close(anim._fig)


def test_save_figures(tmp_path, result):
    fig = viz.plot_trajectory_3d(result)
    paths = viz.save_figures({"traj": fig}, str(tmp_path), dpi=60)
    assert len(paths) == 1
    assert (tmp_path / "traj.png").exists()
    plt.close(fig)