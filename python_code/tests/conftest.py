"""Общие фикстуры pytest."""

from __future__ import annotations

import numpy as np
import pytest

from mag_nav import (
    GaussianNoise,
    RotatingDipole,
    Scenario,
    SensorConfig,
    Solenoid,
    grid3d,
    random_first_octant,
)

B0_REF = 1.0e-7  # Тл*м^3  ~  m0 = 0.796 А*м^2


@pytest.fixture(scope="session")
def b0_ref() -> float:
    return B0_REF


@pytest.fixture(scope="session")
def dipole(b0_ref: float) -> RotatingDipole:
    """Эталонный точечный диполь в начале координат, ось вращения Z."""
    return RotatingDipole(
        origin=np.zeros(3),
        moment_amplitude=b0_ref * 4.0 * np.pi / (4.0e-7 * np.pi),
        b0=b0_ref,
        quadrant_policy="off",
    )


@pytest.fixture(scope="session")
def solenoid(b0_ref: float) -> Solenoid:
    """Конечный соленоид с той же дипольной константой, что и ``dipole``."""
    # m = N*I*pi*Rs^2 = b0*4*pi/mu0  ->  подберём N*I
    n_turns = 120
    radius = 0.05
    m_target = b0_ref * 4.0 * np.pi / (4.0e-7 * np.pi)
    current = m_target / (n_turns * np.pi * radius**2)
    return Solenoid(
        origin=np.zeros(3),
        n_turns=n_turns,
        current=current,
        radius=radius,
        length=0.15,
        n_quad=128,
        quadrant_policy="off",
    )


@pytest.fixture(scope="session")
def octant_points() -> np.ndarray:
    """2000 точек, равномерно распределённых в первом октанте, r ∈ [0.1, 1.5] м."""
    return random_first_octant(2000, 0.1, 1.5, seed=20240607)


@pytest.fixture(scope="session")
def coarse_grid() -> np.ndarray:
    return grid3d((0.15, 1.2), (0.15, 1.2), (0.15, 1.2), (5, 5, 5))


@pytest.fixture(scope="session")
def clean_scenario(dipole: RotatingDipole, octant_points: np.ndarray) -> Scenario:
    return Scenario(source=dipole, trajectory=octant_points, inject_at="none", seed=1)


@pytest.fixture(scope="session")
def noisy_ab_scenario(dipole: RotatingDipole, octant_points: np.ndarray) -> Scenario:
    return Scenario(
        source=dipole,
        trajectory=octant_points,
        noise=GaussianNoise(sigma=1e-3, relative=True),
        inject_at="ab",
        seed=7,
    )


@pytest.fixture(scope="session")
def timeseries_scenario(dipole: RotatingDipole) -> Scenario:
    traj = random_first_octant(120, 0.25, 1.0, seed=99)
    return Scenario(
        source=dipole,
        trajectory=traj,
        sensor=SensorConfig(
            samples_per_period=48,
            n_periods=4,
            noise=GaussianNoise(sigma=1e-4, relative=True),
            demod_method="lstsq",
        ),
        inject_at="samples",
        seed=13,
    )