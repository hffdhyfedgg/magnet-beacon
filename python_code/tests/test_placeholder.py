"""Смоук-тест: пакет импортируется, публичный API на месте."""

from __future__ import annotations

import mag_nav


def test_version():
    assert mag_nav.__version__ == "0.1.0"


def test_public_api_present():
    expected = {
        "RotatingDipole", "Solenoid", "solve_position", "Scenario", "SensorConfig",
        "GaussianNoise", "ab_from_source", "monte_carlo", "RotationMatrix", "visualize",
    }
    missing = expected - set(mag_nav.__all__)
    assert not missing, f"В __all__ отсутствуют: {missing}"
    for name in expected:
        assert hasattr(mag_nav, name), name


def test_submodules_importable():
    import importlib

    for mod in (
        "mag_nav.convention",
        "mag_nav.source.base", "mag_nav.source.dipole", "mag_nav.source.solenoid",
        "mag_nav.signal.ab_params", "mag_nav.signal.recovery",
        "mag_nav.navigation.position", "mag_nav.navigation.orientation",
        "mag_nav.noise.base", "mag_nav.noise.gaussian", "mag_nav.noise.composite",
        "mag_nav.sim.trajectory", "mag_nav.sim.scenario", "mag_nav.sim.visualize",
    ):
        importlib.import_module(mod)