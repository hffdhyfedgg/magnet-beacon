from __future__ import annotations

import numpy as np
import pytest

from mag_nav import (
    ConstantBias,
    GaussianNoise,
    RotatingDipole,
    Scenario,
    SensorConfig,
    Solenoid,
    SumNoise,
    calibration_sweep,
    monte_carlo,
    random_first_octant,
    radial_rays,
)
from mag_nav.sim.scenario import boundary_radius_map, model_discrepancy_map


# --------------------------------------------------------------------------- #
# Идеальный режим
# --------------------------------------------------------------------------- #
def test_clean_run_is_exact(clean_scenario):
    res = clean_scenario.run()
    assert res.valid.all()
    assert np.allclose(res.error_abs, 0.0, atol=1e-9)
    assert res.statistics()["rms_norm"] < 1e-9
    assert res.inject_at == "none"


def test_clean_run_solenoid_has_systematic_error(solenoid):
    """Для соленоида обратное решение (дипольная алгебра) даёт смещение ≠ 0."""
    traj = random_first_octant(400, 0.25, 1.5, seed=3)
    res = Scenario(source=solenoid, trajectory=traj, inject_at="none").run()
    assert res.valid.all()
    rms = res.statistics()["rms_norm"]
    assert rms > 0.0
    # в дальней зоне смещение мало
    far = traj[np.linalg.norm(traj, axis=-1) > 1.2]
    res_far = Scenario(source=solenoid, trajectory=far, inject_at="none").run()
    assert res_far.statistics()["rms_rel"] < rms / max(np.linalg.norm(far, axis=-1).min(), 1e-9) * 5 + 0.05


def test_reproducibility(noisy_ab_scenario):
    r1 = noisy_ab_scenario.run()
    r2 = noisy_ab_scenario.run()
    assert np.allclose(r1.estimated_positions, r2.estimated_positions)


def test_seed_override_changes_result(noisy_ab_scenario):
    r1 = noisy_ab_scenario.run()
    r2 = noisy_ab_scenario.run(seed=noisy_ab_scenario.seed + 1)
    assert not np.allclose(r1.estimated_positions, r2.estimated_positions)


# --------------------------------------------------------------------------- #
# Шум
# --------------------------------------------------------------------------- #
def test_ab_noise_degrades_accuracy(clean_scenario, noisy_ab_scenario):
    assert clean_scenario.run().statistics()["rms_norm"] < noisy_ab_scenario.run().statistics()["rms_norm"]


def test_timeseries_mode_averages_noise(dipole):
    """Усреднение за много периодов должно снижать СКП примерно как sqrt(N)."""
    traj = random_first_octant(60, 0.3, 1.0, seed=21)
    base = dict(source=dipole, trajectory=traj, inject_at="samples")

    def rms(spp, nper):
        sc = Scenario(
            **base,
            sensor=SensorConfig(
                samples_per_period=spp,
                n_periods=nper,
                noise=GaussianNoise(sigma=5e-3, relative=True),
            ),
            seed=5,
        )
        return sc.run().statistics()["rms_norm"]

    short, long_ = rms(32, 1), rms(32, 64)
    assert long_ < short / 3.0, (short, long_)


def test_timeseries_mode_with_dc_bias(dipole):
    """МНК с dc-членом компенсирует постоянное внешнее поле."""
    traj = random_first_octant(80, 0.3, 1.0, seed=22)
    noise = SumNoise([GaussianNoise(sigma=1e-4, relative=True), ConstantBias(bias=np.array([1e-8, -2e-8, 5e-9]))])

    def rms(with_dc):
        sc = Scenario(
            source=dipole,
            trajectory=traj,
            inject_at="samples",
            sensor=SensorConfig(samples_per_period=48, n_periods=4, noise=noise, with_dc=with_dc),
            seed=3,
        )
        return sc.run().statistics()["rms_norm"]

    assert rms(True) < rms(False)


def test_projection_vs_lstsq_both_work(timeseries_scenario):
    sc_proj = Scenario(
        source=timeseries_scenario.source,
        trajectory=timeseries_scenario.trajectory,
        sensor=SensorConfig(
            samples_per_period=48, n_periods=4,
            noise=GaussianNoise(sigma=1e-4, relative=True),
            demod_method="projection",
        ),
        inject_at="samples",
        seed=13,
    )
    assert sc_proj.run().statistics()["rms_norm"] > 0


def test_noise_level_monotonic(dipole):
    traj = random_first_octant(150, 0.3, 1.2, seed=31)
    prev = -1.0
    for s in (0.0, 1e-3, 1e-2, 5e-2):
        sc = Scenario(
            source=dipole, trajectory=traj,
            noise=GaussianNoise(sigma=s, relative=True),
            inject_at="ab" if s > 0 else "none", seed=2,
        )
        rms = sc.run().statistics()["rms_norm"]
        assert rms >= prev - 1e-12
        prev = rms


# --------------------------------------------------------------------------- #
# Монте-Карло
# --------------------------------------------------------------------------- #
def test_monte_carlo_aggregation(noisy_ab_scenario):
    mc = monte_carlo(noisy_ab_scenario, n_trials=8, base_seed=100)
    assert mc.n_trials == 8
    assert mc.error_norm_mean.shape == (noisy_ab_scenario.trajectory.shape[0],)
    assert mc.per_component_rms.shape == (3,)
    assert 0.0 <= mc.validity_ratio <= 1.0
    assert mc.overall_rms > 0
    assert mc.overall_p95 >= mc.overall_rms * 0.5
    assert "MonteCarlo" in mc.summary()


def test_monte_carlo_zero_trials(noisy_ab_scenario):
    with pytest.raises(ValueError):
        monte_carlo(noisy_ab_scenario, n_trials=0)


# --------------------------------------------------------------------------- #
# Калибровка
# --------------------------------------------------------------------------- #
def test_calibration_sweep_range_scaling(dipole):
    traj = random_first_octant(80, 0.3, 1.0, seed=41)
    sc = Scenario(source=dipole, trajectory=traj, inject_at="none")
    out = calibration_sweep(sc, np.array([0.5, 1.0, 2.0]))
    # ошибка дальности для f: |f^(1/3) - 1|
    expected = np.abs(np.array([0.5, 1.0, 2.0]) ** (1 / 3) - 1.0)
    assert np.allclose(out["rms_range_relative"], expected, rtol=1e-6)
    # угловая погрешность не зависит от калибровки
    assert np.allclose(out["angular_error_deg"], 0.0, atol=1e-6)
    assert out["rms_position"][1] == pytest.approx(0.0, abs=1e-9)


def test_b0_override(dipole):
    traj = random_first_octant(40, 0.3, 1.0, seed=42)
    sc = Scenario(source=dipole, trajectory=traj, inject_at="none", b0_override=dipole.dipole_constant * 8.0)
    res = sc.run()
    assert np.allclose(res.r_est / res.r_true, 2.0, rtol=1e-10)
    assert np.allclose(
        res.estimated_positions / np.linalg.norm(res.estimated_positions, axis=-1, keepdims=True),
        res.true_positions / np.linalg.norm(res.true_positions, axis=-1, keepdims=True),
        rtol=1e-10,
    )


# --------------------------------------------------------------------------- #
# Расхождение моделей
# --------------------------------------------------------------------------- #
def test_model_discrepancy_map_shapes(solenoid, dipole):
    out = model_discrepancy_map(
        solenoid, dipole, plane="xz", fixed=0.2,
        u_range=(0.1, 1.0), v_range=(0.1, 1.0), n=(25, 30),
    )
    assert out["error_abs"].shape == (25, 30, 3)
    assert out["error_norm"].shape == (25, 30)
    assert out["valid"].shape == (25, 30)
    assert np.all(out["error_norm"][out["valid"]] >= 0)


def test_model_discrepancy_grows_near_beacon(solenoid, dipole):
    near = model_discrepancy_map(
        solenoid, dipole, plane="xz", fixed=0.05,
        u_range=(0.10, 0.25), v_range=(0.10, 0.25), n=(20, 20),
    )
    far = model_discrepancy_map(
        solenoid, dipole, plane="xz", fixed=0.05,
        u_range=(1.20, 1.60), v_range=(1.20, 1.60), n=(20, 20),
    )
    e_near = np.nanmean(near["error_rel"][near["valid"]])
    e_far = np.nanmean(far["error_rel"][far["valid"]])
    assert e_near > e_far * 5, (e_near, e_far)


def test_boundary_radius_map(solenoid, dipole):
    out = boundary_radius_map(
        solenoid, dipole, threshold=0.02,
        elevation_deg=45, azimuth_deg=45,
        r_grid=np.geomspace(0.08, 3.0, 150),
    )
    assert out["r"].size == 150
    br = float(out["boundary_r"][0])
    assert np.isfinite(br) and br > 0
    # за границей погрешность действительно ниже порога
    idx = np.searchsorted(out["r"], br)
    tail = out["metric"][idx:]
    assert np.all(tail[np.isfinite(tail)] < 0.02)


# --------------------------------------------------------------------------- #
# Разное
# --------------------------------------------------------------------------- #
def test_without_noise_helper(noisy_ab_scenario):
    clean = noisy_ab_scenario.without_noise().run()
    assert clean.inject_at == "none"
    assert clean.statistics()["rms_norm"] < 1e-9


def test_with_trajectory_helper(clean_scenario):
    new_traj = radial_rays(3, np.array([0.5, 1.0]))
    res = clean_scenario.with_trajectory(new_traj).run()
    assert res.true_positions.shape == new_traj.shape


def test_frame_alignment_identity_is_noop(dipole):
    from mag_nav.navigation.orientation import IDENTITY, RotationMatrix

    traj = random_first_octant(50, 0.3, 1.0, seed=51)
    base = Scenario(source=dipole, trajectory=traj, inject_at="none").run()
    rot = Scenario(
        source=dipole, trajectory=traj, inject_at="none", frame_alignment=IDENTITY
    ).run()
    assert np.allclose(base.estimated_positions, rot.estimated_positions)


def test_frame_alignment_rotation_invariance(dipole):
    """Поворот рамки приёмника компенсируется и не меняет оценку позиции."""
    from mag_nav.navigation.orientation import RotationMatrix

    traj = random_first_octant(60, 0.3, 1.0, seed=52)
    R = RotationMatrix.from_euler(0.13, -0.22, 0.31)
    base = Scenario(source=dipole, trajectory=traj, inject_at="none").run()
    rot = Scenario(source=dipole, trajectory=traj, inject_at="none", frame_alignment=R).run()
    assert np.allclose(base.estimated_positions, rot.estimated_positions, rtol=1e-10, atol=1e-12)


def test_result_repr_and_summary(clean_scenario):
    res = clean_scenario.run()
    assert "SimResult" in repr(res)
    assert "RMS" in res.summary()
    assert set(res.ab_error()) == {"abs_a", "abs_b", "rel_a", "rel_b"}