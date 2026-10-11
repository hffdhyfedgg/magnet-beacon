#!/usr/bin/env python3
"""Демонстрация: полный прогон симулятора навигации по магнитному маяку.

Запуск::

    uv run python examples/demo.py            # или: python examples/demo.py
    python examples/demo.py --no-save         # без записи PNG
    python examples/demo.py --animate         # + анимация (matplotlib window)

Все фигуры пишутся в ``examples/out/``.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib

if "--no-window" in sys.argv or os.environ.get("MPLBACKEND") is None:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mag_nav import (  # noqa: E402
    ConstantBias,
    GaussianNoise,
    RotatingDipole,
    Scenario,
    SensorConfig,
    Solenoid,
    SumNoise,
    boundary_radius_map,
    calibration_sweep,
    harmonic_residual,
    monte_carlo,
    radial_rays,
    random_first_octant,
)
from mag_nav.sim import visualize as viz  # noqa: E402

MU0 = 4e-7 * np.pi


def make_beacons():
    """Диполь и соленоид с одинаковой дипольной константой."""
    n_turns, radius, length = 150, 0.05, 0.18
    current = 3.0
    sol = Solenoid(
        origin=np.zeros(3),
        n_turns=n_turns,
        current=current,
        radius=radius,
        length=length,
        n_quad=128,
        quadrant_policy="off",
    )
    dip = RotatingDipole(
        origin=np.zeros(3),
        b0=sol.dipole_constant,
        quadrant_policy="off",
    )
    return dip, sol


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "out"))
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--animate", action="store_true")
    ap.add_argument("--quick", action="store_true", help="меньше точек — быстрее")
    args = ap.parse_args(argv)

    viz.apply_style()
    dip, sol = make_beacons()
    n_pts = 300 if args.quick else 1500

    print("=" * 78)
    print("МАГНИТНАЯ НАВИГАЦИЯ ПО ВРАЩАЮЩЕМУСЯ МАЯКУ — демо")
    print("=" * 78)
    print(f"Дипольная константа B0 = {dip.dipole_constant:.6e} Тл·м³")
    print(f"Магнитный момент  m0 = {sol.magnetic_moment:.6e} А·м²")
    print(f"Соленоид: N={sol.n_turns}, I={sol.current} А, R={sol.radius} м, L={sol.length} м")

    figures: dict[str, plt.Figure] = {}

    # ------------------------------------------------------------------ 1
    print("\n[1] Поле маяка")
    figures["01_field_slice"] = viz.plot_field_slice(
        dip, plane="xz", fixed=0.0, extent=(-0.8, 0.8, 0.0, 1.6), n=180
    )
    figures["02_field_lines"] = viz.plot_field_lines(
        dip,
        seeds=np.array([[0.0, 0.0, z] for z in np.linspace(0.15, 0.6, 5)]),
        plane="xz",
        step=5e-3,
        n_steps=700,
        extent=(-0.8, 0.8, -0.8, 0.8),
    )
    figures["03_polarization"] = viz.plot_polarization(
        dip,
        points=np.zeros((1, 3)),
        plane="xz",
        fixed=0.0,
        extent=(-0.8, 0.8, 0.05, 1.6),
        n=120,
    )

    # ------------------------------------------------------------------ 2
    print("[2] Идеальный режим (без шума)")
    traj = random_first_octant(n_pts, 0.15, 1.6, seed=2024)
    clean = Scenario(source=dip, trajectory=traj, inject_at="none").run()
    print("   ", clean.summary())
    figures["04_trajectory_clean"] = viz.plot_trajectory_3d(clean)

    # ------------------------------------------------------------------ 3
    print("[3] Шум в a/b-параметрах (быстрый режим)")
    noisy_ab = Scenario(
        source=dip,
        trajectory=traj,
        noise=GaussianNoise(sigma=2e-2, relative=True),
        inject_at="ab",
        seed=11,
    ).run()
    print("   ", noisy_ab.summary())
    figures["05_trajectory_noisy"] = viz.plot_trajectory_3d(noisy_ab)
    figures["06_error_components"] = viz.plot_error_components(noisy_ab)
    figures["07_error_hist"] = viz.plot_error_histogram(noisy_ab, bins=70, log=True)
    figures["08_error_vs_range"] = viz.plot_error_vs_range(noisy_ab)
    figures["09_ab_invariants"] = viz.plot_ab_invariants(noisy_ab)

    # ------------------------------------------------------------------ 4
    print("[4] Реалистичный тракт: отсчёты + шум + МНК-демодуляция")
    sensor_noise = SumNoise(
        [
            GaussianNoise(sigma=3e-3, relative=True),
            ConstantBias(bias=np.array([2e-9, -3e-9, 1e-9])),
        ]
    )
    ts = Scenario(
        source=dip,
        trajectory=random_first_octant(400, 0.2, 1.4, seed=55),
        sensor=SensorConfig(
            samples_per_period=48,
            n_periods=6,
            noise=sensor_noise,
            demod_method="lstsq",
            with_dc=True,
        ),
        inject_at="samples",
        seed=12,
    ).run()
    print("   ", ts.summary())
    figures["10_timeseries_errors"] = viz.plot_error_vs_range(ts)

    # ------------------------------------------------------------------ 5
    print("[5] Монте-Карло")
    mc_sc = Scenario(
        source=dip,
        trajectory=random_first_octant(200, 0.2, 1.4, seed=56),
        noise=GaussianNoise(sigma=1.5e-2, relative=True),
        inject_at="ab",
        seed=1,
    )
    mc = monte_carlo(mc_sc, n_trials=25, base_seed=1000)
    print("   ", mc.summary())
    figures["11_monte_carlo"] = viz.plot_monte_carlo(mc)

    # ------------------------------------------------------------------ 6
    print("[6] Ошибка калибровки B0")
    cal = calibration_sweep(
        Scenario(
            source=dip, trajectory=random_first_octant(150, 0.25, 1.3, seed=57), inject_at="none"
        ),
        np.geomspace(0.5, 2.0, 15),
    )
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    axes[0].loglog(cal["factors"], cal["rms_position"], "o-")
    axes[0].set_xlabel("B0_исп / B0_ист")
    axes[0].set_ylabel("RMS |Δr|, м")
    axes[0].set_title("Погрешность позиции")
    axes[1].semilogx(cal["factors"], cal["rms_range_relative"], "o-", color="tab:red")
    axes[1].set_xlabel("B0_исп / B0_ист")
    axes[1].set_ylabel("RMS(Δr/r)")
    axes[1].set_title("Погрешность дальности  (~ f^⅓ − 1)")
    axes[2].semilogx(cal["factors"], cal["angular_error_deg"], "o-", color="tab:green")
    axes[2].set_xlabel("B0_исп / B0_ист")
    axes[2].set_ylabel("угол, °")
    axes[2].set_title("Угловая погрешность (не зависит от B0)")
    for ax in axes:
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    figures["12_calibration"] = fig
    print(
        f"    RMS(Δr/r) при f=2: {cal['rms_range_relative'][-1]:.4f} "
        f"(теория {2 ** (1 / 3) - 1:.4f})"
    )

    # ------------------------------------------------------------------ 7
    print("[7] Расхождение моделей: соленоид → диполь")
    figures["13_discrepancy_map"] = viz.plot_model_discrepancy(
        sol,
        dip,
        plane="xz",
        fixed=0.03,
        u_range=(0.12, 1.2),
        v_range=(0.12, 1.2),
        n=(80, 80),
        relative=True,
        log=True,
        component="norm",
    )
    figures["14_boundary_radius"] = viz.plot_boundary_radius(
        sol,
        dip,
        thresholds=(0.01, 0.05, 0.1),
        elevations=(15.0, 30.0, 45.0, 60.0, 75.0),
        azimuth=45.0,
        r_range=(0.1, 2.5),
        n_r=200,
    )
    for el in (15.0, 45.0, 75.0):
        br = boundary_radius_map(
            sol,
            dip,
            threshold=0.02,
            elevation_deg=el,
            azimuth_deg=45.0,
            r_grid=np.geomspace(0.08, 3.0, 300),
        )
        print(f"    el={el:5.1f}°  →  r*(2%) = {float(br['boundary_r'][0]):.3f} м")
        _, rel = harmonic_residual(
            sol,
            radial_rays(
                1,
                np.array([1.0]),
                elevation_deg=el,
            ),
            n_harmonic=96,
        )  # azimuth_deg=45.0
        print(f"               несинусоидальность на r=1 м: {float(rel.max()):.2e}")

    # ------------------------------------------------------------------ 8
    if args.animate:
        anim = viz.animate_navigation(clean, interval_ms=30, trail=80, source=dip)
        anim._fig.savefig(
            os.path.join(args.out, "anim_preview.png"), dpi=110
        ) if not args.no_save else None
        plt.show()

    if not args.no_save:
        os.makedirs(args.out, exist_ok=True)
        paths = viz.save_figures(figures, args.out, dpi=140)
        print(f"\n[✓] Сохранено {len(paths)} фигур в {args.out}")
    else:
        for f in figures.values():
            plt.close(f)

    print("\nГотово.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
