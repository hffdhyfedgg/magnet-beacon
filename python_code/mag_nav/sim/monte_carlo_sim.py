"""Универсальный класс симуляции Монте-Карло для навигации по магнитному маяку.

Даёт единый интерфейс для прогона ансамбля реализаций с произвольной
конфигурацией источника, траектории, шума и тракта измерения.

Метрики по ансамблю
--------------------
``l2``
    Корень из среднего квадрата отклонения (RMSE) по всем точкам и реализациям:
    ``sqrt( mean( ||est - true||^2 ) )``.

``l_inf``
    Средний по реализациям максимум модуля отклонения по любой из координат:
    ``mean_over_trials( max_over_points_and_axes( |est - true| ) )``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..convention import as_points, norm, safe_norm
from ..navigation.position import PositionSolution, solve_position
from ..noise.base import NoiseModel
from ..signal.ab_params import ABParams, ab_from_source
from ..signal.recovery import (
    DemodulationResult,
    demodulate,
    demodulate_lstsq,
    make_phase_grid,
    synthesize_series,
)
from ..source.base import MagneticSource
from .scenario import Scenario, SensorConfig, SimResult

__all__ = ["MonteCarloSimulator", "MonteCarloResult"]


# --------------------------------------------------------------------------- #
# Результат
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MonteCarloResult:
    """Агрегированный результат ансамбля Монте-Карло."""

    n_trials: int
    n_points: int

    # Метрики
    l2: float
    """RMSE по всем точкам и реализациям, м."""

    l2_per_point: NDArray[np.float64]
    """RMSE в каждой точке траектории (усреднение по реализациям), (n_points,)."""

    l_inf: float
    """Средний по реализациям максимум |отклонения| по любой координате, м."""

    l_inf_per_trial: NDArray[np.float64]
    """Максимум |отклонения| для каждой реализации, (n_trials,)."""

    l2_relative: float
    """RMSE, нормированная на среднюю дальность."""

    # Сырые данные
    errors: NDArray[np.float64]
    """(n_trials, n_points, 3) — вектор ошибки."""

    valid: NDArray[np.bool_]
    """(n_trials, n_points) — маска корректных решений."""

    per_trial_results: tuple[SimResult, ...]
    """Полные результаты каждого пролёта (для детальной диагностики)."""

    meta: dict[str, Any] = field(default_factory=dict)

    # -- удобство ------------------------------------------------------------ #
    @property
    def validity_ratio(self) -> float:
        return float(self.valid.mean()) if self.valid.size else 0.0

    def summary(self) -> str:
        return (
            f"MonteCarlo[{self.n_trials} trials × {self.n_points} pts]: "
            f"L2={self.l2:.4e} м, L∞={self.l_inf:.4e} м, "
            f"L2_rel={self.l2_relative:.3e}, valid={100 * self.validity_ratio:.1f}%"
        )

    def __repr__(self) -> str:
        return self.summary()


# --------------------------------------------------------------------------- #
# Симулятор
# --------------------------------------------------------------------------- #
@dataclass
class MonteCarloSimulator:
    """Универсальный инструмент для ансамблевых симуляций навигации.

    Parameters
    ----------
    source :
        Модель маяка (:class:`RotatingDipole`, :class:`Solenoid`, …).
    trajectory :
        Точки наблюдения ``(n, 3)`` в мировой системе, м.
    sensor :
        Конфигурация тракта измерения. ``None`` → шум в a/b или отсутствует.
    noise :
        Модель шума. Приоритет у ``sensor.noise``, если задан.
    inject_at :
        ``"samples"`` | ``"ab"`` | ``"none"`` — куда вносить шум.
    b0_override :
        Переопределение дипольной константы для обратного решения
        (исследование ошибки калибровки).
    n_harmonic :
        Число фаз для выделения первой гармоники (для не-дипольных источников).
    """

    source: MagneticSource
    trajectory: NDArray[np.float64]
    sensor: SensorConfig | None = None
    noise: NoiseModel | None = None
    inject_at: Literal["samples", "ab", "none"] = "samples"
    b0_override: float | None = None
    n_harmonic: int = 64

    def __post_init__(self) -> None:
        self.trajectory = as_points(self.trajectory)

    # -- конфигурация -------------------------------------------------------- #
    @property
    def active_noise(self) -> NoiseModel | None:
        if self.sensor is not None and self.sensor.noise is not None:
            return self.sensor.noise
        return self.noise

    @property
    def b0_solve(self) -> float:
        if self.b0_override is not None:
            return float(self.b0_override)
        return self.source.dipole_constant

    def with_noise(self, noise: NoiseModel | None) -> "MonteCarloSimulator":
        """Копия с другой моделью шума."""
        return MonteCarloSimulator(
            source=self.source,
            trajectory=self.trajectory.copy(),
            sensor=self.sensor,
            noise=noise,
            inject_at=self.inject_at,
            b0_override=self.b0_override,
            n_harmonic=self.n_harmonic,
        )

    def with_trajectory(self, trajectory: ArrayLike) -> "MonteCarloSimulator":
        """Копия с другой траекторией."""
        return MonteCarloSimulator(
            source=self.source,
            trajectory=as_points(trajectory),
            sensor=self.sensor,
            noise=self.noise,
            inject_at=self.inject_at,
            b0_override=self.b0_override,
            n_harmonic=self.n_harmonic,
        )

    def with_b0_factor(self, factor: float) -> "MonteCarloSimulator":
        """Копия с масштабированной калибровкой."""
        return MonteCarloSimulator(
            source=self.source,
            trajectory=self.trajectory.copy(),
            sensor=self.sensor,
            noise=self.noise,
            inject_at=self.inject_at,
            b0_override=self.source.dipole_constant * factor,
            n_harmonic=self.n_harmonic,
        )

    # -- один пролёт ---------------------------------------------------------- #
    def _single_run(self, rng: np.random.Generator) -> SimResult:
        """Выполнить один прогон без создания нового ГСЧ извне."""
        pts = self.trajectory
        ab_true = ab_from_source(self.source, pts, n_harmonic=self.n_harmonic)
        a_true, b_true = ab_true.a, ab_true.b
        rel_true = self.source.relative_points(pts)

        noise = self.active_noise
        mode = self.inject_at
        if noise is None:
            mode = "none"

        meta: dict[str, Any] = {
            "mode": mode,
            "source": type(self.source).__name__,
            "b0_true": float(self.source.dipole_constant),
            "b0_solve": self.b0_solve,
            "n_harmonic": self.n_harmonic,
        }

        if mode == "none":
            a_meas, b_meas = a_true.copy(), b_true.copy()
        elif mode == "ab":
            a_meas = noise.apply(a_true, rng)
            b_meas = noise.apply(b_true, rng)
            meta["noise"] = noise.describe()
        elif mode == "samples":
            sensor = self.sensor or SensorConfig()
            phases = sensor.phases()
            series = synthesize_series(a_true, b_true, phases)
            if noise is not None:
                series = noise.apply(series, rng)
                meta["noise"] = noise.describe()
            if sensor.demod_method == "lstsq":
                demod = demodulate_lstsq(series, phases, with_dc=sensor.with_dc)
            else:
                demod = demodulate(series, phases)
            a_meas, b_meas = demod.a, demod.b
            meta["demod_method"] = sensor.demod_method
            meta["n_samples"] = sensor.n_samples
        else:
            raise ValueError(f"Неизвестный inject_at: {mode!r}")

        sol = solve_position(a_meas, b_meas, self.b0_solve)
        est_rel = sol.points
        est_abs = est_rel + self.source.origin

        return SimResult(
            true_positions=rel_true,
            estimated_positions=est_rel,
            true_positions_abs=pts,
            estimated_positions_abs=est_abs,
            a_true=a_true,
            b_true=b_true,
            a_meas=a_meas,
            b_meas=b_meas,
            solution=sol,
            b0=self.b0_solve,
            valid=sol.valid,
            inject_at=mode,
            meta=meta,
        )

    # -- ансамбль ------------------------------------------------------------- #
    def run(
        self,
        n_trials: int,
        *,
        base_seed: int = 0,
        keep_results: bool = True,
    ) -> MonteCarloResult:
        """Прогнать ``n_trials`` реализаций и агрегировать метрики.

        Parameters
        ----------
        n_trials :
            Число реализаций шума.
        base_seed :
            Начальный seed; реализация ``k`` получает ``base_seed + k``.
        keep_results :
            Сохранять ли полные :class:`SimResult` (для детальной диагностики).

        Returns
        -------
        MonteCarloResult
        """
        if n_trials < 1:
            raise ValueError("n_trials >= 1")

        n_pts = self.trajectory.shape[0]
        errors = np.full((n_trials, n_pts, 3), np.nan, dtype=np.float64)
        valid = np.zeros((n_trials, n_pts), dtype=np.bool_)
        results: list[SimResult] = []

        for k in range(n_trials):
            rng = np.random.default_rng(base_seed + k)
            res = self._single_run(rng)
            errors[k] = res.error_abs
            valid[k] = res.valid
            if keep_results:
                results.append(res)

        # --- L2 (RMSE) ---
        # По каждой точке: RMSE по реализациям
        err_norm = np.linalg.norm(errors, axis=-1)  # (n_trials, n_pts)
        err_norm_masked = np.where(valid, err_norm, np.nan)

        with np.errstate(invalid="ignore"):
            l2_per_point = np.sqrt(np.nanmean(err_norm_masked**2, axis=0))
            all_valid_norms = err_norm[valid]
            l2_global = (
                float(np.sqrt(np.mean(all_valid_norms**2)))
                if all_valid_norms.size
                else float("nan")
            )

        # Относительная: нормировка на среднюю дальность
        r_true = norm(self.trajectory - self.source.origin, axis=-1)
        mean_r = float(np.mean(r_true[r_true > 0])) if np.any(r_true > 0) else 1.0
        l2_relative = l2_global / mean_r if mean_r > 0 else float("nan")

        # --- L_inf ---
        # Для каждой реализации: максимум |ошибки| по точкам и осям
        abs_errors = np.abs(errors)  # (n_trials, n_pts, 3)
        abs_errors_masked = np.where(valid[..., np.newaxis], abs_errors, 0.0)
        l_inf_per_trial = np.max(abs_errors_masked, axis=(1, 2))  # (n_trials,)
        l_inf_global = float(np.mean(l_inf_per_trial))

        return MonteCarloResult(
            n_trials=n_trials,
            n_points=n_pts,
            l2=l2_global,
            l2_per_point=l2_per_point,
            l_inf=l_inf_global,
            l_inf_per_trial=l_inf_per_trial,
            l2_relative=l2_relative,
            errors=errors,
            valid=valid,
            per_trial_results=tuple(results) if keep_results else (),
            meta={
                "source": type(self.source).__name__,
                "b0_true": float(self.source.dipole_constant),
                "b0_solve": self.b0_solve,
                "inject_at": self.inject_at,
                "noise": self.active_noise.describe() if self.active_noise else "none",
            },
        )
