"""Восстановление a/b-параметров из зашумлённой выборки сигнала.

Аналог MATLAB-скриптов серии ``NoizeSimulation*`` и ``MNKerr.m``.

Два метода:

``demodulate``
    Прямая квадратурная проекция (быстро, точно при целом числе периодов
    и равномерной сетке).

``demodulate_lstsq``
    МНК-оценка по модели ``B_i = a*cos(phi_i) + b*sin(phi_i) + c``
    (устойчива к нецелому числу периодов и к постоянному смещению,
    например к фону геомагнитного поля).

Все функции **чистые**: ГСЧ передаётся извне, побочные эффекты отсутствуют.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..convention import TWO_PI, norm, safe_norm

__all__ = [
    "DemodulationResult",
    "make_phase_grid",
    "synthesize_series",
    "demodulate",
    "demodulate_lstsq",
    "estimate_amplitude",
    "demodulation_error",
]


def make_phase_grid(
    samples_per_period: int,
    n_periods: int = 1,
    phase_offset: float = 0.0,
) -> NDArray[np.float64]:
    """Равномерная сетка фаз ``(n_samples,)`` без дублирования последней точки."""
    spp = int(samples_per_period)
    nper = int(n_periods)
    if spp < 4:
        raise ValueError("samples_per_period >= 4 (теорема Котельникова для 1-й гармоники)")
    if nper < 1:
        raise ValueError("n_periods >= 1")
    n = spp * nper
    return TWO_PI * np.arange(n) / spp + phase_offset


def synthesize_series(a_vec: ArrayLike, b_vec: ArrayLike, phases: ArrayLike) -> NDArray[np.float64]:
    """Синтез ``B(t) = a cos(phi) + b sin(phi)``.

    Parameters
    ----------
    a_vec, b_vec : (..., 3)
    phases : (n_samples,)

    Returns
    -------
    (n_samples, ..., 3)
    """
    a = np.asarray(a_vec, dtype=np.float64)
    b = np.asarray(b_vec, dtype=np.float64)
    ph = np.asarray(phases, dtype=np.float64).reshape(-1)
    shape = (ph.size,) + (1,) * (a.ndim - 1) + (1,)
    c = np.cos(ph).reshape(shape)
    s = np.sin(ph).reshape(shape)
    return c * a[np.newaxis, ...] + s * b[np.newaxis, ...]


@dataclass(frozen=True)
class DemodulationResult:
    """Результат демодуляции."""

    a: NDArray[np.float64]
    b: NDArray[np.float64]
    dc: NDArray[np.float64] | None = None
    residual_rms: NDArray[np.float64] | None = None
    method: str = "projection"

    @property
    def shape(self) -> tuple[int, ...]:
        return self.a.shape[:-1]


def demodulate(samples: ArrayLike, phases: ArrayLike) -> DemodulationResult:
    """Квадратурная проекция (дискретное преобразование на известной частоте).

    Parameters
    ----------
    samples : (n_samples, ..., 3)
        Мгновенные отсчёты магнитометра, ось 0 — время.
    phases : (n_samples,)
        Фазы ``w*t_i``, рад.

    Returns
    -------
    DemodulationResult с ``a``, ``b`` формы ``(..., 3)``.
    """
    s = np.asarray(samples, dtype=np.float64)
    ph = np.asarray(phases, dtype=np.float64).reshape(-1)
    n = ph.size
    if s.shape[0] != n:
        raise ValueError(f"samples.shape[0]={s.shape[0]} != len(phases)={n}")

    shape = (n,) + (1,) * (s.ndim - 2) + (1,)
    a = (2.0 / n) * np.sum(s * np.cos(ph).reshape(shape), axis=0)
    b = (2.0 / n) * np.sum(s * np.sin(ph).reshape(shape), axis=0)
    return DemodulationResult(a=a, b=b, method="projection")


def demodulate_lstsq(
    samples: ArrayLike,
    phases: ArrayLike,
    *,
    with_dc: bool = True,
    rcond: float | None = None,
) -> DemodulationResult:
    """МНК-оценка ``(a, b[, dc])`` — аналог ``MNKerr.m``.

    Модель: ``B_i = a*cos(phi_i) + b*sin(phi_i) + c``. Наличие ``c`` делает
    оценку несмещённой при постоянном фоне (геомагнитное поле, offset АЦП).
    """
    s = np.asarray(samples, dtype=np.float64)
    ph = np.asarray(phases, dtype=np.float64).reshape(-1)
    n = ph.size
    if s.shape[0] != n:
        raise ValueError("Длина оси времени не совпадает с phases")

    cols = [np.cos(ph), np.sin(ph)]
    if with_dc:
        cols.append(np.ones(n))
    design = np.stack(cols, axis=1)  # (n, k)

    trailing = s.shape[1:]
    y = s.reshape(n, -1)  # (n, m)
    coef, *_ = np.linalg.lstsq(design, y, rcond=rcond)  # (k, m)

    a = coef[0].reshape(trailing)
    b = coef[1].reshape(trailing)
    dc = coef[2].reshape(trailing) if with_dc else None

    residual = y - design @ coef  # (n_samples, n_points*3)
    residual_3d = residual.reshape((n,) + trailing)  # (n_samples, n_points, 3)
    # RMS модуля векторного остатка, усреднённый по времени → (n_points,)
    rms = np.sqrt(np.mean(np.sum(residual_3d**2, axis=-1), axis=0))
    return DemodulationResult(a=a, b=b, dc=dc, residual_rms=rms, method="lstsq")


def estimate_amplitude(
    samples: ArrayLike, phases: ArrayLike, *, method: str = "lstsq"
) -> DemodulationResult:
    """Единая точка входа: выбрать метод демодуляции по имени."""
    if method == "projection":
        return demodulate(samples, phases)
    if method == "lstsq":
        return demodulate_lstsq(samples, phases)
    raise ValueError(f"Неизвестный метод демодуляции: {method!r}")


def demodulation_error(
    recovered: DemodulationResult, a_true: ArrayLike, b_true: ArrayLike
) -> dict[str, NDArray[np.float64] | float]:
    """Метрики качества восстановления a/b."""
    a_t = np.asarray(a_true, dtype=np.float64)
    b_t = np.asarray(b_true, dtype=np.float64)

    err_a = recovered.a - a_t
    err_b = recovered.b - b_t
    scale = np.maximum(np.sqrt((np.sum(a_t**2, -1) + np.sum(b_t**2, -1)) / 2.0), 1e-300)
    return {
        "abs_a": norm(err_a),
        "abs_b": norm(err_b),
        "rel_a": norm(err_a) / scale,
        "rel_b": norm(err_b) / scale,
        "rms_abs": float(np.sqrt(np.mean(norm(err_a) ** 2 + norm(err_b) ** 2) / 2.0)),
        "max_rel": float(
            np.max(
                np.concatenate(
                    [np.atleast_1d(norm(err_a) / scale), np.atleast_1d(norm(err_b) / scale)]
                )
            )
        ),
    }
