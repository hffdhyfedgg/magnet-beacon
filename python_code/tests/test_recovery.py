from __future__ import annotations

import numpy as np
import pytest

from mag_nav import (
    ConstantBias,
    GaussianNoise,
    demodulate,
    demodulate_lstsq,
    demodulation_error,
    make_phase_grid,
    synthesize_series,
)


def _truth(n_pts=8, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n_pts, 3))
    b = rng.standard_normal((n_pts, 3))
    return a, b


def test_projection_exact_integer_periods():
    a, b = _truth()
    ph = make_phase_grid(64, 4)
    s = synthesize_series(a, b, ph)
    res = demodulate(s, ph)
    assert np.allclose(res.a, a, atol=1e-12)
    assert np.allclose(res.b, b, atol=1e-12)


def test_lstsq_exact_integer_periods():
    a, b = _truth()
    ph = make_phase_grid(64, 4)
    res = demodulate_lstsq(synthesize_series(a, b, ph), ph)
    assert np.allclose(res.a, a, atol=1e-11)
    assert np.allclose(res.b, b, atol=1e-11)
    assert np.allclose(res.dc, 0.0, atol=1e-11)


def test_lstsq_rejects_dc_offset():
    """Проекция смещена при нецелом числе периодов; МНК с dc — нет."""
    a, b = _truth()
    spp, nper = 50, 1
    ph = make_phase_grid(spp, nper) + np.linspace(0, 0.7, spp * nper)  # «рваная» фаза
    s = synthesize_series(a, b, ph) + 0.35

    lstsq = demodulate_lstsq(s, ph, with_dc=True)
    proj = demodulate(s, ph)

    err_l = demodulation_error(lstsq, a, b)
    err_p = demodulation_error(proj, a, b)
    assert err_l["rms_abs"] < err_p["rms_abs"]
    assert np.allclose(lstsq.dc, 0.35, atol=1e-2)


def test_underdetermined_raises():
    with pytest.raises(ValueError):
        make_phase_grid(2, 1)
    with pytest.raises(ValueError):
        make_phase_grid(16, 0)


def test_length_mismatch_raises():
    a, b = _truth()
    ph = make_phase_grid(32, 1)
    s = synthesize_series(a, b, ph)
    with pytest.raises(ValueError):
        demodulate(s, ph[:-1])
    with pytest.raises(ValueError):
        demodulate_lstsq(s, ph[:-1])


def test_noise_averaging_scales_with_sqrt_n():
    """СКП оценки a падает ~ sqrt(N) с ростом числа отсчётов."""
    a, b = _truth(n_pts=1)
    a = a[0][None, :]
    b = b[0][None, :]
    sig = 0.1

    def rms_for(spp, nper, trials=200):
        ph = make_phase_grid(spp, nper)
        errs = []
        for k in range(trials):
            rng = np.random.default_rng(k)
            s = synthesize_series(a, b, ph) + GaussianNoise(sigma=sig).sample(
                (ph.size, 1, 3), rng
            )
            r = demodulate_lstsq(s, ph)
            errs.append(np.linalg.norm(r.a - a))
        return float(np.sqrt(np.mean(np.square(errs))))

    e1 = rms_for(16, 1)
    e2 = rms_for(16, 16)
    assert e2 < e1 / 2.5, (e1, e2)


def test_demodulation_error_keys():
    a, b = _truth(n_pts=3)
    ph = make_phase_grid(32, 2)
    res = demodulate(synthesize_series(a, b, ph), ph)
    err = demodulation_error(res, a, b)
    assert set(err) >= {"abs_a", "abs_b", "rel_a", "rel_b", "rms_abs", "max_rel"}
    assert err["rms_abs"] < 1e-11