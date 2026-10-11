from __future__ import annotations

import numpy as np
import pytest

from mag_nav import (
    ClipNoise,
    CompositeNoise,
    ConstantBias,
    GaussianNoise,
    NoiseChain,
    QuantizationNoise,
    SumNoise,
)


def _statistics(sigma, n=400_000, seed=0, relative=False):
    rng = np.random.default_rng(seed)
    model = GaussianNoise(sigma=sigma, relative=relative)
    ref = np.ones((n, 3)) * 2.0 if relative else None
    s = model.sample((n, 3), rng, reference=ref)
    return float(s.mean()), float(s.std())


def test_gaussian_zero_mean_unit_sigma():
    mean, std = _statistics(1.0)
    assert abs(mean) < 3e-3
    assert abs(std - 1.0) < 3e-3


def test_gaussian_sigma_scaling():
    _, s1 = _statistics(1.0)
    _, s2 = _statistics(0.25)
    assert abs(s2 / s1 - 0.25) < 0.01


def test_gaussian_per_component_sigma():
    rng = np.random.default_rng(1)
    model = GaussianNoise(sigma=np.array([1.0, 2.0, 4.0]))
    s = model.sample((300_000, 3), rng)
    stds = s.std(axis=0)
    assert np.allclose(stds, [1.0, 2.0, 4.0], rtol=0.02)


def test_gaussian_relative_mode():
    rng = np.random.default_rng(2)
    model = GaussianNoise(sigma=0.05, relative=True)
    ref = np.tile(np.array([3.0, 4.0, 0.0]), (200_000, 1))  # |ref| = 5
    s = model.sample(ref.shape, rng, reference=ref)
    assert abs(s.std() - 0.05 * 5.0) < 0.05


def test_gaussian_relative_requires_reference():
    model = GaussianNoise(sigma=0.05, relative=True)
    with pytest.raises(ValueError):
        model.sample((10, 3), np.random.default_rng(0))


def test_gaussian_cov_correlation():
    cov = np.array([[1.0, 0.8, 0.0], [0.8, 1.0, 0.0], [0.0, 0.0, 0.25]])
    rng = np.random.default_rng(3)
    model = GaussianNoise(cov=cov)
    s = model.sample((400_000, 3), rng)
    emp = np.cov(s.T)
    assert np.allclose(emp, cov, atol=0.02)


def test_clip_sigma_limits_outliers():
    rng = np.random.default_rng(4)
    model = GaussianNoise(sigma=1.0, clip_sigma=2.0)
    s = model.sample((200_000, 3), rng)
    assert np.abs(s).max() <= 2.0 + 1e-12


def test_constant_bias():
    rng = np.random.default_rng(5)
    model = ConstantBias(bias=np.array([1.0, -2.0, 0.5]))
    out = model.apply(np.zeros((100, 3)), rng)
    assert np.allclose(out, [1.0, -2.0, 0.5])
    assert model.effective_sigma == pytest.approx(np.linalg.norm([1.0, -2.0, 0.5]))


def test_sum_noise_variance_addition():
    rng = np.random.default_rng(6)
    m = SumNoise([GaussianNoise(sigma=3.0), GaussianNoise(sigma=4.0)])
    s = m.sample((300_000, 3), rng)
    assert abs(s.std() - 5.0) < 0.05
    assert m.effective_sigma == pytest.approx(5.0)


def test_composite_alias():
    assert CompositeNoise is SumNoise


def test_sum_noise_rejects_non_additive():
    with pytest.raises(TypeError):
        SumNoise([QuantizationNoise(lsb=0.1)])


def test_chain_sequential_application():
    rng = np.random.default_rng(7)
    chain = NoiseChain([GaussianNoise(sigma=0.0), ConstantBias(bias=0.7), QuantizationNoise(lsb=0.5)])
    out = chain.apply(np.zeros((50, 3)), rng)
    assert np.allclose(out, np.round(0.7 / 0.5) * 0.5)


def test_chain_requires_models():
    with pytest.raises(ValueError):
        NoiseChain([])


def test_quantization_noise():
    rng = np.random.default_rng(8)
    q = QuantizationNoise(lsb=0.25)
    out = q.apply(np.array([[0.13, 0.37, -0.61]]), rng)
    assert np.allclose(out, [[0.25, 0.25, -0.5]])
    assert q.effective_sigma == pytest.approx(0.25 / np.sqrt(12))


def test_clip_noise():
    rng = np.random.default_rng(9)
    c = ClipNoise(limit=1.0)
    out = c.apply(np.array([[-5.0, 0.3, 7.0]]), rng)
    assert np.allclose(out, [[-1.0, 0.3, 1.0]])


def test_reproducibility_same_seed():
    model = GaussianNoise(sigma=1.0)
    a = model.sample((1000, 3), np.random.default_rng(42))
    b = model.sample((1000, 3), np.random.default_rng(42))
    assert np.array_equal(a, b)


def test_different_seeds_differ():
    model = GaussianNoise(sigma=1.0)
    a = model.sample((1000, 3), np.random.default_rng(1))
    b = model.sample((1000, 3), np.random.default_rng(2))
    assert not np.array_equal(a, b)


def test_apply_preserves_signal():
    rng = np.random.default_rng(10)
    signal = np.full((100_000, 3), 7.0)
    out = GaussianNoise(sigma=0.5).apply(signal, rng)
    assert abs(out.mean() - 7.0) < 0.01