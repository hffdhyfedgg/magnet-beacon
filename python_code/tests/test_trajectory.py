from __future__ import annotations

import numpy as np
import pytest

from mag_nav.sim import trajectory as trj


def test_line_endpoints():
    p = trj.line([0, 0, 0], [1, 2, 3], 11)
    assert p.shape == (11, 3)
    assert np.allclose(p[0], [0, 0, 0])
    assert np.allclose(p[-1], [1, 2, 3])
    d = np.linalg.norm(np.diff(p, axis=0), axis=-1)
    assert np.allclose(d, d[0])


def test_line_no_endpoint():
    p = trj.line([0, 0, 0], [1, 0, 0], 4, endpoint=False)
    assert np.allclose(p[:, 0], [0, 0.25, 0.5, 0.75])


def test_line_too_few_points():
    with pytest.raises(ValueError):
        trj.line([0, 0, 0], [1, 1, 1], 1)


def test_waypoint_path_continuity():
    wp = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]]
    p = trj.waypoint_path(wp, n_per_segment=20)
    assert p.shape[1] == 3
    assert np.allclose(p[0], [0, 0, 0])
    assert np.allclose(p[-1], [1, 1, 1])


def test_waypoint_path_closed():
    wp = [[0, 0, 0], [1, 0, 0], [1, 1, 0]]
    p = trj.waypoint_path(wp, n_per_segment=10, closed=True)
    assert np.allclose(p[0], p[-1], atol=1e-12) or np.allclose(p[-1], [0, 0, 0])


def test_arc_circle_closed():
    c = trj.circle([0, 0, 0], 0.5, 200, plane="xz")
    r = np.linalg.norm(c, axis=-1)
    assert np.allclose(r, 0.5)
    assert np.allclose(c[0], c[-1], atol=1e-12)


def test_helix_radius_and_pitch():
    h = trj.helix([0, 0, 0], [0, 0, 1], radius=0.3, pitch=0.2, n_turns=3, n=600)
    rho = np.hypot(h[:, 0], h[:, 1])
    assert np.allclose(rho, 0.3, atol=1e-12)
    assert np.isclose(h[-1, 2] - h[0, 2], 3 * 0.2, atol=1e-9)


def test_grid3d_count_and_exclusion():
    g = trj.grid3d((0, 1), (0, 1), (0, 1), (4, 5, 6))
    assert g.shape == (120, 3)
    g2 = trj.grid3d((0, 1), (0, 1), (0, 1), (4, 5, 6), exclude_radius=0.3)
    assert g2.shape[0] < 120
    assert np.all(np.linalg.norm(g2, axis=-1) >= 0.3 - 1e-12)


def test_plane_grid_labels_and_shapes():
    for plane, labels in (("xy", ("x, м", "y, м")), ("xz", ("x, м", "z, м")), ("yz", ("y, м", "z, м"))):
        pts, U, V, lab = trj.plane_grid(plane, 0.25, (0, 1), (0, 1), (6, 7))
        assert pts.shape == (6, 7, 3)
        assert U.shape == (6, 7)
        assert lab == labels
        idx = {"xy": 2, "xz": 1, "yz": 0}[plane]
        assert np.allclose(pts[..., idx], 0.25)


def test_plane_grid_bad_plane():
    with pytest.raises(ValueError):
        trj.plane_grid("ab", 0.0, (0, 1), (0, 1), (3, 3))


def test_random_first_octant_bounds_and_reproducibility():
    p1 = trj.random_first_octant(5000, 0.2, 1.0, seed=7)
    p2 = trj.random_first_octant(5000, 0.2, 1.0, seed=7)
    assert np.array_equal(p1, p2)
    r = np.linalg.norm(p1, axis=-1)
    assert np.all((r >= 0.2 - 1e-12) & (r <= 1.0 + 1e-12))
    assert np.all(p1 >= -1e-15)
    # равномерность по объёму: медиана r ближе к r_max, чем середина отрезка
    assert np.median(r) > 0.5 * (0.2 + 1.0)


def test_random_first_octant_bad_range():
    with pytest.raises(ValueError):
        trj.random_first_octant(10, 1.0, 0.5)


def test_radial_rays_shape():
    p = trj.radial_rays(5, np.geomspace(0.1, 1.0, 7), elevation_deg=30)
    assert p.shape == (35, 3)
    assert np.all(p >= -1e-15)


def test_uniform_moving_observer():
    t = np.array([0.0, 1.0, 2.0])
    p = trj.uniform_moving_observer([0, 0, 0.5], [0.1, 0.2, 0.0], t)
    assert np.allclose(p, [[0, 0, 0.5], [0.1, 0.2, 0.5], [0.2, 0.4, 0.5]])