"""2. Sigma5 geometry relations."""
import math

import numpy as np

from generator.config import THETA, THETA_DEG, TILT_AXIS
from generator.geometry import (cubic_rotation_group, disorientation_angle_deg,
                                orientation_matrices, rotation_angle_axis)


def test_theta_value():
    assert abs(THETA_DEG - 36.86989764584402) < 1e-9
    assert abs(THETA - 2.0 * math.atan(1.0 / 3.0)) < 1e-15


def test_misorientation_angle_and_axis():
    M1, M2 = orientation_matrices(THETA)
    R = M2 @ M1.T
    angle, axis = rotation_angle_axis(R)
    assert abs(angle - THETA_DEG) < 1e-9
    assert abs(axis[2]) > 0.9999
    assert abs(axis[0]) < 1e-9 and abs(axis[1]) < 1e-9
    assert np.allclose(axis, [0.0, 0.0, 1.0]) or np.allclose(axis, [0.0, 0.0, -1.0])


def test_cubic_symmetry_disorientation():
    M1, M2 = orientation_matrices(THETA)
    R = M2 @ M1.T
    group = cubic_rotation_group()
    assert len(group) == 24
    dmin = disorientation_angle_deg(R, group)
    assert abs(dmin - THETA_DEG) < 1e-9


def test_gb_plane_is_310():
    M1, M2 = orientation_matrices(THETA)
    # GB normal in the lab frame is x; in crystal frame 1 it must be [3,-1,0]
    n1 = M1.T @ np.array([1.0, 0.0, 0.0])
    n2 = M2.T @ np.array([1.0, 0.0, 0.0])
    ideal1 = np.array([3.0, -1.0, 0.0]) / math.sqrt(10.0)
    ideal2 = np.array([3.0, 1.0, 0.0]) / math.sqrt(10.0)
    assert abs(np.dot(n1, ideal1)) > 0.999999
    assert abs(np.dot(n2, ideal2)) > 0.999999


def test_orientation_matrices_are_rotations():
    M1, M2 = orientation_matrices(THETA)
    for M in (M1, M2):
        assert np.allclose(M @ M.T, np.eye(3), atol=1e-12)
        assert abs(np.linalg.det(M) - 1.0) < 1e-12
