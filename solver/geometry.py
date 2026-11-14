"""휨 곡면과 TIM 두께 지도.

제어점 3x3 (2차 clamped B-spline 곡면). 배열 P[j][i]: j = y 방향 인덱스, i = x 방향 인덱스.

        y=L  [ 0      w_mid_x   0 ]
             [ w_mid_y w_center w_mid_y ]
        y=0  [ 0      w_mid_x   0 ]
              x=0               x=L

여기에 조립 기울기(tilt_x, tilt_y)를 평면으로 더한다.
TIM 두께 = tim_min + (휨의 최댓값 - 휨)   → 가장 높이 솟은 점에서 TIM이 가장 얇다.
"""
import numpy as np
from scipy.interpolate import BSpline

_KNOTS = np.array([0, 0, 0, 1, 1, 1.0])


def grid_coords(L, n):
    """셀 중심 좌표 (1D)와 2D 격자 X, Y. 배열 모양은 (y, x)."""
    x = (np.arange(n) + 0.5) * L / n
    X, Y = np.meshgrid(x, x, indexing="xy")
    return x, X, Y


def basis(L, n):
    """2차 B-spline 기저함수 3개를 격자점에서 계산 → (3, n)."""
    x, _, _ = grid_coords(L, n)
    return np.stack([BSpline(_KNOTS, np.eye(3)[k], 2)(x / L) for k in range(3)])


def control_net(w_center, w_mid_x, w_mid_y):
    """대칭 휨 변수 3개 → 3x3 제어점 높이 [m]."""
    return np.array([[0.0, w_mid_x, 0.0],
                     [w_mid_y, w_center, w_mid_y],
                     [0.0, w_mid_x, 0.0]])


def warp_surface(P, L, n):
    """3x3 제어점 → 휨 곡면 w(x, y) [m], 모양 (n, n), 행 = y."""
    B = basis(L, n)
    return B.T @ np.asarray(P, float) @ B


def tilt_plane(tilt_x, tilt_y, L, n):
    """양 끝 높이 차가 tilt_x, tilt_y인 평면 [m]."""
    x, X, Y = grid_coords(L, n)
    return tilt_x * (X / L - 0.5) + tilt_y * (Y / L - 0.5)


def tim_thickness(w, tim_min):
    """휨 곡면 → TIM 두께 지도. 가장 높은 점의 두께 = tim_min."""
    return tim_min + (w.max() - w)
