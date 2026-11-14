"""해석기 검증 테스트.  실행: pytest -q"""
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "solver"))
from solver import load_params, solve, solve_fields, tim_map, heat_flux_map, warp_peak_to_valley  # noqa: E402

P = load_params()
F = P["fixed"]
HOT = [55e-6, 25e-6, 25e-6, 0, 0, 3.0, 22e-3, 4e-3, 170e4]   # 모서리 근처 핫스팟 + 볼록 휨


def flat_x(k_tim=3.0, hs_q=0.0):
    return [0, 0, 0, 0, 0, k_tim, F["die_size"] / 2, F["die_size"] / 2, hs_q]


@pytest.mark.parametrize("nz", [1, 2, 3])
def test_uniform_matches_analytic(nz):
    """휨 없음 + 균일 발열 → 층 수와 상관없이 T = T_c + q * R (1D 해석해와 정확히 일치)."""
    _, T = solve(flat_x(), P, nz=nz)
    R = F["die_thickness"] / F["k_si"] + F["tim_min"] / 3.0 + 1 / F["h_coldplate"]
    assert np.allclose(T, F["t_coolant"] + F["q_avg"] * R, rtol=1e-9)


@pytest.mark.parametrize("nz", [1, 2])
def test_energy_balance(nz):
    """발열 총량 = 냉각수로 빠져나간 열 총량 (옆면 단열)."""
    x = [55e-6, 20e-6, 10e-6, 5e-6, -5e-6, 2.5, 20e-3, 6e-3, 200e4]
    q = heat_flux_map(x, P)
    _, out = solve_fields(tim_map(x, P), q, 2.5, P, nz=nz, return_flux=True)
    assert np.isclose(out.sum(), q.sum(), rtol=1e-8)


def test_symmetry():
    """대칭 휨 + 중앙 핫스팟 → 온도 지도도 좌우·상하·대각 대칭."""
    x = [50e-6, 20e-6, 20e-6, 0, 0, 3.0, F["die_size"] / 2, F["die_size"] / 2, 150e4]
    _, T = solve(x, P)
    assert np.allclose(T, T[:, ::-1]) and np.allclose(T, T[::-1, :]) and np.allclose(T, T.T)


def test_dome_tim_thin_center_thick_edges():
    """볼록 휨이면 TIM은 가운데가 가장 얇고 모서리가 가장 두껍다."""
    t = tim_map([55e-6, 25e-6, 25e-6, 0, 0, 3.0, 0.01, 0.01, 1e6], P)
    n = t.shape[0]
    assert np.isclose(t.min(), F["tim_min"])
    assert t[n // 2, n // 2] < t[0, 0]


def test_thicker_tim_is_hotter():
    """같은 조건에서 TIM 열전도율이 낮으면 더 뜨겁다."""
    assert solve(flat_x(k_tim=1.5, hs_q=1e6), P)[0] > solve(flat_x(k_tim=6.0, hs_q=1e6), P)[0]


def test_default_layers_close_to_fine_reference():
    """기본 층 수(params) 결과가 실리콘 12층 기준과 0.5 K 이내 (층 수 결정의 회귀 테스트)."""
    a = solve(HOT, P)[0]; ref = solve(HOT, P, nz=12)[0]
    assert abs(a - ref) < 0.5


def test_more_layers_closer_to_reference():
    """층을 늘리면 12층 기준에 가까워진다 (1층 오차 > 2층 오차)."""
    ref = solve(HOT, P, nz=12)[0]
    assert abs(solve(HOT, P, nz=1)[0] - ref) > abs(solve(HOT, P, nz=2)[0] - ref)


def test_warp_constraint_helper():
    assert warp_peak_to_valley(flat_x(), P) == 0.0
