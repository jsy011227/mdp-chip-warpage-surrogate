"""해석기 검증 테스트.  실행: pytest -q"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "solver"))
from solver import load_params, solve, solve_fields, tim_map, heat_flux_map, warp_peak_to_valley  # noqa: E402

P = load_params()
F = P["fixed"]


def flat_x(k_tim=3.0, hs_q=0.0):
    return [0, 0, 0, 0, 0, k_tim, F["die_size"] / 2, F["die_size"] / 2, hs_q]


def test_uniform_matches_analytic():
    """휨 없음 + 균일 발열 → T = T_c + q * R (손계산과 정확히 일치해야 함)."""
    tmax, T = solve(flat_x(), P)
    R = F["die_thickness"] / F["k_si"] + F["tim_min"] / 3.0 + 1 / F["h_coldplate"]
    exact = F["t_coolant"] + F["q_avg"] * R
    assert np.allclose(T, exact, rtol=1e-9)


def test_energy_balance():
    """발열 총량 = 냉각판으로 빠져나간 열 총량 (가장자리 단열)."""
    x = [55e-6, 20e-6, 10e-6, 5e-6, -5e-6, 2.5, 20e-3, 6e-3, 200e4]
    _, T = solve(x, P)
    t_tim = tim_map(x, P); q = heat_flux_map(x, P)
    R = F["die_thickness"] / F["k_si"] + t_tim / 2.5 + 1 / F["h_coldplate"]
    assert np.isclose(((T - F["t_coolant"]) / R).sum(), q.sum(), rtol=1e-8)


def test_symmetry():
    """대칭 휨 + 중앙 핫스팟 → 온도 지도도 좌우·상하 대칭."""
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


def test_warp_constraint_helper():
    assert warp_peak_to_valley(flat_x(), P) == 0.0
