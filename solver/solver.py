"""2.5D 칩 평면 열 해석기 (정상상태).

지배 방정식 (칩 평면 x-y, 가장자리 단열):
    k_si * t_die * ∇²T - (T - T_c) / R(x, y) + q(x, y) = 0
    R(x, y) = t_die / k_si + t_tim(x, y) / k_tim + 1 / h

사용법:
    from solver import load_params, solve
    p = load_params()
    tmax, T = solve(x, p)          # x: 입력 9개 (params.yaml의 inputs 순서, SI 단위)
"""
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import yaml

from geometry import grid_coords, control_net, warp_surface, tilt_plane, tim_thickness

ROOT = Path(__file__).resolve().parents[1]


def _numify(o):
    """YAML이 '80.0e4' 같은 값을 문자열로 읽는 문제 방지: 숫자처럼 생긴 문자열은 float로 바꾼다."""
    if isinstance(o, dict): return {k: _numify(v) for k, v in o.items()}
    if isinstance(o, list): return [_numify(v) for v in o]
    if isinstance(o, str):
        try: return float(o)
        except ValueError: return o
    return o


def load_params(path=ROOT / "params.yaml"):
    with open(path, encoding="utf-8") as f:
        return _numify(yaml.safe_load(f))


def input_names(p):
    return [d["name"] for d in p["inputs"]]


def unpack(x, p):
    names = input_names(p)
    x = np.asarray(x, float)
    if x.shape != (len(names),):
        raise ValueError(f"x must have {len(names)} values in order {names}")
    return dict(zip(names, x))


def tim_map(x, p, n=None):
    """입력 x → TIM 두께 지도 [m]."""
    f = p["fixed"]; n = n or p["solver"]["grid"]; v = unpack(x, p)
    L = f["die_size"]
    w = warp_surface(control_net(v["w_center"], v["w_mid_x"], v["w_mid_y"]), L, n)
    w = w + tilt_plane(v["tilt_x"], v["tilt_y"], L, n)
    return tim_thickness(w, f["tim_min"])


def heat_flux_map(x, p, n=None):
    """입력 x → 발열 지도 q(x, y) [W/m2]: 평균 + 가우시안 핫스팟."""
    f = p["fixed"]; n = n or p["solver"]["grid"]; v = unpack(x, p)
    _, X, Y = grid_coords(f["die_size"], n)
    r2 = (X - v["hs_x"]) ** 2 + (Y - v["hs_y"]) ** 2
    return f["q_avg"] + v["hs_q"] * np.exp(-r2 / (2 * f["hs_sigma"] ** 2))


def solve_fields(t_tim, q, k_tim, p):
    """TIM 두께 지도와 발열 지도로 온도 지도를 푼다 (k_tim은 숫자 또는 지도)."""
    f = p["fixed"]; n = t_tim.shape[0]
    L, t_die, k_si = f["die_size"], f["die_thickness"], f["k_si"]
    R = t_die / k_si + t_tim / np.broadcast_to(k_tim, t_tim.shape) + 1 / f["h_coldplate"]
    dx = L / n
    c = k_si * t_die / dx ** 2
    idx = np.arange(n * n).reshape(n, n)
    rows, cols, vals = [], [], []
    diag = 1.0 / R.ravel()
    for di, dj in ((0, 1), (0, -1), (1, 0), (-1, 0)):     # 이웃 연결 (가장자리 = 단열)
        a = idx[max(0, -di):n - max(0, di), max(0, -dj):n - max(0, dj)].ravel()
        b = idx[max(0, di):n - max(0, -di), max(0, dj):n - max(0, -dj)].ravel()
        rows.append(a); cols.append(b); vals.append(np.full(a.size, -c))
        diag[a] += c
    A = sp.csr_matrix((np.concatenate(vals + [diag]),
                       (np.concatenate(rows + [idx.ravel()]), np.concatenate(cols + [idx.ravel()]))),
                      shape=(n * n, n * n))
    rhs = q.ravel() + f["t_coolant"] / R.ravel()
    return spla.spsolve(A.tocsc(), rhs).reshape(n, n)


def solve(x, p=None, n=None):
    """입력 x (9개) → (Tj,max [°C], 온도 지도 T [°C], 모양 (n, n), 행 = y)."""
    p = p or load_params()
    v = unpack(x, p)
    T = solve_fields(tim_map(x, p, n), heat_flux_map(x, p, n), v["k_tim"], p)
    return float(T.max()), T


def warp_peak_to_valley(x, p, n=None):
    """기울기를 뺀 휨 곡면의 최고-최저 차 [m] (DOE 제약 확인용)."""
    f = p["fixed"]; n = n or p["solver"]["grid"]; v = unpack(x, p)
    w = warp_surface(control_net(v["w_center"], v["w_mid_x"], v["w_mid_y"]), f["die_size"], n)
    return float(w.max() - w.min())


if __name__ == "__main__":
    import time
    p = load_params()
    x = [55e-6, 25e-6, 25e-6, 0, 0, 3.0, 22e-3, 4e-3, 170e4]
    t = time.time(); tmax, T = solve(x, p)
    print(f"Tj,max = {tmax:.2f} °C  ({time.time() - t:.3f} s, grid {T.shape})")
