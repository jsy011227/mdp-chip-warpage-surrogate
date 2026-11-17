"""칩 평면 열 해석기 (정상상태) — 실리콘을 두께 방향으로 si_layers층 나눈 '다층 2.5D' 모델.

구조 (위 → 아래):
    냉각수 T_c ─ 1/h ─ TIM (위치별 1D 저항 t_tim(x,y)/k_tim) ─ 실리콘 si_layers층 ─ 능동면 (발열 q(x,y))
    - 실리콘: x, y, z 모두 열전도를 푼다 (층마다 x-y 격자, 층 사이 z 연결). 옆면 단열.
    - TIM: 가로 전도는 실리콘의 약 0.2%라 무시 → 위치마다 독립된 1D 저항 (정상상태에서 정확).
    - 냉각판: 열전달계수 h 하나 (바닥판의 가로 퍼짐 무시 → 보수적).
    - 출력 온도 = 능동면(아랫면) 온도 Tj(x, y).

si_layers 선택 근거 (solver/studies/layer_count.py): 실리콘 비오 수 Bi_z = 0.09~0.29 > 0.1 → 1층은
핫스팟 온도를 2.8 K(최대 3.9 K) 낮게 계산. 2층이면 12층 대비 최대 0.21 K(최악 조건 0.44 K), 1회 약 40 ms.

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


def solve_fields(t_tim, q, k_tim, p, nz=None, return_flux=False):
    """TIM 두께 지도와 발열 지도 → 능동면 온도 지도 Tj [°C] (k_tim은 숫자 또는 지도).

    nz: 실리콘 층 수 (기본 params solver.si_layers). return_flux=True면 냉각수로 빠져나간
    위치별 열유속 [W/m2]도 함께 돌려준다 (에너지 보존 확인용).
    """
    f = p["fixed"]; n = t_tim.shape[0]; nz = int(nz or p["solver"].get("si_layers", 1))
    L, H, k = f["die_size"], f["die_thickness"], f["k_si"]
    dx, dz = L / n, H / nz
    idx = np.arange(nz * n * n).reshape(nz, n, n)          # (z, y, x), z = 0 이 능동면(아래)
    rows, cols, vals = [], [], []
    diag = np.zeros(idx.size)

    def link(a, b, g):                                     # 두 셀 사이 열전도 [W/K]
        rows.extend([a, b]); cols.extend([b, a]); vals.extend([-g, -g])
        np.add.at(diag, a, g); np.add.at(diag, b, g)

    gx, gz = k * dz, k * dx * dx / dz
    link(idx[:, :, :-1].ravel(), idx[:, :, 1:].ravel(), np.full(nz * n * (n - 1), gx))   # x 방향
    link(idx[:, :-1, :].ravel(), idx[:, 1:, :].ravel(), np.full(nz * n * (n - 1), gx))   # y 방향
    if nz > 1:
        link(idx[:-1].ravel(), idx[1:].ravel(), np.full((nz - 1) * n * n, gz))           # z 방향
    top = idx[-1].ravel()
    R_top = dz / 2 / k + (t_tim / np.broadcast_to(k_tim, t_tim.shape)).ravel() + 1 / f["h_coldplate"]
    g_top = dx * dx / R_top                                 # 윗셀 중심 → 냉각수 [W/K]
    diag[top] += g_top
    rhs = np.zeros(idx.size)
    rhs[top] += g_top * f["t_coolant"]
    rhs[idx[0].ravel()] += q.ravel() * dx * dx              # 능동면 발열
    allidx = np.arange(idx.size)
    A = sp.csr_matrix((np.concatenate(vals + [diag]), (np.concatenate(rows + [allidx]), np.concatenate(cols + [allidx]))),
                      shape=(idx.size, idx.size))
    T = spla.spsolve(A.tocsc(), rhs).reshape(nz, n, n)
    Tj = T[0] + q * (dz / 2) / k                            # 아랫셀 중심 → 능동면으로 보정
    if return_flux:
        return Tj, ((T[-1].ravel() - f["t_coolant"]) / R_top).reshape(n, n)
    return Tj


def solve(x, p=None, n=None, nz=None):
    """입력 x (9개) → (Tj,max [°C], 능동면 온도 지도 Tj [°C], 모양 (n, n), 행 = y)."""
    p = p or load_params()
    v = unpack(x, p)
    T = solve_fields(tim_map(x, p, n), heat_flux_map(x, p, n), v["k_tim"], p, nz)
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
    print(f"Tj,max = {tmax:.2f} °C  ({time.time() - t:.3f} s, grid {T.shape}, Si {p['solver']['si_layers']} layers)")
