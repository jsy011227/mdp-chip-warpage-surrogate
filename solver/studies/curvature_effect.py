"""실리콘을 평면으로 근사한 오차: 실제로 휜(곡률이 큰) 실리콘으로 풀어 비교한다.

- 두 경우 모두 TIM 두께 지도 t(x,y)와 발열 q(x,y)는 같다 (같은 문제).
- '평면': 실리콘 z ∈ [-H, 0]  (현재 2.5D 해석기의 가정)
- '휨':   실리콘 z ∈ [-H, 0] + s·w(x,y)  (다이 자체가 휜 껍질, s = 휨 배율)
- 3D 유한요소 (육면체 1차, scikit-fem). 같은 격자 연결을 쓰고 절점 높이만 바꾸므로
  이산화 오차는 상쇄되고 '곡률 효과'만 남는다.
- 경계: 아랫면 열유속 q, 윗면 TIM+냉각판 직렬 저항 (투영면적 기준 h_eff = 1/R), 옆면 단열.
  기울어진 면에서는 실제 면적이 √(1+|∇w|²)배라 면적당 값을 그만큼 나눈다 (총 발열·총 컨덕턴스 동일).
결과 → results/curvature_effect.csv, figures/b_curvature_effect_v1.png
"""
import csv, sys, time
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
from scipy.sparse.linalg import spsolve
import skfem as fem
from skfem.helpers import dot, grad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from solver import load_params, input_names, unpack, solve, warp_peak_to_valley                      # noqa: E402
from _common import plt, FIG, OUT, clean, random_inputs                         # noqa: E402

p = load_params(); f = p["fixed"]
L, H, KSI = f["die_size"], f["die_thickness"], f["k_si"]
KN = np.array([0, 0, 0, 1, 1, 1.0])
B = [BSpline(KN, np.eye(3)[k], 2) for k in range(3)]


def warp_fn(v):
    P = np.array([[0, v["w_mid_x"], 0], [v["w_mid_y"], v["w_center"], v["w_mid_y"]], [0, v["w_mid_x"], 0]])
    def w(X, Y):
        bx = np.stack([b(np.clip(X / L, 0, 1)) for b in B]); by = np.stack([b(np.clip(Y / L, 0, 1)) for b in B])
        return np.einsum("j...,ji,i...->...", by, P, bx) + v["tilt_x"] * (X / L - .5) + v["tilt_y"] * (Y / L - .5)
    g = np.linspace(0, L, 401); GX, GY = np.meshgrid(g, g)
    return w, w(GX, GY).max()


def slope2(w, X, Y, e=1e-6):
    wx = (w(X + e, Y) - w(X - e, Y)) / (2 * e); wy = (w(X, Y + e) - w(X, Y - e)) / (2 * e)
    return wx ** 2 + wy ** 2


NXY, NZ = 64, 6
m0 = fem.MeshHex.init_tensor(np.linspace(0, L, NXY + 1), np.linspace(0, L, NXY + 1), np.linspace(-H, 0, NZ + 1))
TOP = m0.facets_satisfying(lambda x: x[2] > -1e-12)
BOT = m0.facets_satisfying(lambda x: x[2] < -H + 1e-12)
bot_nodes = np.unique(m0.facets[:, BOT])
E = fem.ElementHex1()


@fem.BilinearForm
def cond(u, v, w): return KSI * dot(grad(u), grad(v))


@fem.BilinearForm
def robin(u, v, w): return w["hA"] * u * v


@fem.LinearForm
def robin_rhs(v, w): return w["hA"] * f["t_coolant"] * v


@fem.LinearForm
def flux(v, w): return w["qA"] * v


def solve_fem(v, s):
    """s = 0 → 평면 실리콘, s = 1 → 실제 휨, s > 1 → 휨을 과장."""
    w, wmax = warp_fn(v)
    X, Y, Z = m0.p
    m = fem.MeshHex(np.vstack([X, Y, Z + s * w(X, Y)]), m0.t)
    ib = fem.Basis(m, E)
    ft, fb = fem.FacetBasis(m, E, facets=TOP), fem.FacetBasis(m, E, facets=BOT)

    def coeffs(fbasis):
        x, y = fbasis.global_coordinates().value[:2]
        corr = np.sqrt(1 + s ** 2 * slope2(w, x, y))                 # 실제 면적 / 투영 면적
        t = f["tim_min"] + (wmax - w(x, y))                          # TIM 두께 (두 경우 동일)
        heff = 1 / (t / v["k_tim"] + 1 / f["h_coldplate"])
        q = f["q_avg"] + v["hs_q"] * np.exp(-((x - v["hs_x"]) ** 2 + (y - v["hs_y"]) ** 2) / (2 * f["hs_sigma"] ** 2))
        return heff / corr, q / corr

    hA, _ = coeffs(ft); _, qA = coeffs(fb)
    A = cond.assemble(ib) + robin.assemble(ft, hA=hA)
    b = robin_rhs.assemble(ft, hA=hA) + flux.assemble(fb, qA=qA)
    T = spsolve(A.tocsr(), b)
    return T[bot_nodes].max()


if __name__ == "__main__":
    names = input_names(p); rng = np.random.default_rng(5)
    S = [0.0, 1.0, 10.0, 30.0]
    X = random_inputs(400, p, seed=5)
    X = np.array([x for x in X if warp_peak_to_valley(x, p) <= p['constraints']['max_warp_peak_to_valley']])[:16]
    rows = []; t0 = time.time()
    for i, x in enumerate(X):
        v = unpack(x, p)
        T = {s: solve_fem(v, s) for s in S}
        fd = solve(x, p)[0]
        rows.append([i] + [T[s] for s in S] + [fd])
        print(i, " ".join(f"{T[s] - T[0]:+.4f}" for s in S[1:]), f"FD-FEM {fd - T[0]:+.3f}", f"{time.time() - t0:.0f}s", flush=True)
    R = np.array(rows)
    with open(OUT / "curvature_effect.csv", "w", newline="") as fh:
        wtr = csv.writer(fh); wtr.writerow(["case"] + [f"Tjmax_s{s:g}" for s in S] + ["Tjmax_FD"]); wtr.writerows(np.round(R, 5).tolist())

    d = R[:, 2:2 + len(S) - 1] - R[:, [1]]
    print("mean |dT|", np.abs(d).mean(0), "max", np.abs(d).max(0))
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.loglog(S[1:], np.abs(d).mean(0), "o-", color="#2f6fd6", label="평균 |ΔTj,max|")
    ax.loglog(S[1:], np.abs(d).max(0), "s-", color="#d7301f", label="최대 |ΔTj,max|")
    for s, y in zip(S[1:], np.abs(d).max(0)): ax.text(s, y * 1.35, f"{y:.3g} K", ha="center", fontsize=9)
    ax.axvline(1, color="#6b7280", ls=":", lw=1); ax.text(1.06, ax.get_ylim()[0] * 1.6, "실제 휨\n(곡률 반경 ≈ 1.8 m)", fontsize=9, color="#6b7280")
    ax.set_xlabel("휨 배율 (1 = 실제, 30 = 곡률 반경 약 6 cm)"); ax.set_ylabel("|Tj,max(휜 실리콘) − Tj,max(평면)| [K]")
    ax.set_title("실리콘을 평면으로 근사한 오차", fontsize=11.5, fontweight="bold", loc="left"); ax.legend(frameon=False, fontsize=9)
    clean(ax); fig.tight_layout(); fig.savefig(FIG / "b_curvature_effect_v1.png", dpi=160)
