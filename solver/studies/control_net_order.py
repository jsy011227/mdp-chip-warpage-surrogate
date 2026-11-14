"""제어점 3×3 (2차)으로 충분한가?

'실제' 휨 = 4×4 제어점 3차 Bezier 곡면 (돔 + 무작위 고차 성분, 크기 별로 3단계).
이 곡면을 3×3 2차 곡면으로 최소제곱 근사한 뒤, 같은 조건에서 Tj,max 차이를 본다.
실행: python solver/studies/control_net_order.py
"""
import csv
import numpy as np
from scipy.interpolate import BSpline
from _common import load_params, random_inputs, plt, clean, FIG, OUT
from solver import solve_fields, heat_flux_map
from geometry import grid_coords, warp_surface, tim_thickness

p = load_params(); f = p["fixed"]; n = p["solver"]["grid"]; L = f["die_size"]
_, X, Y = grid_coords(L, n); u = (np.arange(n) + 0.5) / n

def basis(deg, m):
    k = np.r_[np.zeros(deg + 1), np.ones(deg + 1)]
    return np.stack([BSpline(k, np.eye(m)[i], deg)(u) for i in range(m)])

B3, B4 = basis(2, 3), basis(3, 4)
A3 = np.einsum("in,jm->nmij", B3, B3).reshape(n * n, 9)      # 3x3 least-squares design matrix

def surf4(P4): return B4.T @ P4 @ B4
def fit3(w):
    c, *_ = np.linalg.lstsq(A3, w.ravel(), rcond=None)
    return warp_surface(c.reshape(3, 3).T, L, n)   # c index = (i_x, j_y) -> P[j][i]

rng = np.random.default_rng(5)
dome4 = np.array([[0, 1, 1, 0], [1, 2.2, 2.2, 1], [1, 2.2, 2.2, 1], [0, 1, 1, 0]], float)
levels = {"작음 (±2 μm)": 2e-6, "보통 (±5 μm)": 5e-6, "큼 (±10 μm)": 10e-6}
cases = random_inputs(60, p, seed=3)
rows, res = [], {}
for lab, amp in levels.items():
    d = []
    for k, x in enumerate(cases):
        base = dome4 / (surf4(dome4).max() - surf4(dome4).min()) * rng.uniform(10e-6, 40e-6)
        w4 = surf4(base + rng.uniform(-amp, amp, (4, 4)))
        w3 = fit3(w4)
        q = heat_flux_map(x, p); kt = x[5]
        T4 = solve_fields(tim_thickness(w4, f["tim_min"]), q, kt, p).max()
        T3 = solve_fields(tim_thickness(w3, f["tim_min"]), q, kt, p).max()
        resid = np.abs(w4 - w3).max()
        d.append(T3 - T4); rows.append([lab, k, f"{resid*1e6:.3f}", f"{T4:.3f}", f"{T3:.3f}", f"{T3-T4:.4f}"])
    res[lab] = np.array(d)

with open(OUT / "control_net_order.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["level", "case", "max_shape_residual_um", "Tmax_4x4", "Tmax_3x3", "diff_K"]); w.writerows(rows)

fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.9))
# example surfaces
base = dome4 / (surf4(dome4).max() - surf4(dome4).min()) * 30e-6
w4 = surf4(base + np.random.default_rng(9).uniform(-10e-6, 10e-6, (4, 4))); w3 = fit3(w4)
j = n // 3
axs[0].plot(X[j] * 1e3, w4[j] * 1e6, color="#d7301f", lw=2, label="'실제' 휨 (4×4, 고차 성분 ±10 μm)")
axs[0].plot(X[j] * 1e3, w3[j] * 1e6, color="#2f6fd6", lw=2, ls="--", label="3×3 근사")
axs[0].set_xlabel("x [mm]  (y = 칩의 1/3 지점 단면)", fontsize=10); axs[0].set_ylabel("휨 높이 [μm]", fontsize=10)
axs[0].set_title("3×3은 큰 모양은 따라가지만 비대칭 성분 일부를 놓친다", fontsize=11.5, fontweight="bold", loc="left"); axs[0].legend(fontsize=8.5, frameon=False)
axs[1].boxplot([np.abs(v) for v in res.values()], tick_labels=list(res.keys()), widths=.5, patch_artist=True,
               boxprops=dict(facecolor="#dbe7fa", color="#2f6fd6"), medianprops=dict(color="#1f2328"))
axs[1].axhline(1.0, color="#d7301f", ls=":", lw=1.2); axs[1].text(3.3, 1.05, "목표 1 K", color="#d7301f", fontsize=9)
axs[1].set_ylabel("|Tj,max(3×3) − Tj,max(4×4)| [K]", fontsize=10); axs[1].set_xlabel("4×4 곡면의 고차 성분 크기", fontsize=10)
axs[1].set_title("고차 성분이 놓친 온도 오차 (조건 60개씩)", fontsize=11.5, fontweight="bold", loc="left")
for a in axs: clean(a)
fig.tight_layout(w_pad=3); fig.savefig(FIG / "b_control_net_3x3_vs_4x4_v1.png", dpi=160)
for lab, v in res.items(): print(lab, "mean", round(np.abs(v).mean(), 3), "max", round(np.abs(v).max(), 3))
