"""격자 수렴 연구: 격자 수에 따라 Tj,max가 얼마나 바뀌나?

기준(참값 대용) = 가장 촘촘한 격자 256. 무작위 입력 30개에서 |Tj,max(n) - Tj,max(256)| 분포를 본다.
실행: python solver/studies/grid_convergence.py
"""
import time, csv
import numpy as np
from _common import load_params, random_inputs, plt, clean, FIG, OUT
from solver import solve

p = load_params()
X = random_inputs(30, p, seed=1)
grids = [16, 24, 32, 48, 64, 96, 128]
ref_n = 256

ref = np.array([solve(x, p, ref_n)[0] for x in X])
err, sec = {}, {}
for n in grids:
    t = time.time(); v = np.array([solve(x, p, n)[0] for x in X]); sec[n] = (time.time() - t) / len(X)
    err[n] = np.abs(v - ref)

with open(OUT / "grid_convergence.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["grid", "mean_abs_err_K", "max_abs_err_K", "sec_per_solve"])
    for n in grids: w.writerow([n, f"{err[n].mean():.4f}", f"{err[n].max():.4f}", f"{sec[n]:.4f}"])

fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
mx = [err[n].max() for n in grids]; mn = [err[n].mean() for n in grids]
axs[0].loglog(grids, mx, "o-", color="#d7301f", label="최대 오차 (30개 중)")
axs[0].loglog(grids, mn, "s--", color="#2f6fd6", label="평균 오차")
axs[0].axhline(0.1, color="#9ca3af", lw=1, ls=":"); axs[0].text(grids[0], 0.11, "0.1 K", fontsize=9, color="#6b7280")
axs[0].axvline(p["solver"]["grid"], color="#1a7f4b", lw=1.2, alpha=.6)
axs[0].text(p["solver"]["grid"] * 1.05, max(mx) * 0.6, f"채택: {p['solver']['grid']}", color="#1a7f4b", fontsize=9.5)
axs[0].set_xlabel("칩 한 변의 격자 수", fontsize=10); axs[0].set_ylabel(f"|Tj,max − Tj,max(격자 {ref_n})| [K]", fontsize=10)
axs[0].set_title("격자를 늘리면 오차가 줄어든다", fontsize=11.5, fontweight="bold", loc="left"); axs[0].legend(fontsize=9, frameon=False)
axs[1].loglog(grids, [sec[n] * 1e3 for n in grids], "o-", color="#6b7280")
axs[1].set_xlabel("칩 한 변의 격자 수", fontsize=10); axs[1].set_ylabel("1회 계산 시간 [ms]", fontsize=10)
axs[1].set_title("대신 계산 시간이 늘어난다", fontsize=11.5, fontweight="bold", loc="left")
for a in axs: clean(a)
fig.tight_layout(w_pad=3); fig.savefig(FIG / "b_grid_convergence_v1.png", dpi=160)

g = p["solver"]["grid"]
print(f"grid {g}: mean {err[g].mean():.3f} K, max {err[g].max():.3f} K, {sec[g]*1e3:.1f} ms/solve")
for n in grids: print(n, round(err[n].mean(), 4), round(err[n].max(), 4), round(sec[n] * 1e3, 1))
