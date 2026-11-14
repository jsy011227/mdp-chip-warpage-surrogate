"""layer_count.csv → 실리콘 층 수 그림 (1층 이상만) → figures/b_silicon_layers_v2.png"""
import csv
import numpy as np
from _common import plt, FIG, OUT, clean

rows = [r for r in csv.DictReader(open(OUT / "layer_count.csv")) if int(r["layers"]) >= 1]
L = [int(r["layers"]) for r in rows]
rnd = [float(r["random_maxabs_dTmax_K"]) for r in rows]; wst = [float(r["worst_maxabs_dTmax_K"]) for r in rows]
ms = [float(r["ms_per_solve"]) for r in rows]
xs = np.arange(len(L)); lab = [f"{n}층\n(채택)" if n == 2 else f"{n}층" for n in L]
fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.9))
axs[0].bar(xs - .2, rnd, .38, color="#2f6fd6", label="무작위 24개 중 최대")
axs[0].bar(xs + .2, wst, .38, color="#d7301f", label="최악 조건 3개 중 최대")
axs[0].set_yscale("log"); axs[0].axhline(.5, color="#6b7280", ls=":", lw=1); axs[0].text(len(L) - .8, .55, "0.5 K", fontsize=9, color="#6b7280")
axs[0].set_xticks(xs, lab); axs[0].set_ylabel("|Tj,max − Tj,max(12층)| [K]", fontsize=10)
axs[0].set_title("층을 나누면 오차가 빠르게 줄어든다", fontsize=11.5, fontweight="bold", loc="left"); axs[0].legend(fontsize=9, frameon=False)
axs[1].bar(xs, ms, .5, color="#9ca3af")
for i, m in enumerate(ms): axs[1].text(i, m * 1.04 + 5, f"{m:.0f} ms", ha="center", fontsize=9)
axs[1].set_xticks(xs, lab); axs[1].set_ylabel("1회 계산 시간 [ms]", fontsize=10)
axs[1].set_title("대신 계산 시간이 늘어난다 (격자 64×64)", fontsize=11.5, fontweight="bold", loc="left")
for a in axs: clean(a)
fig.tight_layout(w_pad=3); fig.savefig(FIG / "b_silicon_layers_v2.png", dpi=160)
