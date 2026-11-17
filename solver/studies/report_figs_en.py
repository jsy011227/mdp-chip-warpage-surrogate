"""English report figures, re-plotted from results/ (no re-computation). → figures/b_*_en.png"""
import csv, json
import numpy as np
from _common import plt, FIG, OUT, clean

C1, C2, G = "#2f6fd6", "#d7301f", "#6b7280"
rd = lambda f: list(csv.DictReader(open(OUT / f)))

# grid convergence
r = rd("grid_convergence.csv"); n = [int(x["grid"]) for x in r]
fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.loglog(n, [float(x["mean_abs_err_K"]) for x in r], "o-", c=C1, label="mean")
ax.loglog(n, [float(x["max_abs_err_K"]) for x in r], "s-", c=C2, label="max")
ax.axvline(64, c=G, ls=":"); ax.text(66, .02, "64 (used)", color=G, fontsize=9)
ax.set_xlabel("grid n (n×n)"); ax.set_ylabel("|Tj,max − Tj,max(256)| [K]"); ax.legend(frameon=False); clean(ax)
fig.tight_layout(); fig.savefig(FIG / "b_grid_en.png", dpi=160)

# silicon layers
r = [x for x in rd("layer_count.csv") if int(x["layers"]) >= 1]; L = [int(x["layers"]) for x in r]; xs = np.arange(len(L))
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.4))
axs[0].bar(xs - .2, [float(x["random_maxabs_dTmax_K"]) for x in r], .38, color=C1, label="random 24 (max)")
axs[0].bar(xs + .2, [float(x["worst_maxabs_dTmax_K"]) for x in r], .38, color=C2, label="worst 3 (max)")
axs[0].set_yscale("log"); axs[0].set_ylabel("|ΔTj,max| vs 12 layers [K]"); axs[0].legend(frameon=False, fontsize=9)
ms = [float(x["ms_per_solve"]) for x in r]; axs[1].bar(xs, ms, .5, color="#9ca3af"); axs[1].set_ylabel("ms / solve")
for a in axs: a.set_xticks(xs, [f"{k}{' (used)' if k == 2 else ''}" for k in L]); a.set_xlabel("silicon layers"); clean(a)
fig.tight_layout(); fig.savefig(FIG / "b_layers_en.png", dpi=160)

# 3x3 fit error vs higher-order amplitude
r = rd("control_net_order.csv"); lv = [("작음 (±2 μm)", "±2 μm"), ("보통 (±5 μm)", "±5 μm"), ("큼 (±10 μm)", "±10 μm")]
fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.boxplot([[abs(float(x["diff_K"])) for x in r if x["level"] == k] for k, _ in lv], tick_labels=[e for _, e in lv], showfliers=True)
ax.set_xlabel("higher-order amplitude (truth = 4×4 cubic)"); ax.set_ylabel("|Tj,max(3×3) − Tj,max(truth)| [K]"); clean(ax)
fig.tight_layout(); fig.savefig(FIG / "b_cn3_fit_en.png", dpi=160)

# 3x3 vs 4x4 vs 5x5
S = json.load(open(OUT / "control_net_4x4_summary.json"))["summary"]; rat = ["0.1", "0.2", "0.3"]
show = [("Q3-sym", "3×3 sym (used)", C1, "o-"), ("Q3", "3×3 full", "#60a5fa", "o--"), ("C4", "4×4 cubic", C2, "s-"),
        ("Q4", "4×4 quad", "#f59e0b", "^-"), ("Q5-sym", "5×5 sym", "#1a7f4b", "D-")]
fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.4), sharey=True)
for ax, fa, t in zip(axs, ["skew", "flat", "random"], ["asymmetric", "symmetric (flat top)", "random"]):
    for k, nm, c, st in show:
        ax.plot([10, 20, 30], [S[fa][q][k]["mean"] for q in rat], st, c=c, label=nm, ms=5)
    ax.set_title(t, loc="left", fontweight="bold"); ax.set_xlabel("higher-order [% of dome]"); ax.set_xticks([10, 20, 30]); clean(ax)
axs[0].set_ylabel("mean |ΔTj,max| [K]"); axs[0].legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig(FIG / "b_cn4_en.png", dpi=160)

# curved vs flat silicon
R = np.array([[float(v) for v in x.values()] for x in rd("curvature_effect.csv")]); s = [1, 10, 30]
d = np.abs(R[:, 2:5] - R[:, [1]])
fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.loglog(s, d.mean(0), "o-", c=C1, label="mean"); ax.loglog(s, d.max(0), "s-", c=C2, label="max")
for a, b in zip(s, d.max(0)): ax.text(a, b * 1.3, f"{b:.3f}", ha="center", fontsize=9)
ax.set_xticks(s, ["1\n(real, R≈1.8 m)", "10\n(R≈18 cm)", "30\n(R≈6 cm)"])
ax.set_xlabel("warp scale"); ax.set_ylabel("|Tj,max(curved) − Tj,max(flat)| [K]"); ax.legend(frameon=False); clean(ax)
fig.tight_layout(); fig.savefig(FIG / "b_curv_en.png", dpi=160)
print("ok")
