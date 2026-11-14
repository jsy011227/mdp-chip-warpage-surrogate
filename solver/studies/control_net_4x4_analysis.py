"""휨 곡면 제어점 3×3 → 4×4 변경이 얼마나 유의미한가.

비교하는 표현 (방향마다 1D 기저의 텐서곱, 상수항은 TIM 두께 식에서 상쇄되므로 항상 포함)
  Q3-sym : 3×3 2차, x·y 좌우대칭          (현재 프로젝트: 휨 변수 3개 + 기울기 2개)
  Q3     : 3×3 2차, 대칭 없음              (변수 8개)
  C4     : 4×4 3차 Bezier, 대칭 없음       (변수 15개)
  Q4     : 4×4 2차 B-spline (knot 0.5)     (변수 15개)
  Q4-sym : 4×4 2차 B-spline, 좌우대칭      (변수 3개 + 기울기 2개)
  Q5-sym : 5×5 2차 B-spline, 좌우대칭      (변수 8개 + 기울기 2개)

'실제' 휨 = 돔(진폭 A) + 고차 성분(진폭 r·A) 네 가족
  dome   : 고차 성분 없음
  skew   : 비대칭 3차 성분 (꼭대기가 한쪽으로 쏠림)
  flat   : 대칭 4차 성분 (꼭대기가 평평하거나 가장자리만 처짐)
  random : 8×8 2차 B-spline 무작위 성분 (일반적인 울퉁불퉁)
각 표현으로 최소제곱 근사 → 같은 조건(무작위 TIM·핫스팟)에서 Tj,max 차이를 본다.
실행: python solver/studies/control_net_4x4_analysis.py
"""
import json, time
import numpy as np
from scipy.interpolate import BSpline
from _common import load_params, random_inputs, plt, clean, FIG, OUT
from solver import solve_fields, heat_flux_map
from geometry import grid_coords

p = load_params(); f = p["fixed"]; n = p["solver"]["grid"]; L = f["die_size"]
u = (np.arange(n) + 0.5) / n; s = 2 * u - 1
_, X, Y = grid_coords(L, n)


def basis1d(deg, m):
    inner = np.linspace(0, 1, m - deg + 1)[1:-1]
    k = np.r_[np.zeros(deg + 1), inner, np.ones(deg + 1)]
    return np.stack([BSpline(k, np.eye(m)[i], deg)(u) for i in range(m)])        # (m, n)


def sym1d(B):
    m = B.shape[0]; cols = []
    for i in range((m + 1) // 2):
        j = m - 1 - i
        cols.append(B[i] + (B[j] if j != i else 0))
    return np.stack(cols)


REPS = {
    "Q3-sym": sym1d(basis1d(2, 3)), "Q3": basis1d(2, 3), "C4": basis1d(3, 4),
    "Q4": basis1d(2, 4), "Q4-sym": sym1d(basis1d(2, 4)), "Q5-sym": sym1d(basis1d(2, 5)),
}
NVAR = {"Q3-sym": "3 + 기울기 2", "Q3": "8", "C4": "15", "Q4": "15", "Q4-sym": "3 + 기울기 2", "Q5-sym": "8 + 기울기 2"}
NIN = {"Q3-sym": 9, "Q3": 12, "C4": 19, "Q4": 19, "Q4-sym": 9, "Q5-sym": 14}
DESIGN = {k: np.einsum("an,bm->mnab", B, B).reshape(n * n, -1) for k, B in REPS.items()}   # rows: (y, x)


def fit(name, w):
    D = DESIGN[name]
    c, *_ = np.linalg.lstsq(D, w.ravel(), rcond=None)
    return (D @ c).reshape(n, n)


B8 = basis1d(2, 8)
def truth(family, A, r, rng):
    dome = A * (2 - s[None, :] ** 2 - s[:, None] ** 2) / 2
    if family == "dome": hi = 0 * dome
    elif family == "skew":
        a, b = rng.uniform(-1, 1, 2); hi = (a * s[None, :] ** 3 + b * s[:, None] ** 3) / 2
    elif family == "flat":
        a, b = rng.uniform(-1, 1, 2); hi = 4 * (a * (s[None, :] ** 4 - s[None, :] ** 2) + b * (s[:, None] ** 4 - s[:, None] ** 2)) / 2
    else:
        hi = B8.T @ rng.uniform(-1, 1, (8, 8)) @ B8
    if np.ptp(hi) > 0: hi = hi / np.ptp(hi)
    return dome + r * A * hi


def tmax(w, x):
    gap = f["tim_min"] + (w.max() - w)
    return solve_fields(gap, heat_flux_map(x, p), x[5], p).max()


if __name__ == "__main__":
    t0 = time.time()
    fams = ["dome", "skew", "flat", "random"]; ratios = [0.1, 0.2, 0.3]; NC = 30
    reps = list(REPS)
    res = {fa: {str(r): {k: {"dT": [], "shape": []} for k in reps} for r in ratios} for fa in fams}
    rng = np.random.default_rng(21)
    for fa in fams:
        for r in (ratios if fa != "dome" else ratios[:1]):
            X_ = random_inputs(NC, p, seed=int(100 * r) + len(fa))
            for x in X_:
                A = rng.uniform(20e-6, 40e-6) * rng.choice([1, 1, 1, -1])
                wt = truth(fa, A, r, rng); Tt = tmax(wt, x)
                for k in reps:
                    wf = fit(k, wt)
                    res[fa][str(r)][k]["dT"].append(float(tmax(wf, x) - Tt))
                    res[fa][str(r)][k]["shape"].append(float(np.abs(wf - wt).max() * 1e6))
    json.dump(res, open(OUT / "control_net_4x4_analysis.json", "w"))
    print(f"done {time.time() - t0:.0f} s")

    def stat(fa, r, k, key="dT", fn=lambda v: np.abs(v).mean()):
        return fn(np.array(res[fa][str(r)][k][key]))

    # ---- Fig 1: basis functions ----
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.3), sharey=True)
    for ax, (nm, B, t) in zip(axs, [("3×3 (2차)", basis1d(2, 3), "현재: 기저 3개, 하나가 전체에 영향"),
                                    ("4×4 (3차 Bezier)", basis1d(3, 4), "기저 4개, 여전히 전체에 영향"),
                                    ("4×4 (2차 B-spline, knot 0.5)", basis1d(2, 4), "기저 4개, 반쪽씩 담당 (국소 제어)")]):
        for i, b in enumerate(B): ax.plot(u * L * 1e3, b, lw=1.8, color=["#2f6fd6", "#d7301f", "#1a7f4b", "#9333ea", "#c2410c"][i])
        if "knot" in nm: ax.axvline(L * 1e3 / 2, color="#9ca3af", ls=":", lw=1)
        ax.set_title(nm + "\n", fontsize=11.5, fontweight="bold", loc="left"); ax.text(0.0, 1.02, t, transform=ax.transAxes, fontsize=9.2, color="#6b7280", ha="left", va="bottom")
        ax.set_xlabel("x [mm]", fontsize=9.5); clean(ax)
    axs[0].set_ylabel("기저함수 값 (제어점의 영향)", fontsize=9.5)
    fig.tight_layout(); fig.savefig(FIG / "b_cn4_basis_v1.png", dpi=160)

    # ---- Fig 2: example profiles for each family ----
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    rr = np.random.default_rng(3); j = n // 2
    for ax, fa, ttl in zip(axs, ["skew", "flat", "random"], ["비대칭 (꼭대기 쏠림)", "대칭 고차 (평평한 꼭대기·가장자리 처짐)", "무작위 고차"]):
        wt = truth(fa, 35e-6, 0.3, rr)
        ax.plot(X[j] * 1e3, wt[j] * 1e6, color="#1f2328", lw=2.4, label="'실제' 휨 (고차 30%)")
        for k, c, ls in (("Q3-sym", "#2f6fd6", "--"), ("C4", "#d7301f", "-."), ("Q5-sym", "#1a7f4b", ":")):
            ax.plot(X[j] * 1e3, fit(k, wt)[j] * 1e6, color=c, lw=1.8, ls=ls, label={"Q3-sym": "3×3 대칭 (현재)", "C4": "4×4 3차", "Q5-sym": "5×5 대칭"}[k])
        ax.set_title(ttl, fontsize=11, fontweight="bold", loc="left"); ax.set_xlabel("x [mm] (칩 가운데 단면)", fontsize=9.5); clean(ax)
    axs[0].set_ylabel("휨 높이 [μm]", fontsize=9.5); axs[0].legend(fontsize=8.5, frameon=False, loc="lower center")
    fig.tight_layout(); fig.savefig(FIG / "b_cn4_profiles_v1.png", dpi=160)

    # ---- Fig 3: mean |dTmax| vs ratio per family and representation ----
    show = [("Q3-sym", "3×3 대칭 (현재)", "#2f6fd6", "o-"), ("Q3", "3×3 전체", "#60a5fa", "o--"),
            ("C4", "4×4 3차", "#d7301f", "s-"), ("Q4", "4×4 2차+knot", "#f59e0b", "^-"), ("Q5-sym", "5×5 대칭", "#1a7f4b", "D-")]
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.7), sharey=True)
    for ax, fa, ttl in zip(axs, ["skew", "flat", "random"], ["비대칭 고차 성분", "대칭 고차 성분", "무작위 고차 성분"]):
        for k, nm, c, st in show:
            ax.plot([r * 100 for r in ratios], [stat(fa, r, k) for r in ratios], st, color=c, label=nm, lw=1.8, ms=5)
        ax.axhline(0.15, color="#9ca3af", ls=":", lw=1); ax.text(30, 0.12, "해석기 자체 오차 0.15 K", fontsize=8.5, color="#6b7280", ha="right", va="top")
        ax.set_title(ttl, fontsize=11.5, fontweight="bold", loc="left"); ax.set_xlabel("고차 성분 크기 [돔 진폭의 %]", fontsize=9.5)
        ax.set_xticks([10, 20, 30]); clean(ax)
    axs[0].set_ylabel("평균 |Tj,max(근사) − Tj,max(실제)| [K]", fontsize=9.5); axs[0].legend(fontsize=8.5, frameon=False, loc="upper left")
    axs[0].annotate("3×3 대칭 = 5×5 대칭\n(대칭 표현은 비대칭 성분을\n전혀 못 잡아서 같은 값)", (20, stat("skew", 0.2, "Q3-sym")), (21, 0.82), fontsize=8.5,
                    arrowprops=dict(arrowstyle="-", color="#9ca3af"))
    axs[0].text(14, 0.05, "4×4 두 가지는 거의 0", fontsize=8.5, ha="left", color="#d7301f")
    axs[1].annotate("3×3 대칭 · 3×3 전체 · 4×4 3차 ·\n4×4 2차+knot 네 개가 모두 같은 값\n(4×4로 늘려도 개선 없음)", (20, stat("flat", 0.2, "Q3-sym")), (11, 0.62), fontsize=8.5,
                    arrowprops=dict(arrowstyle="-", color="#9ca3af"))
    fig.tight_layout(); fig.savefig(FIG / "b_cn4_error_v1.png", dpi=160)

    summary = {fa: {str(r): {k: {"mean": round(stat(fa, r, k), 3), "max": round(stat(fa, r, k, fn=lambda v: np.abs(v).max()), 3),
                                  "shape_um": round(stat(fa, r, k, key="shape", fn=np.mean), 2)} for k in reps}
                    for r in (ratios if fa != "dome" else ratios[:1])} for fa in fams}
    json.dump({"summary": summary, "nvar": NVAR, "ninputs": NIN}, open(OUT / "control_net_4x4_summary.json", "w"), ensure_ascii=False, indent=1)
    for fa in fams:
        for r in (ratios if fa != "dome" else ratios[:1]):
            print(fa, r, {k: (summary[fa][str(r)][k]["mean"], summary[fa][str(r)][k]["max"]) for k in reps})
