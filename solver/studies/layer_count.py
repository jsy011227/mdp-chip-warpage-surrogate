"""실리콘을 두께 방향으로 몇 층 나눠야 하나? (기준 = 12층)

'기존 2.5D' = 처음 해석기 (옆 전도를 능동면 온도로 계산). 1층 이상 = 현재 해석기 (층 중심 온도로 계산).

무작위 조건 24개 + 비오 수가 큰 최악 조건 3개에서 Tj,max 오차, 온도 지도 최대 오차, 계산 시간을 비교한다.
실행: python solver/studies/layer_count.py
"""
import copy, csv, time
import numpy as np
from _common import load_params, random_inputs, plt, clean, FIG, OUT
from solver import solve
from solver import tim_map, heat_flux_map
import scipy.sparse as sp, scipy.sparse.linalg as spla


def legacy_2p5d(x, pp):
    """처음 해석기: k_si*t_die*lap(Tj) - (Tj - Tc)/R + q = 0 (옆 전도를 능동면 온도로 계산)."""
    f = pp['fixed']; t = tim_map(x, pp); q = heat_flux_map(x, pp); n = t.shape[0]
    R = f['die_thickness'] / f['k_si'] + t / x[5] + 1 / f['h_coldplate']; c = f['k_si'] * f['die_thickness'] / (f['die_size'] / n) ** 2
    idx = np.arange(n * n).reshape(n, n); rows, cols, vals = [], [], []; diag = 1 / R.ravel()
    for di, dj in ((0, 1), (0, -1), (1, 0), (-1, 0)):
        a = idx[max(0, -di):n - max(0, di), max(0, -dj):n - max(0, dj)].ravel(); b = idx[max(0, di):n - max(0, -di), max(0, dj):n - max(0, -dj)].ravel()
        rows.append(a); cols.append(b); vals.append(np.full(a.size, -c)); diag[a] += c
    A = sp.csr_matrix((np.concatenate(vals + [diag]), (np.concatenate(rows + [idx.ravel()]), np.concatenate(cols + [idx.ravel()]))), shape=(n * n, n * n))
    T = spla.spsolve(A.tocsc(), q.ravel() + f['t_coolant'] / R.ravel()).reshape(n, n)
    return float(T.max()), T

p = load_params(); L = p["fixed"]["die_size"]
REF = 12
layers = [0, 1, 2, 3, 4, 6]          # 0 = 기존 2.5D

cases = [(x, p) for x in random_inputs(24, p, seed=4)]
worst = []
for sig in (1.0e-3, 1.5e-3, 2.0e-3):                       # thin, good TIM (largest Bi), small hotspot
    pp = copy.deepcopy(p); pp["fixed"]["hs_sigma"] = sig
    x = np.array([0, 0, 0, 0, 0, 6.0, 0.3 * L, 0.7 * L, 240e4])
    worst.append((x, pp))

def run(group):
    out = {nz: {"dT": [], "field": [], "sec": []} for nz in layers}
    for x, pp in group:
        Tr, Fr = solve(x, pp, nz=REF)
        for nz in layers:
            t = time.time()
            Tm, F = legacy_2p5d(x, pp) if nz == 0 else solve(x, pp, nz=nz)
            out[nz]["sec"].append(time.time() - t)
            out[nz]["dT"].append(Tm - Tr); out[nz]["field"].append(np.abs(F - Fr).max())
    return {nz: {k: np.array(v) for k, v in d.items()} for nz, d in out.items()}

R, W = run(cases), run(worst)
with open(OUT / "layer_count.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["layers", "random_mean_dTmax_K", "random_maxabs_dTmax_K", "random_maxabs_field_K", "worst_maxabs_dTmax_K", "worst_maxabs_field_K", "ms_per_solve"])
    for nz in layers:
        w.writerow([nz, f"{R[nz]['dT'].mean():+.3f}", f"{np.abs(R[nz]['dT']).max():.3f}", f"{R[nz]['field'].max():.3f}",
                    f"{np.abs(W[nz]['dT']).max():.3f}", f"{W[nz]['field'].max():.3f}", f"{R[nz]['sec'].mean()*1e3:.0f}"])

fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.9))
xs = np.arange(len(layers)); lab = [("기존\n2.5D" if n == 0 else (f"{n}층\n(채택)" if n == p["solver"]["si_layers"] else f"{n}층")) for n in layers]
axs[0].bar(xs - 0.2, [np.abs(R[n]["dT"]).max() for n in layers], 0.38, color="#2f6fd6", label="무작위 24개 중 최대")
axs[0].bar(xs + 0.2, [np.abs(W[n]["dT"]).max() for n in layers], 0.38, color="#d7301f", label="최악 조건 3개 중 최대")
axs[0].set_yscale("log"); axs[0].axhline(0.5, color="#6b7280", ls=":", lw=1); axs[0].text(5.2, 0.55, "0.5 K", fontsize=9, color="#6b7280")
axs[0].set_xticks(xs, lab); axs[0].set_ylabel(f"|Tj,max − Tj,max({REF}층)| [K]", fontsize=10)
axs[0].set_title("층을 나누면 오차가 빠르게 줄어든다", fontsize=11.5, fontweight="bold", loc="left"); axs[0].legend(fontsize=9, frameon=False)
axs[1].bar(xs, [R[n]["sec"].mean() * 1e3 for n in layers], 0.5, color="#9ca3af")
for i, n in enumerate(layers):
    axs[1].text(i, R[n]["sec"].mean() * 1e3 * 1.05, f"{R[n]['sec'].mean()*1e3:.0f} ms", ha="center", fontsize=9)
axs[1].set_xticks(xs, lab); axs[1].set_ylabel("1회 계산 시간 [ms]", fontsize=10)
axs[1].set_title("대신 계산 시간이 늘어난다 (격자 64×64)", fontsize=11.5, fontweight="bold", loc="left")
for a in axs: clean(a)
fig.tight_layout(w_pad=3); fig.savefig(FIG / "b_silicon_layers_v1.png", dpi=160)

print("layers | random mean dT | random max|dT| | random max field | worst max|dT| | worst max field | ms")
for nz in layers:
    print(nz, f"{R[nz]['dT'].mean():+.3f}", f"{np.abs(R[nz]['dT']).max():.3f}", f"{R[nz]['field'].max():.3f}",
          f"{np.abs(W[nz]['dT']).max():.3f}", f"{W[nz]['field'].max():.3f}", f"{R[nz]['sec'].mean()*1e3:.0f}")
