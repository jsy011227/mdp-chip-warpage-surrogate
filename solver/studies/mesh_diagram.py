"""2.5D(실리콘 2층) 해석기의 격자 구성 설명 그림. 흰 배경 + 반투명 검은 선.
실행: python solver/studies/mesh_diagram.py  → figures/b_mesh_diagram_v1.png
"""
import numpy as np
from _common import plt, FIG

K = (0, 0, 0)
def ln(a): return dict(color=K, alpha=a)

fig = plt.figure(figsize=(14, 9.6), facecolor="white")
fig.text(0.02, 0.965, "2.5D 해석기 격자: 실리콘만 x·y·z로 쪼개고, TIM과 냉각판은 위치별 저항으로 붙인다",
         fontsize=15, fontweight="bold")

# ---------------- (a) exploded oblique view (hand-drawn projection) ----------------
ax = fig.add_axes([0.02, 0.40, 0.50, 0.50]); ax.set_facecolor("white"); ax.set_axis_off()
n, S, ca, sa = 8, 1.0, 0.55 * np.cos(np.radians(30)), 0.55 * np.sin(np.radians(30))
def P(x, y, z): return x + ca * y, z + sa * y            # cabinet projection

def slab(z0, h, grid=True, fill=.04, edge=.6, cols=False):
    top = [P(0, 0, z0 + h), P(n, 0, z0 + h), P(n, n, z0 + h), P(0, n, z0 + h)]
    front = [P(0, 0, z0), P(n, 0, z0), P(n, 0, z0 + h), P(0, 0, z0 + h)]
    side = [P(n, 0, z0), P(n, n, z0), P(n, n, z0 + h), P(n, 0, z0 + h)]
    for poly, f in ((front, fill * 1.6), (side, fill * 2.4), (top, fill)):
        ax.add_patch(plt.Polygon(poly, closed=True, fc=(0, 0, 0, f), ec=(0, 0, 0, edge), lw=1.1))
    if grid:
        for k in range(1, n):
            ax.plot(*zip(P(k, 0, z0 + h), P(k, n, z0 + h)), lw=.7, **ln(.38))
            ax.plot(*zip(P(0, k, z0 + h), P(n, k, z0 + h)), lw=.7, **ln(.38))
            ax.plot(*zip(P(k, 0, z0), P(k, 0, z0 + h)), lw=.7, **ln(.38))
            ax.plot(*zip(P(n, k, z0), P(n, k, z0 + h)), lw=.7, **ln(.38))
    if cols:
        for i in range(n):
            for j in range(n):
                ax.plot(*P(i + .5, j + .5, z0 + h), "o", ms=2.2, color=(0, 0, 0, .55))

Z = {"si1": (0.0, 1.1), "si2": (1.7, 1.1), "tim": (4.1, 0.35), "cp": (5.8, 0.9)}
slab(*Z["si1"]); slab(*Z["si2"])
slab(*Z["tim"], grid=False, fill=.12, cols=True)
slab(*Z["cp"], grid=False, fill=.02, edge=.45)
for i in (1.5, 4.0, 6.5):
    x0, y0 = P(i, 0, -1.3); x1, y1 = P(i, 0, -0.05)
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=(0, 0, 0, .6), lw=1.3))
    x0, y0 = P(i, n, 6.8); x1, y1 = P(i, n, 7.7)
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=(0, 0, 0, .4), lw=1.1))
for (a0, b0) in ((1.1, 1.7), (2.8, 4.1), (4.45, 5.8)):                                 # dashed gaps (exploded)
    for xx in (0, n):
        ax.plot(*zip(P(xx, 0, a0), P(xx, 0, b0)), lw=.8, ls=(0, (3, 3)), **ln(.35))
xr = P(n, n, 0)[0] + .5
def lab(z, t): ax.text(xr, P(n, n / 2, z)[1], t, fontsize=10, va="center")
lab(0.55, "실리콘 1층 (아래 = 능동면 쪽)\n64×64 셀, 두께 0.39 mm")
lab(2.25, "실리콘 2층\n64×64 셀, 두께 0.39 mm")
lab(4.25, "TIM: 격자 없음\n점 하나 = 위치별 1D 저항")
lab(6.25, "냉각판: 열전달계수 h")
ax.text(P(4.0, 0, -1.9)[0], P(4.0, 0, -1.9)[1], "발열 q(x, y) → 능동면으로 들어옴", ha="center", fontsize=10)
ax.text(P(4.0, n, 8.1)[0], P(4.0, n, 8.1)[1], "냉각수 T_c 로 빠져나감", ha="center", fontsize=10)
ax.set_xlim(-0.5, xr + 7.5); ax.set_ylim(-2.4, 10.9); ax.set_aspect("equal")
fig.text(0.03, 0.905, "(a) 층별로 떼어 놓은 모습 · 두께 과장, 그림은 8×8 (실제 64×64)", fontsize=11.5, fontweight="bold")

# ---------------- (b) cross-section with resistor network ----------------
bx = fig.add_axes([0.53, 0.40, 0.45, 0.52]); bx.set_facecolor("white"); bx.set_axis_off()
fig.text(0.54, 0.905, "(b) 단면 한 줄: 셀(점)과 열저항(지그재그)의 연결", fontsize=11.5, fontweight="bold")
nc, W = 5, 1.6
y_lay = [0.8, 2.0]                       # silicon layer centres (cell height 1.2 each, from 0.2)

def zig(x0, y0, x1, y1, a=.6, amp=.09, k=6):
    t = np.linspace(0, 1, 2 * k + 3); px = x0 + (x1 - x0) * t; py = y0 + (y1 - y0) * t
    off = np.r_[0, np.tile([amp, -amp], k), amp, 0][:t.size] * np.r_[0, np.ones(t.size - 2), 0]
    nx, ny = -(y1 - y0), (x1 - x0); nn = np.hypot(nx, ny); nx, ny = nx / nn, ny / nn
    bx.plot(px + off * nx, py + off * ny, lw=1.2, **ln(a))

for i in range(nc):                      # cells
    for j, (yb, yt) in enumerate([(0.2, 1.4), (1.4, 2.6)]):
        bx.add_patch(plt.Rectangle((i * W, yb), W, yt - yb, fill=False, lw=1.0, **ln(.45)))
for i in range(nc):
    xc = i * W + W / 2
    for yc in y_lay: bx.plot(xc, yc, "o", ms=6, color=(0, 0, 0, .8))
    zig(xc, y_lay[0] + .12, xc, y_lay[1] - .12, a=.65)                        # z link
    zig(xc, y_lay[1] + .12, xc, 2.6, a=.65, amp=.06, k=3)                      # half cell to top
    bx.add_patch(plt.Rectangle((xc - .32, 2.6), .64, .55, fill=True, fc=(0, 0, 0, .10), lw=1, ec=(0, 0, 0, .55)))  # TIM resistor
    zig(xc, 3.2, xc, 3.9, a=.65, amp=.07, k=3)                                 # 1/h
    bx.annotate("", (xc, .2), (xc, -.45), arrowprops=dict(arrowstyle="-|>", color=(0, 0, 0, .6), lw=1.2))
for i in range(nc - 1):                  # lateral links
    for yc in y_lay: zig(i * W + W / 2 + .12, yc, (i + 1) * W + W / 2 - .12, yc, a=.55)
bx.plot([-.2, nc * W + .2], [3.9, 3.9], lw=2.2, **ln(.75))
# TIM thickness varies with warpage (exaggerated) : show as varying box height note
bx.text(nc * W + .35, 3.9, "냉각수 T_c", va="center", fontsize=10)
bx.text(nc * W + .35, 3.55, "1/h  (냉각판)", va="center", fontsize=10)
bx.text(nc * W + .35, 2.88, "t_TIM(x,y) / k_TIM\n(휨에 따라 위치마다 다름)", va="center", fontsize=10)
bx.text(nc * W + .35, 2.0, "실리콘 2층 셀", va="center", fontsize=10)
bx.text(nc * W + .35, 0.8, "실리콘 1층 셀", va="center", fontsize=10)
bx.text(nc * W + .35, -.3, "발열 q(x,y)", va="center", fontsize=10)
bx.text(W * 1.0, 2.22, "가로 연결 k·Δz/Δx", fontsize=8.8, ha="center", va="center", color=(0, 0, 0, .8), bbox=dict(fc="white", ec="none", pad=1.5))
bx.text(W / 2 + .2, 1.4, "세로 연결 k·Δx²/Δz", fontsize=8.8, ha="left", va="center", color=(0, 0, 0, .8), bbox=dict(fc="white", ec="none", pad=1.5))
bx.annotate("Tj = 1층 셀 온도 + q·(Δz/2)/k\n(능동면까지 반 셀 보정)", (W / 2, .55), (W * 0.1, -1.35), fontsize=9.3,
            arrowprops=dict(arrowstyle="-", color=(0, 0, 0, .4)))
bx.set_xlim(-.4, nc * W + 3.4); bx.set_ylim(-1.6, 4.3)

# ---------------- (c) numbers ----------------
cx = fig.add_axes([0.03, 0.03, 0.94, 0.33]); cx.set_axis_off()
fig.text(0.03, 0.355, "(c) 크기와 이유", fontsize=11.5, fontweight="bold")
rows = [
    ("층", "어떻게 나눴나", "셀 크기", "이유 (축별 비오 수)"),
    ("실리콘 0.78 mm", "x 64 × y 64 × z 2 = 8,192개 셀", "Δx = Δy = 0.40 mm,  Δz = 0.39 mm", "Bi_z 0.09~0.29 > 0.1 → z 분할 필요, 2층이면 층당 0.04~0.14 · 핫스팟 Bi_xy ≈ 1 → 가로 전도 필요"),
    ("TIM 40~85 μm", "격자 없음, x·y 위치마다 저항 1개 (4,096개)", "—", "발열·가로 흐름 없음 → 1D 저항이 정확 · 가로 Bi 440~800 → 이웃끼리 독립"),
    ("냉각판", "격자 없음, 위치마다 1/h", "—", "구리 바닥판의 가로 퍼짐은 무시 → 온도를 약간 높게 봄 (보수적)"),
    ("결과", "미지수 8,192개 희소 연립방정식 1회 풀이", "1회 약 40 ms", "12층 대비 평균 0.15 K / 최악 0.44 K · 격자 256 대비 평균 0.05 K"),
]
xs = [0.0, 0.13, 0.42, 0.62]
for r, row in enumerate(rows):
    y = 0.90 - r * 0.20
    for c, txt in enumerate(row):
        cx.text(xs[c], y, txt, fontsize=10.5 if r else 10.5, fontweight="bold" if (r == 0 or c == 0) else "normal",
                va="center", color=(0, 0, 0, .9), wrap=True)
    cx.plot([0, 1], [y - 0.10, y - 0.10], lw=.8 if r else 1.2, **ln(.35 if r else .6), transform=cx.transAxes)
cx.set_xlim(0, 1); cx.set_ylim(0, 1)

fig.savefig(FIG / "b_mesh_diagram_v1.png", dpi=170, facecolor="white")
print("saved")
