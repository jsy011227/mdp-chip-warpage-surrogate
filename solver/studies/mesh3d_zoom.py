"""전체 3D 격자를 짰을 때의 형상 확대 그림 → figures/b_mesh3d_zoom_v1.png
왜 데이터 생성에 3D 해석을 쓰지 않는지 보여 준다. 치수는 모두 실제 비율."""
from _common import plt, FIG
import numpy as np
from matplotlib.patches import Polygon, Rectangle, ConnectionPatch

K = (0, 0, 0)
SI, TIM, CU = (0.55, 0.60, 0.68), (0.85, 0.45, 0.10), (0.80, 0.55, 0.35)

# 실제 치수 [mm]
L, H_SI, H_CU = 25.6, 0.78, 2.0          # 냉각판 바닥 두께 2 mm는 그림용 가정
T_MIN, WARP = 0.040, 0.045                # TIM 최소 두께, 휨 (중앙 → 모서리)


def tim_t(x):                              # 돔 휨: 중앙 얇고 가장자리 두꺼움
    u = (x - L / 2) / (L / 2)
    return T_MIN + WARP * u ** 2


fig = plt.figure(figsize=(15, 11.2), facecolor="white")
gs = fig.add_gridspec(3, 3, height_ratios=[0.55, 2.2, 1.25], width_ratios=[1.35, 1, 1], hspace=.32, wspace=.18,
                      left=.04, right=.98, top=.895, bottom=.03)
fig.text(.04, .975, "전체 3D 격자로 짜면: 실제 비율로 본 형상과 격자", fontsize=18, fontweight="bold", va="top")
fig.text(.04, .952, "치수는 모두 실제 비율 (냉각판 바닥 두께 2 mm는 그림용 가정). 휨 45 μm, TIM 40~85 μm.",
         fontsize=10.5, color=(0, 0, 0, .6), va="top")

# ── (a) 칩 전체 단면, 실제 비율 ─────────────────────────────────────────
ax = fig.add_subplot(gs[0, :])
x = np.linspace(0, L, 400); t = tim_t(x)
ax.fill_between(x, 0, H_SI, color=SI, alpha=.35, lw=0)
ax.fill_between(x, H_SI, H_SI + t, color=TIM, alpha=.9, lw=0)
ax.fill_between(x, H_SI + t, H_SI + T_MIN + WARP + H_CU * .35, color=CU, alpha=.35, lw=0)
ax.set_xlim(-3.2, L + 3.6); ax.set_ylim(-.15, H_SI + .9); ax.set_aspect("equal"); ax.axis("off")
ax.text(-.25, H_SI / 2, "실리콘 0.78 mm", ha="right", va="center", fontsize=10)
ax.text(L + .35, H_SI + .07, "TIM (선 한 줄로 보임)", ha="left", va="center", fontsize=10, color=TIM)
ax.text(L + .35, H_SI + .5, "냉각판", ha="left", va="center", fontsize=10)
ax.text(L / 2, -.12, "← 칩 폭 25.6 mm →", ha="center", va="top", fontsize=10)
zx0, zx1 = L - 1.6, L
ax.add_patch(Rectangle((zx0, -.05), zx1 - zx0, H_SI + .7, fill=False, ec=K, lw=1.4, ls="--"))
ax.text(zx0 - .1, H_SI + .72, "(b) 확대 위치", ha="right", fontsize=10)
ax.set_title("(a) 칩 전체 단면 — 실제 비율이면 TIM은 거의 보이지 않는다", fontsize=12, fontweight="bold", loc="left")

# ── (b) 가장자리 1.6 × 1.6 mm 조각의 3D 격자 (cabinet projection) ───────
ax = fig.add_subplot(gs[1, :2])
ca, sa = .42, .30
P = lambda X, Y, Z: (X + ca * Y, Z + sa * Y)
W = D = 1.6; dxy = 0.1
nxy = int(round(W / dxy))
zs_si = np.linspace(0, H_SI, 9)                       # 실리콘 8층
cu_top = H_SI + 0.08 + 0.55                           # 냉각판은 0.55 mm만 보여 줌

def t_at(X):                                          # 조각 안 TIM 두께 (실제 x 위치)
    return tim_t(L - W + X)

def layer_z(X, k, nk=4):                              # TIM k번째 경계면 높이
    return H_SI + t_at(X) * k / nk

def face(pts, fc, a, ec_a=.55, lw=.6):
    ax.add_patch(Polygon([P(*p) for p in pts], closed=True, fc=(*fc, a), ec=(*K, ec_a), lw=lw))

# 앞면 (y = 0): 실리콘 → TIM → 구리, 격자 칸마다
Xs = np.linspace(0, W, nxy + 1)
for i in range(nxy):
    xa, xb = Xs[i], Xs[i + 1]
    for k in range(8):                                 # 실리콘
        face([(xa, 0, zs_si[k]), (xb, 0, zs_si[k]), (xb, 0, zs_si[k + 1]), (xa, 0, zs_si[k + 1])], SI, .22)
    for k in range(4):                                 # TIM 4층
        face([(xa, 0, layer_z(xa, k)), (xb, 0, layer_z(xb, k)), (xb, 0, layer_z(xb, k + 1)), (xa, 0, layer_z(xa, k + 1))],
             TIM, .75, ec_a=.7, lw=.35)
    zc = [layer_z(xa, 4), layer_z(xb, 4)]
    for k in range(5):                                 # 구리 5층
        f0, f1 = k / 5, (k + 1) / 5
        face([(xa, 0, zc[0] + (cu_top - zc[0]) * f0), (xb, 0, zc[1] + (cu_top - zc[1]) * f0),
              (xb, 0, zc[1] + (cu_top - zc[1]) * f1), (xa, 0, zc[0] + (cu_top - zc[0]) * f1)], CU, .18)
# 오른쪽 면 (x = W): 칩 가장자리
Ys = np.linspace(0, D, nxy + 1)
zr = [*zs_si, *[layer_z(W, k) for k in range(1, 5)]]
zr += [layer_z(W, 4) + (cu_top - layer_z(W, 4)) * f for f in np.linspace(.2, 1, 5)]
cols = [SI] * 8 + [TIM] * 4 + [CU] * 5; alph = [.14] * 8 + [.7] * 4 + [.12] * 5
for j in range(nxy):
    for k in range(len(zr) - 1):
        face([(W, Ys[j], zr[k]), (W, Ys[j + 1], zr[k]), (W, Ys[j + 1], zr[k + 1]), (W, Ys[j], zr[k + 1])],
             cols[k], alph[k], ec_a=.45 if cols[k] is not TIM else .6, lw=.35 if cols[k] is TIM else .5)
# 윗면 (구리 꼭대기)
for i in range(nxy):
    for j in range(nxy):
        face([(Xs[i], Ys[j], cu_top), (Xs[i + 1], Ys[j], cu_top), (Xs[i + 1], Ys[j + 1], cu_top), (Xs[i], Ys[j + 1], cu_top)],
             CU, .10, ec_a=.35, lw=.4)
ax.set_xlim(-.55, W + ca * D + 1.25); ax.set_ylim(-.12, cu_top + sa * D + .1); ax.set_aspect("equal"); ax.axis("off")
lx = W + ca * D + .08
ax.annotate("냉각판 (구리)\n격자 0.1 mm", xy=P(W, D * .5, cu_top - .2), xytext=(lx, cu_top + .15), fontsize=10.5,
            arrowprops=dict(arrowstyle="-", color=K, alpha=.5))
ax.annotate("TIM 4층\n한 층 10~21 μm", xy=P(W, D * .55, H_SI + .04), xytext=(lx, H_SI + .05), fontsize=10.5, color=TIM,
            fontweight="bold", arrowprops=dict(arrowstyle="-", color=TIM))
ax.annotate("실리콘 8층\n한 층 0.1 mm", xy=P(W, D * .5, .35), xytext=(lx, .2), fontsize=10.5,
            arrowprops=dict(arrowstyle="-", color=K, alpha=.5))
ax.text(W / 2, -.1, "← 1.6 mm (가장자리 조각, 축척 그대로) →", ha="center", va="top", fontsize=10)
ax.set_title("(b) 가장자리 1.6 × 1.6 mm 조각의 3D 격자 — 이 조각만 16×16×17 = 4,352칸", fontsize=12, fontweight="bold", loc="left")
zoom_anchor = P(W * .5, 0, H_SI + .03)
ax.add_patch(Rectangle((zoom_anchor[0] - .18, zoom_anchor[1] - .09), .36, .18, fill=False, ec=K, lw=1.3, ls="--"))
ax_b = ax

# ── (c) TIM 부근 한 번 더 확대: 칸 모양 비교 ─────────────────────────────
ax = fig.add_subplot(gs[1, 2])
ax.set_xlim(-.02, .40); ax.set_ylim(.66, 1.08); ax.set_aspect("equal"); ax.axis("off")
xx = np.array([0, .1, .2, .3]); X0 = W * .5 - .15
for i in range(3):
    xa, xb = xx[i], xx[i + 1]
    for k in (6, 7):
        ax.add_patch(Polygon([(xa, zs_si[k]), (xb, zs_si[k]), (xb, zs_si[k + 1]), (xa, zs_si[k + 1])], fc=(*SI, .25), ec=K, lw=.8))
    for k in range(4):
        ax.add_patch(Polygon([(xa, layer_z(X0 + xa, k)), (xb, layer_z(X0 + xb, k)), (xb, layer_z(X0 + xb, k + 1)),
                              (xa, layer_z(X0 + xa, k + 1))], fc=(*TIM, .75), ec=K, lw=.6))
    z4a, z4b = layer_z(X0 + xa, 4), layer_z(X0 + xb, 4)
    ax.add_patch(Polygon([(xa, z4a), (xb, z4b), (xb, z4b + .1), (xa, z4a + .1)], fc=(*CU, .25), ec=K, lw=.8))
ax.annotate("", (.305, H_SI), (.305, layer_z(X0 + .3, 4)), arrowprops=dict(arrowstyle="<->", lw=1))
ax.text(.31, H_SI + .03, f"TIM\n{t_at(X0 + .3) * 1e3:.0f} μm", fontsize=10, va="center", color=TIM, fontweight="bold")
ax.text(.15, .668, "← 0.3 mm →", ha="center", fontsize=10)
ax.text(0, 1.075, "TIM 칸: 100 × 100 × 15 μm\n→ 납작한 칸 (가로:세로 약 7:1)\n실리콘 칸: 100 × 100 × 98 μm", fontsize=10, va="top")
ax.set_title("(c) TIM 부근 한 번 더 확대", fontsize=12, fontweight="bold", loc="left")
fig.add_artist(ConnectionPatch((zoom_anchor[0] + .18, zoom_anchor[1]), (-.02, .80), "data", "data", axesA=ax_b, axesB=ax,
                               color=K, lw=1, ls="--", alpha=.6))

# ── (d) 비교표 ──────────────────────────────────────────────────────────
ax = fig.add_subplot(gs[2, :]); ax.axis("off")
rows = [["", "전체 3D 격자 (TIM 해상)", "2.5D (채택)"],
        ["칸 수", "256 × 256 × 20층 ≈ 131만 칸  (TIM에 맞춘 15 μm 정육면체면 약 5.6억 칸)", "64 × 64 × 2층 = 8,192 절점"],
        ["휨이 바뀔 때", "TIM 경계면이 휘므로 표본마다 형상·격자를 다시 만든다", "TIM 두께 지도만 바꾼다 (격자 그대로)"],
        ["1회 계산", "격자 생성 + 풀이, 분 단위 (추정)", "약 40 ms"],
        ["데이터 2,500개", "수 시간~수 일 (추정)", "약 2분"],
        ["정확도", "냉각판 바닥 퍼짐까지 반영 (TIM은 2.5D와 같은 결과)", "실리콘 2층 오차 최대 0.44 K, 냉각판은 보수적"]]
tb = ax.table(cellText=rows, colWidths=[.13, .52, .35], loc="upper left", cellLoc="left")
tb.auto_set_font_size(False); tb.set_fontsize(10.5); tb.scale(1, 2.0)
for (r, c), cell in tb.get_celld().items():
    cell.set_edgecolor((0, 0, 0, .25))
    if r == 0 or c == 0: cell.set_text_props(fontweight="bold"); cell.set_facecolor((.95, .95, .95))
    if c == 2 and r > 0: cell.set_facecolor((.93, .97, .93))
ax.set_title("(d) 그래서 데이터 생성은 2.5D, 3D는 일부 표본의 검증에만 쓴다 (C 담당)", fontsize=12, fontweight="bold", loc="left")

fig.savefig(FIG / "b_mesh3d_zoom_v1.png", dpi=150, facecolor="white")
print("saved")
