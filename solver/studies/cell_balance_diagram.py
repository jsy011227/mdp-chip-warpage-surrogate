"""층별 셀 하나의 열 흐름 방향과 지배식 그림 → figures/b_cell_balance_v1.png"""
from _common import plt, FIG
import numpy as np
from matplotlib.patches import Polygon, FancyArrowPatch

K = (0, 0, 0)
HEAT = (0.80, 0.33, 0.05)      # 발열
COOL = (0.10, 0.35, 0.70)      # 냉각수
ca, sa = 0.45, 0.32            # cabinet projection (y 깊이 방향)


def P(x, y, z): return np.array([x + ca * y, z + sa * y])


def box(ax, x0, z0, w, d, h, fill=.05, edge=.75, lw=1.3, color=K):
    """직육면체 (x0..x0+w, 0..d, z0..z0+h) 투명 솔리드."""
    c = [[(x0, 0, z0), (x0 + w, 0, z0), (x0 + w, 0, z0 + h), (x0, 0, z0 + h)],           # 앞
         [(x0, 0, z0 + h), (x0 + w, 0, z0 + h), (x0 + w, d, z0 + h), (x0, d, z0 + h)],   # 위
         [(x0 + w, 0, z0), (x0 + w, d, z0), (x0 + w, d, z0 + h), (x0 + w, 0, z0 + h)]]   # 오른쪽
    for f in c:
        ax.add_patch(Polygon([P(*v) for v in f], closed=True, fc=(*color, fill), ec=(*K, edge), lw=lw))
    for a, b in [((x0, d, z0), (x0, d, z0 + h)), ((x0, d, z0), (x0 + w, d, z0)), ((x0, d, z0), (x0, 0, z0))]:
        ax.plot(*zip(P(*a), P(*b)), color=K, alpha=.3, lw=1, ls=(0, (3, 2)))


def arrow(ax, a, b, color=K, both=False, lw=1.8, alpha=.85, ms=14):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="<|-|>" if both else "-|>", mutation_scale=ms,
                                 color=color, lw=lw, alpha=alpha, shrinkA=0, shrinkB=0))


def cross(ax, p, s=.13):
    ax.plot([p[0] - s, p[0] + s], [p[1] - s, p[1] + s], color=K, lw=2)
    ax.plot([p[0] - s, p[0] + s], [p[1] + s, p[1] - s], color=K, lw=2)


def zig(ax, x, z0, z1, amp=.13, k=6):
    zs = np.linspace(z0, z1, 2 * k + 3)
    xs = np.r_[x, x + amp * np.tile([1, -1], k + 1)[:2 * k + 1], x]
    ax.plot(xs, zs, color=K, lw=1.5, alpha=.8)


def setup(ax):
    ax.set_xlim(-0.2, 6.6); ax.set_ylim(-1.15, 3.25); ax.set_aspect("equal"); ax.axis("off")


def eqs(ax, title, lines):
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0, .95, title, fontsize=14, fontweight="bold", va="top")
    y = .78
    for kind, t in lines:
        if kind == "h":      # 소제목
            ax.text(0, y, t, fontsize=10.5, color=(0, 0, 0, .6), va="top"); y -= .085
        elif kind == "m":    # 수식
            ax.text(.03, y, t, fontsize=14.5, va="top"); y -= .135
        else:                # 설명
            ax.text(.03, y, t, fontsize=10.5, va="top"); y -= .085


fig = plt.figure(figsize=(15, 20.5), facecolor="white")
gs = fig.add_gridspec(4, 2, width_ratios=[1, 1.35], hspace=.08, wspace=.02, left=.02, right=.98, top=.955, bottom=.03)
fig.text(.02, .985, "층별 셀 하나에서 일어나는 열전달: 방향과 지배식", fontsize=18, fontweight="bold", va="top")
fig.text(.02, .966, "위에서부터 냉각판 → TIM → 실리콘 윗층 → 실리콘 아랫층(능동면).  열은 아래(발열)에서 들어와 위(냉각수)로 나간다.  "
         "양방향 화살표 = 온도차에 따라 방향이 정해지는 교환", fontsize=10.5, color=(0, 0, 0, .65), va="top")
W, D, H = 2.2, 2.0, 1.25
x0, z0 = 1.6, 0.0

# ── ① 냉각판 ────────────────────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 0]); setup(ax)
box(ax, x0, z0, W, D, .45, fill=.10)
xs = np.linspace(0, 1, 60)
for k in range(3):
    ax.plot(x0 + .1 + (W + 1.3) * xs, 2.0 + .25 * k + .06 * np.sin(xs * 10 * np.pi), color=COOL, lw=1.3, alpha=.6)
ax.text(x0 + W + 1.55, 2.25, "냉각수 $T_c$", fontsize=11, color=COOL, va="center")
arrow(ax, P(x0 + W / 2, D / 2, .45) , P(x0 + W / 2, D / 2, 1.85), color=COOL)
ax.text(P(x0 + W / 2, D / 2, 1.2)[0] + .15, P(x0 + W / 2, D / 2, 1.2)[1], "$q''$", fontsize=12, color=COOL)
arrow(ax, P(x0 + W / 2, D / 2, -.9), P(x0 + W / 2, D / 2, 0), lw=1.4, alpha=.6)
ax.text(P(x0 + W / 2, 0, -.6)[0] - 1.65, -.65, "TIM에서 올라옴", fontsize=10, color=(0, 0, 0, .6))
for p in (P(x0 - .45, 0, .22), P(x0 + W + .45, 0, .22)): cross(ax, p)
ax.text(x0 - .2, -.95, "× 옆으로 퍼지는 효과는 계산하지 않음", fontsize=9.5, color=(0, 0, 0, .65))
eqs(fig.add_subplot(gs[0, 1]), "① 냉각판 (베이스 + 냉각수)", [
    ("h", "모델: 위치마다 같은 열전달계수 $h$ 하나 (대류 경계)"),
    ("m", r"$q''(x,y) = h\,\left(T_{cp}(x,y) - T_c\right)$"),
    ("t", r"열저항(단위면적):  1/h = 1 / 70,000 = $1.43\times10^{-5}$ m²·K/W"),
    ("h", "왜 단순화하나"),
    ("t", "베이스 안에서 옆으로 퍼지는 열을 무시 → 핫스팟이 덜 퍼짐"),
    ("t", "→ 온도를 약간 높게 예측 (안전 쪽, 보수적)"),
    ("t", "휨과 무관하게 모든 위치에서 같은 값"),
])

# ── ② TIM ───────────────────────────────────────────────────────────────
ax = fig.add_subplot(gs[1, 0]); setup(ax)
box(ax, x0 + .7, z0, .8, .8, 1.2, fill=.06)                                    # 셀 위의 TIM 기둥
zig(ax, P(x0 + 1.1, .4, 0)[0], P(0, .4, .15)[1], P(0, .4, 1.05)[1])
arrow(ax, P(x0 + 1.1, .4, -.9), P(x0 + 1.1, .4, -.02))
arrow(ax, P(x0 + 1.1, .4, 1.22), P(x0 + 1.1, .4, 2.05))
ax.text(P(x0 + 1.1, .4, 1.6)[0] + .15, P(0, .4, 1.6)[1], "$q''$ (위로만)", fontsize=11)
for p in (P(x0 + .25, .4, .6), P(x0 + 1.95, .4, .6)): cross(ax, p)
ax.annotate("", P(x0 + .65, .4, -.05), P(x0 + .65, .4, 1.25),
            arrowprops=dict(arrowstyle="<->", color=K, alpha=.6, lw=1))
ax.text(P(x0 + .45, .4, .6)[0] - 1.15, P(0, .4, .95)[1], "$t_{TIM}(x,y)$\n40~85 μm", fontsize=10.5, ha="center")
# 오른쪽: 직선 온도 분포
gx, gz0, gz1 = 4.6, 0.05, 1.75
ax.plot([gx, gx], [gz0, gz1 + .15], color=K, lw=1); ax.plot([gx, gx + 1.6], [gz0, gz0], color=K, lw=1)
ax.plot([gx + 1.4, gx + .25], [gz0 + .1, gz1], color=K, lw=2.2, alpha=.8)
ax.text(gx - .1, gz1 + .25, "z", fontsize=10.5); ax.text(gx + 1.65, gz0 - .05, "T", fontsize=10.5, va="center")
ax.text(gx + .55, -.35, "T가 z에 대해 직선", fontsize=9.5, ha="center", color=(0, 0, 0, .65))
eqs(fig.add_subplot(gs[1, 1]), "② TIM (위치마다 독립된 1D 막대)", [
    ("h", "지배식: 발열 없음 + 옆 흐름 없음 → 1차원 전도"),
    ("m", r"$\dfrac{d}{dz}\left(k_{TIM}\dfrac{dT}{dz}\right)=0 \;\Rightarrow\; T(z)$ 는 직선"),
    ("m", r"$q''(x,y) = \dfrac{k_{TIM}}{t_{TIM}(x,y)}\,\left(T_{Si,top} - T_{cp}\right)$"),
    ("h", "왜 격자로 자르지 않나"),
    ("t", "두께 40~85 μm ≪ 셀 폭 400 μm"),
    ("t", "→ 옆 셀과 주고받는 열이 사실상 0, 1D 해가 정확"),
    ("t", "★ 휨은 여기서만 들어온다: 두께 지도 $t_{TIM}(x,y)$"),
])

# ── ③ 실리콘 윗층 셀 ─────────────────────────────────────────────────────
def si_cell(ax, top):
    box(ax, x0, z0, W, D, H, fill=.07)
    c = P(x0 + W / 2, D / 2, H / 2)
    zc = H / 2
    arrow(ax, P(x0 - 1.3, D / 2, zc), P(x0, D / 2, zc), both=True)                       # -x
    arrow(ax, P(x0 + W, D / 2, zc), P(x0 + W + 1.3, D / 2, zc), both=True)               # +x
    arrow(ax, P(x0 + W / 2, -1.3, zc), P(x0 + W / 2, 0, zc), both=True, alpha=.6)        # -y (앞)
    arrow(ax, P(x0 + W / 2, D, zc), P(x0 + W / 2, D + 1.0, zc), both=True, alpha=.6)     # +y (뒤)
    ax.text(*(P(x0 - 1.3, D / 2, zc) + [0, .15]), "$G_{xy}$ (x)", fontsize=11)
    bb = dict(fc="white", ec="none", pad=1.2, alpha=.9)
    ax.text(*(P(x0 + W + 1.0, D / 2, zc) + [0, .15]), "$G_{xy}$ (x)", fontsize=11, bbox=bb)
    ax.text(*(P(x0 + W / 2, -1.3, zc) + [-.1, -.32]), "$G_{xy}$ (y)", fontsize=10, color=(0, 0, 0, .75), bbox=bb)
    ax.text(*(P(x0 + W / 2, D + 1.0, zc) + [.12, .02]), "$G_{xy}$ (y)", fontsize=10, color=(0, 0, 0, .75), bbox=bb)
    o = np.array([-.05, 2.55])                                                          # 좌표축
    for v, t in [((.5, 0), "x"), ((ca * .7, sa * .7), "y"), ((0, .5), "z")]:
        arrow(ax, o, o + v, lw=1, ms=8, alpha=.6); ax.text(*(o + np.array(v) * 1.25), t, fontsize=9, alpha=.7)
    ax.plot(*c, "o", color=K, ms=5)
    ax.text(c[0] + .12, c[1] - .25, "$T_i$", fontsize=11.5)
    if top:
        arrow(ax, P(x0 + W / 2, D / 2, H), P(x0 + W / 2, D / 2, H + 1.0), color=COOL)
        ax.text(P(x0 + W / 2, D / 2, H + .6)[0] + .15, P(0, D / 2, H + .6)[1], "$G_{top}$ → TIM·냉각판", fontsize=11, color=COOL)
        arrow(ax, P(x0 + W / 2, D / 2, -1.05), P(x0 + W / 2, D / 2, 0))
        ax.text(P(x0 + W / 2, D / 2, -.6)[0] + .15, P(0, D / 2, -.6)[1], "$G_z$ (아랫층에서)", fontsize=11)
    else:
        arrow(ax, P(x0 + W / 2, D / 2, H), P(x0 + W / 2, D / 2, H + 1.0))
        ax.text(P(x0 + W / 2, D / 2, H + .6)[0] + .15, P(0, D / 2, H + .6)[1], "$G_z$ (윗층으로)", fontsize=11)
        for xx in np.linspace(x0 + .3, x0 + W - .3, 4):
            arrow(ax, P(xx, 0, -.75), P(xx, 0, -.03), color=HEAT, lw=1.6, ms=11)
        ax.text(x0 + W + .35, -.55, r"발열 $q(x,y)\,\Delta x^2$", fontsize=11, color=HEAT)
    ax.text(x0 + .05, H + .05 if not top else -.3, "", fontsize=9)


ax = fig.add_subplot(gs[2, 0]); setup(ax); si_cell(ax, top=True)
eqs(fig.add_subplot(gs[2, 1]), "③ 실리콘 윗층 셀 (3D 전도, 냉각 쪽)", [
    ("h", "지배식 (연속):"),
    ("m", r"$\nabla\cdot\left(k_{Si}\,\nabla T\right) = 0$"),
    ("h", "셀 하나의 열수지 (유한체적): 나가는 열의 합 = 0"),
    ("m", r"$\sum_{j\in xy} G_{xy}(T_i-T_j) + G_z(T_i-T_i^{\downarrow}) + G_{top}(T_i-T_c) = 0$"),
    ("m", r"$G_{xy}=k\,\Delta z \quad G_z=\dfrac{k\,\Delta x^2}{\Delta z} \quad G_{top}=\dfrac{\Delta x^2}{R_{top}}$"),
    ("m", r"$R_{top}=\dfrac{\Delta z/2}{k_{Si}}+\dfrac{t_{TIM}(x,y)}{k_{TIM}}+\dfrac{1}{h}$"),
    ("t", "①②의 저항이 직렬로 묶여 이 셀의 윗면 경계가 된다"),
])

# ── ④ 실리콘 아랫층 셀 (능동면) ─────────────────────────────────────────
ax = fig.add_subplot(gs[3, 0]); setup(ax); si_cell(ax, top=False)
eqs(fig.add_subplot(gs[3, 1]), "④ 실리콘 아랫층 셀 (능동면, 발열 쪽)", [
    ("h", "경계조건 (능동면, z = 0):"),
    ("m", r"$-\left.k_{Si}\,\dfrac{\partial T}{\partial z}\,\right|_{z=0} = q(x,y)$"),
    ("h", "셀 하나의 열수지: 나가는 열의 합 = 들어온 발열"),
    ("m", r"$\sum_{j\in xy} G_{xy}(T_i-T_j) + G_z(T_i-T_i^{\uparrow}) = q_i\,\Delta x^2$"),
    ("h", "구하는 값: 셀 중심에서 능동면까지 반 칸 보정"),
    ("m", r"$T_j = T_i + q_i\,\dfrac{\Delta z/2}{k_{Si}}$"),
    ("t", "옆면 4개는 단열 → 가장자리 셀은 이웃 수가 적을 뿐"),
])

for y in (.7365, .5005, .2645):
    fig.add_artist(plt.Line2D([.02, .98], [y, y], color=K, alpha=.15, lw=1))
fig.savefig(FIG / "b_cell_balance_v1.png", dpi=150, facecolor="white")
print("saved", FIG / "b_cell_balance_v1.png")
