"""English versions of the drawing figures (string-swapped copies) → figures/b_mesh_en.png, b_cell_en.png"""
import re, runpy, tempfile
from pathlib import Path
H = Path(__file__).resolve().parent
MESH = {
 "2.5D 해석기 격자: 실리콘만 x·y·z로 쪼개고, TIM과 냉각판은 위치별 저항으로 붙인다": "2.5D grid: silicon meshed in x·y·z; TIM and cold plate = per-location resistances",
 "실리콘 1층 (아래 = 능동면 쪽)\\n64×64 셀, 두께 0.39 mm": "Si layer l=1 (active side)\\n64×64 nodes, 0.39 mm",
 "실리콘 2층\\n64×64 셀, 두께 0.39 mm": "Si layer l=2\\n64×64 nodes, 0.39 mm",
 "TIM: 격자 없음\\n점 하나 = 위치별 1D 저항": "TIM: no mesh\\ndot = 1D resistance",
 "냉각판: 열전달계수 h": "cold plate: h",
 "발열 q(x, y) → 능동면으로 들어옴": "heat flux q(x, y) in (active surface)",
 "냉각수 T_c 로 빠져나감": "out to coolant T_c",
 "(a) 층별로 떼어 놓은 모습 · 두께 과장, 그림은 8×8 (실제 64×64)": "(a) exploded layers · thickness exaggerated, 8×8 shown (64×64 used)",
 "(b) 단면 한 줄: 셀(점)과 열저항(지그재그)의 연결": "(b) one row: nodes (dots) and resistances (zigzag)",
 "냉각수 T_c": "coolant T_c", "1/h  (냉각판)": "1/h  (cold plate)",
 "t_TIM(x,y) / k_TIM\\n(휨에 따라 위치마다 다름)": "t_TIM(x,y) / k_TIM\\n(varies with warpage)",
 "실리콘 2층 셀": "Si node l=2", "실리콘 1층 셀": "Si node l=1", "발열 q(x,y)": "heat flux q(x,y)",
 "가로 연결 k·Δz/Δx": "lateral k·Δz·Δy/Δx", "세로 연결 k·Δx²/Δz": "vertical k·Δx²/Δz",
 "Tj = 1층 셀 온도 + q·(Δz/2)/k\\n(능동면까지 반 셀 보정)": "Tj = T(l=1) + q·(Δz/2)/k\\n(half-cell to active surface)",
 "(c) 크기와 이유": "(c) sizes",
 '("층", "어떻게 나눴나", "셀 크기", "이유 (축별 비오 수)")': '("layer", "discretization", "spacing", "note")',
 '"x 64 × y 64 × z 2 = 8,192개 셀"': '"64 × 64 × 2 = 8,192 nodes"',
 '"Bi_z 0.09~0.29 > 0.1 → z 분할 필요, 2층이면 층당 0.04~0.14 · 핫스팟 Bi_xy ≈ 1 → 가로 전도 필요"': '"2 layers: max error 0.44 K vs 12 layers"',
 '"격자 없음, x·y 위치마다 저항 1개 (4,096개)"': '"no mesh, 1 resistance per (x,y) (4,096)"',
 '"발열·가로 흐름 없음 → 1D 저항이 정확 · 가로 Bi 440~800 → 이웃끼리 독립"': '"no generation, t ≪ Δx → 1D plane wall"',
 '"냉각판"': '"cold plate"', '"격자 없음, 위치마다 1/h"': '"no mesh, 1/h"',
 '"구리 바닥판의 가로 퍼짐은 무시 → 온도를 약간 높게 봄 (보수적)"': '"base spreading neglected → conservative"',
 '"결과"': '"result"', '"미지수 8,192개 희소 연립방정식 1회 풀이"': '"one sparse solve, 8,192 unknowns"',
 '"1회 약 40 ms"': '"~40 ms"', '"12층 대비 평균 0.15 K / 최악 0.44 K · 격자 256 대비 평균 0.05 K"': '"mean 0.15 K vs 12 layers · 0.05 K vs grid 256"',
 "b_mesh_diagram_v1.png": "b_mesh_en.png", "\"실리콘 0.78 mm\"": "\"Si 0.78 mm\"",
}
CELL = {
 "층별 셀 하나에서 일어나는 열전달: 방향과 지배식": "Heat transfer at one node per layer: directions and equations",
 "위에서부터 냉각판 → TIM → 실리콘 윗층 → 실리콘 아랫층(능동면).  열은 아래(발열)에서 들어와 위(냉각수)로 나간다.  ": "Top to bottom: cold plate → TIM → Si l=2 → Si l=1 (active).  Heat enters at the bottom, leaves to the coolant.  ",
 "양방향 화살표 = 온도차에 따라 방향이 정해지는 교환": "double arrow = exchange set by ΔT",
 "냉각수 $T_c$": "coolant $T_c$", "TIM에서 올라옴": "from TIM", "× 옆으로 퍼지는 효과는 계산하지 않음": "× lateral spreading not modeled",
 "① 냉각판 (베이스 + 냉각수)": "① Cold plate", "모델: 위치마다 같은 열전달계수 $h$ 하나 (대류 경계)": "convection, single $h$",
 "열저항(단위면적):  1/h": "R'' = 1/h", "왜 단순화하나": "note",
 "베이스 안에서 옆으로 퍼지는 열을 무시 → 핫스팟이 덜 퍼짐": "base spreading neglected",
 "→ 온도를 약간 높게 예측 (안전 쪽, 보수적)": "→ slightly higher T (conservative)", "휨과 무관하게 모든 위치에서 같은 값": "same at all (x,y)",
 "$q''$ → 냉각판": "$q''$ → cold plate", "T가 z에 대해 직선": "linear T(z)",
 "② TIM (위치마다 독립된 1D 막대)": "② TIM (independent 1D column per (x,y))",
 "지배식: 평면벽 전도 (발열 없음, 옆 흐름 없음 → 정확)": "plane wall (no generation, no lateral flow)",
 "(두께 방향 온도는 직선)": "", "왜 격자로 자르지 않나": "why no mesh",
 "두께 40~85 μm ≪ 셀 폭 400 μm": "t = 40~85 μm ≪ Δx = 400 μm", "→ 옆 셀과 주고받는 열이 사실상 0, 1D 해가 정확": "→ lateral exchange ≈ 0",
 "★ 휨은 여기서만 들어온다: 두께 지도 $t_{TIM}(x,y)$": "★ warpage enters only here: $t_{TIM}(x,y)$",
 "$Q_z$ (l=1 절점에서)": "$Q_z$ (from l=1)", "$Q_z$ (l=2 절점으로)": "$Q_z$ (to l=2)",
 "③ 실리콘 윗층 절점 (l = 2, 냉각 쪽 경계)": "③ Si node l = 2 (convection boundary)",
 "에너지 평형:  $\\dot{E}_{in}-\\dot{E}_{out}+\\dot{E}_{gen}=0$  (모든 열을 '들어옴'으로 씀)": "energy balance $\\dot{E}_{in}-\\dot{E}_{out}+\\dot{E}_{gen}=0$ (all heat written as entering)",
 "$k\\Delta z$로 나누고 $\\Delta x=\\Delta y$,  $\\beta=(\\Delta x/\\Delta z)^2$:": "divide by $k\\Delta z$, $\\Delta x=\\Delta y$, $\\beta=(\\Delta x/\\Delta z)^2$:",
 "$T_{m\\pm1}=T_{m\\pm1,n,l}$ 등 약식.  윗면 = 대류 경계 ($T_\\infty = T_c$). ①②가 직렬로 묶인 등가 대류계수 $h_{eq}$": "$T_{m\\pm1}=T_{m\\pm1,n,l}$.  $T_\\infty = T_c$, $h_{eq}$ = ①② in series",
 "④ 실리콘 아랫층 절점 (l = 1, 능동면, 열유속 경계)": "④ Si node l = 1 (active surface, heat-flux boundary)",
 "에너지 평형: 능동면 열유속 $q''$가 들어옴 (발열 $\\dot{E}_{gen}$은 실리콘 내부에 없음)": "energy balance: $q''$ enters ($\\dot{E}_{gen}=0$ inside Si)",
 "구하는 값: 절점(셀 중심) → 능동면, 반 칸 두께 평면벽": "junction: node → active surface (half-cell plane wall)",
 "옆면 단열 절점: 없는 이웃의 항을 빼면 된다": "insulated side node: drop missing neighbor",
 "b_cell_balance_v1.png": "b_cell_en.png",
}
for src, rep in (("mesh_diagram.py", MESH), ("cell_balance_diagram.py", CELL)):
    s = (H / src).read_text()
    for a, b in rep.items():
        if a not in s: print("MISSING", src, a[:40])
        s = s.replace(a, b)
    s = re.sub(r'""".*?"""', '', s, count=1, flags=re.S)
    left = [l for l in s.splitlines() if re.search("[가-힣]", l.split("#")[0])]
    if left: print("KOREAN LEFT", src, *left[:5], sep="\n")
    t = H / ("_en_" + src); t.write_text(s); runpy.run_path(str(t), run_name="__main__"); t.unlink()
