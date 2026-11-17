"""Shared helpers for the B (solver) studies."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "solver"))
from solver import load_params  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

for _f in fm.findSystemFonts():
    if "NotoSansCJK" in _f or "NanumGothic" in _f or "malgun" in _f.lower():
        try: fm.fontManager.addfont(_f)
        except Exception: pass
_names = {f.name for f in fm.fontManager.ttflist}
for _cand in ("Noto Sans CJK JP", "Noto Sans CJK KR", "NanumGothic", "Malgun Gothic"):
    if _cand in _names:
        plt.rcParams["font.family"] = _cand; break
plt.rcParams["axes.unicode_minus"] = False

FIG = ROOT / "figures"
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)


def random_inputs(n, p, seed=0):
    """params.yaml 범위 안에서 균일 무작위 입력 n개 (기울기 제외 휨 제약은 무시)."""
    lo = np.array([d["range"][0] for d in p["inputs"]]); hi = np.array([d["range"][1] for d in p["inputs"]])
    return lo + np.random.default_rng(seed).random((n, len(lo))) * (hi - lo)


def clean(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9)
