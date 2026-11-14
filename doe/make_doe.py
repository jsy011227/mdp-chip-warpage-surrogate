"""[C 담당 시작 코드] LHS로 DOE 입력표를 만들고 해석기로 데이터셋을 생성한다.

예시:
    python doe/make_doe.py --n 100 --out data/data_v0      # 파일럿
    python doe/make_doe.py --n 2500 --out data/data_v1     # 본 데이터

결과 폴더:
    X.npy (n×9, SI 단위)   Tmax.npy (n)   Tfield.npy (n×grid×grid, float32)
    doe_X.csv              split.json (train/val/test 번호)   meta.json
"""
import argparse, json, sys, time, subprocess
from pathlib import Path
import numpy as np
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "solver"))
from solver import load_params, input_names, solve, warp_peak_to_valley  # noqa: E402


def lhs_inputs(n, p):
    lo = np.array([d["range"][0] for d in p["inputs"]]); hi = np.array([d["range"][1] for d in p["inputs"]])
    keep, batch, seed = [], max(2 * n, 50), p["seeds"]["doe"]
    while sum(len(k) for k in keep) < n:          # 제약을 어기는 점은 버리고 다시 뽑는다
        u = qmc.LatinHypercube(d=len(lo), optimization="random-cd", seed=seed).random(batch)
        X = lo + u * (hi - lo)
        ok = [warp_peak_to_valley(x, p) <= p["constraints"]["max_warp_peak_to_valley"] for x in X]
        keep.append(X[np.array(ok)]); seed += 1
    return np.concatenate(keep)[:n]


def split(n, p, frac=(0.7, 0.15, 0.15)):
    idx = np.random.default_rng(p["seeds"]["split"]).permutation(n)
    a, b = int(frac[0] * n), int((frac[0] + frac[1]) * n)
    return {"train": idx[:a].tolist(), "val": idx[a:b].tolist(), "test": idx[b:].tolist()}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=100); ap.add_argument("--out", default="data/data_v0")
    a = ap.parse_args(); p = load_params(); out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True)
    X = lhs_inputs(a.n, p); g = p["solver"]["grid"]
    Tmax = np.zeros(a.n); Tf = np.zeros((a.n, g, g), np.float32); t0 = time.time()
    for k, x in enumerate(X):
        Tmax[k], Tf[k] = solve(x, p)
        if (k + 1) % 100 == 0: print(f"{k + 1}/{a.n}  {time.time() - t0:.1f} s")
    np.save(out / "X.npy", X); np.save(out / "Tmax.npy", Tmax); np.save(out / "Tfield.npy", Tf)
    np.savetxt(out / "doe_X.csv", X, delimiter=",", header=",".join(input_names(p)), comments="")
    (out / "split.json").write_text(json.dumps(split(a.n, p)))
    try: rev = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception: rev = "unknown"
    meta = {"n": a.n, "grid": g, "params_version": p["version"], "git": rev, "inputs": input_names(p),
            "created": time.strftime("%Y-%m-%d %H:%M"), "seconds": round(time.time() - t0, 1)}
    (out / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2))
    print(f"done: {a.n} samples, Tmax {Tmax.min():.1f}–{Tmax.max():.1f} °C, {meta['seconds']} s → {out}")


if __name__ == "__main__":
    main()
