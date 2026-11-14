import numpy as np, sys
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1]))
from solver import load_params, solve_fields, tim_map, heat_flux_map, input_names
from geometry import grid_coords
from numpy.polynomial import legendre as Lg
p=load_params(); f=p["fixed"]; n=p["solver"]["grid"]; L=f["die_size"]
names=input_names(p); rng=np.random.default_rng(0)
_,X,Y=grid_coords(L,n); u=2*X/L-1; v=2*Y/L-1
def P(m,z): c=np.zeros(m+1); c[m]=1; return Lg.legval(z,c)
modes={}
for m in (2,4,6,8):
    s=P(m,u)+P(m,v); s=(s-s.min())/(s.max()-s.min()); modes[m]=s  # 1 μm p-v when ×1e-6
res={m:[] for m in modes}
for i in range(40):
    x=dict(w_center=rng.uniform(10,40)*1e-6,w_mid_x=rng.uniform(0,20)*1e-6,w_mid_y=rng.uniform(0,20)*1e-6,
           tilt_x=0,tilt_y=0,k_tim=rng.uniform(1.5,6),hs_x=rng.uniform(2e-3,23.6e-3),hs_y=rng.uniform(2e-3,23.6e-3),hs_q=rng.uniform(80e4,240e4))
    xv=np.array([x[k] for k in names]); t0=tim_map(xv,p); q=heat_flux_map(xv,p)
    T0=solve_fields(t0,q,x["k_tim"],p).max()
    for m,s in modes.items():
        T1=solve_fields(t0+2e-6*s,q,x["k_tim"],p).max()   # 2 μm p-v 모드 추가
        res[m].append(T1-T0)
for m in res:
    a=np.abs(res[m]); print(f"P{m}: mean {a.mean():.3f} K, max {a.max():.3f} K  (per 2um p-v)")
