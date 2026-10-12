"""Bounded re-tuning of the model-31 parameters (max +/-20 % from BC18 refined) against BC18,
with the extended simulator. Writes tune_bc18_result.json."""
import json, time, dataclasses as dc
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from simulator import MODEL31_REFINED, Strain, bc18_experiment, simulate
from validate_bc18 import load_bc18_offline, load_bc18_online

P31 = ["mumax1", "mumax2", "YXS1", "YPS2", "KS1", "KP", "kd", "KSI2", "KPI2"]
ref = Strain().copy(**MODEL31_REFINED)
x_ref = np.array([MODEL31_REFINED[p] for p in P31])
off = load_bc18_offline()
on = load_bc18_online()
on_pts = on[(on.t > 20) & (on.t < 143)].iloc[::60]          # ~every 10 h, fed phase only
SCALE = dict(X=15.0, Glc=20.0, Xyl=30.0, Xol=150.0, EtOH=6.0, Gly=1.2)
W = dict(X=3.0, Glc=0.5, Xyl=1.0, Xol=1.0, EtOH=0.7, Gly=0.5, pO2=0.5)

def residuals(f, return_sim=False):
    st = ref.copy(**{p: v for p, v in zip(P31, x_ref * f)})
    ex = dc.replace(bc18_experiment(st))
    tr = simulate(ex, dt_out=0.1, rtol=1e-5, atol=1e-7).truth
    res = []
    for k, sc in SCALE.items():
        d = off[["t", k]].dropna()
        sim = np.interp(d.t, tr.t, tr[k])
        res.append(W[k] * (sim - d[k].values) / sc)
    res.append(W["pO2"] * (np.interp(on_pts.t, tr.t, tr.pO2) - on_pts.pO2.values) / 100.0)
    r = np.concatenate(res)
    return (r, tr) if return_sim else r

t0 = time.time()
r0 = residuals(np.ones(9)); print("start cost %.4f" % (0.5 * r0 @ r0))
sol = least_squares(residuals, np.ones(9), bounds=(0.8 * np.ones(9), 1.2 * np.ones(9)),
                    diff_step=0.02, x_scale=0.1, max_nfev=120, verbose=1)
print("end cost %.4f  (%.0f s, %d evals)" % (sol.cost, time.time() - t0, sol.nfev))
out = {p: dict(refined=float(a), tuned=float(a * f), change_pct=float(100 * (f - 1)))
       for p, a, f in zip(P31, x_ref, sol.x)}
json.dump(dict(params=out, cost_start=float(0.5 * r0 @ r0), cost_end=float(sol.cost)),
          open("tune_bc18_result.json", "w"), indent=1)
for p, v in out.items():
    print(f"{p:7s} {v['refined']:9.4f} -> {v['tuned']:9.4f}  ({v['change_pct']:+5.1f} %)")
