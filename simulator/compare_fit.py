import json, numpy as np
from simulator import MODEL31_REFINED, Strain, bc18_experiment, simulate
from validate_bc18 import load_bc18_offline, load_bc18_online
r=json.load(open('tune_bc18_result.json'))['params']
off=load_bc18_offline(); on=load_bc18_online()
for lab,st in [("refined",Strain().copy(**MODEL31_REFINED)),("tuned",Strain())]:
    res=simulate(bc18_experiment(st)); tr=res.truth; cb=res.carbon_balance()
    rm={}
    for k in ['X','Xyl','Xol','EtOH','Gly']:
        d=off[['t',k]].dropna(); rm[k]=np.sqrt(np.mean((np.interp(d.t,tr.t,tr[k])-d[k])**2))
    o=on[(on.t>20)]; rm['pO2_fed']=np.sqrt(np.mean((np.interp(o.t,tr.t,tr.pO2)-o.pO2)**2))
    print(f"== {lab}: RMSE "+"  ".join(f"{k} {v:.2f}" for k,v in rm.items())+f"   xylitol C excess {cb['xylitol_excess_Cmol']:+.3f}")
    for t in [4.3,7.65,23.4,29.1,45.6,73.2,97.7,120,143]:
        s=tr.iloc[(tr.t-t).abs().argmin()]; d=off.iloc[(off.t-t).abs().argmin()]; oo=on.iloc[(on.t-t).abs().argmin()]
        print(f" t={t:6.1f} X {s.X:5.1f}/{d.X:5.1f}  Xyl {s.Xyl:5.1f}/{d.Xyl:5.1f}  Xol {s.Xol:6.1f}/{d.Xol:6.1f}  EtOH {s.EtOH:4.1f}/{d.EtOH:4.1f}  pO2 {s.pO2:5.1f}/{oo.pO2:5.1f}  OUR {s.OUR:4.1f}")
