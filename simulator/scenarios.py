import sys, time; sys.path.insert(0,'.')
from simulator import *
import dataclasses as dc
base = bc18_experiment()
def run(name, ex):
    t0=time.time(); r=simulate(ex); tr=r.truth; e=tr.iloc[-1]
    cb=r.carbon_balance()
    print(f"{name:38s} {time.time()-t0:4.1f}s X={e.X:5.1f} Glc={e.Glc:5.2f} Xyl={e.Xyl:5.1f} Xol={e.Xol:6.1f} EtOH_max={tr.EtOH.max():4.1f} For={e.For:5.2f} Fur={e.Fur:4.2f} FurOH={e.FurOH:4.2f} Trc={e.Trc:4.1f} minpO2={tr.pO2.min():5.1f} base={e.cumBase*1000:5.1f}mL acid={e.cumAcid*1000:4.1f}mL C={cb['closure_without_YPS2_excess']:.4f}")
    return r
run("BC18 reference", base)
run("xylose-only feed (no glucose)", dc.replace(base, feed_composition={"Xyl":317.0}))
run("glucose:xylose feed 1:1", dc.replace(base, feed_composition={"Glc":200.0,"Xyl":200.0}))
run("300 rpm (O2 limited)", dc.replace(base, stirrer=300.0))
run("DO cascade 30%, N 300-1200", dc.replace(base, stirrer=300.0, do_cascade={"setpoint":30,"N_min":300,"N_max":1200,"gain":40}))
run("35 degC", dc.replace(base, T_set=35.0))
run("pH 4.5", dc.replace(base, pH_set=4.5))
cimv = {"Xyl":140.0*0.7+0,"Glc":30.0,"Ara":12.0,"Olig":25.0,"For":14.0,"Ac":5.0,"Fur":0.8,"HMF":0.3,"Trc":70.0}
run("CIMV 70% feed, H1274 (no FDH)", dc.replace(base, feed_composition=cimv))
run("CIMV 70% feed, FDH strain", dc.replace(base, feed_composition=cimv, strain=Strain(name="H4031", fdh=1.0)))
run("CIMV 70% feed, FDH + pH 6.0", dc.replace(base, feed_composition=cimv, pH_set=6.0, strain=Strain(name="H4031", fdh=1.0)))
run("CIMV 70% + enzyme at 24 h, FDH", dc.replace(base, feed_composition=cimv, strain=Strain(name="H4031", fdh=1.0), enzyme_additions=[(24.0, 5.0)]))
pulses=[Pulse(t, 0.010, {"Glc":400.0*0.25,"Xyl":400.0*0.75}, "glucose+xylose pulse") for t in (8, 24, 48, 72)]
run("pulsed batch (4 pulses)", dc.replace(base, feed_segments=[], pulses=pulses))
run("exp feed mu=0.03", dc.replace(base, feed_segments=[FeedSegment(7.66,143,0.001,0.03)]))
