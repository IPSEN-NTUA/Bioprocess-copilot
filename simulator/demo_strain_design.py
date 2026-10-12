"""Engineer a few strains (3 clones each), simulate them on the BC18 recipe and on an O2-limited
variant (300 rpm), and print the outcome. Usage: python demo_strain_design.py"""
import dataclasses as dc
import numpy as np, pandas as pd
from simulator import bc18_experiment, simulate
from strain_design import Modification as M, engineer, parent_record, ground_truth

DESIGNS = {
    "XYL1 2x": [M("XYL1", "overexpression", 2.0, "pTDH3-XYL1")],
    "ZWF1+ GND1+ PGI1- XYL1+": [M("ZWF1", "overexpression", 2.0, "pTEF1-ZWF1"), M("GND1", "overexpression", 2.0, "pTEF1-GND1"),
                                M("PGI1", "downregulation", 0.5, "pRNR2-PGI1"), M("XYL1", "overexpression", 2.0, "pTDH3-XYL1")],
    "GDH1 down": [M("GDH1", "downregulation", 0.5, "pREV1-GDH1")],
    "SDH1 deletion": [M("SDH1", "deletion", label="sdh1")],
    "GPD1 deletion": [M("GPD1", "deletion", label="gpd1")],
    "random mutagenesis (5 hits)": [M(None, "random_mutation") for _ in range(5)],
}
parent = parent_record()
rows = []
for cond, ex0 in [("BC18 recipe", bc18_experiment()), ("O2-limited, 300 rpm", dc.replace(bc18_experiment(), stirrer=300.0, antifoam_additions=[]))]:
    for dname, mods in [("parent H1274", [])] + list(DESIGNS.items()):
        for clone in range(1 if not mods else 3):
            rec = parent if not mods else engineer(parent, mods, f"{dname} c{clone+1}", np.random.default_rng(100 + clone))
            tr = simulate(dc.replace(ex0, strain=rec.strain)).truth
            e = tr.iloc[-1]
            rows.append(dict(condition=cond, design=dname, clone=clone + 1, genotype=rec.genotype()[:70],
                             dXol_NRA=round(100 * (np.exp(rec.dlnJ["xylitol"]) - 1), 1),
                             Xol=round(e.Xol, 1), Xyl=round(e.Xyl, 1), X=round(e.X, 1),
                             EtOH_end=round(e.EtOH, 2), Gly_end=round(e.Gly, 2), Suc_end=round(e.Suc, 3)))
df = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(df.to_string(index=False))
