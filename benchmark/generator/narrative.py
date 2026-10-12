"""
Template-based prose: lab notes (4 personas), meeting minutes, weekly reports, deviation reports,
SOPs, certificates of analysis, strain construction records, maintenance log.

All numbers quoted in prose are taken from the MEASURED data summaries (never from the truth).
"""
from __future__ import annotations


import numpy as np

import measure as ms
import world as W

N = W.SYNTHETIC_NOTICE


def c(x, nd=1, comma=False):
    s = f"{x:.{nd}f}"
    return s.replace(".", ",") if comma else s


def header(title, author, date):
    return f"# {title}\n\n_{N}_\n\n**Author:** {W.PEOPLE[author].name}  \n**Date:** {date}\n\n"


# ---------------------------------------------------------------------------
# Lab notes
# ---------------------------------------------------------------------------
def note_reactor(rec, sm, rng) -> str:
    spec = rec["spec"]
    who = spec.operator
    strain = f"{spec.strain_id} c{rec['clone'] + 1}" if rec["clone"] is not None else spec.strain_id
    fault_codes = {f["code"] for f in spec.faults}
    if who == "EV":
        t = header(f"{spec.run_id} - {spec.purpose}", who, spec.start.strftime("%Y-%m-%d"))
        t += (f"**Vessel:** {spec.vessel} ({W.VESSELS[spec.vessel].model})  \n**Strain:** {strain}  \n"
              f"**Inoculated:** {spec.start:%Y-%m-%d %H:%M}\n\n## Setup\n"
              f"300 mL {W.YP} with 20 g/L glucose, 30 C, pH {spec.pH} (base/acid), {spec.stirrer:.0f} rpm, 0.7 L/min air. "
              f"Inoculated to OD600 {sm['od0']:.1f}. Feed {sm['feed_name']} started at {spec.feed_start_h:.1f} h, "
              f"{spec.feed_rate_mLh:.1f} mL/h.\n\n## Course\n")
        t += (f"- Glucose below detection by the {sm['glc_out_t']:.0f} h sample; ethanol peaked at "
              f"{sm['etoh_peak']:.1f} g/L ({sm['etoh_peak_t']:.0f} h) and was re-consumed by {sm['etoh_gone_t']:.0f} h.\n")
        if sm.get("po2_min") is not None:
            t += (f"- pO2 reached its minimum of {sm['po2_min']:.0f} % around {sm['po2_min_t']:.0f} h and sat near "
                  f"{sm['po2_fed']:.0f} % during the fed phase.\n")
        if sm.get("etoh_late"):
            t += (f"- Ethanol came back in the late fed phase (first seen at {sm['etoh_late_t']:.0f} h, "
                  f"{sm['etoh_late']:.1f} g/L by the end).\n")
        t += f"- Base used: {sm['base_mL']:.1f} mL; acid: {sm['acid_mL']:.1f} mL.\n"
        if sm.get("rerun_samples"):
            t += (f"- Update 2025-04-09: xylitol for {', '.join(sm['rerun_samples'])} replaced by the HPLC re-run after "
                  "DEV-2025-02 (original values were high).\n")

        if "F-FEEDCAL" in fault_codes:
            t += ("- Feed balance shows noticeably more feed consumed than the pump counter implies "
                  f"(about {sm['feed_bal_mL']:.0f} mL by balance vs {sm['feed_cnt_mL']:.0f} mL by pump). Pump calibration to be checked.\n")
        t += (f"\n## Result\nFinal sample ({sm['t_end']:.0f} h): xylitol {sm['xol_end']:.1f} g/L, residual xylose "
              f"{sm['xyl_end']:.1f} g/L, DW {sm['dw_end']:.1f} g/L (max {sm['dw_max']:.1f}).")
        if sm.get("for_end") is not None and spec.lot:
            t += f" Formate {sm['for_end']:.2f} g/L at the end."
        t += "\n"
        return t
    if who == "MD":
        cm = True
        lines = [f"# {spec.run_id.replace('FB-0', 'FB')} {spec.vessel} {spec.strain_id}", "", f"_{N}_", "",
                 f"- {spec.start:%d/%m/%y} inoc {spec.start:%H:%M}, OD0 {c(sm['od0'], 1, cm)}"]
        if "F-INOC" in fault_codes:
            lines.append("- preculture was 2 days old, OD lower than usual - slow start expected")
        lines.append(f"- feed {sm['feed_name']} on at {c(spec.feed_start_h, 1, cm)}h"
                     + (f", lot {spec.lot}" if spec.lot else ""))
        lines.append(f"- EtOH max {c(sm['etoh_peak'], 1, cm)} g/l")
        if "F-POWER" in fault_codes:
            f = next(x for x in spec.faults if x["code"] == "F-POWER")
            lines.append(f"- power cut ~{c(f['t0'], 0)}h (UPS tripped). feed pump + logger off, restarted "
                         f"{c(f['t1'] - f['t0'], 1, cm)}h later. see DEV")
        if "F-ANTIFOAM" in fault_codes:
            f = next(x for x in spec.faults if x["code"] == "F-ANTIFOAM")
            lines.append(f"- {c(f['t'], 0)}h lots of foam, added AF ({c(f['ppm'], 0)} ppm)")
        lines.append(f"- end {c(sm['t_end'], 0)}h: Xol {c(sm['xol_end'], 1, cm)} / Xyl {c(sm['xyl_end'], 1, cm)} / "
                     f"DW {c(sm['dw_end'], 1, cm)}" + (f" / For {c(sm['for_end'], 2, cm)}" if sm.get('for_end') is not None and spec.lot else ""))
        if sm.get("dead_t") is not None:
            lines.append(f"- culture stopped working by ~{c(sm['dead_t'], 0)}h: xylose piling up, xylitol not going up, "
                         f"DW down to {c(sm['dw_end'], 1, cm)}"
                         + (", formate building up" if (sm.get('for_end') or 0) > 1.0 else ""))
        return "\n".join(lines) + "\n"
    # PN
    t = header(f"{spec.run_id}: {spec.purpose}", who, spec.start.strftime("%d %b %Y"))
    t += (f"Ran {strain} in {spec.vessel}. ")
    if spec.lot:
        t += (f"The feed was {int(spec.frac * 100)} % of the C5 fraction (lot {spec.lot}), glucose topped up to 30 g/L and "
              "xylose to 140 g/L. ")
    if spec.enzyme_h:
        t += f"Added the hemicellulase cocktail (5 U/L) at {spec.enzyme_h:.0f} h to see whether we can release xylose from the oligomers. "
    t += (f"\n\nThe batch phase looked normal - glucose was gone by about {sm['glc_out_t']:.0f} h and ethanol went up to "
          f"{sm['etoh_peak']:.1f} g/L. ")
    if sm.get("po2_min") is not None:
        if sm["etoh_peak_t"] <= sm["po2_min_t"] <= sm["etoh_gone_t"] + 2:
            t += f"pO2 dropped to about {sm['po2_min']:.0f} % while the ethanol was being used up, which I think is the respiratory phase. "
        else:
            t += f"pO2 dipped to about {sm['po2_min']:.0f} % at around {sm['po2_min_t']:.0f} h. "
    t += (f"\n\nAt the end ({sm['t_end']:.0f} h) xylitol was {sm['xol_end']:.1f} g/L with {sm['xyl_end']:.1f} g/L xylose left. ")
    if spec.enzyme_h:
        t += "Not sure the enzyme did much - need to compare with the run without it before saying anything. "
    t += "\n"
    return t


def note_flask_screen(phase, runs, sums, decisions, rng, scatter=None, hour=None) -> str:
    title = "Clone screen in shake flasks" if phase == "P2" else "C5 tolerance and FDH clone screen"
    date = runs[0]["spec"].start.strftime("%d %b %Y")
    t = header(f"{title} ({phase})", "PN", date)
    if phase == "P2":
        t += ("50 mL YP + 20 g/L glucose + 50 g/L xylose in 250 mL baffled flasks, CaCO3 buffered, 30 C, 200 rpm, "
              "inoculated to OD 2. Each clone in duplicate, parent in duplicate.\n\n")
    else:
        t += ("Same flask set-up, but with the C5 fraction (lot C5-2501) mixed into the medium at 30, 50 or 70 %, "
              "plus 20 g/L glucose. Two flasks without CaCO3 to see how much buffering matters.\n\n")
    t += f"Xylitol (g/L, HPLC) at {hour:.0f} h, mean of the flasks:\n\n| Strain / clone | Xylitol | n |\n|---|---|---|\n"
    for k, (m, n) in sums.items():
        t += f"| {k} | {'n.d.' if m is None else f'{m:.1f}'} | {n} |\n"
    t += "\n" + ("n.d. = below the HPLC detection limit.\n\n" if any(m is None for m, _ in sums.values()) else "")
    for d in decisions:
        t += f"Picked **{d['strain']} {d['chosen']}** for further work (highest xylitol in the screen). "
    if scatter is not None:
        between, within = scatter
        if between < 2 * within:
            t += (f"\n\nDifferences between clones (sd of clone means {between:.1f} g/L) are small compared with the "
                  f"scatter between duplicate flasks (about {within:.1f} g/L), so I would not over-read them.\n")
        else:
            t += (f"\n\nClone differences (sd of clone means {between:.1f} g/L) are larger than the duplicate-flask "
                  f"scatter (about {within:.1f} g/L).\n")
    else:
        t += "\n"
    return t


# ---------------------------------------------------------------------------
# Meeting minutes (EV), weekly reports (PN)
# ---------------------------------------------------------------------------
def minutes(date, attendees, items, actions) -> str:
    t = header(f"Fermentation team meeting - {date}", "EV", date)
    t += f"**Present:** {', '.join(W.PEOPLE[a].name for a in attendees)}\n\n## Discussion\n"
    for it in items:
        t += f"- {it}\n"
    t += "\n## Decisions and actions\n"
    for a in actions:
        t += f"- {a}\n"
    return t


def weekly_report(date, paragraphs) -> str:
    t = header(f"Weekly report - week of {date}", "PN", date)
    return t + "\n\n".join(paragraphs) + "\n"


# ---------------------------------------------------------------------------
# Deviation reports
# ---------------------------------------------------------------------------
def deviation(dev_id, date, author, run, title, description, cause, impact, action) -> str:
    t = header(f"{dev_id}: {title}", author, date)
    t += (f"**Affected:** {run}\n\n## Description\n{description}\n\n## Root cause\n{cause}\n\n"
          f"## Impact on data\n{impact}\n\n## Corrective action\n{action}\n")
    return t


# ---------------------------------------------------------------------------
# SOPs
# ---------------------------------------------------------------------------
def sop_hplc(version: int) -> str:
    eff = "2024-09-02" if version == 1 else "2025-06-01"
    analytes = ms.HPLC_ANALYTES_V2 if version == 2 else ms.HPLC_ANALYTES_V1
    t = header(f"SOP-ANA-003 HPLC analysis of fermentation samples (version {version})", "TR", eff)
    t += (f"**Effective from:** {eff}\n\n## Scope\nQuantification of "
          f"{', '.join(ms.HPLC_NAMES[a].lower() for a in analytes)} in culture supernatants.\n\n"
          f"## Method\n- Column: organic acid column, 0.6 mL/min 5 mM H2SO4, 55 C.\n"
          f"- Detection: refractive index" + (" + UV 280 nm for furfural and HMF" if version == 2 else "") + ".\n"
          f"- Samples centrifuged, supernatant frozen at -20 C, analysed in batches on Tuesdays and Fridays.\n"
          f"- Dilute so that every analyte is below {ms.HPLC_MAX_VIAL:.0f} g/L in the vial (1:1 to 1:50).\n"
          f"- Five-point external calibration per sequence; check standard every 10 injections.\n"
          f"- Report values below the LOD as '<LOD' (LOD 0.02-0.05 g/L in the vial).\n")
    if version == 2:
        t += "\n## Changes from version 1\n- Succinate, furfural and HMF added. Values for these analytes before June 2025 are not available.\n"
    return t


def reactor_setup() -> str:
    t = header("Bench reactor set-up and data logging (R1-R4, GA-4)", "MD", "2025-01-10")
    p = __import__("measure").PUMP_ML_PER_COUNT
    t += ("| Vessel | Controller | Online export | Pumps logged as |\n|---|---|---|---|\n"
          "| R1, R2 | BioCtrl B1 | tab-separated text, one line every 10 min | pump counts (Acid_Pump, Base_Pump, Feed_Pump) |\n"
          "| R3, R4 | VesselView VV-1 | csv, one line every 10 min | BASE_mL, ACID_mL in mL; feed only via FEED_BAL_g |\n\n"
          "## Pump calibration (BioCtrl, checked 2025-01-10)\n"
          f"- Acid pump: {p['acid']:.2f} mL per count (1 M H2SO4)\n- Base pump: {p['base']:.2f} mL per count (2 M KOH)\n"
          f"- Feed pump: {p['feed']:.2f} mL per count\n\n"
          "## Time stamps\n"
          "- BioCtrl: the header gives 'Log start' and 'Inoculation'; Date/Time are wall-clock (dd.mm.yyyy).\n"
          "- VesselView: Elapsed_h counts from the start of the recipe (when the logger is started), not from inoculation. "
          "The inoculation time is in the run sheets and the LIMS register (the S01 sample is taken ~6 min after inoculation).\n"
          "- Feed: VesselView logs the feed bottle mass (FEED_BAL_g). BioCtrl does not log the balance; the feed bottle is "
          "weighed at the start and end of the run and written in the run workbook.\n\n"
          "## Off-gas\n"
          "GasLab GA-4 is connected to the R1 and R2 exhaust lines only (one stream file per run). Cycle every 30 min; "
          "dry gas, N2/O2/Ar/CO2 in %. Inlet air: 20.95 % O2, 0.04 % CO2, 0.7 L/min. RelTime[s] counts from the analyser "
          "start, which is not the inoculation time.\n")
    return t


def sop_fedbatch() -> str:
    t = header("SOP-FB-001 Xylitol fed-batch in 1 L bench reactors (version 3)", "EV", "2025-01-15")
    t += ("## Batch phase\n300 mL YP + 20 g/L glucose, 30 C, pH 5.5 controlled with 2 M KOH / 1 M H2SO4, 900 rpm, 0.7 L/min air. "
          "Inoculate from an overnight YPD preculture to OD600 2.\n\n## Fed phase\nStart the feed at 7.5 h (normally just after glucose "
          "depletion) at 2.0 mL/h. Feed A: 100 g/L glucose + 300 g/L xylose. C5 feeds: C5 fraction diluted with water, glucose "
          "supplemented to 30 g/L, xylose topped up to 140 g/L. Feed B (P6 only): 150 g/L glucose + 300 g/L xylose.\n\n"
          "Feed bottles stand on the feed balance; feed added = balance loss / density. Measured densities: Feed A "
          f"{W.FEED_DENSITY['A']:.2f} g/mL, Feed B {W.FEED_DENSITY['B']:.2f} g/mL, C5 feeds {W.FEED_DENSITY['C5']:.2f} g/mL.\n\n## Sampling\n7 mL at 0, 4, 8, 23.5, 30, 47.5, 72, 96, 120 and 143.5 h "
          "(167.5 h for C5 runs). OD600, dry weight (3 x 2 mL) and HPLC (SOP-ANA-003).\n\n## Antifoam\nOnly if foam reaches the "
          "headplate; record the amount.\n")
    return t


# ---------------------------------------------------------------------------
# Certificates of analysis (supplier measurement, own noise)
# ---------------------------------------------------------------------------
def coa(lot, rng, store=None) -> str:
    L = W.LOTS[lot]
    t = f"# Certificate of analysis - C5 sugar fraction, lot {lot}\n\n_{N}_\n\n**Supplier:** {L['supplier']}  \n" \
        f"**Date shipped:** {L['received']}\n\n| Component | g/L |\n|---|---|\n"
    names = dict(Xyl="Xylose (monomeric)", Glc="Glucose", Ara="Arabinose", Olig="Xylo-oligosaccharides (xylose eq.)",
                 For="Formic acid", Ac="Acetic acid", Fur="Furfural", HMF="5-HMF")
    for k, nm in names.items():
        v = L[k] * (1 + rng.normal(0, 0.03))
        t += f"| {nm} | {v:.2f} |\n"
        if store is not None:
            store.setdefault(lot, {})[k] = round(v, 2)
    t += f"| pH | {L['pH']:.1f} |\n\nStore at 4 C. Not sterile.\n"
    return t


# ---------------------------------------------------------------------------
# Strain construction records
# ---------------------------------------------------------------------------
def construction_record(plan, sdata, decision) -> str:
    t = header(f"Construction record {plan.strain_id}", plan.constructed_by, plan.date)
    t += f"**Parent:** {sdata['built_on'] or 'lab host'}\n\n**Change:** {plan.description}\n\n"
    if plan.mods:
        t += "**Parts:**\n" + "".join(f"- {m.label} ({m.kind})\n" for m in plan.mods) + "\n"
    if plan.fdh_nominal:
        t += "**Parts:**\n- pRS-FDH, heterologous formate dehydrogenase, TEF1 promoter\n\n"
    n = len(sdata["clones"])
    if n > 1:
        t += (f"Transformed, {n} clones confirmed by colony PCR and kept as c1-c{n}. ")
        if decision:
            t += f"Clone {decision['chosen']} selected after the flask screen on {decision['date']}.\n"
        else:
            t += "Clone selection pending.\n"
    if plan.strain_id in ("YX-003", "YX-004"):
        t += ("\nNote: built on YX-002 c1, which was the first PCR-confirmed clone; the YX-002 clone screen had not been "
              "run at that point.\n")
    return t
