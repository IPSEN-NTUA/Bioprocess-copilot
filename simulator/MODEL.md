# Xylitol bioreactor simulator: model description

`simulator.py` produces the true state of a yeast xylitol cultivation, which the benchmark's document renderers (HPLC sheets, online logs, off-gas, lab notes …) will later read from. This file lists every equation and parameter, and where each value comes from.

**Core:** VTT kinetic model 31, starting from the BC18 *refined* parameters in `realdata_summary.json` and re-tuned against BC18 within ±20% of each value (`tune_bc18.py`). **Extensions:** physically motivated terms that leave model 31 as it is at the BC18 reference conditions (30 °C, pH 5.5, aerobic, glucose-limited feed, no inhibitors).

## State vector

| State | Unit | Notes |
|---|---|---|
| V | L | broth volume |
| X | g DW/L | viable biomass |
| Glc, Xyl, Xol | g/L | glucose (S1), xylose (S2), xylitol (P) |
| Ara, Olig | g/L | arabinose (not consumed), xylo-oligomers (xylose eq.) |
| EtOH, Gly, Ac, Suc, Pyr | g/L | by-products |
| For, Fur, HMF, FurOH, HMFOH | g/L | hydrolysate inhibitors and furan alcohols |
| Trc | % | tracer = % CIMV C5 fraction in the reactor |
| Enz | U/L | added hydrolytic enzyme |
| DO | mmol/L | dissolved O2 (pO2 = DO / DO\*) |
| Ind | – | induction of respiratory / gluconeogenic enzymes (0–1) |
| cum* | L, mmol, g | cumulative feed, base, acid, evaporation, O2, CO2, stripped ethanol, lysed biomass, xylitol formed |

## Model 31 (kept exactly)

```
mu1  = mumax1 · S1/(KS1+S1)
qP31 = mumax2 · S2 / ((KP + S2 + S2²/KSI2) · (1 + P/KPI2))
dX/dt = mu1·X − kd·X − D·X        dS1/dt = −mu1·X/YXS1 + D(S1f − S1)
dS2/dt = −qP31·X/YPS2 + D(S2f−S2) dP/dt  = qP31·X − D·P
```
| Parameter | BC18 refined | Tuned (default H1274) | Change |
|---|---|---|---|
| mumax1 (h⁻¹) | 0.331 | 0.395 | +19.5% |
| mumax2 (g/g/h) | 3.997 | 4.795 | +20.0% |
| YXS1 (g/g) | 0.690 | 0.769 | +11.5% |
| YPS2 (g/g) | 1.240 | 1.084 | −12.6% |
| KS1 (g/L) | 0.102 | 0.081 | −20.0% |
| KP (g/L) | 40.0 | 32.0 | −20.0% |
| kd (h⁻¹) | 0.0210 | 0.0252 | +20.0% |
| KSI2 (g/L) | 33.8 | 29.0 | −14.3% |
| KPI2 (g/L) | 11.43 | 13.70 | +19.9% |

The originals are kept in `MODEL31_REFINED`; use `Strain().copy(**MODEL31_REFINED)` to run with them. The tuning objective was the HPLC data (biomass, glucose, xylose, xylitol, ethanol, glycerol) plus fed-phase online pO2 from BC18, all normalised. Six of the nine parameters ended at or near the ±20% bound, so the data would pull further if allowed.

## Extensions

**Environment factors.** These multiply every uptake rate; each one equals 1 at the reference conditions.
- f_T, f_pH: Rosso cardinal models, T (5, 32, 40 °C) and pH (2.5, 5.0, 8.5), normalised at 30 °C and pH 5.5.
- f_inh = 1/(1+(For_u/0.03)²) · 1/(1+Ac_u/1.0) · 1/(1+Fur/1.0) · 1/(1+HMF/2.0) · (1 − EtOH/90). The subscript u means the undissociated acid, `C·1/(1+10^(pH−pKa))` (pKa formate 3.75, acetate 4.76). So pH and buffering change toxicity, as VTT saw in the flask screens.

**1. Glucose: oxidative channel (model 31), Pasteur effect, Crabtree overflow.**
- `qG31 = (mumax1/YXS1)·S1/(KS1+S1)·f_env` is the model-31 uptake.
- The respired part is `qG31·fO2`, with `fO2 = DO/(0.0015+DO)`. Growth on it is exactly `YXS1·qG31·fO2`.
- The non-respired part `qG31·(1−fO2)` (the Pasteur effect) and the overflow `qGof = 1.6·S1²/(1.5²+S1²)·f_env` go to a fermentative channel. Its C-mol split is X 0.12, EtOH 0.53, CO2 0.287, Gly 0.045, Ac 0.010, Suc 0.004, Pyr 0.004. The Hill form keeps overflow at about 0 under glucose limitation (S1 ≈ 0.01 g/L in the fed phase), so model 31 is undisturbed there.

**2. Re-consumption after glucose (diauxie).** Ethanol, glycerol, acetate and pyruvate are taken up at `q_i,max·C/(K+C)·Ind·DO/(0.04+DO)·f_env`. Their O2 affinity (0.04 mmol/L) is lower than that of glucose respiration (0.0015), so under O2 limitation glucose gets the oxygen first. `dInd/dt = (Ki²/(Ki²+S1²) − Ind)/1.5 h` gives glucose repression and the diauxic lag (Ki = 0.1 g/L). Biomass yields in C-mol per C-mol: EtOH 0.60 (with 0.05 going to acetate), Gly 0.40, Ac 0.30, Pyr 0.30. The remainder becomes CO2.

**3. Cofactor coupling of XR.** XR needs 1 NAD(P)H per xylitol. The available NAD(P)H is computed with degree-of-reduction balances: the γ of each consumed substrate minus the γ of its products.
- `A` = supply from external co-substrates (respired glucose, ethanol, glycerol, acetate, pyruvate, and formate via FDH at 2 e⁻/mol). `E_tot = A` + fermentative channel + endogenous respiration of decayed biomass (f_endo = 0.6 of kd·X).
- The XR rate is `qP = qP31 · f_env · φ · ψ`, where `φ = A_C/(A_C + 0.49 C-mmol/g/h)` and A_C is the carbon uptake of external co-substrates (XR needs active co-substrate catabolism; carbon-based so that small biomass-yield changes don't swing XR) and `ψ = (1 + (D/0.9·E_tot)^4)^−¼` is a soft cap on the reductase demand D (xylitol + furfural + HMF reduction).
- Result: no glucose in the feed means little xylitol, and the glucose:xylose feed ratio becomes a real design variable. In the BC18 fed phase, φ ≈ 0.71.

**4. Oxygen and CO2.** `qO2 = (E_tot − electrons used by reductases)/4`. CO2 is computed by carbon closure of every flux. `dDO/dt = kLa(DO* − DO) − qO2·X`, with `kLa = 88·(N/900)^1.2·(V/0.3 L)^−0.4·(vvm/2.33)^0.5 / (1 + AF/100 ppm) h⁻¹` (Labfors 1). The V term is the fall in power per volume as the broth grows at fixed stirrer speed. AF is antifoam, a dose event that is diluted by volume changes. kLa_ref = 88 h⁻¹ reproduces the end of the BC18 O2 crash at about 17 h. DO* comes from the Benson–Krause polynomial. An optional DO→stirrer P-cascade is available. Off-gas composition comes from an inert-gas balance (`Result.offgas()`).

**5. Hydrolysate (CIMV C5).**
- Formate is oxidised by native FDH (0.002 g/g/h) or heterologous FDH (0.25 g/g/h at fdh = 1).
- Furfural and HMF are reduced to their alcohols (0.10 and 0.04 g/g/h), consuming NAD(P)H in competition with XR.
- Undissociated formic acid and furfural also kill cells: `z = kd·((For_u/0.05)² + Fur/2)`, saturated as `kd_tox = z/(1 + z/0.30)` (at most 0.30 1/h). Basal decay `kd` stays as in model 31 (60 % respired, 40 % lysed). Inhibitor-killed cells go to a non-viable pool `Xd` that is not respired: it still counts in dry weight and OD and lyses slowly (0.005 1/h, half-life ~6 d). So a culture killed by formate keeps a measurable dry weight instead of dropping to zero. Without inhibitors (BC18 conditions) `Xd` stays 0 and the model is unchanged.
- Added enzyme releases xylose from oligomers at `k·Enz·Olig/(5+Olig)`, and the enzyme deactivates at 0.01 h⁻¹.

**6. Reactor physics.**
- Every liquid state follows `dC/dt = r + (F_feed·C_f − F_in·C + F_evap·C)/V`, with F_in = feed + base + acid. Water evaporates through a 4 °C condenser (about 0.25 mL/h at 0.7 L/min), and ethanol is stripped at 0.002 h⁻¹ at the reference aeration.
- **pH control:** the base and acid flows cancel the net proton production: secreted acetate, succinate and pyruvate, NH4⁺ assimilation (all of it with an NH4 N-source, 15% with YP), and protons consumed by formate oxidation. That's why FDH strains on hydrolysate need acid addition, as in the VTT sheets.
- Sampling, pulses and enzyme additions are discrete events. Feed is blocked above the working-volume limit (0.8 L).

## Validation against BC18, `validation_bc18.png`

RMSE vs BC18, refined → tuned (same extended model, final O2 setup):

| | Biomass | Xylose | Xylitol | Ethanol | Glycerol | Fed-phase pO2 |
|---|---|---|---|---|---|---|
| Refined | 1.49 g/L | 10.0 g/L | 5.7 g/L | 0.47 g/L | 0.46 g/L | 14.1 % pts |
| Tuned | 1.19 g/L | 0.94 g/L | 5.6 g/L | 0.57 g/L | 0.43 g/L | 5.4 % pts |

- **Residual xylose fixed:** 31.8 vs 32.9 g/L at 143 h (was 51). The main move was YPS2 1.24 → 1.08 g/g. Once sampling, evaporation and base are in the volume balance, that is needed to consume enough xylose.
- **Fed-phase pO2 fixed by antifoam, not kinetics.** No single kLa fits BC18. Ending the O2 crash at 17 h needs about 88 h⁻¹; about 50% pO2 in the fed phase at VTT's own measured OUR (~3.4 mmol/L/h) needs about 30 h⁻¹. The reference recipe therefore includes an **inferred** 200 ppm antifoam dose at 16.7 h, which is **not in the VTT records**. A slow pO2 decline of about 15% points over 120 h remains; that looks like probe drift and belongs in the measurement layer.
- **Still off:**
  - Biomass at 23–29 h: 12.3–12.9 vs 15.1–14.9 g/L. Model-31 parameters within ±20% can't reach it, because the gain is set by the ethanol phase.
  - The pO2 minimum during the ethanol phase is about 10% vs about 1% measured.
  - Xylitol at 45.6 h: 66 vs 79 g/L. That point looks high relative to its neighbours.
- **Carbon balance** closes to 1.00000 once the xylitol carbon implied by YPS2 = 1.084 > 1.013 g/g is set aside. The excess is 4% of input carbon, down from 11% with the refined YPS2. `Result.carbon_balance()` reports it.

## Strain design (`strain_design.py`)

**Ground truth:** kinetic model **978,2** (index 922 of the 939-model A_oxidative ORACLE ensemble), saved standalone in `strain_design/ground_truth_978_2.npz` (FCC 283×283, CCC 190×283, reference fluxes and concentrations of state 978).
- It is one of the most representative models: weighted distance of its boundary-flux control coefficients to the ensemble median is 0.35, against an ensemble median of 0.76.
- All of its top-10 effects per boundary flux have the ensemble sign.
- Its xylitol control pattern correlates 0.95 with the ensemble median (XR 0.37; ZWF1 +0.08, GND +0.04, NADK +0.03, GDH1 −0.05).

**1. Modification → enzyme activity.** The realised fold is `F_nominal · u`, with `u ~ U(0.5, 1.5)` drawn per transformant. The nominal fold must be in (1, 2] for overexpression and [0.5, 1) for downregulation. u is truncated so that:
- **the direction never flips:** overexpression ≥ 1×, downregulation ≤ 1×;
- **the expression change never exceeds 2-fold:** overexpression lands in 1–2×, downregulation in 0.5–1×.

Examples: 1.5× nominal → 1.0–2.0× (mean 1.5); 2× nominal → 1.0–2.0× (mean 1.5); 0.5× nominal → 0.5–0.75×.

Other cases:
- **Targeted point mutations** (known intent) use the same rule.
- **Deletions** are deterministic at 0.10× residual activity.
- **Random mutations** (mutagenesis, evolution, or an untargeted point mutation): the enzyme is random if not given. The fold is drawn from a loss-skewed distribution of effects: 45% neutral (0.97–1.03×), 30% mild loss (0.5–0.97×), 15% strong loss (0.05–0.5×), 5% near-null (0.005–0.05×), 5% gain (1.03–2×). No extra u is applied.

**2. Activity → fluxes (NRA, linear).** `dln e = ln(fold)`, summed along the lineage, and `dln J = FCC · dln e`. Large responses are compressed with `L·tanh(dlnJ/L)`, L = ln 3: about 1:1 below |dln J| ≈ 0.3, saturating at 3×. Children are always rebuilt from the root strain with the cumulative `dln e`, so compression is applied once.

**3. Fluxes → production-phase parameters.** Each multiplier m acts as `m**Ind`, so the batch phase is untouched:

| Kinetic flux | Macro parameter |
|---|---|
| XYLTt | `mumax2 × exp(dlnJ)` |
| biomass | `prod_mumax1_mult = exp(dlnJ_growth)` |
| biomass − GLCt1 | `prod_YXS1_mult = exp(dlnJ_growth − dlnJ_GLC)` |
| ETOHt, GLYCt, SUCCt, ACtr | `prod_div_mult` (share of O2-limited glucose fermented) + `prod_split`, solved so each by-product flux changes by exactly `exp(dlnJ)` with carbon and electrons closed |
| SUCCt | `prod_suc_mult` (oxidative succinate leak) |
| O2t, CO2t | not mapped; the macro balances compute them (use for consistency checks) |

FDH, inhibitor tolerance and the temperature optimum stay macro traits, set on `Strain` directly.

**Behaviour:** run `demo_strain_design.py` (results depend on the clone draws).

## Usage

```python
from simulator import Experiment, FeedSegment, Pulse, Strain, simulate, bc18_experiment
res = simulate(bc18_experiment())
res.truth        # 1-min grid: all states, rates, OUR/CER/RQ, pO2, stirrer, flows
res.samples      # true state at each sampling event (before sample removal)
res.events       # inoculation, feed start/stop, samples, pulses, enzyme, max volume …
res.offgas()     # O2 % / CO2 % / OUR / CER / RQ
res.carbon_balance()
```
Strain design:
```python
from strain_design import Modification as M, engineer, parent_record
import numpy as np
parent = parent_record()                                    # H1274
s1 = engineer(parent, [M("XYL1", "overexpression", 2.0, "pTDH3-XYL1")], "H1274-XR2", np.random.default_rng(1))
s2 = engineer(s1, [M("ZWF1", "overexpression", 2.0, "pTEF1-ZWF1")], "H1274-XR2-ZWF1", np.random.default_rng(2))
s2.summary(); s2.genotype()                                  # NRA flux changes, realised genotype
simulate(dc.replace(bc18_experiment(), strain=s2.strain))
```
A run takes about 1.5 s. `validate_bc18.py` regenerates both figures, `tune_bc18.py` reruns the ±20% tuning (writes `tune_bc18_result.json`), and `compare_fit.py` prints refined vs tuned errors.

## Open points for Stefanos

1. Confirm the inferred antifoam dose in BC18, or replace it with what VTT actually did.
2. Every extension parameter is a literature-typical guess for *S. cerevisiae* (marked in the `Strain` dataclass). These are what the NRA strain-design step will perturb.
3. Measurement noise, analytical artefacts and faults are *not* in this module on purpose. They belong in the measurement and rendering layer, so the truth stays clean for answer keys.
