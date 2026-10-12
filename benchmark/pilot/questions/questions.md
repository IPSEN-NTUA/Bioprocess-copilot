# Pilot benchmark questions

Gold answers are computed from the corpus files by `generator/questions.py` and `generator/questions_extra.py`. Difficulty tiers: L1 lookup, L2 single document + calculation, L3 joins, L4 provenance / faults / harder calculation, L5 synthesis, U unanswerable or false premise.

Interoperability tags: `id_resolution` - the same entity is named differently across sources (FB02 / FB-002 / FB002-S01, YX-002-c1 / 'YX-002 (clone c3)'); `lod_encoding` - below-detection results are written as '<0.05', '<LOD', 'n.d.' or a blank cell; `number_date_format` - comma decimals and dd/mm/yy or dd.mm.yyyy dates in some sources; `time_alignment` - process time must be reconstructed from wall-clock stamps, logger start or analyser start; `vocabulary` - different names for the same quantity (Xol / Xylitol / 'Xylitol (g/L)', formic acid / formate); `units_conversion` - pump counts, balance grams and mL must be converted with documented factors; `versioning` - superseded results, method versions and corrected values; `source_authority` - sources disagree; the answer must pick the authoritative one and say why; `implicit_context` - a fact is not written in the record and must be inferred from another document.


## Core set (difficulty tiers L1-L5, U)


### L1-01 (L1, lookup, 1 hop)

**Q:** Which strain was used in bioreactor run FB-013?

**A:** YX-005


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/notes/FB-013_MD.md  
_Grading:_ exact  
_Interop:_ id_resolution


### L1-02 (L1, lookup, 1 hop)

**Q:** In which bioreactor vessel was FB-022 run?

**A:** R1


_Evidence:_ corpus/notes/FB-022_PN.md  
_Grading:_ exact


### L1-03 (L1, lookup, 1 hop)

**Q:** What formic acid concentration does the certificate of analysis report for C5 lot C5-2503?

**A:** 23.74 g/L


_Evidence:_ corpus/materials/CoA_C5-2503.md  
_Grading:_ numeric


### L1-04 (L1, lookup, 1 hop)

**Q:** Who constructed strain YX-004?

**A:** Mark Doherty


_Evidence:_ corpus/strains/strain_registry.xlsx, corpus/strains/construction_YX-004.md  
_Grading:_ exact


### L1-05 (L1, lookup, 1 hop)

**Q:** From which date is version 2 of the HPLC method SOP-ANA-003 effective?

**A:** 2025-06-01


_Evidence:_ corpus/sops/SOP-ANA-003_v2.md  
_Grading:_ exact


### L1-06 (L1, lookup, 1 hop)

**Q:** According to the fed-batch SOP, when is the feed started and at what rate?

**A:** At 7.5 h, at 2.0 mL/h


_Evidence:_ corpus/sops/SOP-FB-001_v3.md  
_Grading:_ key_points (key points: 7.5 h, 2.0 mL/h)


### L1-07 (L1, lookup, 1 hop)

**Q:** Which clone of YX-006 is recorded as selected in the strain registry?

**A:** c3


_Evidence:_ corpus/strains/strain_registry.xlsx  
_Grading:_ exact


### L1-08 (L1, lookup, 1 hop)

**Q:** What final xylitol concentration did Elena Varga's workbook report for FB-001?

**A:** 148.6 g/L


_Evidence:_ corpus/spreadsheets/EV_FB-001.xlsx  
_Grading:_ numeric


### L1-09 (L1, lookup, 1 hop)

**Q:** What dilution factor was used for the HPLC analysis of sample FB018-S06?

**A:** 1:10


_Evidence:_ corpus/hplc/processed/HPLC_results_2025-06.xlsx, corpus/hplc/raw/  
_Grading:_ exact


### L2-01 (L2, aggregation, 1 hop)

**Q:** In FB-003, at which sampling time was glucose first below the HPLC detection limit?

**A:** 8 h


_Evidence:_ corpus/spreadsheets/EV_FB-003.xlsx, corpus/hplc/processed/HPLC_results_2025-02.xlsx  
_Grading:_ numeric  
_Interop:_ lod_encoding


### L2-02 (L2, aggregation, 1 hop)

**Q:** What was the peak ethanol concentration in FB-001, and at which sample time?

**A:** 5.78 g/L at 4 h


_Evidence:_ corpus/spreadsheets/EV_FB-001.xlsx  
_Grading:_ key_points (key points: 5.78 g/L, 4 h)


### L2-03 (L2, lookup, 1 hop)

**Q:** How much base was used in FB-008 according to the lab note?

**A:** 15.2 mL


_Evidence:_ corpus/notes/FB-008_EV.md  
_Grading:_ numeric


### L2-04 (L2, aggregation, 1 hop)

**Q:** What was the maximum dry weight reached in FB-018?

**A:** 14.03 g/L


_Evidence:_ corpus/spreadsheets/EV_FB-018.xlsx  
_Grading:_ numeric


### L2-05 (L2, aggregation, 1 hop)

**Q:** In the March shake-flask clone screen, what mean xylitol did YX-002 clone c3 reach at 48 h, and from how many flasks?

**A:** 51.2 g/L from 2 flasks


_Evidence:_ corpus/notes/SF_screen_P2_PN.md, corpus/spreadsheets/PN_clone_screen_Mar2025.xlsx  
_Grading:_ key_points (key points: 51.2 g/L, 2 flasks)


### L2-06 (L2, aggregation, 1 hop)

**Q:** What was the lowest dissolved-oxygen reading in FB-005 between 2 h and 40 h after inoculation, and when did it occur?

**A:** 8.7 % at about 15.0 h


_Evidence:_ corpus/online/FB-005_R2.txt  
_Grading:_ key_points (key points: 8.7 %, ~15 h)  
_Interop:_ time_alignment


### L2-07 (L2, aggregation, 1 hop)

**Q:** How many HPLC sequences were run in May 2025?

**A:** 6


_Evidence:_ corpus/hplc/raw/  
_Grading:_ numeric


### L2-08 (L2, lookup, 1 hop)

**Q:** Which analytes were added to the HPLC method in version 2 of SOP-ANA-003?

**A:** Succinate, furfural and HMF


_Evidence:_ corpus/sops/SOP-ANA-003_v2.md  
_Grading:_ key_points (key points: succinate, furfural, HMF)


### L2-09 (L2, aggregation, 1 hop)

**Q:** How many shake flasks were in the March clone screen?

**A:** 20


_Evidence:_ corpus/spreadsheets/PN_clone_screen_Mar2025.xlsx  
_Grading:_ numeric


### L3-01 (L3, join, 3 hops)

**Q:** Which bioreactor runs used strain YX-006, and which clone of YX-006 was it?

**A:** FB-014, FB-015, FB-017, FB-018, FB-019, FB-020, FB-021, FB-022, FB-023, FB-024, FB-025 (11 runs); clone c3


_Evidence:_ corpus/notes/, corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/spreadsheets/, corpus/strains/strain_registry.xlsx  
_Grading:_ key_points (key points: FB-014, FB-015, FB-017, FB-018, FB-019, FB-020, FB-021, FB-022, FB-023, FB-024, FB-025, c3)  
_Interop:_ id_resolution, implicit_context


### L3-02 (L3, join, 2 hops)

**Q:** Rank the six C5-fraction fed-batches FB-012 to FB-017 by final xylitol (highest first).

**A:** FB-013 (67.6) > FB-014 (67.3) > FB-017 (65.9) > FB-015 (64.3) > FB-016 (62.2) > FB-012 (19.6)


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ ranking (key points: FB-013, FB-014, FB-017, FB-015, FB-016, FB-012)  
_Interop:_ id_resolution


### L3-03 (L3, join, 3 hops)

**Q:** What was the mean and standard deviation of final xylitol in the three baseline runs FB-001 to FB-003?

**A:** 145.9 +/- 2.7 g/L (FB-001 148.6, FB-002 145.8, FB-003 143.2)


_Evidence:_ corpus/spreadsheets/EV_FB-001.xlsx, corpus/spreadsheets/EV_FB-003.xlsx, corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/meetings/minutes_2025-02-14.md  
_Grading:_ numeric  
_Interop:_ id_resolution


### L3-04 (L3, join, 2 hops)

**Q:** Which person operated the most bioreactor runs, and how many?

**A:** Mark Doherty (MD), 13 runs (EV 10, PN 2)


_Evidence:_ corpus/notes/, corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ key_points (key points: Mark Doherty, 13)


### L3-05 (L3, join, 2 hops)

**Q:** For which bioreactor runs is off-gas (O2/CO2) data available, and what do those runs have in common?

**A:** FB-001, FB-002, FB-004, FB-005, FB-008, FB-011, FB-012, FB-013, FB-016, FB-017, FB-018, FB-019, FB-022, FB-023 (14 runs); they are exactly the runs in vessels R1 and R2 - the gas analyser only measures the R1/R2 exhaust streams, so no run in R3 or R4 has off-gas data.


_Evidence:_ corpus/offgas/, corpus/notes/, corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ key_points (key points: FB-001, FB-002, FB-004, FB-005, FB-008, FB-011, FB-012, FB-013, FB-016, FB-017, FB-018, FB-019, FB-022, FB-023, R1, R2)


### L3-06 (L3, provenance, 3 hops)

**Q:** Which strains were built on a YX-002 clone other than the one finally selected from the YX-002 screen?

**A:** YX-003 and YX-004 were both built on YX-002 c1 (the first PCR-confirmed clone, before the screen); the screen selected YX-002 c3, and the team decided to keep c1 as the background.


_Evidence:_ corpus/strains/construction_YX-003.md, corpus/strains/construction_YX-004.md, corpus/strains/strain_registry.xlsx, corpus/meetings/minutes_2025-03-14.md  
_Grading:_ key_points (key points: YX-003, YX-004, c1, c3)  
_Interop:_ id_resolution


### L3-07 (L3, join, 2 hops)

**Q:** Compare the final xylitol of FB-022 (no enzyme) and FB-023 (enzyme added). What was the difference?

**A:** FB-022 71.5 g/L, FB-023 74.7 g/L; difference +3.2 g/L (FB-023 higher)


_Evidence:_ corpus/spreadsheets/PN_FB022-023_lot_comparison.xlsx, corpus/notes/FB-023_PN.md  
_Grading:_ numeric


### L3-08 (L3, join, 2 hops)

**Q:** Which C5 lot and what C5 fraction were used in the feed of FB-016?

**A:** C5-2501, C5 70%


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/notes/FB-016_MD.md  
_Grading:_ key_points (key points: C5-2501, 70)  
_Interop:_ number_date_format


### L3-09 (L3, rationale, 2 hops)

**Q:** Which YX-005 clone was selected, and on what basis?

**A:** c3: highest xylitol at 72 h in the 50 % C5 flask screen (c1 69.9, c2 66.5, c3 71.6 g/L; one flask per clone).


_Evidence:_ corpus/strains/strain_registry.xlsx, corpus/notes/SF_screen_P4_PN.md, corpus/spreadsheets/PN_C5_tolerance_Apr2025.xlsx  
_Grading:_ key_points (key points: c3, 72 h, 50 %, 71.6)


### L3-10 (L3, join, 2 hops)

**Q:** Which bioreactor runs did Priya Nair operate?

**A:** FB-022, FB-023


_Evidence:_ corpus/notes/  
_Grading:_ key_points (key points: FB-022, FB-023)


### L4-01 (L4, fault_diagnosis, 3 hops)

**Q:** Why was FB-005 excluded from the strain comparison, and what evidence supports it?

**A:** The R2 feed pump over-delivered by about 15 % (worn tubing): the feed balance showed 314 mL consumed vs 273 mL by the pump counter (DEV-2025-01); on 2025-03-26 the pump was found to deliver 2.30 mL/h at 2.00 mL/h setpoint and was recalibrated. FB-005 therefore received more glucose and xylose than planned.


_Evidence:_ corpus/deviations/DEV-2025-01.md, corpus/notes/FB-005_EV.md, corpus/equipment/maintenance_log.csv, corpus/meetings/minutes_2025-04-11.md  
_Grading:_ key_points (key points: feed pump, 15 %, 314, 273, DEV-2025-01)


### L4-02 (L4, fault_diagnosis, 2 hops)

**Q:** What happened in FB-010 around 50 h, and how long is the gap in its online data?

**A:** The lab UPS tripped (power failure): the R4 controller and feed pump were off and the feed was restarted by hand (DEV-2025-03, ~3.5 h without feed). The online log has no records between 49.3 h and 53.0 h (~3.7 h gap).


_Evidence:_ corpus/online/FB010_R4.csv, corpus/notes/FB-010_MD.md, corpus/deviations/DEV-2025-03.md, corpus/equipment/maintenance_log.csv  
_Grading:_ key_points (key points: power, 3.7 h, feed)  
_Interop:_ time_alignment


### L4-03 (L4, fault_diagnosis, 3 hops)

**Q:** Which P3 run has an unreliable dissolved-oxygen signal in the second half, and what is the evidence?

**A:** FB-009 (R3, YX-003): pO2 falls steadily from about 73 % at 40 h to 38 % at 140 h, whereas its replicate FB-006 (same strain, same vessel) stays at about 81-67 % and the offline data of FB-009 are normal. The R3 pO2 probe showed slow response and low readings and its membrane was replaced on 2025-04-07 (maintenance log; minutes 2025-04-11).


_Evidence:_ corpus/online/FB009_R3.csv, corpus/online/FB006_R3.csv, corpus/equipment/maintenance_log.csv, corpus/meetings/minutes_2025-04-11.md  
_Grading:_ key_points (key points: FB-009, drift, membrane, 2025-04-07)  
_Interop:_ time_alignment


### L4-04 (L4, provenance, 4 hops)

**Q:** Which C5 lot was FB-024 planned with, which lot was actually used, and how can this be shown from the data?

**A:** Planned C5-2503 (run sheet and MD note); actually C5-2502 (MD 'Media & feeds' sheet, DEV-2025-04). The FB-024 feed sample contained 7.7 g/L formate - like the C5-2502 feed of FB-022 (7.7) and the 60 % of the C5-2502 CoA value (7.9 g/L), and unlike the C5-2503 feed of FB-025 (14.8; CoA-based 14.2 g/L).


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/hplc/processed/HPLC_results_2025-06.xlsx, corpus/materials/CoA_C5-2502.md, corpus/materials/CoA_C5-2503.md, corpus/deviations/DEV-2025-04.md, corpus/reports/weekly_PN_2025-07-04.md  
_Grading:_ key_points (key points: C5-2503, C5-2502, 7.7, formate)  
_Interop:_ vocabulary, source_authority


### L4-05 (L4, provenance, 3 hops)

**Q:** For sample FB005-S07, what xylitol value was first reported, what is the corrected value, and why did it change?

**A:** First 103.88 g/L (SEQ-20250321), corrected 97.82 g/L from the re-run SEQ-20250409: the xylitol standard used in SEQ-20250321 was prepared from a partially hydrated batch, so xylitol results of that sequence were biased high (DEV-2025-02).


_Evidence:_ corpus/hplc/processed/HPLC_results_2025-03.xlsx, corpus/hplc/raw/, corpus/deviations/DEV-2025-02.md  
_Grading:_ key_points (key points: 103.88, 97.82, DEV-2025-02)  
_Interop:_ versioning


### L4-06 (L4, data_quality, 3 hops)

**Q:** Is any xylitol value in Mark Doherty's master log inconsistent with the HPLC results? If so, which and why?

**A:** Yes. (1) FB-015 S08 is logged as 490.7 g/L, but the HPLC result is 49.07 g/L - a transcription error (decimal point shifted, factor 10) that also breaks the time course of the run. (2) 8 values (FB-006 S04, S05, S06, S07; FB-007 S04, S05, S06, S07) still show the original results of HPLC sequence SEQ-20250321; they were superseded by the re-run SEQ-20250409 after DEV-2025-02 (biased xylitol standard), but the master log was not updated.


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/hplc/processed/, corpus/deviations/DEV-2025-02.md  
_Grading:_ key_points (key points: FB-015 S08, 490.7, 49.07, SEQ-20250321, not updated, FB-006, FB-007)  
_Interop:_ versioning, source_authority, number_date_format


### L4-07 (L4, fault_diagnosis, 2 hops)

**Q:** Did the antifoam added to FB-014 affect dissolved oxygen? Quantify.

**A:** Yes: antifoam was added at ~30 h (MD note); pO2 fell from about 79 % (27-30 h) to about 53 % (31-34 h), consistent with antifoam lowering the oxygen transfer (kLa).


_Evidence:_ corpus/notes/FB-014_MD.md, corpus/online/FB014_R3.csv  
_Grading:_ key_points (key points: 30 h, 79, 53)  
_Interop:_ time_alignment


### L4-08 (L4, rationale, 3 hops)

**Q:** Why did FB-002 start more slowly than FB-001 and FB-003?

**A:** The preculture was two days old (MD note), giving a lower start OD600 (1.36 vs 2.04 and 2.14) and lower dry weight at 4 h (3.77 vs 5.80 and 5.62 g/L); the final titre was nevertheless close to the others.


_Evidence:_ corpus/notes/FB-002_MD.md, corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/spreadsheets/EV_FB-001.xlsx, corpus/spreadsheets/EV_FB-003.xlsx  
_Grading:_ key_points (key points: preculture, 1.36, 3.77)  
_Interop:_ number_date_format


### L4-09 (L4, provenance, 4 hops)

**Q:** Trace the lineage of YX-006 back to the parent strain, including the clones used and the genetic change at each step.

**A:** YX-006 = YX-003-c3 + heterologous FDH (plasmid); YX-003 = YX-002-c1 + pTEF1-ZWF1 and pTEF1-GND1 (oxidative PPP for NADPH); YX-002 = YX-001 + second, integrated XYL1 copy (pTDH3); YX-001 = the lab host with plasmid-borne xylose reductase (XYL1).


_Evidence:_ corpus/strains/construction_YX-006.md, corpus/strains/construction_YX-003.md, corpus/strains/construction_YX-002.md, corpus/strains/construction_YX-001.md, corpus/strains/strain_registry.xlsx  
_Grading:_ key_points (key points: YX-003-c3, YX-002-c1, FDH, ZWF1, GND1, XYL1)


### L5-01 (L5, synthesis, 5 hops)

**Q:** In the P3 bioreactor comparison (FB-004 to FB-011, excluding FB-005, using corrected HPLC values), which engineered strain gave the highest mean final xylitol, by how much compared with the parent control, and how reliable is that difference?

**A:** YX-004: 147.4 g/L vs control 139.6 g/L (+5.6 %); YX-003 143.0, YX-002 139.7 (FB-008 only, because FB-005 is excluded). The evidence is weak: only 1-2 runs per strain, a difference of 7.9 g/L against a run-to-run sd of about 2.7 g/L (baseline FB-001 to FB-003), and the control mean includes FB-010 (136.4 g/L), which had ~3.5 h without feed in the power cut (DEV-2025-03). Against FB-004 alone (142.8 g/L) the gain shrinks to +3.3 %. Round-1 values also had to be corrected for the HPLC standard error (DEV-2025-02).


_Evidence:_ corpus/hplc/processed/, corpus/lims/sample_register.csv, corpus/deviations/DEV-2025-01.md, corpus/deviations/DEV-2025-02.md, corpus/deviations/DEV-2025-03.md, corpus/meetings/minutes_2025-02-14.md, corpus/meetings/minutes_2025-04-11.md  
_Grading:_ key_points (key points: YX-004, 147.4, 139.6, FB-010, few runs)  
_Interop:_ versioning


### L5-02 (L5, synthesis, 5 hops)

**Q:** Is the heterologous FDH needed for xylitol production at 70 % C5 fraction? Support the answer with flask and bioreactor data.

**A:** Yes. In flasks at 70 % C5 the parent YX-001 made 1.5 g/L xylitol and did not grow (OD 2.8 at 72 h from OD 2.4), while the FDH strains grew and made xylitol. In the bioreactor, FB-012 (YX-003, no FDH, 70 % C5) stopped working by ~96 h (xylose accumulating, formate 7.2 g/L at the end) and finished at only 19.6 g/L xylitol, DW 3.0 g/L; with FDH at 70 % C5 the runs reached FB-015 64.3, FB-016 62.2, FB-017 65.9 g/L - FB-015 and FB-017 with YX-006, which is YX-003 plus FDH, i.e. the same background as FB-012.


_Evidence:_ corpus/spreadsheets/PN_C5_tolerance_Apr2025.xlsx, corpus/meetings/minutes_2025-05-02.md, corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/notes/FB-012_MD.md, corpus/strains/construction_YX-006.md, corpus/meetings/minutes_2025-05-30.md  
_Grading:_ key_points (key points: yes, formate, FB-012, 19.6, FB-015)


### L5-03 (L5, synthesis, 5 hops)

**Q:** In the P6 study with YX-006 (900 vs 500 rpm, feed A vs the more glucose-rich feed B), which condition was oxygen-limited, what is the evidence, and did it measurably reduce the final xylitol titre?

**A:** FB-021 (500 rpm, feed B): its fed-phase pO2 median was 1 % (vs 74 / 41 / 52 % for FB-018/019/020), and ethanol re-appeared from 96 h (3.2 g/L at the end) while it stayed below detection in the other three runs - overflow/fermentative metabolism under O2 limitation. Final xylitol was 143.6 / 150.0 / 151.4 / 146.1 g/L; the 5.3 g/L lower titre of FB-021 vs FB-020 is a single-run comparison within the run-to-run scatter (baseline sd 2.7 g/L), so a titre loss is not demonstrated.


_Evidence:_ corpus/spreadsheets/EV_FB-018.xlsx, corpus/spreadsheets/EV_FB-019.xlsx, corpus/spreadsheets/EV_FB-020.xlsx, corpus/spreadsheets/EV_FB-021.xlsx, corpus/notes/FB-021_EV.md, corpus/online/, corpus/meetings/minutes_2025-06-20.md, corpus/meetings/minutes_2025-02-14.md  
_Grading:_ key_points (key points: FB-021, 1 %, ethanol, not significant)


### L5-04 (L5, synthesis, 4 hops)

**Q:** Taking the lot mix-up into account, what can be concluded about the effect of the enzyme on lot C5-2503?

**A:** Nothing: the only real C5-2503 run is FB-025 (with enzyme, 69.87 g/L); FB-024, meant as the no-enzyme C5-2503 control, was fed with C5-2502 (71.14 g/L), so there is no C5-2503 run without enzyme to compare against (a repeat was planned for August). For C5-2502 the enzyme changed final xylitol by +3.2 g/L (FB-023 vs FB-022, single runs).


_Evidence:_ corpus/deviations/DEV-2025-04.md, corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/spreadsheets/PN_FB022-023_lot_comparison.xlsx, corpus/meetings/minutes_2025-07-11.md  
_Grading:_ key_points (key points: cannot, FB-025, FB-024, C5-2502)


### L5-05 (L5, data_quality, 3 hops)

**Q:** Do any HPLC results of FB-012 or FB-013 break the time course in a way that suggests a sample mix-up? Explain.

**A:** Yes, sample S07 (72 h) of both runs. FB-012 S07 shows xylitol 39.87 g/L (S06 23.29, S08 25.3) and no formate although formate is present in all its other fed-phase samples; FB-013 S07 shows xylitol 26.96 g/L (S06 27.6, S08 46.69) and formate 4.32 g/L although all other FB-013 samples have at most 0.08 g/L formate. Each S07 looks like the other run's sample, so the two tubes were most likely swapped before HPLC (both sampled the same morning, 2025-05-08 09:30 and 2025-05-08 10:05). OD and dry weight of S07 fit their own runs.


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/hplc/processed/HPLC_results_2025-05.xlsx, corpus/lims/sample_register.csv  
_Grading:_ key_points (key points: S07, swap, FB-012, FB-013, formate)  
_Interop:_ id_resolution


### L5-06 (L5, synthesis, 4 hops)

**Q:** Which bioreactor run reached the highest final xylitol in the whole campaign, and is that result trustworthy?

**A:** FB-005 (153.7 g/L), but it is not trustworthy: its feed pump over-delivered by ~15 % (DEV-2025-01), so it received more substrate and was excluded. The highest valid run is FB-020 (151.4 g/L), followed by FB-019 (150.0), FB-001 (148.6) - differences between these are within run-to-run scatter.


_Evidence:_ corpus/hplc/processed/, corpus/lims/sample_register.csv, corpus/deviations/, corpus/spreadsheets/  
_Grading:_ key_points (key points: FB-005, 153.7, DEV-2025-01, FB-020)  
_Interop:_ versioning, source_authority


### U-01 (U, unanswerable, 1 hop)

**Q:** What was the dissolved CO2 concentration in the broth of FB-019?

**A:** Not available: dissolved CO2 was not measured. Only the exhaust-gas CO2 (%) was recorded for FB-019 (off-gas file).


_Evidence:_ corpus/offgas/GA4_fb019_R2.xlsx, corpus/sops/  
_Grading:_ unanswerable (key points: not measured)


### U-02 (U, unanswerable, 1 hop)

**Q:** What succinate concentration was measured in FB-006?

**A:** Not available: FB-006 (March 2025) samples were analysed with SOP-ANA-003 v1, which does not include succinate; succinate was only added in v2, effective 2025-06-01.


_Evidence:_ corpus/sops/SOP-ANA-003_v1.md, corpus/sops/SOP-ANA-003_v2.md, corpus/hplc/processed/HPLC_results_2025-03.xlsx  
_Grading:_ unanswerable (key points: not measured, v1, 2025-06-01)  
_Interop:_ versioning


### U-03 (U, unanswerable, 1 hop)

**Q:** What was the respiratory quotient (RQ) of FB-006 during the fed phase?

**A:** Not available: there is no off-gas data for FB-006 (vessel R3); off-gas was only recorded for R1 and R2 runs.


_Evidence:_ corpus/offgas/  
_Grading:_ unanswerable (key points: no off-gas, R3)


### U-04 (U, unanswerable, 1 hop)

**Q:** What final xylitol titre did run FB-026 reach?

**A:** There is no run FB-026 in the records (runs go from FB-001 to FB-025); the minutes of 2025-07-11 only plan a repeat for August.


_Evidence:_ corpus/lims/sample_register.csv, corpus/meetings/minutes_2025-07-11.md  
_Grading:_ unanswerable (key points: no such run)


### U-05 (U, false_premise, 1 hop)

**Q:** Why was strain YX-004 dropped after the P3 comparison?

**A:** False premise: YX-004 was not dropped. The registry lists it as 'on hold', and the 2025-04-11 minutes decided to keep it on hold.


_Evidence:_ corpus/strains/strain_registry.xlsx, corpus/meetings/minutes_2025-04-11.md  
_Grading:_ unanswerable (key points: not dropped, on hold)


### U-06 (U, unanswerable, 1 hop)

**Q:** How much did lot C5-2502 of the C5 fraction cost?

**A:** Not available: no document gives prices or costs for the C5 lots.


_Evidence:_ corpus/materials/  
_Grading:_ unanswerable (key points: not documented)


### U-07 (U, false_premise, 1 hop)

**Q:** What final xylitol did the C5-2503 run without enzyme reach?

**A:** There is no such run: FB-024 was planned as the no-enzyme C5-2503 run but was fed with C5-2502 (DEV-2025-04). The only C5-2503 run, FB-025, had enzyme.


_Evidence:_ corpus/deviations/DEV-2025-04.md, corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ unanswerable (key points: no such run, FB-024, C5-2502)


## Semantic interoperability (IO)


### IO-01 (L2, lookup, 2 hops)

**Q:** What xylitol concentration was measured in FB-006 at 96 h?

**A:** 112.63 g/L (sample FB006-S08; written as run 'FB06', sample 'S08', column 'Xol' in the MD master log and as 'FB006-S08' / 'Xylitol' in the HPLC results).


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/hplc/processed/HPLC_results_2025-03.xlsx, corpus/lims/sample_register.csv  
_Grading:_ numeric  
_Interop:_ id_resolution, vocabulary


### IO-02 (L2, lookup, 2 hops)

**Q:** What glycerol concentration was measured in FB-002 at 23.5 h?

**A:** Below the detection limit: the MD master log shows '<0.05' and the HPLC results flag it '<LOD' (LOD 0.02-0.05 g/L in the vial). It was measured, so it is not missing, and it should not be reported as 0 or as a number.


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/hplc/processed/HPLC_results_2025-02.xlsx, corpus/sops/SOP-ANA-003_v1.md  
_Grading:_ key_points (key points: below detection, <0.05)  
_Interop:_ lod_encoding, id_resolution  
_Trap:_ 0 g/L, or 'not measured'


### IO-03 (L2, lookup, 2 hops)

**Q:** On which date and at what time was FB-002 inoculated?

**A:** 03 February 2025, 10:15 (2025-02-03). MD writes dates as dd/mm/yy ('03/02/25' in his note and master log), which must not be read as month/day; the LIMS register confirms the date.


_Evidence:_ corpus/notes/FB-002_MD.md, corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/lims/sample_register.csv  
_Grading:_ key_points (key points: 2025-02-03, 10:15)  
_Interop:_ number_date_format  
_Trap:_ 02 March 2025 (month/day read the wrong way round)


### IO-04 (L3, time_series, 3 hops)

**Q:** What was the dissolved-oxygen (pO2) reading in FB-021 at 47.5 h after inoculation?

**A:** About 32.6 % (record 2025-06-11T10:30:00). In the VesselView export Elapsed_h counts from the logger start, which was 15.8 h before inoculation, so the row with Elapsed_h = 47.5 is actually 31.7 h after inoculation (pO2 2.4 %). Inoculation time from the LIMS S01 sample / run sheet.


_Evidence:_ corpus/online/FB021_R4.csv, corpus/lims/sample_register.csv, corpus/spreadsheets/EV_FB-021.xlsx, corpus/equipment/reactor_setup.md  
_Grading:_ numeric  
_Interop:_ time_alignment  
_Trap:_ 2.4 % (row with Elapsed_h = 47.5)


### IO-05 (L3, time_series, 3 hops)

**Q:** How many hours after inoculation did the exhaust CO2 of FB-001 peak, and what was the peak value?

**A:** At about 3.5 h after inoculation (03.02.2025 13:00), CO2 1.17 %, during the glucose batch phase. Inoculation is in the BioCtrl header (03 Feb 2025 09:30); the analyser RelTime counts from the analyser start (2025-02-02 17:17), so RelTime/3600 = 19.7 h is not process time.


_Evidence:_ corpus/offgas/GA4_fb001_R1.xlsx, corpus/online/FB-001_R1.txt, corpus/equipment/reactor_setup.md  
_Grading:_ key_points (key points: 3.5 h, 1.17 %)  
_Interop:_ time_alignment  
_Trap:_ 19.7 h (analyser RelTime)


### IO-06 (L3, lookup, 2 hops)

**Q:** What was the arabinose concentration at the end of FB-016?

**A:** 6.39 g/L (FB016-S11, HPLC results). The MD master log has no arabinose column, so the value is only in the HPLC results workbook; arabinose comes in with the C5 fraction and is not consumed.


_Evidence:_ corpus/hplc/processed/HPLC_results_2025-05.xlsx, corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ numeric  
_Interop:_ source_authority, id_resolution  
_Trap:_ 'not measured' (because the MD log lacks the column)


### IO-07 (L3, provenance, 3 hops)

**Q:** What is the valid xylitol concentration of sample FB006-S05?

**A:** 37.54 g/L from the re-run sequence SEQ-20250409. The first result (41.63 g/L, SEQ-20250321) was superseded after DEV-2025-02 (biased xylitol standard); the MD master log still shows the old value (41.63).


_Evidence:_ corpus/hplc/processed/HPLC_results_2025-03.xlsx, corpus/deviations/DEV-2025-02.md, corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ numeric  
_Interop:_ versioning, source_authority  
_Trap:_ 41.63 g/L (superseded)


### IO-08 (L3, provenance, 3 hops)

**Q:** Was YX-003 built on the same YX-002 clone that was tested in the bioreactor runs FB-005 and FB-008?

**A:** No. YX-003 was built on YX-002-c1 (construction record, registry parent), while FB-005 and FB-008 used YX-002 clone c3 (EV workbooks: 'YX-002 (clone c3)'). The same clone is written 'YX-002-c1', 'YX-002 c1' and 'clone c1' in different documents.


_Evidence:_ corpus/strains/construction_YX-003.md, corpus/strains/strain_registry.xlsx, corpus/spreadsheets/EV_FB-005.xlsx, corpus/spreadsheets/EV_FB-008.xlsx  
_Grading:_ key_points (key points: no, c1, c3)  
_Interop:_ id_resolution


### IO-09 (L2, lookup, 2 hops)

**Q:** What final dry weight did Mark Doherty write in his FB-016 note, and does it agree with his master log?

**A:** 3.3 g/L - written '3,3' with a decimal comma in the note; the master log gives 3.28 g/L, so they agree.


_Evidence:_ corpus/notes/FB-016_MD.md, corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ key_points (key points: 3.3, agree)  
_Interop:_ number_date_format  
_Trap:_ 33 g/L (comma read as a thousands separator)


### IO-10 (L3, calculation, 3 hops)

**Q:** How much base (2 M KOH) was added in FB-002 and in FB-006?

**A:** FB-002 (R2, BioCtrl): 300 base-pump counts x 0.05 mL/count = 15.0 mL. FB-006 (R3, VesselView): 15.45 mL, logged directly in mL (BASE_mL). The two controllers log base differently; the pump calibration is in the reactor set-up document.


_Evidence:_ corpus/online/FB-002_R2.txt, corpus/online/FB006_R3.csv, corpus/equipment/reactor_setup.md  
_Grading:_ key_points (key points: 15.0, 15.45)  
_Interop:_ units_conversion, vocabulary  
_Trap:_ 300 mL (counts taken as mL)


### IO-11 (L4, provenance, 3 hops)

**Q:** The MD master log has a row 'FB15 / S08'. Which LIMS sample is this, when was it taken, and in which HPLC sequence was it analysed?

**A:** LIMS sample FB015-S08 (run FB-015, 96 h), sampled 2025-05-09 11:15, analysed on 2025-05-13 in sequence SEQ-20250513 (method v1). This is also the row with the decimal-point slip in the master log.


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/lims/sample_register.csv, corpus/hplc/processed/HPLC_results_2025-05.xlsx, corpus/hplc/raw/  
_Grading:_ key_points (key points: FB015-S08, 2025-05-09, SEQ-20250513)  
_Interop:_ id_resolution


### IO-12 (L3, provenance, 2 hops)

**Q:** PN's April flask note lists 'YX-006 c3, 70 % C5' with 64.9 g/L xylitol. Which flask is that, and what OD600 did it reach at 72 h?

**A:** Flask SF-033 (sheet YX-006 of PN_C5_tolerance_Apr2025.xlsx, 'Tolerance, 70 % C5'); OD600 18.46 at 72 h. The minutes of 2025-05-02 refer to the same flask by its ID.


_Evidence:_ corpus/notes/SF_screen_P4_PN.md, corpus/spreadsheets/PN_C5_tolerance_Apr2025.xlsx, corpus/meetings/minutes_2025-05-02.md  
_Grading:_ key_points (key points: SF-033, 18.46)  
_Interop:_ id_resolution


### IO-13 (L4, inference, 3 hops)

**Q:** Which clone of YX-005 was used in FB-013?

**A:** Clone c3. MD's records only say 'YX-005'; the clone follows from the strain registry / construction record (clone c3 selected after the flask screen on 2025-04-29), which was before FB-013 started on 2025-05-05.


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/notes/FB-013_MD.md, corpus/strains/strain_registry.xlsx, corpus/strains/construction_YX-005.md  
_Grading:_ key_points (key points: c3)  
_Interop:_ implicit_context, number_date_format


## Calculations (CA)


### CA-01 (L2, calculation, 1 hop)

**Q:** Estimate the specific growth rate of FB-001 in the batch phase from the dry weight between the first two samples.

**A:** mu = ln(5.80/0.82) / (4 - 0.1 h) = 0.50 1/h (OD600 gives a similar value).


_Evidence:_ corpus/spreadsheets/EV_FB-001.xlsx  
_Grading:_ numeric


### CA-02 (L2, calculation, 2 hops)

**Q:** Recompute the xylitol concentration of sample FB018-S06 from the raw HPLC sequence file.

**A:** 6.1226 g/L in the vial x dilution 10 = 61.23 g/L (SEQ-20250613.csv); this matches the processed result.


_Evidence:_ corpus/hplc/raw/SEQ-20250613.csv, corpus/hplc/processed/HPLC_results_2025-06.xlsx  
_Grading:_ numeric  
_Interop:_ units_conversion


### CA-03 (L3, calculation, 3 hops)

**Q:** How much feed (mL) was added to FB-006?

**A:** Balance loss 805.42 - 484.60 = 320.8 g; divided by the Feed A density 1.15 g/mL (SOP-FB-001) = 279 mL. VesselView has no feed pump counter, so the balance is the only record.


_Evidence:_ corpus/online/FB006_R3.csv, corpus/sops/SOP-FB-001_v3.md, corpus/equipment/reactor_setup.md, corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ numeric  
_Interop:_ units_conversion


### CA-04 (L3, calculation, 4 hops)

**Q:** Check the claim in DEV-2025-01 that the R2 feed pump over-delivered by about 15 % in FB-005, using the raw records.

**A:** Pump: 2728 counts x 0.10 mL = 272.8 mL. Balance: (769.3 - 408.5) g / 1.15 g/mL = 313.7 mL. Ratio 1.150, i.e. 15 % more delivered than counted - the claim holds.


_Evidence:_ corpus/online/FB-005_R2.txt, corpus/spreadsheets/EV_FB-005.xlsx, corpus/sops/SOP-FB-001_v3.md, corpus/equipment/reactor_setup.md, corpus/deviations/DEV-2025-01.md  
_Grading:_ numeric  
_Interop:_ units_conversion


### CA-05 (L3, calculation, 3 hops)

**Q:** At most, how much extra xylose (g per litre of feed) could the enzyme release from the FB-023 feed, if it hydrolysed all xylo-oligosaccharides? Compare it with the xylose in that feed.

**A:** The feed is 60 % of lot C5-2502, whose CoA lists 42.82 g/L xylo-oligosaccharides (xylose equivalents): 0.6 x 42.82 = 25.7 g/L of feed, about 18 % on top of the 140 g/L xylose the feed is topped up to (SOP). Oligomers are not measured by the HPLC method, so the actual release was not verified.


_Evidence:_ corpus/materials/CoA_C5-2502.md, corpus/sops/SOP-FB-001_v3.md, corpus/notes/FB-023_PN.md, corpus/sops/SOP-ANA-003_v2.md  
_Grading:_ numeric  
_Interop:_ vocabulary


### CA-06 (L4, calculation, 3 hops)

**Q:** Calculate the respiratory quotient (RQ) of FB-018 during the fed phase (after 30 h) from the off-gas data.

**A:** With an inert-gas balance (inlet 20.95 % O2 / 0.04 % CO2; outlet dry gas): RQ = (yCO2,out*f - yCO2,in) / (yO2,in - yO2,out*f), f = (1 - yO2,in - yCO2,in)/(1 - yO2,out - yCO2,out). Median over the fed phase = 1.80 (interquartile range 1.60-2.10). The O2 difference over the reactor is only ~0.1 %, so single cycles are noisy. RQ > 1 is expected: xylose reductase takes NAD(P)H that would otherwise be reoxidised with O2.


_Evidence:_ corpus/offgas/GA4_fb018_R1.xlsx, corpus/online/FB-018_R1.txt, corpus/equipment/reactor_setup.md  
_Grading:_ numeric  
_Interop:_ time_alignment


### CA-07 (L4, calculation, 2 hops)

**Q:** What concentration of undissociated formic acid was present at the end of FB-012?

**A:** Formate at the end was 7.23 g/L (MD master log); at the controlled pH 5.5 (SOP) and pKa 3.75 the undissociated fraction is 1/(1 + 10^(5.5 - 3.75)) = 0.0175, so about 0.13 g/L. The pKa of formic acid is general chemistry, not in the corpus.


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/sops/SOP-FB-001_v3.md  
_Grading:_ numeric  
_Interop:_ vocabulary


### CA-08 (L5, calculation, 4 hops)

**Q:** Estimate the xylitol yield on consumed xylose (g/g) for FB-018 with a mass balance, and say whether it is physically plausible.

**A:** Xylose fed: 275 mL feed (balance loss / 1.15 g/mL) x 307.9 g/L (feed sample FB018-F) = 84.8 g. Xylose left: 12.55 g/L x 500 mL in the reactor plus 0.30 g removed in 7 mL samples = 6.6 g. Xylitol formed: 143.6 g/L x 500 mL + 3.23 g in samples = 75.0 g. Yield = 75.0/78.2 = 0.96 g/g. Plausible: the stoichiometric maximum is 1.01 g/g (152.15/150.13), and the estimate is within the measurement error of it (volumes are recorded to 10 mL, HPLC ~3 %).


_Evidence:_ corpus/spreadsheets/EV_FB-018.xlsx, corpus/hplc/processed/HPLC_results_2025-06.xlsx, corpus/sops/SOP-FB-001_v3.md  
_Grading:_ numeric  
_Interop:_ units_conversion


### CA-09 (L3, calculation, 3 hops)

**Q:** In the P3 comparison (excluding FB-005), which strain left the least residual xylose at the end (mean of its runs)?

**A:** YX-003 (14.0 g/L); YX-001 16.1, YX-002 14.7, YX-004 14.7 g/L. Residual xylose reflects how much of the fed xylose a strain converted, so it is a more direct readout of conversion capacity than the final titre.


_Evidence:_ corpus/hplc/processed/, corpus/lims/sample_register.csv, corpus/deviations/DEV-2025-01.md  
_Grading:_ key_points (key points: YX-003, 14.0)  
_Interop:_ id_resolution


### CA-10 (L4, statistics, 3 hops)

**Q:** Using the run-to-run scatter of the three baseline runs, how many runs per strain would be needed to detect a 5 g/L difference in final xylitol between two strains (two-sided alpha 0.05, 80 % power)?

**A:** sd = 2.7 g/L (FB-001 to FB-003, n = 3). n = 2((1.96 + 0.84) sd / 5)^2 = 4.6, i.e. about 5 runs per strain. With only 3 baseline runs the sd itself is very uncertain, so this is a lower bound; the P3 comparison used 1-2 runs per strain.


_Evidence:_ corpus/spreadsheets/EV_FB-001.xlsx, corpus/spreadsheets/EV_FB-003.xlsx, corpus/spreadsheets/MD_fermentation_master_2025.xlsx  
_Grading:_ numeric  
_Interop:_ id_resolution


## Experiment design (XD)


### XD-01 (L5, experiment_design, 6 hops)

**Q:** Based on all the data so far, what should the team try next? Propose up to five experiments, and for each give the evidence that motivates it and what result would change a decision.

**A:** Any well-argued subset of the rubric items; at least three grounded in the corpus evidence, with no penalised claim.


| Item | Suggestion | Evidence in the corpus | Simulator check (answer key) |
|---|---|---|---|
| R1 | Repeat the C5-2503 run without enzyme (ideally with replicates). | FB-024 was fed C5-2502 by mistake (DEV-2025-04), so lot C5-2503 has no no-enzyme run; the 2025-07-11 minutes plan a repeat. | Expected enzyme effect at 5 U/L on C5-2503 only +1.4 g/L, below single-run scatter: one repeat run will not settle it. |
| R2 | Re-run the strain comparison with more replicates per strain (balanced over vessels/days) and a more sensitive readout (residual xylose, xylitol at 48-72 h, or a higher feed rate). | Differences between strains (control 139.6, best 147.4 g/L) are of the order of the run-to-run sd (2.7 g/L); only 1-2 runs per strain; YX-002 has a single valid run (FB-008); the control includes the power-cut run FB-010. | True YX-003 vs YX-001 difference 6.04 g/L with run-to-run sd 4.82 g/L: about 10 runs per strain needed. |
| R3 | Increase the xylose supply (higher feed rate or more concentrated xylose) - titre is limited by the xylose fed. | Feed A carries 300 g/L xylose and gives ~140-150 g/L xylitol; the C5 feeds are topped up to only 140 g/L xylose and give ~62-75 g/L, even with FDH strains whose formate is fully removed; residual xylose stays low (~7-18 g/L), so most fed xylose is converted. | Feed A at 3 instead of 2 mL/h: YX-006 145 -> 168 g/L, and the YX-003 vs YX-001 difference grows from 5.6 to 8.9 g/L. |
| R4 | With the glucose-rich feed B keep oxygen non-limiting: >= 700-900 rpm or a pO2-controlled stirrer cascade. | FB-021 (500 rpm, feed B) had a fed-phase pO2 median of 1 % and ethanol re-appeared late; FB-020 (900 rpm) stayed at 52 %. | Feed B: 500 rpm 143 g/L with 6.3 g/L ethanol; 700 rpm 148, 900 rpm 148, pO2 cascade at 20 % 148 g/L. |
| R5 | Push the C5 fraction above 70 % (up to undiluted) with YX-006, and test lot C5-2503 (more formic acid) at 70 %. | With FDH, 70 % C5 (FB-015 64.3, FB-017 65.9 g/L) performs like 50 % (FB-014 67.3) and formate is removed (end formate below detection (<0.05 g/L) in FB-014, FB-015 and FB-017); C5-2503 has the highest formic acid in its CoA. | YX-006 on C5-2501: 70 % 68, 80 % 67, 90 % 66, 100 % 68 g/L - no penalty up to 100 %. |
| R6 | Do an enzyme dose-response and measure the oligomers directly (e.g. acid hydrolysis), since the HPLC method does not see them. | The +3.2 g/L of FB-023 vs FB-022 is a single-run difference within scatter; the CoAs list 29-43 g/L oligomers, i.e. up to ~26 g/L extra xylose in a 60 % feed. | C5-2502: no enzyme 70, 5 U/L 72, 20 U/L 76 g/L. |
| R7 | Test PN's idea that the PPP strain (YX-003) needs more glucose to show its NADPH advantage: YX-001 vs YX-003 on feed B. | PN weekly report 2025-03-14; the P6 feed-B runs used only YX-006. | Informative but negative: YX-003 advantage on feed B +2.8 g/L vs on feed A +5.6 g/L. |
| R8 | Close analytical gaps: analyse C5 runs with HPLC v2 (succinate, furfural, HMF), add a viability measurement, and repeat the single YX-002 run. | P5 C5 runs were analysed with v1 (no furans); FB-012's death is only inferred from xylose accumulation; YX-002 rests on FB-008 alone. | Not simulated (measurement gap, not a process change). |


Penalised: (1) Uses FB-005 as valid data (excluded, DEV-2025-01). (2) Uses superseded SEQ-20250321 xylitol values or the 490.7 g/L slip in the MD log. (3) States that a strain is clearly better than the parent from the P3 data without addressing the scatter. (4) States that the enzyme clearly increases xylitol, or treats FB-024 as a C5-2503 run. (5) Recommends unbuffered flasks or a strain without FDH for >= 50 % C5.


_Evidence:_ corpus/  
_Grading:_ rubric (key points: R1, R2, R3, R4, R5, R6, R7, R8)  
_Interop:_ source_authority, versioning


### XD-02 (L5, experiment_design, 5 hops)

**Q:** The team wants to show in the bioreactor that the PPP strain YX-003 makes more xylitol than the parent YX-001. Design that experiment (conditions, readout, number of runs, controls) and justify each choice from the data.

**A:** Key elements: (1) replicates - the measured difference so far is +3.4 g/L against a run-to-run sd of 2.7 g/L, so ~10 runs per strain (power calculation), not 1-2; (2) make xylose conversion, not xylose supply, limiting - higher feed rate or xylose concentration - and use residual xylose and xylitol at 48-72 h as readouts besides final titre; (3) balance strains over vessels R1-R4 and days, same feed batch, and keep pO2 non-limiting (900 rpm); (4) exclude or repeat runs with deviations (feed pump, power cut, probe drift) and use re-run-corrected HPLC values.


| Item | Suggestion | Evidence in the corpus | Simulator check (answer key) |
|---|---|---|---|
| D1 | replicates from a power calculation | sd 2.7 g/L, difference +3.4 g/L | 10 runs per strain at standard conditions |
| D2 | higher xylose supply so that conversion capacity limits | titre follows the xylose fed; residual xylose low | difference grows from 5.6 to 8.9 g/L at 3 mL/h |
| D3 | residual xylose / intermediate xylitol as readouts | P3 residual xylose: lowest for YX-003 | residual xylose at 3 mL/h: YX-001 25.6 vs YX-003 16.6 g/L |
| D4 | balanced vessel/day assignment, non-limiting O2, exclude deviation runs | vessel kLa differ, FB-005/FB-010/FB-009 deviations | Not simulated. |


Penalised: (1) Relies on final titre from 1-2 runs per strain. (2) Uses FB-005 or superseded HPLC values.


_Evidence:_ corpus/  
_Grading:_ rubric (key points: replicates, feed rate, residual xylose, balanced vessels)  
_Interop:_ source_authority


### XD-03 (L4, experiment_design, 5 hops)

**Q:** Is it worth running YX-006 on undiluted (100 %) C5 fraction? What does the existing data say, and what risk remains untested?

**A:** Probably yes: with FDH, 70 % C5 gave the same titre as 50 % (FB-015 64.3, FB-017 65.9 vs FB-014 67.3 g/L) and formate was removed (end formate below detection (<0.05 g/L) in FB-014, FB-015 and FB-017 in the master log), while the strain without FDH (FB-012) failed at 70 %. Untested: fractions above 70 %, furfural/HMF and acetate at higher fractions (P5 runs were analysed with HPLC v1, which has no furans), and lot C5-2503, which has more formic acid. Titre will stay limited by the xylose in the feed unless the xylose top-up is raised.


| Item | Suggestion | Evidence in the corpus | Simulator check (answer key) |
|---|---|---|---|
| U1 | worth it: no penalty 50 -> 70 % with FDH | FB-013..017 | YX-006: 70 % 68, 80 % 67, 90 % 66, 100 % 68 g/L |
| U2 | risk: furans/acetate not measured in P5 | v1 method | Not simulated. |
| U3 | risk: C5-2503 formic acid | CoA | YX-006 on C5-2503: 60 % 68, 70 % 67, 80 % 66 g/L |
| U4 | titre stays xylose-limited | C5 feeds topped up to 140 g/L | See design_checks: titre flat with C5 fraction at constant xylose top-up. |


Penalised: (1) Claims the C5 fraction is the main limit on titre. (2) Ignores that FB-012 lacked FDH.


_Evidence:_ corpus/spreadsheets/MD_fermentation_master_2025.xlsx, corpus/meetings/minutes_2025-05-30.md, corpus/sops/SOP-ANA-003_v1.md, corpus/materials/CoA_C5-2503.md, corpus/sops/SOP-FB-001_v3.md  
_Grading:_ rubric (key points: 70 % = 50 %, formate removed, furans not measured, C5-2503)  
_Interop:_ versioning, vocabulary
