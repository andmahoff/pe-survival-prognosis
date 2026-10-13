# Registration log

The analyses were developed in a local git repository that was never pushed
to a public server. The table lists its 87 commits in order, so a reader can
see when each registration document, addendum (A1 to A38) and deviation (D1
to D19) was committed relative to the analyses it covers. The dates are local
commit times (UK time, UTC+1) and cannot be verified independently; they are
given as a record, not as proof.

The pre-registration (PREREGISTRATION.md) was the first commit, on
30 September 2026 at 10:21, before any model was fitted to a P2 analysis
frame. Addenda were committed with the words "before running" when they
preceded the analysis they specify.

This public repository was created afresh from the final state of that
local repository. Compared with the local copy, the scripts differ only in
the following ways:

- personal folder paths were replaced by the `PROJECTS_DIR` environment
  variable, and the BigQuery billing project by `GCP_PROJECT`;
- a few printed messages and account hints were reworded;
- header and other comments were rewritten for readability.

A token-level comparison (R parse tokens without comments, Python syntax
trees without the module docstring) confirmed that nothing else in the code
changed. The registration documents in this folder are copied unchanged.

| # | Committed (UK time) | Commit | Message |
|---|---|---|---|
| 1 | 2026-09-30 10:21 | `2b72c64` | Pre-registration of P2 analyses, before any model fitting |
| 2 | 2026-09-30 10:36 | `8df37a2` | D1: censor event-free patients at the end of day tau, before any model fitting |
| 3 | 2026-09-30 11:02 | `1c48dea` | Stage 5 gates: all 27 pass, before any analysis |
| 4 | 2026-09-30 13:10 | `16b9987` | D2: pooled out-of-fold discrimination bias; stage 6a scripts, before any evaluation |
| 5 | 2026-09-30 13:19 | `39b2b05` | Stage 7: pre-registered primary contrast (H1 supported) |
| 6 | 2026-09-30 13:25 | `2ecd251` | Stage 7b: primary robust to excluding hospice discharges |
| 7 | 2026-09-30 13:34 | `098dd14` | Stage 8: late fusion of linear scores (S1, S3, S5) |
| 8 | 2026-09-30 13:44 | `825ab3e` | Stage 7c: S5 primary in three-modality cohort |
| 9 | 2026-09-30 14:13 | `b7a15f1` | Stage 11a: bootstrap intervals for S1 and S3 |
| 10 | 2026-09-30 15:14 | `468ae39` | D3 + S2: competing-risk formulations; primary robust to extreme scores |
| 11 | 2026-09-30 15:46 | `814d82a` | Stage 11b: calibration of death and composite predictions (D4) |
| 12 | 2026-09-30 16:04 | `0bc88e3` | Addendum A1: time-varying ECG hazard ratio, specified before running |
| 13 | 2026-09-30 16:19 | `90a5595` | Addendum A1: Royston-Parmar curve, shape comparison, landmark supermodel |
| 14 | 2026-09-30 16:24 | `83eaff7` | Addendum A1 results: Royston-Parmar curve, shape comparison, landmark supermodel |
| 15 | 2026-09-30 16:36 | `907fe0a` | Addendum A2-A4: decision curves, stability, ECG features, before running |
| 16 | 2026-09-30 17:30 | `d3fb7b6` | Addendum A2-A4: decision curves, stability and subgroups, ECG feature timing |
| 17 | 2026-09-30 17:43 | `271c92b` | Addendum A5: conformal lower bounds on survival time, before running |
| 18 | 2026-09-30 18:00 | `4c48669` | D6 and addendum A6-A8: sPESI benchmark, S6 with Debray framework, missing-data strategies, before running |
| 19 | 2026-09-30 18:35 | `6ac22c4` | D7: sPESI inputs from presentation vitals (stage 3b pull), before any sPESI result |
| 20 | 2026-09-30 19:23 | `723bf98` | D8: CDR score of survivors corrected, before re-running stage 16a |
| 21 | 2026-09-30 19:47 | `e645ae4` | Addendum A9: fusion gap mechanism along follow-up time, before running |
| 22 | 2026-09-30 20:25 | `f96ad0d` | Addendum A10: fusion architecture by time-varying effects, before running |
| 23 | 2026-09-30 21:37 | `d15c047` | D9: stage 18a feature winsorising, overflow-proof Breslow, wider penalty grid |
| 24 | 2026-09-30 22:23 | `02c638c` | D10: registered block-standardised early arms (18b); 18a early arms relabelled as plain concatenation |
| 25 | 2026-09-30 22:52 | `6bfa4d3` | Addendum A11: horizon-tuned PH models, before running |
| 26 | 2026-09-30 23:22 | `15e1ff3` | Stage 6b forests (D5), before running |
| 27 | 2026-10-01 10:31 | `fd6f232` | D11: tuning ranges at their edges, forest sensitivity grid, S4 implementation, before running |
| 28 | 2026-10-01 11:36 | `20a1877` | D11 reruns: longer penalty path, edge reporting, 18a/18b split cleaned, before running |
| 29 | 2026-10-01 13:42 | `664f1ac` | D12 (S6 path and winsorising) and A12 (zero-shot forests and horizons), before running |
| 30 | 2026-10-01 14:24 | `9d995b8` | A13 (time-varying-coefficient stack) and scripts for A12 and A13, before running |
| 31 | 2026-10-01 14:47 | `e448039` | 16d: restore R's core setting after the probe tasks; record warnings |
| 32 | 2026-10-01 18:38 | `f7f437b` | A14 (super learner) and A15 (CTPA time pattern and tuning), before running |
| 33 | 2026-10-01 20:53 | `b068785` | D13 and stage 6c (wide-grid forests, S4 sensitivity), before running |
| 34 | 2026-10-02 01:15 | `71ea94a` | 6c: save.memory forests, release memory, checkpoint and resume (first run ran out of memory) |
| 35 | 2026-10-02 09:41 | `ded007a` | 11d: S4 sensitivity evaluation with the wide-grid forests |
| 36 | 2026-10-02 10:23 | `5e97606` | A16 (retuned ridge-Cox, oblique forests, CoxBoost, super learner v2), before running |
| 37 | 2026-10-02 11:05 | `74ea2a8` | 6d: grow oblique forests in batches of 50 trees (first probe ran out of memory) |
| 38 | 2026-10-02 15:39 | `b07ecc4` | 6d: record failed forest fits, refit affected tasks on resume, FAIL if any remain |
| 39 | 2026-10-02 17:50 | `79721fa` | A17 (registered analyses on retuned ridge-Cox), A18 (Analysis B), A19 (spline Cox), before running |
| 40 | 2026-10-02 19:26 | `079dfb8` | 13r: compare decisions and TRUE/FALSE flags instead of skipping files where they changed |
| 41 | 2026-10-02 19:28 | `3989220` | 13r: handle results files with no number columns |
| 42 | 2026-10-02 21:48 | `092856e` | A20-A22: binary retune, spline and tuning changes, Super learner v4, before running |
| 43 | 2026-10-03 14:04 | `73c5f0e` | A23 (interpretability) and A24 (temporal validation), before running |
| 44 | 2026-10-03 20:50 | `a1d7b4c` | 21a v2: timing probe drops columns constant in the development years |
| 45 | 2026-10-03 21:56 | `a9bc199` | D14: 21a v3 records failed fits and applies a fixed fallback rule |
| 46 | 2026-10-03 23:22 | `fd998eb` | D15: 21a v4 cross-fits test-year scores (oblique-forest scale artefact) |
| 47 | 2026-10-04 00:43 | `822a82e` | A25 (stage 12 sensitivity) and A26 (stage 10 neural learners, A10 arms), before running |
| 48 | 2026-10-04 00:47 | `0bed133` | 12 v2: final-year part uses the composite models present in stage 11b's file |
| 49 | 2026-10-04 01:20 | `be2a5f6` | A26 part 2: stage 11e (S4 neural evaluation), 18i/18j (A10 forest and neural arms), 18k (their evaluation), before running |
| 50 | 2026-10-04 11:28 | `08c8e80` | 11e v2: fix the sign of the log-likelihood skill score |
| 51 | 2026-10-04 13:02 | `85f652e` | A27: block forest optimised (18l), 18k v2, Super learner v5 specified, before running |
| 52 | 2026-10-04 17:32 | `760cdf0` | A27: stage 18m, Super learner v5 with the optimised block forest, before running |
| 53 | 2026-10-04 19:50 | `d1fff59` | A28: EHR 79 level-1 models with edge-extending tuning (6i), CTPA blocks over time (16f), Super learner v6 specified, before running |
| 54 | 2026-10-05 00:16 | `c677f26` | A28: the final-model rule uses AUROC (clarified before stage 18n) |
| 55 | 2026-10-05 00:49 | `54f2b3f` | A28: stage 18n, Super learner v6 and per-horizon weights, before running |
| 56 | 2026-10-05 00:58 | `7adaa09` | A29: subgroup evaluation (3c demographics pull, 22a), with 18n, before running |
| 57 | 2026-10-05 11:15 | `a7db8d6` | A30: EHR+ (diagnosis codes, destination, unused measurements) and Super learner v7 specified, before running |
| 58 | 2026-10-05 11:15 | `f6bc552` | A30 scoped as the discharge model, with 3d, 6j and the inventory for the at-diagnosis model (A31) |
| 59 | 2026-10-05 11:15 | `5ced7f6` | Stage 16c (version that was run) and stage 17a, previously uncommitted |
| 60 | 2026-10-05 14:09 | `aca71f8` | D16: ridge penalty on the EHR+ spline model's linear terms, with a spline-only refit mode |
| 61 | 2026-10-05 14:31 | `e1ab6bd` | A30: stage 18o, Super learner v7 (discharge model) with the hospice-excluded evaluation, before running |
| 62 | 2026-10-05 14:54 | `ea1320c` | 22a v2: subgroup evaluation rerun on Super learner v7 (A30 adopted it) |
| 63 | 2026-10-05 15:08 | `34a12d5` | A31: what the hospice indicator adds (6j --no-hospice, 18p), before running |
| 64 | 2026-10-05 16:51 | `bb61fe4` | A32: at-diagnosis model specified, with the presentation pull (3e) and ECG/CTPA inventory (0p), before running |
| 65 | 2026-10-05 17:00 | `d4a27ae` | 3e v2 (item IDs as numbers or text), inventory script restored, comorbidity install script |
| 66 | 2026-10-05 17:10 | `6534288` | D17 (presentation vital signs 73.5%) and stage 4p, the at-diagnosis builder, before running |
| 67 | 2026-10-05 17:18 | `b0fea23` | 0q: ECG time resolution and unmatched-column check for A32 |
| 68 | 2026-10-05 17:46 | `e0dc91a` | D18: first ECG by date (day before or day of admission), axis columns mapped, CTPA hours from the notes' own admission times |
| 69 | 2026-10-05 19:23 | `c4e3735` | A32: stage 6p, at-diagnosis level-1 fits (EHR, first ECG, presentation CTPA), before running |
| 70 | 2026-10-06 08:56 | `e3822e5` | A32: stage 18q, the at-diagnosis super learner and its sPESI, PESI-10, EHR-alone and Frame A comparisons, before running |
| 71 | 2026-10-06 09:39 | `62fffa1` | A33 (Option B: composite and CV readmission, both pipelines, matched binary models) and stage 4r, before running |
| 72 | 2026-10-06 09:40 | `ca35fff` | A33: secondary dissertation-style CV readmission analysis (death-first excluded per horizon), before any Option B model |
| 73 | 2026-10-06 10:15 | `ab86bc8` | A33: stage 6q, survival level-1 fits for the composite and CV readmission in both pipelines, before running |
| 74 | 2026-10-07 09:18 | `c6b83b5` | A33: stage 18r, composite and CV super learners in both pipelines, v7-core check and dissertation-style CV secondary, before running |
| 75 | 2026-10-07 10:13 | `5fc8ecc` | A34 (no-hospice composite/CV at discharge, component-based composite) and stage 6q v2 with --no-hospice, before running |
| 76 | 2026-10-07 16:45 | `5f51046` | A34: stage 18r v2, no-hospice composite/CV at discharge, hospice-excluded evaluation, component-based composite, before running |
| 77 | 2026-10-07 17:44 | `df1ecb3` | A33: stage 6r, matched binary level-1 models (ridge logistic, probability forests), both pipelines, before running |
| 78 | 2026-10-09 08:27 | `1a0b215` | A33: stage 18s, binary super learners and survival-versus-binary comparison, before running |
| 79 | 2026-10-09 14:19 | `0bd406f` | A35: stage 0r diagnosis-window audit, stage 13a Table 1, stage 21b temporal validation, before running |
| 80 | 2026-10-09 23:50 | `f05b107` | A36 (history known on arrival as the main at-diagnosis model) and D19 (first ECG by the clock within 24 h): stage 4p v3, 18r v3, 13a v2, before running |
| 81 | 2026-10-11 19:05 | `fc524b5` | Stage 13a v3: fixed decimals per row, complementary 30-day small-cell suppression, Portuguese grouped as White (stage 22a) |
| 82 | 2026-10-11 19:50 | `3ca57fc` | A37: interpretability (20p, 20q) and CTPA block split (16p) for the A36 at-diagnosis model, before running |
| 83 | 2026-10-12 12:32 | `5e25362` | Stage 20q v2: clearer lab names, identical Charlson/Elixhauser indicators merged in importance figures |
| 84 | 2026-10-12 12:42 | `c70120a` | Stage 20q v3: labelled log-scale ticks, share figures from 0 |
| 85 | 2026-10-12 13:21 | `56c4d99` | Stage 13b-d: final-model tables, calibration and decision curves (A38) |
| 86 | 2026-10-12 13:54 | `8402c77` | Stage 13b-d: final-model tables, calibration and decision curves (A38) |
| 87 | 2026-10-12 21:33 | `18825c3` | Stage 13e: verified results file (A38) |
