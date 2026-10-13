# Pre-registration addendum (P2)

Analyses specified after results were seen. Each entry is committed before
any of its own analyses is run. All are exploratory: intervals only, no
decision rules, no claims of significance. Changes after commit go in
DEVIATIONS.md.

## A1. 30 September 2026: the ECG hazard ratio as a curve over time

Specified after: the primary result (ECG hazard-ratio ratio, days 0-30 over
days 31-365, 1.260, 95% CI 1.061 to 1.511), its hospice and extreme-score
sensitivities (D3: winsorised per-period hazard ratios 1.53 and 1.20 per SD),
S1, S2, S3, S5 and the calibration stage. Committed before any analysis below
was run.

Question: after discharge, does the ECG's prognostic weight for death decline
smoothly, or change at one point?

Data: Frame B, death to 365 days after discharge. Seed-averaged out-of-fold
ridge-Cox death scores from stage 6a. Main analyses use the scores winsorised
at their 0.5th and 99.5th percentiles and then standardised, because the
registered per-SD scale is inflated by one corrupted ECG recording (D3). The
registered scaling is a sensitivity analysis. Main cohort: EHR+ECG
(n = 2,909, 620 deaths).

A1.1 Flexible parametric (Royston-Parmar) model, R package rstpm2
- Proportional-hazards link; baseline log cumulative hazard a natural spline
  in log time with 4 df; a time-varying effect with 3 df for each score.
- Output: hazard ratio per 1 SD of each score at days 1 to 365, for a patient
  with average scores (all scores 0), with 95% pointwise delta-method bands.
- Summary, main model only: hazard ratio at 7, 14, 30, 90, 198 and 365 days,
  and the ratio HR(14 d) / HR(198 d), the midpoints of the primary's two
  periods, with 1,000-sample patient bootstrap percentile intervals
  (seed 42). Likelihood-ratio test of each score's time-varying effect (3 df).
- Sensitivity: baseline df 3, 4 or 5 crossed with time-varying df 1 to 4
  (AIC and key hazard ratios); registered scaling; hospice discharges
  excluded; three-modality cohort with the CTPA score added.

A1.2 Shape of the change (Cox models, data split at every event time)
- Compared by AIC: constant effects; a step at day 30 (the primary's model);
  effects linear in log time; effects as a natural spline in log time (3 df).
- Change-point profile: the ECG effect allowed to step at day 7, 14, 21, 30,
  45, 60, 90, 120, 180 or 270 (EHR effect constant), reporting the partial
  log-likelihood and the hazard ratios either side of each step. A sharp peak
  would point to a step; a flat profile to a gradual change.

A1.3 Landmark supermodel (van Houwelingen 2007)
- Landmarks every 15 days from 0 to 330 days after discharge. At each
  landmark s, the patients still alive are followed for deaths on days s to
  s + 30 and administratively censored at s + 31 (as in D1).
- One stacked Cox model stratified by landmark, with each score's log hazard
  ratio a quadratic in landmark time,
  theta0 + theta1 (s / 330) + theta2 (s / 330)^2, and robust variance
  clustered by patient.
- Output: the hazard ratio per 1 SD for death in the next 30 days as a
  function of time since discharge, with 95% bands; landmark-by-landmark
  estimates for comparison; the Wald test of theta1 = theta2 = 0 (2 df); and
  the ratio HR(s = 0) / HR(s = 180).
- Sensitivity: a 60-day window (landmarks 0 to 300); landmarks every 30 days;
  registered scaling; hospice discharges excluded; three-modality cohort with
  CTPA.

## A2 to A4. 30 September 2026: clinical utility, stability, ECG features

Specified after everything listed under A1 and the A1 results (the ECG's
hazard ratio falls mostly within two to three weeks of discharge; a smooth
decline in log time fits better than a step at day 30). Committed before
any analysis below was run.

A2. Decision curve analysis
- Net benefit at threshold pt, treating everyone whose predicted risk is at
  least pt: TP/n - FP/n x pt / (1 - pt). Follow-up is complete, so the
  outcome is whether the event happened by the horizon; for CV readmission a
  death first counts as no event (van Geloven et al., BMJ 2022). References:
  treat all and treat none.
- Death, Frame A 30 d and Frame B 30 and 365 d, both cohorts: EHR alone
  (level-1 ridge-Cox absolute risk), Cox stack, binary stack and late fusion
  recalibrated (stage 11b predictions, out of fold, seed-averaged).
- CV readmission, same frames, horizons and cohorts: the four S2
  formulations.
- Curves at thresholds 1% to 50%, with net reduction in interventions per
  100 patients. At 5%, 10% and 20%, each model's net benefit and these
  differences, with 2,000-sample bootstrap intervals (patient-clustered in
  Frame A), seed 42: Cox stack - EHR alone, binary stack - EHR alone, Cox
  stack - binary stack; F2 - F1, F3 - F1, F4 - F1.
- Expected: F1's overprediction lowers its net benefit at thresholds near
  and above the observed risk.

A3. Stability (Riley and Collins, Biometrical Journal 2023) and subgroups
- Seed instability: across the six seeds' out-of-fold Cox-stack and
  binary-stack death predictions, each patient's mean absolute deviation
  from the seed average (overall and by risk decile), and the share of
  patients whose classification at 10% or 20% differs between seeds.
- Bootstrap instability of the fusion stage: the Cox stack and binary stack
  refitted, including winsorising and standardising, to the seed-averaged
  out-of-fold level-1 scores in the full cohort (original model) and in
  1,000 bootstrap samples (seed 42; patient-clustered in Frame A), each
  applied to the original sample. Reported: mean absolute prediction error,
  classification instability at 10% and 20%, ICI of each bootstrap model,
  net benefit at 5%, 10% and 20%, and instability-plot data. Death, Frame A
  30 d and Frame B 30 and 365 d, both cohorts. Level-1 models are held fixed,
  so this measures the instability of the fusion stage only.
- Subgroups: Frame A death 30 d and Frame B death 365 d, EHR+ECG, Cox stack
  and binary stack: observed and mean predicted risk, O/E, AUROC, calibration
  slope and ICI by age band (<65, 65-79, 80+) and by sex if an EHR feature
  records it, with 1,000-sample bootstrap intervals. Recorded race is not in
  the local data and is not analysed here.

A4. Which ECG features carry the early signal
- Frame B EHR+ECG, death to 365 d. Each of the 264 ECG features (71
  SCP-statement logits, 193 measurements) winsorised at its 0.5th and 99.5th
  percentiles and standardised; complete cases per feature.
- Main: a Cox model with the feature's coefficient split at day 30 (the
  primary's periods), adjusted for the EHR score (also split). Reported:
  hazard ratio per SD early and late, their ratio with Wald 95% CI and p,
  and Benjamini-Hochberg q across the 264 features; also the whole-year
  hazard ratio. Each feature is oriented so that higher values mean higher
  one-year risk, so a ratio above 1 always means front-loaded.
- Sensitivity: the effect linear in log time (A1.2's best-fitting shape),
  slope with p and q; the main model without EHR adjustment.
- Summaries: counts of front-loaded and late-loaded features (q < 0.05) by
  block (SCP statements, measurements), and the 20 most front-loaded.
- Expectation, not formally tested: statements linked to right-ventricular
  strain (sinus tachycardia, right bundle branch block, right ventricular
  hypertrophy, right atrial overload, T-wave inversion) will be among the
  most front-loaded.

## A5. 30 September 2026: conformal lower bounds on survival time

Specified after everything listed under A1 to A4 and their results.
Committed before any analysis below was run.

Method: Candes, Lei and Ren (JRSSB 2023), Algorithm 1. Frame B follows
every patient to day 366, a censoring time observed for everyone and
independent of outcome and covariates, so the method applies with threshold
c0 = 366 and all weights equal to 1; the target is Y = min(T, 366).
- Data: Frame B, death to one year after discharge, EHR+ECG cohort (main)
  and three-modality cohort; one row per patient.
- Base models: the Cox stack (level-2 Cox on the winsorised ridge-Cox scores
  of every modality) and an EHR-only Cox model, each fitted on the four
  training folds of each seed; predicted survival from the Breslow baseline.
- Calibration: within each held-out fold (whose predictions come from models
  that never saw it), 10 random half splits, each half calibrating the other:
  600 calibration/test splits per cohort and level.
- Scores: CDR (main), V = alpha - F(Y | x), bound = earliest time at which
  predicted risk reaches alpha - eta; CQR (sensitivity), V = q_alpha(x) - Y,
  bound = q_alpha(x) - eta; eta = the ceil((1 - alpha)(n + 1))-th smallest
  calibration score. Bounds clipped to [0, 366].
- Comparators: the Cox model's own alpha-quantile without calibration; a
  covariate-free conformal bound. Group-conditional (Mondrian) CDR within
  tertiles of predicted one-year risk (cut-points from the training folds)
  and within age bands (<65, 65-79, 80+).
- Levels: alpha = 0.10 (main), 0.05, 0.20.
- Reported: coverage (mean, 5th and 95th percentiles across splits); mean
  and median bound; share bounded at 366 days and below 30 days; coverage
  and mean bound by risk decile, risk tertile and age band (alpha 0.10); the
  per-patient difference in mean bound, fused minus EHR only, with a
  2,000-sample patient bootstrap interval.
- Expected: marginal coverage at or above 1 - alpha for every conformal
  method (by construction); the uncalibrated Cox quantile may miss it;
  unconditional CDR may under-cover the highest-risk tertile, which Mondrian
  calibration corrects; the fused model gives slightly higher bounds than
  EHR alone.