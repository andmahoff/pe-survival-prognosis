# Pre-registration: time-to-event re-analysis of multimodal PE prognosis (P2)

Author: Andrzej Machowski
Registered: 30 September 2026, before any model was fitted to a P2 analysis frame.
Record: the first git commit containing this file.

## 1. What had been seen before registration

- Outputs of scripts 01 to 04b: cohort sizes, event counts at each horizon,
  data-quality checks and the INSPECT censoring pattern. No model had been
  fitted to any P2 frame.
- From earlier work on overlapping data: the DISSERTATION and SIDE PROJECT
  binary 30-day results for these MIMIC cohorts, and the SIDE PROJECT
  INSPECT/MIMIC EHR+CTPA transfer results. The hypotheses below were written
  knowing those results.

## 2. Analysis frames

| Frame | Data | Clock starts | Cohort | Horizon | n | Events |
|---|---|---|---|---|---|---|
| A | MIMIC-IV | index admission | all index admissions, EHR+ECG (three-modality subset) | 30 d | 3,501 (1,636) | death 397, CV readmission first 172, composite 547 |
| B | MIMIC-IV | index discharge | first admission per patient, alive at discharge | 365 d | 2,909 (1,405) | death 620, CV readmission first 431, composite 899 |
| M | MIMIC-IV | index admission | first admission per patient, three-modality | 365 d | 1,490 | death 369 |
| I | INSPECT | index CTPA scan | PE-positive persons | 365 d | 3,300 | death 539, 475 censored before one year |

- Frame B excludes in-hospital deaths and patients with an in-window ECG
  recorded after discharge, so every feature is known when the clock starts.
- Three-modality subsets are nested: the patient's first admission must itself
  have a CTPA.
- CV readmission: a later admission with a harmonised MACE code (MI, heart
  failure, cardiac arrest, AF, ischaemic stroke; DISSERTATION script 80 codes).
- Competing risks: 1 = CV readmission first, 2 = death first. A CV readmission
  and a death on the same day count as CV.
- Time: MIMIC is day-level. An event on day d sits at d + 0.5, censoring at
  tau + 0.5, and horizon h is read as time <= h + 0.5. INSPECT is minute-level
  and horizon h is read as time <= h, reproducing its native labels.
- Horizons: 10, 30, 90, 180 and 365 days (6 months = 180 days).

## 3. Modalities

- EHR: 28 features.
- ECG: 264 features (71 SCP-ECG statement logits, 193 waveform measurements,
  strict Tp-e only).
- CTPA: 46 report features; 38 are shared with INSPECT (the device block is
  MIMIC-only).
- Feature lists are fixed in data/processed/feature_whitelist*.csv. No label
  or timing column is a feature.

## 4. Models and validation

- Level 1, per modality, out-of-fold: ridge-penalised Cox, random survival
  forest, and their rank average. Cause-specific versions for CV readmission.
- Binary twins with matched tuning budgets at each horizon: L2 logistic
  regression, random forest, and their rank average.
- Neural comparators (secondary): DeepSurv, DeepHit and an MLP binary twin.
- Validation: 5-fold StratifiedGroupKFold x 6 seeds (42, 7, 13, 1, 2, 3),
  grouped by patient, folds fixed in data/processed/folds*.csv. Tuning by inner
  3-fold CV inside each training fold. Imputation and scaling fitted inside
  training folds only. Predictions averaged over seeds before evaluation.
- Fusion: late fusion by weighted rank average with weights tuned inside
  training folds; stacked Cox, cause-specific and Fine-Gray models on
  out-of-fold modality scores; early fusion by ridge-Cox on block-standardised
  concatenated features and by a random survival forest on the concatenation.

## 5. Primary contrast

- Frame B, EHR+ECG cohort, all-cause death to 365 days after discharge.
- A Cox model on the standardised, seed-averaged, out-of-fold ridge-Cox risk
  scores of the EHR and ECG modalities, with each score's coefficient allowed
  to differ between days 0 to 30 and days 31 to 365 after discharge (split at
  30.5).
- Estimand: the ratio of the ECG hazard ratio in days 0 to 30 to the ECG hazard
  ratio in days 31 to 365.
- Uncertainty: 95% percentile interval from 2,000 patient-level bootstrap
  refits of that model.
- H1: the ratio is above 1 (the ECG's prognostic information is concentrated
  in the first month).
- Decision rule: H1 is supported if the lower bound exceeds 1. The reverse is
  reported if the upper bound is below 1. Otherwise the result is inconclusive.
  The EHR ratio is reported alongside with no decision rule.

## 6. Secondary analyses (exploratory: intervals only, no significance claims)

- S1, formulation factorial (Frame A, all three endpoints at 30 days):
  per-modality objective (cross-entropy vs Cox partial likelihood) crossed with
  fusion criterion (AUROC at 30 days vs Uno's C to 30 days). Outcomes: fused
  time-dependent AUC, Uno's C and modality weights. Expected: the survival and
  binary formulations agree within noise at 30 days.
- S2, CV readmission formulations (Frame A to 30 days; Frame B to 365 days):
  conditional binary (death-first excluded, the DISSERTATION definition),
  unconditional binary, cause-specific Cox and Fine-Gray. Judged on absolute
  CV readmission risk in the full cohort and on the conditional estimand.
  Expected: the conditional binary model is worst on the full-cohort estimand.
- S3, horizon curve: in-fold late-fusion ECG weight at each horizon, Frame B
  death. Expected: declining with horizon.
- S4, learner families including DeepSurv and DeepHit. Expected: no gain over
  ridge-Cox.
- S5, three-modality cohort: the primary contrast and S1 repeated with CTPA.
- S6, INSPECT and MIMIC (Frames I and M, EHR plus shared CTPA, death to 365
  days): transfer in both directions, discrimination and calibration reported
  separately, and the bias from dropping censored patients quantified against
  Kaplan-Meier.

## 7. Evaluation

- Harrell's C; Uno's C truncated at the horizon; cumulative/dynamic
  time-dependent AUC at each horizon; Brier score and IPA; calibration at 30 and
  365 days (observed/expected, slope, ICI, E50, E90).
- Paired contrasts on the same patients: influence-function intervals
  (riskRegression) in Frames B, M and I; a 2,000-resample patient-clustered
  bootstrap in Frame A.
- Late-fusion outputs are ranks, so their calibration is assessed only after
  in-fold recalibration.

## 8. Sensitivity analyses

- Hospice discharges excluded (Frame B).
- Index possibly in the final data year excluded (Frame B CV readmission and
  composite).
- All admissions with a patient-clustered bootstrap (Frame B).
- 48-hour landmark (Frame A).
- INSPECT one-year mortality bounds (all censored survive; all die).
- sPESI as a single-covariate Cox benchmark, if its inputs are available.

## 9. Gates before any analysis (stage 5)

- No label or timing column among the features.
- Binary per-modality models reproduce the SIDE PROJECT MIMIC results within
  0.005 where a MIMIC-trained reference exists.
- In Frame A, time-dependent AUC at 30.5 days equals the ordinary death_30d
  AUROC to three decimals.
- Every risk score has C above 0.5; a permuted score gives about 0.5.
- Harrell's and Uno's C agree in Frame A.

## 10. Deviations

Any change made after this commit is logged in DEVIATIONS.md with its date,
what changed and why, and is reported in the write-up.