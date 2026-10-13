# Deviations from PREREGISTRATION.md

## D1. 30 September 2026: censoring time of event-free patients

- Registered: event-free patients censored at tau + 0.5 in the MIMIC frames;
  INSPECT frame truncated at 365 days.
- Changed to: censored at the end of day tau, i.e. tau + 1 (31 and 366 days);
  INSPECT frame truncated at 366 days.
- Why: deaths on day tau are also placed at tau + 0.5, so the censoring time
  coincided with them. Tested with scikit-survival 0.22.2 and scikit-learn
  1.3.2 on simulated data shaped like Frame A: the time-dependent AUC at the
  horizon could not be computed (the evaluation time equalled the end of
  follow-up, and just below it the censoring survival function was zero), and
  Uno's C, which truncates strictly below its limit, dropped the day-tau
  deaths. With censoring at tau + 1, the time-dependent AUC at 30.5 days equals
  the binary AUROC exactly and Uno's C equals Harrell's C.
- Effect: no event indicator and no count in section 2 changes. Horizons are
  still read as time <= h + 0.5 (MIMIC) or time <= h (INSPECT); C-index
  truncation is at h + 1 (MIMIC) or just past h (INSPECT). The INSPECT frame
  now also holds deaths between days 365 and 366, which fall outside every
  horizon.
- Made before any model was fitted.

## D2. 30 September 2026: bias in pooled out-of-fold discrimination

- Found: in the stage 6a sanity check, two logistic CTPA scores for CV
  readmission by day 10 in the Frame A three-modality cohort (21 events) had
  pooled out-of-fold AUROCs of 0.28 and 0.34. The follow-up check
  (06a_foldbias.py) showed this is the pooling bias described by Parker,
  Gunter and Bedo (BMC Bioinformatics 2007;8:326), not a sign error: the
  average of the per-fold AUROCs was 0.43 and 0.45 (per-fold SD 0.16), and a
  test fold's event rate correlated -0.78 and -0.75 with its mean prediction.
  The matching Cox linear predictor, which has no intercept, showed no gap
  (pooled 0.579, per-fold 0.583, correlation -0.07). The shared folds are
  stratified on the event code at each frame's horizon, so rare short-horizon
  labels are unevenly spread across folds.
- Changed: (1) the stage 6a sanity check "scores ranking no better than
  chance" is replaced by a per-fold rule: a score fails only if the average of
  its per-fold estimates is below 0.40. (2) Every discrimination comparison in
  S1, S2 and S4 is also reported as the average of per-fold estimates, as a
  sensitivity analysis, because pooled evaluation penalises models with an
  intercept (the binary twins, and Cox absolute risks, which carry a
  fold-specific baseline hazard) relative to Cox linear predictors when
  signal is weak.
- Not changed: the primary contrast, which uses Cox linear predictors in
  Frame B EHR+ECG, and the pooled evaluation, which remains the main
  analysis.
- Made before any performance result was examined.

## D3. 30 September 2026: winsorised level-2 inputs in S2, and a robustness
## check of the primary

- Found: the first S2 run failed its hard checks. The cause-specific Cox
  absolute risks (F3) overflowed for a few patients: risks far outside 0 to 1,
  missing values for two patients in Frame B, and risks that fell as the
  horizon lengthened. The cause is a small number of extreme level-1 ridge-Cox
  linear predictors, which enter the cause-specific model through exp(). The
  same extreme scores are the likely source of the long upper tails of the
  primary's per-period bootstrap hazard ratios.
- Changed: in every S2 level-2 model (all four formulations, so the
  comparison stays fair), each level-1 input is winsorised at its
  training-fold 0.5th and 99.5th percentiles before standardising. The
  pre-registration did not specify level-2 input scaling.
- Added: a sensitivity analysis of the primary with the scores winsorised the
  same way, and with rank-based normal scores. The registered primary result
  is unchanged and remains the primary.
- Timing: made after the failed S2 run, whose F1, F2 and F4 results were
  visible. The change applies identically to all four formulations and was
  triggered by the numerical failure, not by those results.
- Outcome (09a_lp_diagnostic.R): no level-1 score reached the +-50 cap. The
  extreme ECG scores come from one patient whose ECG gives extreme values on
  over a dozen unrelated diagnostic statements at once, almost certainly a
  corrupted recording. The primary ratio was unchanged: registered 1.260
  (1.061 to 1.511), winsorised 1.272 (1.056 to 1.548), rank-based 1.241
  (1.038 to 1.499).

## D4. 30 September 2026: intervals in Frame B, and two evaluation choices

- Registered: influence-function intervals (riskRegression) in Frames B, M
  and I.
- Changed: from stage 11a on, Frame B intervals come from a 2,000-sample
  patient bootstrap (seed 42), matching Frame A's clustered bootstrap. The S3
  horizon curve re-chooses the fusion weights inside every resample, which
  influence functions cannot do, and one method across frames keeps the
  contrasts comparable. With complete follow-up in Frame B, both target the
  same quantities. Adopted when stage 11a was written, before any of its
  results existed; recorded here late.
- Clarified: the D3 winsorisation applies to every level-2 model fitted after
  S2, including the stage 11b calibration stacks.
- Clarified: IPA uses the observed event proportion at the horizon as its
  reference model. With complete follow-up this equals the Kaplan-Meier
  reference.

## D5. 30 September 2026: forest implementation and tuning (stage 6b)

- Registered: a random survival forest and a random forest binary twin at
  level 1, tuned by inner 3-fold CV inside each training fold. No package was
  named.
- Changed: both twins are fitted with randomForestSRC in R, the reference
  implementation of random survival forests: a survival forest with log-rank
  splitting for each of the four stage 6a ridge-Cox targets (death,
  composite, cause-specific CV readmission and cause-specific death), and a
  classification forest with Gini splitting for each binary label. Both twins
  share one package, 500 trees, nsplit = 10 random split points, in-fold
  median imputation and the same 6-point grid (node size 5, 15 or 40; mtry
  ceiling(sqrt(p)) or ceiling(3 sqrt(p))). They are tuned on out-of-bag
  performance: Harrell's C for survival forests, AUROC for classification
  forests. Survival ensembles use a grid of 60 training event times plus, for
  each horizon, the last training event time at or before it, so horizon
  risks are exact while memory stays bounded.
- Why: out-of-bag error is the standard internal validation for bagged
  forests and approximates cross-validation without refitting. Inner 3-fold
  CV plus a final refit would need about four times as many forest fits. One
  package for both twins keeps the splitting engine, tree count and grid
  identical, so the survival-versus-binary comparison within the forest
  family stays matched.
- Made before any forest was fitted.

## D6. 30 September 2026: S6 implementation, and calibration under censoring

- Registered: S6, transfer between INSPECT (Frame I) and MIMIC (Frame M) on
  EHR plus shared CTPA features, evaluated as in section 7 (calibration as
  observed/expected, slope, ICI, E50 and E90).
- Changed: (1) S6 uses ridge-penalised Cox fitted with glmnet in R (alpha 0,
  lambda by cross-validated C) rather than scikit-survival, so that fitting,
  IPCW evaluation (riskRegression::Score) and model updating run in one
  script; ridge-Cox remains the learner, stacked into EHR, CTPA and EHR +
  CTPA models as elsewhere. Forests are not used in S6. (2) In Frames I and
  M, calibration is reported as Kaplan-Meier observed / expected, calibration
  slope, and decile tables of Kaplan-Meier risk against mean predicted risk.
  ICI, E50 and E90 are not computed, because under INSPECT's censoring they
  need a censoring-aware smoother used nowhere else in the pipeline; the same
  measures are used in Frame M for comparability.
- Made before any S6 model was fitted.

## D7. 30 September 2026: sPESI inputs (addendum 2, A6)

- Registered (A6): the six sPESI components taken from the 28 EHR features by
  name, a missing value counting as absent.
- Found: the 28 EHR features contain no oxygen saturation, and their heart
  rate and systolic BP are averages over the patient's history up to index
  discharge, whereas sPESI is scored at presentation. Found by the first
  sPESI run, which stopped at column matching before any outcome was
  analysed.
- Changed: SpO2 and presentation vital signs were pulled from MIMIC-IV on
  BigQuery (stage 3b, 03b_bq_ehr_extra.py), after that script reproduced the
  28 features' own vital signs (heart rate 99.97%, systolic BP 100%) and labs
  (median 99.65%). Main sPESI: heart rate, systolic BP and SpO2 at
  presentation (triage of the ED visit that led to the index admission, else
  that visit's first ED vital signs, else the first ICU value within 24 h),
  falling back to the averaged value where no presentation value exists
  (about 27% of admissions); a value still missing counts as absent. Age,
  cancer, heart failure and COPD come from the 28 features as registered.
- Added sensitivity analyses: averaged vital signs for everyone; and
  presentation values restricted to admissions with all three (complete
  cases, about 73%), with every model re-evaluated on those admissions.
- Made before any sPESI outcome result existed.

## D8. 30 September 2026: the CDR score of survivors (addendum, A5)

- Registered (A5): the CDR score V = alpha - F(Y | x), with Y = min(T, 366).
- Found: the first stage 16a run set F(Y | x) for patients alive at one year
  to the model's predicted one-year risk. Y has a point mass at 366, so its
  cumulative distribution there is 1 and their score should be alpha - 1.
  The bounds remained valid, since conformal coverage holds for any score,
  but were far too conservative: coverage 0.98 against a 0.90 target and a
  median bound of 32 days in the EHR+ECG cohort.
- Changed: survivors' score corrected to alpha - 1 in CDR and in both
  Mondrian CDR methods. Nothing else changes. The first run's outputs are
  archived in results/conformal_v1; CQR, the uncalibrated Cox quantile and
  the covariate-free bound use the same random splits and must reproduce
  them exactly, which the corrected script checks.
- Timing: found from the first run's results (CDR coverage far above
  target), so made after they were seen. The correction implements the
  registered score rather than changing it.

## D9. 30 September 2026: stage 18a preprocessing and penalty grid
## (addendum 2, A10)

- Registered (A10): median imputation and standardisation inside each
  training fold; block penalty factors 1, 4, 16 or 64; Breslow risks for the
  PH arms.
- Found: the first stage 18a run gave one or two patients risks of exactly
  0 or 1 in most arms fitted directly on features, EHR alone included, and
  in one fold of one seed the EHR-unpenalised arm returned no risk for 294
  patients (Frame A, three-modality cohort). The late arms, whose inputs are
  winsorised modality scores (D3), were unaffected. Cause: a few extreme
  feature values, in the EHR as well as the ECG block, gave linear
  predictors large enough to overflow the Breslow calculation. Separately,
  the inner cross-validation chose the largest block penalty (64) in a
  third to a half of folds in the EHR+ECG cohort.
- Changed: (1) each non-binary feature is winsorised at its training-fold
  0.5th and 99.5th percentiles before standardising, as D3 did for level-2
  inputs; binary features are left untouched. (2) The Breslow baseline is
  computed on linear predictors centred at their training maximum, with
  risks as -expm1(-H): algebraically identical, but it cannot overflow.
  (3) The block penalty grid gains 256. The stage 18a checks now FAIL on
  missing risks and WARN on risks of exactly 0 or 1.
- Timing: made after the first run's checks and a diagnostic of missing and
  saturated risks, and after seeing which block penalties were chosen;
  before any architecture was evaluated or compared. The only performance
  figure seen was the stage 18a sanity check that its EHR-alone model
  matched stage 6a.

## D10. 30 September 2026: stage 18a block scaling, and the
## plain-concatenation arms (addendum 2, A10)

- Registered (A10): Early PH and Early TV on the block-standardised
  concatenation, each block's standardised columns divided by the square
  root of the block's size.
- Found: glmnet standardises every column internally unless told not to,
  which undid the division by the block size. Stage 18a's "Early PH" and
  "Early TV" are therefore ridge models on the plain concatenation of
  standardised features. The other stage 18a arms are unaffected: their
  columns were already standardised, and glmnet applies penalty factors
  after its own standardisation either way.
- Timing: found while the second stage 18a run was still going, before any
  of its results were seen. That run then finished and was evaluated by
  stage 18c before the fix, so results with the plain-concatenation arms
  under the registered labels have been seen.
- Changed: stage 18b fits the registered arms with standardize = FALSE,
  using stage 18a's folds, preprocessing and inner folds, after checking
  that refitting stage 18a's early model the way 18a did reproduces its
  predictions. It then relabels stage 18a's early arms "Early PH plain
  concatenation" and "Early TV plain concatenation". Stage 18a itself stays
  as run.
- Kept as additional arms: the two plain-concatenation arms, already fitted
  and matching the SIDE PROJECT's comparison of plain against
  block-standardised early fusion. They are reported under their own names
  and contrasted with the registered arms, which remain A10's early arms.

## D11. 1 October 2026: tuning ranges at their edges, and how S4 is run

- Found, glmnet (stages 18a, 18b, 18d): in stage 18d, 46% of the
  horizon-specific penalty choices sat at an end of glmnet's default penalty
  path (lambda.min.ratio 1e-4), and the evidence shows it was the
  least-penalised end. For the block-standardised early model in the
  EHR+ECG cohort, every horizon-tuned choice equalled the whole-follow-up
  choice, both at that end, so that contrast was zero by construction. Every
  penalised model whose risks are used directly had calibration slopes
  above 1 (EHR alone 1.08 to 1.23; registered early 1.33 to 1.55; early TV
  1.21 to 1.50), while the unpenalised late stacks were near 1 (Late TV
  0.98 to 0.99): the pattern of shrinkage forced by a truncated path.
- Changed, glmnet: every glmnet model in stages 18a, 18b and 18d uses a
  longer path, lambda.min.ratio 1e-7 with nlambda 80 (about 11 values per
  decade, as before), and each stage reports the share of choices at either
  end. The first runs' predictions and evaluation are archived
  (data/processed/fusion_arch/v1, results/arch_v1) and reported as
  superseded; the reruns are A10's and A11's results.
- Found, forests (stage 6b): out-of-bag tuning chose the largest node size
  (40) in 58% to 75% of folds for the ECG and CTPA survival forests and all
  three classification forests, and the smaller mtry in 78% to 91% of folds
  for every forest: the most smoothed corner of the D5 grid. The checks
  file holds no performance numbers, so this was seen blind.
- Added, forests: a sensitivity run of stage 6b with node sizes 5, 15, 40,
  80 and 160 and mtry ceiling(sqrt(p) / 2), ceiling(sqrt(p)) and
  ceiling(3 sqrt(p)), everything else as D5. The registered forests remain
  S4's primary forests; the extended-grid forests are a sensitivity
  analysis of S4. If its timing probe projects more than about 10 hours,
  the grid is trimmed before any fit is kept, and recorded here.
- S4 implementation (stage 11c), registered as "learner families ...
  Expected: no gain over ridge-Cox": per frame, cohort and modality, the
  survival family (ridge-Cox, random survival forest, their rank average)
  and the binary family (L2 logistic, random forest, their rank average),
  for death, composite and CV readmission at every horizon. Measures: AUROC
  at each horizon (the cumulative/dynamic AUC under complete follow-up);
  Harrell's C over the follow-up (equal there to Uno's C); and IPA,
  log-likelihood skill score, observed/expected and calibration slope for
  every model with absolute risks. ICI, E50 and E90 are left to stage 11b.
  Fused, death: for each survival family, a Cox stack per seed and fold on
  rank-based normal scores of the modality scores, the same transform for
  every family, because forest ensemble mortality is not on a log-hazard
  scale. Contrasts: forest and rank average minus the linear learner, with
  the 2,000-sample bootstrap (D4) and the D2 per-fold sensitivity analysis.
  The neural learners (stage 10) are added when they exist. Also reported:
  the share of stage 6a and 6b tuning choices at each end of their grids.
- Timing: made after stage 18d's results and stage 6b's blinded checks;
  before any forest performance was seen, before S4 was run, and before
  any rerun.
- Outcome (glmnet reruns): the share of choices at the least-penalised end
  fell from 46% (stage 18d) to 0% for the early models; EHR alone, 19%, sits
  at an end that is now effectively unpenalised. The early models'
  calibration slopes moved towards 1 (registered early, Frame B EHR+ECG:
  1.33 to 1.55 before, 1.14 to 1.35 after). In Frame A's three-modality
  cohort, for the full set, plain concatenation, block penalties and EHR
  unpenalised became unstable: with 338 features and about 120 training
  deaths, near-unpenalised fits often failed to converge (glmnet warnings,
  including 11 empty models), 40% to 43% of choices sat at the
  least-penalised end, and calibration slopes fell to 0.70 to 0.92. Those
  three models in that cell are reported as numerically unstable and not
  interpreted; the registered and late models in the cell were unaffected.

## D12. 1 October 2026: S6 penalty path and feature winsorising
## (addendum 2, A7, A8; stage 16b)

- Registered (D6, A7, A8): ridge-Cox per modality with glmnet (alpha 0,
  lambda by cross-validated C) on imputed features, with glmnet's default
  penalty path.
- Found: in MIMIC, D9 showed that a few extreme feature values distort
  ridge models (winsorising alone raised the EHR model's 30-day AUROC from
  0.855 to 0.870 in Frame A), and D11 showed that glmnet's default path
  (lambda.min.ratio 1e-4) can stop short of the penalty the inner
  cross-validation prefers. Stage 16b uses that default path and
  unwinsorised features, and never recorded where its choices fell.
- Changed: (1) every glmnet model in stage 16b uses the D11 path (80 values
  down to 1e-7 of the largest); (2) after imputation by each A8 strategy,
  every non-binary feature is winsorised at the source (training) data's
  0.5th and 99.5th percentiles, applied to source and target alike;
  (3) the share of penalty choices at each end of the path is recorded and
  checked (WARN above 20% at the least-penalised end). Lambda is still
  chosen by cross-validated C, as D6 registered.
- The first S6 run is archived in results/s6_v1 and reported as superseded.
- Timing: made after S6's first results were seen, prompted by D9 and D11
  rather than by S6's own numbers; before the rerun and before anything in
  A12 was run.
  ## D13. 1 October 2026: how the wide-grid forest rerun is cut if too slow
## (D11)

- Registered (D11): if the wide-grid rerun's timing probe projects more
  than about 10 hours, the grid is trimmed before any fit is kept.
- Changed: the grid is never trimmed, because its new ends (node sizes 80
  and 160, mtry ceiling(sqrt(p) / 2)) are the point of the rerun. Instead,
  if the full run is projected above 10 hours, its scope is reduced to the
  largest of: (1) all four survival targets and every binary label; (2)
  all four survival targets and the death labels; (3) the death survival
  target and the death labels; (4) the death survival target only. The
  scope used is printed, written to the checks file and reported. Death is
  S4's main endpoint and the only one in its fused models.
- Made before the rerun's probe was run.
## D14. 3 October 2026: failed fits in stage 21a (A24)

The first full run of stage 21a stopped before producing any results
because one level-1 fit failed inside the development years, and the
script did not record which. Failures are now handled by a fixed rule,
decided before any result was seen: a forest tuning candidate that fails
is skipped and the best remaining candidate used (as in stages 6c-6g); a
spline Cox fit that fails is retried with smaller bases (k = 4, then 3)
and, if every attempt fails, that fold uses the ridge-Cox score. Every
failure, its message and every substitution are listed in
temporal_failures.csv and counted in the checks; the run stops only if a
model family has no scores at all.
## D15. 3 October 2026: cross-fitted test-year scores in stage 21a (A24)

Investigating stage 21a's calibration showed that only the stacks built on
the oblique forest over-predicted test-year risk (observed/expected 0.48
at 10 and 30 days, against 0.80-1.15 for every other model). The cause is
a property of the software, not the data: aorsf's predicted mortality
grows with the number of distinct event times in the training data, so
learners refitted on all development patients scored the test years on a
larger scale than the out-of-fold development scores the stacks were
fitted to. Test-year patients are now scored with the mean of the five
development-fold models (cross-fitting), so every score comes from a
model trained on the same share of the data; nothing else changes. Both
sets of results are reported; the third version's are kept in
results/temporal_v3/. The main analyses are unaffected, since all their
scores come from fold models.
## D16. Spline Cox on EHR+ (stage 6j): ridge penalty on the linear terms

In stage 6j's first run, 25 of 60 spline Cox fits failed at every basis
size with "indefinite penalized likelihood" (24 in the three-modality
cohort) and fell back to ridge-Cox under D14. A simulation reproduced the
error when a rare indicator had no deaths, or only deaths, in a training
fold (quasi-separation; for example discharge to hospice), which leaves
its unpenalised coefficient without a finite estimate. Change, made
before any EHR+ result was evaluated: in the EHR+ spline Cox family only,
the terms entered linearly carry a ridge penalty (mgcv paraPen, strength
estimated by REML with the smooths). Only that family is refitted
(--refit-spline); ridge-Cox and the forests from the first run are reused
unchanged.
## D17. Presentation vital signs available for 73.5% of admissions (A32)

Stage 3e's check expected at least 90% of admissions to have a first vital
sign from 24 hours before to 24 hours after admission; 73.5% did. MIMIC-IV
holds no ward observations: vital signs come only from the emergency
department (MIMIC-IV-ED, 2011-2019) or the ICU, so direct, transferred and
pre-2011 admissions without early ICU care have none. First labs were
available for 99.8%. Decided before any at-diagnosis model was fitted:
proceed as A32 specified (missing vital signs imputed for the models, with
a "no presentation vital signs" indicator; scored 0 in sPESI and PESI),
and add a sensitivity analysis repeating the sPESI and PESI comparison in
the admissions with presentation vital signs.
## D18. The first ECG is identified by date (A32)

Stage 0q found every ECG time in ecg_record_index.csv at exactly midnight,
and the admission times in index_admission_times.csv likewise, so ECG
order within a day cannot be recovered and A32's 24-hour ECG window cannot
be applied as specified. Changed before any at-diagnosis model was fitted:
the ECG at presentation is the ECG(s) dated the day before or the day of
admission (an ECG on the admission day is at most about 24 hours after
admission); where several share the earliest such date, their features are
averaged (axes as circular means). The 95 "_worst" columns of ecg_opt.csv
are left out (for a single ECG they duplicate their "_mean" columns); the
3 "_circ" axis columns are mapped to the per-ECG axes. CTPA report hours
are measured from ctpa_notes_index.csv's own admission times, which carry
the clock time. Blood tests and vital signs are unaffected (stage 3e used
MIMIC's admission times).
D19 (replaces D18). The stage 0r audit recovered real ECG acquisition times from MIMIC-IV-ECG on BigQuery. Therefore, the at-diagnosis ECG is now the first ECG recorded from 24 hours before to 24 hours after admission, the same window as the labs and vital signs. D18’s date rule had dropped 1,806 ECGs dated the day after admission, and 55% of these were taken within 24 hours. ECG coverage rises from 90.9% to 96.6% of admissions. Where D18 averaged several ECGs from the first date, D19 uses the first ECG only.