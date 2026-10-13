# Pre-registration addendum 2 (P2)

Continues PREREGISTRATION_ADDENDUM.md. Same rules: each entry is committed
before any of its own analyses is run; all are exploratory (intervals only,
no decision rules); changes after commit go in DEVIATIONS.md.

## A6. 30 September 2026: the sPESI benchmark in detail

Specified after A1 to A5 and their results. The pre-registration lists
"sPESI as a single-covariate Cox benchmark" among the sensitivity analyses;
this entry fixes how it is run.
- sPESI (0 to 6): age over 80; cancer; chronic cardiopulmonary disease
  (heart failure or chronic lung disease, one point); heart rate 110 or
  more; systolic blood pressure under 100; oxygen saturation under 90%.
  Components are taken from the 28 EHR features by name; a missing value
  counts as absent. The chosen columns are printed and saved before any
  outcome is analysed.
- Cells: Frame A death 30 d; Frame B death 30 and 365 d; both cohorts.
- sPESI as a single-covariate Cox model fitted in the stage 4 training folds
  (out of fold, seed-averaged), against EHR alone (level-1 ridge-Cox risk),
  the Cox stack and the binary stack (stage 11b predictions).
- Measures: AUROC, Brier, IPA, observed/expected, calibration slope, net
  benefit at 5%, 10% and 20%. As a rule (sPESI 1 or more = not low risk):
  sensitivity, specificity, PPV, NPV, the share classed low risk and their
  death rate, and its net benefit at the same thresholds.
- Matched comparisons, for each model: the death rate in its lowest-risk
  group of the same size as sPESI's low-risk group, and its specificity at
  the sensitivity sPESI achieves.
- Intervals: 2,000-sample bootstrap, patient-clustered in Frame A, seed 42,
  for every measure and for the differences model minus sPESI.

## A7. 30 September 2026: S6 read through the Debray framework, with updating

Specified after A1 to A5 and their results, before any S6 model was fitted
(implementation choices in D6).
- S6 as registered: ridge-Cox per modality and the EHR + CTPA stack, trained
  on all of one dataset and applied to all of the other, both directions;
  internal 5-fold cross-validation in each dataset as the reference; death
  at 30 and 365 days; discrimination and calibration reported separately;
  the bias from dropping censored INSPECT patients quantified against
  Kaplan-Meier.
- Step 1, relatedness (Debray et al., J Clin Epidemiol 2015): out-of-fold
  AUROC of ridge-logistic membership models telling the datasets apart (all
  features; all features with missing indicators; EHR only; CTPA only);
  per-feature standardised mean differences and missingness; mean and SD of
  the transferred fused score in source and target; Kaplan-Meier mortality.
- Step 2, performance in the target: time-dependent AUC and Brier / IPA with
  IPCW and influence-function intervals and paired contrasts; Uno's C;
  observed / expected; calibration slope; decile calibration tables.
- Step 3, updating, fused model, main imputation, by 5-fold cross-validation
  in the target: baseline hazard re-estimated; baseline and slope
  re-estimated; full refit (the internal model).
- Expected: the datasets are clearly distinguishable (membership AUROC well
  above 0.5); ranking transfers better than absolute risk; updating the
  baseline recovers most of the calibration.

## A8. 30 September 2026: missing-data strategies for the transfer

- Four strategies, each fitted on the source only and applied to both: source
  medians (the pipeline's method, main); source medians with missing
  indicators; regression imputation fitted on the source without the outcome
  (Sisk et al., Stat Methods Med Res 2023); the target's own medians (a
  deployment site imputing with its own data).
- Each strategy is carried through the whole transfer and evaluated as in A7,
  with paired contrasts against the main strategy.
- Expected: small differences overall; missing indicators may help or harm
  depending on whether missingness relates to the outcome, and target
  medians may improve calibration in the target.

## A9. 30 September 2026: the fusion gap mechanism along follow-up time

Specified after A1 to A8 and all their results, including the S3 horizon
curve, and before any analysis below was run.
- Question: does the relation between the unimodal gap and the fusion gain,
  found across the SIDE PROJECT's 12 transfer and outcome cells (r = -0.902)
  and in the DISSERTATION's script 140, also hold along follow-up time; and
  is it what the theory of combining markers predicts (Su and Liu, JASA
  1993; Pencina et al., Am J Epidemiol 2012)?
- Data: Frame B, death after discharge, both cohorts; the seed-averaged
  out-of-fold level-1 Cox scores and the stage 6a folds.
- Cells: modality pair x horizon. Pairs EHR+ECG (both cohorts), EHR+CTPA and
  ECG+CTPA (three-modality cohort); horizons 10, 14, 21, 30, 45, 60, 90,
  120, 180, 270 and 365 days.
- Models, fitted in the training folds of every seed and averaged over
  seeds: at each horizon, a logistic stack on the winsorised scores of each
  single modality and each pair (main); one Cox stack per model with
  Breslow risks at every horizon (sensitivity).
- Per cell: gap = AUROC of the better single modality minus the other;
  fusion gain of the pair over the better single modality, by AUROC, IPA, a
  log-likelihood skill score (1 - LL / LL of the event-rate model) and the
  continuous NRI; the within-class correlation of the two scores (normal
  scores); the theoretical gain of the optimal linear combination under the
  binormal model.
- Tests: slope and correlation of gain on gap along time within each pair,
  and pooled over all cells; observed against theoretical gain (slope and
  ratio); the same regression for the IPA, log-likelihood and NRI gains.
  Incident/dynamic AUROC (Heagerty and Zheng, Biometrics 2005) for each
  single modality and pair (Cox-stack score): for every death, the share of
  patients still at risk that it outranks, kernel-smoothed on log time
  (bandwidth 0.4) at the same 11 times; the incident gap and gain, and the
  regression of one on the other.
- Intervals: 2,000-sample patient bootstrap, seed 42; the three-modality
  cohort is taken as the subset of each resample.
- Expected: negative slopes, as the ECG's weight declines with time while
  the EHR's does not; observed gains close to, and slightly below, the
  theoretical gains, since out-of-fold stacks lose a little; all four gain
  measures falling as the gap widens, and any that does not would show the
  pattern is partly specific to the AUROC scale. Adjacent horizons share
  patients and deaths, so the pooled correlations are descriptive; the
  bootstrap intervals carry the uncertainty.

## A10. 30 September 2026: fusion architecture by time-varying effects

Specified after A1 to A9 and all their results, before any analysis below
was run. Extends the registered early fusion (section 4: ridge-Cox on
block-standardised concatenated features, and a random survival forest on
the concatenation).
- Question: does the SIDE PROJECT's ranking of fusion families (late at
  least as good as early, both better than intermediate, across 14
  architectures and 12 transfer and outcome cells) hold for time-to-event
  prediction, and does letting each modality's effect change with time
  alter it?
- Design: fusion family (late, early, intermediate) x effects over time:
  proportional hazards (PH: one model, constant effects, Breslow baseline)
  or time-varying (TV: a separate model at each horizon).
- Frames and outcome: Frame A death at 10 and 30 days; Frame B death at 10,
  30, 90, 180 and 365 days; both cohorts. Modality sets: EHR+ECG; in the
  three-modality cohort also EHR+CTPA and EHR+ECG+CTPA.
- Linear arms (stage 18a): glmnet in R, ridge penalty, lambda by inner
  3-fold subject-grouped cross-validation on the deviance (partial
  likelihood for Cox, binomial for logistic); median imputation and
  standardisation inside each training fold; every stage 4 seed and fold,
  predictions averaged over seeds.
  - EHR alone, PH (ridge-Cox) and TV (ridge logistic per horizon).
  - Early PH (registered): ridge-Cox on the block-standardised
    concatenation, each block's standardised columns divided by the square
    root of its size. Early TV: ridge logistic per horizon on the same.
  - Full set only: Early PH with one penalty factor per block (EHR 1,
    each other block 1, 4, 16 or 64, the combination chosen by inner
    cross-validation; Boulesteix et al. 2017), and Early PH with the EHR
    block unpenalised (Herrmann et al. 2021).
  - Late PH and Late TV, nested: per-modality ridge-Cox scores out of fold
    within the training fold (inner 3-fold), winsorised as in D3, then a
    Cox stack (PH) or a logistic stack per horizon (TV). These nested arms,
    not the stage 8 fusion, are the late comparators, so every arm sees the
    same training data.
- Forest arms (the registered random survival forest on the concatenation,
  and a survival block forest; Hornung and Wright 2019) and intermediate
  arms (neural networks with one encoder per modality and a fusion layer,
  with a Cox output for PH or a discrete-time hazard output for TV, with
  multimodal dropout; covering the registered DeepSurv, DeepHit and MLP
  twin) will be specified in a later entry before they are run.
- Evaluation (stage 18c), every arm on the same prediction type (risk at
  each horizon): AUROC, IPA, log-likelihood skill score, observed/expected,
  calibration slope; 2,000-sample bootstrap, patient-clustered in Frame A,
  seed 42. Contrasts: each arm minus the matching late arm; TV minus PH
  within each family; block-aware minus registered early; each fused arm
  minus EHR alone; dense versus sparse, (Early - Late) with EHR+ECG minus
  (Early - Late) with EHR+CTPA, in the three-modality cohort.
- Expected, from the SIDE PROJECT and the DISSERTATION: late at least as
  good as the registered early arm at every horizon; the early deficit
  larger with the dense ECG block than with the sparse CTPA block
  (DISSERTATION script 141); block-aware early fusion closing much of that
  deficit; TV minus PH largest at 10 and 30 days, where the ECG's weight is
  front-loaded.

## A11. 30 September 2026: horizon-tuned proportional-hazards models

Specified after A10's linear results (stages 18a to 18c, with D9 and D10)
had been seen, and before anything below was run.
- Motivation: in A10 every PH model was tuned once, on the whole follow-up.
  In Frame B that choice shrank the ECG block strongly (64 times the EHR in
  80% of folds), and the block-penalty model lost to plain concatenation at
  10 days but beat it at 180 and 365 days. The TV models, fitted separately
  at each horizon, lost at 10 days. That loss mixes two things: targeting
  the horizon, and learning only from the deaths before it.
- Models (stage 18d): EHR alone; Early PH (block-standardised, each set);
  Early PH block penalties (full set, the A10 grid with 256 added as in
  D9); Late PH (nested, each set: every modality model tuned at the
  horizon, then a Cox stack per horizon on their inner out-of-fold scores,
  winsorised as in D3). Each is a ridge-Cox model fitted on every training
  death, as its A10 counterpart, but with its penalty, and for the
  block-penalty model its block weights, chosen separately at each horizon:
  the choice with the lowest inner-cross-validated binomial log-likelihood
  loss of the model's Breslow risk at that horizon, the rule that tunes the
  TV models. Same 50-value penalty path as cv.glmnet, same preprocessing,
  inner folds and seeds as stage 18a; glmnet with standardize = FALSE.
- Frames, cohorts, sets, horizons, measures and bootstrap as in A10
  (evaluated by stage 18c).
- Contrasts: horizon-tuned minus whole-follow-up PH, for each model (the
  gain from targeting the horizon); TV minus horizon-tuned, for EHR alone,
  early and late (the effect of fitting on the horizon's label instead of
  every death); horizon-tuned early minus horizon-tuned late; horizon-tuned
  block penalties minus horizon-tuned Early PH; each horizon-tuned fused
  model minus horizon-tuned EHR alone; dense versus sparse on the
  horizon-tuned models. Also reported: the block weights chosen at each
  horizon, and each horizon's penalty against the whole-follow-up choice.
- Expected: in Frame B the ECG block is shrunk less at 10 and 30 days than
  at 180 and 365 days; horizon tuning helps most at 10 and 30 days and
  little at 365 days, where the whole-follow-up tuning is already dominated
  by later deaths; TV minus horizon-tuned is negative at 10 days, because
  the horizon-tuned model learns from every death, and near zero at 365
  days; late remains at least as good as early.

## A12. 1 October 2026: zero-shot survival transfer with forests, along
## follow-up time

Specified after S6's first run (stage 16b), S4 (stage 11c) and the A10 and
A11 reruns had been seen, and before anything below was run. Builds on S6
as rerun under D12.
- Questions: (1) Does the forest advantage that S4 found inside MIMIC for
  the EHR survive zero-shot transfer between INSPECT and MIMIC, or do
  forests transfer worse than ridge-Cox? (2) How does zero-shot transfer
  performance change along follow-up time, and does the time pattern of the
  EHR and CTPA modalities' discrimination seen in MIMIC also appear in
  INSPECT?
- Data: Frames I and M; the 28 EHR and 38 CTPA features shared by both;
  source medians and source-percentile winsorising (S6's main strategy
  under D12); the stage 4b folds and their six seeds, predictions averaged
  over seeds.
- Learners, per modality (EHR, CTPA): ridge-Cox exactly as S6 under D12; a
  random survival forest (randomForestSRC, log-rank splitting, 500 trees,
  nsplit 10) tuned on out-of-bag Harrell's C over node sizes 5, 15, 40, 80
  and 160 and mtry ceiling(sqrt(p) / 2), ceiling(sqrt(p)) and
  ceiling(3 sqrt(p)) (the D11 grid); and the rank average of the two scores
  (ranking only).
- Fused (EHR + CTPA), for each of the three: a Cox stack on rank-based
  normal scores of the modality scores, the transform fitted on the source
  and the same for every family (as S4), with stack inputs out of sample
  within the training data (inner cross-validation for ridge-Cox, as in S6;
  out-of-bag predictions for forests); Breslow risks.
- Designs: internal (5-fold cross-validation within each dataset) and
  transfer (trained on all of one dataset, applied unchanged to all of the
  other), both directions. No target outcome is used in transfer.
- Horizons: 10, 30, 90, 180 and 365 days.
- Measures in the target: IPCW time-dependent AUC (Uno et al. 2007), IPCW
  Brier score and IPA (Kaplan-Meier censoring model in INSPECT; complete
  follow-up in MIMIC), Kaplan-Meier observed/expected, calibration slope
  (point estimate); for the single modalities, the AUC of each modality
  score at each horizon, internal and transferred.
- Contrasts: forest minus ridge-Cox and rank average minus ridge-Cox,
  internal and transferred, in each direction and at each horizon; the
  transfer loss (transferred minus internal) of each learner; forest minus
  ridge-Cox in transfer loss (a difference of differences); the EHR minus
  CTPA gap in AUC along the horizons, in each dataset.
- Intervals: 2,000-sample bootstrap of the target patients, seed 42, with
  the censoring weights re-estimated in every resample and the predictions
  held fixed. S6 used influence-function intervals; the bootstrap is used
  here because the difference of differences needs it.
- Expected: (1) forests beat ridge-Cox for the EHR internally in both
  datasets but lose more in transfer, because the datasets are almost
  perfectly distinguishable (membership AUROC 0.98) and trees cannot
  extrapolate beyond the source data's range; the forest-minus-ridge
  difference in transfer loss is negative, and the forest's advantage after
  transfer is small or absent. (2) Transfer loss in AUC is largest at 10
  and 30 days, where the two clocks differ most (the CTPA scan in INSPECT,
  admission in MIMIC), and smallest at 365 days; observed/expected is off
  in the same direction at every horizon. (3) In both datasets the CTPA's
  discrimination falls with horizon faster than the EHR's, as the ECG's
  does in MIMIC, because it carries acute-phase information.

## A13. 1 October 2026: a late-fusion Cox stack with time-varying
## coefficients

Specified after A1, the A10 and A11 reruns (D11), S4 and the S6 rerun (D12)
had been seen, and before anything below was run. The registered primary
analysis, with its step at day 30, is unchanged.
- Question: can one Cox stack whose modality weights change with follow-up
  time combine the strengths of the PH stack (one model fitted on every
  death) and of the per-horizon stacks (weights that follow the horizon,
  and better calibration at long horizons), while giving each patient one
  coherent survival curve?
- Data and inputs (stage 18e): Frames A and B, both cohorts and the A10
  modality sets; the nested late-fusion modality scores exactly as stage
  18a computes them (D9, D10, D11), with 18a's Late PH and Late TV rebuilt
  from them as a check.
- Model: a Cox stack on the winsorised modality scores Z with hazard
  h0(t) exp(Z b + Z c g(t)): each score has a constant coefficient b and a
  coefficient c on a time function g shared by all scores. Candidate g:
  none (the PH stack); a step at day 7; a step at day 30 (Frame B only);
  log time, constant within intervals with boundaries at 1, 2, 4, 7, 14,
  30, 60, 120 and 240 days (Frame A: up to 14) and valued at each
  interval's geometric midpoint. Fitted with coxph on data split at those
  boundaries (Efron ties); risks at each horizon from the Breslow baseline
  with the time-varying linear predictor.
- Arms: Late TVC, with g chosen inside each training fold by the
  cross-validated partial log-likelihood (Verweij and van Houwelingen 1993)
  over 18a's inner folds; Late TVC log-time, with g fixed at log time (a
  sensitivity analysis); EHR alone TVC, the chosen-g model on the EHR score
  alone (for the fusion gain).
- Horizons, measures and bootstrap as in A10 (stage 18c).
- Contrasts: Late TVC minus Late PH, minus Late TV and minus Late PH
  horizon-tuned; Late TVC log-time minus Late TVC; Late TVC minus EHR alone
  TVC (fusion gain). Also reported: the share of folds choosing each g, and
  the log-time model's b and c for each modality.
- Expected: a time function chosen in most Frame B folds, most often log
  time (A1: a decline linear in log time fitted best); the ECG's c
  negative (its weight fades) and the EHR's smaller in magnitude; Late TVC
  at least as good as Late PH on AUROC, with log-likelihood and calibration
  slope better than Late PH at 10 and 30 days and close to Late TV at 90 to
  365 days.
  ## A14. 1 October 2026: a survival super learner as the second layer

Specified after A10 to A13 (with D11), S4 and A12 had been seen, and before
anything below was run.
- Question: does a super learner over late-fusion stacks and learner
  families beat each single stack, in particular S4's ridge-Cox and forest
  Cox stacks?
- Inputs (stage 18f): Frames A and B, both cohorts, the A10 modality sets;
  the stage 6a ridge-Cox and stage 6b random survival forest (RSF)
  out-of-fold death scores, per seed and outer fold, as rank-based normal
  scores on the rows each stack is fitted to (as S4).
- Candidates: Kaplan-Meier (no covariates); EHR alone, ridge-Cox and RSF (a
  Cox stack on the one score); and for each set, the Cox stack (PH), the
  per-horizon logistic stacks (TV) and A13's log-time-varying Cox stack
  (TVC), each on ridge-Cox scores and on RSF scores, plus a Cox stack on
  both families' scores. Ten candidates per set.
- Super learner (Polley and van der Laan 2011; Westling et al. 2023):
  inside each outer training fold, every candidate is fitted on two of three
  inner subject-grouped folds and predicts the third. Weights, non-negative
  and summing to one, minimise the squared error of those predicted risks
  over the frame's horizons; follow-up is complete, so no censoring weights
  are needed. The candidates are then refitted on the whole training fold
  and combined. Arms: Super learner (one weight vector for all horizons);
  Super learner, per-horizon weights (one vector per horizon, the
  time-dependent analogue of TD-RMSurv; Flack et al. 2025); Super learner,
  best single (the candidate with the lowest inner-cross-validated loss).
  Six seeds, predictions averaged.
- Measures: AUROC, IPA, log-likelihood skill score, observed/expected,
  calibration slope (point estimate); 2,000-sample bootstrap, seed 42,
  patient-clustered in Frame A.
- Contrasts: Super learner minus the ridge-Cox, RSF and both-families Cox
  stacks, the best single, EHR alone (ridge-Cox), and stage 18a's Late PH
  and Late TV; per-horizon minus common weights. The comparisons with Late
  PH and Late TV mix two layer-1 pipelines (stage 6a/6b scores here, nested
  ridge-Cox scores in 18a) and are reported as such. Also reported: the
  mean weights and how often each candidate is the best single.
- Checks: the ridge-Cox and RSF Cox stacks on each cohort's full set must
  reproduce S4's fused AUROCs.
- Expected: the super learner at least as good as every candidate on IPA
  and log-likelihood (its oracle property), with most weight on RSF-based
  candidates; its AUROC gain over the RSF Cox stack below 0.01;
  per-horizon weights no better than common ones at these event counts.

## A15. 1 October 2026: the CTPA's time pattern like for like, and its
## level-1 tuning

Specified after S5, stage 8, A12 and A13 had been seen, and before anything
below was run. Motivation: in MIMIC the CTPA's hazard ratio was
front-loaded (S5: early/late ratio 1.194 [1.026, 1.450]) and its stack
weight faded (stage 8, A13), whereas in A12 its cumulative AUROC rose with
horizon in INSPECT. Those analyses differed in measure, clock, dataset and
feature set (the 46 MIMIC features, with the device block, against the 38
shared with INSPECT).
- Scores (stage 16e): ridge-Cox death scores fitted exactly as stage 18a's
  EHR alone (glmnet, D11 path, deviance tuning, winsorised features), out
  of fold over the stage 4 / 4b folds and every seed, seed-averaged. Frame
  B (three-modality cohort, discharge clock): CTPA 46, CTPA 38, EHR. Frame
  M (admission clock): CTPA 46, CTPA 38, EHR. Frame I (INSPECT, scan
  clock): CTPA 38, EHR. Frame A (three-modality cohort): CTPA 46, for part
  2 only.
- Part 1, time pattern (Frames B, M, I; scores winsorised and
  standardised): per-SD hazard ratios for death in days 0-30 and 31-365
  and their ratio, as in S5, for the CTPA alone and with the EHR score;
  incident/dynamic AUROC as in A9 (kernel-smoothed on log time, bandwidth
  0.4) at 10, 30, 90, 180 and 365 days for the CTPA and EHR scores, their
  gap, and the change from 10 to 365 days; cumulative IPCW AUROC at the
  same horizons, as in A12, for reconciliation (point estimates).
  Contrasts of the early/late ratio: 46 over 38 features in B and in M
  (paired: the device block); M over I with 38 features (independent: the
  hospital and clock).
- Part 2, CTPA tuning (Frames A and B, CTPA 46): these ridge-Cox risks
  against stage 6a's CTPA ridge-Cox (tuned on Harrell's C over five
  penalties): AUROC, IPA, log-likelihood skill score, observed/expected and
  calibration slope at every horizon, with intervals for the differences;
  the share of penalty choices at either end of the path.
- Intervals: 2,000-sample bootstrap, seed 42, scores held fixed,
  patient-clustered where patients have several admissions.
- Checks: Frame B's EHR scores must reproduce stage 18a's EHR alone PH
  risks.
- Expected: (1) a front-loaded CTPA hazard ratio (above 1) in Frames B and
  M with the 46 features, smaller with the 38 shared features, because the
  device block carries in-hospital deterioration; in INSPECT a ratio closer
  to 1. (2) The CTPA's incident AUROC falls with time wherever its ratio is
  above 1, while its cumulative AUROC need not, reconciling S5 and A12.
  (3) Retuning leaves the CTPA's AUROC almost unchanged but moves its
  calibration slope from well above 1 (S4: about 2.6) towards 1.
  ## A16. 2 October 2026: retuned ridge-Cox, oblique forests and CoxBoost
## at layer 1, and a larger super learner

Specified after A14, A15, A12 and the stage 6c sensitivity analysis of S4
had been seen, and before anything below was run. Motivated by A15
(retuning the CTPA's stage 6a ridge-Cox with the D11 path and deviance
tuning fixed its calibration) and by a neutral benchmark of survival
learners on low-dimensional data (Burk et al., Bioinformatics 2026), in
which oblique random survival forests and likelihood-based boosting had
the best average ranks.
- Level 1 (stage 6d): death, every frame, cohort, modality and seed, out
  of fold over the stage 4 folds, with stage 18a's preprocessing (median
  imputation, winsorising, standardisation) and inner folds.
  - Ridge-Cox retuned: exactly stage 18a's EHR alone (glmnet, D11 path,
    deviance tuning), for every modality; its EHR risks must reproduce
    stage 18a's EHR alone PH.
  - Oblique random survival forest (aorsf; Jaeger et al. 2019, 2024): 500
    trees, Cox-based linear combinations, tuned on out-of-bag Harrell's C
    over minimum leaf size 5, 20, 60 x mtry ceiling(sqrt(p) / 2),
    ceiling(sqrt(p)), ceiling(2 sqrt(p)); score = mortality.
  - CoxBoost (component-wise likelihood-based boosting; Binder and
    Schumacher 2008): default penalty (9 x events), number of steps (0 to
    300) by the inner cross-validated partial log-likelihood; score = the
    linear predictor.
  The stage is blinded: it writes completeness, the D2 per-fold rule, the
  reproduction check and where the tuning choices fell.
- Layer 2 (stage 18g): A14's super learner with five families: ridge-Cox
  (stage 6a), ridge-Cox retuned, RSF (stage 6b), oblique RSF and CoxBoost.
  Candidates: Kaplan-Meier; EHR alone for each family; for each set and
  family the Cox stack (PH), the per-horizon logistic stacks (TV) and the
  log-time-varying Cox stack (TVC); A14's Cox stack on ridge-Cox and RSF
  scores together. A14's eleven candidates are a subset and must
  reproduce A14's predictions. Arms: Super learner v2 (one weight vector)
  and its best single candidate; per-horizon weights are dropped, as A14
  found they added nothing.
- Layer-1 comparison: each modality's ridge-Cox risks, retuned against
  stage 6a (A15's comparison, extended to the EHR and ECG).
- Measures and bootstrap as in A14.
- Contrasts: Super learner v2 minus A14's super learner, minus Late PH and
  minus its best single; ridge-Cox retuned PH minus ridge-Cox PH (whether
  retuning layer 1 changes the fused stack); oblique RSF minus RSF, for PH
  and TV; CoxBoost PH minus ridge-Cox retuned PH; for each modality,
  retuned minus stage 6a ridge-Cox at layer 1.
- Expected: (1) at layer 1, retuning moves calibration slopes that are
  above 1 towards 1, with AUROC changes under 0.01; (2) in the fused
  stacks, retuned ridge-Cox within 0.005 AUROC of stage 6a ridge-Cox,
  because stacks rescale their inputs; (3) oblique forests at least as good
  as standard forests for the EHR, and no better for the ECG; CoxBoost
  close to ridge-Cox; (4) Super learner v2 at least as good as A14's,
  gaining under 0.005 AUROC.
  ## A17. 2 October 2026: the registered analyses rerun on retuned
## ridge-Cox scores

Specified after A15 and A16 showed that stage 6a's ridge-Cox (tuned on
Harrell's C over five penalties) was sub-optimal, and before any
registered script was rerun on retuned scores.
- Level 1 (stage 6e): every stage 6a ridge-Cox target, modality, frame,
  cohort and seed refitted exactly as stage 18a's EHR alone (D11 path,
  deviance tuning, winsorised features). The binary-family (L2 logistic)
  scores are unchanged. The death scores must reproduce stage 6d exactly.
- The registered scripts (stages 7, 7b, 7c, 8, 9, 11a and 11b) are rerun
  unchanged in a working copy of the project whose level-1 ridge-Cox file
  is stage 6e's. Registered outputs are not touched. Every registered
  estimate is reported beside its retuned counterpart (stage 13r).
- Status: a sensitivity analysis. The registered results remain the
  registered results.
- Expected: the primary decision unchanged (the ECG's early/late hazard
  ratio above 1), with estimates moving by less than their interval
  widths; better calibration in stage 11b.

## A18. 2 October 2026: Analysis B, a stronger EHR model

Proposed on 30 September; specified here before anything was run.
- Question: does the ECG's added value, and its front-loaded weight (the
  primary estimand), survive a much stronger EHR model? This is the
  DISSERTATION's script 140 test in survival form.
- Frame B only: the extra labs run to discharge, so in Frame A they could
  leak in-hospital deaths. Both cohorts.
- EHR versions, all ridge-Cox exactly as stage 18a (D11): EHR 28 (stage
  6d's scores); EHR 79, the 28 plus the 51 extra index-admission labs of
  ehr_extra.csv (vancomycin level and blood-gas temperature excluded, as
  recorded at the pull); EHR 79 + MI, adding a missing-value indicator for
  every extra lab with missing values (EHR+ECG cohort only). ECG and CTPA
  scores from stage 6d.
- Fused models, per seed and fold on rank-based normal scores: Cox stack
  (PH) and per-horizon logistic stacks (TV), for EHR alone, EHR + ECG,
  and in the three-modality cohort EHR + CTPA and EHR + ECG + CTPA.
- The primary's estimand (the ECG's hazard-ratio ratio, days 0-30 over
  31-365) with EHR 28 and with EHR 79 as the adjusting score, and their
  ratio.
- Contrasts (AUROC, IPA, log-likelihood, 10-365 days; 2,000-sample
  bootstrap): EHR 79 minus EHR 28, alone and fused; EHR 79 + MI minus EHR
  79; the fusion gain of each modality set with each EHR version, and the
  difference of those gains.
- Expected: EHR 79 beats EHR 28 alone by 0.01-0.03 AUROC; the ECG's
  fusion gain shrinks with EHR 79 but stays positive at 10-30 days; the
  ECG's early/late ratio stays above 1; missing indicators add little.

## A19. 2 October 2026: a penalised-spline Cox model for the EHR

Specified after S4, the stage 6c sensitivity analysis and A16 showed a
consistent forest advantage for the EHR, and before anything below was
run.
- Question: is that advantage smooth non-linearity that an additive Cox
  model can capture?
- Level 1 (stage 19b): EHR death scores, every frame, cohort and seed,
  out of fold, with stage 18a's preprocessing. Additive Cox model
  (mgcv::gam, family cox.ph): features with at least 10 distinct values
  as cubic regression splines (k = 5), the rest linear; smoothness by
  REML, select = TRUE.
- Compared, through the same Cox (PH) and per-horizon logistic (TV)
  stacks, with the EHR scores of ridge-Cox retuned, RSF and oblique RSF,
  for EHR alone, EHR + ECG and (three-modality cohort) EHR + ECG + CTPA.
  Also the spline model's own risks against ridge-Cox retuned's.
- Contrasts: spline minus ridge-Cox retuned; oblique RSF minus spline;
  RSF minus spline. 2,000-sample bootstrap, patient-clustered in Frame A.
- Expected: the spline model closes part of the gap between ridge-Cox and
  the oblique forest for the EHR, and calibrates better than the forests;
  in fused stacks, differences under 0.005 AUROC.
  ## A20. 2 October 2026: binary models retuned like the ridge-Cox, and the
## Cox-versus-binary analyses rerun

Specified after A17 showed that S1's registered binary-over-Cox advantage
disappeared when only the ridge-Cox was retuned, and before anything below
was run.
- Level 1 (stage 6f): every stage 6a binary model (L2 logistic: death,
  composite, and CV readmission with and without deaths in training, at
  every horizon), modality, frame, cohort and seed refitted with stage
  6e's approach: winsorised features; glmnet, alpha 0, 80 penalties down
  to 1e-7 of the largest; lambda by inner 3-fold cross-validated log-loss,
  with inner folds grouped by patient and stratified by the label. A label
  with fewer than 3 events in a training fold gets that fold's event rate
  (counted in the checks). The ridge-Cox columns are stage 6e's.
- The registered scripts for stages 8, 9, 11a and 11b are rerun unchanged
  in a second working copy (SurvivalProject_retuned2) whose level-1 file
  is stage 6f's, and compared with the registered results and with A17's
  copy (the latter isolates what retuning the binary models changed).
- Status: a sensitivity analysis.
- Expected: the binary models gain calibration and a little
  discrimination; at 30 days (S1) Cox and binary end within 0.005 AUROC
  of each other; at 90 days and beyond (S3) Cox keeps an advantage,
  because one survival model borrows strength across horizons.

## A21. 2 October 2026: spline Cox as a sixth super-learner family
## (Super learner v3)

Specified after A19 showed that the penalised-spline Cox model matched the
random survival forest for the EHR, and before anything below was run.
- Families: A16's five plus spline Cox, whose EHR scores are stage 19b's
  (per seed) and whose ECG and CTPA scores are ridge-Cox retuned, as in
  stage 19b's fused models.
- Fitted in stage 18h, which also refits A16's library (Super learner v2,
  which must reproduce A16 exactly) and A22's (Super learner v4) from the
  same candidate predictions. Bootstrap as in A16.
- Contrasts: v3 minus v2; spline Cox minus ridge-Cox retuned (PH and TV
  stacks); oblique RSF minus spline Cox (PH and TV); and A16's contrasts
  again.
- Expected: v3 within 0.003 AUROC of v2, with spline-based stacks taking
  weight mainly from the ridge-Cox stacks.

## A22. 2 October 2026: four tuning changes from the literature, and
## Super learner v4

Specified after the stage 6c, A16 and A19 results and a literature check
on tuning survival models, and before anything below was run.
1. Proper-score tuning of the forests (stage 6g): the random survival
   forest (stage 6c's wide grid) and the oblique random survival forest
   also chosen by out-of-bag Brier score summed over the frame's horizons,
   not only by out-of-bag Harrell's C. Node size mainly affects
   calibration, which C cannot see. The C choices must reproduce stages
   6c and 6d exactly.
2. Elastic-net Cox (stage 6h): alpha 0, 0.25, 0.5, 0.75 or 1 chosen with
   lambda on inner cross-validated deviance, for every modality; alpha 0
   must reproduce stage 6d's ridge-Cox exactly.
3. Per-horizon recalibration of the super learner (stage 18h): a logistic
   intercept and slope on the logit of the super learner's risk, fitted on
   its inner out-of-fold predictions in each training fold and applied to
   the test fold. Recalibration is monotone, so it cannot change AUROC.
4. Extended grids: the oblique forest's mtry up to ceiling(4 sqrt(p))
   (stage 6g); CoxBoost refitted with up to 500 steps in the folds where
   stage 6d's choice hit 300 (stage 6h). The other folds keep stage 6d's
   fit, on the assumption that an optimum found below the cap is global
   (the cross-validated curve up to 300 steps is identical).
- Super learner v4: A21's library plus five families (RSF Brier-tuned,
  oblique RSF wide grid, oblique RSF Brier-tuned, elastic-net Cox,
  CoxBoost 500), with its best single candidate and a recalibrated
  version; v3 is also recalibrated.
- Contrasts: v4 minus v3, v2, Late PH, Late TV and its best single;
  recalibrated minus not (v3, v4); recalibrated v4 minus Late PH; each new
  family's stack minus its predecessor; and, on each modality's own
  level-1 risks, elastic net minus ridge-Cox retuned, Brier-tuned minus
  C-tuned (both forests), and oblique forest wide minus original grid.
- Expected: (1) Brier tuning chooses larger node sizes and improves the
  forests' own calibration (IPA, log-likelihood, slopes nearer 1), with
  AUROC changes under 0.005; (2) elastic net close to ridge-Cox for the
  EHR and CTPA, with at most a small gain for the ECG; (3) recalibration
  brings super-learner slopes close to 1 and improves log-likelihood
  slightly; (4) little or no gain from the extended grids; Super learner
  v4 within 0.003 AUROC of v3.
  ## A23. 3 October 2026: interpretability of the final model's components

Specified after A21 and A22 made the penalised-spline Cox model the largest
super-learner family, and before anything below was run. Descriptive only;
no hypothesis tests.
- Spline Cox partial effects (stage 20a): the stage 19b model refitted once
  on every patient of each frame and cohort, for interpretation rather than
  prediction. For each EHR feature: effective degrees of freedom, Wald
  chi-square and p-value; the hazard ratio across its 1st to 99th
  percentile relative to its median, with 95% intervals from the model's
  covariance; and the SD of its term across patients as an importance score.
- SurvSHAP(t) (Krzyzinski et al. 2023) for the Brier-tuned RSF and oblique
  RSF on the EHR alone, EHR+ECG cohort, Frames B and A, each refitted once
  at the settings stage 6g chose most often. Interventional Shapley values
  of predicted risk at every horizon, estimated by 32 antithetic
  permutations against 100 background patients, for 300 randomly chosen
  patients. Each patient's values must add up exactly to their risk minus
  the background mean.
- Summaries: mean absolute SurvSHAP(t) per feature and horizon; each
  feature's share of the total at each horizon; Spearman agreement between
  forest importance at the last horizon and the spline term SD.
- Figures (stage 20b), Frame B: the 8 strongest partial effects; importance
  share over time for each forest; importance agreement.
- Expected: the forests and the spline model agree on the strongest
  features (Spearman at least 0.7), led by age, cancer and acute
  physiology; acute physiology loses share of importance between 10 and
  365 days while chronic conditions gain, mirroring the front-loaded ECG
  signal.

## A24. 3 October 2026: temporal validation of the final pipeline

Specified before anything below was run. This is not the registered check
for truncation at the end of MIMIC's collection window, which remains in
the sensitivity list.
- Frame B. Each index admission's approximate real year is the midpoint of
  its MIMIC anchor window (approx_year_hi minus 1). Development: 2013 or
  earlier. Test: 2014 onwards, pooled and by period (2014-2016; 2017
  onwards).
- Core library: ridge-Cox retuned, spline Cox (EHR; ridge-Cox for ECG and
  CTPA), RSF Brier-tuned and oblique RSF Brier-tuned, which together carry
  nearly all of Super learner v4's weight. Every tuning step is replayed
  inside the development data: level-1 scores out of fold (5 folds) for
  development patients and from refits on all development patients for
  the test patients; then stage 18h's candidates, super learner weights
  (inner 3-fold) and per-horizon recalibration, all fitted on development
  only and applied unchanged to the test years.
- Evaluation at 10-365 days: AUROC, IPA, log-likelihood skill score,
  observed/expected, calibration slope, in the test years (2,000-sample
  bootstrap) and, for reference, cross-validated within the development
  years. Contrasts in the test years: recalibrated super learner minus the
  late-fusion analogue (ridge-Cox retuned Cox stack), minus EHR alone
  (spline Cox), and minus the uncalibrated super learner. Drift table:
  patients, deaths and Kaplan-Meier one-year risk per period.
- Expected: test-year AUROC within 0.02 of the development estimate;
  observed/expected drifting further from 1 than discrimination does; the
  super learner's advantage over the late-fusion analogue persisting at
  most horizons.
  ## A25. 3 October 2026: the remaining registered sensitivity analyses (stage 12)

Specified before running; the pre-registration (section 8) named these
analyses without the details below. Hospice discharges (stage 7b), sPESI
(stage 16c) and the INSPECT one-year mortality bounds (stages 4b and 16b)
are already done and not repeated.
- Final data year (Frame B CV readmission and composite): patients whose
  index admission may fall in MIMIC's final data year (stage 4's
  final_year_possible) are excluded and the registered out-of-fold
  predictions re-scored, models unchanged: S2's four formulations (AUROC
  in the full cohort and among patients not dying first, Brier,
  observed/expected; F2, F3 and F4 minus F1, F4 minus F3) and stage 11b's
  composite models (Cox stack, binary stack, recalibrated late fusion:
  AUROC, IPA, observed/expected), beside the full cohort.
- All admissions (Frame B, the primary): Frame B rebuilt with every
  admission of the EHR+ECG cohort under stage 4's exclusions and time
  coding (its first admissions must reproduce frame_B.csv exactly). EHR and
  ECG death scores from ridge-Cox retuned (D11), out of fold over 5 folds
  grouped by patient, the stage 4 seeds. The primary estimand with a
  2,000-sample patient-clustered bootstrap, beside the first-admission
  frame scored by the same learner (stage 6d), so only the cohort differs.
- 48-hour landmark (Frame A): admissions with no death and no CV
  readmission on calendar days 0 or 1 (event times above 2). S1's scores
  re-scored for outcomes by day 30 (AUROC at 30.5 days, Uno's C to 31
  days) with S1's contrasts, beside the full Frame A cohort;
  patient-clustered bootstrap.
- Expected: no conclusion changes. The ECG early/late ratio stays above 1
  with its interval excluding 1 in all admissions; S2's ordering of
  formulations is unchanged without the final data year; at the landmark
  S1's binary-minus-Cox contrasts move by less than 0.01 AUROC, while the
  fusion gain shrinks, consistent with the ECG's front-loaded signal.

## A26. 3 October 2026: the registered neural learners (stage 10) and A10's forest and intermediate arms (stages 18i and 18j)

Specified before running, completing S4 and the arms A10 deferred.
- Stage 10 (S4's neural family), per frame, cohort, modality and seed, out
  of fold over the stage 4 folds: DeepSurv (Cox partial likelihood) for
  death, composite and cause-specific CV readmission; DeepHit (competing
  risks, CV readmission first and death first; 20 intervals including
  every horizon; alpha 0.5, sigma 0.1); the MLP binary twin (one output per
  horizon, binary cross-entropy) for death, composite and unconditional CV
  readmission. Two hidden layers (64 and 32 units), Adam (1e-3), batches of
  256, weight decay 1e-4. Matched tuning budget: two candidates per learner
  (dropout 0.1 or 0.5), chosen on a 20% validation split of the training
  fold (grouped by patient) with early stopping (patience 15, at most 300
  epochs), then refitted on the whole training fold for the best number of
  epochs. Stage 18a's preprocessing.
- Stage 11e: DeepSurv minus ridge-Cox (stage 6a, and retuned for death),
  MLP twin minus L2 logistic, DeepHit minus cause-specific ridge-Cox, and
  fused DeepSurv minus fused ridge-Cox (Cox stack on normal scores);
  2,000-sample bootstrap, patient-clustered in Frame A.
- Stage 18i, forest arms on the concatenated modalities, every A10 frame,
  cohort and set, every seed and fold: Early RSF (the registered random
  survival forest on the concatenation; node size 10 or 40 by mtry
  ceiling(sqrt(p)) or ceiling(2 sqrt(p)), chosen on out-of-bag Brier) and a
  survival block forest (Hornung and Wright 2019; each feature offered with
  probability proportional to its modality's weight divided by the
  modality's size; EHR weight 1, other weights from 0.25, 0.5 or 1 chosen
  on out-of-bag Brier).
- Stage 18j, intermediate arms: one encoder per modality (EHR 16 units, ECG
  32, CTPA 16), a fusion layer (32 units), multimodal dropout 0.25. Heads:
  Intermediate PH (Cox partial likelihood, Breslow risks) and Intermediate
  TV (discrete-time logistic hazard over the intervals between horizons).
  Two candidates per head (dropout 0.1 or 0.4), tuned as in stage 10.
- Stage 18k: each new arm minus Late PH, minus Late TV and minus EHR alone
  PH; Block forest minus Early RSF; Early RSF minus Early PH (plain
  concatenation); Intermediate TV minus Intermediate PH; 2,000-sample
  bootstrap, patient-clustered in Frame A.
- Expected: neural learners no better than ridge-Cox (within 0.01 AUROC or
  below), the MLP twin no better than L2 logistic, DeepHit no better than
  cause-specific ridge-Cox and worse calibrated. No new arm beats Late PH
  at most horizons, as intermediate fusion never beat late fusion in the
  side project; the block forest at least matches Early RSF where the ECG
  (the dense block) is included.
  ## A27. 3 October 2026: the block forest, optimised, and a super learner that can use it

Specified after stage 18k, before anything below was run. Exploratory:
A10's matched-budget comparison (stage 18k, first version) remains the
registered result.
- Why: in stage 18i the block forest was the only alternative architecture
  to match late fusion, but its weight choices sat at the lowest value of
  a small grid in 53-63% of folds, and its node size and mtry were fixed.
- Stage 18l: the stage 18i block forest tuned in two steps on every
  training fold, by out-of-bag Brier summed over the frame's horizons:
  (1) the other modality's weight in {0.05, 0.1, 0.25, 0.5, 1}, or (ECG,
  CTPA) in {0.05, 0.1, 0.25} x {0.1, 0.25, 0.5}, at node size 10 and mtry
  ceiling(sqrt(p)); (2) at those weights, node size {5, 10, 20} x mtry
  ceiling(f sqrt(p)), f in {1, 2, 3}. Every frame, cohort and set of A10,
  every stage 4 seed and fold; out-of-fold risks saved per seed.
- Stage 18k, second version: adds "Block forest, optimised", with the
  contrasts against Late PH, Late TV and EHR alone PH, and against stage
  18i's block forest.
- Stage 18m, Super learner v5: Super learner v4's library plus the
  optimised block forest as a candidate in every set (its own risk at
  every horizon; for the inner weights, its out-of-fold risks on the
  training folds), fitted per seed and outer fold and recalibrated per
  horizon exactly as v4. It must reproduce v4 exactly before the block
  forest is added. Contrasts: v5 minus v4 (both recalibrated), v5 minus
  Late PH; the block forest's weight.
- Expected: the optimised block forest gains at most 0.005 AUROC over
  stage 18i's and stays within 0.005 of late fusion; Super learner v5 is
  within 0.003 AUROC of v4, with the block forest taking a weight below
  0.2, since v4 already contains stacks of per-modality forests.
  ## A28. 3 October 2026: the 79-feature EHR in the super learner, per-horizon weights, and the CTPA's content over time

Specified before anything below was run. The subgroup evaluation by age,
sex and ethnicity follows once the final models are fixed by the rule in
item 3.

1. Level-1 models on the 79-feature EHR (stage 6i). Frame B only (the
   extra labs run to discharge and could leak in-hospital deaths in Frame
   A), both cohorts, the stage 4 folds and six seeds; EHR 79 built exactly
   as stage 19a. Families: ridge-Cox retuned (D11 path, extended 100-fold
   if the smallest penalty is chosen, at most twice); spline Cox (stage
   19b's model with mgcv's optimizer = "efs", which with 79 features is
   about 40 times faster and gave the same fit in simulation; D14's
   fallback rule); RSF Brier-tuned (starting grid node {5, 15, 40, 80,
   160} x mtry f sqrt(p), f in {0.5, 1, 3}); oblique RSF Brier-tuned
   (leaf {5, 20, 60} x mtry f sqrt(p), f in {0.5, 1, 2, 4}). Forest grids
   are extended by tune_expand.R inside each training fold: while the best
   value of a hyperparameter is at an edge of the values tried, the next
   value beyond it is tried with the others at their best (node size and
   leaf halved down to 1 or doubled; mtry halved down to 1 or x1.75 up to
   p), kept only if it improves the out-of-bag Brier score by more than
   0.1%, at most three extensions per direction.
2. Super learner v6 and per-horizon weights (stage 18n, Frame B). Libraries
   v4, v5 (v4 + optimised block forest) and v6 (v5 + four EHR 79
   families: each family above, with the EHR score from stage 6i and the
   ECG and CTPA scores of the same family as in v4); candidates built as
   stage 18h. Each library fitted with pooled weights (one set across
   horizons, as v4) and with per-horizon weights (fitted separately at
   each horizon), every version recalibrated per horizon as v4. v4 and v5
   with pooled weights must reproduce stage 18m exactly. Contrasts: v6
   minus v5 and minus v4; per-horizon minus pooled for each library; v6
   minus Late PH; each EHR 79 candidate minus its EHR 28 counterpart (EHR
   alone and PH stacks). 2,000-sample bootstrap.
3. Choosing the final models before the subgroup evaluation: in the
   primary cohort (Frame B, EHR+ECG; recalibrated versions), v6 replaces
   v4 only if it is significantly better at two or more horizons and
   significantly worse at none; per-horizon weights replace pooled ones
   under the same rule. Otherwise Super learner v4 with pooled weights
   remains the final model. Frame A's final model stays v4.
4. The CTPA's content over time (stage 16f). Frame B, three-modality
   cohort. The 46 CTPA features in five blocks (PE descriptors 13,
   incidental thoracic findings 7, comorbidity mentions 16, support devices
   8, report metadata 2); for each block and for all 46, a ridge-Cox death
   score as stage 16e, out of fold, seed-averaged. Per block: hazard
   ratios per SD for death in days 0-30 and 31-365 and their ratio, alone
   and with the EHR score; AUROC for death by day 30 and for death in days
   31-365 among day-30 survivors. Contrasts: ratio of early/late ratios of
   each block over all 46 and over the PE descriptors; landmark AUROC
   (days 31-365) of each block minus the PE descriptors. 2,000-sample
   bootstrap.
- Expected: (1-2) the EHR 79 families improve EHR alone by 0.01-0.02
  AUROC, as Analysis B did with ridge-Cox, and v6 gains 0.005-0.015 over
  v4 at 30-90 days; per-horizon weights are within 0.003 of pooled
  weights, with no gain. (4) Incidental findings and comorbidity mentions
  are late-loaded (early/late ratio below 1, and below all 46 features'),
  support devices are early-loaded (ratio above 1), the PE descriptors
  carry little signal at any time, and incidental findings have the
  highest AUROC for deaths in days 31-365.
  - Clarification before stage 18n was run: the final-model rule in item 3 uses AUROC, the measure used for every earlier choice of final model.
  ## A29. 3 October 2026: performance of the final models by age, sex and recorded race

Specified before stage 18n's result and before any subgroup data were
pulled.
- Models: Frame B, EHR+ECG cohort, the final model chosen by A28's rule
  (read automatically from stage 18n's decision file); Frame A, EHR+ECG
  cohort, Super learner v4 recalibrated. Comparators: Late PH and EHR alone
  (spline Cox stack).
- Subgroups: age at the index admission (under 65, 65-79, 80 and over);
  sex (MIMIC patients.gender); recorded race (MIMIC admissions.race) in six
  groups: White, Black, Hispanic or Latino, Asian, Other, Unknown. Race is
  not a model input and is used only for evaluation.
- Measures per subgroup and horizon: patients, deaths, AUROC, IPA,
  observed/expected and calibration slope of the final model; the final
  model minus Late PH and minus EHR alone (AUROC); each subgroup minus the
  rest of the cohort (AUROC, observed/expected). Estimates only where the
  subgroup has at least 10 deaths and 10 survivors at that horizon.
  2,000-sample bootstrap of the whole cohort, subgroups recomputed in each
  resample (patient-clustered in Frame A).
- Expected: AUROC lower within age bands than overall (age is a strong
  predictor, so a narrower age range leaves less to separate); no
  significant difference between women and men; observed/expected between
  0.8 and 1.25 in every group with at least 50 deaths; the final model's
  advantage over Late PH positive in every group with at least 50 deaths.
  No direction is predicted for the race groups; the smaller ones will
  have wide intervals.
  ## A30. 3 October 2026: information known at discharge (EHR+) and Super learner v7

Specified before any of it was run. Frame B only (every addition is known
by discharge; in Frame A the diagnosis codes would leak the outcome).
1. Stage 3d: the diagnosis codes of every Frame B index admission
   (MIMIC-IV diagnoses_icd, BigQuery).
2. Stage 6j: "EHR+" = EHR 79 + measurements pulled by stage 3b but never
   used (spo2, nlr, hr_first, sbp_first, spo2_first, sex_f, n_prior,
   n_codes) + discharge destination indicators (home health care, skilled
   nursing, rehabilitation, long-term acute care, hospice, other; home the
   reference) + comorbidity from the index admission's codes (R package
   comorbidity, Quan's ICD-9-CM and ICD-10 maps: the 17 Charlson and 31
   Elixhauser categories, each kept if its prevalence is 1%-99%, the
   Charlson score and the van Walraven Elixhauser score). The four
   families of stage 6i (ridge-Cox retuned, spline Cox, RSF Brier-tuned,
   oblique RSF Brier-tuned), tuned exactly as stage 6i.
3. Stage 18o: Super learner v7 = v6's library + the four EHR+ families
   (each with its EHR+ score and the same family's ECG and CTPA scores),
   pooled weights, recalibrated per horizon; v6 must be reproduced exactly.
   Contrasts: v7 minus v6 and minus Late PH; each EHR+ candidate minus its
   EHR 79 counterpart; the same contrasts with patients discharged to
   hospice removed from the evaluation.
4. Rule: v7 replaces v6 as the final Frame B model only if, in the primary
   cohort (Frame B, EHR+ECG; recalibrated; AUROC), it is significantly
   better at two or more horizons and significantly worse at none, AND its
   AUROC gain stays positive at those horizons with hospice discharges
   removed. If v7 is adopted, the subgroup evaluation (A29) is rerun on it.
- Expected: the comorbidity and destination blocks help most at 90-365
  days; v7 gains 0.005-0.02 AUROC over v6 there, most in patients aged 80
  and over; the gain shrinks but stays positive without hospice discharges.
  - Scope, stated before stage 6j was run: A30 builds the DISCHARGE model
  (prediction at discharge of death in the following year). Every EHR+
  feature is known by discharge, so it is legitimate for that question,
  but not for prognosis at diagnosis, which the at-diagnosis model
  (A31, specified separately) addresses with information available at
  presentation only.
  ## A31. 3 October 2026: how much the hospice destination adds to the discharge model

Specified before it was run. Descriptive; Super learner v7 remains the
final discharge model.
- Stage 6j (--no-hospice): EHR+ rebuilt without the hospice destination
  indicator (hospice discharges carry no destination indicator, like home
  discharges); the four EHR+ families refitted exactly as stage 6j's
  second version.
- Stage 18p: Super learner v7nh = v6 + these four families, fitted,
  weighted and recalibrated exactly as v7; v2-v7 must reproduce stage 18o
  exactly. Contrasts (Frame B, recalibrated): v7 minus v7nh (what the
  hospice indicator adds) and v7nh minus v6 (what EHR+ adds without it),
  on every patient and with hospice discharges removed; AUROC and IPA.
- Expected: the hospice indicator accounts for most of v7's gain at 10
  days and little at 365 days; v7nh is significantly better than v6 at
  90-365 days.
- The at-diagnosis model, previously called A31, becomes A32.
## A32. 3 October 2026: the at-diagnosis model, benchmarked against sPESI and PESI

Specified before any of it was run.
- Question: risk of death within 10, 30, 90, 180 and 365 days of
  ADMISSION (in-hospital deaths included), using only information
  available at presentation. Cohort and folds: Frame A's (3,501
  admissions; the three-modality subset for CTPA analyses).
- Presentation window: from 24 hours before to 24 hours after the
  admission time for blood tests, vital signs and ECG; from 48 hours
  before to 24 hours after for the CTPA report (the dissertation's
  window).
- EHR at presentation: age; sex; the FIRST value in the window of heart
  rate, systolic and diastolic blood pressure, respiratory rate, oxygen
  saturation and temperature (implausible values discarded); the FIRST
  value of every lab analyte of the 79-feature EHR; history from EARLIER
  admissions only (Charlson and Elixhauser categories with 1%-99%
  prevalence, both scores, number of earlier admissions and codes).
  Excluded: anything aggregated over the stay, the index admission's
  codes, discharge destination, and the 28-feature EHR's four condition
  flags (their timing is not recorded).
- ECG at presentation: the features of the FIRST ECG in the window, taken
  from the per-ECG files with the same columns as the ECG modality; for
  admissions with a single ECG they must reproduce ecg_opt.csv exactly
  (FAIL otherwise). No ECG in the window: features missing (training-fold
  medians) plus an indicator.
- CTPA at presentation: the 46 CTPA features only where the report falls
  in its window; otherwise the admission counts as without a CTPA.
- Models: EHR at presentation, the four families of stage 6j's second
  version (D16 included) with tune_expand.R; ECG and CTPA at presentation,
  ridge-Cox retuned, RSF Brier-tuned and oblique RSF Brier-tuned (spline
  Cox uses ridge-Cox for them, as in v4). Super learner (pooled weights,
  recalibrated per horizon), candidates built as stage 18h.
- Benchmarks from presentation values: sPESI (age over 80, cancer, chronic
  cardiopulmonary disease, heart rate 110 or more, systolic pressure below
  100, saturation below 90) and a 10-item PESI (age, male sex +10, cancer
  +30, heart failure +10, chronic lung disease +10, heart rate 110 or more
  +20, systolic below 100 +30, respiratory rate 30 or more +20,
  temperature below 36 +20, saturation below 90 +20; altered mental
  status is not recorded in MIMIC). Missing vital signs score 0;
  comorbidities from earlier admissions, as for the model.
- Comparisons: AUROC at 30 days (primary) and every horizon; the super
  learner against sPESI and against PESI-10; its calibration; its
  specificity at sPESI's sensitivity (bootstrap recomputing sPESI's
  operating point); the super learner against EHR at presentation alone
  (what ECG and CTPA add at diagnosis); against Frame A's Super learner
  v4 at 10 and 30 days on the same patients (what later information added
  to the earlier, mixed models). 2,000-sample bootstrap.
- Sensitivity: cancer, heart failure and chronic lung disease taken from
  earlier OR index-admission codes (history as taken at presentation), for
  the model and both scores alike, since about half of patients have no
  earlier admission.
- Expected: the at-diagnosis super learner scores 0.02-0.05 below Frame
  A's v4 at 30 days; it beats sPESI and PESI-10 by at least 0.07 at 30
  days; ECG and CTPA add more over the EHR at presentation than they did
  in the mixed and discharge models.
  ## A33. Option B: the composite and CV readmission through both final pipelines, with matched binary models

Specified before any of it was run.
- Pipelines: at diagnosis (Frame P; EHR at presentation, first ECG,
  presentation CTPA) and at discharge (Frame B; EHR+, ECG, CTPA), in the
  EHR+ECG and three-modality cohorts, folds as for death.
- Outcomes, horizons 10, 30, 90, 180 and 365 days from each frame's clock:
  the composite (first of CV readmission or death) and CV readmission
  first, death competing. Frame P's outcomes from stage 4r (same-day CV
  and death count as CV, as stage 4; must reproduce Frame A to day 30);
  Frame B's as stage 4.
- Survival level 1: ridge-Cox retuned, spline Cox (D16; EHR only), RSF
  and oblique RSF (Brier-tuned, tune_expand.R) on every modality. The
  composite as a single event. CV readmission cause-specific (deaths
  censored); the RSF as a competing-risks forest (Brier on CV first by
  each horizon); the oblique RSF cause-specific.
- Survival super learners (pooled weights, recalibrated per horizon):
  the composite exactly as death (EHR-alone and PH/TV/TVC candidates); CV
  readmission targets its cumulative incidence through per-horizon
  logistic stacks on "CV readmission first by the horizon" (death first
  counts as 0), with each modality's CV and death scores of the same
  family as inputs, plus EHR-alone candidates; hazard stacks are not used
  for CV because they ignore the competing risk.
- Discharge library: v7-core (the four EHR+ families with the same
  families' ECG and CTPA scores). Checked on death: a v7-core super learner
  from existing scores is representative of v7 if within 0.005 AUROC at
  every horizon in the primary cohort; otherwise a stated limitation.
- Matched binary models (primary cohort, EHR+ECG, both pipelines; death,
  composite and CV readmission): per horizon, ridge logistic regression
  and probability forests (Brier-tuned, out-of-bag) on each modality;
  binary super learner per horizon (EHR-alone and per-set logistic-stack
  candidates, simplex weights, recalibrated). Same features, folds and
  targets as the survival models (CV: CV readmission first by the
  horizon).
- Comparisons: survival minus binary super learner per pipeline, outcome
  and horizon (AUROC, IPA; calibration slope reported); each outcome's
  performance at diagnosis and at discharge. 2,000-sample bootstrap
  (clustered by patient in Frame P).
- Expected: composite AUROC below death and CV readmission lower still;
  every outcome better predicted at discharge than at diagnosis; survival
  and binary super learners within 0.01 AUROC at most horizons.
  - Secondary, dissertation-style CV analysis (added before any Option B
  model was run): at each horizon, patients who died before any CV
  readmission by that horizon are excluded, and the outcome is CV
  readmission by the horizon among the rest (at 30 days from admission
  this is the dissertation's and side project's cv_first definition).
  Binary models (primary cohort, both pipelines) are trained and
  evaluated in that per-horizon cohort, as the dissertation did; the
  survival models' CV predictions are evaluated in the same cohort on
  AUROC only, since their probabilities refer to readmission before death
  rather than to readmission given no earlier death. Reported beside the
  primary competing-risk analysis, not in place of it.
  ## A34. Option B extensions before the binary comparison: no-hospice and component-based composite

Specified before any of it was run.
- No-hospice at discharge (as A31 for death): the EHR+ families refitted
  for the composite and CV readmission without the hospice destination
  indicator (stage 6q --no-hospice; ECG and CTPA reused); super learners
  built from them; contrasts with and without the indicator, and every
  discharge super learner (death, composite, CV) also evaluated with
  hospice discharges removed. Descriptive.
- Component-based composite: in both pipelines, extra composite
  candidates whose inputs are each modality's death score and CV score
  from the same family (PH and TV stacks, per family and set), added to
  the composite library. The extended super learner replaces the
  composite super learner of stage 18r only if significantly better on
  AUROC at two or more horizons and significantly worse at none, in the
  primary cohort of that pipeline.
- Not pursued (literature review): new competing-risks learners (classical
  cause-specific models match machine-learning ones, and the CV pipeline
  already stacks cause-specific Cox, splines and competing-risks forests
  on the cumulative incidence); readmission-specific features (length of
  stay, recent utilisation, procedures, discharge medications), recorded
  as future work.
- Expected: the hospice indicator matters for the composite much as for
  death and little for CV readmission; component-based stacking changes
  composite AUROC by less than 0.01.
  A35. Temporal validation of the final models (stage 21b). Both final pipelines, at diagnosis and at discharge, are refitted on admissions up to 2013 and tested on admissions from 2014 onwards. This covers death, the composite and CV readmission. Every tuning step, super learner weight and recalibration is fitted on the development years only. Test admissions of patients who also have a development admission are excluded, so no patient appears in both sets. Results are reported for the pooled test years and for 2014 to 2016 and 2017 onwards. The main cross-validated model, scored on the same test patients, is the reference. The at-diagnosis run follows the stage 0r audit, because that audit may change the at-diagnosis features.
  A36. History known on arrival in the main at-diagnosis model (stage 4p, third version). The stage 0r audit showed that 47% of admissions had no earlier admission, so the strict rule missed much of the history a clinician would know. Therefore, a chronic condition now counts when it is coded on an earlier admission or on the index admission. However, categories that the stay can cause, or that describe the PE itself, still count from earlier admissions only. These are myocardial infarction, cerebrovascular disease, paralysis, peptic ulcer, arrhythmia, valvular disease, pulmonary circulation disorders, other neurological disorders, coagulopathy, fluid and electrolyte disorders, blood loss anaemia and weight loss. Index codes for acute right heart failure are not counted as heart failure, since the PE usually causes it. sPESI and PESI-10 use the same history. The strict version is kept as a sensitivity analysis, and from A36 onwards the Psens and “(sens)” labels refer to it. Every at-diagnosis model is refitted with these inputs, and A35’s temporal validation uses them.
  ## A37. 11 October 2026: interpretability and the CTPA block split for the A36 at-diagnosis model

Specified before anything below was run. Descriptive, extending A23 and
A28 item 4 to the at-diagnosis model (Frame P, death within 365 days of
admission).
1. Interpretability (stage 20p), as A23 with these changes. Spline Cox
   partial effects use stage 6p's EHR spline model on ehrP.csv (linear
   terms ridge-penalised as D16, optimizer "efs", D14's fallback),
   refitted once on every admission of each cohort. Features are ranked
   by the SD of their term, because with select = TRUE mgcv's chi-square
   for a heavily shrunk term is unstable. SurvSHAP(t) uses stage 6p's
   Brier-tuned RSF and oblique RSF on the EHR alone, EHR+ECG cohort, at
   stage 6p's modal settings, with A23's design (300 explained and 100
   background admissions, 32 antithetic permutation orders, 10 to 365
   days). Added summary: the share of importance held by each EHR domain
   (demographics, vital signs, laboratory tests, history, earlier care),
   where a domain's value per admission is the sum of its features'
   values.
2. Figures (stage 20q): the 8 largest partial effects; feature and domain
   shares over time for each forest; importance in all three models for
   features in any model's top 10; CTPA block hazard ratios.
3. CTPA block split (stage 16p), as A28 item 4 in Frame P's three-modality
   cohort, with two sub-blocks of the PE descriptors: right-heart strain
   signs (rv_strain, septal_bow, reflux, mpa_enlarge) and clot location
   and extent (the other 9). Block scores follow stage 6p's CTPA ridge-Cox
   including its path extension, and all 46 features together must
   reproduce stage 6p's CTPA score. Added contrasts: AUROC for death by day
   30 of each block minus the PE descriptors; right-heart strain signs
   against clot location and extent. The bootstrap resamples patients, so
   repeat admissions stay together.
- Expected: (1) cancer and age lead in all three models, with acute
  presentation measurements (urea, lactate, respiratory rate, oxygen
  saturation) close behind; vital signs and laboratory tests hold a larger
  share of forest importance at 10 days than at 365 days, and history a
  smaller one; Spearman agreement of at least 0.7 between each forest and
  the spline model. (3) At diagnosis the PE descriptors carry more early
  signal than after discharge, driven by the right-heart strain signs
  (early/late ratio above 1 and above clot location and extent); the
  comorbidity mentions remain the most durable block (highest AUROC for
  deaths in days 31-365).
  - Clarification after stages 16p and 20p ran, for presentation only: in
  the importance figures, Charlson and Elixhauser indicators that are
  identical columns are shown as one feature (SurvSHAP(t) values summed
  per admission, spline term SDs summed), and ambiguous MIMIC item names
  are clarified. The results CSVs are unchanged.
 ## A38. Stage 13 results tables and figures (12 October 2026)

Stage 13 reports the final models. It fits no new model and changes no earlier result.

**Final models.** At diagnosis, death uses Super learner P (A32, A36). At discharge, death uses Super learner v7 (A30). The composite and CV readmission use the stage 18r super learners in both pipelines (A33, A34). All predictions are out of fold, averaged over seeds and recalibrated per horizon.

**Stage 13b (calibration and decision curves).** The script reads the local out-of-fold predictions and writes aggregates only. It first reproduces the AUROC, IPA, observed/expected and calibration slope of stages 18q and 18r (FAIL beyond 1e-6). It adds the calibration intercept and slope with intervals, and ICI, E50 and E90 from a loess curve (span 0.75, degree 2, as in stage 11b). Table 2 shows ICI without an interval. Resampling noise inflates the bootstrap ICI, so its percentile interval sat above the estimate in 4 of 60 cells. Intervals use 2,000 patient-clustered bootstrap samples with seed 42. For the plots, deciles of predicted risk are merged from the lowest up until each group holds at least 10 events. No plotted group can then imply a count from 1 to 9.

**Decision curves (exploratory).** Net benefit is computed at thresholds from 0.5% to 50% for the model, treat all and treat none. At diagnosis, for death, two rules are added: sPESI 1 or more, and PESI-10 class III to V (86 points or more). PESI-10 lacks altered mental status, so its classes may understate risk. Differences in net benefit at 5%, 10% and 20% get bootstrap intervals.

**Stage 13c (Tables 2, 3, S2 and S3).** Table 2 reports events, AUROC, IPA, observed/expected, calibration slope and ICI at each horizon. AUROC, IPA and observed/expected keep their stage 18q and 18r intervals. Table 3 reports the gain over sPESI, PESI-10 and the EHR alone, with bold marking an interval that excludes zero. At discharge, the death gain over the EHR alone uses v7-core (A33), because only v7-core has a matching EHR-alone super learner. v7-core is within 0.002 AUROC of v7 at every horizon. Tables S2 and S3 repeat both tables for the three-modality cohort.

**Stage 13d (figures).** The ECG figures redraw the registered primary analysis (stage 7, H1) and its sensitivity analyses, so no estimate is new. The calibration and decision-curve figures come from stage 13b. Titles that make a claim are written from the numbers by fixed rules, so a title cannot overstate a result.

**Stage 13e (verified results).** The script reads the saved aggregate CSVs and writes one reference file of every headline number, each tagged with its source file and rows. The write-up quotes numbers from this file only. It stops if death counts disagree between stage 4r and Table 1, or if any final stage's checks file holds a FAIL.