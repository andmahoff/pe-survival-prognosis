# P2 verified results

Generated 2026-10-12 21:33 by `scripts/13e_verified_results.py` from the saved aggregate CSVs. Every number below is read by code from the file named under its table; `verified_results.csv` holds each one with its exact source rows. Intervals are 95% patient-clustered bootstrap intervals unless stated. Bold marks an interval that excludes zero; a bound printed as 0.000 then lies just beyond zero. Counts of 1 to 9 are shown as <10.


## 1. Cohorts and events

| Pipeline | Cohort | Size |
|---|---|---|
| At diagnosis | EHR+ECG | 3,501 admissions |
| At diagnosis | 3-mod | 1,636 admissions |
| At discharge | EHR+ECG | 2,909 discharges |
| At discharge | 3-mod | 1,405 discharges |

At diagnosis, the 3,501 admissions come from 3,160 patients (Table 1).
*Source: `results/optionB_outcome_events.csv` (frame, cohort, horizon).*

**At diagnosis, EHR+ECG cohort: events by each horizon, n (%)**

| Outcome | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Death | 178 (5.1%) | 397 (11.3%) | 646 (18.5%) | 801 (22.9%) | 971 (27.7%) |
| Death or CV readmission | 216 (6.2%) | 547 (15.6%) | 898 (25.6%) | 1,100 (31.4%) | 1,291 (36.9%) |
| CV readmission (death competing) | 39 (1.1%) | 172 (4.9%) | 325 (9.3%) | 418 (11.9%) | 501 (14.3%) |

**At diagnosis, 3-mod cohort: events by each horizon, n (%)**

| Outcome | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Death | 68 (4.2%) | 152 (9.3%) | 281 (17.2%) | 344 (21.0%) | 420 (25.7%) |
| Death or CV readmission | 88 (5.4%) | 219 (13.4%) | 377 (23.0%) | 469 (28.7%) | 556 (34.0%) |
| CV readmission (death competing) | 21 (1.3%) | 76 (4.6%) | 134 (8.2%) | 177 (10.8%) | 212 (13.0%) |

**At discharge, EHR+ECG cohort: events by each horizon, n (%)**

| Outcome | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Death | 81 (2.8%) | 177 (6.1%) | 351 (12.1%) | 469 (16.1%) | 620 (21.3%) |
| Death or CV readmission | 175 (6.0%) | 352 (12.1%) | 572 (19.7%) | 728 (25.0%) | 899 (30.9%) |
| CV readmission (death competing) | 101 (3.5%) | 198 (6.8%) | 289 (9.9%) | 360 (12.4%) | 431 (14.8%) |

**At discharge, 3-mod cohort: events by each horizon, n (%)**

| Outcome | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Death | 28 (2.0%) | 79 (5.6%) | 173 (12.3%) | 226 (16.1%) | 292 (20.8%) |
| Death or CV readmission | 66 (4.7%) | 152 (10.8%) | 262 (18.6%) | 342 (24.3%) | 414 (29.5%) |
| CV readmission (death competing) | 42 (3.0%) | 83 (5.9%) | 123 (8.8%) | 161 (11.5%) | 188 (13.4%) |

**Presentation data (at diagnosis)**

| Block | Item | Value |
|---|---|---|
| frame | admissions / deaths by 30 / by 365 days | 3501 / 397 / 971 |
| EHR | admissions with presentation vital signs | 73.6% |
| EHR | first labs: median coverage (registered 18 / extra) | 99.4% / 32.1% |
| EHR | features | 128 |
| EHR | cancer / heart failure / chronic lung disease (main) | 30.7% / 19.6% / 29.4% |
| EHR | the same, strict (earlier admissions only) | 18.9% / 9.4% / 15.8% |
| EHR | Charlson score median, main / strict | 2 / 0 |
| scores | sPESI 0 / 1 / 2+ (main) | 29.3% / 39.8% / 30.9% |
| scores | sPESI 0 / 1 / 2+ (strict) | 44.8% / 35.9% / 19.3% |
| scores | PESI-10 median (IQR), main | 88 (69-109) |
| ECG | admissions with an ECG within 24 h (D18 rule by dates) | 96.6% (90.9%) |
| ECG | first ECG, hours from admission: median (IQR) | -3.1 (-5.3 to -0.8) |
| ECG | first ECG before admission / after admission | 80.4% / 19.6% |
| ECG | admissions averaging ECGs with the same time | 0.0% |
| CTPA | admissions with a presentation CTPA (deaths by 365 days) | 1636 (420) |

*Source: `results/presentation_summary.csv` (all rows).*


## 2. Registered primary analysis (H1), Frame B, EHR+ECG cohort

Deaths are split at day 30 after discharge. The estimand is the ECG hazard ratio in days 0 to 30 over that in days 31 to 365.

| Quantity | Value (95% bootstrap interval) |
|---|---|
| patients | 2,909 |
| deaths days 0-30 | 177 |
| deaths days 31-365 | 443 |
| HR per SD, ECG, days 0-30 | 1.36 (1.12 to 1.83) |
| HR per SD, ECG, days 31-365 | 1.08 (0.94 to 1.37) |
| HR per SD, EHR, days 0-30 | 2.37 (2.01 to 3.57) |
| HR per SD, EHR, days 31-365 | 2.51 (2.23 to 3.36) |
| ECG ratio, early / late (primary) | 1.26 (1.06 to 1.51) |
| Decision | H1 SUPPORTED: ECG prognostic signal concentrated in days 0-30 |
| EHR ratio, early / late (no decision rule) | 0.95 (0.83 to 1.13) |

*Source: `results/primary.csv` (bootstrap percentile rows).*

**The ECG's early / late ratio across sensitivity analyses**

| Analysis | Source | ECG ratio (95% interval) |
|---|---|---|
| Registered primary analysis | primary.csv | 1.26 (1.06 to 1.51) |
| Scores winsorised | primary_robust.csv | 1.27 (1.06 to 1.55) |
| Scores rank-normalised | primary_robust.csv | 1.24 (1.04 to 1.50) |
| Hospice discharges excluded | primary_hospice.csv | 1.32 (1.07 to 1.65) |
| Three-modality cohort | s5_primary.csv | 1.25 (0.96 to 1.63) |
| Three-modality cohort, adjusted for CTPA | s5_primary.csv | 1.23 (0.93 to 1.64) |
| Retuned EHR score, 28 features | ab_ecg_ratio.csv | 1.24 (1.04 to 1.48) |
| Larger EHR score, 79 features | ab_ecg_ratio.csv | 1.23 (1.03 to 1.49) |
| 79 features plus missing indicators | ab_ecg_ratio.csv | 1.23 (1.03 to 1.49) |


## 3. Final models: discrimination and calibration (Tables 2 and S2)

At diagnosis, death uses Super learner P. At discharge, death uses Super learner v7. The composite and CV readmission use the stage 18r super learners. ICI has no interval (stage 13c note).

### Table 2: EHR+ECG cohort

**At diagnosis: Death**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 178 | 397 | 646 | 801 | 971 |
| AUROC (95% CI) | 0.855 (0.832 to 0.878) | 0.867 (0.850 to 0.884) | 0.882 (0.869 to 0.895) | 0.885 (0.873 to 0.897) | 0.889 (0.877 to 0.901) |
| IPA (95% CI) | 0.123 (0.082 to 0.162) | 0.237 (0.197 to 0.276) | 0.331 (0.296 to 0.366) | 0.368 (0.334 to 0.401) | 0.405 (0.373 to 0.436) |
| Observed / expected (95% CI) | 1.00 (0.86 to 1.13) | 1.00 (0.92 to 1.08) | 1.00 (0.94 to 1.06) | 1.00 (0.95 to 1.05) | 1.00 (0.95 to 1.05) |
| Calibration slope (95% CI) | 1.02 (0.91 to 1.16) | 1.02 (0.93 to 1.12) | 1.02 (0.94 to 1.10) | 1.02 (0.95 to 1.10) | 1.02 (0.95 to 1.10) |
| ICI | 0.008 | 0.008 | 0.012 | 0.015 | 0.018 |
| Calibration intercept (95% CI) | −0.01 (−0.18 to 0.15) | 0.00 (−0.12 to 0.11) | 0.00 (−0.11 to 0.11) | 0.00 (−0.11 to 0.11) | 0.00 (−0.11 to 0.11) |

**At diagnosis: Death or CV readmission**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 216 | 547 | 898 | 1,100 | 1,291 |
| AUROC (95% CI) | 0.830 (0.807 to 0.854) | 0.845 (0.828 to 0.861) | 0.867 (0.853 to 0.880) | 0.869 (0.856 to 0.882) | 0.874 (0.861 to 0.887) |
| IPA (95% CI) | 0.100 (0.067 to 0.130) | 0.223 (0.188 to 0.257) | 0.343 (0.309 to 0.377) | 0.374 (0.342 to 0.406) | 0.408 (0.377 to 0.440) |
| Observed / expected (95% CI) | 1.00 (0.88 to 1.12) | 1.00 (0.94 to 1.07) | 1.00 (0.95 to 1.05) | 1.00 (0.96 to 1.04) | 1.00 (0.96 to 1.04) |
| Calibration slope (95% CI) | 1.02 (0.91 to 1.15) | 1.02 (0.93 to 1.11) | 1.02 (0.94 to 1.10) | 1.02 (0.95 to 1.10) | 1.02 (0.94 to 1.10) |
| ICI | 0.007 | 0.016 | 0.015 | 0.018 | 0.023 |
| Calibration intercept (95% CI) | 0.00 (−0.16 to 0.14) | 0.00 (−0.10 to 0.10) | 0.00 (−0.10 to 0.10) | 0.00 (−0.10 to 0.10) | 0.00 (−0.09 to 0.10) |

**At diagnosis: CV readmission, with death as a competing event**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 39 | 172 | 325 | 418 | 501 |
| AUROC (95% CI) | 0.737 (0.650 to 0.819) | 0.790 (0.758 to 0.821) | 0.805 (0.779 to 0.830) | 0.805 (0.782 to 0.826) | 0.799 (0.777 to 0.820) |
| AUROC excluding earlier deaths (95% CI) | 0.743 (0.657 to 0.822) | 0.810 (0.781 to 0.839) | 0.838 (0.814 to 0.860) | 0.842 (0.821 to 0.862) | 0.843 (0.823 to 0.863) |
| IPA (95% CI) | 0.011 (−0.002 to 0.022) | 0.059 (0.030 to 0.085) | 0.129 (0.091 to 0.165) | 0.167 (0.130 to 0.203) | 0.174 (0.141 to 0.208) |
| Observed / expected (95% CI) | 0.99 (0.69 to 1.29) | 0.99 (0.85 to 1.14) | 0.99 (0.89 to 1.09) | 0.99 (0.91 to 1.08) | 0.99 (0.92 to 1.07) |
| Calibration slope (95% CI) | 1.06 (0.70 to 1.43) | 1.01 (0.87 to 1.16) | 1.01 (0.90 to 1.13) | 1.01 (0.91 to 1.12) | 1.02 (0.93 to 1.12) |
| ICI | 0.004 | 0.006 | 0.007 | 0.007 | 0.006 |
| Calibration intercept (95% CI) | −0.01 (−0.38 to 0.26) | −0.02 (−0.19 to 0.15) | −0.01 (−0.14 to 0.11) | −0.01 (−0.13 to 0.10) | −0.01 (−0.12 to 0.10) |

**At discharge: Death**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 81 | 177 | 351 | 469 | 620 |
| AUROC (95% CI) | 0.936 (0.913 to 0.957) | 0.925 (0.907 to 0.940) | 0.910 (0.895 to 0.924) | 0.907 (0.894 to 0.920) | 0.900 (0.888 to 0.913) |
| IPA (95% CI) | 0.249 (0.160 to 0.324) | 0.314 (0.250 to 0.372) | 0.362 (0.315 to 0.407) | 0.387 (0.346 to 0.427) | 0.409 (0.375 to 0.446) |
| Observed / expected (95% CI) | 0.98 (0.79 to 1.16) | 0.99 (0.87 to 1.11) | 1.00 (0.92 to 1.07) | 1.00 (0.93 to 1.06) | 1.00 (0.94 to 1.05) |
| Calibration slope (95% CI) | 1.02 (0.88 to 1.21) | 1.01 (0.91 to 1.14) | 1.01 (0.92 to 1.11) | 1.01 (0.93 to 1.11) | 1.01 (0.94 to 1.11) |
| ICI | 0.005 | 0.004 | 0.004 | 0.008 | 0.010 |
| Calibration intercept (95% CI) | −0.03 (−0.31 to 0.20) | −0.01 (−0.20 to 0.16) | −0.01 (−0.14 to 0.13) | −0.01 (−0.14 to 0.12) | −0.01 (−0.13 to 0.11) |

**At discharge: Death or CV readmission**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 175 | 352 | 572 | 728 | 899 |
| AUROC (95% CI) | 0.865 (0.842 to 0.888) | 0.870 (0.851 to 0.887) | 0.877 (0.862 to 0.890) | 0.880 (0.867 to 0.893) | 0.879 (0.866 to 0.892) |
| IPA (95% CI) | 0.161 (0.110 to 0.207) | 0.253 (0.209 to 0.295) | 0.329 (0.292 to 0.366) | 0.370 (0.335 to 0.404) | 0.399 (0.367 to 0.431) |
| Observed / expected (95% CI) | 0.99 (0.87 to 1.12) | 1.00 (0.91 to 1.08) | 1.00 (0.94 to 1.06) | 1.00 (0.95 to 1.05) | 1.00 (0.96 to 1.04) |
| Calibration slope (95% CI) | 1.01 (0.90 to 1.14) | 1.01 (0.92 to 1.12) | 1.01 (0.93 to 1.10) | 1.01 (0.94 to 1.10) | 1.01 (0.94 to 1.09) |
| ICI | 0.009 | 0.012 | 0.017 | 0.023 | 0.024 |
| Calibration intercept (95% CI) | −0.01 (−0.18 to 0.15) | −0.01 (−0.14 to 0.12) | 0.00 (−0.12 to 0.11) | 0.00 (−0.12 to 0.10) | 0.00 (−0.11 to 0.09) |

**At discharge: CV readmission, with death as a competing event**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 101 | 198 | 289 | 360 | 431 |
| AUROC (95% CI) | 0.837 (0.804 to 0.867) | 0.837 (0.811 to 0.862) | 0.845 (0.824 to 0.866) | 0.846 (0.826 to 0.866) | 0.839 (0.819 to 0.857) |
| AUROC excluding earlier deaths (95% CI) | 0.842 (0.808 to 0.872) | 0.846 (0.820 to 0.869) | 0.860 (0.840 to 0.879) | 0.865 (0.845 to 0.884) | 0.864 (0.846 to 0.881) |
| IPA (95% CI) | 0.050 (0.013 to 0.083) | 0.101 (0.062 to 0.136) | 0.164 (0.125 to 0.201) | 0.206 (0.165 to 0.246) | 0.225 (0.186 to 0.260) |
| Observed / expected (95% CI) | 0.99 (0.81 to 1.17) | 0.99 (0.86 to 1.12) | 0.99 (0.89 to 1.10) | 1.00 (0.91 to 1.09) | 1.00 (0.92 to 1.07) |
| Calibration slope (95% CI) | 1.00 (0.85 to 1.17) | 1.00 (0.88 to 1.13) | 1.00 (0.90 to 1.12) | 1.01 (0.91 to 1.12) | 1.01 (0.91 to 1.11) |
| ICI | 0.010 | 0.015 | 0.015 | 0.016 | 0.016 |
| Calibration intercept (95% CI) | −0.01 (−0.24 to 0.18) | −0.01 (−0.18 to 0.14) | −0.01 (−0.15 to 0.13) | −0.01 (−0.14 to 0.12) | −0.01 (−0.12 to 0.11) |

*Source: `results/tables23_long.csv` (table 2; original stage named in the CSV. Calibration intercept: `results/fin_metrics.csv`, measure citl).*


### Table S2: 3-mod cohort

**At diagnosis: Death**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 68 | 152 | 281 | 344 | 420 |
| AUROC (95% CI) | 0.877 (0.843 to 0.908) | 0.887 (0.865 to 0.908) | 0.900 (0.883 to 0.916) | 0.896 (0.879 to 0.913) | 0.905 (0.889 to 0.920) |
| IPA (95% CI) | 0.106 (0.037 to 0.167) | 0.219 (0.155 to 0.275) | 0.352 (0.298 to 0.402) | 0.373 (0.320 to 0.422) | 0.432 (0.383 to 0.477) |
| Observed / expected (95% CI) | 0.97 (0.76 to 1.19) | 0.98 (0.85 to 1.11) | 0.99 (0.91 to 1.08) | 0.99 (0.92 to 1.07) | 0.99 (0.93 to 1.06) |
| Calibration slope (95% CI) | 1.05 (0.87 to 1.28) | 1.03 (0.91 to 1.18) | 1.02 (0.92 to 1.14) | 1.02 (0.92 to 1.14) | 1.02 (0.92 to 1.14) |
| ICI | 0.006 | 0.013 | 0.020 | 0.023 | 0.026 |
| Calibration intercept (95% CI) | −0.03 (−0.32 to 0.21) | −0.02 (−0.23 to 0.15) | −0.02 (−0.18 to 0.14) | −0.02 (−0.17 to 0.13) | −0.01 (−0.17 to 0.13) |

**At diagnosis: Death or CV readmission**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 88 | 219 | 377 | 469 | 556 |
| AUROC (95% CI) | 0.832 (0.792 to 0.868) | 0.853 (0.831 to 0.875) | 0.874 (0.856 to 0.892) | 0.878 (0.859 to 0.894) | 0.889 (0.872 to 0.904) |
| IPA (95% CI) | 0.093 (0.043 to 0.139) | 0.209 (0.159 to 0.255) | 0.337 (0.288 to 0.383) | 0.376 (0.330 to 0.419) | 0.434 (0.388 to 0.475) |
| Observed / expected (95% CI) | 0.98 (0.79 to 1.18) | 0.99 (0.87 to 1.10) | 0.99 (0.92 to 1.07) | 0.99 (0.94 to 1.06) | 1.00 (0.94 to 1.05) |
| Calibration slope (95% CI) | 1.04 (0.87 to 1.26) | 1.02 (0.91 to 1.16) | 1.02 (0.92 to 1.14) | 1.02 (0.92 to 1.12) | 1.02 (0.92 to 1.13) |
| ICI | 0.008 | 0.020 | 0.024 | 0.026 | 0.028 |
| Calibration intercept (95% CI) | −0.03 (−0.27 to 0.19) | −0.02 (−0.20 to 0.14) | −0.01 (−0.16 to 0.13) | −0.01 (−0.15 to 0.12) | −0.01 (−0.15 to 0.13) |

**At diagnosis: CV readmission, with death as a competing event**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 21 | 76 | 134 | 177 | 212 |
| AUROC (95% CI) | 0.612 (0.470 to 0.744) | 0.782 (0.734 to 0.827) | 0.807 (0.771 to 0.841) | 0.814 (0.781 to 0.842) | 0.812 (0.781 to 0.840) |
| AUROC excluding earlier deaths (95% CI) | 0.620 (0.478 to 0.751) | 0.803 (0.757 to 0.845) | 0.843 (0.811 to 0.874) | 0.853 (0.825 to 0.879) | 0.856 (0.828 to 0.883) |
| IPA (95% CI) | 0.003 (−0.006 to 0.012) | 0.039 (0.008 to 0.066) | 0.095 (0.049 to 0.139) | 0.147 (0.097 to 0.192) | 0.173 (0.121 to 0.222) |
| Observed / expected (95% CI) | 0.99 (0.61 to 1.43) | 0.99 (0.79 to 1.21) | 1.00 (0.85 to 1.16) | 0.99 (0.87 to 1.13) | 0.99 (0.88 to 1.11) |
| Calibration slope (95% CI) | 1.06 (−0.20 to 2.17) | 1.04 (0.83 to 1.29) | 1.04 (0.88 to 1.24) | 1.03 (0.89 to 1.18) | 1.03 (0.89 to 1.19) |
| ICI | 0.004 | 0.009 | 0.012 | 0.016 | 0.011 |
| Calibration intercept (95% CI) | −0.01 (−0.50 to 0.37) | −0.01 (−0.26 to 0.21) | −0.01 (−0.20 to 0.18) | −0.01 (−0.18 to 0.17) | −0.01 (−0.17 to 0.15) |

**At discharge: Death**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 28 | 79 | 173 | 226 | 292 |
| AUROC (95% CI) | 0.941 (0.908 to 0.971) | 0.915 (0.886 to 0.940) | 0.906 (0.886 to 0.926) | 0.904 (0.885 to 0.922) | 0.903 (0.885 to 0.920) |
| IPA (95% CI) | 0.181 (0.040 to 0.290) | 0.253 (0.165 to 0.334) | 0.339 (0.270 to 0.400) | 0.365 (0.303 to 0.420) | 0.405 (0.349 to 0.458) |
| Observed / expected (95% CI) | 0.97 (0.66 to 1.29) | 0.99 (0.81 to 1.19) | 1.00 (0.89 to 1.12) | 1.00 (0.91 to 1.10) | 1.00 (0.93 to 1.08) |
| Calibration slope (95% CI) | 1.06 (0.83 to 1.46) | 1.02 (0.86 to 1.24) | 1.02 (0.90 to 1.17) | 1.02 (0.90 to 1.15) | 1.02 (0.91 to 1.15) |
| ICI | 0.005 | 0.004 | 0.008 | 0.013 | 0.016 |
| Calibration intercept (95% CI) | −0.04 (−0.50 to 0.32) | −0.01 (−0.29 to 0.25) | 0.00 (−0.20 to 0.20) | 0.00 (−0.19 to 0.18) | 0.01 (−0.16 to 0.17) |

**At discharge: Death or CV readmission**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 66 | 152 | 262 | 342 | 414 |
| AUROC (95% CI) | 0.855 (0.813 to 0.893) | 0.864 (0.839 to 0.890) | 0.880 (0.859 to 0.900) | 0.886 (0.868 to 0.904) | 0.891 (0.873 to 0.908) |
| IPA (95% CI) | 0.116 (0.051 to 0.173) | 0.212 (0.152 to 0.268) | 0.322 (0.266 to 0.378) | 0.375 (0.326 to 0.422) | 0.423 (0.375 to 0.471) |
| Observed / expected (95% CI) | 0.99 (0.77 to 1.21) | 0.99 (0.86 to 1.13) | 1.00 (0.91 to 1.09) | 1.00 (0.93 to 1.07) | 1.00 (0.94 to 1.06) |
| Calibration slope (95% CI) | 1.04 (0.86 to 1.25) | 1.02 (0.90 to 1.19) | 1.02 (0.90 to 1.15) | 1.02 (0.92 to 1.15) | 1.02 (0.91 to 1.15) |
| ICI | 0.008 | 0.013 | 0.020 | 0.026 | 0.027 |
| Calibration intercept (95% CI) | −0.02 (−0.30 to 0.23) | −0.01 (−0.21 to 0.18) | −0.01 (−0.17 to 0.15) | −0.01 (−0.16 to 0.14) | −0.01 (−0.16 to 0.14) |

**At discharge: CV readmission, with death as a competing event**

| Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Events | 42 | 83 | 123 | 161 | 188 |
| AUROC (95% CI) | 0.804 (0.742 to 0.859) | 0.826 (0.785 to 0.866) | 0.841 (0.809 to 0.872) | 0.850 (0.820 to 0.878) | 0.847 (0.818 to 0.875) |
| AUROC excluding earlier deaths (95% CI) | 0.811 (0.749 to 0.865) | 0.836 (0.796 to 0.875) | 0.862 (0.832 to 0.891) | 0.875 (0.848 to 0.900) | 0.877 (0.851 to 0.903) |
| IPA (95% CI) | 0.036 (−0.002 to 0.065) | 0.078 (0.024 to 0.130) | 0.130 (0.072 to 0.185) | 0.184 (0.123 to 0.239) | 0.213 (0.154 to 0.268) |
| Observed / expected (95% CI) | 0.99 (0.71 to 1.28) | 0.99 (0.78 to 1.19) | 0.99 (0.83 to 1.15) | 0.99 (0.87 to 1.12) | 0.99 (0.88 to 1.11) |
| Calibration slope (95% CI) | 1.06 (0.80 to 1.38) | 1.01 (0.84 to 1.22) | 1.00 (0.86 to 1.17) | 1.02 (0.88 to 1.18) | 1.02 (0.88 to 1.19) |
| ICI | 0.006 | 0.017 | 0.022 | 0.026 | 0.028 |
| Calibration intercept (95% CI) | −0.01 (−0.37 to 0.27) | −0.01 (−0.29 to 0.21) | −0.01 (−0.24 to 0.18) | −0.01 (−0.20 to 0.17) | −0.01 (−0.19 to 0.16) |

*Source: `results/tables23_long.csv` (table S2; original stage named in the CSV. Calibration intercept: `results/fin_metrics.csv`, measure citl).*


## 4. Added value over clinical scores and the EHR alone (Tables 3, S3)

Each cell is the final model minus the comparator. At discharge, the death rows use v7-core (A33).

### Table 3: EHR+ECG cohort

**At diagnosis: Death**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over sPESI | **0.150 (0.119 to 0.183)** | **0.147 (0.126 to 0.169)** | **0.142 (0.125 to 0.160)** | **0.144 (0.128 to 0.160)** | **0.141 (0.125 to 0.157)** |
| AUROC gain over PESI-10 | **0.101 (0.075 to 0.130)** | **0.102 (0.084 to 0.121)** | **0.088 (0.072 to 0.103)** | **0.087 (0.073 to 0.102)** | **0.082 (0.069 to 0.096)** |
| Specificity gain at sPESI's sensitivity | **0.242 (0.151 to 0.303)** | **0.246 (0.151 to 0.296)** | **0.232 (0.170 to 0.287)** | **0.231 (0.187 to 0.283)** | **0.211 (0.167 to 0.277)** |
| AUROC gain over the EHR alone | 0.007 (−0.002 to 0.016) | **0.006 (0.001 to 0.011)** | 0.003 (0.000 to 0.006) | 0.002 (−0.001 to 0.004) | **0.003 (0.000 to 0.006)** |
| IPA gain over the EHR alone | 0.006 (−0.009 to 0.020) | **0.013 (0.002 to 0.025)** | 0.008 (0.000 to 0.016) | **0.008 (0.001 to 0.014)** | **0.012 (0.004 to 0.020)** |

**At diagnosis: Death or CV readmission**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | **0.011 (0.002 to 0.020)** | **0.010 (0.004 to 0.016)** | **0.006 (0.002 to 0.010)** | **0.004 (0.001 to 0.007)** | **0.005 (0.001 to 0.008)** |
| IPA gain over the EHR alone | 0.003 (−0.011 to 0.016) | **0.015 (0.002 to 0.027)** | **0.016 (0.006 to 0.025)** | **0.015 (0.007 to 0.023)** | **0.018 (0.009 to 0.026)** |

**At diagnosis: CV readmission, with death as a competing event**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | 0.016 (−0.013 to 0.052) | 0.010 (0.000 to 0.022) | **0.010 (0.003 to 0.017)** | **0.008 (0.002 to 0.013)** | **0.008 (0.002 to 0.014)** |
| IPA gain over the EHR alone | 0.003 (−0.002 to 0.008) | 0.007 (−0.001 to 0.016) | 0.006 (−0.003 to 0.014) | 0.003 (−0.004 to 0.011) | 0.005 (−0.002 to 0.012) |

**At discharge: Death**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | **0.008 (0.002 to 0.016)** | **0.005 (0.001 to 0.009)** | 0.002 (0.000 to 0.004) | 0.001 (−0.001 to 0.002) | 0.002 (0.000 to 0.004) |
| IPA gain over the EHR alone | 0.006 (−0.020 to 0.030) | 0.006 (−0.007 to 0.019) | 0.003 (−0.005 to 0.011) | 0.004 (−0.002 to 0.009) | **0.007 (0.001 to 0.014)** |

**At discharge: Death or CV readmission**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | **0.015 (0.007 to 0.023)** | **0.007 (0.002 to 0.012)** | **0.004 (0.000 to 0.007)** | 0.002 (0.000 to 0.005) | **0.003 (0.000 to 0.005)** |
| IPA gain over the EHR alone | −0.002 (−0.018 to 0.013) | 0.005 (−0.005 to 0.017) | **0.010 (0.001 to 0.019)** | **0.009 (0.002 to 0.016)** | **0.012 (0.005 to 0.020)** |

**At discharge: CV readmission, with death as a competing event**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | 0.005 (−0.002 to 0.012) | 0.002 (−0.003 to 0.007) | 0.002 (−0.002 to 0.005) | 0.002 (−0.001 to 0.004) | 0.003 (0.000 to 0.006) |
| IPA gain over the EHR alone | 0.001 (−0.006 to 0.009) | 0.002 (−0.005 to 0.009) | 0.002 (−0.004 to 0.009) | 0.001 (−0.004 to 0.006) | 0.003 (−0.002 to 0.009) |

Bold: the 95% interval excludes zero.
*Source: `results/tables23_long.csv` (table 3).*


### Table S3: 3-mod cohort

**At diagnosis: Death**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over sPESI | **0.129 (0.081 to 0.175)** | **0.119 (0.090 to 0.152)** | **0.130 (0.107 to 0.154)** | **0.130 (0.108 to 0.154)** | **0.131 (0.111 to 0.151)** |
| AUROC gain over PESI-10 | **0.078 (0.041 to 0.116)** | **0.076 (0.049 to 0.105)** | **0.075 (0.055 to 0.097)** | **0.076 (0.055 to 0.096)** | **0.077 (0.058 to 0.096)** |
| Specificity gain at sPESI's sensitivity | **0.249 (0.203 to 0.346)** | **0.268 (0.074 to 0.339)** | **0.227 (0.138 to 0.325)** | **0.166 (0.125 to 0.273)** | **0.168 (0.124 to 0.259)** |
| AUROC gain over the EHR alone | 0.005 (−0.010 to 0.020) | 0.011 (0.000 to 0.023) | 0.004 (−0.002 to 0.010) | 0.003 (−0.002 to 0.009) | 0.005 (0.000 to 0.010) |
| IPA gain over the EHR alone | 0.010 (−0.017 to 0.039) | **0.030 (0.004 to 0.055)** | **0.019 (0.000 to 0.036)** | **0.017 (0.001 to 0.032)** | **0.022 (0.006 to 0.038)** |

**At diagnosis: Death or CV readmission**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | 0.012 (−0.004 to 0.028) | **0.017 (0.006 to 0.028)** | **0.008 (0.001 to 0.015)** | **0.010 (0.003 to 0.017)** | **0.011 (0.005 to 0.018)** |
| IPA gain over the EHR alone | 0.006 (−0.017 to 0.028) | 0.021 (−0.002 to 0.042) | **0.022 (0.005 to 0.041)** | **0.030 (0.013 to 0.047)** | **0.037 (0.020 to 0.054)** |

**At diagnosis: CV readmission, with death as a competing event**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | 0.067 (−0.007 to 0.156) | 0.015 (−0.004 to 0.035) | **0.014 (0.001 to 0.029)** | **0.016 (0.005 to 0.028)** | **0.019 (0.008 to 0.032)** |
| IPA gain over the EHR alone | 0.003 (−0.001 to 0.007) | 0.004 (−0.007 to 0.015) | 0.008 (−0.006 to 0.023) | 0.012 (−0.002 to 0.026) | **0.018 (0.003 to 0.033)** |

**At discharge: Death**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | 0.007 (−0.006 to 0.028) | 0.003 (−0.004 to 0.011) | 0.001 (−0.002 to 0.004) | 0.001 (−0.001 to 0.003) | 0.003 (0.000 to 0.005) |
| IPA gain over the EHR alone | 0.002 (−0.025 to 0.029) | 0.013 (−0.005 to 0.033) | 0.002 (−0.007 to 0.011) | 0.004 (−0.002 to 0.011) | 0.009 (0.000 to 0.018) |

**At discharge: Death or CV readmission**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | 0.008 (−0.002 to 0.021) | 0.007 (−0.001 to 0.014) | 0.003 (−0.001 to 0.008) | 0.002 (−0.002 to 0.006) | 0.004 (0.000 to 0.008) |
| IPA gain over the EHR alone | −0.009 (−0.025 to 0.007) | 0.003 (−0.013 to 0.018) | 0.011 (−0.001 to 0.023) | **0.011 (0.001 to 0.021)** | **0.017 (0.005 to 0.030)** |

**At discharge: CV readmission, with death as a competing event**

| Comparison | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| AUROC gain over the EHR alone | 0.006 (−0.006 to 0.018) | 0.003 (−0.003 to 0.009) | 0.002 (−0.003 to 0.008) | 0.002 (−0.002 to 0.005) | **0.004 (0.000 to 0.008)** |
| IPA gain over the EHR alone | 0.003 (−0.003 to 0.009) | −0.005 (−0.012 to 0.002) | 0.000 (−0.007 to 0.009) | −0.001 (−0.007 to 0.006) | 0.005 (−0.003 to 0.012) |

Bold: the 95% interval excludes zero.
*Source: `results/tables23_long.csv` (table S3).*


### The scores' and the EHR-alone model's own AUROC at diagnosis (death)

**EHR+ECG cohort**

| Model or score (AUROC) | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Frame A v4, recalibrated | 0.904 (0.882 to 0.925) | 0.901 (0.885 to 0.915) |  |  |  |
| PESI-10 | 0.754 (0.721 to 0.786) | 0.765 (0.743 to 0.787) | 0.795 (0.777 to 0.813) | 0.798 (0.781 to 0.814) | 0.807 (0.791 to 0.822) |
| Super learner P, recalibrated | 0.855 (0.832 to 0.878) | 0.867 (0.850 to 0.884) | 0.882 (0.869 to 0.895) | 0.885 (0.873 to 0.897) | 0.889 (0.877 to 0.901) |
| Super learner P-EHR, recalibrated | 0.849 (0.823 to 0.873) | 0.861 (0.843 to 0.878) | 0.880 (0.866 to 0.893) | 0.884 (0.871 to 0.895) | 0.886 (0.874 to 0.898) |
| sPESI | 0.705 (0.671 to 0.738) | 0.720 (0.696 to 0.743) | 0.740 (0.721 to 0.759) | 0.741 (0.723 to 0.759) | 0.748 (0.731 to 0.765) |

Frame A's mixed-time model (v4) minus the at-diagnosis model:

| Contrast | 10 d | 30 d |
|---|---|---|
| AUROC: Frame A v4 minus Super learner P | 0.049 (0.024 to 0.075) | 0.033 (0.018 to 0.050) |
| IPA: Frame A v4 minus Super learner P | 0.188 (0.132 to 0.245) | 0.115 (0.073 to 0.158) |

**3-mod cohort**

| Model or score (AUROC) | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Frame A v4, recalibrated | 0.891 (0.855 to 0.924) | 0.903 (0.882 to 0.924) |  |  |  |
| PESI-10 | 0.799 (0.757 to 0.842) | 0.811 (0.781 to 0.840) | 0.824 (0.801 to 0.847) | 0.820 (0.797 to 0.843) | 0.828 (0.805 to 0.850) |
| Super learner P, recalibrated | 0.877 (0.843 to 0.908) | 0.887 (0.865 to 0.908) | 0.900 (0.883 to 0.916) | 0.896 (0.879 to 0.913) | 0.905 (0.889 to 0.920) |
| Super learner P-EHR, recalibrated | 0.872 (0.838 to 0.904) | 0.875 (0.853 to 0.898) | 0.895 (0.878 to 0.912) | 0.892 (0.875 to 0.909) | 0.900 (0.883 to 0.915) |
| sPESI | 0.749 (0.695 to 0.800) | 0.768 (0.734 to 0.800) | 0.769 (0.743 to 0.795) | 0.766 (0.741 to 0.790) | 0.774 (0.752 to 0.797) |

Frame A's mixed-time model (v4) minus the at-diagnosis model:

| Contrast | 10 d | 30 d |
|---|---|---|
| AUROC: Frame A v4 minus Super learner P | 0.014 (−0.019 to 0.046) | 0.017 (−0.004 to 0.035) |
| IPA: Frame A v4 minus Super learner P | 0.095 (0.025 to 0.171) | 0.086 (0.030 to 0.142) |

*Source: `results/slP_summary.csv` (evaluation 'all').*


## 5. Decision curves (EHR+ECG cohort; exploratory)

Net benefit counts correctly flagged patients per patient, minus a penalty for false alarms. Cells are the model's net benefit minus each strategy's, with 95% intervals.

**At diagnosis**

| Outcome | Horizon | Threshold | Model net benefit | Minus treat all | Minus sPESI 1+ | Minus PESI-10 III to V |
|---|---|---|---|---|---|---|
| Death | 30 d | 5% | 0.087 (0.077 to 0.098) | 0.021 (0.018 to 0.023) | 0.011 (0.008 to 0.014) | 0.011 (0.008 to 0.015) |
| Death | 30 d | 10% | 0.068 (0.058 to 0.078) | 0.053 (0.048 to 0.058) | 0.027 (0.022 to 0.032) | 0.018 (0.013 to 0.023) |
| Death | 30 d | 20% | 0.044 (0.035 to 0.053) | 0.152 (0.142 to 0.160) | 0.086 (0.077 to 0.094) | 0.054 (0.045 to 0.062) |
| Death | 365 d | 5% | 0.250 (0.233 to 0.267) | 0.010 (0.008 to 0.013) | 0.010 (0.006 to 0.014) | 0.023 (0.018 to 0.030) |
| Death | 365 d | 10% | 0.231 (0.214 to 0.248) | 0.034 (0.030 to 0.038) | 0.017 (0.013 to 0.022) | 0.022 (0.016 to 0.028) |
| Death | 365 d | 20% | 0.198 (0.181 to 0.216) | 0.101 (0.093 to 0.109) | 0.046 (0.038 to 0.053) | 0.030 (0.022 to 0.038) |
| Death or CV readmission | 30 d | 5% | 0.127 (0.115 to 0.139) | 0.015 (0.012 to 0.017) |  |  |
| Death or CV readmission | 30 d | 10% | 0.105 (0.093 to 0.118) | 0.042 (0.038 to 0.047) |  |  |
| Death or CV readmission | 30 d | 20% | 0.070 (0.059 to 0.082) | 0.125 (0.115 to 0.134) |  |  |
| Death or CV readmission | 365 d | 5% | 0.339 (0.321 to 0.358) | 0.004 (0.001 to 0.006) |  |  |
| Death or CV readmission | 365 d | 10% | 0.317 (0.298 to 0.336) | 0.019 (0.015 to 0.023) |  |  |
| Death or CV readmission | 365 d | 20% | 0.280 (0.259 to 0.299) | 0.069 (0.061 to 0.076) |  |  |
| CV readmission (death competing) | 30 d | 5% | 0.021 (0.014 to 0.028) | 0.022 (0.017 to 0.025) |  |  |
| CV readmission (death competing) | 30 d | 10% | 0.012 (0.006 to 0.017) | 0.068 (0.062 to 0.074) |  |  |
| CV readmission (death competing) | 30 d | 20% | 0.001 (−0.002 to 0.005) | 0.190 (0.181 to 0.198) |  |  |
| CV readmission (death competing) | 365 d | 5% | 0.105 (0.093 to 0.118) | 0.007 (0.004 to 0.010) |  |  |
| CV readmission (death competing) | 365 d | 10% | 0.079 (0.067 to 0.092) | 0.031 (0.025 to 0.038) |  |  |
| CV readmission (death competing) | 365 d | 20% | 0.050 (0.039 to 0.061) | 0.121 (0.111 to 0.131) |  |  |


**At discharge**

| Outcome | Horizon | Threshold | Model net benefit | Minus treat all |
|---|---|---|---|---|
| Death | 30 d | 5% | 0.044 (0.035 to 0.052) | 0.032 (0.029 to 0.035) |
| Death | 30 d | 10% | 0.035 (0.027 to 0.044) | 0.079 (0.074 to 0.084) |
| Death | 30 d | 20% | 0.025 (0.018 to 0.033) | 0.199 (0.192 to 0.206) |
| Death | 365 d | 5% | 0.189 (0.173 to 0.203) | 0.017 (0.014 to 0.019) |
| Death | 365 d | 10% | 0.170 (0.154 to 0.184) | 0.044 (0.039 to 0.048) |
| Death | 365 d | 20% | 0.141 (0.126 to 0.156) | 0.124 (0.116 to 0.133) |
| Death or CV readmission | 30 d | 5% | 0.095 (0.084 to 0.106) | 0.020 (0.017 to 0.023) |
| Death or CV readmission | 30 d | 10% | 0.076 (0.065 to 0.087) | 0.053 (0.047 to 0.058) |
| Death or CV readmission | 30 d | 20% | 0.049 (0.038 to 0.060) | 0.148 (0.138 to 0.157) |
| Death or CV readmission | 365 d | 5% | 0.279 (0.261 to 0.295) | 0.006 (0.003 to 0.009) |
| Death or CV readmission | 365 d | 10% | 0.258 (0.240 to 0.275) | 0.025 (0.021 to 0.030) |
| Death or CV readmission | 365 d | 20% | 0.223 (0.205 to 0.240) | 0.087 (0.078 to 0.095) |
| CV readmission (death competing) | 30 d | 5% | 0.041 (0.033 to 0.050) | 0.022 (0.018 to 0.026) |
| CV readmission (death competing) | 30 d | 10% | 0.030 (0.022 to 0.038) | 0.065 (0.059 to 0.071) |
| CV readmission (death competing) | 30 d | 20% | 0.008 (0.002 to 0.014) | 0.173 (0.163 to 0.182) |
| CV readmission (death competing) | 365 d | 5% | 0.114 (0.101 to 0.127) | 0.010 (0.007 to 0.013) |
| CV readmission (death competing) | 365 d | 10% | 0.094 (0.082 to 0.107) | 0.041 (0.035 to 0.047) |
| CV readmission (death competing) | 365 d | 20% | 0.068 (0.056 to 0.080) | 0.132 (0.123 to 0.141) |

*Source: `results/fin_dca_key.csv` (cohort EHR+ECG; 10, 90 and 180 days are in the CSV).*


## 6. Survival super learner minus matched binary super learner (AUROC)

Positive values favour the survival model. CV excl: patients who died first are left out, as in the dissertation.

| Pipeline | Outcome | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|---|
| At diagnosis | Death | 0.007 (−0.010 to 0.023) | −0.004 (−0.011 to 0.004) | 0.001 (−0.002 to 0.004) | 0.002 (−0.001 to 0.004) | 0.000 (−0.003 to 0.002) |
| At diagnosis | Death or CV readmission | 0.009 (−0.006 to 0.024) | 0.004 (−0.001 to 0.010) | 0.002 (−0.001 to 0.004) | 0.001 (−0.002 to 0.003) | −0.001 (−0.003 to 0.002) |
| At diagnosis | CV readmission (death competing) | −0.036 (−0.088 to 0.011) | **0.016 (0.005 to 0.027)** | **0.007 (0.001 to 0.013)** | **0.005 (0.000 to 0.009)** | 0.002 (−0.002 to 0.006) |
| At diagnosis | CV readmission (CV excl) | −0.036 (−0.088 to 0.014) | **0.013 (0.002 to 0.022)** | −0.003 (−0.009 to 0.003) | **−0.010 (−0.016 to −0.003)** | **−0.015 (−0.022 to −0.008)** |
| At discharge | Death | 0.007 (−0.011 to 0.027) | 0.003 (−0.003 to 0.010) | 0.004 (0.000 to 0.008) | 0.000 (−0.002 to 0.003) | −0.001 (−0.003 to 0.001) |
| At discharge | Death or CV readmission | **0.012 (0.001 to 0.023)** | 0.005 (−0.001 to 0.011) | **0.003 (0.000 to 0.006)** | −0.001 (−0.004 to 0.002) | −0.002 (−0.005 to 0.000) |
| At discharge | CV readmission (death competing) | 0.013 (−0.004 to 0.030) | 0.005 (−0.004 to 0.014) | 0.006 (0.000 to 0.012) | **0.006 (0.001 to 0.011)** | **0.005 (0.000 to 0.010)** |
| At discharge | CV readmission (CV excl) | 0.010 (−0.005 to 0.026) | 0.001 (−0.008 to 0.009) | 0.001 (−0.004 to 0.006) | −0.003 (−0.007 to 0.001) | **−0.009 (−0.013 to −0.005)** |

Within 0.01 AUROC of each other: 32 of 40 cells (A33's rule). Bold: the interval excludes zero.
*Source: `results/binsurv_summary.csv` (all rows).*


## 7. Temporal validation (development 2010 to 2013, test 2014 onwards)

The super learner is refitted on the development years and tested on later years. The internal CV model is the main cross-validated model scored on the same test admissions.

### At diagnosis

| Outcome | Period | Admissions | Events by 30 d | Events by 365 d | Risk by 365 d |
|---|---|---|---|---|---|
| Death | development | 1,479 | 159 | 381 | 25.8% |
| Death | 2014-2016 | 965 | 108 | 271 | 28.1% |
| Death | 2017 onwards | 1,011 | 124 | 304 | 30.1% |
| Death or CV readmission | development | 1,479 | 226 | 532 | 36.0% |
| Death or CV readmission | 2014-2016 | 965 | 148 | 362 | 37.5% |
| Death or CV readmission | 2017 onwards | 1,011 | 163 | 376 | 37.2% |
| CV readmission (death competing) | development | 1,479 | 79 | 220 | 17.1% |
| CV readmission (death competing) | 2014-2016 | 965 | 43 | 134 | 16.3% |
| CV readmission (death competing) | 2017 onwards | 1,011 | 46 | 137 | 16.1% |

*Source: `results/tv2_P_drift.csv` (cohort EHR+ECG).*

| Outcome | Model or contrast | Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|---|---|
| Death | Super learner | AUROC | 0.868 (0.838 to 0.896) | 0.861 (0.839 to 0.883) | 0.879 (0.862 to 0.897) | 0.883 (0.867 to 0.900) | 0.887 (0.872 to 0.903) |
| Death | Super learner | O/E | 0.82 (0.67 to 0.96) | 0.89 (0.80 to 0.99) | 0.99 (0.92 to 1.07) | 1.03 (0.97 to 1.10) | 0.97 (0.92 to 1.03) |
| Death | Super learner | Slope | 1.11 | 1.00 | 1.07 | 1.13 | 1.07 |
| Death | Internal CV model | AUROC | 0.868 (0.838 to 0.896) | 0.874 (0.853 to 0.894) | 0.889 (0.872 to 0.905) | 0.892 (0.877 to 0.908) | 0.895 (0.880 to 0.909) |
| Death | Internal CV model | O/E | 0.90 (0.74 to 1.06) | 0.94 (0.84 to 1.04) | 0.98 (0.91 to 1.05) | 1.00 (0.94 to 1.06) | 0.98 (0.93 to 1.04) |
| Death | Internal CV model | Slope | 1.07 | 1.05 | 1.05 | 1.06 | 1.05 |
| Death | SL minus Internal CV model | AUROC diff | −0.001 (−0.012 to 0.011) | **−0.013 (−0.021 to −0.005)** | **−0.009 (−0.015 to −0.004)** | **−0.009 (−0.015 to −0.003)** | **−0.008 (−0.013 to −0.002)** |
| Death | SL minus SL-EHR | AUROC diff | 0.005 (−0.002 to 0.013) | **0.007 (0.003 to 0.011)** | **0.004 (0.001 to 0.007)** | 0.003 (0.000 to 0.005) | **0.004 (0.001 to 0.007)** |
| Death | SL minus sPESI | AUROC diff | **0.159 (0.114 to 0.206)** | **0.146 (0.114 to 0.178)** | **0.139 (0.117 to 0.163)** | **0.143 (0.121 to 0.165)** | **0.145 (0.124 to 0.166)** |
| Death | SL minus PESI-10 | AUROC diff | **0.121 (0.081 to 0.164)** | **0.108 (0.082 to 0.134)** | **0.092 (0.072 to 0.112)** | **0.087 (0.068 to 0.107)** | **0.083 (0.065 to 0.101)** |
| Death or CV readmission | Super learner | AUROC | 0.835 (0.802 to 0.865) | 0.841 (0.820 to 0.862) | 0.868 (0.852 to 0.885) | 0.869 (0.852 to 0.885) | 0.874 (0.858 to 0.889) |
| Death or CV readmission | Super learner | O/E | 0.72 (0.60 to 0.84) | 0.85 (0.77 to 0.93) | 0.90 (0.85 to 0.96) | 0.93 (0.88 to 0.97) | 0.92 (0.88 to 0.97) |
| Death or CV readmission | Super learner | Slope | 1.05 | 0.97 | 1.08 | 1.03 | 1.02 |
| Death or CV readmission | Internal CV model | AUROC | 0.844 (0.812 to 0.873) | 0.853 (0.832 to 0.873) | 0.879 (0.863 to 0.896) | 0.876 (0.860 to 0.892) | 0.881 (0.865 to 0.896) |
| Death or CV readmission | Internal CV model | O/E | 0.89 (0.74 to 1.03) | 0.95 (0.87 to 1.04) | 0.97 (0.92 to 1.03) | 0.98 (0.93 to 1.03) | 0.98 (0.94 to 1.02) |
| Death or CV readmission | Internal CV model | Slope | 1.07 | 1.04 | 1.08 | 1.04 | 1.04 |
| Death or CV readmission | SL minus Internal CV model | AUROC diff | −0.009 (−0.020 to 0.000) | **−0.012 (−0.019 to −0.006)** | **−0.011 (−0.016 to −0.006)** | **−0.007 (−0.012 to −0.002)** | **−0.007 (−0.012 to −0.002)** |
| Death or CV readmission | SL minus SL-EHR | AUROC diff | 0.006 (−0.005 to 0.017) | **0.011 (0.005 to 0.018)** | **0.007 (0.002 to 0.012)** | **0.006 (0.001 to 0.010)** | **0.006 (0.002 to 0.010)** |
| CV readmission (death competing) | Super learner | AUROC | 0.683 (0.535 to 0.823) | 0.780 (0.734 to 0.820) | 0.799 (0.762 to 0.830) | 0.788 (0.752 to 0.818) | 0.788 (0.756 to 0.818) |
| CV readmission (death competing) | Super learner | O/E | 0.44 (0.23 to 0.65) | 0.69 (0.56 to 0.84) | 0.76 (0.66 to 0.87) | 0.77 (0.69 to 0.86) | 0.82 (0.74 to 0.91) |
| CV readmission (death competing) | Super learner | Slope | 0.67 | 0.80 | 0.90 | 0.85 | 0.90 |
| CV readmission (death competing) | Internal CV model | AUROC | 0.648 (0.492 to 0.803) | 0.782 (0.737 to 0.822) | 0.805 (0.767 to 0.837) | 0.800 (0.765 to 0.830) | 0.799 (0.768 to 0.827) |
| CV readmission (death competing) | Internal CV model | O/E | 0.71 (0.38 to 1.06) | 0.88 (0.71 to 1.07) | 0.92 (0.80 to 1.05) | 0.93 (0.82 to 1.03) | 0.95 (0.85 to 1.05) |
| CV readmission (death competing) | Internal CV model | Slope | 0.71 | 0.96 | 1.00 | 0.98 | 1.01 |
| CV readmission (death competing) | SL minus Internal CV model | AUROC diff | **0.034 (0.013 to 0.055)** | −0.002 (−0.016 to 0.010) | −0.006 (−0.019 to 0.007) | **−0.012 (−0.024 to −0.001)** | **−0.011 (−0.022 to −0.001)** |
| CV readmission (death competing) | SL minus SL-EHR | AUROC diff | 0.015 (−0.020 to 0.055) | 0.007 (−0.013 to 0.027) | **0.012 (0.001 to 0.023)** | **0.010 (0.000 to 0.019)** | **0.011 (0.002 to 0.019)** |

*Source: `results/tv2_P_summary.csv and tv2_P_contrasts.csv` (EHR+ECG cohort, test pooled; per-period rows are in the CSVs).*


### At discharge

| Outcome | Period | Admissions | Events by 30 d | Events by 365 d | Risk by 365 d |
|---|---|---|---|---|---|
| Death | development | 1,253 | 69 | 234 | 18.7% |
| Death | 2014-2016 | 822 | 49 | 189 | 23.0% |
| Death | 2017 onwards | 834 | 59 | 197 | 23.6% |
| Death or CV readmission | development | 1,253 | 147 | 371 | 29.6% |
| Death or CV readmission | 2014-2016 | 822 | 96 | 269 | 32.7% |
| Death or CV readmission | 2017 onwards | 834 | 109 | 259 | 31.1% |
| CV readmission (death competing) | development | 1,253 | 88 | 191 | 16.3% |
| CV readmission (death competing) | 2014-2016 | 822 | 52 | 119 | 15.8% |
| CV readmission (death competing) | 2017 onwards | 834 | 58 | 121 | 15.9% |

*Source: `results/tv2_B_drift.csv` (cohort EHR+ECG).*

| Outcome | Model or contrast | Measure | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|---|---|
| Death | Super learner | AUROC | 0.938 (0.908 to 0.963) | 0.922 (0.900 to 0.942) | 0.910 (0.892 to 0.926) | 0.909 (0.892 to 0.924) | 0.901 (0.884 to 0.917) |
| Death | Super learner | O/E | 0.99 (0.76 to 1.22) | 0.81 (0.69 to 0.94) | 1.02 (0.92 to 1.12) | 1.03 (0.95 to 1.11) | 1.01 (0.94 to 1.08) |
| Death | Super learner | Slope | 1.30 | 0.99 | 1.03 | 1.07 | 1.08 |
| Death | Internal CV model | AUROC | 0.948 (0.919 to 0.971) | 0.929 (0.909 to 0.949) | 0.913 (0.895 to 0.929) | 0.912 (0.895 to 0.926) | 0.902 (0.885 to 0.917) |
| Death | Internal CV model | O/E | 0.96 (0.74 to 1.18) | 0.91 (0.78 to 1.05) | 1.00 (0.91 to 1.10) | 1.00 (0.92 to 1.08) | 1.00 (0.93 to 1.07) |
| Death | Internal CV model | Slope | 1.16 | 1.05 | 1.03 | 1.05 | 1.03 |
| Death | SL minus Internal CV model | AUROC diff | **−0.010 (−0.018 to −0.003)** | **−0.007 (−0.012 to −0.002)** | −0.003 (−0.007 to 0.002) | −0.003 (−0.007 to 0.002) | −0.001 (−0.005 to 0.003) |
| Death | SL minus SL-EHR | AUROC diff | **0.006 (0.002 to 0.012)** | 0.003 (−0.002 to 0.008) | 0.003 (0.000 to 0.006) | 0.001 (−0.001 to 0.004) | **0.003 (0.000 to 0.005)** |
| Death or CV readmission | Super learner | AUROC | 0.873 (0.845 to 0.899) | 0.874 (0.853 to 0.895) | 0.882 (0.865 to 0.900) | 0.882 (0.865 to 0.900) | 0.881 (0.865 to 0.898) |
| Death or CV readmission | Super learner | O/E | 0.82 (0.68 to 0.96) | 0.79 (0.70 to 0.88) | 0.86 (0.79 to 0.92) | 0.87 (0.81 to 0.93) | 0.90 (0.85 to 0.95) |
| Death or CV readmission | Super learner | Slope | 1.15 | 0.99 | 1.01 | 0.95 | 0.98 |
| Death or CV readmission | Internal CV model | AUROC | 0.877 (0.847 to 0.905) | 0.878 (0.856 to 0.899) | 0.887 (0.870 to 0.904) | 0.887 (0.871 to 0.904) | 0.885 (0.869 to 0.901) |
| Death or CV readmission | Internal CV model | O/E | 0.94 (0.78 to 1.10) | 0.93 (0.83 to 1.04) | 0.96 (0.89 to 1.03) | 0.97 (0.90 to 1.03) | 0.98 (0.92 to 1.03) |
| Death or CV readmission | Internal CV model | Slope | 1.09 | 1.04 | 1.05 | 1.03 | 1.03 |
| Death or CV readmission | SL minus Internal CV model | AUROC diff | −0.004 (−0.012 to 0.005) | −0.004 (−0.009 to 0.002) | **−0.005 (−0.009 to 0.000)** | **−0.005 (−0.009 to −0.001)** | −0.004 (−0.008 to 0.001) |
| Death or CV readmission | SL minus SL-EHR | AUROC diff | 0.003 (−0.008 to 0.014) | 0.003 (−0.004 to 0.011) | 0.003 (−0.001 to 0.007) | 0.000 (−0.004 to 0.004) | 0.001 (−0.003 to 0.005) |
| CV readmission (death competing) | Super learner | AUROC | 0.828 (0.779 to 0.874) | 0.836 (0.802 to 0.868) | 0.842 (0.810 to 0.870) | 0.832 (0.801 to 0.859) | 0.832 (0.805 to 0.857) |
| CV readmission (death competing) | Super learner | O/E | 0.67 (0.50 to 0.85) | 0.75 (0.63 to 0.88) | 0.74 (0.64 to 0.84) | 0.76 (0.67 to 0.85) | 0.80 (0.72 to 0.88) |
| CV readmission (death competing) | Super learner | Slope | 0.93 | 0.91 | 0.91 | 0.84 | 0.89 |
| CV readmission (death competing) | Internal CV model | AUROC | 0.838 (0.795 to 0.879) | 0.848 (0.818 to 0.877) | 0.855 (0.828 to 0.880) | 0.851 (0.824 to 0.875) | 0.849 (0.824 to 0.874) |
| CV readmission (death competing) | Internal CV model | O/E | 0.86 (0.65 to 1.10) | 0.91 (0.76 to 1.07) | 0.90 (0.78 to 1.03) | 0.91 (0.81 to 1.02) | 0.93 (0.84 to 1.03) |
| CV readmission (death competing) | Internal CV model | Slope | 0.97 | 1.02 | 1.02 | 0.99 | 1.03 |
| CV readmission (death competing) | SL minus Internal CV model | AUROC diff | −0.010 (−0.030 to 0.009) | −0.012 (−0.029 to 0.004) | −0.013 (−0.027 to 0.000) | **−0.018 (−0.030 to −0.007)** | **−0.017 (−0.028 to −0.006)** |
| CV readmission (death competing) | SL minus SL-EHR | AUROC diff | 0.007 (−0.005 to 0.022) | 0.003 (−0.005 to 0.011) | 0.005 (0.000 to 0.010) | 0.002 (−0.002 to 0.007) | 0.003 (−0.001 to 0.008) |

*Source: `results/tv2_B_summary.csv and tv2_B_contrasts.csv` (EHR+ECG cohort, test pooled; per-period rows are in the CSVs).*


## 8. Model choices and hospice discharges

| Choice | Result | Source |
|---|---|---|
| Discharge death model (A30) | Super learner v7, recalibrated; better at 5 of 5 horizons, worse at 0 | sl7_decision.csv |
| Composite, at diagnosis (A34 SL+ rule) | Super learner, recalibrated (SL+ better at 0, worse at 0 horizons) | optB2_decision.csv |
| Composite, at discharge (A34 SL+ rule) | Super learner, recalibrated (SL+ better at 0, worse at 4 horizons) | optB2_decision.csv |
| v7-core against v7 (A33) | largest AUROC gap 0.0016; representative at every horizon: True | optB2_core_check.csv |

**The hospice indicator's contribution at discharge (death): v7 minus v7 without it, AUROC**

| Evaluation | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| all | **0.021 (0.010 to 0.033)** | **0.016 (0.009 to 0.024)** | **0.007 (0.002 to 0.011)** | **0.005 (0.002 to 0.008)** | **0.003 (0.001 to 0.005)** |
| hospice excluded | −0.002 (−0.007 to 0.003) | −0.001 (−0.003 to 0.002) | **−0.002 (−0.003 to 0.000)** | −0.001 (−0.002 to 0.001) | 0.000 (−0.002 to 0.001) |

Bold: the 95% interval excludes zero.
*Source: `results/sl7nh_summary.csv` (measure auc, contrast v7 minus v7nh).*

**Final discharge models with hospice discharges left out (AUROC, EHR+ECG cohort)**

| Outcome | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|
| Death | 0.901 (0.862 to 0.935) | 0.897 (0.875 to 0.919) | 0.892 (0.876 to 0.909) | 0.893 (0.879 to 0.907) | 0.890 (0.876 to 0.903) |
| Death or CV readmission | 0.843 (0.816 to 0.869) | 0.851 (0.830 to 0.870) | 0.864 (0.848 to 0.879) | 0.870 (0.857 to 0.884) | 0.871 (0.857 to 0.885) |
| CV readmission (death competing) | 0.839 (0.807 to 0.870) | 0.838 (0.812 to 0.863) | 0.846 (0.826 to 0.866) | 0.846 (0.826 to 0.866) | 0.839 (0.819 to 0.858) |

*Source: `results/optB2_sl_summary.csv` (evaluation 'hospice excluded').*


## 9. At-diagnosis sensitivity analyses (death, EHR+ECG cohort, AUROC)

Strict history counts chronic conditions from earlier admissions only (A36). The second block keeps only admissions with presentation vital signs (D17).

| Patients | Model or score | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|---|
| all | Super learner P | 0.855 (0.832 to 0.878) | 0.867 (0.850 to 0.884) | 0.882 (0.869 to 0.895) | 0.885 (0.873 to 0.897) | 0.889 (0.877 to 0.901) |
| all | Super learner Psens | 0.860 (0.834 to 0.885) | 0.849 (0.830 to 0.867) | 0.856 (0.841 to 0.871) | 0.854 (0.840 to 0.868) | 0.851 (0.837 to 0.865) |
| all | sPESI | 0.705 (0.671 to 0.738) | 0.720 (0.696 to 0.743) | 0.740 (0.721 to 0.759) | 0.741 (0.723 to 0.759) | 0.748 (0.731 to 0.765) |
| all | sPESI (sens) | 0.698 (0.662 to 0.733) | 0.675 (0.650 to 0.702) | 0.702 (0.681 to 0.723) | 0.704 (0.684 to 0.724) | 0.706 (0.688 to 0.725) |
| all | PESI-10 | 0.754 (0.721 to 0.786) | 0.765 (0.743 to 0.787) | 0.795 (0.777 to 0.813) | 0.798 (0.781 to 0.814) | 0.807 (0.791 to 0.822) |
| all | PESI-10 (sens) | 0.750 (0.716 to 0.783) | 0.724 (0.700 to 0.747) | 0.757 (0.736 to 0.776) | 0.759 (0.741 to 0.777) | 0.764 (0.746 to 0.782) |
| with presentation vital signs | Super learner P | 0.844 (0.816 to 0.870) | 0.864 (0.845 to 0.881) | 0.884 (0.869 to 0.898) | 0.886 (0.872 to 0.899) | 0.892 (0.879 to 0.905) |
| with presentation vital signs | Super learner Psens | 0.853 (0.826 to 0.879) | 0.847 (0.827 to 0.866) | 0.859 (0.842 to 0.874) | 0.857 (0.840 to 0.872) | 0.854 (0.837 to 0.868) |
| with presentation vital signs | sPESI | 0.698 (0.661 to 0.734) | 0.721 (0.695 to 0.747) | 0.742 (0.721 to 0.765) | 0.742 (0.720 to 0.763) | 0.750 (0.730 to 0.770) |
| with presentation vital signs | sPESI (sens) | 0.702 (0.665 to 0.740) | 0.682 (0.654 to 0.710) | 0.708 (0.685 to 0.731) | 0.708 (0.685 to 0.730) | 0.708 (0.687 to 0.729) |
| with presentation vital signs | PESI-10 | 0.739 (0.702 to 0.775) | 0.759 (0.734 to 0.784) | 0.792 (0.772 to 0.812) | 0.795 (0.775 to 0.813) | 0.806 (0.787 to 0.824) |
| with presentation vital signs | PESI-10 (sens) | 0.744 (0.706 to 0.780) | 0.725 (0.698 to 0.751) | 0.759 (0.737 to 0.781) | 0.760 (0.739 to 0.781) | 0.764 (0.744 to 0.784) |

*Source: `results/slP_metrics.csv` (cohort EHR+ECG, measure auc).*


## 10. Interpretability at diagnosis (death, EHR+ECG cohort)

**Penalised spline Cox: the 12 largest terms (spread of each term's log hazard ratio)**. Names follow the figures (stage 20q): Blood pH, ALT and CK-MB are 'pH', 'Alanine aminotransferase alt' and 'Creatine kinase mb isoenzyme' in the CSV.

| Feature | Domain | Type | Term SD | edf |
|---|---|---|---|---|
| Charlson score | History | spline | 0.581 | 2.80 |
| Age | Demographics | spline | 0.513 | 0.99 |
| Elixhauser score | History | spline | 0.253 | 2.32 |
| Creatinine | Laboratory tests | spline | 0.213 | 3.13 |
| Urea nitrogen | Laboratory tests | spline | 0.192 | 2.05 |
| Blood pH | Laboratory tests | spline | 0.185 | 3.14 |
| ALT | Laboratory tests | spline | 0.165 | 2.36 |
| Lymphocytes | Laboratory tests | spline | 0.159 | 1.60 |
| CK-MB | Laboratory tests | spline | 0.157 | 1.50 |
| Glucose | Laboratory tests | spline | 0.138 | 1.70 |
| Bicarbonate | Laboratory tests | spline | 0.132 | 0.93 |
| Alkaline phosphatase | Laboratory tests | spline | 0.131 | 0.96 |

*Source: `results/interpP_spline_terms.csv` (cohort EHR+ECG, sorted by term_sd).*

**RSF, Brier-tuned: the 10 largest SurvSHAP(t) shares at 365 days**

| Rank | Feature | Domain | Share |
|---|---|---|---|
| 1 | Metastatic cancer (Elixhauser) | History | 10.9% |
| 2 | Metastatic cancer (Charlson) | History | 10.1% |
| 3 | Charlson score | History | 8.3% |
| 4 | Cancer (Charlson) | History | 5.2% |
| 5 | Age | Demographics | 4.9% |
| 6 | Solid tumour (Elixhauser) | History | 4.2% |
| 7 | Elixhauser score | History | 3.4% |
| 8 | Lymphocytes | Laboratory tests | 2.2% |
| 9 | Alkaline phosphatase | Laboratory tests | 1.7% |
| 10 | Albumin | Laboratory tests | 1.6% |

**oblique RSF, Brier-tuned: the 10 largest SurvSHAP(t) shares at 365 days**

| Rank | Feature | Domain | Share |
|---|---|---|---|
| 1 | Metastatic cancer (Charlson) | History | 7.7% |
| 2 | Metastatic cancer (Elixhauser) | History | 6.9% |
| 3 | Charlson score | History | 5.5% |
| 4 | Age | Demographics | 4.6% |
| 5 | Cancer (Charlson) | History | 3.8% |
| 6 | Elixhauser score | History | 3.6% |
| 7 | Solid tumour (Elixhauser) | History | 3.0% |
| 8 | Haemoglobin | Laboratory tests | 2.5% |
| 9 | MCV | Laboratory tests | 2.1% |
| 10 | Alkaline phosphatase | Laboratory tests | 2.1% |

Charlson and Elixhauser twins are identical columns, so each model splits their credit: Metastatic cancer (Charlson, Elixhauser); Chronic lung disease (Charlson, Elixhauser); Heart failure (Charlson, Elixhauser); Peripheral vascular disease (Charlson, Elixhauser); Hemiplegia or paralysis (Charlson, Elixhauser).
*Source: `results/interpP_shap_importance.csv` (horizon 365).*

**Share of importance by domain, SurvSHAP(t)**

| Model | Domain | 10 d | 30 d | 90 d | 180 d | 365 d |
|---|---|---|---|---|---|---|
| RSF, Brier-tuned | Demographics | 6.2% | 5.1% | 6.0% | 6.2% | 7.4% |
| RSF, Brier-tuned | Vital signs | 6.5% | 5.3% | 4.3% | 4.1% | 3.7% |
| RSF, Brier-tuned | Laboratory tests | 51.9% | 45.4% | 38.3% | 35.6% | 33.2% |
| RSF, Brier-tuned | History | 34.8% | 43.6% | 50.7% | 53.3% | 54.9% |
| RSF, Brier-tuned | Earlier care | 0.7% | 0.5% | 0.7% | 0.8% | 0.7% |
| oblique RSF, Brier-tuned | Demographics | 8.3% | 8.0% | 8.3% | 8.4% | 8.6% |
| oblique RSF, Brier-tuned | Vital signs | 7.1% | 6.6% | 6.0% | 5.8% | 5.3% |
| oblique RSF, Brier-tuned | Laboratory tests | 47.1% | 41.4% | 37.6% | 35.9% | 34.6% |
| oblique RSF, Brier-tuned | History | 36.2% | 42.5% | 46.9% | 48.9% | 50.4% |
| oblique RSF, Brier-tuned | Earlier care | 1.4% | 1.4% | 1.1% | 1.1% | 1.1% |

*Source: `results/interpP_shap_domains.csv` (all rows).*

**Agreement in feature ranking (Spearman)**

| Model | Comparison | Spearman |
|---|---|---|
| RSF, Brier-tuned | importance at 365 days vs spline term SD | 0.45 |
| oblique RSF, Brier-tuned | importance at 365 days vs spline term SD | 0.39 |
| RSF, Brier-tuned vs oblique RSF, Brier-tuned | importance at 365 days, forest vs forest | 0.67 |

*Source: `results/interpP_agreement.csv` (all rows).*


## 11. CTPA blocks at diagnosis (three-modality cohort, death)

**Block score on its own**

| Block | Features | HR 0 to 30 d | HR 31 to 365 d | Early / late | AUROC by 30 d | AUROC 31 to 365 d |
|---|---|---|---|---|---|---|
| All CTPA | 46 | 2.19 (1.99 to 2.45) | 1.89 (1.73 to 2.10) | 1.16 (1.03 to 1.30) | 0.799 (0.763 to 0.836) | 0.743 (0.709 to 0.778) |
| PE descriptors | 13 | 1.26 (1.08 to 1.51) | 1.26 (1.12 to 1.44) | 1.00 (0.82 to 1.24) | 0.563 (0.519 to 0.607) | 0.569 (0.532 to 0.607) |
| Incidental findings | 7 | 1.86 (1.65 to 2.12) | 1.77 (1.61 to 1.97) | 1.05 (0.91 to 1.21) | 0.713 (0.666 to 0.759) | 0.676 (0.637 to 0.715) |
| Comorbidity | 16 | 1.74 (1.59 to 1.92) | 1.70 (1.57 to 1.87) | 1.02 (0.91 to 1.15) | 0.683 (0.634 to 0.732) | 0.679 (0.640 to 0.720) |
| Support devices | 8 | 1.28 (1.20 to 1.53) | 1.09 (0.99 to 1.21) | 1.17 (1.09 to 1.41) | 0.606 (0.550 to 0.658) | 0.500 (0.457 to 0.541) |
| Report metadata | 2 | 1.31 (1.13 to 1.52) | 1.33 (1.16 to 1.51) | 0.99 (0.81 to 1.22) | 0.563 (0.513 to 0.616) | 0.562 (0.519 to 0.605) |
| PE: right-heart strain signs | 4 | 1.06 (0.96 to 1.27) | 1.20 (1.03 to 1.52) | 0.89 (0.70 to 1.10) | 0.495 (0.451 to 0.543) | 0.518 (0.480 to 0.557) |
| PE: clot location and extent | 9 | 1.25 (1.06 to 1.48) | 1.24 (1.10 to 1.42) | 1.01 (0.82 to 1.25) | 0.567 (0.521 to 0.611) | 0.567 (0.530 to 0.605) |

**Block score with the EHR score in the model**

| Block | Features | HR 0 to 30 d | HR 31 to 365 d | Early / late | AUROC by 30 d | AUROC 31 to 365 d |
|---|---|---|---|---|---|---|
| All CTPA | 46 | 1.46 (1.28 to 1.69) | 1.20 (1.04 to 1.37) | 1.22 (1.04 to 1.45) |  |  |
| PE descriptors | 13 | 1.11 (0.94 to 1.34) | 1.17 (1.01 to 1.35) | 0.95 (0.77 to 1.19) |  |  |
| Incidental findings | 7 | 1.22 (1.05 to 1.41) | 1.17 (1.02 to 1.33) | 1.05 (0.88 to 1.26) |  |  |
| Comorbidity | 16 | 1.21 (1.07 to 1.35) | 1.17 (1.04 to 1.31) | 1.03 (0.89 to 1.18) |  |  |
| Support devices | 8 | 1.18 (1.06 to 1.30) | 0.85 (0.69 to 0.99) | 1.38 (1.17 to 1.72) |  |  |
| Report metadata | 2 | 0.97 (0.83 to 1.12) | 0.98 (0.84 to 1.12) | 0.99 (0.82 to 1.20) |  |  |
| PE: right-heart strain signs | 4 | 1.02 (0.92 to 1.27) | 1.12 (0.98 to 1.48) | 0.91 (0.69 to 1.15) |  |  |
| PE: clot location and extent | 9 | 1.09 (0.93 to 1.30) | 1.16 (1.01 to 1.34) | 0.94 (0.77 to 1.17) |  |  |

*Source: `results/ctpa_blocks_P_time.csv` (all rows).*

**Block contrasts** (bold: the interval excludes no difference)

| Contrast | Model | Quantity | Estimate |
|---|---|---|---|
| PE descriptors over All CTPA | alone | ratio of early/late ratios | 0.86 (0.69 to 1.09) |
| PE descriptors over All CTPA | with EHR | ratio of early/late ratios | 0.78 (0.61 to 1.00) |
| Incidental findings over All CTPA | alone | ratio of early/late ratios | 0.91 (0.82 to 1.00) |
| Incidental findings over All CTPA | with EHR | ratio of early/late ratios | **0.86 (0.75 to 0.97)** |
| Incidental findings over PE descriptors | alone | ratio of early/late ratios | 1.05 (0.81 to 1.36) |
| Incidental findings over PE descriptors | with EHR | ratio of early/late ratios | 1.10 (0.82 to 1.48) |
| Incidental findings minus PE descriptors | alone | AUROC, deaths in days 31-365 among day-30 survivors | **0.107 (0.052 to 0.161)** |
| Incidental findings minus PE descriptors | alone | AUROC, death by day 30 | **0.149 (0.080 to 0.219)** |
| Comorbidity over All CTPA | alone | ratio of early/late ratios | **0.88 (0.81 to 0.96)** |
| Comorbidity over All CTPA | with EHR | ratio of early/late ratios | **0.84 (0.75 to 0.94)** |
| Comorbidity over PE descriptors | alone | ratio of early/late ratios | 1.02 (0.79 to 1.30) |
| Comorbidity over PE descriptors | with EHR | ratio of early/late ratios | 1.08 (0.80 to 1.42) |
| Comorbidity minus PE descriptors | alone | AUROC, deaths in days 31-365 among day-30 survivors | **0.109 (0.050 to 0.168)** |
| Comorbidity minus PE descriptors | alone | AUROC, death by day 30 | **0.120 (0.047 to 0.192)** |
| Support devices over All CTPA | alone | ratio of early/late ratios | 1.01 (0.90 to 1.25) |
| Support devices over All CTPA | with EHR | ratio of early/late ratios | 1.13 (0.91 to 1.45) |
| Support devices over PE descriptors | alone | ratio of early/late ratios | 1.17 (0.94 to 1.55) |
| Support devices over PE descriptors | with EHR | ratio of early/late ratios | **1.45 (1.09 to 1.98)** |
| Support devices minus PE descriptors | alone | AUROC, deaths in days 31-365 among day-30 survivors | **−0.069 (−0.123 to −0.016)** |
| Support devices minus PE descriptors | alone | AUROC, death by day 30 | 0.043 (−0.029 to 0.113) |
| Report metadata over All CTPA | alone | ratio of early/late ratios | 0.86 (0.70 to 1.05) |
| Report metadata over All CTPA | with EHR | ratio of early/late ratios | **0.81 (0.65 to 1.00)** |
| Report metadata over PE descriptors | alone | ratio of early/late ratios | 0.99 (0.72 to 1.36) |
| Report metadata over PE descriptors | with EHR | ratio of early/late ratios | 1.04 (0.76 to 1.42) |
| Report metadata minus PE descriptors | alone | AUROC, deaths in days 31-365 among day-30 survivors | −0.007 (−0.072 to 0.055) |
| Report metadata minus PE descriptors | alone | AUROC, death by day 30 | 0.000 (−0.071 to 0.077) |
| PE: right-heart strain signs over PE: clot location and extent | alone | ratio of early/late ratios | 0.88 (0.65 to 1.15) |
| PE: right-heart strain signs over PE: clot location and extent | with EHR | ratio of early/late ratios | 0.97 (0.69 to 1.29) |
| PE: right-heart strain signs minus PE: clot location and extent | alone | AUROC, death by day 30 | **−0.072 (−0.126 to −0.016)** |
| PE: right-heart strain signs minus PE: clot location and extent | alone | AUROC, deaths in days 31-365 among day-30 survivors | −0.050 (−0.101 to 0.004) |

*Source: `results/ctpa_blocks_P_contrasts.csv` (all rows).*


## 12. Subgroups at discharge (Super learner v7, death, AUROC)

Cells with fewer than 10 deaths have no AUROC.

| Variable | Group | Discharges | Deaths 30 d | AUROC 30 d | Deaths 365 d | AUROC 365 d |
|---|---|---|---|---|---|---|
| age | under 65 | 1,471 | 58 | 0.953 (0.925 to 0.974) | 200 | 0.925 (0.903 to 0.944) |
| age | 65-79 | 950 | 75 | 0.916 (0.886 to 0.945) | 244 | 0.891 (0.869 to 0.912) |
| age | 80 and over | 488 | 44 | 0.842 (0.779 to 0.895) | 176 | 0.797 (0.756 to 0.837) |
| sex | female | 1,501 | 81 | 0.928 (0.904 to 0.949) | 309 | 0.898 (0.880 to 0.915) |
| sex | male | 1,408 | 96 | 0.922 (0.895 to 0.945) | 311 | 0.903 (0.885 to 0.922) |
| race | White | 2,020 | 136 | 0.911 (0.889 to 0.931) | 471 | 0.893 (0.878 to 0.909) |
| race | Black | 499 | 20 | 0.955 (0.916 to 0.983) | 77 | 0.914 (0.879 to 0.945) |
| race | Hispanic or Latino | 107 | <10 | not estimated | 10 | 0.994 (0.978 to 1.000) |
| race | Asian | 59 | <10 | not estimated | 24 | 0.888 (0.797 to 0.958) |
| race | Other | 108 | <10 | not estimated | 23 | 0.855 (0.776 to 0.923) |
| race | Unknown | 116 | <10 | not estimated | 15 | 0.906 (0.831 to 0.964) |

*Source: `results/subgroup_metrics.csv` (frame B, measure auc).*


## 13. Checks across the final stages

| Checks file | PASS | WARN | FAIL | INFO |
|---|---|---|---|---|
| presentation_checks.csv | 14 | 0 | 0 | 8 |
| optionB_outcome_checks.csv | 9 | 0 | 0 | 3 |
| table1_checks.csv | 28 | 0 | 0 | 11 |
| slP_checks.csv | 9 | 0 | 0 | 2 |
| optB2_sl_checks.csv | 24 | 0 | 0 | 8 |
| sl7_checks.csv | 20 | 0 | 0 | 3 |
| binsurv_checks.csv | 17 | 0 | 0 | 1 |
| tv2_P_checks.csv | 24 | 0 | 0 | 2 |
| tv2_B_checks.csv | 24 | 0 | 0 | 2 |
| interpP_checks.csv | 9 | 0 | 0 | 1 |
| ctpa_blocks_P_checks.csv | 11 | 0 | 0 | 1 |
| subgroup_checks.csv | 5 | 0 | 0 | 7 |
| fin_checks.csv | 28 | 0 | 0 | 2 |
| tables23_checks.csv | 4 | 0 | 0 | 1 |


## Source files

| File | Last modified |
|---|---|
| ab_ecg_ratio.csv | 2026-10-02 18:23 |
| binsurv_checks.csv | 2026-10-11 00:00 |
| binsurv_summary.csv | 2026-10-11 00:00 |
| ctpa_blocks_P_checks.csv | 2026-10-11 19:52 |
| ctpa_blocks_P_contrasts.csv | 2026-10-11 19:52 |
| ctpa_blocks_P_time.csv | 2026-10-11 19:52 |
| fin_checks.csv | 2026-10-12 13:20 |
| fin_dca_key.csv | 2026-10-12 13:20 |
| fin_metrics.csv | 2026-10-12 13:20 |
| interpP_agreement.csv | 2026-10-12 07:57 |
| interpP_checks.csv | 2026-10-12 07:57 |
| interpP_fig_twins.csv | 2026-10-12 12:42 |
| interpP_shap_domains.csv | 2026-10-12 07:57 |
| interpP_shap_importance.csv | 2026-10-12 07:57 |
| interpP_spline_terms.csv | 2026-10-11 19:53 |
| optB2_core_check.csv | 2026-10-10 23:56 |
| optB2_decision.csv | 2026-10-10 23:56 |
| optB2_sl_checks.csv | 2026-10-10 23:56 |
| optB2_sl_summary.csv | 2026-10-10 23:56 |
| optionB_outcome_checks.csv | 2026-10-06 09:39 |
| optionB_outcome_events.csv | 2026-10-06 09:39 |
| presentation_checks.csv | 2026-10-09 23:52 |
| presentation_summary.csv | 2026-10-09 23:52 |
| primary.csv | 2026-09-30 13:11 |
| primary_hospice.csv | 2026-09-30 13:20 |
| primary_robust.csv | 2026-09-30 15:04 |
| s5_primary.csv | 2026-09-30 13:35 |
| sl7_checks.csv | 2026-10-05 14:44 |
| sl7_decision.csv | 2026-10-05 14:44 |
| sl7nh_summary.csv | 2026-10-05 16:39 |
| slP_checks.csv | 2026-10-10 23:44 |
| slP_metrics.csv | 2026-10-10 23:44 |
| slP_summary.csv | 2026-10-10 23:44 |
| subgroup_checks.csv | 2026-10-05 14:55 |
| subgroup_metrics.csv | 2026-10-05 14:55 |
| table1.csv | 2026-10-11 19:05 |
| table1_checks.csv | 2026-10-11 19:05 |
| tables23_checks.csv | 2026-10-12 13:54 |
| tables23_long.csv | 2026-10-12 13:54 |
| tv2_B_checks.csv | 2026-10-11 16:42 |
| tv2_B_contrasts.csv | 2026-10-11 16:42 |
| tv2_B_drift.csv | 2026-10-11 00:12 |
| tv2_B_summary.csv | 2026-10-11 16:42 |
| tv2_P_checks.csv | 2026-10-11 18:23 |
| tv2_P_contrasts.csv | 2026-10-11 18:23 |
| tv2_P_drift.csv | 2026-10-11 16:42 |
| tv2_P_summary.csv | 2026-10-11 18:23 |

