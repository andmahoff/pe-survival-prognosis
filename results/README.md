# Results

This folder holds the aggregate tables written by the scripts in `scripts/`, and no file here holds a row per patient.

The best place to start is [`verified_results.md`](verified_results.md), which gathers the headline numbers that `scripts/13e_verified_results.py` reads from the CSV named under each table, and `verified_results.csv` holds every number with its source rows.

## Conventions

- `*_checks.csv`: the checks a script ran, each marked PASS, WARN, FAIL or INFO. However, three FAIL rows remain in the published files: two in `level1_linear_checks.csv` (below-chance CTPA scores, followed up in `06a_diag.py`, `06a_foldbias.py` and DEVIATIONS.md D2) and one in `pres_pull_summary.csv` (vital-sign coverage at presentation, DEVIATIONS.md D17).
- `*_metrics.csv`, `*_summary.csv`: performance by frame, cohort, model, horizon and measure, in long format, with 95% bootstrap intervals in `lo` and `hi` where computed.
- `*_contrasts.csv`: paired differences between models, with intervals.
- `*_weights.csv`, `*_tuning*.csv`, `*_lambda.csv`: fitted super learner weights and chosen tuning values.
- Frames A and P are at diagnosis (clock at admission), B is at discharge, I is INSPECT and M is its MIMIC-IV partner, and the cohorts are `EHR+ECG` and `3-mod` (with a CTPA report).
- `table1.xlsx` and `tables_2_3.xlsx` are formatted versions of `table1.csv` and of `table2.csv`, `table3.csv`, `tableS2.csv` and `tableS3.csv` respectively.

## Small cells and files left out

[`tools/suppress_small_cells.py`](../tools/suppress_small_cells.py) applied the following rules to the tables in this folder, and [`docs/suppressed_cells.csv`](../docs/suppressed_cells.csv) lists every changed cell. Firstly, counts of 1 to 9 patients or events are shown as `<10`. Secondly, where a set of cells adds up to a known total, a second cell is shown as `suppressed` (or `†` in Table 1) so the small one cannot be worked out by subtraction. Finally, rates from those counts are blank, and shares that would give 1 to 9 patients are shown as `suppressed`.

Three tables are not published, because small counts could be worked out from many of their rows:

- `cal_curves.csv` (`11b_calibration.R`): observed risk in bins of predicted risk. However, the final calibration in `fin_cal_bins.csv` merges bins until each holds at least 10 events.
- `s6_deciles.csv` (`16b_transfer_s6.R`): observed and predicted risk by decile in the transfer analysis.
- `lm_points.csv` (`14b_landmark_supermodel.R`): numbers at risk and deaths at each landmark, although the fitted curve is in `lm_curve.csv`.

In addition, terminal logs (`*_log.txt`), archived first runs of a few stages and the output of the one-off diagnostic `18z_diag_preds.R` are not included.

## Which script wrote which file

| Script | Files |
|---|---|
| `00n_inventory_presentation.R` | `inventory.csv` |
| `00p_inventory_ecg_ctpa.R` | `inventory_ecg_ctpa.csv` |
| `00q_ecg_time_check.R` | `ecg_time_check.csv` |
| `00r_audit_diagnosis_window.py` | `diagwin_ecg.csv`, `diagwin_history.csv`, `diagwin_labs.csv`, `diagwin_summary.csv` |
| `02_audit.py` | `audit_checks.csv`, `audit_horizons.csv` |
| `03_bq_followup.py` | `bq_followup_checks.csv`, `bq_followup_horizons.csv` |
| `03b_bq_ehr_extra.py` | `ehr_extra_checks.csv`, `ehr_extra_dictionary.csv`, `ehr_extra_validation.csv` |
| `03c_bq_demographics.py` | `demographics_summary.csv` |
| `03d_bq_diagnoses.py` | `dx_index_summary.csv` |
| `03e_bq_presentation.py` | `pres_pull_summary.csv` |
| `04_frames.py` | `frames_events.csv`, `frames_summary.csv` |
| `04b_inspect_frame.py` | `frames_im_events.csv`, `frames_im_summary.csv` |
| `04p_build_presentation.R` | `presentation_checks.csv`, `presentation_ecg_sources.csv`, `presentation_history.csv`, `presentation_summary.csv` |
| `04r_optionB_outcomes.R` | `optionB_outcome_checks.csv`, `optionB_outcome_events.csv` |
| `05_gates.py` | `gates.csv` |
| `05b_gates.R` | `gates_r.csv` |
| `06a_diag.py` | `level1_diag_fits.csv`, `level1_diag_scores.csv` |
| `06a_foldbias.py` | `level1_foldbias.csv` |
| `06a_level1_linear.py` | `level1_linear_checks.csv` |
| `06b_level1_forests.R` | `level1_forest_checks.csv` |
| `06c_level1_forests_wide.R` | `level1_forest_wide_checks.csv` |
| `06d_level1_extra.R` | `level1_extra_checks.csv` |
| `06e_level1_linear_retuned.R` | `level1_retuned_checks.csv` |
| `06f_level1_binary_retuned.R` | `level1_binary_retuned_checks.csv` |
| `06g_level1_forests_proper.R` | `level1_proper_checks.csv` |
| `06h_level1_enet_cboost.R` | `level1_enet_checks.csv` |
| `06i_level1_ehr79.R` | `ehr79_checks.csv`, `ehr79_tuning_summary.csv` |
| `06j_level1_ehrplus.R` (`_nh`: run with `--no-hospice`) | `ehrplus_checks.csv`, `ehrplus_features.csv`, `ehrplus_nh_checks.csv`, `ehrplus_nh_features.csv`, `ehrplus_nh_tuning_summary.csv`, `ehrplus_tuning_summary.csv` |
| `06p_level1_presentation.R` (`Psens`: run with `--sens`) | `P_level1_checks.csv`, `P_level1_tuning_summary.csv`, `Psens_level1_checks.csv`, `Psens_level1_tuning_summary.csv` |
| `06q_level1_optionB.R` (`_nh`: run with `--no-hospice`) | `Bcomp_level1_checks.csv`, `Bcomp_level1_tuning_summary.csv`, `Bcomp_nh_level1_checks.csv`, `Bcomp_nh_level1_tuning_summary.csv`, `Bcv_level1_checks.csv`, `Bcv_level1_tuning_summary.csv`, `Bcv_nh_level1_checks.csv`, `Bcv_nh_level1_tuning_summary.csv`, `Pcomp_level1_checks.csv`, `Pcomp_level1_tuning_summary.csv`, `Pcv_level1_checks.csv`, `Pcv_level1_tuning_summary.csv` |
| `06r_level1_binary.R` | `bin_B_level1_checks.csv`, `bin_B_tuning_summary.csv`, `bin_P_level1_checks.csv`, `bin_P_tuning_summary.csv` |
| `07_primary.R` | `primary.csv`, `primary_boot.csv` |
| `07b_primary_hospice.R` | `primary_hospice.csv`, `primary_hospice_boot.csv` |
| `07c_s5_three_modality.R` | `s5_primary.csv`, `s5_primary_boot.csv` |
| `08_late_fusion_linear.py` | `fusion_s1.csv`, `fusion_s3.csv`, `fusion_weights.csv` |
| `09_s2_competing.R` | `s2_checks.csv`, `s2_contrasts.csv`, `s2_metrics.csv` |
| `09a_lp_diagnostic.R` | `lp_diag.csv`, `lp_extreme_features.csv`, `primary_robust.csv` |
| `10_deep.py` | `level1_deep_checks.csv` |
| `11a_eval_s1_s3.R` | `eval_s1.csv`, `eval_s3.csv`, `eval_s3_weights.csv` |
| `11b_calibration.R` | `cal_checks.csv`, `cal_contrasts.csv`, `cal_metrics.csv` |
| `11c_eval_s4.R` | `s4_checks.csv`, `s4_contrasts.csv`, `s4_metrics.csv`, `s4_perfold.csv`, `s4_tuning.csv` |
| `11d_eval_s4_wide.R` | `s4w_checks.csv`, `s4w_contrasts.csv`, `s4w_metrics.csv` |
| `11e_eval_s4_neural.R` | `s4n_checks.csv`, `s4n_contrasts.csv`, `s4n_metrics.csv` |
| `12_sensitivity.R` | `sens_all_admissions.csv`, `sens_checks.csv`, `sens_final_year.csv`, `sens_landmark.csv` |
| `13a_table1.py` | `table1.csv`, `table1.xlsx`, `table1_checks.csv` |
| `13b_calibration_dca_final.R` | `fin_cal_bins.csv`, `fin_cal_curves.csv`, `fin_checks.csv`, `fin_dca_curves.csv`, `fin_dca_key.csv`, `fin_metrics.csv` |
| `13c_results_tables.py` | `table2.csv`, `table3.csv`, `tableS2.csv`, `tables23_checks.csv`, `tables23_long.csv`, `tableS3.csv`, `tables_2_3.xlsx` |
| `13d_fig_results.R` | `fin_fig_captions.csv` |
| `13e_verified_results.py` | `verified_checks.csv`, `verified_results.csv`, `verified_results.md` |
| `13r_compare_retuned.R` (one set per working copy compared) | `retuned2_vs_registered.csv`, `retuned2_vs_registered_summary.csv`, `retuned2_vs_registered_text.csv`, `retuned2_vs_retuned.csv`, `retuned2_vs_retuned_summary.csv`, `retuned2_vs_retuned_text.csv`, `retuned_vs_registered.csv`, `retuned_vs_registered_summary.csv`, `retuned_vs_registered_text.csv` |
| `14a_rp_curve.R` | `rp_curve.csv`, `rp_models.csv`, `rp_summary.csv` |
| `14b_landmark_supermodel.R` | `lm_curve.csv`, `lm_tests.csv` |
| `15a_dca.R` | `dca_checks.csv`, `dca_contrasts.csv`, `dca_curves.csv` |
| `15b_stability.R` | `stab_bins.csv`, `stab_seed.csv`, `stab_summary.csv`, `subgroups.csv` |
| `15c_ecg_features.R` | `ecg_features_blocks.csv`, `ecg_features_tv.csv` |
| `16a_conformal.R` | `conformal_checks.csv`, `conformal_contrast.csv`, `conformal_groups.csv`, `conformal_summary.csv` |
| `16b_transfer_s6.R` | `s6_censoring.csv`, `s6_checks.csv`, `s6_contrasts.csv`, `s6_features.csv`, `s6_lambda.csv`, `s6_performance.csv`, `s6_relatedness.csv` |
| `16c_spesi.R` | `spesi_checks.csv`, `spesi_contrasts.csv`, `spesi_distribution.csv`, `spesi_items.csv`, `spesi_metrics.csv` |
| `16d_transfer_forests.R` | `a12_checks.csv`, `a12_contrasts.csv`, `a12_gap.csv`, `a12_performance.csv`, `a12_tuning.csv`, `a12_warnings.csv` |
| `16e_ctpa_time.R` | `ctpa_checks.csv`, `ctpa_contrasts.csv`, `ctpa_cumulative.csv`, `ctpa_hr.csv`, `ctpa_incident.csv`, `ctpa_retune.csv` |
| `16f_ctpa_blocks.R` | `ctpa_blocks_checks.csv`, `ctpa_blocks_contrasts.csv`, `ctpa_blocks_time.csv` |
| `16p_ctpa_blocks_P.R` | `ctpa_blocks_P_checks.csv`, `ctpa_blocks_P_contrasts.csv`, `ctpa_blocks_P_time.csv`, `ctpa_blocks_P_vs_B.csv` |
| `17a_gap_time.R` | `gap_cells.csv`, `gap_checks.csv`, `gap_incident.csv`, `gap_regressions.csv` |
| `18a_early_linear.R` | `early_checks.csv`, `early_lambda.csv`, `early_tuning.csv`, `early_warnings.csv` |
| `18b_early_blockstd.R` | `blockstd_checks.csv`, `blockstd_lambda.csv`, `blockstd_warnings.csv` |
| `18c_arch_eval.R` | `arch_checks.csv`, `arch_contrasts.csv`, `arch_metrics.csv` |
| `18d_horizon_tuned.R` | `htune_blocks.csv`, `htune_checks.csv`, `htune_lambda.csv`, `htune_warnings.csv` |
| `18e_late_tvc.R` | `tvc_checks.csv`, `tvc_coefs.csv`, `tvc_forms.csv`, `tvc_warnings.csv` |
| `18f_super_learner.R` | `sl_checks.csv`, `sl_contrasts.csv`, `sl_metrics.csv`, `sl_weights.csv` |
| `18g_super_learner_v2.R` | `sl2_checks.csv`, `sl2_contrasts.csv`, `sl2_metrics.csv`, `sl2_weights.csv` |
| `18h_super_learner_v4.R` | `sl4_checks.csv`, `sl4_contrasts.csv`, `sl4_metrics.csv`, `sl4_weights.csv` |
| `18i_arch_forest.R` | `arch_forest_checks.csv` |
| `18j_arch_neural.py` | `arch_neural_checks.csv` |
| `18k_arch_eval2.R` | `arch2_checks.csv`, `arch2_contrasts.csv`, `arch2_metrics.csv` |
| `18l_block_forest_opt.R` | `blockopt_checks.csv`, `blockopt_params_summary.csv` |
| `18m_super_learner_v5.R` | `sl5_checks.csv`, `sl5_contrasts.csv`, `sl5_metrics.csv`, `sl5_weights.csv` |
| `18n_super_learner_v6.R` | `sl6_checks.csv`, `sl6_contrasts.csv`, `sl6_decision.csv`, `sl6_metrics.csv`, `sl6_weights.csv` |
| `18o_super_learner_v7.R` | `sl7_checks.csv`, `sl7_contrasts.csv`, `sl7_decision.csv`, `sl7_metrics.csv`, `sl7_weights.csv` |
| `18p_v7_no_hospice.R` | `sl7nh_checks.csv`, `sl7nh_contrasts.csv`, `sl7nh_metrics.csv`, `sl7nh_summary.csv`, `sl7nh_weights.csv` |
| `18q_super_learner_P.R` | `slP_checks.csv`, `slP_contrasts.csv`, `slP_metrics.csv`, `slP_summary.csv`, `slP_weights.csv` |
| `18r_optionB_super_learners.R` (`optB_`: first version, A33; `optB2_`: final version, A34) | `optB2_core_check.csv`, `optB2_decision.csv`, `optB2_sl_checks.csv`, `optB2_sl_contrasts.csv`, `optB2_sl_metrics.csv`, `optB2_sl_summary.csv`, `optB2_sl_weights.csv`, `optB_core_check.csv`, `optB_sl_checks.csv`, `optB_sl_contrasts.csv`, `optB_sl_metrics.csv`, `optB_sl_summary.csv`, `optB_sl_weights.csv` |
| `18s_binary_vs_survival.R` | `binsurv_checks.csv`, `binsurv_contrasts.csv`, `binsurv_metrics.csv`, `binsurv_summary.csv`, `binsurv_weights.csv` |
| `19a_analysis_b.R` | `ab_checks.csv`, `ab_contrasts.csv`, `ab_ecg_ratio.csv`, `ab_metrics.csv` |
| `19b_spline_cox.R` | `spline_checks.csv`, `spline_contrasts.csv`, `spline_metrics.csv` |
| `20a_interpretability.R` | `interp_agreement.csv`, `interp_checks.csv`, `interp_shap_importance.csv`, `interp_spline_curves.csv`, `interp_spline_terms.csv` |
| `20p_interpretability_P.R` | `interpP_agreement.csv`, `interpP_checks.csv`, `interpP_shap_domains.csv`, `interpP_shap_importance.csv`, `interpP_spline_curves.csv`, `interpP_spline_terms.csv` |
| `20q_fig_interpretability_P.R` | `interpP_fig_captions.csv`, `interpP_fig_twins.csv` |
| `21a_temporal.R` | `temporal_checks.csv`, `temporal_contrasts.csv`, `temporal_drift.csv`, `temporal_failures.csv`, `temporal_metrics.csv`, `temporal_weights.csv` |
| `21b_temporal_final.R` | `tv2_B_checks.csv`, `tv2_B_contrasts.csv`, `tv2_B_drift.csv`, `tv2_B_metrics.csv`, `tv2_B_summary.csv`, `tv2_P_checks.csv`, `tv2_P_contrasts.csv`, `tv2_P_drift.csv`, `tv2_P_metrics.csv`, `tv2_P_summary.csv` |
| `22a_subgroups.R` | `subgroup_checks.csv`, `subgroup_contrasts.csv`, `subgroup_metrics.csv` |
