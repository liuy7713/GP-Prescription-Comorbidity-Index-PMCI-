# GP-Prescription-Comorbidity-Index-PMCI-
Yihao Liu Undergraduate Thesis Project
 
This repository contains the analysis code for the Prescription-Mediated Comorbidity Index (PMCI), developed using UK Biobank GP prescription data. The PMCI is constructed via SHARP-calibrated stability selection applied to an elastic net Cox proportional hazards model, with Random Survival Forests included as a non-linear benchmark.
 
---
 
## Data Access
 
The data used in this project are from the UK Biobank (Application ID: 677583) and are not publicly available. Access must be obtained directly through the [UK Biobank Access Management System](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access). Scripts 07 onwards require HPC submission and cannot be run locally.
 
---
 
## Main Analysis Pipeline
 
`01_build_analysis_dataset_gp_complete.R` builds the modelling-ready dataset from raw GP prescription records. It standardises BNF codes to paragraph level, restricts to chapters 1–12, retains categories prescribed to at least 500 participants, and constructs a binary ever-prescribed matrix across 124 categories.
 
`02_merge_baseline_sex_age.R` merges the GP prescription dataset with baseline sex and age variables.
 
`03_elastic_net_cox_gp.R` fits a sex-stratified elastic net Cox model as a benchmark. Alpha is selected via grid search and age is included as an unpenalised covariate.
 
`04_stability_selection.R` implements naive stability selection. For each sex, it subsamples 50% of the training set 100 times, fits a CV-tuned elastic net Cox model on each subsample, and declares features with selection proportion above 0.8 as stable.
 
`05_sharp.R` implements SHARP-calibrated stability selection with an elastic net Cox base learner, performing joint calibration over a grid of B, π, and λ to identify a data-justified operating point.
 
`06_derive_PMCI_SHARP.R` derives the PMCI from the SHARP-selected features as both an integer point score and a continuous linear predictor, and evaluates incremental predictive value over age alone using C-index with 95% bootstrap confidence intervals.
 
`07_rsf_tune_one.R` fits a single RSF model for one hyperparameter combination from `rsf_grid.csv`. Intended for HPC job array submission.
 
`07_rsf_aggregate_refit_CI.R` aggregates tuning results, selects the best hyperparameters, and refits the final RSF model with permutation variable importance and bootstrap confidence intervals.
 
---
 
## Appendix and Diagnostic Scripts
 
`appendix_A_bnf_feature_table.R` outputs the full table of all 124 BNF prescription categories.
 
`appendix_B_bnf_feature_table.R` outputs the SHARP-selected features with their corresponding PMCI point weights.
 
`data_quality_sensitivity_check.R` produces plots of GP prescription and participant counts by year to justify the early censoring date of 1990-01-01.
 
`km_pmci_pts.R` produces Kaplan-Meier survival curves stratified by PMCI score for each sex.
 
`selection_proportion_distributions_plot.R` plots selection proportion distributions for naive stability selection and SHARP side by side.
 
`sensitivity_plot_early_death_CI.R` runs a sensitivity analysis excluding early deaths at varying thresholds, with bootstrap confidence intervals on the C-index.
 
`sharp_calibration_heatmap.R` produces SHARP calibration heatmaps showing the joint optimal operating point for each sex.
 
`union_heatmap.R` produces a correlation heatmap over the union of naive and SHARP selected features for males.
 
`upset_plot.R` produces the UpSet plot comparing feature selections across SHARP and RSF for both sexes.
