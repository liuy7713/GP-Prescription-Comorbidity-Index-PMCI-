# 07_rsf_aggregate_refit_CI.R

# Same as 07_rsf_aggregate_refit.R but adds bootstrap CI computation
# for test C-index while the model is still in memory.

# This script is essentially the same as 07_rsf_aggregate_refit.R, with the 
# addition of Bootstrapping for 95% confidence interval on the test C-Index 
# while the model is still in memory. 

# This utilises the same Bootstrapping idea as the PMCI 95% Confidence interval 
# bootstrapping method. 

# This script follows from 07_rsf_tune_one.R. 
# It picks the best parameters based on the tuning results from 07_rsf_tune_one.R 
# and refits the best parameters again but includes variable importance 'permute' 
# OOB C-Index produced in refit should be identical (or very similar) to corresponding 
# OOB C-Index obtained during tuning. 

# This script allows to only refit the final RSF model for one sex. It allows 
# the user to submit two jobs to HPC each running the final refit for one sex to 
# speed up the procedure. 

# This file can only be run in the HPC, does not run on RStudio. 


library(data.table) 
library(randomForestSRC) 
library(survival) 


# Configurations 
DATA_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"
EN_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/en_cox"
OUT_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/rsf"

SEED <- 20260119
BLOCK_SIZE <- 100

# Bootstrap parameters (matching PMCI code)
B_BOOT <- 500
CI_LEVEL <- 0.95
BOOT_STRATIFY_EVENT <- TRUE

# Refit only one sex (This comes from export RSF_SEX="Male" or "Female" in PBS script)
ONLY_SEX <- Sys.getenv("RSF_SEX", unset = "")
if (ONLY_SEX != "") {
  cat("RSF_SEX set -> refitting only:", ONLY_SEX, "\n")
}

# Threads from PBS
NTHREAD <- as.integer(Sys.getenv("PBS_NCPUS", unset = "8"))
data.table::setDTthreads(NTHREAD)
Sys.setenv(OMP_NUM_THREADS = NTHREAD,
           MKL_NUM_THREADS = 1,
           OPENBLAS_NUM_THREADS = 1)

# Run tag must match the tuning run
RUN_TAG <- Sys.getenv("RUN_TAG")
if (RUN_TAG == "") stop("RUN_TAG is not set. Export RUN_TAG to point to the tuning run folder.")

OUT_DIR_RUN <- file.path(OUT_DIR, RUN_TAG)
CHUNK_DIR <- file.path(OUT_DIR_RUN, "tuning_chunks")
if (!dir.exists(CHUNK_DIR)) stop("Missing chunk dir: ", CHUNK_DIR)

# Harrell's C-Index Function (Same as previous models) 
cindex_risk <- function(time, event, risk_score) {
  cc <- survival::concordance(survival::Surv(time, event) ~ I(-risk_score))
  as.numeric(cc$concordance)
}

# Bootstrap CI function (matching PMCI code)
bootstrap_cindex_ci <- function(time, event, lp, B = 500, level = 0.95,
                                seed = 1, stratify_event = TRUE) {
  set.seed(seed)
  n <- length(time)
  stopifnot(length(event) == n, length(lp) == n)
  
  idx_event <- which(event == 1)
  idx_cens  <- which(event == 0)
  
  cvals <- numeric(B)
  
  for (b in seq_len(B)) {
    if (stratify_event) {
      # Sample within event strata to avoid degenerate samples
      ib <- c(
        sample(idx_event, length(idx_event), replace = TRUE),
        sample(idx_cens,  length(idx_cens),  replace = TRUE)
      )
    } else {
      ib <- sample.int(n, n, replace = TRUE)
    }
    
    cc <- survival::concordance(survival::Surv(time[ib], event[ib]) ~ I(-lp[ib]))
    cvals[b] <- as.numeric(cc$concordance)
  }
  
  alpha <- (1 - level) / 2
  ci <- as.numeric(stats::quantile(cvals, probs = c(alpha, 1 - alpha), na.rm = TRUE))
  list(cvals = cvals, lo = ci[1], hi = ci[2])
}

# Loading Train Test Split (Same as previous models) 
load_split_idx <- function(sex_name, n) {
  sfile <- file.path(EN_DIR, paste0("split_idx_", sex_name, ".rds"))
  if (!file.exists(sfile)) stop("Missing split index file: ", sfile)
  sp <- readRDS(sfile)
  tr <- sp$tr; te <- sp$te
  if (any(tr < 1 | tr > n) || any(te < 1 | te > n)) stop("split indices out of range for ", sex_name)
  list(tr = tr, te = te, seed = sp$seed)
}

# Obtaining feature columns (Same as previous models) 
get_feature_cols <- function(df_sex, sex_name) {
  bnf_cols <- grep("^bnf_", names(df_sex), value = TRUE)
  
  tmp <- df_sex[, ..bnf_cols]
  tmp <- as.data.frame(lapply(tmp, as.numeric))
  nonzero <- colSums(tmp, na.rm = TRUE) > 0
  bnf_cols <- bnf_cols[nonzero]
  if (length(bnf_cols) == 0) stop("No nonzero bnf_ columns found for ", sex_name)
  
  if (sex_name == "Male") {
    drop_male <- c(
      "bnf_Female_sex_hormones_and_their_modulators",
      "bnf_Preparations_for_vaginal_and_vulval_changes",
      "bnf_Vaginal_and_vulval_infections",
      "bnf_Combined_hormonal_contraceptives_and_systems",
      "bnf_Progestogen_only_contraceptives",
      "bnf_Spermicidal_contraceptives",
      "bnf_Emergency_contraception"
    )
    drop_male <- intersect(drop_male, bnf_cols)
    bnf_cols <- setdiff(bnf_cols, drop_male)
  }
  bnf_cols
}

# Out of Bag C-Index (Same as Code/07_rsf_tune_one.R)
oob_cindex_rsf <- function(fit) {
  1 - tail(fit$err.rate, 1)
}


# Read the Chunk files and aggregate 
chunk_files <- list.files(CHUNK_DIR, pattern = "^tune_task.*\\.csv$", full.names = TRUE)
if (length(chunk_files) == 0) stop("No tuning chunk CSVs found in: ", CHUNK_DIR)

tune_all <- rbindlist(lapply(chunk_files, fread), fill = TRUE)

fwrite(tune_all, file.path(OUT_DIR_RUN, "rsf_tuning_all.csv"))

# Pick the best per sex by Out of Bag C-Index 
# If there is tie break, then choose faster elapsed 
setorder(tune_all, sex, -oob_cindex, elapsed_sec)
best_by_sex <- tune_all[, .SD[1], by = sex]
fwrite(best_by_sex, file.path(OUT_DIR_RUN, "rsf_best_params_by_sex.csv"))

print(best_by_sex)

cat("\nSelected best params (used for refit):\n")
print(best_by_sex[, .(sex, ntree, mtry, nodesize, nsplit, oob_cindex, elapsed_sec)])


# Load the full dataset again 
df <- as.data.table(readRDS(DATA_RDS))
stopifnot(all(c("sex", "time", "event", "age_at_entry") %in% names(df)))
df[, event := as.integer(event)]
df[, time := as.numeric(time)]
df[, age_at_entry := as.numeric(age_at_entry)]

# Refit the final models for each sex 
all_metrics <- list()

# Section here added for running with one single sex 
# This can be avoided by exporting both "Male" and "Female" in the PBS script 
sexes_to_run <- unique(best_by_sex$sex) 
if (ONLY_SEX != "") sexes_to_run <- intersect(sexes_to_run, ONLY_SEX) 

for (sex_name in sexes_to_run) {
  
  cat("\n============================\n")
  cat("Final RSF refit for:", sex_name, "\n")
  cat("============================\n")
  
  bp <- best_by_sex[sex == sex_name]
  
  model_path <- file.path(OUT_DIR_RUN, paste0("model_rsf_", sex_name, "_tuned_CI.rds"))
  vimp_path <- file.path(OUT_DIR_RUN, paste0("rsf_vimp_", sex_name, "_tuned_CI.csv"))
  met_path <- file.path(OUT_DIR_RUN, paste0("metrics_", sex_name, "_rsf_CI.csv"))
  
  if (file.exists(model_path) && file.exists(vimp_path) && file.exists(met_path)) {
    cat("Outputs already exist for ", sex_name, " -> skipping refit.\n", sep = "")
    next
  }
  
  # Loading the best parameters 
  mtry <- as.integer(bp$mtry[1])
  nodesize <- as.integer(bp$nodesize[1])
  nsplit <- as.integer(bp$nsplit[1])
  ntree_best <- as.integer(bp$ntree[1])
  
  df_sex <- df[sex == sex_name]
  n <- nrow(df_sex)
  
  # Features
  bnf_cols <- get_feature_cols(df_sex, sex_name)
  rhs <- paste(c("age_at_entry", bnf_cols), collapse = " + ")
  fml <- as.formula(paste0("Surv(time, event) ~ ", rhs))
  
  # Split
  sp <- load_split_idx(sex_name, n)
  dtr <- df_sex[sp$tr]
  dte <- df_sex[sp$te]
  
  cat("Split: train=", nrow(dtr), " test=", nrow(dte), "\n", sep="")
  cat("Best params: mtry=", mtry, " nodesize=", nodesize, " nsplit=", nsplit, "\n", sep="")
  
  set.seed(SEED)
  t0 <- Sys.time()
  
  # Final fit function here, using importance = "permute", variable importance included 
  final_fit <- rfsrc(
    formula = fml,
    data = dtr,
    ntree = ntree_best,
    mtry = mtry,
    nodesize = nodesize,
    nsplit = nsplit,
    importance = "permute",
    nthread = NTHREAD,
    block.size = BLOCK_SIZE
  )
  
  t1 <- Sys.time()
  fit_min <- as.numeric(difftime(t1, t0, units = "mins"))
  
  # Evaluate
  pred_te <- predict(final_fit, newdata = dte)
  c_te <- cindex_risk(dte$time, dte$event, pred_te$predicted)
  c_oob_final <- oob_cindex_rsf(final_fit)
  
  cat(sprintf("OOB C-index (final)  = %.4f\n", c_oob_final))
  cat(sprintf("Test C-index (final) = %.4f\n", c_te))
  cat(sprintf("Final fit time (min) = %.2f\n", fit_min))
  
  # This section here is for calculating the bootstrap C-Index 
  cat("\n--- Bootstrapping test C-index ---\n")
  
  boot_result <- bootstrap_cindex_ci(
    time = dte$time,
    event = dte$event,
    lp = pred_te$predicted,
    B = B_BOOT,
    level = CI_LEVEL,
    seed = SEED + ifelse(sex_name == "Female", 1001, 2001),
    stratify_event = BOOT_STRATIFY_EVENT
  )
  
  cat(sprintf("Bootstrap %.0f%% CI: [%.4f, %.4f]\n", 
              100 * CI_LEVEL, boot_result$lo, boot_result$hi))
  
  # Save bootstrap distribution
  boot_dist_file <- file.path(OUT_DIR_RUN, paste0("rsf_bootstrap_dist_", sex_name, ".rds"))
  saveRDS(boot_result$cvals, boot_dist_file)
  cat("Saved bootstrap distribution to:", boot_dist_file, "\n")
  
  # Save model
  saveRDS(final_fit, file.path(OUT_DIR_RUN, paste0("model_rsf_", sex_name, "_tuned_CI.rds")))
  
  # Save vimp (Variable Importance) 
  vimp_dt <- data.table(
    feature = names(final_fit$importance),
    vimp = as.numeric(final_fit$importance)
  )
  setorder(vimp_dt, -vimp)
  fwrite(vimp_dt, file.path(OUT_DIR_RUN, paste0("rsf_vimp_", sex_name, "_tuned_CI.csv")))
  
  # Metrics (with bootstrap CI columns added)
  metrics_dt <- data.table(
    run_tag = RUN_TAG,
    sex = sex_name,
    method = "RSF",
    n = n,
    n_train = nrow(dtr),
    n_test = nrow(dte),
    ntree_final = ntree_best,
    best_mtry = mtry,
    best_nodesize = nodesize,
    best_nsplit = nsplit,
    oob_cindex_best = as.numeric(bp$oob_cindex[1]),
    oob_cindex_final = c_oob_final,
    cindex_test = c_te,
    cindex_test_lo = boot_result$lo,      # For bootstrap 
    cindex_test_hi = boot_result$hi,      # For bootstrap 
    boot_B = B_BOOT,                      # For bootstrap 
    boot_level = CI_LEVEL,                # For bootstrap 
    boot_stratify_event = BOOT_STRATIFY_EVENT,  # For bootstrap 
    final_fit_minutes = fit_min,
    nthread = NTHREAD,
    block_size = BLOCK_SIZE
  )
  
  fwrite(metrics_dt, file.path(OUT_DIR_RUN, paste0("metrics_", sex_name, "_rsf_CI.csv")))
  all_metrics[[sex_name]] <- metrics_dt
}

metrics_all <- rbindlist(all_metrics, fill = TRUE)
fwrite(metrics_all, file.path(OUT_DIR_RUN, "metrics_rsf_all_CI.csv"))

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), con = file.path(OUT_DIR_RUN, "sessionInfo_rsf.txt"))

cat("\nAll done. Outputs in:\n", OUT_DIR_RUN, "\n", sep="")

