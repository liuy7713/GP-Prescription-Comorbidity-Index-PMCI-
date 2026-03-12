# 06_derive_PMCI_SHARP.R 

# The purpose of the file is to derive the Prescription Mediated Comorbidity 
# Index (PMCI) using the stable features that I obtained in the fifth file 
# 05_sharp.R. 

# The derivation is achieved via two methods. First method follows the pipeline 
# outlined by Wanqing, using a integer based point system. The second method is 
# by employing a continuous linear predictor to produce a potentially marginally 
# better PMCI. 

# This file also produces a series of C-Index tests to evaluate models with only 
# age, age and PMCI_PTS, and age and PMCI_LP. 

# This fie now incorporates bootstrapping to calculate the 95% confidence interval 
# for the three models. 


library(data.table) 
library(survival) 

# Configurations 
DATA_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"
EN_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/en_cox"
SHARP_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/sharp"
OUT_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/sharp/pmci_sharp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE) 

SEED <- 20260120 
SCALE_FACTOR <- 5 

# Bootstrapping Configurations 
B_BOOT <- 500 
CI_LEVEL <- 0.95 
BOOT_STRATIFY_EVENT <- TRUE 

bootstrap_cindex_ci <- function(time, event, lp, B = 1000, level = 0.95, seed = 1, stratify_event = TRUE) {
  set.seed(seed)
  n <- length(time)
  stopifnot(length(event) == n, length(lp) == n)
  
  idx_event <- which(event == 1)
  idx_cens <- which(event == 0)
  
  cvals <- numeric(B)
  
  for (b in seq_len(B)) {
    if (stratify_event) {
      # sample within event strata to avoid degenerate samples with few/no events
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

# Harrell's C-Index function (same as previous implementations) 
cindex_cox <- function(time, event, lp) { 
  cc <- survival::concordance(survival::Surv(time, event) ~ I(-lp)) 
  as.numeric(cc$concordance) 
}

# Loading Train Test Split (same as previous models) 
load_split_idx <- function(sex_name, n) {
  sfile <- file.path(EN_DIR, paste0("split_idx_", sex_name, ".rds"))
  if (!file.exists(sfile)) stop("Missing split index file: ", sfile)
  sp <- readRDS(sfile)
  tr <- sp$tr 
  te <- sp$te
  if (any(tr < 1 | tr > n) || any(te < 1 | te > n)) stop("split indices out of range for ", sex_name)
  list(tr = tr, te = te, seed = sp$seed)
}

# Feature filtering function (same as used in 03_elastic_net_cox_gp.R and 04_stability_selection.R) 
get_feature_cols <- function(df_sex, sex_name) { 
  bnf_cols <- grep("^bnf_", names(df_sex), value = TRUE) 
  # Drop all zero columns in this sex 
  tmp <- df_sex[, ..bnf_cols] 
  tmp <- as.data.frame(lapply(tmp, as.numeric)) 
  nonzero <- colSums(tmp, na.rm = TRUE) > 0 
  bnf_cols <- bnf_cols[nonzero] 
  
  if (length(bnf_cols) == 0) stop("No nonzero bnf_ columns found for ", sex_name) 
  
  # Drop female specific drugs for males 
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

# Function to read the stable features from the SHARP results 
read_stable_vars <- function(sex_name) { 
  f <- file.path(SHARP_DIR, paste0("stable_variables_", sex_name, "_sharp.csv")) 
  if (!file.exists(f)) stop("Missing stable features file: ", f) 
  dt <- fread(f) 
  if ("feature" %in% names(dt)) feats <- dt$feature else feats <-  dt[[1]] 
  feats <- unique(as.character(feats)) 
  feats <- feats[!is.na(feats) & nzchar(feats)]
  feats 
}

# Function to load the refit cox model 
refit_cox <- function(df_train, stable_feats, sex_name) { 
  rds_path <- file.path(SHARP_DIR, paste0("cox_refit_", sex_name, "_sharp.rds")) 
  fit <- readRDS(rds_path) 
  return(fit) 
} 

coef_to_weights <- function(beta, scale_factor = 5) { 
  # Integer points conversion 
  round(beta * scale_factor) 
}

# Function to derive the PMCI for a given sex 
derive_pmci_one_sex <- function(df_sex, sex_name) { 
  cat("\n============================\n")
  cat("Deriving PMCI from SHARP for:", sex_name, "\n")
  cat("============================\n")
  
  set.seed(SEED) 
  req <- c("time", "event", "age_at_entry") 
  df_sex[, event := as.integer(event)] 
  df_sex[, time := as.numeric(time)] 
  df_sex[, age_at_entry := as.numeric(age_at_entry)]
  
  allowed_bnf <- get_feature_cols(df_sex, sex_name) 
  
  # Read SHARP stable variables and find them with allowed and existing columns 
  stable_feats_raw <- read_stable_vars(sex_name) 
  stable_feats <- intersect(stable_feats_raw, allowed_bnf) 
  stable_feats <- intersect(stable_feats, names(df_sex)) 
  stable_feats <- setdiff(stable_feats, "age_at_entry") 
  cat("Stable feats in file: ", length(stable_feats_raw), "\n", sep = "") 
  cat("Stable feats after sex filters and existence checks: ", length(stable_feats), "\n", sep = "") 
  
  n <- nrow(df_sex) 
  sp <- load_split_idx(sex_name, n) 
  tr <- sp$tr 
  te <- sp$te 
  dtr <- df_sex[tr] 
  dte <- df_sex[te] 
  
  # Load the unpenalised cox refit on training set 
  cox_refit <- refit_cox(dtr, stable_feats, sex_name) 
  
  # Extracting coefficients 
  b <- coef(cox_refit) 
  b <- b[!is.na(b)] 
  b_dt <- data.table(feature = names(b), beta = as.numeric(b)) 
  beta_age <- b_dt[feature == "age_at_entry", beta]
  if (length(beta_age) == 0) { 
    beta_age <- NA_real_ 
  }
  drug_dt <- b_dt[feature != "age_at_entry"] 
  # Only keep the stable features for scoring 
  drug_dt <- drug_dt[feature %in% stable_feats]
  
  # Construct the weights for the drugs 
  drug_dt[, weight := coef_to_weights(beta, SCALE_FACTOR)] 
  drug_dt[, abs_beta := abs(beta)] 
  setorder(drug_dt, -abs_beta)
  drug_dt[, abs_beta := NULL] 
  cat("Cox refit coefficients: age present = ", !is.na(beta_age), " | drug coefficient count = ", nrow(drug_dt), "\n", sep = "") 
  
  # Build the PMCI score on full sex specific dataset 
  # Define PMCI_lp as drug only 
  # Define PMCI_pts as drug only integer points 
  if (nrow(drug_dt) == 0) { 
    df_sex[, pmci_lp := 0.0] 
    df_sex[, pmci_pts := 0L] 
  } else { 
    cols <- drug_dt$feature 
    stopifnot(all(cols %in% names(df_sex))) 
    X_dt <- df_sex[, lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = cols] 
    X_drug <- as.matrix(X_dt) 
    # Ensure that X_drug is numeric 
    storage.mode(X_drug) <- "double" 
    beta_vec <- drug_dt$beta 
    names(beta_vec) <- drug_dt$feature 
    w_vec <- drug_dt$weight 
    names(w_vec) <- drug_dt$feature 
    df_sex[, pmci_lp := as.numeric(X_drug %*% beta_vec)] 
    df_sex[, pmci_pts := as.integer(X_drug %*% w_vec)] 
  }
  
  # Now only evaluate the models on the test split 
  # First start with age only 
  fit_age <- coxph(Surv(time, event) ~ age_at_entry, data = dtr, ties = "efron") 
  # Now try with age and PMCI_PTS 
  fit_pts <- coxph(Surv(time, event) ~ age_at_entry + pmci_pts, data = df_sex[tr], ties = "efron") 
  # Now try with age and PMCI_LP 
  fit_lp <- coxph(Surv(time, event) ~ age_at_entry + pmci_lp, data = df_sex[tr], ties = "efron") 
  
  # Now compute the C index on test 
  lp_age_te <- predict(fit_age, newdata = df_sex[te], type = "lp") 
  lp_pts_te <- predict(fit_pts, newdata = df_sex[te], type = "lp") 
  lp_lp_te <- predict(fit_lp, newdata = df_sex[te], type = "lp") 
  c_age <- cindex_cox(df_sex[te]$time, df_sex[te]$event, lp_age_te) 
  c_pts <- cindex_cox(df_sex[te]$time, df_sex[te]$event, lp_pts_te) 
  c_lp <- cindex_cox(df_sex[te]$time, df_sex[te]$event, lp_lp_te) 
  cat("Test C index | age only: ", round(c_age, 4), " | age + PMCI_pts: ", round(c_pts, 4), " | age + PMCI_lp: ", round(c_lp, 4), "\n", sep = "") 
  
  # Bootstrap CIs on test set (fixed fitted models)
  boot_age <- bootstrap_cindex_ci(
    time = df_sex[te]$time, event = df_sex[te]$event, lp = lp_age_te,
    B = B_BOOT, level = CI_LEVEL, seed = SEED + 101,
    stratify_event = BOOT_STRATIFY_EVENT
  )
  
  boot_pts <- bootstrap_cindex_ci(
    time = df_sex[te]$time, event = df_sex[te]$event, lp = lp_pts_te,
    B = B_BOOT, level = CI_LEVEL, seed = SEED + 202,
    stratify_event = BOOT_STRATIFY_EVENT
  )
  
  boot_lp <- bootstrap_cindex_ci(
    time = df_sex[te]$time, event = df_sex[te]$event, lp = lp_lp_te,
    B = B_BOOT, level = CI_LEVEL, seed = SEED + 303,
    stratify_event = BOOT_STRATIFY_EVENT
  )
  
  cat(sprintf("Bootstrap %.0f%% CI (test C-index)\n", 100 * CI_LEVEL))
  cat(sprintf("  age only:      %.4f [%.4f, %.4f]\n", c_age, boot_age$lo, boot_age$hi))
  cat(sprintf("  age + PMCI_pts:%.4f [%.4f, %.4f]\n", c_pts, boot_pts$lo, boot_pts$hi))
  cat(sprintf("  age + PMCI_lp: %.4f [%.4f, %.4f]\n", c_lp,  boot_lp$lo,  boot_lp$hi))
  
  # Save the weights 
  weights_out <- copy(drug_dt) 
  weights_out[, sex := sex_name] 
  weights_out[, scale_factor := SCALE_FACTOR] 
  fwrite(weights_out, file.path(OUT_DIR, paste0("pmci_weights_", sex_name, "_sharp.csv"))) 
  
  # Save the scores and keep the ID 
  id_col <- NULL 
  for (cand in c("eid", "ID", "id", "participant_id")) { 
    if (cand %in% names(df_sex)) { 
      id_col <- cand; 
      break 
    }
  }
  keep_cols <- c(id_col, "time", "event", "age_at_entry", "pmci_lp", "pmci_pts") 
  keep_cols <- keep_cols[!is.na(keep_cols)] 
  scores_dt <- df_sex[, ..keep_cols] 
  fwrite(scores_dt, file.path(OUT_DIR, paste0("pmci_scores_", sex_name, "_sharp.csv"))) 
  saveRDS(scores_dt, file.path(OUT_DIR, paste0("pmci_scores_", sex_name, "_sharp.rds")))
  
  # Save all the metrics 
  metrics_dt <- data.table(
    sex = sex_name, 
    method = "SHARP_PMCI", 
    seed = SEED, 
    n = nrow(df_sex), 
    n_train = length(tr), 
    n_test = length(te), 
    n_stable = length(stable_feats), 
    n_drug_coef = nrow(drug_dt), 
    scale_factor = SCALE_FACTOR, 
    cindex_test_age_only = c_age, 
    cindex_test_age_pmci_pts = c_pts, 
    cindex_test_age_pmci_lp = c_lp, 
    cindex_test_age_only_lo = boot_age$lo,
    cindex_test_age_only_hi = boot_age$hi,
    cindex_test_age_pmci_pts_lo = boot_pts$lo,
    cindex_test_age_pmci_pts_hi = boot_pts$hi,
    cindex_test_age_pmci_lp_lo = boot_lp$lo,
    cindex_test_age_pmci_lp_hi = boot_lp$hi,
    boot_B = B_BOOT,
    boot_level = CI_LEVEL,
    boot_stratify_event = BOOT_STRATIFY_EVENT
  )
  fwrite(metrics_dt, file.path(OUT_DIR, paste0("pmci_metrics_", sex_name, "_sharp.csv"))) 
  
  # Now save the models 
  models <- list(
    cox_refit_stable = cox_refit, 
    fit_age_only = fit_age, 
    fit_age_pmci_pts = fit_pts, 
    fit_age_pmci_lp = fit_lp, 
    stable_feats_raw = stable_feats_raw, 
    stable_feats = stable_feats
  )
  saveRDS(models, file.path(OUT_DIR, paste0("pmci_models_", sex_name, "_sharp.rds"))) 
  
  invisible(list(scores = scores_dt, weights = weights_out, metrics = metrics_dt, models = models)) 
}


# Main Data Processing 
cat("Loading data: ", DATA_RDS, "\n", sep = "") 
df <- as.data.table(readRDS(DATA_RDS)) 
stopifnot(all(c("sex", "time", "event", "age_at_entry") %in% names(df))) 

df_f <- df[sex == "Female"] 
df_m <- df[sex == "Male"] 
cat("N Female: ", nrow(df_f), " | events: ", sum(as.integer(df_f$event)), "\n", sep = "") 
cat("N Male: ", nrow(df_m), " | events: ", sum(as.integer(df_m$event)), "\n", sep = "") 

res_f <- derive_pmci_one_sex(df_f, "Female") 
res_m <- derive_pmci_one_sex(df_m, "Male") 

cat("\nAll done. Outputs to: ", OUT_DIR, "\n", sep = "")






