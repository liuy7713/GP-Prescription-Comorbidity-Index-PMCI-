# 04_stability_selection.R 
# Purpose of this file is to move on from the standard Elastic Net Cox model in 03_elastic_net_cox_gp.R 
# Here implementing a form of stability selection (coded, not implementing SHARP yet) 

# For each sex, I repeatedly subsample 50% of the training set B times 
# On each subsample, fit a CV tuned elastic net cox model using glmnet and record 
# which penalised predictors have non zero coefficients at the chosen lambda. 
# Obtain each predictors selection proportion and predictors with proportion 
# greater equal to PI_THRESH are considered stable. 
# Then refit an unpenalised Cox model on the stable predictors with age to get 
# coefficients used to build the final index. 


library(data.table) 
library(glmnet) 
library(survival) 
library(Matrix) 

# Configurations 
DATA_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"
EN_DIR   <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/en_cox"
OUT_DIR  <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/stab_naive"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE) 

SEED <- 20260117 

# Stability Selection Settings 
B <- 100  # Previously set to 50 for speed 
SUBSAMPLE_FRAC <- 0.5  # Takes 50% of the subsample every time 
N_FOLDS <- 10  # Previously set to 5 for speed 
PI_THRESH <- 0.8 

# Choice of lambda inside each sample (lambda.min OR lambda.1se) 
LAMBDA_CHOICE <- "lambda.1se" 

# Choice of Alpha 
# Either direct from the results of 03_elastic_net_cox_gp.R 
# Or a fixed alpha for replication e.g. 0.5 or 1.0 
ALPHA_MODE <- "from_03"  # Either "fixed" or "from_03" 
ALPHA_FIXED <- 0.5 

# C-Index function 
cindex_cox <- function(time, event, lp) { 
  cc <- survival::concordance(survival::Surv(time, event) ~ I(-lp)) 
  as.numeric(cc$concordance) 
}

# Use the same feature filtering rules from 03_elastic_net_cox_gp.R 
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

# Read the alphas for sex 
read_alpha_for_sex <- function(sex_name) { 
  if (ALPHA_MODE == "fixed") return(ALPHA_FIXED) 
  mfile <- file.path(EN_DIR, paste0("metrics_", sex_name, ".csv")) 
  if (!file.exists(mfile)) { 
    warning("metrics file not found for ", sex_name, ". Falling back to fixed alpha.") 
    return (ALPHA_FIXED) 
  }
  m <- fread(mfile) 
  a <- as.numeric(m$alpha[1]) 
  if (!is.finite(a)) { 
    warning("alpha in metrics not readable for ", sex_name, ". Falling back to fixed alpha.")
    return (ALPHA_FIXED) 
  }
  a 
}

load_split_idx <- function(sex_name, n) { 
  sfile <- file.path(EN_DIR, paste0("split_idx_", sex_name, ".rds")) 
  if (!file.exists(sfile)) { 
    stop("Missing split index file: ", sfile, "\n", "Patch and rerun 03_elastic_net_cox_gp.R to save split indices.") 
  }
  sp <- readRDS(sfile) 
  tr <- sp$tr 
  te <- sp$te 
  if (any(tr < 1 | tr > n) || any(te < 1 | te > n)) stop("split indices out of range for ", sex_name) 
  list(tr = tr, te = te, seed = sp$seed) 
}

# Stability Selection fitting function 
fit_stability_one_sex <- function(df_sex, sex_name) { 
  cat("\n============================\n")
  cat("Stability selection for:", sex_name, "\n")
  cat("============================\n")
  set.seed(SEED) 
  bnf_cols <- get_feature_cols(df_sex, sex_name) 
  alpha <- read_alpha_for_sex(sex_name) 
  cat("Using alpha = ", alpha, " | BNF features = ", length(bnf_cols), "\n")
  
  # Build matrix 
  X <- as.matrix(df_sex[, c("age_at_entry", bnf_cols), with = FALSE]) 
  y <- survival::Surv(df_sex$time, df_sex$event) 
  
  # penalty.factor = 0 for age unpenalised 
  pf <- c(0, rep(1, length(bnf_cols))) 
  
  # Reuse the same train split from 03_elastic_net_cox_gp.R 
  n <- nrow(df_sex) 
  sp <- load_split_idx(sex_name, n) 
  tr <- sp$tr 
  te <- sp$te 
  X_tr <- X[tr, , drop = FALSE] 
  y_tr <- y[tr] 
  X_te <- X[te, , drop = FALSE] 
  y_te <- y[te] 
  
  feats <- colnames(X_tr) 
  
  # Only count stability for the penalised predictors, always keep age, but do not count it 
  penalised_feats <- setdiff(feats, "age_at_entry") 
  sel_counts <- setNames(integer(length(penalised_feats)), penalised_feats) 
  
  # Stability section 
  for (b in seq_len(B)) { 
    ss <- sample(seq_len(nrow(X_tr)), size = floor(SUBSAMPLE_FRAC * nrow(X_tr)), replace = FALSE) 
    cvfit <- cv.glmnet(x = X_tr[ss, , drop = FALSE], y = y_tr[ss], family = "cox", alpha = alpha, nfolds = N_FOLDS, penalty.factor = pf, standardize = TRUE) 
    lam <- if (LAMBDA_CHOICE == "lambda.1se") cvfit$lambda.1se else cvfit$lambda.min 
    beta <- as.matrix(coef(cvfit, s = lam)) 
    
    selected <- rownames(beta)[beta[, 1] != 0] 
    selected <- intersect(selected, penalised_feats) 
    sel_counts[selected] <- sel_counts[selected] + 1L 
  }
  
  sel_prop <- sel_counts / B 
  sel_dt <- data.table(feature = names(sel_prop), sel_prop = as.numeric(sel_prop)) 
  setorder(sel_dt, -sel_prop) 
  # Stable features are selected more than PI_THRESH probability 
  stable_feats <- sel_dt[sel_prop >= PI_THRESH, feature]
  stable_feats <- stable_feats[!is.na(stable_feats)]
  cat("Stable features (>= ", PI_THRESH, "): ", length(stable_feats), "\n", sep = "") 
  
  # Tag to distinguish lambda results 
  tag <- if (LAMBDA_CHOICE == "lambda.1se") "l1se" else "lmin"
  
  # Save selection proportions and stable list 
  fwrite(sel_dt, file.path(OUT_DIR, paste0("selection_proportions_", sex_name, "_", tag, ".csv"))) 
  fwrite(data.table(feature = stable_feats), file.path(OUT_DIR, paste0("stable_variables_", sex_name, "_", tag, ".csv"))) 
  
  # Refit the unpenalised cox on stable set (with age) 
  rhs <- paste(c("age_at_entry", stable_feats), collapse = " + ") 
  fml <- as.formula(paste0("Surv(time, event) ~ ", rhs)) 
  
  dtr <- df_sex[tr] 
  dte <- df_sex[te] 
  
  # Function to fit an unpenalised Cox proportional hazards model 
  cox_refit <- coxph(fml, data = dtr, ties = "efron") 
  
  lp_te <- predict(cox_refit, newdata = dte, type = "lp") 
  c_te <- cindex_cox(dte$time, dte$event, lp_te) 
  
  # Save refit and coefficients 
  saveRDS(cox_refit, file.path(OUT_DIR, paste0("cox_refit_", sex_name, "_", tag, ".rds"))) 
  coefs <- data.table(feature = names(coef(cox_refit)), coef = as.numeric(coef(cox_refit))) 
  coefs[, abs_coef := abs(coef)] 
  setorder(coefs, -abs_coef) 
  coefs[, abs_coef := NULL] 
  fwrite(coefs, file.path(OUT_DIR, paste0("cox_refit_coef_", sex_name, "_", tag, ".csv"))) 
  
  # Metrics to store 
  fwrite(data.table(sex = sex_name, B = B, subsample_frac = SUBSAMPLE_FRAC, alpha = alpha, lambda_choice = LAMBDA_CHOICE, pi_thresh = PI_THRESH, n = n, n_train = length(tr), n_test = length(te), n_stable = length(stable_feats), cindex_test = c_te), file.path(OUT_DIR, paste0("metrics_", sex_name, "_", tag, ".csv"))) 
  cat("Test C-index (refit Cox) = ", round(c_te, 4), "\n") 
}


# Main Data Processing 
cat("Loading data: ", DATA_RDS, "\n") 
df <- as.data.table(readRDS(DATA_RDS)) 
stopifnot(all(c("sex", "time", "event", "age_at_entry") %in% names(df))) 
df[, event := as.integer(event)] 
df[, time := as.numeric(time)] 
df[, age_at_entry := as.numeric(age_at_entry)] 

df_f <- df[sex == "Female"] 
df_m <- df[sex == "Male"] 

cat("N Female: ", nrow(df_f), " | events: ", sum(df_f$event), "\n")
cat("N Male: ", nrow(df_m), " | events: ", sum(df_m$event), "\n") 

fit_stability_one_sex(df_f, "Female") 
fit_stability_one_sex(df_m, "Male") 

cat("\nAll done. Outputs in: \n", OUT_DIR, "\n", sep = "") 



