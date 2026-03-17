# 03_elastic_net_cox_gp.R 
# Purpose to construct a one iteration elastic net cox model 
# Simple Elastic Net Cox benchmark model with sex stratified with GP prescription features 

# Input file: data_full_gp_cov.rds (derived from 02_merge_baseline_sex_age.R) 
# Output files: 
# cv_fit_<sex>.rds 
# alpha_grid_summary_<sex>.csv 
# coef_<sex>_lambda_min.csv 
# coef_<sex>_lambda_1se.csv 
# metrics_<sex>.csv 

# Age is unpenalised here by using penalty.factor = 0 for age_at_entry 
# Predictors here is BNF paragraph indicators in binary and age 

library(data.table) 
library(survival) 
library(glmnet) 
library(Matrix) 

# Configurations 
DATA_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"
OUT_DIR  <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/en_cox"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE) 

SEED <- 20260108 
TRAIN_FRAC <- 0.8 
N_FOLDS <- 10 
# alpha grid configuration 
ALPHAS <- c(0.2, 0.5, 0.8, 1.0) 

# Harrell's C-index function 
cindex_cox <- function(time, event, lp) { 
  cc <- survival::concordance(survival::Surv(time, event) ~ I(-lp)) 
  as.numeric(cc$concordance) 
}

# Extract all non zero coefficients into a sorted data.table 
extract_nz_coefs <- function(cvfit, s) { 
  beta <- as.matrix(coef(cvfit, s = s)) 
  nz <- which(beta != 0) 
  out <- data.table(feature = rownames(beta)[nz], coef = as.numeric(beta[nz, 1])) 
  out[, abs_coef := abs(coef)] 
  setorder(out, -abs_coef) 
  out[, abs_coef := NULL] 
  out 
}

# Fitting function for a given sex 
fit_one_sex <- function(df_sex, sex_name, alphas, out_dir) { 
  cat("\n============================\n")
  cat("Fitting Elastic Net Cox for:", sex_name, "\n")
  cat("============================\n")
  set.seed(SEED) 
  stopifnot(all(c("time", "event", "age_at_entry") %in% names(df_sex))) 
  # Feature columns are all bnf_* + age_at_entry 
  bnf_cols <- grep("^bnf_", names(df_sex), value = TRUE) 
  # Drop all columns with zero variance or are all zero in this sex 
  # This is an attempt to remove female only drugs in male database 
  tmp <- df_sex[, ..bnf_cols] 
  tmp <- as.data.frame(lapply(tmp, as.numeric)) 
  nonzero <- colSums(tmp, na.rm = TRUE) > 0 
  bnf_cols <- bnf_cols[nonzero] 
  if (length(bnf_cols) == 0) stop("No bnf_ columns found in dataset. ")
  # Drop specified female only drugs in male database to confirm  
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
    if (length(drop_male) > 0) { 
      cat("Dropping female specific categories for males: ", length(drop_male), "\n")
      cat("-", paste(drop_male, collapse = "\n-"), "\n")
      bnf_cols <- setdiff(bnf_cols, drop_male) 
    }
  }
  
  cat("Using", length(bnf_cols), "BNF features for", sex_name, "\n")
  
  # Building the matrix 
  X <- as.matrix(df_sex[, c("age_at_entry", bnf_cols), with = FALSE]) 
  y <- survival::Surv(df_sex$time, df_sex$event) 
  # penalty.factor = 0 for age and 1 for all BNF features 
  pf <- c(0, rep(1, length(bnf_cols))) 
  # Train Test split 
  n <- nrow(df_sex) 
  idx <- sample.int(n) 
  n_train <- floor(TRAIN_FRAC * n) 
  tr <- idx[1:n_train] 
  te <- idx[(n_train + 1):n] 
  
  X_tr <- X[tr, , drop = FALSE] 
  y_tr <- y[tr] 
  X_te <- X[te, , drop = FALSE] 
  y_te <- y[te] 
  
  set.seed(SEED) 
  foldid <- sample(rep(1:N_FOLDS, length.out = length(tr))) 
  saveRDS(list(tr = tr, te = te, foldid = foldid, seed = SEED), file.path(out_dir, paste0("split_idx_", sex_name, ".rds")))
  
  # Track alpha grid results 
  alpha_rows <- list() 
  
  # Best is chosen by CV partial likelihood deviance 
  best <- list(
    alpha = NA_real_,
    cvfit = NULL,
    lambda_min = NA_real_,
    lambda_1se = NA_real_,
    cv_score = Inf,
    cindex_train_min = NA_real_,
    cindex_test_min  = NA_real_,
    cindex_train_1se = NA_real_,
    cindex_test_1se  = NA_real_,
    nnz_min = NA_integer_,
    nnz_1se = NA_integer_
  )
  
  for (a in alphas) { 
    cat("CV for alpha = ", a, "\n")
    cvfit <- cv.glmnet(x = X_tr, y = y_tr, family = "cox", alpha = a, foldid = foldid, penalty.factor = pf, standardize = TRUE) 
    lam_min <- cvfit$lambda.min 
    lam_1se <- cvfit$lambda.1se 
    # lambda.min processing 
    lp_tr_min <- as.numeric(predict(cvfit, newx = X_tr, s = lam_min, type = "link")) 
    lp_te_min <- as.numeric(predict(cvfit, newx = X_te, s = lam_min, type = "link")) 
    c_tr_min <- cindex_cox(df_sex$time[tr], df_sex$event[tr], lp_tr_min) 
    c_te_min <- cindex_cox(df_sex$time[te], df_sex$event[te], lp_te_min) 
    nnz_min <- sum(as.matrix(coef(cvfit, s = lam_min)) != 0)
    # lambda.1se processing 
    lp_tr_1se <- as.numeric(predict(cvfit, newx = X_tr, s = lam_1se, type = "link")) 
    lp_te_1se <- as.numeric(predict(cvfit, newx = X_te, s = lam_1se, type = "link")) 
    c_tr_1se <- cindex_cox(df_sex$time[tr], df_sex$event[tr], lp_tr_1se) 
    c_te_1se <- cindex_cox(df_sex$time[te], df_sex$event[te], lp_te_1se) 
    nnz_1se <- sum(as.matrix(coef(cvfit, s = lam_1se)) != 0) 
    
    alpha_rows[[length(alpha_rows) + 1]] <- data.table(sex = sex_name, alpha = a, lambda_min = lam_min, cindex_train_min = c_tr_min, cindex_test_min = c_te_min, nnz_min = nnz_min, lambda_1se = lam_1se, cindex_train_1se = c_tr_1se, cindex_test_1se = c_te_1se, cv_score = min(cvfit$cvm), nnz_1se = nnz_1se, n = n, n_train = length(tr), n_test = length(te), events = sum(df_sex$event)) 
    
    cv_score <- min(cvfit$cvm)
    
    if (cv_score < best$cv_score - 1e-12 || (abs(cv_score - best$cv_score) <= 1e-12 && a > best$alpha)) {
      best$alpha <- a
      best$cvfit <- cvfit
      best$lambda_min <- lam_min
      best$lambda_1se <- lam_1se
      best$cv_score <- cv_score
      
      # store metrics for the chosen alpha 
      best$cindex_train_min <- c_tr_min
      best$cindex_test_min  <- c_te_min
      best$cindex_train_1se <- c_tr_1se
      best$cindex_test_1se  <- c_te_1se
      best$nnz_min <- nnz_min
      best$nnz_1se <- nnz_1se
    }
  }
  
  alpha_summary <- rbindlist(alpha_rows)
  fwrite(alpha_summary, file.path(out_dir, paste0("alpha_grid_summary_", sex_name, ".csv"))) 
  
  cat(
    "Best alpha for ", sex_name,
    " = ", best$alpha,
    " | CV score (min cvm) = ", signif(best$cv_score, 5),
    " | test C-index (lambda.min) = ", round(best$cindex_test_min, 4),
    " | lambda.min = ", signif(best$lambda_min, 4),
    " | lambda.1se = ", signif(best$lambda_1se, 4),
    "\n", sep = ""
  )
  
  # Save best cvfit 
  saveRDS(best$cvfit, file.path(out_dir, paste0("cv_fit_", sex_name, ".rds"))) 
  
  # Coefficients at lambda.min and lambda.1se 
  coef_min <- extract_nz_coefs(best$cvfit, best$lambda_min) 
  coef_1se <- extract_nz_coefs(best$cvfit, best$lambda_1se) 
  
  fwrite(coef_min, file.path(out_dir, paste0("coef_", sex_name, "_lambda_min.csv"))) 
  fwrite(coef_1se, file.path(out_dir, paste0("coef_", sex_name, "_lambda_1se.csv"))) 
  
  # Metrics for both lambda choices on the same train test split 
  lp_tr_1se_best <- as.numeric(predict(best$cvfit, newx = X_tr, s = best$lambda_1se, type = "link")) 
  lp_te_1se_best <- as.numeric(predict(best$cvfit, newx = X_te, s = best$lambda_1se, type = "link")) 
  c_tr_1se_best <- cindex_cox(df_sex$time[tr], df_sex$event[tr], lp_tr_1se_best) 
  c_te_1se_best <- cindex_cox(df_sex$time[te], df_sex$event[te], lp_te_1se_best) 
  
  metrics <- data.table(sex = sex_name, n = n, events = sum(df_sex$event), alpha = best$alpha, lambda_min = best$lambda_min, cindex_train_min = best$cindex_train_min, cindex_test_min = best$cindex_test_min, nnz_min = nrow(coef_min), lambda_1se = best$lambda_1se, cindex_train_1se = c_tr_1se_best, cindex_test_1se = c_te_1se_best, nnz_1se = nrow(coef_1se)) 
  fwrite(metrics, file.path(out_dir, paste0("metrics_", sex_name, ".csv"))) 
  
  invisible(list(alpha_summary = alpha_summary, metrics = metrics, coef_min = coef_min, coef_1se = coef_1se)) 
}


# Main Code Processing 
cat("Loading data: ", DATA_RDS, "\n")
df <- readRDS(DATA_RDS) 
df <- as.data.table(df) 

stopifnot(all(c("sex", "time", "event", "age_at_entry") %in% names(df))) 

# Ensures that all types are correct 
df[, event := as.integer(event)] 
df[, time := as.numeric(time)] 
df[, age_at_entry := as.numeric(age_at_entry)]

# Split by sex 
df_f <- df[sex == "Female"] 
df_m <- df[sex == "Male"] 

cat("N Female: ", nrow(df_f), " | events: ", sum(df_f$event), "\n")
cat("N Male: ", nrow(df_m), " | events: ", sum(df_m$event), "\n")

# Fitting by sex 
res_f <- fit_one_sex(df_f, "Female", ALPHAS, OUT_DIR) 
res_m <- fit_one_sex(df_m, "Male", ALPHAS, OUT_DIR) 

cat("\nAll Done. Outputs in: \n", OUT_DIR, "\n", sep = "") 