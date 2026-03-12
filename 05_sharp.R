# 05_sharp.R 

# Purpose of this file here is to fully implement the SHARP package into stability selection 

# The idea here is to implement SHARP on an elastic net cox base learner validated 
# with B, pi grid and lambda window to obtain a stable set of prescription features 
# that I can use to construct a reliable (relatively high C-Index) comorbidity index 


library(data.table) 
library(survival) 
library(glmnet) 
library(sharp) 

# Configurations 
DATA_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"
EN_DIR   <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/en_cox"
OUT_DIR  <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/sharp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE) 

SEED <- 20260119 

# SHARP Specific Settings 
# This choice of fixed alpha value is a specific modelling choice which is 
# explained in the progress report 
ALPHA_FIXED <- 0.8 
# B is defined as the number of subsamples chosen 
# This is defined to be the same as 04_stability_selection.R for consistency 
B <- 100 
# Subsample Fraction (By convention defined to be 50%) 
SUBSAMPLE_FRAC <- 0.5 
# Instantiating the PI value grid 
PI_LIST <- seq(0.5, 0.99, by = 0.01) 

# C Index Evaluation function 
cindex_cox <- function(time, event, lp) { 
  cc <- survival::concordance(survival::Surv(time, event) ~ I(-lp)) 
  as.numeric(cc$concordance)
}

# Define the get_feature_cols function 
# Same feature filtering rules as used in 03_elastic_net_cox_gp.R 
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

# Loading the split index function 
# Same split index function as used in 04_stability_selection.R 
load_split_idx <- function(sex_name, n) { 
  sfile <- file.path(EN_DIR, paste0("split_idx_", sex_name, ".rds")) 
  if (!file.exists(sfile)) { 
    stop("Missing split index file: ", sfile, "\n")
  }
  sp <- readRDS(sfile) 
  tr <- sp$tr 
  te <- sp$te 
  if (any(tr < 1 | tr > n) || any (te < 1 | te > n)) stop("split indices out of range for ", sex_name) 
  list(tr = tr, te = te, seed = sp$seed)
}

# Defining the actual fitting function with SHARP wrapped 
fit_sharp_one_sex <- function(df_sex, sex_name) { 
  cat("\n============================\n")
  cat("SHARP for:", sex_name, "\n")
  cat("============================\n")
  set.seed(SEED) 
  
  bnf_cols <- get_feature_cols(df_sex, sex_name) 
  alpha <- ALPHA_FIXED 
  cat("Using alpha = ", alpha, " | BNF features = ", length(bnf_cols), "\n") 
  
  # Build X and y for Matrix construction 
  X <- as.matrix(df_sex[, c("age_at_entry", bnf_cols), with = FALSE]) 
  y <- Surv(df_sex$time, df_sex$event) 
  
  # Initialising penalty.factor 
  # Keeping age unpenalised (penalty.factor = 0) and bnf columns penalised (1) 
  pf <- c(0, rep(1, length(bnf_cols))) 
  
  # Reuse the same train test split from 03_elastic_net_cox_gp.R 
  n <- nrow(df_sex) 
  sp <- load_split_idx(sex_name, n) 
  tr <- sp$tr 
  te <- sp$te 
  cat("Split: train=", length(tr), " test=", length(te), "\n", sep = "")
  
  X_tr <- X[tr, , drop = FALSE] 
  y_tr <- y[tr] 
  dtr <- df_sex[tr] 
  dte <- df_sex[te] 
  
  # Build a lambda grid on the training set via glmnet 
  gfit <- glmnet(x = X_tr, y = y_tr, family = "cox", alpha = alpha, penalty.factor = pf, standardize = TRUE) 
  
  # Restricting lambda grid to a sensible range 
  # This is done by restricting to models which will provide a sensible number of features 
  beta <- as.matrix(coef(gfit)) 
  nnz <- colSums(beta[-1, , drop = FALSE] != 0) 
  keep <- which(nnz >= 0 & nnz <= 124) 
  cat("nnz range across lambda path: ", min(nnz), " to ", max(nnz), "\n", sep = "") 
  cat("Keeping ", length(keep), " lambdas \n", sep = "") 
  stopifnot(length(keep) > 0) 
  Lambda <- matrix(gfit$lambda[keep], ncol = 1) 
  colnames(Lambda) <- "lambda"
  
  # Implementing SHARP stability selection package here 
  # Here I am using the VariableSelection() function from SHARP 
  stab <- sharp::VariableSelection( 
    xdata = X_tr, 
    ydata = y_tr, 
    #Lambda = Lambda, 
    pi_list = PI_LIST, 
    K = B, 
    tau = SUBSAMPLE_FRAC, 
    seed = SEED, 
    family = "cox", 
    implementation = sharp::PenalisedRegression, 
    resampling = "subsampling", 
    n_cat = 3, 
    # These are used to DISABLE PFER error control 
    PFER_method = "MB", 
    PFER_thr = Inf, 
    FDP_thr = Inf, 
    #n_cores = 6, 
    penalty.factor = pf, 
    alpha = alpha 
    #standardize = TRUE 
  ) 
  cat("Calibrated (lambda*, pi*): \n")
  print(sharp::Argmax(stab)) 
  cat("Calibrated pi per lambda (stab$P) and stable count per lambda (stab$Q_s): \n")
  print(head(stab$P)) 
  print(head(stab$Q_s)) 
  
  # Determining the height of the selection proportions 
  sel_prop_mat <- sharp::SelectionProportions(stab) 
  cat("Max selection proportion: ", max(sel_prop_mat), "\n")
  
  # Saving calibration plot 
  pdf(file.path(OUT_DIR, paste0("calibration_", sex_name, "_sharp.pdf"))) 
  sharp::CalibrationPlot(stab) 
  dev.off() 
  saveRDS(stab, file.path(OUT_DIR, paste0("sharp_object_", sex_name, ".rds"))) 
  
  # Extracting selection proportions and the calibrated stable feature set 
  selprop <- sharp::SelectionProportions(stab) 
  if (is.null(dim(selprop))) { 
    sel_dt <- data.table(feature = names(selprop), sel_prop = as.numeric(selprop)) 
  } else { 
    sel_dt <- data.table(feature = colnames(selprop), sel_prop = apply(selprop, 2, max)) 
  }
  setorder(sel_dt, -sel_prop) 
  # Extract the stable variables at the calibrated tuple 
  sv <- sharp::SelectedVariables(stab) 
  if (is.null(dim(sv))) { 
    sel_idx <- which(as.integer(sv) == 1L) 
    stable_feats <- names(sv)[sel_idx] 
  } else { 
    if (nrow(sv) == 1) { 
      sel_idx <- which(as.integer(sv[1, ]) == 1L) 
      stable_feats <- colnames(sv)[sel_idx] 
    } else if (ncol(sv) == 1) { 
      sel_idx <- which(as.integer(sv[, 1]) == 1L) 
      stable_feats <- rownames(sv)[sel_idx] 
    } else { 
      sel_idx <- which(as.integer(sv[1, ]) == 1L)
      stable_feats <- colnames(sv)[sel_idx] 
    }
  }
  stable_feats <- setdiff(stable_feats, "age_at_entry") 
  stable_feats <- intersect(stable_feats, colnames(df_sex)) 
  # Just to confirm that all feature names exists 
  stopifnot(all(stable_feats %in% colnames(X_tr))) 
  cat("Stable Features: ", paste(stable_feats, collapse = ", "), "\n", sep = "") 
  
  rhs <- if (length(stable_feats) == 0) "age_at_entry" else paste(c("age_at_entry", stable_feats), collapse = " + ") 
  fml <- as.formula(paste0("Surv(time, event) ~ ", rhs)) 
  
  fwrite(sel_dt, file.path(OUT_DIR, paste0("selection_proportions_", sex_name, "_sharp.csv"))) 
  fwrite(data.table(feature = stable_feats), file.path(OUT_DIR, paste0("stable_variables_", sex_name, "_sharp.csv"))) 
  cat("Stable features count: ", length(stable_feats), "\n", sep = "") 
  
  # Refit unpenalised Cox on the stable set with age 
  cox_refit <- coxph(fml, data = dtr, ties = "efron") 
  lp_te <- predict(cox_refit, newdata = dte, type = "lp") 
  c_te <- cindex_cox(dte$time, dte$event, lp_te) 
  saveRDS(cox_refit, file.path(OUT_DIR, paste0("cox_refit_", sex_name, "_sharp.rds"))) 
  fwrite(data.table(feature = names(coef(cox_refit)), coef = as.numeric(coef(cox_refit))), file.path(OUT_DIR, paste0("cox_refit_coef_", sex_name, "_sharp.csv")))
  fwrite(data.table(sex = sex_name, method = "SHARP", B = B, subsample_frac = SUBSAMPLE_FRAC, alpha = alpha, n = n, n_train = length(tr), n_test = length(te), n_stable = length(stable_feats), cindex_test = c_te), file.path(OUT_DIR, paste0("metrics_", sex_name, "_sharp.csv"))) 
  cat("Test C-index = ", round(c_te, 4), "\n") 
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

fit_sharp_one_sex(df_f, "Female") 
fit_sharp_one_sex(df_m, "Male") 

cat("\nAll done. Outputs in: \n", OUT_DIR, "\n", sep = "") 


