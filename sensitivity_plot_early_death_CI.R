# sensitivity_plot_early_death_CI.R 

# This script is a replacement for the script sensitivity_plot_early_death.R. 
# They both perform the same procedures, but this file includes the ability to 
# produce a 95% confidence interval via bootstrapping and has included this 
# result in the form of error bars on the sensitivity plot. 

# The resulting diagram can be found in the Appendix in the paper. 


library(data.table)
library(survival)
library(ggplot2)

# Configurations 
ROOT <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp"
DATA_RDS  <- file.path(ROOT, "data_full_gp_cov.rds")
EN_DIR    <- file.path(ROOT, "en_cox")
SHARP_DIR <- file.path(ROOT, "sharp")
OUT_DIR   <- file.path(ROOT, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Output Paths 
OUT_PDF <- file.path(OUT_DIR, "fig_sensitivity_early_death_exclusion_ci.pdf")
OUT_PNG <- file.path(OUT_DIR, "fig_sensitivity_early_death_exclusion_ci.png")

set.seed(20260211)
windows_years <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
B <- 500   # bootstrap replicates 

TIME_COL  <- "time"
EVENT_COL <- "event"

# Harrell's C-Index Function (same as all previous implementations) 
cindex_lp <- function(time, event, lp) {
  cc <- survival::concordance(survival::Surv(time, event) ~ I(-lp))
  as.numeric(cc$concordance)
}

# Bootstrap C-index for a fixed restricted test set
boot_cindex <- function(dte, lp, B = 500) {
  n <- nrow(dte)
  if (n < 50) return(list(mean = NA_real_, lo = NA_real_, hi = NA_real_))
  idx_mat <- replicate(B, sample.int(n, n, replace = TRUE))
  vals <- apply(idx_mat, 2, function(ii) cindex_lp(dte[[TIME_COL]][ii], dte[[EVENT_COL]][ii], lp[ii]))
  list(mean = mean(vals, na.rm = TRUE),
       lo   = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
       hi   = as.numeric(quantile(vals, 0.975, na.rm = TRUE)))
}

# Loading Train Test Split (same as all previous implementations) 
load_split_idx <- function(sex_name) {
  sfile <- file.path(EN_DIR, paste0("split_idx_", sex_name, ".rds"))
  if (!file.exists(sfile)) stop("Missing split index file: ", sfile)
  readRDS(sfile)
}

# Loading the SHARP refit model 
load_sharp_refit <- function(sex_name) {
  rfile <- file.path(SHARP_DIR, paste0("cox_refit_", sex_name, "_sharp.rds"))
  if (!file.exists(rfile)) stop("Missing SHARP refit Cox: ", rfile)
  readRDS(rfile)
}


# Loading Dataset 
df <- as.data.table(readRDS(DATA_RDS))
df[, (EVENT_COL) := as.integer(get(EVENT_COL))]
df[, (TIME_COL)  := as.numeric(get(TIME_COL))]
stopifnot(all(c("sex", TIME_COL, EVENT_COL) %in% names(df)))

# Bootstrapping Function for one sex 
run_one_sex <- function(sex_name) {
  dsex <- df[sex == sex_name]
  sp <- load_split_idx(sex_name)
  te <- sp$te
  dte_full <- dsex[te]
  
  fit <- load_sharp_refit(sex_name)
  
  # compute LP once on the full test set 
  lp_full <- as.numeric(predict(fit, newdata = dte_full, type = "lp"))
  
  out <- rbindlist(lapply(windows_years, function(wy) {
    cut_t <- wy * 365.25   # Used as time is in days 
    keep <- which(dte_full[[TIME_COL]] >= cut_t)
    
    dte <- dte_full[keep]
    lp  <- lp_full[keep]
    
    bci <- boot_cindex(dte, lp, B = B)
    
    data.table(
      years = wy,
      Sex = sex_name,
      n_test = nrow(dte),
      events = sum(dte[[EVENT_COL]]),
      cindex = cindex_lp(dte[[TIME_COL]], dte[[EVENT_COL]], lp),
      ci_lo = bci$lo,
      ci_hi = bci$hi
    )
  }))
  
  out
}

res <- rbindlist(list(
  run_one_sex("Female"),
  run_one_sex("Male")
))

# Plotting Functions 
res[, Sex := factor(Sex, levels = c("Female", "Male"))]

p <- ggplot(res, aes(x = years, y = cindex, colour = Sex, group = Sex)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.08, linewidth = 0.6) +
  scale_x_continuous(breaks = windows_years, name = "Exclusion window for early deaths (years)") +
  scale_y_continuous(limits = c(0.6, 1.0), name = "C-index") +
  theme_bw(base_size = 12) +
  theme(legend.title = element_blank())

ggsave(OUT_PDF, p, width = 6.8, height = 4.6, useDingbats = FALSE)
ggsave(OUT_PNG, p, width = 6.8, height = 4.6, dpi = 300)

fwrite(res, file.path(OUT_DIR, "sensitivity_early_death_exclusion_bootstrap_ci.csv"))
cat("Saved:\n", OUT_PDF, "\n", OUT_PNG, "\n",
    file.path(OUT_DIR, "sensitivity_early_death_exclusion_bootstrap_ci.csv"), "\n")

