# sharp_calibration_heatmap.R 

# Purpose of this file is to utilise the internal SHARP function CalibrationPlot
# to create two separate SHARP calibration heatmap plots. 

# The two output plots can be found in the Appendix in the paper. 

# The two figures illustrate the optimal calibration point for SHARP and the 
# lambda^* and pi^* values. 


library(sharp)

# Configurations 
ROOT <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp"
SHARP_DIR <- file.path(ROOT, "sharp")
OUT_DIR   <- file.path(ROOT, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

plot_one <- function(sex) {
  stab_path <- file.path(SHARP_DIR, paste0("sharp_object_", sex, ".rds"))
  if (!file.exists(stab_path)) stop("Missing: ", stab_path)
  
  stab <- readRDS(stab_path)
  
  a <- sharp::Argmax(stab)
  lam_star <- as.numeric(a[1, "lambda"])
  pi_star  <- as.numeric(a[1, "pi"])
  
  lambda_vals <- as.numeric(stab$Lambda[, 1])
  nL <- length(lambda_vals)
  
  out_pdf <- file.path(OUT_DIR, paste0("fig_sharp_calibration_", sex, ".pdf"))
  pdf(out_pdf, width = 10.5, height = 6.5, useDingbats = FALSE)
  
  op <- par(no.readonly = TRUE)
  on.exit({par(op); dev.off()}, add = TRUE)
  
  par(mar = c(5, 7, 6, 8))
  
  par(xaxt = "n")
  print(sharp::Argmax(stab)) 
  sharp::CalibrationPlot(stab)
  
  
  par(xaxt = "s")
  
  mtext(expression(lambda), side = 1)  # x-axis label a bit lower
  
  # title + annotation (kept compact)
  title(main = paste0("SHARP calibration plot (", sex, ")"), font.main = 2, cex.main = 1.5)
  mtext(sprintf("Calibrated point:  lambda* = %.5f,  pi* = %.2f", lam_star, pi_star),
        side = 3, line = 0.6, adj = 0, cex = 0.95)
  
  cat("Saved:", out_pdf, "\n")
}

plot_one("Female")
plot_one("Male")

