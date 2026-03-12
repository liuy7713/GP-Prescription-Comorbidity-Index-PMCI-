# km_pmci_pts.R 

# The purpose of this file is to create the two Kaplan-Meier plots for the PMCI 
# for both sexes. 

# These two diagrams examine the differences between user with PMCI < 0 and 
# PMCI >= 0. 

# Both figures are placed in the Appendix in the paper. 


library(data.table) 
library(survival) 

# Configurations 
ROOT <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp"
PMCI_DIR <- file.path(ROOT, "sharp", "pmci_sharp") 
OUT_DIR  <- file.path(ROOT, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TIME_COL <- "time" 
EVENT_COL <- "event" 
SCORE_COL <- "pmci_pts" 

TIME_IN_DAYS <- 1 

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-4) return("< 1e-4")
  formatC(p, format = "g", digits = 3)
}

# Kaplan-Meier Plotting Function 
plot_km_pmci_pts_one_sex <- function(sex) {
  rds_path <- file.path(PMCI_DIR, paste0("pmci_scores_", sex, "_sharp.rds"))
  if (!file.exists(rds_path)) stop("Missing: ", rds_path)
  
  d <- as.data.table(readRDS(rds_path))
  
  # Basic checks
  req <- c(TIME_COL, EVENT_COL, SCORE_COL)
  miss <- setdiff(req, names(d))
  if (length(miss) > 0) stop("Missing columns in ", rds_path, ": ", paste(miss, collapse = ", "))
  
  d <- d[!is.na(get(TIME_COL)) & !is.na(get(EVENT_COL)) & !is.na(get(SCORE_COL))]
  d[, (EVENT_COL) := as.integer(get(EVENT_COL))]
  d[, (TIME_COL)  := as.numeric(get(TIME_COL)) * TIME_IN_DAYS]
  d[, (SCORE_COL) := as.numeric(get(SCORE_COL))]
  
  # PMCI group: >0 vs <=0
  d[, pmci_group := ifelse(get(SCORE_COL) > 0, "PMCI > 0", "PMCI <= 0")]
  d[, pmci_group := factor(pmci_group, levels = c("PMCI <= 0", "PMCI > 0"))]
  
  # Survival fit (time in DAYS)
  fml <- as.formula(paste0("Surv(", TIME_COL, ", ", EVENT_COL, ") ~ pmci_group"))
  fit <- survfit(fml, data = d)
  
  # Log-rank test 
  lr <- survdiff(fml, data = d)
  p_lr <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  # Output paths
  out_pdf <- file.path(OUT_DIR, paste0("fig_km_pmci_pts_", sex, ".pdf"))
  out_png <- file.path(OUT_DIR, paste0("fig_km_pmci_pts_", sex, ".png"))
  
  # Try survminer if available 
  has_survminer <- requireNamespace("survminer", quietly = TRUE)
  
  if (has_survminer) {
    # ggsurvplot expects time in whatever units you give it; we give DAYS but label years.
    library(ggplot2)
    g <- survminer::ggsurvplot(
      fit, data = d,
      risk.table = TRUE,
      pval = paste0("Log-rank p = ", fmt_p(p_lr)),
      conf.int = FALSE,
      censor = TRUE,
      legend.title = NULL,
      legend.labs = paste0(levels(d$pmci_group), " (n=", as.integer(table(d$pmci_group)[levels(d$pmci_group)]), ")"),
      xlab = "Time since entry (years)",
      ylab = "Survival probability",
      break.time.by = 365.25 * 2,  # ticks every 2 years (in DAYS)
      ggtheme = theme_bw(base_size = 12)
    )
    
    # Convert x-axis labels from days -> years
    g$plot <- g$plot + scale_x_continuous(
      labels = function(x) sprintf("%.0f", x / 365.25)
    )
    g$table <- g$table + scale_x_continuous(
      labels = function(x) sprintf("%.0f", x / 365.25)
    )
    
    g$plot <- g$plot + ggtitle(paste0("Kaplan Meier curves by PMCI (", sex, ")")) +
      theme(plot.title = element_text(face = "bold"))
    
    # Save output 
    grDevices::pdf(out_pdf, width = 8.2, height = 7.2, useDingbats = FALSE)
    print(g)
    dev.off()
    
    grDevices::png(out_png, width = 1400, height = 1200, res = 180)
    print(g)
    dev.off()
    
  } else {
    # Base survival plotting 
    grDevices::pdf(out_pdf, width = 7.8, height = 6.2, useDingbats = FALSE)
    op <- par(no.readonly = TRUE)
    on.exit({par(op); dev.off()}, add = TRUE)
    
    par(mar = c(5, 5, 4.2, 2))
    cols <- c("PMCI <= 0" = "black", "PMCI > 0" = "red3")
    
    # Plot in DAYS but label YEARS on axis
    plot(
      fit,
      col = cols[levels(d$pmci_group)],
      lwd = 2,
      xlab = "Time since entry (years)",
      ylab = "Survival probability",
      mark.time = TRUE,
      xaxt = "n"
    )
    
    # custom x-axis: show years
    xmax <- max(d[[TIME_COL]], na.rm = TRUE)
    yr_breaks <- pretty(c(0, xmax / 365.25))
    axis(1, at = yr_breaks * 365.25, labels = sprintf("%.0f", yr_breaks))
    
    title(main = paste0("Kaplan Meier curves by PMCI (", sex, ")"), font.main = 2)
    
    legend(
      "bottomleft",
      legend = paste0(levels(d$pmci_group), " (n=", as.integer(table(d$pmci_group)[levels(d$pmci_group)]), ")"),
      col = cols[levels(d$pmci_group)],
      lwd = 2,
      bty = "n"
    )
    
    mtext(paste0("Log-rank p = ", fmt_p(p_lr)), side = 3, line = 0.3, adj = 0, cex = 0.95)
    dev.off()
    
    # PNG
    grDevices::png(out_png, width = 1200, height = 900, res = 160)
    par(mar = c(5, 5, 4.2, 2))
    plot(
      fit,
      col = cols[levels(d$pmci_group)],
      lwd = 2,
      xlab = "Time since entry (years)",
      ylab = "Survival probability",
      mark.time = TRUE,
      xaxt = "n"
    )
    axis(1, at = yr_breaks * 365.25, labels = sprintf("%.0f", yr_breaks))
    title(main = paste0("Kaplan Meier curves by PMCI (", sex, ")"), font.main = 2)
    legend(
      "bottomleft",
      legend = paste0(levels(d$pmci_group), " (n=", as.integer(table(d$pmci_group)[levels(d$pmci_group)]), ")"),
      col = cols[levels(d$pmci_group)],
      lwd = 2,
      bty = "n"
    )
    mtext(paste0("Log-rank p = ", fmt_p(p_lr)), side = 3, line = 0.3, adj = 0, cex = 0.95)
    dev.off()
  }
  
  cat("Saved:\n", out_pdf, "\n", out_png, "\n")
  invisible(list(pdf = out_pdf, png = out_png, logrank_p = p_lr))
}


plot_km_pmci_pts_one_sex("Female") 
plot_km_pmci_pts_one_sex("Male") 
