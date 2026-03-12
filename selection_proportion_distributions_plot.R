# selection_proportion_distributions_plot.R 

# Purpose of this file is to create a combined plot of the selection proportion 
# distributions of naive stability selection (both lambda_1se and lambda_min) 
# as well as SHARP calibrated results. 

# The finalised diagram can be found in the Appendix. 
# Calibration results combined with this diagram show that SHARP only considers 
# features as "stable" if they consistently appear in almost every subsample 
# in the current models. 


library(data.table)
library(ggplot2)

# Configurations 
ROOT <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp"
NAIVE_DIR <- file.path(ROOT, "stab_naive")
SHARP_DIR <- file.path(ROOT, "sharp")
OUT_DIR <- file.path(ROOT, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PDF <- file.path(OUT_DIR, "fig_selection_proportion_distributions.pdf")
OUT_PNG <- file.path(OUT_DIR, "fig_selection_proportion_distributions.png")

read_selprop <- function(path) {
  dt <- fread(path)
  prop_col <- intersect(names(dt), c("sel_prop", "selection_proportion", "pi", "prob", "proportion"))
  if (length(prop_col) == 0) {
    stop("Cannot find selection proportion column in: ", path,
         "\nColumns: ", paste(names(dt), collapse = ", "))
  }
  setnames(dt, prop_col[1], "sel_prop")
  dt[, sel_prop := as.numeric(sel_prop)]
  dt
}

find_naive_file <- function(sex, tag) {
  file.path(NAIVE_DIR, paste0("selection_proportions_", sex, "_", tag, ".csv"))
}

find_sharp_file <- function(sex) {
  file.path(SHARP_DIR, paste0("selection_proportions_", sex, "_sharp.csv"))
}

sexes <- c("Female", "Male")
all_dt <- rbindlist(lapply(sexes, function(sex) {
  
  f_lmin <- find_naive_file(sex, "lmin")
  dt_lmin <- read_selprop(f_lmin)[, .(sel_prop)]
  dt_lmin[, `:=`(sex = sex, method = "Naive (lambda-min)")]
  
  f_l1se <- find_naive_file(sex, "l1se")
  dt_l1se <- read_selprop(f_l1se)[, .(sel_prop)]
  dt_l1se[, `:=`(sex = sex, method = "Naive (lambda-1se)")]
  
  f_sharp <- find_sharp_file(sex)
  dt_sharp <- read_selprop(f_sharp)[, .(sel_prop)]
  dt_sharp[, `:=`(sex = sex, method = "SHARP (calibrated)")]
  
  rbindlist(list(dt_lmin, dt_l1se, dt_sharp), use.names = TRUE)
}), use.names = TRUE)

all_dt[, sex := factor(sex, levels = c("Female", "Male"))]
all_dt[, method := factor(method, levels = c("Naive (lambda-min)", "Naive (lambda-1se)", "SHARP (calibrated)"))]

# Using dashed lines for naive panels (they are fixed at arbitrary pi^* = 0.8) 
vline_naive <- CJ(
  sex = factor(c("Female", "Male"), levels = c("Female", "Male")),
  method = factor(c("Naive (lambda-min)", "Naive (lambda-1se)"),
                  levels = levels(all_dt$method))
)
vline_naive[, x := 0.8]

# Using dotted line for SHARP panels (these are calibrated results, here pi^* = 0.99) 
vline_sharp <- data.table(
  sex = factor(c("Female", "Male"), levels = c("Female", "Male")),
  method = factor(rep("SHARP (calibrated)", 2), levels = levels(all_dt$method)),
  x = c(0.99, 0.99)
)

p <- ggplot(all_dt, aes(x = sel_prop)) +
  geom_histogram(binwidth = 0.02, boundary = 0, closed = "left") +
  facet_grid(sex ~ method) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "Selection proportion",
    y = "Number of features",
    caption = "Dashed: fixed pi = 0.8 (naive). Dotted: SHARP calibrated pi* (sex-specific)." 
  ) +
  geom_vline(
    data = vline_naive, 
    aes(xintercept = x), 
    linetype = "dashed", 
    linewidth = 0.4, 
    inherit.aes = FALSE 
  ) + 
  geom_vline(
    data = vline_sharp, 
    aes(xintercept = x), 
    linetype = "dotted", 
    linewidth = 0.6, 
    inherit.aes = FALSE 
  ) + 
  theme_bw(base_size = 12) +
  theme(
    strip.text.y = element_text(face = "bold"),  # makes Female/Male more obvious
    strip.text.x = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

ggsave(OUT_PDF, p, width = 11, height = 6)
ggsave(OUT_PNG, p, width = 11, height = 6, dpi = 300)

cat("Saved:\n", OUT_PDF, "\n", OUT_PNG, "\n")



