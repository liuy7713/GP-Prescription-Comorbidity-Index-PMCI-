# upset_plot.R

# The script creates the UpSet plot comparing feature selections across SHARP 
# and RSF for both sexes. 

# The resulting UpSet plot can be found in the Appendix in the paper. 

library(data.table)
library(UpSetR)

# Configurations 
BASE_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp"
SHARP_DIR <- file.path(BASE_DIR, "sharp")
RSF_DIR   <- file.path(BASE_DIR, "rsf", "rsf_1683030")
FIG_DIR   <- file.path(BASE_DIR, "figures")

# Number of top RSF features to retain (Does not include age at entry) 
RSF_TOP_N <- 30

# Load the stable features from SHARP 
sharp_f <- fread(file.path(SHARP_DIR, "stable_variables_Female_sharp.csv"))$feature
sharp_m <- fread(file.path(SHARP_DIR, "stable_variables_Male_sharp.csv"))$feature

cat("SHARP Female:", length(sharp_f), "features\n")
cat("SHARP Male:  ", length(sharp_m), "features\n")

# Load the top 30 RSF features by permutation importance 
rsf_vimp_f <- fread(file.path(RSF_DIR, "rsf_vimp_Female_tuned_CI.csv"))
rsf_vimp_m <- fread(file.path(RSF_DIR, "rsf_vimp_Male_tuned_CI.csv"))

# Take top N, exclude age_at_entry (it is not a prescription feature)
rsf_f <- rsf_vimp_f[feature != "age_at_entry"][order(-vimp)][1:RSF_TOP_N]$feature
rsf_m <- rsf_vimp_m[feature != "age_at_entry"][order(-vimp)][1:RSF_TOP_N]$feature

cat("RSF Female:  ", length(rsf_f), "features\n")
cat("RSF Male:    ", length(rsf_m), "features\n")

# Build the binary membership matrix 
# Union of all features across the four sets
all_feats <- sort(unique(c(sharp_f, sharp_m, rsf_f, rsf_m)))
cat("\nTotal unique features across all sets:", length(all_feats), "\n")

# Clean feature names for display 
clean_name <- function(x) gsub("_", " ", sub("^bnf_", "", x))

membership <- data.frame(
  Feature         = clean_name(all_feats),
  `SHARP Female`  = as.integer(all_feats %in% sharp_f),
  `SHARP Male`    = as.integer(all_feats %in% sharp_m),
  `RSF Female`    = as.integer(all_feats %in% rsf_f),
  `RSF Male`      = as.integer(all_feats %in% rsf_m),
  check.names     = FALSE
)

# Print summary
cat("\nSet sizes:\n")
cat("  SHARP Female:", sum(membership$`SHARP Female`), "\n")
cat("  SHARP Male:  ", sum(membership$`SHARP Male`), "\n")
cat("  RSF Female:  ", sum(membership$`RSF Female`), "\n")
cat("  RSF Male:    ", sum(membership$`RSF Male`), "\n")

# Save membership table for reference
fwrite(membership, file.path(FIG_DIR, "upset_membership_table.csv"))

# UpSet Plotting 
# PDF 
pdf(file.path(FIG_DIR, "fig_upset_sharp_rsf.pdf"), width = 12, height = 7)

print(upset(
  membership,
  sets = c("RSF Male", "RSF Female", "SHARP Male", "SHARP Female"),
  keep.order = TRUE,
  order.by = "freq",
  nintersects = 30,
  sets.bar.color = c("blue", "lightblue", "orange", "red"),
  main.bar.color = "grey30",
  matrix.color = "grey30",
  point.size = 3,
  line.size = 1,
  text.scale = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.1),
  mb.ratio = c(0.6, 0.4),
  set_size.show = TRUE,
  set_size.scale_max = max(
    sum(membership$`SHARP Female`),
    sum(membership$`SHARP Male`),
    sum(membership$`RSF Female`),
    sum(membership$`RSF Male`)
  ) + 5
))

dev.off()

# PNG 
png(file.path(FIG_DIR, "fig_upset_sharp_rsf.png"), width = 2400, height = 1400, res = 200)

print(upset(
  membership,
  sets = c("RSF Male", "RSF Female", "SHARP Male", "SHARP Female"),
  keep.order = TRUE,
  order.by = "freq",
  nintersects = 30,
  sets.bar.color = c("blue", "lightblue", "orange", "red"),
  main.bar.color = "grey30",
  matrix.color = "grey30",
  point.size = 3,
  line.size = 1,
  text.scale = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.1),
  mb.ratio = c(0.6, 0.4),
  set_size.show = TRUE,
  set_size.scale_max = max(
    sum(membership$`SHARP Female`),
    sum(membership$`SHARP Male`),
    sum(membership$`RSF Female`),
    sum(membership$`RSF Male`)
  ) + 5
))

dev.off()

cat("\nDone. Outputs saved to:\n")
cat("  ", file.path(FIG_DIR, "fig_upset_sharp_rsf.pdf"), "\n")
cat("  ", file.path(FIG_DIR, "fig_upset_sharp_rsf.png"), "\n")
cat("  ", file.path(FIG_DIR, "upset_membership_table.csv"), "\n")

