# union_heatmap.R 

# This script outputs the correlation heatmap included in the Appendix of the paper. 
# This heatmap was created as an attempt to discover any close proxies of the two 
# medication categores Iron Deficiency Anaemias and Cardiac Glycosides in the 
# SHARP selected prescription categories for males. 

# The set of features being evaluated is essentially the union of naive (1se) 
# and sharp selected features for males. 


library(data.table)
library(ggplot2)

# Configurations 
ROOT <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp"
DATA_RDS <- file.path(ROOT, "data_full_gp_cov.rds")
NAIVE_DIR <- file.path(ROOT, "stab_naive")
SHARP_DIR <- file.path(ROOT, "sharp")
OUT_DIR   <- file.path(ROOT, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Naive and SHARP Males selected variables 
naive_file <- file.path(NAIVE_DIR, "stable_variables_Male_l1se.csv")
sharp_file <- file.path(SHARP_DIR, "stable_variables_Male_sharp.csv")

# Feature loading functions 
read_feats <- function(path) {
  dt <- fread(path)
  if ("feature" %in% names(dt)) as.character(dt$feature) else as.character(dt[[1]])
}

naive_feats <- read_feats(naive_file)
sharp_feats <- read_feats(sharp_file)

union_feats <- sort(unique(c(naive_feats, sharp_feats)))

# Loading the Male dataset 
df <- as.data.table(readRDS(DATA_RDS))
dfm <- df[sex == "Male"]

union_feats <- intersect(union_feats, names(dfm))

# Building binary numeric matrix 
X <- dfm[, lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = union_feats]
X <- as.matrix(X)
X[is.na(X)] <- 0

# Correlation calculation 
C <- cor(X, use = "pairwise.complete.obs")

# Hierarchical clustering to get a better structure 
hc <- hclust(as.dist(1 - abs(C)))
ord <- hc$order
C_ord <- C[ord, ord]

# Use data.table for ggplot 
dt_hm <- as.data.table(as.table(C_ord))


setnames(
  dt_hm,
  old = names(dt_hm),
  new = c("Feature1", "Feature2", "Corr")
)

dt_hm[, Feature1 := factor(Feature1, levels = rownames(C_ord))]
dt_hm[, Feature2 := factor(Feature2, levels = rownames(C_ord))]

# Plotting 
p <- ggplot(dt_hm, aes(Feature1, Feature2, fill = Corr)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, limits = c(-1, 1), 
    name = "Correlation"
  ) + 
  coord_fixed() +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank()
  ) +
  labs(
    title = "Male: Correlation heatmap of union(SHARP, naive) features"
  )

out_pdf <- file.path(OUT_DIR, "fig_corr_union_male_sharp_vs_naive_l1se.pdf")
ggsave(out_pdf, p, width = 9, height = 8, useDingbats = FALSE)

cat("Saved:", out_pdf, "\n")


