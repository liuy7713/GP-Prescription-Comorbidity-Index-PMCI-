# 07_rsf_tune_one.R

# Reads tuning metrics from Code/rsf_grid.csv
# Run ONE RSF tuning fit for ONE (sex, nodesize, mtry, nsplit) 
# Writes one-row CSV into OUT_DIR/RUN_TAG/tuning_chunks/.

# This file is the first file for RSF Hyper-Parameter Tuning. 
# It is run for every row of Code/rsf_grid.csv 
# This script does not conduct any metrics comparison and only outputs tuning 
# fit results for the given tuning parameters. 

# These initial tuning runs do not incorporate variable importance during fitting 
# Variable importance is currently set to 'none' as it is not necessary during 
# tuning. Variable importance will become relevant only in the final refitting 
# conducted in the next rsf file 07_rsf_aggregate_refit_CI.R 

# This file can only be run in the HPC, does not run on RStudio. 


library(data.table) 
library(randomForestSRC) 
library(survival) 


# Configurations 
DATA_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"
EN_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/en_cox"
OUT_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/rsf"

SEED <- 20260119
BLOCK_SIZE <- 100

# Threads from PBS (Loads the allocated number of threads from HPC) 
NTHREAD <- as.integer(Sys.getenv("PBS_NCPUS", unset = "8"))
data.table::setDTthreads(NTHREAD)
Sys.setenv(OMP_NUM_THREADS = NTHREAD,
           MKL_NUM_THREADS = 1,
           OPENBLAS_NUM_THREADS = 1)

# Run tag to avoid overwriting (set by PBS script)
RUN_TAG <- Sys.getenv("RUN_TAG", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
OUT_DIR_RUN <- file.path(OUT_DIR, RUN_TAG)
CHUNK_DIR <- file.path(OUT_DIR_RUN, "tuning_chunks")
dir.create(CHUNK_DIR, recursive = TRUE, showWarnings = FALSE)

# Loading Train Test Split (Same as previous models) 
load_split_idx <- function(sex_name, n) {
  sfile <- file.path(EN_DIR, paste0("split_idx_", sex_name, ".rds"))
  if (!file.exists(sfile)) stop("Missing split index file: ", sfile)
  sp <- readRDS(sfile)
  tr <- sp$tr; te <- sp$te
  if (any(tr < 1 | tr > n) || any(te < 1 | te > n)) stop("split indices out of range for ", sex_name)
  list(tr = tr, te = te, seed = sp$seed)
}

# Function to obtain feature columns (Same as previous models) 
get_feature_cols <- function(df_sex, sex_name) {
  bnf_cols <- grep("^bnf_", names(df_sex), value = TRUE)
  
  tmp <- df_sex[, ..bnf_cols]
  tmp <- as.data.frame(lapply(tmp, as.numeric))
  nonzero <- colSums(tmp, na.rm = TRUE) > 0
  bnf_cols <- bnf_cols[nonzero]
  if (length(bnf_cols) == 0) stop("No nonzero bnf_ columns found for ", sex_name)
  
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

# Out of Bag (OOB) C-Index Function for tuning only 
# OOB evaluates the function on a subset of the training data using the portion 
# of the ensemble which was not trained on this subset to avoid information 
# leakage and unreliable internal reference C-Index 
oob_cindex_rsf <- function(fit) {
  1 - tail(fit$err.rate, 1)
}

# Getting the arguments such as rsf_grid.csv and task_id within it. 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 07_rsf_tune_one.R <grid_csv> <task_id>")
}
grid_csv <- args[1]
task_id_raw <- args[2]
task_id <- suppressWarnings(as.integer(trimws(task_id_raw)))

grid <- fread(grid_csv)
stopifnot(all(c("task_id","sex","nodesize","mtry","nsplit","ntree") %in% names(grid)))

# Force task_id column value into integers 
grid[, task_id := suppressWarnings(as.integer(trimws(as.character(task_id))))]

if (is.na(task_id)) {
  cat("DEBUG: task_id_raw='", task_id_raw, "' parsed NA\n", sep="")
  quit(status = 2)
}

# Debug code (Not used if task_id was allocated and parsed correctly) 
cat("DEBUG: grid_csv=", grid_csv, "\n", sep="")
cat("DEBUG: task_id_raw='", task_id_raw, "' task_id=", task_id, "\n", sep="")
cat("DEBUG: grid task_id class=", class(grid$task_id), "\n", sep="")
cat("DEBUG: first task_ids=", paste(head(grid$task_id, 20), collapse=","), "\n", sep="")

hits <- which(grid$task_id == task_id)
cat("DEBUG: hits=", length(hits), "\n", sep="")

row <- grid[hits]
if (nrow(row) != 1) {
  cat("DEBUG: rows matching task_id:\n")
  print(grid[task_id == task_id])
  stop("Could not find unique row for task_id=", task_id)
}

sex_name <- as.character(row$sex[1])
nodesize <- as.integer(row$nodesize[1])
mtry <- as.integer(row$mtry[1])
nsplit <- as.integer(row$nsplit[1])
ntree <- as.integer(row$ntree[1])

cat("RUN_TAG=", RUN_TAG, "\n", sep="")
cat("Task ", task_id, ": sex=", sex_name,
    " nodesize=", nodesize, " mtry=", mtry, " nsplit=", nsplit, " ntree=", ntree, "\n", sep="")

# Loading Data 
df <- as.data.table(readRDS(DATA_RDS))
stopifnot(all(c("sex", "time", "event", "age_at_entry") %in% names(df)))
df[, event := as.integer(event)]
df[, time := as.numeric(time)]
df[, age_at_entry := as.numeric(age_at_entry)]

df_sex <- df[sex == sex_name]
n <- nrow(df_sex)
if (n == 0) stop("No rows for sex=", sex_name)

# Features
bnf_cols <- get_feature_cols(df_sex, sex_name)
rhs <- paste(c("age_at_entry", bnf_cols), collapse = " + ")
fml <- as.formula(paste0("Surv(time, event) ~ ", rhs))

# Split
sp <- load_split_idx(sex_name, n)
dtr <- df_sex[sp$tr]


# Fitting only one tune function 
set.seed(SEED)
t0 <- Sys.time()

# rfsrc function used solely to run one specific tune 
# Variable importance set to 'none' 
fit <- rfsrc(
  formula = fml,
  data = dtr,
  ntree = ntree,
  mtry = mtry,
  nodesize = nodesize,
  nsplit = nsplit,
  importance = "none",
  nthread = NTHREAD,
  block.size = BLOCK_SIZE
)

t1 <- Sys.time()
oob_c <- oob_cindex_rsf(fit)
elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))

cat(sprintf("DONE: OOB C=%.4f | %.1fs\n", oob_c, elapsed_sec))

# Write the results into the Chunk csv files 
out <- data.table(
  run_tag = RUN_TAG,
  task_id = task_id,
  sex = sex_name,
  ntree = ntree,
  mtry = mtry,
  nodesize = nodesize,
  nsplit = nsplit,
  oob_cindex = oob_c,
  elapsed_sec = elapsed_sec,
  nthread = NTHREAD,
  block_size = BLOCK_SIZE
)

chunk_file <- file.path(CHUNK_DIR, sprintf("tune_task%03d_%s.csv", task_id, sex_name))
fwrite(out, chunk_file)

cat("Wrote: ", chunk_file, "\n", sep="")


