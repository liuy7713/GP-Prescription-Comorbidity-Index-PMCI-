# 02_merge_baseline_sex_age.R 
# Purpose of this file is to merge the dataset with baseline sex and age variables
# The baseline sex and age variables come from the file Baseline.rds
# Age is computed as true age at first GP prescription date
# using Year-of-birth and Month-of-birth from Baseline.rds
# Day of birth is not available so the 15th of the birth month is used
# as a standard midpoint approximation

library(data.table)

BASELINE_RDS <- "/rds/general/user/rw1317/projects/chadeau_ukbb_folder/live/projects/Comorbidity-index/Data/Baseline.rds"
GP_RDS       <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp.rds"
OUT_RDS      <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"

cat("Loading Baseline ...\n")
baseline <- readRDS(BASELINE_RDS)
baseline <- as.data.table(baseline, keep.rownames = "eid")
baseline[, eid := as.numeric(eid)]

SEX_COL         <- "Sex"
YOB_COL         <- "Year-of-birth"
MOB_COL         <- "Month-of-birth"

# Extract sex, year of birth and month of birth
base_small <- baseline[, .(
  eid,
  sex_raw        = get(SEX_COL),
  year_of_birth  = as.integer(get(YOB_COL)),
  month_of_birth = as.integer(get(MOB_COL))
)]

base_small[, sex := fifelse(
  sex_raw %in% c("Female", "F", 0, "0", 2, "2"), "Female",
  fifelse(sex_raw %in% c("Male", "M", 1, "1"), "Male", NA_character_)
)]

base_small <- base_small[!is.na(sex) & !is.na(year_of_birth) & !is.na(month_of_birth)]

cat("Baseline rows after filtering: ", nrow(base_small), "\n")
cat("Sex distribution in baseline:\n")
print(base_small[, .N, by = sex])

cat("Loading GP modelling data ...\n")
df <- readRDS(GP_RDS)
df <- as.data.table(df)
df[, .rowid := .I]

cat("Merging ...\n")
df2 <- merge(df, base_small[, .(eid, sex, year_of_birth, month_of_birth)], by = "eid", all.x = TRUE)
setorder(df2, .rowid)
df2[, .rowid := NULL]

cat("Rows before merge: ", nrow(df), "\n")
cat("Rows after merge:  ", nrow(df2), "\n")
cat("Missing sex:           ", sum(is.na(df2$sex)), "\n")
cat("Missing year of birth: ", sum(is.na(df2$year_of_birth)), "\n")

# Compute age at first prescription
# Day of birth is unknown so the 15th of the birth month is used as a midpoint
df2[, dob_approx := as.Date(paste(year_of_birth, sprintf("%02d", month_of_birth), "15", sep = "-"))]
df2[, age_at_entry := as.numeric(as.Date(first_rx_date) - dob_approx) / 365.25]

# Clean up intermediate columns
df2[, c("year_of_birth", "month_of_birth", "dob_approx") := NULL]

cat("Missing age_at_entry: ", sum(is.na(df2$age_at_entry)), "\n")

# Sanity check on age range
cat("Age at entry summary:\n")
print(summary(df2$age_at_entry))

# Drop rows with missing sex or age
df2 <- df2[!is.na(sex) & !is.na(age_at_entry)]
cat("Rows after dropping missing sex or age: ", nrow(df2), "\n")

cat("Sex distribution:\n")
print(df2[, .N, by = sex])

saveRDS(df2, OUT_RDS)
cat("Saved: ", OUT_RDS, "\n")