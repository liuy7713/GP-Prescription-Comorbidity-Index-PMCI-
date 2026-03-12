# Purpose of this file is to merge the dataset with baseline sex and age variables 
# The baseline sex and age variables come from the file Baseline.rds 

library(data.table) 

BASELINE_RDS <- "/rds/general/user/rw1317/projects/chadeau_ukbb_folder/live/projects/Comorbidity-index/Data/Baseline.rds"
GP_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp.rds"
OUT_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"

cat("Loading Baseline ... \n")
baseline <- readRDS(BASELINE_RDS) 
baseline <- as.data.table(baseline, keep.rownames = "eid") 
baseline[, eid := as.numeric(eid)] 

SEX_COL <- "Sex" 
AGE_COL <- if ("Age-at-recruitment" %in% names(baseline)) "Age-at-recruitment" else "Age-when-attended-assessment-centre"

base_small <- baseline[, .(eid, sex_raw = get(SEX_COL), age_at_entry = as.numeric(get(AGE_COL)))] 
base_small[, sex := fifelse(sex_raw %in% c("Female", "F", 0, "0", 2, "2"), "Female", fifelse(sex_raw %in% c("Male", "M", 1, "1"), "Male", NA_character_))] 
base_small <- base_small[!is.na(sex) & !is.na(age_at_entry)]

cat("Loading GP modelling data ... \n")
df <- readRDS(GP_RDS) 
df <- as.data.table(df) 

cat("Merging ... \n")
df2 <- merge(df, base_small[, .(eid, sex, age_at_entry)], by = "eid", all.x = TRUE) 

cat("Rows after dropping missing sex or age: ", nrow(df2), "\n") 
cat("Sex distribution: \n")
print(df2[, .N, by = sex]) 

saveRDS(df2, OUT_RDS) 
cat("Saved: ", OUT_RDS, "\n")

