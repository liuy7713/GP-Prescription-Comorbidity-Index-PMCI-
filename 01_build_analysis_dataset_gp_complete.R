# 01_build_analysis_dataset_gp_complete.R 
# Purpose is to build the complete modelling ready dataset from GP prescriptions + death registry 
# Standardise BNF codes into 6 digit paragraph 
# Create dot code "1.1.1" and section "1.1" with no leading zeros 
# Restrict to chapters 1 to 12 and the whitelist (keep sections) 
# Keep common categories (>= 500 people) 
# Manual exclusions and regrouping to reach 124 categories 
# Map dot codes to human readable names 
# Build binary matrix with ever prescribed (1 / 0) for each category 
# Survival t0 is the first prescription date in retained GP data 

# Output: 
# derived_gp/data_full_gp.rds 
# derived_gp/bnf_category_prevalence.csv 

# I am not using the age_at_entry = 2006 - YOB proxy 
# Age and Sex are merged later from Baseline.rds 

library(data.table) 
library(stringr) 
library(lubridate) 
library(Matrix) 

# Configurations 
ROOT <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583" 
GP_FILE <- file.path(ROOT, "gp_scripts.txt")
DEATH_FILE <- file.path(ROOT, "death.txt")
OUT_DIR <- file.path(ROOT, "derived_gp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE) 
OUT_RDS <- file.path(OUT_DIR, "data_full_gp.rds") 
OUT_META <- file.path(OUT_DIR, "bnf_category_prevalence.csv") 

# Remove UKB artefact dates as data before 1990 is not reliable 
MIN_GPDATA_DATE <- as.IDate("1990-01-01") 

# Prevalence filter 
MIN_PATIENTS_PER_BNF <- 500 

FIXED_CENSOR_DATE <- as.IDate(NA) 

parse_date <- function(x) as.IDate(lubridate::dmy(x)) 

# Convert BNF codes from paragraph to dot 
para6_to_dot <- function(p6) { 
  chapter <- as.integer(substr(p6, 1, 2)) 
  section <- as.integer(substr(p6, 3, 4)) 
  paragraph <- as.integer(substr(p6, 5, 6)) 
  paste(chapter, section, paragraph, sep = ".") 
}

dot_to_section <- function(dot) sub("^((\\d+\\.\\d+)).*$", "\\1", dot) 

safe_colname <- function(x) { 
  x <- str_replace_all(x, "[^A-Za-z0-9]+", "_") 
  x <- str_replace_all(x, "_+", "_") 
  x <- str_replace_all(x, "^_|_$", "") 
  x 
}

# Dedicated whitelist 
keep_sections <- c( 
  paste0("1.", 1:9), 
  paste0("2.", 1:13), 
  paste0("3.", 1:11), 
  paste0("4.", 1:4), 
  paste0("4.", 7:11), 
  paste0("5.", 1:3), 
  paste0("6.", 1:6), 
  paste0("7.", 1:4), 
  paste0("8.", 1:2), 
  paste0("9.", 1:2), 
  paste0("10.", 1:2), 
  paste0("11.", 6:6), 
  paste0("12.", 1:3)
)

# Dot code to name mapping 
bnf_code_to_name <- c(
  "1.1.1" = "Antacids and simeticone", 
  "1.1.2" = "Compound alginates and proprietary indigestion preparations", 
  "1.2.0" = "Antispasmodics and other drugs altering gut motility", 
  "1.3.1" = "H2-receptor antagonists", 
  "1.3.3" = "Chelates and complexes", 
  "1.3.5" = "Proton pump inhibitors", 
  "1.3.7" = "Other Antisecretory drugs and mucosal protectants", 
  "1.4.2" = "Antimotility drugs", 
  "1.5.1" = "Aminosalicylates", 
  "1.5.2" = "Corticosteroids(for IBD)", 
  "1.5.3" = "Drugs affecting immune response", 
  "1.6.1" = "Bulk-forming laxatives", 
  "1.6.2" = "Stimulant laxatives", 
  "1.6.3" = "Faecal softeners", 
  "1.6.4" = "Osmotic laxatives", 
  "2.1.1" = "Cardiac glycosides", 
  "2.2.1" = "Thiazides and related diuretics", 
  "2.2.2" = "Loop diuretics", 
  "2.2.3" = "Potassium-sparing diuretics and aldosterone antagonists", 
  "2.2.4" = "Potassium sparing diuretics and compounds", 
  "2.3.2" = "Drugs for arrhythmias", 
  "2.4.0" = "Beta-adrenoceptor blocking drugs", 
  "2.5.1" = "Vasodilator antihypertensive drugs", 
  "2.5.2" = "Centrally-acting antihypertensive drugs", 
  "2.5.4" = "Alpha-adrenoceptor blocking drugs", 
  "2.5.5" = "Renin-angiotensin system drugs", 
  "2.6.1" = "Nitrates", 
  "2.6.2" = "Calcium-channel blockers", 
  "2.6.3" = "Other antianginal drugs", 
  "2.6.4" = "Peripheral vasodilators and related drugs", 
  "2.8.1" = "Parenteral anticoagulants", 
  "2.8.2" = "Oral anticoagulants", 
  "2.9.0" = "Antiplatelet drugs", 
  "2.11.0" = "Antifibrinolytic drugs and haemostatics", 
  "2.12.0" = "Lipid-regulating drugs", 
  "3.1.1" = "Adrenoceptor agonists", 
  "3.1.2" = "Antimuscarinic bronchodilators",
  "3.1.3" = "Theophylline",
  "3.1.5" = "Other Bronchodilators",
  "3.2.0" = "Corticosteroids (respiratory)",
  "3.3.2" = "Leukotriene receptor antagonists",
  "3.4.1" = "Antihistamines",
  "3.4.3" = "Allergic emergencies",
  "3.7.0" = "Mucolytics",
  "3.8.0" = "Aromatic inhalations",
  "3.9.1" = "Cough suppressants",
  "3.9.2" = "Expectorant and demulcent cough preparations",
  "3.10.0" = "Systemic nasal decongestants",
  "4.1.1" = "Hypnotics",
  "4.1.2" = "Anxiolytics",
  "4.2.1" = "Antipsychotic drugs",
  "4.2.3" = "Drugs used for mania and hypomania",
  "4.3.1" = "Tricyclic and related antidepressant drugs",
  "4.3.3" = "Selective serotonin re-uptake inhibitors",
  "4.3.4" = "Other antidepressant drugs",
  "4.7.1" = "Non-opioid analgesics and compound preparations",
  "4.7.2" = "Opioid analgesics",
  "4.7.3" = "Neuropathic pain",
  "4.7.4" = "Antimigraine drugs",
  "4.8.1" = "Control of epilepsy",
  "4.9.1" = "Dopaminergic drugs used in parkinsonism",
  "4.9.2" = "Antimuscarinic drugs used in parkinsonism",
  "4.9.3" = "Essential tremor, chorea, tics and related disorders",
  "4.10.2" = "Nicotine dependence",
  "4.11.0" = "Drugs for dementia",
  "5.1.1" = "Penicillins",
  "5.1.2" = "Cephalosporins and other beta-lactams",
  "5.1.3" = "Tetracyclines",
  "5.1.5" = "Macrolides",
  "5.1.6" = "Clindamycin and lincomycin",
  "5.1.8" = "Sulfonamides and trimethoprim",
  "5.1.11" = "Metronidazole, tinidazole and ornidazole",
  "5.1.12" = "Quinolones",
  "5.1.13" = "Urinary-tract infections",
  "5.2.1" = "Triazole antifungals",
  "5.2.3" = "Polyene antifungals",
  "5.2.5" = "Other antifungals",
  "5.3.0" = "Other Antiviral drugs",
  "5.3.2" = "Herpesvirus infections",
  "6.1.1" = "Insulin",
  "6.1.2" = "Antidiabetic drugs",
  "6.1.4" = "Treatment of hypoglycaemia",
  "6.1.6" = "Diabetic diagnostic and monitoring agents",
  "6.2.1" = "Thyroid hormones",
  "6.2.2" = "Antithyroid drugs",
  "6.3.2" = "Glucocorticoid therapy",
  "6.4.1" = "Female sex hormones and their modulators",
  "6.4.2" = "Male sex hormones and antagonists",
  "6.5.1" = "Hypothalamic & anterior pituitary hormone & antioestrogens",
  "6.6.2" = "Bisphosphonates and other drugs",
  "7.2.1" = "Preparations for vaginal and vulval changes",
  "7.2.2" = "Vaginal and vulval infections",
  "7.3.1" = "Combined hormonal contraceptives and systems",
  "7.3.2" = "Progestogen-only contraceptives",
  "7.3.3" = "Spermicidal contraceptives",
  "7.3.5" = "Emergency contraception",
  "7.4.1" = "Drugs for urinary retention",
  "7.4.2" = "Drugs for urinary frequency enuresis and incontinence",
  "7.4.3" = "Drugs used in urological pain",
  "7.4.5" = "Drugs for erectile dysfunction",
  "8.1.3" = "Antimetabolites",
  "8.2.2" = "Corticosteroids and other immunosuppressants",
  "9.1.1" = "Iron-deficiency anaemias",
  "9.1.2" = "Drugs used in megaloblastic anaemias",
  "9.1.3" = "Hypoplastic, haemolytic and renal anaemias",
  "9.2.1" = "Oral preparation for fluid and electrolyte imbalance",
  "9.2.2" = "Parent prepn for fluid and electrolyte imb",
  "10.1.1" = "Non-steroidal anti-inflammatory drugs",
  "10.1.2" = "Corticosteroids(for Musculoskeletal use)",
  "10.1.3" = "Rheumatic disease suppressant drugs",
  "10.1.4" = "Gout and cytotoxic induced hyperuicaemia",
  "10.1.5" = "Other drugs for rheumatic diseases",
  "10.2.2" = "Skeletal muscle relaxants",
  "11.6.0" = "Treatment of glaucoma",
  "12.1.1" = "Otitis externa",
  "12.1.3" = "Removal of ear wax and other substances",
  "12.2.1" = "Drugs used in nasal allergy",
  "12.2.2" = "Topical nasal decongestants",
  "12.2.3" = "Nasal preparations for infection",
  "12.3.1" = "Drugs for oral ulceration and inflammation",
  "12.3.2" = "Oropharyngeal anti-infective drugs",
  "12.3.3" = "Lozenges and sprays",
  "12.3.4" = "Mouth-washes, gargles and dentifrices",
  "12.3.5" = "Treatment of dry mouth"
)

# Loading Death Registry 
cat("Loading death file: ", DEATH_FILE, "\n")
death <- fread(DEATH_FILE) 

if (!all(c("eid", "date_of_death") %in% names(death))) { 
  stop("death.txt does not have expected columns (eid, date_of_death). Found: ", paste(names(death), collapse = ", ")) 
}

death[, eid := as.numeric(eid)]
death[, death_date := parse_date(date_of_death)] 
death <- death[!is.na(eid)]

censor_date <- if (!is.na(FIXED_CENSOR_DATE)) FIXED_CENSOR_DATE else max(death$death_date, na.rm = TRUE) 
cat("Censor date: ", as.character(censor_date), "\n\n")

# Load GP Prescriptions 
cat("Loading GP scripts: ", GP_FILE, "\n")
gp <- fread(GP_FILE) 

needed <- c("eid", "issue_date", "bnf_code") 
missing <- setdiff(needed, names(gp)) 
if (length(missing) > 0) { 
  stop("gp_scripts.txt missing columns: ", paste(missing, collapse = ", "), "\nFound: ", paste(names(gp), collapse = ", ")) 
}

gp[, eid := as.numeric(eid)]
gp[, issue_date := parse_date(issue_date)] 
gp <- gp[!is.na(eid) & !is.na(issue_date)] 
MAX_GPDATA_DATE <- as.IDate("2024-12-31") 
gp <- gp[issue_date <= MAX_GPDATA_DATE]
# Remove very early dates 
gp <- gp[issue_date >= MIN_GPDATA_DATE]

# Standardise bnf_code 
gp[, bnf_code := as.character(bnf_code)] 
gp[, bnf_digits := str_replace_all(bnf_code, "[^0-9]", "")] 
gp[, bnf_para6 := substr(bnf_digits, 1, 6)] 
gp <- gp[!is.na(bnf_para6) & nchar(bnf_para6) == 6]

# Dot code and section 
gp[, bnf_dot_code := para6_to_dot(bnf_para6)] 
gp[, bnf_section := dot_to_section(bnf_dot_code)] 

# Restrict to whitelist 
before <- nrow(gp) 
gp <- gp[bnf_section %in% keep_sections] 
cat("After whitelist filter: \n")
cat("Rows before: ", before, "\n")
cat("Rows after: ", nrow(gp), "\n")
cat("Unique eids: ", uniqueN(gp$eid), "\n") 
cat("Unique dot: ", uniqueN(gp$bnf_dot_code), "\n\n")

# Keep one record per (eid, bnf_dot_code) with earliest issue_date 
setorder(gp, eid, bnf_dot_code, issue_date) 
gp_uniq <- gp[, .SD[1L], by = .(eid, bnf_dot_code)] 
cat("Rows after de duplication: ", nrow(gp_uniq), "\n\n")

# Merge in the prevalence filter 
bnf_counts <- gp_uniq[, .(n_persons = uniqueN(eid)), by = bnf_dot_code] 
setorder(bnf_counts, -n_persons) 
common_bnf <- bnf_counts[n_persons >= MIN_PATIENTS_PER_BNF, bnf_dot_code] 
cat("Common dot codes (>= ", MIN_PATIENTS_PER_BNF, " persons): ", length(common_bnf), "\n\n", sep = "") 
gp_common <- gp_uniq[bnf_dot_code %in% common_bnf] 

# Manual exclusions to reach 124 categories 
# Remove unclear and unsure categories 
gp_common <- gp_common[!(bnf_dot_code %in% c("1.1.3", "2.2.7", "5.2.0", "6.1.5", "7.3.4"))]
# Remove categories 1.7 and 1.9 
gp_common <- gp_common[!str_starts(bnf_dot_code, "1.7") & !str_starts(bnf_dot_code, "1.9")]
# Group unspecific drugs into broader category 
gp_common[, bnf_dot_code := fifelse(str_detect(bnf_dot_code, "^1\\.2\\.[0-9]+"), "1.2.0", bnf_dot_code)]
gp_common[, bnf_dot_code := fifelse(str_detect(bnf_dot_code, "^2\\.12\\.[0-9]+"), "2.12.0", bnf_dot_code)]
gp_common[, bnf_dot_code := fifelse(str_detect(bnf_dot_code, "^2\\.11\\.[0-9]+"), "2.11.0", bnf_dot_code)]
gp_common[, bnf_dot_code := fifelse(str_detect(bnf_dot_code, "^2\\.9\\.[0-9]+"),  "2.9.0",  bnf_dot_code)]
gp_common[, bnf_dot_code := fifelse(str_detect(bnf_dot_code, "^11\\.6\\.[0-9]+"), "11.6.0", bnf_dot_code)]

cat("Final unique dot codes after: ", uniqueN(gp_common$bnf_dot_code), "\n\n")

setorder(gp_common, eid, bnf_dot_code, issue_date) 
gp_common <- gp_common[, .SD[1L], by = .(eid, bnf_dot_code)] 

# Map dot codes to category names 
gp_common[, bnf_name := bnf_code_to_name[as.character(bnf_dot_code)]]
gp_common[is.na(bnf_name) | bnf_name == "", bnf_name := bnf_dot_code] 
gp_common[, bnf_col := paste0("bnf_", safe_colname(bnf_name))] 

# Define survival 
first_rx <- gp_common[, .(first_rx_date = min(issue_date, na.rm = TRUE)), by = eid] 
cat("First Rx date range: ", as.character(min(first_rx$first_rx_date)), "to", as.character(max(first_rx$first_rx_date)), "\n\n") 
surv <- merge(first_rx, death[, .(eid, death_date)], by = "eid", all.x = TRUE) 
surv[, event := fifelse(!is.na(death_date) & death_date <= censor_date, 1L, 0L)] 
surv[, end_date := fifelse(event == 1L, death_date, censor_date)] 
surv[, time := as.numeric(as.IDate(end_date) - as.IDate(first_rx_date))] 

bad <- surv[is.na(time) | time <= 0, .N] 
cat("Dropping rows with time<=0 or NA: ", bad, "\n\n")
surv <- surv[!is.na(time) & time > 0] 

cat("Survival summary: \n")
cat("n: ", nrow(surv), "\n")
cat("events: ", sum(surv$event), "\n\n")

gp_common <- gp_common[eid %in% surv$eid]

# Category prevalence filter 
cat_prev <- gp_common[, .(n_patients = uniqueN(eid)), by = bnf_col] 
setorder(cat_prev, -n_patients) 
keep_cols <- cat_prev[n_patients >= MIN_PATIENTS_PER_BNF, bnf_col] 
cat("Keeping categories with >= ", MIN_PATIENTS_PER_BNF, " patients: ", length(keep_cols), "\n\n", sep = "") 
gp_keep <- gp_common[bnf_col %in% keep_cols] 

# Build sparse binary matrix 
eids <- sort(unique(surv$eid)) 
cols <- sort(unique(keep_cols)) 
eid_index <- setNames(seq_along(eids), as.character(eids)) 
col_index <- setNames(seq_along(cols), cols) 
ij <- gp_keep[, .(i = eid_index[as.character(eid)], j = col_index[bnf_col])]
ij <- unique(ij) 

X <- sparseMatrix(i = ij$i, j = ij$j, x = 1L, dims = c(length(eids), length(cols)), dimnames = list(eid = eids, feature = cols)) 

cat("Sparse design matrix: \n")
cat("n_eids: ", nrow(X), "\n")
cat("n_features: ", ncol(X), "\n")
cat("nnz: ", length(X@x), "\n\n")

# Convert to dense data.table format 
X_dt <- as.data.table(as.matrix(X)) 
X_dt[, eid := as.numeric(rownames(as.matrix(X)))] 

# Assemble final dataset 
data_full <- merge(surv[, .(eid, first_rx_date, time, event)], X_dt, by = "eid", all.x = TRUE) 
feat_cols <- setdiff(names(data_full), c("eid", "first_rx_date", "time", "event")) 
for (cc in feat_cols) { 
  set(data_full, which(is.na(data_full[[cc]])), cc, 0L) 
}

cat("Final dataset: \n")
cat("Rows: ", nrow(data_full), "\n")
cat("Cols: ", ncol(data_full), "\n")
cat("Features: ", length(feat_cols), "\n\n")

saveRDS(data_full, OUT_RDS) 
fwrite(cat_prev, OUT_META) 

cat("Saved: \n")
cat(OUT_RDS, "\n", sep = "") 
cat(OUT_META, "\n", sep = "") 
cat("Done. \n")
