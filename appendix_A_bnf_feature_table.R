# appendix_A_bnf_feature_table.R 

# This script is designed to output the table content for the complete 124 BNF 
# features table in the Appendix in the paper 


library(data.table)
library(stringr)

# Function matching preprocessing in 01_build_analysis_dataset_gp_complete.R 
safe_colname <- function(x) {
  x <- str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
  x
}

# Identical function to 01_build_analysis_dataset_gp_complete.R 
# Function used to correspond BNF codes to drug category name 
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

# Load modelling dataset from derived_gp 
DATA_RDS <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/data_full_gp_cov.rds"
DF <- as.data.table(readRDS(DATA_RDS))

feature_cols <- grep("^bnf_", names(DF), value = TRUE)

# Function to build dot code to name 
map_dt <- data.table(
  BNF_Code = names(bnf_code_to_name),
  BNF_Category = unname(bnf_code_to_name)
)
map_dt[, Feature_Name := paste0("bnf_", safe_colname(BNF_Category))]

# Keep only features that exist in DF
map_dt <- map_dt[Feature_Name %in% feature_cols]

# Include the prevalence and sex exclusion criteria 
X <- as.matrix(DF[, ..feature_cols])

n_patients <- colSums(X == 1)

is_female <- DF$sex == "Female"
is_male   <- DF$sex == "Male"

female_any <- colSums(X[is_female, , drop = FALSE] == 1) > 0
male_any   <- colSums(X[is_male,   , drop = FALSE] == 1) > 0

sex_dt <- data.table(
  Feature_Name = feature_cols,
  Sex = ifelse(female_any & male_any, "Female, Male",
               ifelse(female_any, "Female only",
                      ifelse(male_any, "Male only", "None"))),
  N_Patients = as.integer(n_patients)
)

# Merge and export the file 
appendix_table <- merge(map_dt, sex_dt, by = "Feature_Name", all.x = TRUE)
setorder(appendix_table, BNF_Code)

OUT_CSV <- "appendix_bnf_features_with_codes.csv"
fwrite(appendix_table, OUT_CSV)

cat("Saved:", OUT_CSV, "\n")
print(head(appendix_table, 10))

