# appendix_B_bnf_feature_table.R 

# This script is used to output the table content for the SHARP selected features 
# in the Appendix in the paper. 

# The output will also include the corresponding weights PMCI_PTS allocated for 
# each SHARP selected feature. 
# PMCI_PTS weights are essentially round(PMCI_LP weights * c) for c is defined 
# to be 5 for the moment. 

library(data.table)
library(stringr)


# SHARP results data path 
SHARP_DIR <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583/derived_gp/sharp"

safe_colname <- function(x) {
  x <- str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
  x
}

# BNF code mapping function, identical to used in 01_build_analysis_dataset_gp_complete.R 
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

# Load outputs from SHARP 
stable_f <- fread(file.path(SHARP_DIR, "stable_variables_Female_sharp.csv"))
stable_m <- fread(file.path(SHARP_DIR, "stable_variables_Male_sharp.csv"))

coef_f <- fread(file.path(SHARP_DIR, "cox_refit_coef_Female_sharp.csv"))
coef_m <- fread(file.path(SHARP_DIR, "cox_refit_coef_Male_sharp.csv"))

coef_f <- coef_f[feature != "age_at_entry"]
coef_m <- coef_m[feature != "age_at_entry"]

coef_f[, PMCI_F := round(coef * 5)]
coef_m[, PMCI_M := round(coef * 5)]

# Function to build BNF code mapping 
map_dt <- data.table(
  BNF_Code = names(bnf_code_to_name),
  BNF_Category = unname(bnf_code_to_name)
)
map_dt[, Feature := paste0("bnf_", safe_colname(BNF_Category))]

# Selection flags 
map_dt[, Female := Feature %in% stable_f$feature]
map_dt[, Male   := Feature %in% stable_m$feature]

# Merging PMCI weights 
map_dt <- merge(map_dt, coef_f[, .(Feature = feature, PMCI_F)], by = "Feature", all.x = TRUE)
map_dt <- merge(map_dt, coef_m[, .(Feature = feature, PMCI_M)], by = "Feature", all.x = TRUE)

# Retain only selected features 
final_dt <- map_dt[Female | Male]

# Export file formatting 
final_dt[, Female := ifelse(Female, "Yes", "--")]
final_dt[, Male   := ifelse(Male,   "Yes", "--")]
final_dt[, PMCI_F := ifelse(is.na(PMCI_F), "--", PMCI_F)]
final_dt[, PMCI_M := ifelse(is.na(PMCI_M), "--", PMCI_M)]

final_dt <- final_dt[
  , .(
    BNF_Code,
    BNF_Category,
    Female,
    Male,
    PMCI_F,
    PMCI_M
  )
][order(BNF_Code)]

print(final_dt)

# Export 
fwrite(final_dt, "appendix_sharp_pmci_table.csv")

latex_rows <- final_dt[
  , sprintf(
    "%s & %s & %s & %s & %s & %s \\\\",
    BNF_Code,
    BNF_Category,
    Female,
    Male,
    PMCI_F,
    PMCI_M
  )
]

writeLines(latex_rows, "appendix_sharp_pmci_rows.tex")

cat("Files written:\n")
cat(" - appendix_sharp_pmci_table.csv\n")
cat(" - appendix_sharp_pmci_rows.tex\n")


