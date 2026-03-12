# data_quality_sensitivity_check.R 

# Purpose of this file is to conduct a sensitivity analysis to answer any questions 
# regarding the early censoring date (1990-01-01) applied to all UKB data used 
# in this paper. 

# The aim here is to plot the number of GP prescriptions and UKB participants 
# against years. 

# The reason the date (1990-01-01) was chosen as the early censoring date is due 
# to the fact that electronic medical records were only introduced in the UK 
# towards the end of the 1980s. There exists a large number of false and inaccurate
# GP prescription data before 1990. 

# The two plots indicate that there were almost no UKB participants before 1990 
# and the number of GP prescriptions were negligible. 


library(data.table)
library(lubridate)

# Configurations 
ROOT <- "/rds/general/user/yl226/projects/chadeau_ukbb_folder/live/data/project_data/UKB_677583"
GP_FILE <- file.path(ROOT, "gp_scripts.txt")

parse_date <- function(x) as.IDate(lubridate::dmy(x))

gp <- fread(GP_FILE, select = c("eid","issue_date"))
gp[, issue_date := parse_date(issue_date)]
gp <- gp[!is.na(eid) & !is.na(issue_date)]

gp[, year := year(issue_date)]

year_stats <- gp[, .(
  n_records = .N,
  n_participants = uniqueN(eid)
), by = year][order(year)]

print(year_stats[year <= 1995])

pdf("gp_records_by_year.pdf", width = 8, height = 5)
plot(year_stats$year, year_stats$n_records, type = "l",
     xlim = c(1910, 2020),
     xlab = "Year",
     ylab = "Number of GP prescription records")
abline(v = 1990, lty = 2)
text(1990, max(year_stats$n_records, na.rm = TRUE) * 0.05, "1990", pos = 4)
dev.off()

pdf("gp_participants_by_year.pdf", width = 8, height = 5)
plot(year_stats$year, year_stats$n_participants, type = "l",
     xlim = c(1910, 2020),
     xlab = "Year",
     ylab = "Number of participants")
abline(v = 1990, lty = 2)
text(1990, max(year_stats$n_participants, na.rm = TRUE) * 0.05, "1990", pos = 4)
dev.off()

fwrite(year_stats, "gp_year_stats.csv")
cat("Saved: gp_year_stats.csv, gp_records_by_year.png, gp_participants_by_year.png\n")


