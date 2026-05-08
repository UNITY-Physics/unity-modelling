######################################################################
## Create dataset of latest outcome available per subject for FlyWheel
## Brain Health Metrics 
## IHME
## August 2025
#####################################################################


## Set up ---------------------------------------------------------------------------------------------------
rm(list=ls())
username <- Sys.getenv("USER")

library(data.table)
library(dplyr)

#System configuration
l_root <- "/FILEPATH/"

#Set directories
read_dir_flywheel <- paste0(l_root, "FILEPATH/")

save_dir <- paste0(l_root, "FILEPATH/")


## Data manipulation ---------------------------------------------------------------------------------------------------
orig_dt <- as.data.table(read.csv(paste0(read_dir_flywheel, "UNITY-dataTemplate_MRI_clean.csv")))

orig_dt[is.na(childTimepointAge_days) & !is.na(childTimepointAge_months), childTimepointAge_days := round(childTimepointAge_months*30.44, digits = 0)]

# Step 0: Select relevant columns from original data
cols_to_keep <- c(
  "CohortName", "CohortLocation_city", "CohortLocation_country", 
  "StudyID", "UniqueStudyID", "childTimepointAge_days", 
  "gsed_LongForm_DAZScore", "bsid_CCS", "bsid_LCS", "bsid_MCS"
)
dt <- subset(orig_dt, select = cols_to_keep)

# Step 1: Define a function to process each domain
process_domain <- function(dt, score_name, score_col) {
  id_cols <- c("CohortName", "CohortLocation_city", "CohortLocation_country", 
               "StudyID", "UniqueStudyID")
  keep_cols <- c(id_cols, "childTimepointAge_days", score_col)
  
  # Subset to relevant columns and remove NAs
  domain_dt <- dt[, ..keep_cols]
  domain_dt <- domain_dt[!is.na(domain_dt[[score_col]])]
  
  # Print number of UniqueStudyIDs before filtering
  n_before <- length(unique(domain_dt$UniqueStudyID))
  message(score_col, ": UniqueStudyID count before filtering = ", n_before)
  
  # Keep only the row with the maximum age per UniqueStudyID
  domain_dt <- domain_dt[domain_dt[, .I[which.max(childTimepointAge_days)], by = UniqueStudyID]$V1]
  
  # Print number of UniqueStudyIDs after filtering
  n_after <- length(unique(domain_dt$UniqueStudyID))
  message(score_col, ": UniqueStudyID count after filtering = ", n_after)
  
  # Rename the columns
  setnames(domain_dt, 
           old = c("childTimepointAge_days", score_col),
           new = c(paste0("terminal_", score_name, "_Age_days"), paste0("terminal_", score_name)))
  
  return(domain_dt)
}

# Step 2: Process each domain
dt_gsed <- process_domain(dt, "gsed_LongForm_DAZScore", "gsed_LongForm_DAZScore")
dt_ccs  <- process_domain(dt, "bsid_CCS", "bsid_CCS")
dt_lcs  <- process_domain(dt, "bsid_LCS", "bsid_LCS")
dt_mcs  <- process_domain(dt, "bsid_MCS", "bsid_MCS")

# Step 3: Merge all processed datasets by key identifiers
id_cols <- c("CohortName", "CohortLocation_city", "CohortLocation_country", 
             "StudyID", "UniqueStudyID")

merged_dt <- Reduce(function(x, y) merge(x, y, by = id_cols, all = TRUE), 
                    list(dt_gsed, dt_ccs, dt_lcs, dt_mcs))

write.csv(merged_dt, paste0(save_dir, "flywheel_terminal_gsed_bsid.csv"), row.names = FALSE)
