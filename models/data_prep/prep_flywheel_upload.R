##########################################################
## Preparing cohort data for Flywheel upload
## Brain Health Metrics IHME
## April 2026
##########################################################

## Set up ----------------------------------------------------------------------
rm(list=ls())
username <- Sys.getenv("USER")
library(data.table)

# System configuration
l_root <- "/FILEPATH/"

# Directories
read_dir <- paste0(l_root, "FILEPATH/")
save_dir <- paste0(l_root, "FILEPATH/")

# UPDATE THESE FOR EACH RUN
extract_country <- "Uganda" 
extract_cohort <- "LONISAC" 

## Read Data and Mapping -------------------------------------------------------

# Read the fully processed cohort data
file_name <- paste0(extract_cohort, "_", extract_country, "_harmonized_data.csv")
dt_processed <- fread(paste0(read_dir, file_name))

# Read the new Flywheel mapping spreadsheet
fw_map_path <- paste0(l_root, "FILEPATH/flywheel_variable_mapping.csv")
fw_map <- fread(fw_map_path)


## Subset and Rename -----------------------------------------------------------

# 1. Check which mapped variables are missing (for logging purposes)
expected_cols <- fw_map$standard_name
actual_cols <- names(dt_processed)
missing_vars <- setdiff(expected_cols, actual_cols)

if (length(missing_vars) > 0) {
  message("-------------------------------------------------------")
  message(paste("NOTE:", length(missing_vars), "variables from the Flywheel mapping sheet were NOT found in the dataset."))
  message("These will be created as NA columns in the output to maintain schema.")
  message(paste(missing_vars, collapse = ", "))
  message("-------------------------------------------------------")
}

# 2. Initialize the new Flywheel data.table with the requested beginning columns
dt_fw <- data.table(
  CohortName = rep(extract_cohort, nrow(dt_processed)),
  CohortLocation_country = rep(extract_country, nrow(dt_processed))
)

# 3. Build the dataset in the exact order of the mapping sheet
for (i in 1:nrow(fw_map)) {
  s_name <- fw_map$standard_name[i]
  fw_name <- fw_map$flywheel_standard_name[i]
  
  # If the standard_name exists, copy the data. If not, fill with NA.
  if (s_name %in% names(dt_processed)) {
    dt_fw[, (fw_name) := dt_processed[[s_name]]]
  } else {
    dt_fw[, (fw_name) := NA]
  }
}

# 4. Standardize any blank strings to true NAs
for (col in names(dt_fw)) {
  if (is.character(dt_fw[[col]])) {
    # Replace empty strings or strings that are just spaces with NA
    dt_fw[trimws(get(col)) == "", (col) := NA_character_]
  }
}

## Calculate Coverage Report ---------------------------------------------------

# Define the subject ID column name in the Flywheel schema
subj_col <- "StudyID" 
if (!subj_col %in% names(dt_fw)) {
  subj_col <- "subject_id" # Fallback if standard_name wasn't changed
}

if (subj_col %in% names(dt_fw)) {
  n_total_rows <- nrow(dt_fw)
  n_total_subjects <- uniqueN(dt_fw[[subj_col]])
  
  coverage_list <- list()
  
  for (col in names(dt_fw)) {
    # Calculate Row Coverage
    non_na_rows <- sum(!is.na(dt_fw[[col]]))
    pct_rows <- (non_na_rows / n_total_rows) * 100
    
    # Calculate Subject Coverage (Does this subject have data in ANY of their visits?)
    # dt_fw[!is.na(get(col))] creates a subset of rows with data, then we count unique IDs
    non_na_subjects <- uniqueN(dt_fw[!is.na(get(col))][[subj_col]])
    pct_subjects <- (non_na_subjects / n_total_subjects) * 100
    
    # Save to list
    coverage_list[[col]] <- data.table(
      Column_Name = col,
      Percent_Non_NA_Rows = round(pct_rows, 2),
      Percent_Non_NA_Subjects = round(pct_subjects, 2)
    )
  }
  
  # Combine into final coverage dataset
  dt_coverage <- rbindlist(coverage_list)
  
} else {
  message("Warning: Could not find a subject ID column (StudyID) to calculate subject-level coverage.")
}


## Save Output -----------------------------------------------------------------
# Ensure save directory exists
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Save main Flywheel data
out_file_name <- paste0(extract_cohort, "_", extract_country, "_harmonized_demographics_April2026.csv")
fwrite(dt_fw, paste0(save_dir, out_file_name), na = "NA")

# Save coverage report
coverage_file_name <- paste0(extract_cohort, "_", extract_country, "_flywheel_coverage.csv")
if (exists("dt_coverage")) {
  fwrite(dt_coverage, paste0(save_dir, coverage_file_name))
  message("Success! Coverage report saved to: ", paste0(save_dir, coverage_file_name))
}

message("Success! Flywheel upload file saved to: ", paste0(save_dir, out_file_name))

