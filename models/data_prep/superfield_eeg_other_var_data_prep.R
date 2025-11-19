######################################################################
## Create dataset of subcortical MRI, demographics, EEG for predictive analysis
## Brain Health Metrics 
## IHME
## August 2025
#####################################################################


## Set up ---------------------------------------------------------------------------------------------------
rm(list=ls())
username <- Sys.getenv("USER")

library(readxl)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)

#System configuration
l_root <- "/FILEPATH/"

#Set directories
code_dir <- paste0("/FILEPATH/", username, "/FILEPATH/")
read_dir_flywheel <- paste0(l_root, "FILEPATH/")
read_dir_centiles <- paste0(l_root, "FILEPATH/")

save_dir <- paste0(l_root, "FILEPATH/")

#Source functions
source(paste0(code_dir, "flywheel_cohorts/flywheel_functions.R"))

#Set dataset type
# dataset_type <- "all_cohorts_mri"
# dataset_type <- "khula_space_eeg"
# dataset_type <- "khula_space_mri_eeg"
dataset_type <- "khula_mri_eeg"


if (dataset_type == "all_cohorts_mri") {
  ##ALL COHORTS - ORIG, NORM -----------------------------------------------------
  #Set objects
  other_vars <- c('childBirthWeight_kgs', 'childBirthLength_inches', 
                  'MaternalEducation_schoolingYears', 'childBiologicalSex')
  outfile_name <- "all_cohorts_sf_oth_gsed_outcome.csv"
  
  # PREP MRI -- 
  mri_dt <- as.data.table(read.csv(paste0(read_dir_flywheel, "UNITY-dataTemplate_MRI_clean.csv")))
  
  #Merging on SPACE superfield data - temporary until gets added onto main sheet
  dt_space <- as.data.table(read.csv(paste0(read_dir_flywheel, "UCT_SPACE_Superfield_LoMINA_GSED_subset_updated.csv")))
  
  dt_space[, StudyID := sub("_(?!.*_).*", "", subject_id, perl = TRUE)]
  
  dt_space_orig <- mri_dt[CohortName == "SPACE_Explore", ]
  dt_space_orig <- subset(dt_space_orig, select = -c(sf_csf, sf_gray_matter, sf_white_matter, sf_total_icv, sf_left_caudate, 
                                                     sf_left_lentiform, sf_left_hippocampus, sf_right_caudate, sf_right_lentiform, sf_right_hippocampus))
  
  setnames(dt_space, c("CSF.Volume..mm3.", "Gray.Matter.Volume..mm3.", "White.Matter.Volume..mm3.", "Total.Intracranial.Volume..mm3.", "LoMINA.SC.Left.Caudate.mm3", 
                       "LoMINA.SC.Left.Lentiform.mm3", "LoMINA.SC.Left.Hippocampus.mm3", "LoMINA.SC.Right.Caudate.mm3", "LoMINA.SC.Right.Lentiform.mm3", "LoMINA.SC.Right.Hippocampus.mm3"), 
           c("sf_csf", "sf_gray_matter", "sf_white_matter", "sf_total_icv", "sf_left_caudate", 
             "sf_left_lentiform", "sf_left_hippocampus", "sf_right_caudate", "sf_right_lentiform", "sf_right_hippocampus"))
  dt_space <- subset(dt_space, select = c(StudyID, studyTimepoint, sf_csf, sf_gray_matter, sf_white_matter, sf_total_icv, sf_left_caudate, 
                                          sf_left_lentiform, sf_left_hippocampus, sf_right_caudate, sf_right_lentiform, sf_right_hippocampus))
  
  dt_space_orig <- merge(dt_space_orig, dt_space, by = c("StudyID", "studyTimepoint"), all.x = TRUE)
  
  mri_dt <- mri_dt[!CohortName == "SPACE_Explore", ]
  mri_dt <- rbind(mri_dt, dt_space_orig)
  
  mri_dt[, grey_white_ratio := sf_gray_matter/sf_white_matter] 
  mri_dt[, sf_caudate := (sf_left_caudate+sf_right_caudate)/2]
  mri_dt[, sf_lentiform := (sf_left_lentiform+sf_right_lentiform)/2]
  mri_dt[, sf_hippocampus := (sf_left_hippocampus+sf_right_hippocampus)/2]
  
  mri_dt <- mri_dt[,c("CohortName", "StudyID", "childTimepointAge_days", "sf_csf", "sf_gray_matter", 
                      "sf_white_matter", "sf_total_icv", "sf_caudate", "sf_lentiform", "sf_hippocampus")]
  mri_dt <- unique(mri_dt)
  
  #Normalize MRI
  cols_to_normalize <- c("sf_csf", "sf_gray_matter", "sf_white_matter", "sf_caudate", "sf_lentiform", "sf_hippocampus")
  normalize_by_col <- c("sf_total_icv")
  
  mri_dt <- mri_dt[!is.na(get(normalize_by_col)), ]
  
  mri_dt <- normalize_mri(mri_dt, cols_to_normalize, normalize_by_col)
  
  mri_dt[, paste0(normalize_by_col, "_norm") := get(normalize_by_col)]
  
  setnames(mri_dt, c("sf_csf", "sf_gray_matter", "sf_white_matter", "sf_total_icv", "sf_caudate", "sf_lentiform", "sf_hippocampus"), 
           c("sf_csf_orig_value", "sf_gray_matter_orig_value", "sf_white_matter_orig_value", "sf_total_icv_orig_value", "sf_caudate_orig_value", "sf_lentiform_orig_value", "sf_hippocampus_orig_value"))
  
  setnames(mri_dt, c("sf_csf_norm", "sf_gray_matter_norm", "sf_white_matter_norm", "sf_total_icv_norm", "sf_caudate_norm", "sf_lentiform_norm", "sf_hippocampus_norm"), 
           c("sf_csf_norm_value", "sf_gray_matter_norm_value", "sf_white_matter_norm_value", "sf_total_icv_norm_value", "sf_caudate_norm_value", "sf_lentiform_norm_value", "sf_hippocampus_norm_value"))
  
  setnames(mri_dt, c("StudyID", "childTimepointAge_days"), c("subject_id", "mri_age_days"))
  
  mri_dt <- mri_dt[!is.na(mri_age_days),]
  mri_dt[, mri_age_days := as.numeric(mri_age_days)]
  
  mri_wide <- as.data.table(mri_dt)
  
  # PREP OUTCOME -- 
  outcome_dt <- as.data.table(read.csv(paste0(read_dir_flywheel, "UNITY-dataTemplate_MRI_clean.csv")))
  
  outcome_dt <- outcome_dt[,c("CohortName", "StudyID", "childTimepointAge_days", "gsed_LongForm_DAZScore")]
  outcome_dt <- outcome_dt[!is.na(gsed_LongForm_DAZScore),]
  outcome_dt <- unique(outcome_dt)
  
  setnames(outcome_dt, c("StudyID", "childTimepointAge_days"), c("subject_id", "outcome_age_days"))
  
  outcome_dt <- outcome_dt[!is.na(outcome_age_days),]
  outcome_dt[, outcome_age_days := as.numeric(outcome_age_days)]
  
  outcome_dt[,outcome_binary := ifelse(gsed_LongForm_DAZScore < -1, 0,1)]
  
  outcome_dt <- as.data.table(outcome_dt)
  
  ## SELECT EARLIEST MRI + EEG, LATEST OUTCOME AND MERGE--------------------------------------------------------------------------
  #Subset to earliest time point - MRI
  length(unique(mri_wide$subject_id))
  mri_wide <- mri_wide[mri_wide[, .I[which.min(mri_age_days)], by = subject_id]$V1]
  length(unique(mri_wide$subject_id))
  
  #Subset to letest time point - Outcome
  length(unique(outcome_dt$subject_id))
  outcome_dt <- outcome_dt[outcome_dt[, .I[which.max(outcome_age_days)], by = subject_id]$V1]
  length(unique(outcome_dt$subject_id))
  
  #Merge
  dt_merged <- merge(mri_wide, outcome_dt, by=c("subject_id", "CohortName"))
  
  #Ensure MRI/EEG is earlier than outcome
  dt_merged <- dt_merged[mri_age_days<=outcome_age_days, ]
}

if (dataset_type == "khula_mri_eeg") {
  ##KHULA - ORIG, NORM, CENTILES -----------------------------------------------------
  #Set objects
  other_vars <- c('timepointFamilySize', 
                  'childBirthWeight_kgs', 'childBirthLength_inches', 
                  'MaternalEducation_schoolingYears', 'maternalAgeAtChildBirth_years', 'EverStunted', 'EverWasted')
  # outfile_name <- "Khula_eeg_sf_oth_gsed_outcome.csv"
  # outfile_name <- "Khula_eeg_sf_oth_bsid_CCS_outcome.csv"
  # outfile_name <- "Khula_eeg_sf_oth_bsid_LCS_outcome.csv"
  outfile_name <- "Khula_eeg_sf_oth_bsid_MCS_outcome.csv"
  
  # outcome_var_name <- "gsed_LongForm_DAZScore"
  # outcome_var_name <- "bsid_CCS"
  # outcome_var_name <- "bsid_LCS"
  outcome_var_name <- "bsid_MCS"
  
  # PREP MRI -- 
  mri_dt <- as.data.table(read.csv(paste0(read_dir_centiles, "mri_centiles_sexstratified.csv")))
  
  #Convert to wide
  mri_dt <- subset(mri_dt, select = c(age_days, measure, subject_id, orig_value, child_sex, value, centile))
  mri_long <- pivot_longer(
    data = mri_dt,
    cols = c(orig_value, value, centile),
    names_to = "metric",
    values_to = "metric_value"
  )
  
  mri_long$metric[mri_long$metric == "value"] <- "norm_value"
  mri_long$measure_metric <- paste0(mri_long$measure, "_", mri_long$metric)
  
  mri_selected <- mri_long[, c("subject_id", "age_days", "child_sex", "measure_metric", "metric_value")]
  
  mri_wide <- pivot_wider(
    data = mri_selected,
    names_from = measure_metric,
    values_from = metric_value
  )
  
  mri_wide <- as.data.table(mri_wide)
  
  setnames(mri_wide, c("age_days"), c("mri_age_days"))
  
  # PREP EEG -- 
  eeg_power_dt <- as.data.table(read.csv(paste0(read_dir_centiles, "eeg_power_centiles.csv")))
  eeg_total_dt <- as.data.table(read.csv(paste0(read_dir_centiles, "eeg_total_centiles.csv")))
  eeg_micro_occurrence_dt <- as.data.table(read.csv(paste0(read_dir_centiles, "eeg_micro_occur_centiles.csv")))
  
  eeg_power_dt <- eeg_power_dt[!brain_region == "Temporal", ]
  eeg_power_dt <- eeg_power_dt[!measure == "Theta", ]
  eeg_micro_occurrence_dt <- eeg_micro_occurrence_dt[measure == "occurenceD"]
  eeg_total_dt <- eeg_total_dt[measure == "avg_alpha_by_chan"]
  
  eeg_dt <- rbind(eeg_power_dt, eeg_total_dt)
  eeg_dt <- rbind(eeg_dt, eeg_micro_occurrence_dt)
  
  #Convert to wide
  eeg_dt <- subset(eeg_dt, select = c(eeg_age_days, measure, brain_region, subject_id, value, centile))
  
  eeg_long <- pivot_longer(
    data = eeg_dt,
    cols = c(value, centile),
    names_to = "metric",
    values_to = "metric_value"
  )
  
  eeg_long$metric[eeg_long$metric == "value"] <- "orig_value"
  
  eeg_long$measure_metric <- paste0(eeg_long$brain_region, "_", eeg_long$measure, "_", eeg_long$metric)
  
  eeg_selected <- eeg_long[, c("subject_id", "eeg_age_days", "measure_metric", "metric_value")]
  
  eeg_selected <- summarise(
    group_by(
      eeg_selected,
      subject_id,
      eeg_age_days,
      measure_metric
    ),
    metric_value = mean(metric_value, na.rm = TRUE)
  )
  
  eeg_wide <- pivot_wider(
    data = eeg_selected,
    names_from = measure_metric,
    values_from = metric_value
  )
  
  eeg_wide <- as.data.table(eeg_wide) 

  #Create binary alpha frequency variable
  eeg_wide[, alpha_peak_binary := ifelse(is.na(All_avg_alpha_by_chan_orig_value), 0, 1)]
  eeg_wide <- subset(eeg_wide, select = -c(All_avg_alpha_by_chan_centile, All_avg_alpha_by_chan_orig_value))

  eeg_wide <- eeg_wide[complete.cases(eeg_wide)]
  
  # PREP OUTCOME --  
  outcome_dt <- as.data.table(read.csv(paste0(read_dir_flywheel, "UNITY-dataTemplate_MRI_clean.csv")))
  outcome_dt <- outcome_dt[CohortName == "Khula", ]
  
  outcome_dt <- outcome_dt[, .SD, .SDcols = c("StudyID", "childTimepointAge_days", outcome_var_name)]
  outcome_dt <- outcome_dt[!is.na(get(outcome_var_name)),]
  outcome_dt <- unique(outcome_dt)
  
  setnames(outcome_dt, c("StudyID", "childTimepointAge_days"), c("subject_id", "outcome_age_days"))
  
  outcome_dt <- outcome_dt[!is.na(outcome_age_days),]
  outcome_dt[, outcome_age_days := as.numeric(outcome_age_days)]
  
  if (outcome_var_name %like% "bsid") {
    outcome_dt[,outcome_binary := ifelse(get(outcome_var_name) < 85, 0,1)]
  } else {
    outcome_dt[,outcome_binary := ifelse(get(outcome_var_name) < -1, 0,1)]
  }
  
  outcome_dt$subject_id <- gsub("_", "-", outcome_dt$subject_id)
  outcome_dt <- as.data.table(outcome_dt)

  ## SELECT EARLIEST MRI + EEG, LATEST OUTCOME AND MERGE--------------------------------------------------------------------------
  #Subset to earliest time point - MRI
  length(unique(mri_wide$subject_id))
  mri_wide <- mri_wide[mri_wide[, .I[which.min(mri_age_days)], by = subject_id]$V1]
  length(unique(mri_wide$subject_id))
  
  #Subset to earliest time point - EEG
  length(unique(eeg_wide$subject_id))
  eeg_wide <- eeg_wide[eeg_wide[, .I[which.min(eeg_age_days)], by = subject_id]$V1]
  length(unique(eeg_wide$subject_id))
  
  #Subset to letest time point - Outcome
  length(unique(outcome_dt$subject_id))
  outcome_dt <- outcome_dt[outcome_dt[, .I[which.max(outcome_age_days)], by = subject_id]$V1]
  length(unique(outcome_dt$subject_id))
  
  #Merge
  dt_merged <- merge(mri_wide, outcome_dt, by=c("subject_id"))
  dt_merged <- merge(dt_merged, eeg_wide, by=c("subject_id"))
  
  #Ensure MRI/EEG is earlier than outcome
  dt_merged <- dt_merged[mri_age_days<=outcome_age_days&eeg_age_days<=outcome_age_days, ]
  
  dt_merged[, CohortName := "Khula"]
  dt_merged$subject_id <- gsub("-", "_", dt_merged$subject_id)
}


## PREP OTHER FLYWHEEL VARS -------------------------------------------------------------
other_var_dt <- as.data.table(read.csv(paste0(read_dir_flywheel, "UNITY-dataTemplate_MRI_clean.csv")))

other_var_dt <- process_unknowns(other_var_dt)
other_var_dt <- other_var_dt[, .SD, .SDcols = c('CohortName', 'StudyID', other_vars)]
other_var_dt <- other_var_dt[CohortName %in% dt_merged$CohortName,]

other_var_dt <- other_var_dt[complete.cases(other_var_dt)]
other_var_dt <- unique(other_var_dt)

setnames(other_var_dt, c("StudyID"), c("subject_id"))


## MERGE AND READ OUT -------------------------------------------------------------
dt_merged_all <- merge(dt_merged, other_var_dt, by=c("subject_id", "CohortName"))

write.csv(dt_merged_all, paste0(save_dir, outfile_name), row.names=FALSE)


