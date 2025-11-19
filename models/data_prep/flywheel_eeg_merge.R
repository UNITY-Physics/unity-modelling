######################################################################
## Merge Khula and SPACE EEG data to MRI FlyWheel dataset
## Brain Health Metrics 
## IHME
## September 2025
#####################################################################


## Set up ---------------------------------------------------------------------------------------------------
rm(list=ls())
username <- Sys.getenv("USER")

library(data.table)
library(dplyr)
library(readxl)

#System configuration
l_root <- "/FILEPATH/"

#Set directories
read_dir_flywheel_gf <- paste0(l_root, "FILEPATH/")
read_dir_flywheel_london <- paste0(l_root, "FILEPATH/")
read_dir <- paste0(l_root, "FILEPATH/")

save_dir <- paste0(l_root, "FILEPATH/")

#Source functions
code_dir <- paste0("/FILEPATH/", username, "/FILEPATH/")
source(paste0(code_dir, "khula_functions.R"))


## Data manipulation for Khula ---------------------------------------------------------------------------------------------------
#read in data
khula_eeg <- as.data.table(read.csv(paste0(save_dir, "khula_subset_all_eeg_selected_vars.csv")))

khula_dt <- as.data.table(read_xlsx(paste0(read_dir, "CHARACTERIZING_DEVELOPING_EXECUTIVE_FUNCTIONS_IN_THE_FIRST_1000_DAYS_IN_SOUTH_AFRICA_AND_MALAWI_THE_KHULA_STUDY_AGGREGATED_Y2025M008D08.XLSX")))

flywheel_dt_gf <- as.data.table(read.csv(paste0(read_dir_flywheel_gf, "UNITY-dataTemplate_MRIEEG_clean.csv")))

#separating Flywheel dt
flywheel_dt_to_merge <- flywheel_dt[CohortName == "Khula", ]
flywheel_dt_to_keep <- flywheel_dt[!CohortName == "Khula", ]

#deal with duplicate MRI rows by updating age
duplicates1 <- flywheel_dt_to_merge[duplicated(flywheel_dt_to_merge[, .(StudyID, studyTimepoint)]) | duplicated(flywheel_dt_to_merge[, .(StudyID, studyTimepoint)], fromLast = TRUE)]
duplicates2 <- flywheel_dt_to_merge[duplicated(flywheel_dt_to_merge[, .(StudyID, childTimepointAge_days)]) | duplicated(flywheel_dt_to_merge[, .(StudyID, childTimepointAge_days)], fromLast = TRUE)]
write.csv(duplicates2, paste0(save_dir, "flywheel_khula_duplicates_StudyID_age.csv"))

#if duplicates are within a week, retain later scan and copy over all time-invariant variables - many StudyIDs w duplicates also have other timepoints without duplicates

#if duplicates are more than a week apart, wait for clarification

#merging on timepoint from Khula file
khula_dt_mri <- clean_mri(khula_dt)
khula_dt_mri <- subset(khula_dt_mri, select = c(subject_id, hyp_age, time))
setnames(khula_dt_mri, c("hyp_age", "subject_id"), c("childTimepointAge_days", "StudyID"))
khula_dt_mri$StudyID <- gsub("-", "_", khula_dt_mri$StudyID)

#flywheel_dt_to_merge[, childTimepointAge_days := round(childTimepointAge_days)]

#flywheel_dt_to_merge_with_timepoint <- merge(flywheel_dt_to_merge[!is.na(childTimepointAge_days),], khula_dt_mri, by = c("childTimepointAge_days", "StudyID"))

## Data manipulation for SPACE ---------------------------------------------------------------------------------------------------
#open time point mapping files
time_map_mri <- read.csv(paste0(read_dir_flywheel_london, "SPACEEXPLORE-HyperfineDates_08.09.csv"))
time_map_eeg_orig <- read.xlsx(paste0(read_dir_flywheel_london, "SPACE_longitudinal_EEG_withtimepoint.xlsx"))
time_map_eeg_orig <- unique(time_map_eeg_orig)

#subset master file to space and delete EEG variables
space_dt <- unique(flywheel_dt[CohortName == "SPACE_Explore", ])
cols_to_remove <- grep("power|alpha|eeg", names(space_dt), value = TRUE)
space_dt[, (cols_to_remove) := NULL]

#if age in days is NA, calculate mri age in days from age in month variable for master file
space_dt[is.na(childTimepointAge_days) & !is.na(childTimepointAge_months), childTimepointAge_days := childTimepointAge_months*30.436875]

#if age in days is still NA, pull all data for those IDs and see if possible to calculate age from a different timepoint
space_dt_keep <- space_dt[!is.na(childTimepointAge_days), ]
space_dt_fix <- space_dt[is.na(childTimepointAge_days), ]
ids_prob <- unique(space_dt_fix$StudyID) 
all_rows_fix <- space_dt[StudyID %in% ids_prob,] #these are unique

#harmonize map file variable name and structure in map to master file
time_map_mri$subject_id <- gsub("-", "_", time_map_mri$subject_id)
setnames(time_map_mri, c("subject_id", "hyperfine_date", "hyperfine_child_age"), c("StudyID", "studyTimepoint", "childTimepointAge_days"))
date_obj <- as.Date(time_map_mri$studyTimepoint, format = "%Y-%m-%d")
time_map_mri$studyTimepoint <- format(date_obj, "%d/%m/%Y")

#remove duplicate row in mri map file
time_map_mri <- time_map_mri[!(time_map_mri$StudyID=="XXX" & time_map_mri$timepoint=="t2"), ]

#merge MRI time point variable to master - calculate timepoint if missing
space_dt2 <- merge(space_dt, time_map_mri, by = c("StudyID", "studyTimepoint", "childTimepointAge_days"), all.x=TRUE)
space_dt2[is.na(timepoint) & childTimepointAge_days<182, timepoint := "t1"]
space_dt2[is.na(timepoint) & childTimepointAge_days>=182, timepoint := "t2"]

#Deal with NA timepoints if childTimepointAge_days is NA by remerging only on StudyID and studyTimepoint
space_dt2_keep <- space_dt2[!is.na(timepoint), ]
space_dt2_fix <- space_dt2[is.na(timepoint), ]
ids_prob <- unique(space_dt2_fix$StudyID) 
all_rows_fix <- space_dt2[StudyID %in% ids_prob,] 

#manually add in age in days based on other rows of data with age and study timepoint
space_dt2_fix[StudyID=="XXX", childTimepointAge_days := 85]
space_dt2_fix[StudyID=="XXX", childTimepointAge_days := 80]
space_dt2_fix[is.na(timepoint) & childTimepointAge_days<182, timepoint := "t1"]
space_dt2_fix[is.na(timepoint) & childTimepointAge_days>=182, timepoint := "t2"]

#manually add in timepoint  based on only row of eeg data that does have age and timepoint (age = 86)
space_dt2_fix[StudyID=="XXX", timepoint := "t1"]

#bind fixed rows back on
space_dt3 <- rbind(space_dt2_keep, space_dt2_fix)

#merge eeg map time point to most up-to-date SPACE EEG file Lauryn created
space_up_to_date <- read.csv(paste0(save_dir, "space_subset_all_eeg_selected_vars.csv"))
time_map_eeg <- merge(space_up_to_date, time_map_eeg_orig[,c("subject_id", "eeg_age_days", "timepoint")], by = c("subject_id", "eeg_age_days"))

#merge EEG variables to master including timepoint variable - first harmonize variable name and structure in map to master file
setnames(time_map_eeg, "subject_id", "StudyID")
time_map_eeg$timepoint <- gsub("T", "t", time_map_eeg$timepoint)
time_map_eeg$StudyID <- gsub("-", "_", time_map_eeg$StudyID)
setDT(time_map_eeg)

#check NAs in timepoint variable in merge - confirm whether these are unique individuals or a merge problem and update timepoint based on age
na_dt = time_map_eeg[is.na(timepoint), ]
ids <- unique(na_dt$StudyID)
na_dt_subj <- time_map_eeg[StudyID %in% ids,] #seems to be a lack of timepoint, not a merge problem

#create timepoint variable where NA based on logic of under 6 months or above 6 months and confirm this doesn't create duplicate timepoints for a subject
time_map_eeg[is.na(timepoint) & eeg_age_days<182, timepoint := "t1"]
time_map_eeg[is.na(timepoint) & eeg_age_days>=182, timepoint := "t2"]

#check no missingness in timpoint
na_dt = time_map_eeg[is.na(timepoint), ]
ids <- unique(na_dt$StudyID)
na_dt_subj <- time_map_eeg[StudyID %in% ids,] #confirm 0

#one triplicate row of data for one subject - XXX - for now retain only the row that has complete data that follows age logic
time_map_eeg <- time_map_eeg[!(StudyID=="XXX" & (eeg_age_days==272 | eeg_age_days==171)), ]

#XXX: no MRI for T1 but two EEG timepoints that are very different from each other, one with total and one with spectral power
#Solution - tag as T2 bc almost identical timing to MRI data tagged as T2
time_map_eeg[StudyID=="XXX" & eeg_age_days==178, timepoint=="t2"]

#for the merge, retain EEG rows where there is no MRI for the timepoint
space_dt4 <- merge(space_dt3, time_map_eeg, by = c("StudyID", "timepoint"), all=TRUE)
space_dt4[StudyID=="XXX" & is.na(timepoint) & is.na(childTimepointAge_days), childTimepointAge_days := 163]

#need to add time invariant variables to the added EEG rows that didn't merge to MRI
merged_eeg_only <- space_dt4[is.na(childTimepointAge_days), ]
merged_rest <- space_dt4[!is.na(childTimepointAge_days), ]

vars_to_fill <- c( "CohortName", "CohortLocation_city", "CohortLocation_country",
                  "UniqueStudyID", "childBiologicalSex", "childGestation_days",
                  "childGestation_weeks", "childBirthWeight_lbs", "ChildBirthWeight_kgs",
                  "childBirthLength_inches", "maternalAgeAtChildBirth_years", "MaternalEducation_HHS",
                  "MaternalEducation_schoolingYears", "childBirthOrdinal", "childCountryOfBirth",
                  "childCityOfBirth", "childRace", "childEthnicity", "childReligion", "childCaste",
                  "EverStunted", "EverWasted")

#determine which variables are in SPACE and which are fully NA
all_na_cols <- names(merged_rest)[colSums(is.na(merged_rest)) == nrow(merged_rest)]

vars_to_fill_space <- c("CohortName", "CohortLocation_city", "CohortLocation_country",
                      "UniqueStudyID", "childBiologicalSex", "childGestation_days",
                      "childGestation_weeks", "childBirthWeight_lbs", "childBirthWeight_kgs",
                      "childBirthLength_inches", "maternalAgeAtChildBirth_years",
                      "MaternalEducation_schoolingYears", "childCountryOfBirth",
                      "childCityOfBirth", "EverStunted", "EverWasted")

#Remove specified columns from EEG only dt, retain for filled dt and make unique to StudyID and timepoint
merged_eeg_only = merged_eeg_only[, -..vars_to_fill_space]
ids <- unique(merged_eeg_only$StudyID)
oth_var_dt <- merged_rest[StudyID %in% ids, ]
vars_to_fill_space2 <- c(vars_to_fill_space, "StudyID")
oth_var_dt <- unique(oth_var_dt[, ..vars_to_fill_space2])

#Merge additional variables and see how many remaining as all NA
merged_eeg_only2 <- merge(oth_var_dt, merged_eeg_only, by=c("StudyID"), all=TRUE)
merged_rest_add <- merged_eeg_only2[!is.na(CohortName), ]
merged_eeg_only2 <- merged_eeg_only2[is.na(CohortName), ]

merged_rest2 <- rbind(merged_rest, merged_rest_add, fill = TRUE)
merged_interim <- rbind(merged_rest2, merged_eeg_only2, fill = TRUE)

#Add a few additional invariant variables
merged_interim[is.na(mriCollectedAtTimepoint), mriCollectedAtTimepoint := "No"]
merged_interim[is.na(CohortName), CohortName := "SPACE_Explore"]
merged_interim[is.na(CohortLocation_country), CohortLocation_country := "South_Africa"]
merged_interim[is.na(CohortLocation_city), CohortLocation_city := "Cape_Town"]
merged_interim[is.na(UniqueStudyID), UniqueStudyID := paste("SPACE_Explore_",StudyID)]

#Remove subject with multiple rows with incorrect logic - retain correct row based on discussion
merged_interim <- merged_interim[!(StudyID=="XXX" & (eeg_age_days==172 | eeg_age_days==272)), ]

#Update remaining duplicates (2 StudyIDs)
#XXX: appears to have two MRI timpoints about 2 months apart but both time1 - one EEG matched to both but close to the younger age
#Solution - only match EEG to the younger age and retain other row of MRI

#XXX: two MRI timepoints and no EEG; one of the MRI rows has very very different volume values than the other row but has demographics other row doesn't have
#Solution - flag to Sean




write.csv(merged_interim, paste0(save_dir, "space_updated_to_merge_eeg_interim.csv"), row.names=FALSE)
          
#vars to manual: mriCollectedAtTimepoint

## PART 2

######################################################################
## Merge Khula and SPACE EEG data to MRI FlyWheel dataset
## Brain Health Metrics 
## IHME
## September 2025
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
read_dir <- paste0(l_root, "FILEPATH/")

save_dir <- paste0(l_root, "FILEPATH/")

#Source functions
code_dir <- paste0("/FILEPATH/", username, "/FILEPATH/")
source(paste0(code_dir, "khula_functions.R"))

#Read in FlyWheel master dataset
flywheel_dt <- as.data.table(read.csv(paste0(read_dir_flywheel, "UNITY-dataTemplate_MRIEEG_clean.csv")))

## Data manipulation for Khula ---------------------------------------------------------------------------------------------------
#read in data
khula_eeg <- as.data.table(read.csv(paste0(save_dir, "khula_subset_all_eeg_selected_vars.csv")))

khula_dt <- as.data.table(read_xlsx(paste0(read_dir, "CHARACTERIZING_DEVELOPING_EXECUTIVE_FUNCTIONS_IN_THE_FIRST_1000_DAYS_IN_SOUTH_AFRICA_AND_MALAWI_THE_KHULA_STUDY_AGGREGATED_Y2025M008D08.XLSX")))

#separating Flywheel dt
flywheel_dt_to_merge <- flywheel_dt[CohortName == "Khula", ]
flywheel_dt_to_keep <- flywheel_dt[!CohortName == "Khula", ]

flywheel_dt_to_merge[StudyID == "XXX", childBiologicalSex := "Female"]
flywheel_dt_to_merge[StudyID == "XXX", childBiologicalSex := "Female"]
flywheel_dt_to_merge[StudyID == "XXX", childBiologicalSex := "Female"]
flywheel_dt_to_merge[StudyID == "XXX", childBiologicalSex := "Female"]

flywheel_dt_to_merge <- unique(flywheel_dt_to_merge)

#deal with duplicate MRI rows by updating age
duplicates <- flywheel_dt_to_merge[duplicated(flywheel_dt_to_merge[, .(StudyID, childTimepointAge_days)]) | duplicated(flywheel_dt_to_merge[, .(StudyID, childTimepointAge_days)], fromLast = TRUE)]
duplicates <- duplicates[!is.na(childTimepointAge_days),]
write.csv(duplicates, paste0(save_dir, "flywheel_khula_duplicates_StudyID_age.csv"))

flywheel_dt_to_merge_dup_fix <- setdiff(flywheel_dt_to_merge, duplicates)

#if duplicates are within a week, retain later scan and copy over all time-invariant variables - many StudyIDs w duplicates also have other timepoints without duplicates
#if duplicates are more than a week apart, keep earlier timepoint only
duplicates_fixed <- as.data.table(read.csv(paste0(save_dir, "flywheel_khula_duplicates_StudyID_age_fixed.csv")))

flywheel_dt_to_merge_dup_fix <- rbind(flywheel_dt_to_merge_dup_fix, duplicates_fixed)

#merging on MRI timepoint from Khula file
khula_dt_mri <- clean_mri(khula_dt)
khula_dt_mri <- subset(khula_dt_mri, select = c(subject_id, hyp_age, time))
setnames(khula_dt_mri, c("hyp_age", "subject_id"), c("childTimepointAge_days_round", "StudyID"))
khula_dt_mri$StudyID <- gsub("-", "_", khula_dt_mri$StudyID)
khula_dt_mri <- unique(khula_dt_mri)
khula_dt_mri <- as.data.table(khula_dt_mri)
khula_dt_mri[, childTimepointAge_days_round := round(childTimepointAge_days_round / 10) * 10]

flywheel_dt_to_merge_dup_fix[, childTimepointAge_days_round := round(childTimepointAge_days / 10) * 10]

flywheel_dt_to_merge_dup_fix_with_timepoint <- merge(flywheel_dt_to_merge_dup_fix, khula_dt_mri, by = c("childTimepointAge_days_round", "StudyID"), all.x = TRUE)

test <- flywheel_dt_to_merge_dup_fix_with_timepoint[is.na(time) & !is.na(childTimepointAge_days), ]

flywheel_dt_to_merge_dup_fix_with_timepoint[is.na(time) & !is.na(childTimepointAge_days) & childTimepointAge_months < 6, time := '2-5 mo']
flywheel_dt_to_merge_dup_fix_with_timepoint[is.na(time) & !is.na(childTimepointAge_days) & childTimepointAge_months >= 6 & childTimepointAge_months < 12, time := '6-11 mo']
flywheel_dt_to_merge_dup_fix_with_timepoint[is.na(time) & !is.na(childTimepointAge_days) & childTimepointAge_months >= 12 & childTimepointAge_months < 16, time := '12-15 mo']
flywheel_dt_to_merge_dup_fix_with_timepoint[is.na(time) & !is.na(childTimepointAge_days) & childTimepointAge_months >= 16 & childTimepointAge_months < 20, time := '16-19 mo']
flywheel_dt_to_merge_dup_fix_with_timepoint[is.na(time) & !is.na(childTimepointAge_days) & childTimepointAge_months >= 20, time := '20-24 mo']

test <- flywheel_dt_to_merge_dup_fix_with_timepoint[is.na(time) & !is.na(childTimepointAge_days), ]
setcolorder(flywheel_dt_to_merge_dup_fix_with_timepoint, c("time", setdiff(names(flywheel_dt_to_merge_dup_fix_with_timepoint), "time")))

duplicates2 <- flywheel_dt_to_merge_dup_fix_with_timepoint[duplicated(flywheel_dt_to_merge_dup_fix_with_timepoint[, .(StudyID, time)]) | duplicated(flywheel_dt_to_merge_dup_fix_with_timepoint[, .(StudyID, time)], fromLast = TRUE)]
duplicates2 <- duplicates2[!is.na(childTimepointAge_days),]
setcolorder(duplicates2, c("time", setdiff(names(duplicates2), "time")))
unique(duplicates2$StudyID)

flywheel_dt_to_merge_dup_fix_with_timepoint[StudyID == "XXX" & childTimepointAge_days_round == 180, `:=` (time = '6-11 mo', childTimepointAge_days = 330, childTimepointAge_months = 330/30.44)]
flywheel_dt_to_merge_dup_fix_with_timepoint[StudyID == "XXX" & childTimepointAge_days_round == 480, `:=` (time = '16-19 mo')]
flywheel_dt_to_merge_dup_fix_with_timepoint[StudyID == "XXX" & childTimepointAge_days_round == 480, `:=` (time = '16-19 mo')]
flywheel_dt_to_merge_dup_fix_with_timepoint[StudyID == "XXX" & childTimepointAge_days_round == 490, `:=` (time = '12-15 mo')]
flywheel_dt_to_merge_dup_fix_with_timepoint[StudyID == "XXX" & childTimepointAge_days_round == 490, `:=` (time = '12-15 mo')]
flywheel_dt_to_merge_dup_fix_with_timepoint[StudyID == "XXX" & childTimepointAge_days_round == 550, `:=` (time = '12-15 mo')]

duplicates3 <- flywheel_dt_to_merge_dup_fix_with_timepoint[duplicated(flywheel_dt_to_merge_dup_fix_with_timepoint[, .(StudyID, time)]) | duplicated(flywheel_dt_to_merge_dup_fix_with_timepoint[, .(StudyID, time)], fromLast = TRUE)]
duplicates3 <- duplicates3[!is.na(childTimepointAge_days),]
unique(duplicates3$StudyID)

#fixing EEG duplicates
khula_eeg <- khula_eeg[!(subject_id == "XXX" & time == '2-5 mo' & eeg_age_days == 144),]
khula_eeg <- khula_eeg[!(subject_id == "XXX" & time == '6-11 mo' & eeg_age_days == 213),]

duplicates4 <- khula_eeg[duplicated(khula_eeg[, .(subject_id, time)]) | duplicated(khula_eeg[, .(subject_id, time)], fromLast = TRUE)]
khula_eeg_no_dup <- setdiff(khula_eeg, duplicates4)
duplicates4_keep <- duplicates4[is.na(occurenceD), ]
duplicates4_drop <- duplicates4[!is.na(occurenceD), ]
duplicates4_drop <- subset(duplicates4_drop, select = c(subject_id, time, eeg_age_days, occurenceD))
duplicates4_keep <- subset(duplicates4_keep, select = -c(occurenceD))
duplicates4_keep <- merge(duplicates4_keep, duplicates4_drop, by = c("subject_id", "time"))
duplicates4_keep[, eeg_age_days := round((eeg_age_days.x+eeg_age_days.y)/2)]
duplicates4_keep <- subset(duplicates4_keep, select = -c(eeg_age_days.x, eeg_age_days.y))

khula_eeg_no_dup <- rbind(khula_eeg_no_dup, duplicates4_keep)

#for the merge, retain EEG rows where there is no MRI for the timepoint
setnames(khula_eeg_no_dup, c("subject_id", "time"), c("StudyID", "timepoint"))
khula_eeg_no_dup$StudyID <- gsub("-", "_", khula_eeg_no_dup$StudyID)
setnames(flywheel_dt_to_merge_dup_fix_with_timepoint, c("time"), c("timepoint"))
flywheel_dt_to_merge_dup_fix_with_timepoint <- subset(flywheel_dt_to_merge_dup_fix_with_timepoint, select = -c(childTimepointAge_days_round))

cols_to_remove <- grep("power|alpha|eeg", names(flywheel_dt_to_merge_dup_fix_with_timepoint), value = TRUE)
flywheel_dt_to_merge_dup_fix_with_timepoint[, (cols_to_remove) := NULL]
flywheel_dt_to_merge_with_eeg <- merge(flywheel_dt_to_merge_dup_fix_with_timepoint, khula_eeg_no_dup, by = c("StudyID", "timepoint"), all=TRUE)

#need to add time invariant variables to the added EEG rows that didn't merge to MRI
merged_eeg_only <- flywheel_dt_to_merge_with_eeg[is.na(CohortName), ]
merged_rest <- flywheel_dt_to_merge_with_eeg[!is.na(CohortName), ]

vars_to_fill <- c( "CohortName", "CohortLocation_city", "CohortLocation_country",
                   "UniqueStudyID", "childBiologicalSex", "childGestation_days",
                   "childGestation_weeks", "childBirthWeight_lbs", "childBirthWeight_kgs",
                   "childBirthLength_inches", "maternalAgeAtChildBirth_years", "MaternalEducation_HHS",
                   "MaternalEducation_schoolingYears", "childBirthOrdinal", "childCountryOfBirth",
                   "childCityOfBirth", "childRace", "childEthnicity", "childReligion", "childCaste",
                   "EverStunted", "EverWasted")
"childTimepointHeight_inches"             
#determine which variables are in Khula and which are fully NA
all_na_cols <- names(merged_rest)[colSums(is.na(merged_rest)) == nrow(merged_rest)]

vars_to_fill_khula <- c( "CohortName", "CohortLocation_city", "CohortLocation_country",
                         "UniqueStudyID", "childBiologicalSex", "childGestation_days",
                         "childGestation_weeks", "childBirthWeight_lbs", "childBirthWeight_kgs",
                         "childBirthLength_inches", "maternalAgeAtChildBirth_years",
                         "MaternalEducation_schoolingYears", "childBirthOrdinal", "childCountryOfBirth",
                         "childCityOfBirth", "EverStunted", "EverWasted")

#Remove specified columns from EEG only dt, retain for filled dt and make unique to StudyID and timepoint
merged_eeg_only = merged_eeg_only[, -..vars_to_fill_khula]
ids <- unique(merged_eeg_only$StudyID)
oth_var_dt <- merged_rest[StudyID %in% ids, ]
vars_to_fill_khula2 <- c(vars_to_fill_khula, "StudyID")
oth_var_dt <- unique(oth_var_dt[, ..vars_to_fill_khula2])

duplicates5 <- oth_var_dt[duplicated(oth_var_dt[, .(StudyID)]) | duplicated(oth_var_dt[, .(StudyID)], fromLast = TRUE)]
duplicates5_fixed <- duplicates5[!is.na(childGestation_days),]

oth_var_dt <- oth_var_dt[!StudyID %in% duplicates5$StudyID, ]
oth_var_dt <- rbind(oth_var_dt, duplicates5_fixed)

#Merge additional variables and see how many remaining as all NA
merged_eeg_only2 <- merge(oth_var_dt, merged_eeg_only, by=c("StudyID"), all=TRUE)
merged_rest_add <- merged_eeg_only2[!is.na(CohortName), ]
merged_eeg_only2 <- merged_eeg_only2[is.na(CohortName), ]
merged_eeg_only2[, `:=` (CohortName = "Khula", CohortLocation_city = "Cape_Town", CohortLocation_country = "South_Africa", UniqueStudyID = paste0("Khula_", StudyID))]

merged_rest2 <- rbind(merged_rest, merged_rest_add)

merged_interim <- rbind(merged_rest2, merged_eeg_only2)
write.csv(merged_interim, paste0(save_dir, "khula_updated_to_merge_eeg_interim.csv"), row.names=FALSE)


## Combine with SPACE and rest of FlyWheel file ---------------------------------------------------------------------------------------------------
space_dt <- as.data.table(read.csv(paste0(save_dir, "space_updated_to_merge_eeg_interim.csv")))

flywheel_dt_to_keep <- flywheel_dt[!CohortName %in% c("Khula", "SPACE_Explore"), ]
cols_to_remove <- grep("power|alpha|eeg", names(flywheel_dt_to_keep), value = TRUE)
flywheel_dt_to_keep[, (cols_to_remove) := NULL]

flywheel_dt_to_keep <- rbind(flywheel_dt_to_keep, space_dt, fill = TRUE)
flywheel_dt_to_keep <- rbind(flywheel_dt_to_keep, merged_interim, fill = TRUE)

## fix terminal GSED, Bayley ---------------------------------------------------------------------------------------------------
terminal_dt <- as.data.table(read.csv(paste0(save_dir, "flywheel_terminal_gsed_bsid.csv")))
setnames(terminal_dt, c("terminal_gsed_LongForm_DAZScore_Age_days", "terminal_gsed_LongForm_DAZScore", "terminal_bsid_CCS_Age_days",
                        "terminal_bsid_CCS", "terminal_bsid_LCS_Age_days", "terminal_bsid_LCS", "terminal_bsid_MCS_Age_days", "terminal_bsid_MCS"),
         c("MAX_gsed_LongForm_DAZScore_Age_days", "MAX_gsed_LongForm_DAZScore", "MAX_bsid_CCS_Age_days",
           "MAX_bsid_CCS", "MAX_bsid_LCS_Age_days", "MAX_bsid_LCS", "MAX_bsid_MCS_Age_days", "MAX_bsid_MCS"))

flywheel_dt_to_keep <- subset(flywheel_dt_to_keep, select = -c(MAX_gsed_LongForm_DAZScore, MAX_bsid_CCS, MAX_bsid_LCS, MAX_bsid_MCS))
flywheel_dt_to_keep <- merge(flywheel_dt_to_keep, terminal_dt, by = c("CohortName", "CohortLocation_city", "CohortLocation_country", "StudyID","UniqueStudyID"), all.x=TRUE)

## write out new master FlyWheel file ---------------------------------------------------------------------------------------------------
flywheel_dt_to_keep[timepoint == "2-5 mo", timepoint := 't1']
flywheel_dt_to_keep[timepoint == "6-11 mo", timepoint := 't2']
flywheel_dt_to_keep[timepoint == "12-15 mo", timepoint := 't3']
flywheel_dt_to_keep[timepoint == "16-19 mo", timepoint := 't4']
flywheel_dt_to_keep[timepoint == "20-24 mo", timepoint := 't5']

write.csv(flywheel_dt_to_keep, paste0(save_dir, "UNITY-dataTemplate_MRIEEG_clean.csv"), row.names=FALSE)

