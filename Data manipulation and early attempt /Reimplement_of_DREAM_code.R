setwd("~/Desktop/GuanLab/")
# source data file
core <- read.csv("rawdata/prostate_cancer_challenge_data_training/CoreTable_training.csv",
                 header = T, stringsAsFactors = F)


# re-implement creat_DEATH_target.pl
death <- core[, c("RPT", "DEATH")]
death$DEATH <- vapply(death$DEATH, function(x) ifelse(x == "YES", 1, 0), numeric(1))

# re-implement prepare_habini.pl
# Count the average value of ECOG_C
# It is used to fill the missing value of this column
core[, 26] <- as.numeric(core[, 26])
ECOG_C <- sum(core[, 26], na.rm = T)
ECOG_C_count <- sum(!is.na(core[, 26]))
ECOG_C <- ECOG_C/ECOG_C_count

# Calculate the average values of these columns, also for filling the missing value
core[, 30:50] <- lapply(core[, 30:50], function(x) as.numeric(x))
ALP <- colSums(core[, 30:50], na.rm = T)
ALP_count <- unlist(lapply(core[, 30:50], function(x) sum(!is.na(x))))
ALP_avg <- ALP / ALP_count

# Columns for Halabi
core_habini <- core[, c(3, 83, 36, 57, 26, 47, 35, 39, 30, 59)]
core_habini$ANALGESICS <- ifelse(core_habini$ANALGESICS == "YES", 1, 0)
core_habini$LDH <- vapply(core_habini$LDH, 
                          function(x) {
                            if (is.na(x)) {return(1)} # avg for LDH is 263.9051, return 1
                            ifelse(x > 250, 1, 0)
                          }, numeric(1))
# for records in core[, c(57, 58, 60, 61, 62, 63, 66, 67, 69, 72, 73, 75)]
# if BONE is empty & other columns are also empty: 0
# if BONE have record & other columns are empty: 1
# if other columns have records: 2
df_BONE <- core[, c(57, 58, 60, 61, 62, 63, 66, 67, 69, 72, 73, 75)]
df_BONE[is.na(df_BONE)] <- ""
BONEvec <- vector("numeric", length = nrow(core_habini))
for (i in 1:nrow(df_BONE)) {
  if (df_BONE[i, "BONE"] != "Y" & all(df_BONE[i, -1] != "Y")) BONEvec[i] <- 0
  if (df_BONE[i, "BONE"] == "Y" & all(df_BONE[i, -1] != "Y")) BONEvec[i] <- 1
  if (any(df_BONE[i, -1] == "Y")) BONEvec[i] <- 2
}
core_habini$BONE <- BONEvec
core_habini$ECOG_C <- ifelse(is.na(core_habini$ECOG_C), ECOG_C, core_habini$ECOG_C)
core_habini$ALB <- ifelse(is.na(core_habini$ALB), ALP_avg["ALB"], core_habini$ALB)
core_habini$HB <- ifelse(is.na(core_habini$HB), ALP_avg["HB"], core_habini$HB)
core_habini$PSA <- ifelse(is.na(core_habini$PSA), log(ALP_avg["PSA"]+1), log(core_habini$PSA+1))
core_habini$ALP <- ifelse(is.na(core_habini$ALP), log(ALP_avg["ALP"]+1), log(core_habini$ALP+1))
core_habini$LYMPH_NODES <- ifelse(core_habini$LYMPH_NODES == "Y", 1, 0)

# re-implement prepare_habini_validation.pl
# source data file is CoreTable_validation.csv and CoreTable_leaderboard.csv
# Choose the same columns as core_habini and manipulate following the same way
core <- read.csv("rawdata/final/CoreTable_validation.csv", header = T, stringsAsFactors = F)
core[, 26] <- as.numeric(core[, 26])
ECOG_C <- sum(core[, 26], na.rm = T)
ECOG_C_count <- sum(!is.na(core[, 26]))
ECOG_C <- ECOG_C/ECOG_C_count

core[, 30:50] <- lapply(core[, 30:50], function(x) as.numeric(x))
ALP <- colSums(core[, 30:50], na.rm = T)
ALP_count <- unlist(lapply(core[, 30:50], function(x) sum(!is.na(x))))
ALP_avg <- ALP / ALP_count

core_habini_validation1 <- core[, c(3, 83, 36, 57, 26, 47, 35, 39, 30, 59)]
core_habini_validation1$ANALGESICS <- ifelse(core_habini_validation1$ANALGESICS == "YES", 1, 0)
core_habini_validation1$LDH <- vapply(core_habini_validation1$LDH, 
                          function(x) {
                            if (is.na(x)) {return(1)} # avg for LDH is 294.187, return 1
                            ifelse(x > 250, 1, 0)
                          }, numeric(1))
# SAME TO ABOVE
df_BONE <- core[, c(57, 58, 60, 61, 62, 63, 66, 67, 69, 72, 73, 75)]
df_BONE[is.na(df_BONE)] <- ""
BONEvec <- vector("numeric", length = nrow(core_habini_validation1))
for (i in 1:nrow(df_BONE)) {
  if (df_BONE[i, "BONE"] != "Y" & all(df_BONE[i, -1] != "Y")) BONEvec[i] <- 0
  if (df_BONE[i, "BONE"] == "Y" & all(df_BONE[i, -1] != "Y")) BONEvec[i] <- 1
  if (any(df_BONE[i, -1] == "Y")) BONEvec[i] <- 2
}
core_habini_validation1$BONE <- BONEvec
core_habini_validation1$ECOG_C <- ifelse(is.na(core_habini_validation1$ECOG_C), ECOG_C, core_habini_validation1$ECOG_C)
core_habini_validation1$ALB <- ifelse(is.na(core_habini_validation1$ALB), ALP_avg["ALB"], core_habini_validation1$ALB)
core_habini_validation1$HB <- ifelse(is.na(core_habini_validation1$HB), ALP_avg["HB"], core_habini_validation1$HB)
core_habini_validation1$PSA <- ifelse(is.na(core_habini_validation1$PSA), log(ALP_avg["PSA"]+1), log(core_habini_validation1$PSA+1))
core_habini_validation1$ALP <- ifelse(is.na(core_habini_validation1$ALP), log(ALP_avg["ALP"]+1), log(core_habini_validation1$ALP+1))
core_habini_validation1$LYMPH_NODES <- ifelse(core_habini_validation1$LYMPH_NODES == "Y", 1, 0)

core <- read.csv("rawdata/prostate_cancer_challenge_data_leaderboard/CoreTable_leaderboard.csv", 
                 header = T, stringsAsFactors = F)
core[, 26] <- as.numeric(core[, 26])
ECOG_C <- sum(core[, 26], na.rm = T)
ECOG_C_count <- sum(!is.na(core[, 26]))
ECOG_C <- ECOG_C/ECOG_C_count

core[, 30:50] <- lapply(core[, 30:50], function(x) as.numeric(x))
ALP <- colSums(core[, 30:50], na.rm = T)
ALP_count <- unlist(lapply(core[, 30:50], function(x) sum(!is.na(x))))
ALP_avg <- ALP / ALP_count

core_habini_validation2 <- core[, c(3, 83, 36, 57, 26, 47, 35, 39, 30, 59)]
core_habini_validation2$ANALGESICS <- ifelse(core_habini_validation2$ANALGESICS == "YES", 1, 0)
core_habini_validation2$LDH <- vapply(core_habini_validation2$LDH, 
                                      function(x) {
                                        if (is.na(x)) {return(1)} # avg for LDH is 294.187, return 1
                                        ifelse(x > 250, 1, 0)
                                      }, numeric(1))
# SAME TO ABOVE
df_BONE <- core[, c(57, 58, 60, 61, 62, 63, 66, 67, 69, 72, 73, 75)]
df_BONE[is.na(df_BONE)] <- ""
BONEvec <- vector("numeric", length = nrow(core_habini_validation2))
for (i in 1:nrow(df_BONE)) {
  if (df_BONE[i, "BONE"] != "Y" & all(df_BONE[i, -1] != "Y")) BONEvec[i] <- 0
  if (df_BONE[i, "BONE"] == "Y" & all(df_BONE[i, -1] != "Y")) BONEvec[i] <- 1
  if (any(df_BONE[i, -1] == "Y")) BONEvec[i] <- 2
}
core_habini_validation2$BONE <- BONEvec
core_habini_validation2$ECOG_C <- ifelse(is.na(core_habini_validation2$ECOG_C), ECOG_C, core_habini_validation2$ECOG_C)
core_habini_validation2$ALB <- ifelse(is.na(core_habini_validation2$ALB), ALP_avg["ALB"], core_habini_validation2$ALB)
core_habini_validation2$HB <- ifelse(is.na(core_habini_validation2$HB), ALP_avg["HB"], core_habini_validation2$HB)
core_habini_validation2$PSA <- ifelse(is.na(core_habini_validation2$PSA), log(ALP_avg["PSA"]+1), log(core_habini_validation2$PSA+1))
core_habini_validation2$ALP <- ifelse(is.na(core_habini_validation2$ALP), log(ALP_avg["ALP"]+1), log(core_habini_validation2$ALP+1))
core_habini_validation2$LYMPH_NODES <- ifelse(core_habini_validation2$LYMPH_NODES == "Y", 1, 0)
core_habini_validation <- rbind(core_habini_validation1, core_habini_validation2)
rm(list = c("core_habini_validation1", "core_habini_validation2"))

# re-implement prepare_core_table_exclude.pl
# Source data CoreTable_training.csv
core <- read.csv("rawdata/prostate_cancer_challenge_data_training/CoreTable_training.csv",
                 header = T, stringsAsFactors = F)
# before getting the selected features, re-ran the perl file
# There're too many features. So I just ran the perl program and retrieve the features
selected.features <- readLines("code_sub2/v2/selected.features.txt")
core_exclude <- core[, c("RPT", selected.features[c(1, 6:69)])]
core_exclude$AGEGRP2 <- sapply(core_exclude$AGEGRP2, 
                               function(x) {
                                 if (grepl("18-64", x)) return(1)
                                 if (grepl("65-74", x)) return(2)
                                 if (grepl(">=75", x)) return(3)
                               })
core_exclude$BMI <- as.numeric(core_exclude$BMI)
BMIavg <- mean(core_exclude$BMI, na.rm = T)
core_exclude$BMI[is.na(core_exclude$BMI)] <- BMIavg
core_exclude[, 4:18] <- lapply(core_exclude[, 4:18], 
                               function(x) {
                                 x <- as.numeric(x)
                                 avg <- mean(x, na.rm = T)
                                 if (avg > 100) {
                                   ifelse(is.na(x), log(avg+1), log(x+1))
                                 } else {
                                   ifelse(is.na(x), avg, x)
                                 }
                               })
core_exclude[, 19:66] <- lapply(core_exclude[, 19:66], 
                                function(x) {
                                  ifelse(x == "YES" | x == "Y", 1, 0)
                                })
# Give race a dummy code
race <- table(core$RACE_C)
race["Other"] <- race["Other"] + race["Hispanic"]
race_avg <- race/(1600 - 55)
core_exclude$White <- NA
core_exclude$Asian <- NA
core_exclude$Other <- NA
core_exclude$Black <- NA
for (i in 1:1600) {
  if (core[i, "RACE_C"] == "White") 
    core_exclude[i, c("White", "Asian", "Other", "Black")] <- c(1L, 0L, 0L, 0L)
  if (core[i, "RACE_C"] == "Asian")
    core_exclude[i, c("White", "Asian", "Other", "Black")] <- c(0L, 1L, 0L, 0L)
  if (core[i, "RACE_C"] == "Other")
    core_exclude[i, c("White", "Asian", "Other", "Black")] <- c(0L, 0L, 1L, 0L)
  if (core[i, "RACE_C"] == "Black")
    core_exclude[i, c("White", "Asian", "Other", "Black")] <- c(0L, 0L, 0L, 1L)
  if (core[i, "RACE_C"] == "Missing" | core[i, "RACE_C"] == "Hispanic")
    core_exclude[i, c("White", "Asian", "Other", "Black")] <- race_avg[c("White", "Asian", "Other", "Black")]
}

# re-implement prepare_core_table_exclude_validation.pl
core <- read.csv("rawdata/final/CoreTable_validation.csv", header = T, stringsAsFactors = F)
core_exclude_validation1 <- core[, c("RPT", selected.features[c(1, 6:69)])]
core_exclude_validation1$AGEGRP2 <- sapply(core_exclude_validation1$AGEGRP2, 
                               function(x) {
                                 if (grepl("18-64", x)) return(1)
                                 if (grepl("65-74", x)) return(2)
                                 if (grepl(">=75", x)) return(3)
                               })
core_exclude_validation1$BMI <- as.numeric(core_exclude_validation1$BMI)
BMIavg <- mean(core_exclude_validation1$BMI, na.rm = T)
core_exclude_validation1$BMI[is.na(core_exclude_validation1$BMI)] <- BMIavg
core_exclude_validation1[, 4:18] <- lapply(core_exclude_validation1[, 4:18], 
                               function(x) {
                                 x <- as.numeric(x)
                                 avg <- mean(x, na.rm = T)
                                 if (avg > 100) {
                                   ifelse(is.na(x), log(avg+1), log(x+1))
                                 } else {
                                   ifelse(is.na(x), avg, x)
                                 }
                               })
core_exclude_validation1[, 19:66] <- lapply(core_exclude_validation1[, 19:66], 
                                function(x) {
                                  ifelse(x == "YES" | x == "Y", 1, 0)
                                })
# Give race a dummy code
race <- table(core$RACE_C)
race["Other"] <- race["Other"] + race["Hispanic"]
race_avg <- race/(1600 - 55)
core_exclude_validation1$White <- NA
core_exclude_validation1$Asian <- NA
core_exclude_validation1$Other <- NA
core_exclude_validation1$Black <- NA
for (i in 1:313) {
  if (core[i, "RACE_C"] == "White") 
    core_exclude_validation1[i, c("White", "Asian", "Other", "Black")] <- c(1L, 0L, 0L, 0L)
  if (core[i, "RACE_C"] == "Asian")
    core_exclude_validation1[i, c("White", "Asian", "Other", "Black")] <- c(0L, 1L, 0L, 0L)
  if (core[i, "RACE_C"] == "Other")
    core_exclude_validation1[i, c("White", "Asian", "Other", "Black")] <- c(0L, 0L, 1L, 0L)
  if (core[i, "RACE_C"] == "Black")
    core_exclude_validation1[i, c("White", "Asian", "Other", "Black")] <- c(0L, 0L, 0L, 1L)
  if (core[i, "RACE_C"] == "Missing" | core[i, "RACE_C"] == "Hispanic")
    core_exclude_validation1[i, c("White", "Asian", "Other", "Black")] <- race_avg[c("White", "Asian", "Other", "Black")]
}

core <- read.csv("rawdata/prostate_cancer_challenge_data_leaderboard/CoreTable_leaderboard.csv",
                 header = T, stringsAsFactors = F)
core_exclude_validation2 <- core[, c("RPT", selected.features[c(1, 6:69)])]
core_exclude_validation2$AGEGRP2 <- sapply(core_exclude_validation2$AGEGRP2, 
                                           function(x) {
                                             if (grepl("18-64", x)) return(1)
                                             if (grepl("65-74", x)) return(2)
                                             if (grepl(">=75", x)) return(3)
                                           })
core_exclude_validation2$BMI <- as.numeric(core_exclude_validation2$BMI)
BMIavg <- mean(core_exclude_validation2$BMI, na.rm = T)
core_exclude_validation2$BMI[is.na(core_exclude_validation2$BMI)] <- BMIavg
core_exclude_validation2[, 4:18] <- lapply(core_exclude_validation2[, 4:18], 
                                           function(x) {
                                             x <- as.numeric(x)
                                             avg <- mean(x, na.rm = T)
                                             if (avg > 100) {
                                               ifelse(is.na(x), log(avg+1), log(x+1))
                                             } else {
                                               ifelse(is.na(x), avg, x)
                                             }
                                           })
core_exclude_validation2[, 19:66] <- lapply(core_exclude_validation2[, 19:66], 
                                            function(x) {
                                              ifelse(x == "YES" | x == "Y", 1, 0)
                                            })
# Give race a dummy code
race <- table(core$RACE_C)
race["Other"] <- race["Other"] + race["Hispanic"]
race_avg <- race/(1600 - 55)
core_exclude_validation2$White <- NA
core_exclude_validation2$Asian <- NA
core_exclude_validation2$Other <- NA
core_exclude_validation2$Black <- NA
for (i in 1:157) {
  if (core[i, "RACE_C"] == "White") 
    core_exclude_validation2[i, c("White", "Asian", "Other", "Black")] <- c(1L, 0L, 0L, 0L)
  if (core[i, "RACE_C"] == "Asian")
    core_exclude_validation2[i, c("White", "Asian", "Other", "Black")] <- c(0L, 1L, 0L, 0L)
  if (core[i, "RACE_C"] == "Other")
    core_exclude_validation2[i, c("White", "Asian", "Other", "Black")] <- c(0L, 0L, 1L, 0L)
  if (core[i, "RACE_C"] == "Black")
    core_exclude_validation2[i, c("White", "Asian", "Other", "Black")] <- c(0L, 0L, 0L, 1L)
  if (core[i, "RACE_C"] == "Missing" | core[i, "RACE_C"] == "Hispanic")
    core_exclude_validation2[i, c("White", "Asian", "Other", "Black")] <- race_avg[c("White", "Asian", "Other", "Black")]
}
core_exclude_validation <- rbind(core_exclude_validation1, core_exclude_validation2)
rm(list = c("core_exclude_validation1", "core_exclude_validation2"))

# re-implement initialize_habini.pl
# and use GuanRank as DEATH
core <- read.csv("rawdata/prostate_cancer_challenge_data_training/CoreTable_training.csv",
                 header = T, stringsAsFactors = F)
habini_ASC_death <- cbind(DEATH = death[grepl("^ASC", death$RPT), "DEATH"], 
                          DISCONT = core[grepl("^ASC", core$RPT), "DISCONT"],
                          core_habini[grepl("^ASC", core_habini$RPT), -1])
habini_CEL_death <- cbind(DEATH = death[grepl("^CEL", death$RPT), "DEATH"],
                          DISCONT = core[grepl("^CEL", core$RPT), "DISCONT"],
                          core_habini[grepl("^CEL", core_habini$RPT), -1])
habini_VEN_death <- cbind(DEATH = death[grepl("^VEN", death$RPT), "DEATH"], 
                          DISCONT = core[grepl("^VEN", core$RPT), "DISCONT"],
                          core_habini[grepl("^VEN", core_habini$RPT), -1])
habini_ASC_death_GR <- cbind(DEATH = GRtable_ASC[, "rankValue"], 
                          DISCONT = core[grepl("^ASC", core$RPT), "DISCONT"],
                          core_habini[grepl("^ASC", core_habini$RPT), -1])
habini_CEL_death_GR <- cbind(DEATH = GRtable_CEL[, "rankValue"],
                          DISCONT = core[grepl("^CEL", core$RPT), "DISCONT"],
                          core_habini[grepl("^CEL", core_habini$RPT), -1])
habini_VEN_death_GR <- cbind(DEATH = GRtable_VEN[, "rankValue"], 
                          DISCONT = core[grepl("^VEN", core$RPT), "DISCONT"],
                          core_habini[grepl("^VEN", core_habini$RPT), -1])



# re-implement initialize.pl
core <- read.csv("rawdata/prostate_cancer_challenge_data_training/CoreTable_training.csv",
                 header = T, stringsAsFactors = F)
ASC_death <- cbind(DEATH = death[grepl("^ASC", death$RPT), "DEATH"], 
                   DISCONT = core[grepl("^ASC", core$RPT), "DISCONT"],
                   core_exclude[grepl("^ASC", core_exclude$RPT), -1])
CEL_death <- cbind(DEATH = death[grepl("^CEL", death$RPT), "DEATH"],
                   DISCONT = core[grepl("^CEL", core$RPT), "DISCONT"],
                   core_exclude[grepl("^CEL", core_exclude$RPT), -1])
VEN_death <- cbind(DEATH = death[grepl("^VEN", death$RPT), "DEATH"], 
                   DISCONT = core[grepl("^VEN", core$RPT), "DISCONT"],
                   core_exclude[grepl("^VEN", core_exclude$RPT), -1])
ASC_death_GR <- cbind(DEATH = GRtable_ASC[, "rankValue"], 
                   DISCONT = core[grepl("^ASC", core$RPT), "DISCONT"],
                   core_exclude[grepl("^ASC", core_exclude$RPT), -1])
CEL_death_GR <- cbind(DEATH = GRtable_CEL[, "rankValue"],
                   DISCONT = core[grepl("^CEL", core$RPT), "DISCONT"],
                   core_exclude[grepl("^CEL", core_exclude$RPT), -1])
VEN_death_GR <- cbind(DEATH = GRtable_VEN[, "rankValue"], 
                   DISCONT = core[grepl("^VEN", core$RPT), "DISCONT"],
                   core_exclude[grepl("^VEN", core_exclude$RPT), -1])

# re-implement start_eva.pl
final_test <- cbind(ZERO = 0, core_exclude_validation[, -1])

# re-implement start_eva_habini.pl
final_test_habini <- cbind(ZERO = 0, core_habini_validation[, -1])


