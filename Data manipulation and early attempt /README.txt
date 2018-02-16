# Data manipulation
Reimplement_of_DREAM_code.R
Input(from file):
CoreTable_training.csv, CoreTable_validation.csv, CoreTable_leaderboard.csv,
[GRtable_ASC, GRtable_VEN, GRtable_CEL](come from GuanRank.R)
Output(df in program):
death, core_habini, core_habini_validation, core_exclude, core_exclude_validation,
habini_ASC_death, habini_CEL_death, habini_VEN_death, ASC_death, CEL_death,
VEN_death, final_test, final_test_habini, habini_ASC_death_GR,
habini_CEL_death_GR, habini_VEN_death_GR, ASC_death_GR, CEL_death_GR, habini_VEN_death_GR

# GuanRank Implementation
GuanRank.R
Input(from file):
CoreTable_training.csv
Output(df in program):
GRtable_ASC, GRtable_VEN, GRtable_CEL

# Cross-Validation across cohorts
Halabi_ASC+CEL_CV_AUC&AUPRC.R
Input(df in program):
habini_VEN_death, habini_ASC_death, habini_CEL_death
output(figure):
halabi_ASC+CEL_ALB+HB+PSA+ALP_lm_bagging_rf.tiff,
AUPRC_halabi_ASC+CEL_ALB+HB+PSA+ALP_lm_bagging_rf.tiff

Halabi_ASC+VEN_CV_AUC&AUPRC.R
Input(df in program):
habini_VEN_death, habini_ASC_death, habini_CEL_death
output(figure):
halabi_ASC+VEN_ALB+HB+PSA+ALP_lm_bagging_rf.tiff,
AUPRC_halabi_ASC+VEN_ALB+HB+PSA+ALP_lm_bagging_rf.tiff

Halabi_VEN+CEL_CV_AUC&AUPRC.R
Input(df in program):
habini_VEN_death, habini_ASC_death, habini_CEL_death
output(figure):
halabi_CEL+VEN_ALB+HB+PSA+ALP_lm_bagging_rf.tiff,
AUPRC_halabi_CEL+VEN_ALB+HB+PSA+ALP_lm_bagging_rf.tiff

# Cross-Validation in cohorts
InCohort_AUC&AUPRC.R
Input(df in program):
habini_VEN_death, habini_ASC_death, habini_CEL_death,
habini_CEL_death_GR, habini_VEN_death_GR, ASC_death_GR
Output(figure):
AUC_InCohort_GR.tiff, AUC_InCohort.tiff, AUPRC_InCohort.tiff, AUPRC_InCohort_GR.tiff
