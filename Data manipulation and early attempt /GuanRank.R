setwd("~/Desktop/GuanLab/")
core <- read.csv("rawdata/prostate_cancer_challenge_data_training/CoreTable_training.csv",
                   header = T, stringsAsFactors = F)
GRtable <- core[, c("RPT", "LKADT_P", "DEATH")]
GRtable$LKADT_P <- as.numeric(GRtable$LKADT_P)
GRtable$DEATH <- ifelse(GRtable$DEATH == "YES", 1, 0)

GuanRankR <- function(GRtable) {
  
  Survfit <- survfit(Surv(GRtable$LKADT_P, GRtable$DEATH) ~ 1)
  Survdf <- data.frame(time = Survfit$time, prob = Survfit$surv)
  prob <- vapply(GRtable$LKADT_P, 
                 function(x) return(Survdf[Survdf$time == x, "prob"]), numeric(1))
  
  GuanRank <- vector("numeric", length = nrow(GRtable))
  
  for (i in 1:nrow(GRtable)) {
    
    tA <- GRtable[i, "LKADT_P"]
    rA <- prob[i]
    
    if (GRtable[i, "DEATH"] == 1) {
      
      tBgttA <- prob[GRtable$LKADT_P > tA]
      tBletA_sBeq0 <- prob[GRtable$LKADT_P <= tA & GRtable$DEATH == 0]
      tBeqtA_sBeq1 <- prob[GRtable$LKADT_P == tA & GRtable$DEATH == 1]
      GuanRank[i] <- ifelse(length(tBgttA) == 0, 0, 1 * length(tBgttA)) +
        ifelse(length(tBletA_sBeq0) == 0, 0, sum(rA/tBletA_sBeq0)) + 
        ifelse(length(tBeqtA_sBeq1) == 0, 0, 0.5 * length(tBeqtA_sBeq1))
      
    }
    
    if (GRtable[i, "DEATH"] == 0) {
      
      tBgetA_sBeq0 <- prob[GRtable$LKADT_P >= tA & GRtable$DEATH == 0]
      tBgetA_sBeq1 <- prob[GRtable$LKADT_P >= tA & GRtable$DEATH == 1]
      tBlttA_sBeq0 <- prob[GRtable$LKADT_P < tA & GRtable$DEATH == 0]
      GuanRank[i] <- ifelse(length(tBgetA_sBeq0) == 0, 0, sum(1 - 0.5*tBgetA_sBeq0/rA)) + 
        ifelse(length(tBgetA_sBeq1) == 0, 0, sum(1 - tBgetA_sBeq1/rA)) + 
        ifelse(length(tBlttA_sBeq0) == 0, 0, sum(0.5*rA/tBlttA_sBeq0))
      
    }
  }
  
  return(list(prob = prob, rankValue = GuanRank/max(GuanRank)))
}

GRtable_ASC <- GRtable[grepl("^ASC", GRtable$RPT), ]
GRtable_CEL <- GRtable[grepl("^CEL", GRtable$RPT), ]
GRtable_VEN <- GRtable[grepl("^VEN", GRtable$RPT), ]

GRtable_ASC <- cbind(GRtable_ASC, GuanRankR(GRtable_ASC))
GRtable_CEL <- cbind(GRtable_CEL, GuanRankR(GRtable_CEL))
GRtable_VEN <- cbind(GRtable_VEN, GuanRankR(GRtable_VEN))

# test
ASC_test <- read.table("ASCENT2_train_target.txt", sep = "\t", header = F, stringsAsFactors = F)
ASC_test <- ASC_test[order(ASC_test$V1), ]
cor(ASC_test$V2, GRtable_ASC$rankValue)

CEL_test <- read.table("CELGENE_train_target.txt", sep = "\t", header = F, stringsAsFactors = F)
CEL_test <- CEL_test[order(CEL_test$V1), ]
cor(CEL_test$V2, GRtable_CEL$rankValue)

VEN_test <- read.table("EFC6546_train_target.txt", sep = "\t", header = F, stringsAsFactors = F)
VEN_test <- VEN_test[order(VEN_test$V1), ]
cor(VEN_test$V2, GRtable_VEN$rankValue)
