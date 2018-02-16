library(ggplot2)
library(randomForest)
labels <- core[, c("RPT", "LKADT_P", "DEATH", "DISCONT", "ENDTRS_C", "ENTRT_PC")]
labels$DISCONT <- as.integer(labels$DISCONT)
labels$ENTRT_PC <- as.integer(labels$ENTRT_PC)
labels <- na.omit(labels)

labels$Early_Death <- ifelse(labels$DEATH == 'YES' & labels$LKADT_P <= 93, 1, 0)
labels$Early_Discont <- ifelse(grepl("AE", labels$ENDTRS_C) & labels$ENTRT_PC <= 93, 1, 0)
labels$new_label <- ifelse(labels$Early_Death == 0 & labels$Early_Discont == 0 & labels$DISCONT == 0, 0, 1)
# labels$new_label <- labels$Early_Death

table(labels[grepl("^VEN", labels$RPT), "new_label"])
cor(labels$DISCONT, labels$new_label)

coreData <- cbind(labels[, c("RPT", "new_label")], 
                  DEATH = death[death$RPT %in% labels$RPT, "DEATH"],
                  DEATH_day = core[core$RPT %in% labels$RPT, "LKADT_P"],
                  DISCONT = core[core$RPT %in% labels$RPT, "DISCONT"],
                  DISCONT_day = core[core$RPT %in% labels$RPT, "ENTRT_PC"],
                  core_habini[core_habini$RPT %in% labels$RPT, -1],
                  core_exclude[core_exclude$RPT %in% labels$RPT, -1])
coreData_Normal <- coreData
coreData_Normal[, c(7:84)] <- lapply(coreData_Normal[, c(7:84)], function(x) (x-min(x))/(max(x)-min(x)))

ASC <- coreData[grepl("^ASC", coreData$RPT), ]
CEL <- coreData[grepl("^CEL", coreData$RPT), ]
VEN <- coreData[grepl("^VEN", coreData$RPT), ]

ASC_Normal <- coreData_Normal[grepl("^ASC", coreData_Normal$RPT), ]
CEL_Normal <- coreData_Normal[grepl("^CEL", coreData_Normal$RPT), ]
VEN_Normal <- coreData_Normal[grepl("^VEN", coreData_Normal$RPT), ]

# ------------ Basic learner comparison (78 features) -----------------
kFoldCV <- function(data, fold = 5, seed) {
  set.seed(seed)
  data <- data[sample(nrow(data)), ]
  folds <- cut(seq(1, nrow(data)), breaks = fold, labels = F)
  
  result <- data.frame()
  for (i in 1:fold) {
    testIndexes <- which(folds == i, arr.ind = T)
    test <- data[testIndexes, ]
    train <- data[-testIndexes, ]
    
    train$DISCONT <- as.numeric(as.character(train$DISCONT))
    train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
    train <- na.omit(train)
    test$DISCONT <- as.numeric(as.character(test$DISCONT))
    test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
    test <- na.omit(test)
    
    model_bag_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "bag",
                                     seed = 1234, target_column = "DISCONT")
    model_cox_78_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 7:84)], test[, c(5, 7:84)],
                                     model = "cox", seed = 1234, target_column = "DISCONT")
    model_lm_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "lm",
                                    seed = 1234, target_column = "DISCONT")
    model_logit_78_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                       seed = 1234, target_column = "DISCONT")
    model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                    seed = 1234, target_column = "DISCONT")
    
    res <- data.frame(identity = rep(c("model_bag_78_discont", 
                                       "model_lm_78_discont",
                                       "model_cox_78_discont",
                                       "model_logit_78_discont",
                                       "model_rf_78_discont"),
                                     each = 2),
                      value = c(model_bag_78_discont$roc$auc, model_bag_78_discont$pr$auc.integral,
                                model_lm_78_discont$roc$auc, model_lm_78_discont$pr$auc.integral,
                                model_cox_78_discont$roc$auc, model_cox_78_discont$pr$auc.integral,
                                model_logit_78_discont$roc$auc, model_logit_78_discont$pr$auc.integral,
                                model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral),
                      stringsAsFactors = F)
    
    res$model <- sapply(res$identity, function(x) {
      
      if (grepl("bag", x)) return("BAG-CART")
      if (grepl("cox", x)) return("Cox")
      if (grepl("lm", x)) return("Linear Regression")
      if (grepl("logit", x)) return("Logistic Regression")
      if (grepl("rf", x)) return("Random Forest")
      
    })
    res$gold_standard <- unlist(lapply(strsplit(res$identity, "_"), "[", 4))
    res$curve <- rep(c("AUC", "AUPRC"), 5)
    
    result <- rbind(result, res)
  }
  return(result)
}
ASCinCohortCV_Normal <- data.frame()
CELinCohortCV_Normal <- data.frame()
VENinCohortCV_Normal <- data.frame()
for (i in 1:7) {
  ASCinCohortCV_Normal <- rbind(ASCinCohortCV_Normal, kFoldCV(ASC_Normal, 5, i))
  CELinCohortCV_Normal <- rbind(CELinCohortCV_Normal, kFoldCV(CEL_Normal, 5, i))
  VENinCohortCV_Normal <- rbind(VENinCohortCV_Normal, kFoldCV(VEN_Normal, 5, i))
}
ASCinCohortCV_Normal$Cohort <- "ASC"
CELinCohortCV_Normal$Cohort <- "CEL"
VENinCohortCV_Normal$Cohort <- "VEN"
TinCohortCV_Normal <- rbind(ASCinCohortCV_Normal, CELinCohortCV_Normal, VENinCohortCV_Normal)
p_AUC <- ggplot(subset(TinCohortCV_Normal, curve == "AUC"))
p_AUC <- p_AUC + geom_boxplot(aes(Cohort, value, fill = model), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1), legend.position = "bottom") + scale_fill_discrete("Models") + 
  labs(x = "Cohorts", y = "Area Under Curve", title = "Area Under ROC Curve") + scale_y_continuous(limits = c(0.2, 0.8)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + guides(fill = guide_legend(nrow=2,byrow=TRUE))
p_AUPRC <- ggplot(subset(TinCohortCV_Normal, curve == "AUPRC"))
p_AUPRC <- p_AUPRC + geom_boxplot(aes(Cohort, value, fill = model), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1), legend.position = "bottom") + scale_fill_discrete("Models") + 
  labs(x = "Cohorts", y = "Area Under Curve", title = "Area Under PR Curve") + scale_y_continuous(limits = c(0, 0.5)) + 
  guides(fill = guide_legend(nrow=2,byrow=TRUE))
# fix random forest

# -------------------- Part2 predictor choose --------------------
kFoldCV <- function(data, fold = 5, seed) {
  set.seed(seed)
  data <- data[sample(nrow(data)), ]
  folds <- cut(seq(1, nrow(data)), breaks = fold, labels = F)
  
  result <- data.frame()
  for (i in 1:fold) {
    testIndexes <- which(folds == i, arr.ind = T)
    test <- data[testIndexes, ]
    train <- data[-testIndexes, ]
    
    train$DISCONT <- as.numeric(as.character(train$DISCONT))
    train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
    train <- na.omit(train)
    test$DISCONT <- as.numeric(as.character(test$DISCONT))
    test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
    test <- na.omit(test)
    
    model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                    seed = 1234, target_column = "DISCONT")
    model_rf_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
                                    seed = 1234, target_column = "DISCONT")
    
    res <- data.frame(identity = rep(c("model_rf_78_discont", "model_rf_78_new"), each = 2),
                      value = c(model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral,
                                model_rf_78_new$roc$auc, model_rf_78_new$pr$auc.integral),
                      stringsAsFactors = F)

    res$gold_standard <- unlist(lapply(strsplit(res$identity, "_"), "[", 4))
    res$curve <- rep(c("AUC", "AUPRC"), 2)
    
    result <- rbind(result, res)
  }
  return(result)
}
ASCinCohortCV_Normal <- data.frame()
CELinCohortCV_Normal <- data.frame()
VENinCohortCV_Normal <- data.frame()
for (i in 1:7) {
  ASCinCohortCV_Normal <- rbind(ASCinCohortCV_Normal, kFoldCV(ASC_Normal, 5, i))
  CELinCohortCV_Normal <- rbind(CELinCohortCV_Normal, kFoldCV(CEL_Normal, 5, i))
  VENinCohortCV_Normal <- rbind(VENinCohortCV_Normal, kFoldCV(VEN_Normal, 5, i))
}
ASCinCohortCV_Normal$Cohort <- "ASC"
CELinCohortCV_Normal$Cohort <- "CEL"
VENinCohortCV_Normal$Cohort <- "VEN"
TinCohortCV_Normal <- rbind(ASCinCohortCV_Normal, CELinCohortCV_Normal, VENinCohortCV_Normal)

p_AUC <- ggplot(subset(TinCohortCV_Normal, curve == "AUC"))
p_AUC <- p_AUC + geom_boxplot(aes(Cohort, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Cohort", y = "Area Under Curve", title = "Area Under ROC Curve") + scale_y_continuous(limits = c(0.2, 0.8)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed")
p_AUPRC <- ggplot(subset(TinCohortCV_Normal, curve == "AUPRC"))
p_AUPRC <- p_AUPRC + geom_boxplot(aes(Cohort, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Cohort", y = "Area Under Curve", title = "Area Under PR Curve") + scale_y_continuous(limits = c(0, 0.5))










# -------------------- Part3 feature choose --------------------------
# 78 still gives the best performance
kFoldCV <- function(data, fold = 5, seed) {
  set.seed(seed)
  data <- data[sample(nrow(data)), ]
  folds <- cut(seq(1, nrow(data)), breaks = fold, labels = F)
  
  result <- data.frame()
  for (i in 1:fold) {
    testIndexes <- which(folds == i, arr.ind = T)
    test <- data[testIndexes, ]
    train <- data[-testIndexes, ]
    
    train$DISCONT <- as.numeric(as.character(train$DISCONT))
    train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
    train <- na.omit(train)
    test$DISCONT <- as.numeric(as.character(test$DISCONT))
    test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
    test <- na.omit(test)
    
    features <- randomForest(new_label ~ ., train[, c(2, 7:84)])$importance
    result <- rbind(result, features)
  }
  result$feature <- row.names(result)
  row.names(result) <- 1:nrow(result)
  return(result)
}
# ASCfeature <- data.frame()
# CELfeature <- data.frame()
# VENfeature <- data.frame()
coreFeatures <- data.frame()
for (i in 1:7) {
  # ASCfeature <- rbind(ASCfeature, kFoldCV(ASC_Normal, 5, i))
  # CELfeature <- rbind(CELfeature, kFoldCV(CEL_Normal, 5, i))
  # VENfeature <- rbind(VENfeature, kFoldCV(VEN_Normal, 5, i))
  coreFeatures <- rbind(coreFeatures, kFoldCV(coreData_Normal, 5, i))
}
ASCfeature$cohort <- "ASC"
CELfeature$cohort <- "CEL"
VENfeature$cohort <- "VEN"
Tfeature <- rbind(ASCfeature, CELfeature, VENfeature)
# Tfeature$feature <- gsub("[0-9]+", "", Tfeature$feature)
# Tfeature$IncNodePurity <- log(Tfeature$IncNodePurity + 1e-21)
coreFeatures$feature <- gsub("[0-9]+", "", coreFeatures$feature)
coreFeatures$IncNodePurity <- log(coreFeatures$IncNodePurity + 1e-20)
color <- rep("black", 78)
color[c(4, 11, 15, 23, 32, 64, 68, 69, 77)] <- "red"
p_feature <- ggplot(coreFeatures, aes(feature, IncNodePurity))
p_feature + geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Features", y = "log-importance", title = "Importance of features")

# ALP, AST, BMI, CREAT, HB, NEU, PLT, PSA, WBC
# selected_features <- c("ALP", "AST", "BMI", "CREAT", "HB", "NEU", "PLT", "PSA", "WBC", "AGEGRP2", "ALB",
#                        "ALT", "CA", "CREACL", "ECOG_C", "MG", "NA.", "PHOS", "RBC", "TBILI", "TESTO", "TPRO")
# selected_features <- c("AGEGRP2", "ALP", "ALT", "ALB", "AST", "BMI", "CA", "CREAT", "HB", "MG", "NA.", 
#                        "NEU", "PHOS", "PLT", "PSA", "RBC", "CREACAL", "ECOG_C", "TBILI", "TESTO", "TPRO", "WBC")

# Another way to rank the feature importance
library(caret)
control <- trainControl(method = "repeatedcv", number = 5, repeats = 7)
# model_ASC <- train(new_label ~ ., data = ASC_Normal[, c(2, 7:84)], method = "rf", trControl = control, importance = TRUE)
# importance_ASC <- varImp(model_ASC, scale=FALSE)
# model_CEL <- train(new_label ~ ., data = CEL_Normal[, c(2, 7:84)], method = "rf", trControl = control, importance = TRUE)
# importance_CEL <- varImp(model_CEL, scale=FALSE)
# model_VEN <- train(new_label ~ ., data = VEN_Normal[, c(2, 7:84)], method = "rf", trControl = control, importance = TRUE)
# importance_VEN <- varImp(model_VEN, scale=FALSE)
model <- train(new_label ~ ., data = coreData_Normal[, c(2, 7:84)], method = "rf", trControl = control, importance = TRUE)
importance <- varImp(model, scale = F)
plot(importance)
# plot(importance_ASC)
# plot(importance_CEL)
# plot(importance_VEN)

selected_features <- row.names(importance$importance)[which(importance$importance > 1)]

control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats = 10)
results <- rfe(CEL_Normal[, 7:84], CEL_Normal[, 2], sizes=c(1:78), rfeControl=control)




kFoldCV <- function(data, fold = 5, seed) {
  set.seed(seed)
  data <- data[sample(nrow(data)), ]
  folds <- cut(seq(1, nrow(data)), breaks = fold, labels = F)
  
  result <- data.frame()
  for (i in 1:fold) {
    testIndexes <- which(folds == i, arr.ind = T)
    test <- data[testIndexes, ]
    train <- data[-testIndexes, ]
    
    train$DISCONT <- as.numeric(as.character(train$DISCONT))
    train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
    train <- na.omit(train)
    test$DISCONT <- as.numeric(as.character(test$DISCONT))
    test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
    test <- na.omit(test)
    
    model_rf_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
                                seed = 1234, target_column = "DISCONT")
    model_rf_X_new <- modelFit(new_label ~ ., train[, c("new_label", selected_features)], 
                               test[, c("DISCONT", selected_features)], 
                               model = "rf",
                               seed = 1234, target_column = "DISCONT")
    model_rf_9_halabi <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "rf",
                                  seed = 1234, target_column = "DISCONT")
    
    res <- data.frame(identity = rep(c("model_rf_78_new", "model_rf_X_new", "model_rf_9_halabi"), each = 2),
                      value = c(model_rf_78_new$roc$auc, model_rf_78_new$pr$auc.integral,
                                model_rf_X_new$roc$auc, model_rf_X_new$pr$auc.integral,
                                model_rf_9_halabi$roc$auc, model_rf_9_halabi$pr$auc.integral),
                      stringsAsFactors = F)
    
    res$feature_num <- unlist(lapply(strsplit(res$identity, "_"), "[", 3))
    res$gold_standard <- unlist(lapply(strsplit(res$identity, "_"), "[", 4))
    res$category <- paste0(res$gold_standard, " (", res$feature_num, ")")
    res$curve <- rep(c("AUC", "AUPRC"), 3)
    
    result <- rbind(result, res)
  }
  return(result)
}
ASCinCohortCV_Normal <- data.frame()
CELinCohortCV_Normal <- data.frame()
VENinCohortCV_Normal <- data.frame()
for (i in 1:7) {
  ASCinCohortCV_Normal <- rbind(ASCinCohortCV_Normal, kFoldCV(ASC_Normal, 5, i))
  CELinCohortCV_Normal <- rbind(CELinCohortCV_Normal, kFoldCV(CEL_Normal, 5, i))
  VENinCohortCV_Normal <- rbind(VENinCohortCV_Normal, kFoldCV(VEN_Normal, 5, i))
}
ASCinCohortCV_Normal$Cohort <- "ASC"
CELinCohortCV_Normal$Cohort <- "CEL"
VENinCohortCV_Normal$Cohort <- "VEN"
TinCohortCV_Normal <- rbind(ASCinCohortCV_Normal, CELinCohortCV_Normal, VENinCohortCV_Normal)

p_AUC <- ggplot(subset(TinCohortCV_Normal, curve == "AUC"))
p_AUC <- p_AUC + geom_boxplot(aes(Cohort, value, fill = category), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Cohort", y = "Area Under Curve", title = "Area Under ROC Curve") + scale_y_continuous(limits = c(0.2, 0.8)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed")
p_AUPRC <- ggplot(subset(TinCohortCV_Normal, curve == "AUPRC"))
p_AUPRC <- p_AUPRC + geom_boxplot(aes(Cohort, value, fill = category), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Cohort", y = "Area Under Curve", title = "Area Under PR Curve") + scale_y_continuous(limits = c(0, 0.5))


# -------------------- Part4 Ensemble study -----------------------
kFoldCV <- function(data, fold = 5, seed) {
  set.seed(seed)
  data <- data[sample(nrow(data)), ]
  folds <- cut(seq(1, nrow(data)), breaks = fold, labels = F)
  
  result <- data.frame()
  for (i in 1:fold) {
    testIndexes <- which(folds == i, arr.ind = T)
    test <- data[testIndexes, ]
    train <- data[-testIndexes, ]
    
    train$DISCONT <- as.numeric(as.character(train$DISCONT))
    train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
    train <- na.omit(train)
    test$DISCONT <- as.numeric(as.character(test$DISCONT))
    test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
    test <- na.omit(test)
    
    model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                    seed = 1234, target_column = "DISCONT")
    model_rf_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
                                seed = 1234, target_column = "DISCONT")
    model_rf_78_ensemble <- modelFitEnsemble(train[, c(2, 5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                             seed = 1234, target_column = "DISCONT")
    
    res <- data.frame(identity = rep(c("model_rf_78_discont", "model_rf_78_new", "model_rf_78_ensemble"), each = 2),
                      value = c(model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral,
                                model_rf_78_new$roc$auc, model_rf_78_new$pr$auc.integral,
                                model_rf_78_ensemble$roc$auc, model_rf_78_ensemble$pr$auc.integral),
                      stringsAsFactors = F)
    
    res$gold_standard <- unlist(lapply(strsplit(res$identity, "_"), "[", 4))
    res$curve <- rep(c("AUC", "AUPRC"), 3)
    
    result <- rbind(result, res)
  }
  return(result)
}
modelFitEnsemble <- function(train.data, test.data, model = c("bag", "cox", "lm", "logit", "rf"), 
                             seed, target_column) {
  
  require(ipred)
  require(survival)
  require(randomForest)
  require(PRROC)
  
  set.seed(seed)
  
  fit_model <- switch(
    model,
    "bag" = {
      death <- bagging(formula = new_label ~ ., data = train.data[, -which(colnames(train.data) == "DISCONT")])
      discont <- bagging(formula = DISCONT ~ ., data = train.data[, -which(colnames(train.data) == "new_label")])
      list(death = death, discont = discont)
    },
    "lm" = {
      death <- lm(formula = new_label ~ ., data = train.data[, -which(colnames(train.data) == "DISCONT")])
      discont <- lm(formula = DISCONT ~ ., data = train.data[, -which(colnames(train.data) == "new_label")])
      list(death = death, discont = discont)
    },
    "logit" = {
      death <- glm(formula = as.factor(new_label) ~ ., data = train.data[, -which(colnames(train.data) == "DISCONT")], family = binomial(link = "logit"))
      discont <- glm(formula = as.factor(DISCONT) ~ ., data = train.data[, -which(colnames(train.data) == "new_label")], family = binomial(link = "logit"))
      list(death = death, discont = discont)
    },
    "rf" = {
      death <- randomForest(formula = new_label ~ ., data = train.data[, -which(colnames(train.data) == "DISCONT")])
      discont <- randomForest(formula = DISCONT ~ ., data = train.data[, -which(colnames(train.data) == "new_label")])
      list(death = death, discont = discont)
    }
    # "cox" = {
    #   death <- coxph(formula = Surv(DEATH_day, DEATH) ~ ., data = train.data[, -which(colnames(train.data) %in% c("DISCONT", "DISCONT_day"))], ties = "breslow")
    #   discont <- coxph(formula = Surv(DISCONT_day, DISCONT) ~ ., data = train.data[, -which(colnames(train.data) %in% c("DEATH", "DEATH_day"))], ties = "breslow")
    #   list(death = death, discont = discont)
    # }
  )
  
  pred_Death <- predict(fit_model$death, test.data[, -which(colnames(test.data) %in% target_column)])
  pred_Discont <- predict(fit_model$discont, test.data[, -which(colnames(test.data) %in% target_column)])
  pred_ensemble <- (pred_Death + pred_Discont) / 2
  roc <- roc.curve(scores.class0 = pred_ensemble, weights.class0 = test.data[, target_column])
  pr <- pr.curve(scores.class0 = pred_ensemble, weights.class0 = test.data[, target_column], rand.compute = T)
  
  return(list(roc = roc, pr = pr))
}
ASCinCohortCV_Normal <- data.frame()
CELinCohortCV_Normal <- data.frame()
VENinCohortCV_Normal <- data.frame()
for (i in 1:15) {
  ASCinCohortCV_Normal <- rbind(ASCinCohortCV_Normal, kFoldCV(ASC_Normal, 5, i))
  CELinCohortCV_Normal <- rbind(CELinCohortCV_Normal, kFoldCV(CEL_Normal, 5, i))
  VENinCohortCV_Normal <- rbind(VENinCohortCV_Normal, kFoldCV(VEN_Normal, 5, i))
}
ASCinCohortCV_Normal$Cohort <- "ASC"
CELinCohortCV_Normal$Cohort <- "CEL"
VENinCohortCV_Normal$Cohort <- "VEN"
TinCohortCV_Normal <- rbind(ASCinCohortCV_Normal, CELinCohortCV_Normal, VENinCohortCV_Normal)

p_AUC <- ggplot(subset(TinCohortCV_Normal, curve == "AUC"))
p_AUC <- p_AUC + geom_boxplot(aes(Cohort, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Cohort", y = "Area Under Curve", title = "Area Under ROC Curve") + scale_y_continuous(limits = c(0.2, 0.8)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed")
p_AUPRC <- ggplot(subset(TinCohortCV_Normal, curve == "AUPRC"))
p_AUPRC <- p_AUPRC + geom_boxplot(aes(Cohort, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Cohort", y = "Area Under Curve", title = "Area Under PR Curve") + scale_y_continuous(limits = c(0, 0.5))







## ASC + CEL --> VEN --------------
train <- rbind(ASC_Normal, CEL_Normal)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
train <- na.omit(train)
test <- VEN_Normal
test$DISCONT <- as.numeric(as.character(test$DISCONT))
test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
test <- na.omit(test)

# model_bag_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "bag",
#                               seed = 1234, target_column = "DISCONT")
model_bag_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "bag",
                                seed = 1234, target_column = "DISCONT")
model_bag_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "bag",
                            seed = 1234, target_column = "DISCONT")
# model_bag_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "bag",
#                               seed = 1234, target_column = "DISCONT")
model_bag_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "bag",
                                seed = 1234, target_column = "DISCONT")
model_bag_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "bag",
                            seed = 1234, target_column = "DISCONT")
# model_lm_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "lm",
#                              seed = 1234, target_column = "DISCONT")
model_lm_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "lm",
                               seed = 1234, target_column = "DISCONT")
model_lm_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "lm",
                           seed = 1234, target_column = "DISCONT")
# model_lm_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "lm",
#                              seed = 1234, target_column = "DISCONT")
model_lm_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "lm",
                               seed = 1234, target_column = "DISCONT")
model_lm_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "lm",
                           seed = 1234, target_column = "DISCONT")
# model_logit_4_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "logit",
#                                 seed = 1234, target_column = "DISCONT")
model_logit_4_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "logit",
                                  seed = 1234, target_column = "DISCONT")
model_logit_4_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "logit",
                              seed = 1234, target_column = "DISCONT")
# model_logit_9_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "logit",
#                                 seed = 1234, target_column = "DISCONT")
model_logit_9_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "logit",
                                  seed = 1234, target_column = "DISCONT")
model_logit_9_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "logit",
                              seed = 1234, target_column = "DISCONT")
# model_rf_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "rf",
#                               seed = 1234, target_column = "DISCONT")
model_rf_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "rf",
                               seed = 1234, target_column = "DISCONT")
model_rf_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "rf",
                           seed = 1234, target_column = "DISCONT")
# model_rf_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "rf",
#                               seed = 1234, target_column = "DISCONT")
model_rf_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "rf",
                               seed = 1234, target_column = "DISCONT")
model_rf_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "rf",
                           seed = 1234, target_column = "DISCONT")
# model_bag_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "bag",
#                                seed = 1234, target_column = "DISCONT")
model_bag_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "bag",
                                 seed = 1234, target_column = "DISCONT")
model_bag_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "bag",
                             seed = 1234, target_column = "DISCONT")
# model_lm_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "lm",
#                               seed = 1234, target_column = "DISCONT")
model_lm_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "lm",
                                seed = 1234, target_column = "DISCONT")
model_lm_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "lm",
                            seed = 1234, target_column = "DISCONT")
# model_logit_78_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "logit",
#                                  seed = 1234, target_column = "DISCONT")
model_logit_78_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                   seed = 1234, target_column = "DISCONT")
model_logit_78_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "logit",
                               seed = 1234, target_column = "DISCONT")
# model_rf_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "rf",
#                               seed = 1234, target_column = "DISCONT")
model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                seed = 1234, target_column = "DISCONT")
model_rf_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
                            seed = 1234, target_column = "DISCONT")

res <- data.frame(identity = rep(c("model_lm_4_discont", "model_lm_4_new",
                                   "model_lm_9_discont", "model_lm_9_new",
                                   "model_lm_78_discont", "model_lm_78_new",

                                   "model_bag_4_discont", "model_bag_4_new",
                                   "model_bag_9_discont", "model_bag_9_new",
                                   "model_bag_78_discont", "model_bag_78_new",

                                   "model_rf_4_discont", "model_rf_4_new",
                                   "model_rf_9_discont", "model_rf_9_new",
                                   "model_rf_78_discont", "model_rf_78_new",

                                   "model_logit_4_discont", "model_logit_4_new",
                                   "model_logit_9_discont", "model_logit_9_new",
                                   "model_logit_78_discont", "model_logit_78_new"),
                                 each = 2),
                  value = c(
                    model_lm_4_discont$roc$auc, model_lm_4_discont$pr$auc.integral,
                    model_lm_4_new$roc$auc, model_lm_4_new$pr$auc.integral,

                    model_lm_9_discont$roc$auc, model_lm_9_discont$pr$auc.integral,
                    model_lm_9_new$roc$auc, model_lm_9_new$pr$auc.integral,

                    model_lm_78_discont$roc$auc, model_lm_78_discont$pr$auc.integral,
                    model_lm_78_new$roc$auc, model_lm_78_new$pr$auc.integral,


                    model_bag_4_discont$roc$auc, model_bag_4_discont$pr$auc.integral,
                    model_bag_4_new$roc$auc, model_bag_4_new$pr$auc.integral,

                    model_bag_9_discont$roc$auc, model_bag_9_discont$pr$auc.integral,
                    model_bag_9_new$roc$auc, model_bag_9_new$pr$auc.integral,

                    model_bag_78_discont$roc$auc, model_bag_78_discont$pr$auc.integral,
                    model_bag_78_new$roc$auc, model_bag_78_new$pr$auc.integral,


                    model_rf_4_discont$roc$auc, model_rf_4_discont$pr$auc.integral,
                    model_rf_4_new$roc$auc, model_rf_4_new$pr$auc.integral,

                    model_rf_9_discont$roc$auc, model_rf_9_discont$pr$auc.integral,
                    model_rf_9_new$roc$auc, model_rf_9_new$pr$auc.integral,

                    model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral,
                    model_rf_78_new$roc$auc, model_rf_78_new$pr$auc.integral,


                    model_logit_4_discont$roc$auc, model_logit_4_discont$pr$auc.integral,
                    model_logit_4_new$roc$auc, model_logit_4_new$pr$auc.integral,

                    model_logit_9_discont$roc$auc, model_logit_9_discont$pr$auc.integral,
                    model_logit_9_new$roc$auc, model_logit_9_new$pr$auc.integral,

                    model_logit_78_discont$roc$auc, model_logit_78_discont$pr$auc.integral,
                    model_logit_78_new$roc$auc, model_logit_78_new$pr$auc.integral),
                  stringsAsFactors = F)

res$model <- sapply(res$identity, function(x) {

  if (grepl("bag", x)) return("BAG-CART")
  if (grepl("rf", x)) return("RF")
  if (grepl("lm", x)) return("Linear")
  if (grepl("logit", x)) return("Logistic")

})
res$feature_num <- unlist(lapply(strsplit(res$identity, "_"), "[", 3))
res$gold_standard <- unlist(lapply(strsplit(res$identity, "_"), "[", 4))
res$model_feature <- paste0(res$model, " (", res$feature_num, ")")
res$curve <- rep(c("AUC", "AUPRC"), 24)
p_AUC <- ggplot(subset(res, curve == "AUC"))
p_AUC <- p_AUC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("ASC + CEL --> VEN") +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_color_discrete("Predictor") +
  scale_y_continuous(limits = c(0.4, 0.7))
p_AUPRC <- ggplot(subset(res, curve == "AUPRC"))
p_AUPRC <- p_AUPRC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("ASC + CEL --> VEN") +
  scale_color_discrete("Predictor") + scale_y_continuous(limits = c(0, 0.3)) +
  geom_hline(yintercept = model_lm_4_discont$pr$rand$auc.integral, linetype = "dashed")

# ## ASC + VEN --> CEL ---------------------
# train <- rbind(ASC_Normal, VEN_Normal)
# train$DISCONT <- as.numeric(as.character(train$DISCONT))
# train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
# train <- na.omit(train)
# test <- CEL_Normal
# test$DISCONT <- as.numeric(as.character(test$DISCONT))
# test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
# test <- na.omit(test)
# 
# # model_bag_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "bag",
# #                               seed = 1234, target_column = "DISCONT")
# model_bag_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "bag",
#                                 seed = 1234, target_column = "DISCONT")
# model_bag_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "bag",
#                             seed = 1234, target_column = "DISCONT")
# # model_bag_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "bag",
# #                               seed = 1234, target_column = "DISCONT")
# model_bag_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "bag",
#                                 seed = 1234, target_column = "DISCONT")
# model_bag_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "bag",
#                             seed = 1234, target_column = "DISCONT")
# # model_lm_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "lm",
# #                              seed = 1234, target_column = "DISCONT")
# model_lm_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "lm",
#                                seed = 1234, target_column = "DISCONT")
# model_lm_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "lm",
#                            seed = 1234, target_column = "DISCONT")
# # model_lm_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "lm",
# #                              seed = 1234, target_column = "DISCONT")
# model_lm_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "lm",
#                                seed = 1234, target_column = "DISCONT")
# model_lm_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "lm",
#                            seed = 1234, target_column = "DISCONT")
# # model_logit_4_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "logit",
# #                                 seed = 1234, target_column = "DISCONT")
# model_logit_4_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "logit",
#                                   seed = 1234, target_column = "DISCONT")
# model_logit_4_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "logit",
#                               seed = 1234, target_column = "DISCONT")
# # model_logit_9_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "logit",
# #                                 seed = 1234, target_column = "DISCONT")
# model_logit_9_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "logit",
#                                   seed = 1234, target_column = "DISCONT")
# model_logit_9_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "logit",
#                               seed = 1234, target_column = "DISCONT")
# # model_rf_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "rf",
# #                               seed = 1234, target_column = "DISCONT")
# model_rf_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "rf",
#                                seed = 1234, target_column = "DISCONT")
# model_rf_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "rf",
#                            seed = 1234, target_column = "DISCONT")
# # model_rf_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "rf",
# #                               seed = 1234, target_column = "DISCONT")
# model_rf_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "rf",
#                                seed = 1234, target_column = "DISCONT")
# model_rf_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "rf",
#                            seed = 1234, target_column = "DISCONT")
# # model_bag_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "bag",
# #                                seed = 1234, target_column = "DISCONT")
# model_bag_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "bag",
#                                  seed = 1234, target_column = "DISCONT")
# model_bag_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "bag",
#                              seed = 1234, target_column = "DISCONT")
# # model_lm_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "lm",
# #                               seed = 1234, target_column = "DISCONT")
# model_lm_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "lm",
#                                 seed = 1234, target_column = "DISCONT")
# model_lm_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "lm",
#                             seed = 1234, target_column = "DISCONT")
# # model_logit_78_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "logit",
# #                                  seed = 1234, target_column = "DISCONT")
# model_logit_78_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "logit",
#                                    seed = 1234, target_column = "DISCONT")
# model_logit_78_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "logit",
#                                seed = 1234, target_column = "DISCONT")
# # model_rf_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "rf",
# #                               seed = 1234, target_column = "DISCONT")
# model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
#                                 seed = 1234, target_column = "DISCONT")
# model_rf_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
#                             seed = 1234, target_column = "DISCONT")
# 
# res1 <- data.frame(identity = rep(c("model_lm_4_discont", "model_lm_4_new",
#                                    "model_lm_9_discont", "model_lm_9_new",
#                                    "model_lm_78_discont", "model_lm_78_new",
#                                    
#                                    "model_bag_4_discont", "model_bag_4_new",
#                                    "model_bag_9_discont", "model_bag_9_new",
#                                    "model_bag_78_discont", "model_bag_78_new",
#                                    
#                                    "model_rf_4_discont", "model_rf_4_new",
#                                    "model_rf_9_discont", "model_rf_9_new",
#                                    "model_rf_78_discont", "model_rf_78_new",
#                                    
#                                    "model_logit_4_discont", "model_logit_4_new",
#                                    "model_logit_9_discont", "model_logit_9_new",
#                                    "model_logit_78_discont", "model_logit_78_new"),
#                                  each = 2),
#                   value = c(
#                     model_lm_4_discont$roc$auc, model_lm_4_discont$pr$auc.integral,
#                     model_lm_4_new$roc$auc, model_lm_4_new$pr$auc.integral,
#                     
#                     model_lm_9_discont$roc$auc, model_lm_9_discont$pr$auc.integral,
#                     model_lm_9_new$roc$auc, model_lm_9_new$pr$auc.integral,
#                     
#                     model_lm_78_discont$roc$auc, model_lm_78_discont$pr$auc.integral,
#                     model_lm_78_new$roc$auc, model_lm_78_new$pr$auc.integral,
#                     
#                     
#                     model_bag_4_discont$roc$auc, model_bag_4_discont$pr$auc.integral,
#                     model_bag_4_new$roc$auc, model_bag_4_new$pr$auc.integral,
#                     
#                     model_bag_9_discont$roc$auc, model_bag_9_discont$pr$auc.integral,
#                     model_bag_9_new$roc$auc, model_bag_9_new$pr$auc.integral,
#                     
#                     model_bag_78_discont$roc$auc, model_bag_78_discont$pr$auc.integral,
#                     model_bag_78_new$roc$auc, model_bag_78_new$pr$auc.integral,
#                     
#                     
#                     model_rf_4_discont$roc$auc, model_rf_4_discont$pr$auc.integral,
#                     model_rf_4_new$roc$auc, model_rf_4_new$pr$auc.integral,
#                     
#                     model_rf_9_discont$roc$auc, model_rf_9_discont$pr$auc.integral,
#                     model_rf_9_new$roc$auc, model_rf_9_new$pr$auc.integral,
#                     
#                     model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral,
#                     model_rf_78_new$roc$auc, model_rf_78_new$pr$auc.integral,
#                     
#                     
#                     model_logit_4_discont$roc$auc, model_logit_4_discont$pr$auc.integral,
#                     model_logit_4_new$roc$auc, model_logit_4_new$pr$auc.integral,
#                     
#                     model_logit_9_discont$roc$auc, model_logit_9_discont$pr$auc.integral,
#                     model_logit_9_new$roc$auc, model_logit_9_new$pr$auc.integral,
#                     
#                     model_logit_78_discont$roc$auc, model_logit_78_discont$pr$auc.integral,
#                     model_logit_78_new$roc$auc, model_logit_78_new$pr$auc.integral),
#                   stringsAsFactors = F)
# 
# res1$model <- sapply(res1$identity, function(x) {
#   
#   if (grepl("bag", x)) return("BAG-CART")
#   if (grepl("rf", x)) return("RF")
#   if (grepl("lm", x)) return("Linear")
#   if (grepl("logit", x)) return("Logistic")
#   
# })
# res1$feature_num <- unlist(lapply(strsplit(res1$identity, "_"), "[", 3))
# res1$gold_standard <- unlist(lapply(strsplit(res1$identity, "_"), "[", 4))
# res1$model_feature <- paste0(res1$model, " (", res1$feature_num, ")")
# res1$curve <- rep(c("AUC", "AUPRC"), 24)
# 
# p1_AUC <- ggplot(subset(res1, curve == "AUC"))
# p1_AUC <- p1_AUC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("ASC + VEN --> CEL") + 
#   geom_hline(yintercept = 0.5, linetype = "dashed") + scale_color_discrete("Predictor") +
#   scale_y_continuous(limits = c(0.4, 0.7))
# p1_AUPRC <- ggplot(subset(res1, curve == "AUPRC"))
# p1_AUPRC <- p1_AUPRC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("ASC + VEN --> CEL") + 
#   scale_color_discrete("Predictor") + scale_y_continuous(limits = c(0, 0.3)) + 
#   geom_hline(yintercept = model_bag_4_discont$pr$rand$auc.integral, linetype = "dashed")
# 
# 
# ## CEL + VEN --> ASC -----------------------
# train <- rbind(VEN_Normal, CEL_Normal)
# train$DISCONT <- as.numeric(as.character(train$DISCONT))
# train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
# train <- na.omit(train)
# test <- ASC_Normal
# test$DISCONT <- as.numeric(as.character(test$DISCONT))
# test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
# test <- na.omit(test)
# 
# # model_bag_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "bag",
# #                               seed = 1234, target_column = "DISCONT")
# model_bag_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "bag",
#                                 seed = 1234, target_column = "DISCONT")
# model_bag_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "bag",
#                             seed = 1234, target_column = "DISCONT")
# # model_bag_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "bag",
# #                               seed = 1234, target_column = "DISCONT")
# model_bag_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "bag",
#                                 seed = 1234, target_column = "DISCONT")
# model_bag_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "bag",
#                             seed = 1234, target_column = "DISCONT")
# # model_lm_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "lm",
# #                              seed = 1234, target_column = "DISCONT")
# model_lm_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "lm",
#                                seed = 1234, target_column = "DISCONT")
# model_lm_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "lm",
#                            seed = 1234, target_column = "DISCONT")
# # model_lm_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "lm",
# #                              seed = 1234, target_column = "DISCONT")
# model_lm_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "lm",
#                                seed = 1234, target_column = "DISCONT")
# model_lm_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "lm",
#                            seed = 1234, target_column = "DISCONT")
# # model_logit_4_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "logit",
# #                                 seed = 1234, target_column = "DISCONT")
# model_logit_4_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "logit",
#                                   seed = 1234, target_column = "DISCONT")
# model_logit_4_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "logit",
#                               seed = 1234, target_column = "DISCONT")
# # model_logit_9_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "logit",
# #                                 seed = 1234, target_column = "DISCONT")
# model_logit_9_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "logit",
#                                   seed = 1234, target_column = "DISCONT")
# model_logit_9_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "logit",
#                               seed = 1234, target_column = "DISCONT")
# # model_rf_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "rf",
# #                               seed = 1234, target_column = "DISCONT")
# model_rf_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "rf",
#                                seed = 1234, target_column = "DISCONT")
# model_rf_4_new <- modelFit(new_label ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "rf",
#                            seed = 1234, target_column = "DISCONT")
# # model_rf_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "rf",
# #                               seed = 1234, target_column = "DISCONT")
# model_rf_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "rf",
#                                seed = 1234, target_column = "DISCONT")
# model_rf_9_new <- modelFit(new_label ~ ., train[, c(2, 7:15)], test[, c(5, 7:15)], model = "rf",
#                            seed = 1234, target_column = "DISCONT")
# # model_bag_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "bag",
# #                                seed = 1234, target_column = "DISCONT")
# model_bag_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "bag",
#                                  seed = 1234, target_column = "DISCONT")
# model_bag_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "bag",
#                              seed = 1234, target_column = "DISCONT")
# # model_lm_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "lm",
# #                               seed = 1234, target_column = "DISCONT")
# model_lm_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "lm",
#                                 seed = 1234, target_column = "DISCONT")
# model_lm_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "lm",
#                             seed = 1234, target_column = "DISCONT")
# # model_logit_78_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "logit",
# #                                  seed = 1234, target_column = "DISCONT")
# model_logit_78_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "logit",
#                                    seed = 1234, target_column = "DISCONT")
# model_logit_78_new <- modelFit(as.factor(new_label) ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "logit",
#                                seed = 1234, target_column = "DISCONT")
# # model_rf_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "rf",
# #                               seed = 1234, target_column = "DISCONT")
# model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
#                                 seed = 1234, target_column = "DISCONT")
# model_rf_78_new <- modelFit(new_label ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
#                             seed = 1234, target_column = "DISCONT")
# 
# res2 <- data.frame(identity = rep(c("model_lm_4_discont", "model_lm_4_new",
#                                    "model_lm_9_discont", "model_lm_9_new",
#                                    "model_lm_78_discont", "model_lm_78_new",
#                                    
#                                    "model_bag_4_discont", "model_bag_4_new",
#                                    "model_bag_9_discont", "model_bag_9_new",
#                                    "model_bag_78_discont", "model_bag_78_new",
#                                    
#                                    "model_rf_4_discont", "model_rf_4_new",
#                                    "model_rf_9_discont", "model_rf_9_new",
#                                    "model_rf_78_discont", "model_rf_78_new",
#                                    
#                                    "model_logit_4_discont", "model_logit_4_new",
#                                    "model_logit_9_discont", "model_logit_9_new",
#                                    "model_logit_78_discont", "model_logit_78_new"),
#                                  each = 2),
#                   value = c(
#                     model_lm_4_discont$roc$auc, model_lm_4_discont$pr$auc.integral,
#                     model_lm_4_new$roc$auc, model_lm_4_new$pr$auc.integral,
#                     
#                     model_lm_9_discont$roc$auc, model_lm_9_discont$pr$auc.integral,
#                     model_lm_9_new$roc$auc, model_lm_9_new$pr$auc.integral,
#                     
#                     model_lm_78_discont$roc$auc, model_lm_78_discont$pr$auc.integral,
#                     model_lm_78_new$roc$auc, model_lm_78_new$pr$auc.integral,
#                     
#                     
#                     model_bag_4_discont$roc$auc, model_bag_4_discont$pr$auc.integral,
#                     model_bag_4_new$roc$auc, model_bag_4_new$pr$auc.integral,
#                     
#                     model_bag_9_discont$roc$auc, model_bag_9_discont$pr$auc.integral,
#                     model_bag_9_new$roc$auc, model_bag_9_new$pr$auc.integral,
#                     
#                     model_bag_78_discont$roc$auc, model_bag_78_discont$pr$auc.integral,
#                     model_bag_78_new$roc$auc, model_bag_78_new$pr$auc.integral,
#                     
#                     
#                     model_rf_4_discont$roc$auc, model_rf_4_discont$pr$auc.integral,
#                     model_rf_4_new$roc$auc, model_rf_4_new$pr$auc.integral,
#                     
#                     model_rf_9_discont$roc$auc, model_rf_9_discont$pr$auc.integral,
#                     model_rf_9_new$roc$auc, model_rf_9_new$pr$auc.integral,
#                     
#                     model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral,
#                     model_rf_78_new$roc$auc, model_rf_78_new$pr$auc.integral,
#                     
#                     
#                     model_logit_4_discont$roc$auc, model_logit_4_discont$pr$auc.integral,
#                     model_logit_4_new$roc$auc, model_logit_4_new$pr$auc.integral,
#                     
#                     model_logit_9_discont$roc$auc, model_logit_9_discont$pr$auc.integral,
#                     model_logit_9_new$roc$auc, model_logit_9_new$pr$auc.integral,
#                     
#                     model_logit_78_discont$roc$auc, model_logit_78_discont$pr$auc.integral,
#                     model_logit_78_new$roc$auc, model_logit_78_new$pr$auc.integral),
#                   stringsAsFactors = F)
# 
# res2$model <- sapply(res2$identity, function(x) {
#   
#   if (grepl("bag", x)) return("BAG-CART")
#   if (grepl("rf", x)) return("RF")
#   if (grepl("lm", x)) return("Linear")
#   if (grepl("logit", x)) return("Logistic")
#   
# })
# res2$feature_num <- unlist(lapply(strsplit(res2$identity, "_"), "[", 3))
# res2$gold_standard <- unlist(lapply(strsplit(res2$identity, "_"), "[", 4))
# res2$model_feature <- paste0(res2$model, " (", res2$feature_num, ")")
# res2$curve <- rep(c("AUC", "AUPRC"), 24)
# p2_AUC <- ggplot(subset(res2, curve == "AUC"))
# p2_AUC <- p2_AUC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("VEN + CEL --> ASC") + 
#   geom_hline(yintercept = 0.5, linetype = "dashed") + scale_color_discrete("Predictor") +
#   scale_y_continuous(limits = c(0.4, 0.7))
# p2_AUPRC <- ggplot(subset(res2, curve == "AUPRC"))
# p2_AUPRC <- p2_AUPRC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("VEN + CEL --> ASC") + 
#   scale_color_discrete("Predictor") + scale_y_continuous(limits = c(0, 0.3)) + 
#   geom_hline(yintercept = model_bag_4_death$pr$rand$auc.integral, linetype = "dashed")
# 
# 
# 
# 
# # ----------------
# grid_arrange_share_legend(p_AUC, p1_AUC, p2_AUC, nrow = 1, ncol = 3, position = "right")
# grid_arrange_share_legend(p_AUPRC, p1_AUPRC, p2_AUPRC, nrow = 1, ncol = 3, position = "right")

library(dendextend)
hc.matrix <- as.matrix(coreData_Normal[, 7:84])
hc.labels <- unlist(lapply(strsplit(coreData_Normal$RPT, "-"), "[", 1))
row.names(hc.matrix) <- hc.labels
colorCodes <- c(ASC = "red", CELG = "green", VEN = "blue")

dend <- as.dendrogram(hclust(dist(hc.matrix)))
labels_colors(dend) <- colorCodes[hc.labels][order.dendrogram(dend)]

par(cex = 0.2)
plot(dend[[1]])

# right!
library(Rtsne)
RPT <- unlist(lapply(strsplit(coreData_Normal$RPT, split = "-"), "[", 1))
RPT <- as.factor(RPT)
color <- rainbow(length(unique(RPT)))
names(color) <- unique(RPT)
tsne <- Rtsne(coreData_Normal[, c(7:84)])
plot(tsne$Y, t = 'n')
text(tsne$Y, labels=RPT, col=color[RPT])
