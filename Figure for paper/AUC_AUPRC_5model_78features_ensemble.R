library(gridExtra)
library(grid)
# data 
# DEATH, DEATH-DAY, DISCONT, DISCONT-DAY, New Discont, ALL possible features
ASC <- cbind(GRtable_ASC[, c("RPT", "rankValue")], 
             DEATH = death[death$RPT %in% GRtable_ASC$RPT, "DEATH"],
             DEATH_day = core[core$RPT %in% GRtable_ASC$RPT, "LKADT_P"],
             DISCONT = core[core$RPT %in% GRtable_ASC$RPT, "DISCONT"],
             DISCONT_day = core[core$RPT %in% GRtable_ASC$RPT, "ENTRT_PC"],
             core_habini[core_habini$RPT %in% GRtable_ASC$RPT, -1],
             core_exclude[core_exclude$RPT %in% GRtable_ASC$RPT, -1])
CEL <- cbind(GRtable_CEL[, c("RPT", "rankValue")], 
             DEATH = death[death$RPT %in% GRtable_CEL$RPT, "DEATH"],
             DEATH_day = core[core$RPT %in% GRtable_CEL$RPT, "LKADT_P"],
             DISCONT = core[core$RPT %in% GRtable_CEL$RPT, "DISCONT"],
             DISCONT_day = core[core$RPT %in% GRtable_CEL$RPT, "ENTRT_PC"],
             core_habini[core_habini$RPT %in% GRtable_CEL$RPT, -1],
             core_exclude[core_exclude$RPT %in% GRtable_CEL$RPT, -1])
VEN <- cbind(GRtable_VEN[, c("RPT", "rankValue")], 
             DEATH = death[death$RPT %in% GRtable_VEN$RPT, "DEATH"],
             DEATH_day = core[core$RPT %in% GRtable_VEN$RPT, "LKADT_P"],
             DISCONT = core[core$RPT %in% GRtable_VEN$RPT, "DISCONT"],
             DISCONT_day = core[core$RPT %in% GRtable_VEN$RPT, "ENTRT_PC"],
             core_habini[core_habini$RPT %in% GRtable_VEN$RPT, -1],
             core_exclude[core_exclude$RPT %in% GRtable_VEN$RPT, -1])
ASC$DEATH <- ifelse(ASC$DEATH == "YES", 1, 0)
CEL$DEATH <- ifelse(CEL$DEATH == "YES", 1, 0)
VEN$DEATH <- ifelse(VEN$DEATH == "YES", 1, 0)

modelFit <- function(formula, train.data, test.data, model = c("bag", "cox", "lm", "logit", "rf"), 
                     seed, target_column) {
  
  require(ipred)
  require(survival)
  require(randomForest)
  require(PRROC)
  
  set.seed(seed)
  
  fit_model <- switch(
    model,
    "bag" = bagging(formula = formula, data = train.data),
    "lm" = lm(formula = formula, data = train.data),
    "logit" = glm(formula = formula, data = train.data, family = binomial(link = "logit")),
    "rf" = randomForest(formula = formula, data = train.data),
    "cox" = coxph(formula = formula, data = train.data, ties = "breslow")
  )
  
  pred <- predict(fit_model, test.data[, -which(names(test.data) %in% target_column)])
  roc <- roc.curve(scores.class0 = pred, weights.class0 = test.data[, target_column])
  pr <- pr.curve(scores.class0 = pred, weights.class0 = test.data[, target_column], rand.compute = T)
  
  return(list(roc = roc, pr = pr))
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
      death <- bagging(formula = DEATH ~ ., data = train.data[, -which(colnames(train.data) == "DISCONT")])
      discont <- bagging(formula = DISCONT ~ ., data = train.data[, -which(colnames(train.data) == "DEATH")])
      list(death = death, discont = discont)
    },
    "lm" = {
      death <- lm(formula = DEATH ~ ., data = train.data[, -which(colnames(train.data) == "DISCONT")])
      discont <- lm(formula = DISCONT ~ ., data = train.data[, -which(colnames(train.data) == "DEATH")])
      list(death = death, discont = discont)
    },
    "logit" = {
      death <- glm(formula = as.factor(DEATH) ~ ., data = train.data[, -which(colnames(train.data) == "DISCONT")], family = binomial(link = "logit"))
      discont <- glm(formula = as.factor(DISCONT) ~ ., data = train.data[, -which(colnames(train.data) == "DEATH")], family = binomial(link = "logit"))
      list(death = death, discont = discont)
    },
    "rf" = {
      death <- randomForest(formula = DEATH ~ ., data = train.data[, -which(colnames(train.data) == "DISCONT")])
      discont <- randomForest(formula = DISCONT ~ ., data = train.data[, -which(colnames(train.data) == "DEATH")])
      list(death = death, discont = discont)
    },
    "cox" = {
      death <- coxph(formula = Surv(DEATH_day, DEATH) ~ ., data = train.data[, -which(colnames(train.data) %in% c("DISCONT", "DISCONT_day"))], ties = "breslow")
      discont <- coxph(formula = Surv(DISCONT_day, DISCONT) ~ ., data = train.data[, -which(colnames(train.data) %in% c("DEATH", "DEATH_day"))], ties = "breslow")
      list(death = death, discont = discont)
    }
  )
  
  pred_Death <- predict(fit_model$death, test.data[, -which(colnames(test.data) %in% target_column)])
  pred_Discont <- predict(fit_model$discont, test.data[, -which(colnames(test.data) %in% target_column)])
  pred_ensemble <- (pred_Death + pred_Discont) / 2
  roc <- roc.curve(scores.class0 = pred_ensemble, weights.class0 = test.data[, target_column])
  pr <- pr.curve(scores.class0 = pred_ensemble, weights.class0 = test.data[, target_column], rand.compute = T)
  
  return(list(roc = roc, pr = pr))
}


# Cross validation
## ASC + CEL --> VEN --------------
train <- rbind(ASC, CEL)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
train <- na.omit(train)
test <- VEN
test$DISCONT <- as.numeric(as.character(test$DISCONT))
test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
test <- na.omit(test)

# (bag, cox, lm, logit, rf) + 4 features + (Death, Discontinuation, ensemble)
# randomForest Warning:
# Warning message:
# In randomForest.default(m, y, ...) :
#   The response has five or fewer unique values.  Are you sure you want to do regression?
model_bag_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "bag",
                              seed = 1234, target_column = "DISCONT")
model_bag_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "bag",
                                seed = 1234, target_column = "DISCONT")
model_bag_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "bag",
                                seed = 1234, target_column = "DISCONT")
model_cox_4_death <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 11:14)], test[, c(5, 11:14)],
                              model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_4_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 11:14)], test[, c(5, 11:14)],
                                model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_4_ensemble <- modelFitEnsemble(train[, c(3:6, 11:14)], test[, c(5, 11:14)], model = "cox", 
                                        seed = 1234, target_column = "DISCONT")
model_lm_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "lm",
                              seed = 1234, target_column = "DISCONT")
model_lm_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "lm",
                                seed = 1234, target_column = "DISCONT")
model_lm_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "lm",
                               seed = 1234, target_column = "DISCONT")
model_logit_4_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "logit",
                             seed = 1234, target_column = "DISCONT")
model_logit_4_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "logit",
                               seed = 1234, target_column = "DISCONT")
model_logit_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "logit",
                                  seed = 1234, target_column = "DISCONT")
model_rf_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "rf",
                             seed = 1234, target_column = "DISCONT")
model_rf_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "rf",
                               seed = 1234, target_column = "DISCONT")
model_rf_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "rf",
                               seed = 1234, target_column = "DISCONT")

# (bag, cox, lm, logit, rf) + 78 features + (Death, Discontinuation, ensemble)
#
# cox warning:
# Warning messages:
# 1: In fitter(X, Y, strats, offset, init, control, weights = weights,  :
#   Loglik converged before variable  37,38 ; beta may be infinite. 
# Usually you can ignore the message, but the Wald test of significance beta/se(beta) is not valid
# 2: In coxph(formula = formula, data = train.data, ties = "breslow") :
#   X matrix deemed to be singular; variable 21 78
# data in these columns have strong correlation
#
# linear regression warning:
# Warning message:
# In predict.lm(fit_model, test.data[, -which(names(test.data) %in%  :
#   prediction from a rank-deficient fit may be misleading
# This may come from: there is insufficient information contained in your data to estimate the model you desire
#
# randomForest Warning:
# Warning message:
# In randomForest.default(m, y, ...) :
#   The response has five or fewer unique values.  Are you sure you want to do regression?

model_bag_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "bag",
                              seed = 1234, target_column = "DISCONT")
model_bag_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "bag",
                                seed = 1234, target_column = "DISCONT")
model_bag_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "bag",
                                         seed = 1234, target_column = "DISCONT")
model_cox_78_death <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 7:84)], test[, c(5, 7:84)],
                              model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_78_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 7:84)], test[, c(5, 7:84)],
                                model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_78_ensemble <- modelFitEnsemble(train[, c(3:6, 7:84)], test[, c(5, 7:84)],
                                 model = "cox", seed = 1234, target_column = "DISCONT")
model_lm_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "lm",
                             seed = 1234, target_column = "DISCONT")
model_lm_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "lm",
                               seed = 1234, target_column = "DISCONT")
model_lm_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "lm",
                                seed = 1234, target_column = "DISCONT")
model_logit_78_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "logit",
                                seed = 1234, target_column = "DISCONT")
model_logit_78_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                  seed = 1234, target_column = "DISCONT")
model_logit_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                   seed = 1234, target_column = "DISCONT")
model_rf_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "rf",
                             seed = 1234, target_column = "DISCONT")
model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
                               seed = 1234, target_column = "DISCONT")
model_rf_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                seed = 1234, target_column = "DISCONT")

# Performance failed to fullfill the expectation
# # (BAG, linear, rf) + 4 features + new discont value (generated with GuanRank)
# model_bag_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "bag",
#                                   seed = 1234, target_column = "DISCONT")
# model_lm_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "lm",
#                                  seed = 1234, target_column = "DISCONT")
# model_rf_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "rf",
#                                   seed = 1234, target_column = "DISCONT")
# 
# # (BAG, linear, rf) + 78 features + new discont value (generated with GuanRank)
# model_bag_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "bag",
#                                    seed = 1234, target_column = "DISCONT")
# model_lm_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "lm",
#                                    seed = 1234, target_column = "DISCONT")
# model_rf_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
#                                    seed = 1234, target_column = "DISCONT")

res <- data.frame(identity = rep(c("model_bag_4_death", "model_bag_4_discont", "model_bag_4_ensemble",
                               "model_bag_78_death", "model_bag_78_discont", "model_bag_78_ensemble",
                               
                               "model_lm_4_death", "model_lm_4_discont", "model_lm_4_ensemble",
                               "model_lm_78_death", "model_lm_78_discont", "model_lm_78_ensemble",
                               
                               "model_rf_4_death", "model_rf_4_discont", "model_rf_4_ensemble",
                               "model_rf_78_death", "model_rf_78_discont", "model_rf_78_ensemble",
                               
                               "model_cox_4_death", "model_cox_4_discont", "model_cox_4_ensemble",
                               "model_cox_78_death", "model_cox_78_discont", "model_cox_78_ensemble",
                               
                               "model_logit_4_death", "model_logit_4_discont", "model_logit_4_ensemble",
                               "model_logit_78_death", "model_logit_78_discont", "model_logit_78_ensemble"),
                               each = 2),
                  value = c(model_bag_4_death$roc$auc, model_bag_4_death$pr$auc.integral,
                            model_bag_4_discont$roc$auc, model_bag_4_discont$pr$auc.integral,
                            model_bag_4_ensemble$roc$auc, model_bag_4_ensemble$pr$auc.integral,
                            model_bag_78_death$roc$auc, model_bag_78_death$pr$auc.integral,
                            model_bag_78_discont$roc$auc, model_bag_78_discont$pr$auc.integral,
                            model_bag_78_ensemble$roc$auc, model_bag_78_ensemble$pr$auc.integral,
                            
                            model_lm_4_death$roc$auc, model_lm_4_death$pr$auc.integral,
                            model_lm_4_discont$roc$auc, model_lm_4_discont$pr$auc.integral,
                            model_lm_4_ensemble$roc$auc, model_lm_4_ensemble$pr$auc.integral,
                            model_lm_78_death$roc$auc, model_lm_78_death$pr$auc.integral,
                            model_lm_78_discont$roc$auc, model_lm_78_discont$pr$auc.integral,
                            model_lm_78_ensemble$roc$auc, model_lm_78_ensemble$pr$auc.integral,
                            
                            model_rf_4_death$roc$auc, model_rf_4_death$pr$auc.integral,
                            model_rf_4_discont$roc$auc, model_rf_4_discont$pr$auc.integral,
                            model_rf_4_ensemble$roc$auc, model_rf_4_ensemble$pr$auc.integral,
                            model_rf_78_death$roc$auc, model_rf_78_death$pr$auc.integral,
                            model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral,
                            model_rf_78_ensemble$roc$auc, model_rf_78_ensemble$pr$auc.integral,
                            
                            model_cox_4_death$roc$auc, model_cox_4_death$pr$auc.integral,
                            model_cox_4_discont$roc$auc, model_cox_4_discont$pr$auc.integral,
                            model_cox_4_ensemble$roc$auc, model_cox_4_ensemble$pr$auc.integral,
                            model_cox_78_death$roc$auc, model_cox_78_death$pr$auc.integral,
                            model_cox_78_discont$roc$auc, model_cox_78_discont$pr$auc.integral,
                            model_cox_78_ensemble$roc$auc, model_cox_78_ensemble$pr$auc.integral,
                            
                            model_logit_4_death$roc$auc, model_logit_4_death$pr$auc.integral,
                            model_logit_4_discont$roc$auc, model_logit_4_discont$pr$auc.integral,
                            model_logit_4_ensemble$roc$auc, model_logit_4_discont$pr$auc.integral,
                            model_logit_78_death$roc$auc, model_logit_78_death$pr$auc.integral,
                            model_logit_78_discont$roc$auc, model_logit_78_discont$pr$auc.integral,
                            model_logit_78_ensemble$roc$auc, model_logit_78_ensemble$pr$auc.integral),
                  stringsAsFactors = F)

res$model <- sapply(res$identity, function(x) {
  if (grepl("bag", x)) return("BAG-CART")
  if (grepl("cox", x)) return("Cox")
  if (grepl("lm", x)) return("Linear")
  if (grepl("logit", x)) return("Logistic")
  if (grepl("rf", x)) return("RF")
})
res$feature_num <- unlist(lapply(strsplit(res$identity, "_"), "[", 3))
res$gold_standard <- unlist(lapply(strsplit(res$identity, "_"), "[", 4))
res$model_feature <- paste0(res$model, " (", res$feature_num, ")")
res$curve <- rep(c("AUC", "AUPRC"), 30)
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
  geom_hline(yintercept = model_bag_4_death$pr$rand$auc.integral, linetype = "dashed")

## ASC + VEN --> CEL ----------------
train <- rbind(ASC, VEN)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
train <- na.omit(train)
test <- CEL
test$DISCONT <- as.numeric(as.character(test$DISCONT))
test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
test <- na.omit(test)

# (bag, cox, lm, logit, rf) + 4 features + (Death, Discontinuation)
# randomForest Warning:
# Warning message:
# In randomForest.default(m, y, ...) :
#   The response has five or fewer unique values.  Are you sure you want to do regression?
model_bag_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "bag",
                              seed = 1234, target_column = "DISCONT")
model_bag_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "bag",
                                seed = 1234, target_column = "DISCONT")
model_bag_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "bag",
                                         seed = 1234, target_column = "DISCONT")
model_cox_4_death <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 11:14)], test[, c(5, 11:14)],
                              model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_4_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 11:14)], test[, c(5, 11:14)],
                                model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_4_ensemble <- modelFitEnsemble(train[, c(3:6, 11:14)], test[, c(5, 11:14)],
                                model = "cox", seed = 1234, target_column = "DISCONT")
model_lm_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "lm",
                             seed = 1234, target_column = "DISCONT")
model_lm_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "lm",
                               seed = 1234, target_column = "DISCONT")
model_lm_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "lm",
                                seed = 1234, target_column = "DISCONT")
model_logit_4_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "logit",
                                seed = 1234, target_column = "DISCONT")
model_logit_4_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "logit",
                                  seed = 1234, target_column = "DISCONT")
model_logit_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "logit",
                                  seed = 1234, target_column = "DISCONT")
model_rf_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "rf",
                             seed = 1234, target_column = "DISCONT")
model_rf_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "rf",
                               seed = 1234, target_column = "DISCONT")
model_rf_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "rf",
                               seed = 1234, target_column = "DISCONT")

# (bag, cox, lm, logit, rf) + 78 features + (Death, Discontinuation)
#
# cox warning:
# Warning messages:
# 1: In fitter(X, Y, strats, offset, init, control, weights = weights,  :
#   Loglik converged before variable  37,38 ; beta may be infinite. 
# Usually you can ignore the message, but the Wald test of significance beta/se(beta) is not valid
# 2: In coxph(formula = formula, data = train.data, ties = "breslow") :
#   X matrix deemed to be singular; variable 21 78
# data in these columns have strong correlation
#
# linear regression warning:
# Warning message:
# In predict.lm(fit_model, test.data[, -which(names(test.data) %in%  :
#   prediction from a rank-deficient fit may be misleading
# This may come from: there is insufficient information contained in your data to estimate the model you desire
#
# randomForest Warning:
# Warning message:
# In randomForest.default(m, y, ...) :
#   The response has five or fewer unique values.  Are you sure you want to do regression?
#
# logistic regression warning:
# Warning messages:
#   1: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# labels occure at a certain side of values


model_bag_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "bag",
                               seed = 1234, target_column = "DISCONT")
model_bag_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "bag",
                                 seed = 1234, target_column = "DISCONT")
model_bag_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "bag",
                                 seed = 1234, target_column = "DISCONT")
model_cox_78_death <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 7:84)], test[, c(5, 7:84)],
                               model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_78_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 7:84)], test[, c(5, 7:84)],
                                 model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_78_ensemble <- modelFitEnsemble(train[, c(3:6, 7:84)], test[, c(5, 7:84)],
                                 model = "cox", seed = 1234, target_column = "DISCONT")
model_lm_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "lm",
                              seed = 1234, target_column = "DISCONT")
model_lm_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "lm",
                                seed = 1234, target_column = "DISCONT")
model_lm_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "lm",
                                seed = 1234, target_column = "DISCONT")
model_logit_78_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "logit",
                                 seed = 1234, target_column = "DISCONT")
model_logit_78_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                   seed = 1234, target_column = "DISCONT")
model_logit_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                   seed = 1234, target_column = "DISCONT")
model_rf_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "rf",
                              seed = 1234, target_column = "DISCONT")
model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                seed = 1234, target_column = "DISCONT")
model_rf_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                seed = 1234, target_column = "DISCONT")

# # (BAG, linear, rf) + 4 features + new discont value (generated with GuanRank)
# model_bag_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "bag",
#                                   seed = 1234, target_column = "DISCONT")
# model_lm_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "lm",
#                                  seed = 1234, target_column = "DISCONT")
# model_rf_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "rf",
#                                  seed = 1234, target_column = "DISCONT")

# # (BAG, linear, rf) + 78 features + new discont value (generated with GuanRank)
# model_bag_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "bag",
#                                    seed = 1234, target_column = "DISCONT")
# model_lm_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "lm",
#                                   seed = 1234, target_column = "DISCONT")
# model_rf_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
#                                   seed = 1234, target_column = "DISCONT")

res1 <- data.frame(identity = rep(c("model_bag_4_death", "model_bag_4_discont", "model_bag_4_ensemble",
                                   "model_bag_78_death", "model_bag_78_discont", "model_bag_78_ensemble",
                                   
                                   "model_lm_4_death", "model_lm_4_discont", "model_lm_4_ensemble",
                                   "model_lm_78_death", "model_lm_78_discont", "model_lm_78_ensemble",
                                   
                                   "model_rf_4_death", "model_rf_4_discont", "model_rf_4_ensemble",
                                   "model_rf_78_death", "model_rf_78_discont", "model_rf_78_ensemble",
                                   
                                   "model_cox_4_death", "model_cox_4_discont", "model_cox_4_ensemble",
                                   "model_cox_78_death", "model_cox_78_discont", "model_cox_78_ensemble",
                                   
                                   "model_logit_4_death", "model_logit_4_discont", "model_logit_4_ensemble",
                                   "model_logit_78_death", "model_logit_78_discont", "model_logit_78_ensemble"),
                                 each = 2),
                  value = c(model_bag_4_death$roc$auc, model_bag_4_death$pr$auc.integral,
                            model_bag_4_discont$roc$auc, model_bag_4_discont$pr$auc.integral,
                            model_bag_4_ensemble$roc$auc, model_bag_4_ensemble$pr$auc.integral,
                            model_bag_78_death$roc$auc, model_bag_78_death$pr$auc.integral,
                            model_bag_78_discont$roc$auc, model_bag_78_discont$pr$auc.integral,
                            model_bag_78_ensemble$roc$auc, model_bag_78_ensemble$pr$auc.integral,
                            
                            model_lm_4_death$roc$auc, model_lm_4_death$pr$auc.integral,
                            model_lm_4_discont$roc$auc, model_lm_4_discont$pr$auc.integral,
                            model_lm_4_ensemble$roc$auc, model_lm_4_ensemble$pr$auc.integral,
                            model_lm_78_death$roc$auc, model_lm_78_death$pr$auc.integral,
                            model_lm_78_discont$roc$auc, model_lm_78_discont$pr$auc.integral,
                            model_lm_78_ensemble$roc$auc, model_lm_78_ensemble$pr$auc.integral,
                            
                            model_rf_4_death$roc$auc, model_rf_4_death$pr$auc.integral,
                            model_rf_4_discont$roc$auc, model_rf_4_discont$pr$auc.integral,
                            model_rf_4_ensemble$roc$auc, model_rf_4_ensemble$pr$auc.integral,
                            model_rf_78_death$roc$auc, model_rf_78_death$pr$auc.integral,
                            model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral,
                            model_rf_78_ensemble$roc$auc, model_rf_78_ensemble$pr$auc.integral,
                            
                            model_cox_4_death$roc$auc, model_cox_4_death$pr$auc.integral,
                            model_cox_4_discont$roc$auc, model_cox_4_discont$pr$auc.integral,
                            model_cox_4_ensemble$roc$auc, model_cox_4_ensemble$pr$auc.integral,
                            model_cox_78_death$roc$auc, model_cox_78_death$pr$auc.integral,
                            model_cox_78_discont$roc$auc, model_cox_78_discont$pr$auc.integral,
                            model_cox_78_ensemble$roc$auc, model_cox_78_ensemble$pr$auc.integral,
                            
                            model_logit_4_death$roc$auc, model_logit_4_death$pr$auc.integral,
                            model_logit_4_discont$roc$auc, model_logit_4_discont$pr$auc.integral,
                            model_logit_4_ensemble$roc$auc, model_logit_4_discont$pr$auc.integral,
                            model_logit_78_death$roc$auc, model_logit_78_death$pr$auc.integral,
                            model_logit_78_discont$roc$auc, model_logit_78_discont$pr$auc.integral,
                            model_logit_78_ensemble$roc$auc, model_logit_78_ensemble$pr$auc.integral),
                  stringsAsFactors = F)

res1$model <- sapply(res1$identity, function(x) {
  if (grepl("bag", x)) return("BAG-CART")
  if (grepl("cox", x)) return("Cox")
  if (grepl("lm", x)) return("Linear")
  if (grepl("logit", x)) return("Logistic")
  if (grepl("rf", x)) return("RF")
})
res1$feature_num <- unlist(lapply(strsplit(res1$identity, "_"), "[", 3))
res1$gold_standard <- unlist(lapply(strsplit(res1$identity, "_"), "[", 4))
res1$model_feature <- paste0(res1$model, " (", res1$feature_num, ")")
res1$curve <- rep(c("AUC", "AUPRC"), 30)
p1_AUC <- ggplot(subset(res1, curve == "AUC"))
p1_AUC <- p1_AUC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("ASC + VEN --> CEL") + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_color_discrete("Predictor") + 
  scale_y_continuous(limits = c(0.4, 0.7))
p1_AUPRC <- ggplot(subset(res1, curve == "AUPRC"))
p1_AUPRC <- p1_AUPRC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("ASC + VEN --> CEL") + 
  scale_color_discrete("Predictor") + scale_y_continuous(limits = c(0, 0.3)) +
  geom_hline(yintercept = model_bag_4_death$pr$rand$auc.integral, linetype = "dashed")


## VEN + CEL --> ASC --------------------
train <- rbind(CEL, VEN)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train$DISCONT_day <- as.numeric(as.character(train$DISCONT_day))
train <- na.omit(train)
test <- ASC
test$DISCONT <- as.numeric(as.character(test$DISCONT))
test$DISCONT_day <- as.numeric(as.character(test$DISCONT_day))
test <- na.omit(test)

# (bag, cox, lm, logit, rf) + 4 features + (Death, Discontinuation)
# randomForest Warning:
# Warning message:
# In randomForest.default(m, y, ...) :
#   The response has five or fewer unique values.  Are you sure you want to do regression?
model_bag_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "bag",
                              seed = 1234, target_column = "DISCONT")
model_bag_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "bag",
                                seed = 1234, target_column = "DISCONT")
model_bag_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "bag",
                                seed = 1234, target_column = "DISCONT")
model_cox_4_death <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 11:14)], test[, c(5, 11:14)],
                              model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_4_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 11:14)], test[, c(5, 11:14)],
                                model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_4_ensemble <- modelFitEnsemble(train[, c(3:6, 11:14)], test[, c(5, 11:14)],
                                model = "cox", seed = 1234, target_column = "DISCONT")
model_lm_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "lm",
                             seed = 1234, target_column = "DISCONT")
model_lm_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "lm",
                               seed = 1234, target_column = "DISCONT")
model_lm_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "lm",
                               seed = 1234, target_column = "DISCONT")
model_logit_4_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "logit",
                                seed = 1234, target_column = "DISCONT")
model_logit_4_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "logit",
                                  seed = 1234, target_column = "DISCONT")
model_logit_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "logit",
                                  seed = 1234, target_column = "DISCONT")
model_rf_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "rf",
                             seed = 1234, target_column = "DISCONT")
model_rf_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "rf",
                               seed = 1234, target_column = "DISCONT")
model_rf_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "rf",
                               seed = 1234, target_column = "DISCONT")
# (bag, cox, lm, logit, rf) + 78 features + (Death, Discontinuation)
#
# cox warning:
# Warning messages:
# 1: In fitter(X, Y, strats, offset, init, control, weights = weights,  :
#   Loglik converged before variable  37,38 ; beta may be infinite. 
# Usually you can ignore the message, but the Wald test of significance beta/se(beta) is not valid
# 2: In coxph(formula = formula, data = train.data, ties = "breslow") :
#   X matrix deemed to be singular; variable 21 78
# data in these columns have strong correlation
#
# linear regression warning:
# Warning message:
# In predict.lm(fit_model, test.data[, -which(names(test.data) %in%  :
#   prediction from a rank-deficient fit may be misleading
# This may come from: there is insufficient information contained in your data to estimate the model you desire
#
# randomForest Warning:
# Warning message:
# In randomForest.default(m, y, ...) :
#   The response has five or fewer unique values.  Are you sure you want to do regression?
#
# logistic regression warning:
# Warning messages:
#   1: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# labels occure at a certain side of values


model_bag_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "bag",
                               seed = 1234, target_column = "DISCONT")
model_bag_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "bag",
                                 seed = 1234, target_column = "DISCONT")
model_bag_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "bag",
                                 seed = 1234, target_column = "DISCONT")
model_cox_78_death <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 7:84)], test[, c(5, 7:84)],
                               model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_78_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 7:84)], test[, c(5, 7:84)],
                                 model = "cox", seed = 1234, target_column = "DISCONT")
model_cox_78_ensemble <- modelFitEnsemble(train[, c(3:6, 7:84)], test[, c(5, 7:84)],
                                 model = "cox", seed = 1234, target_column = "DISCONT")
model_lm_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "lm",
                              seed = 1234, target_column = "DISCONT")
model_lm_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "lm",
                                seed = 1234, target_column = "DISCONT")
model_lm_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "lm",
                                seed = 1234, target_column = "DISCONT")
model_logit_78_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "logit",
                                 seed = 1234, target_column = "DISCONT")
model_logit_78_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                   seed = 1234, target_column = "DISCONT")
model_logit_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                   seed = 1234, target_column = "DISCONT")
model_rf_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "rf",
                              seed = 1234, target_column = "DISCONT")
model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                seed = 1234, target_column = "DISCONT")
model_rf_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "rf",
                                seed = 1234, target_column = "DISCONT")
# # (BAG, linear, rf) + 4 features + new discont value (generated with GuanRank)
# model_bag_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "bag",
#                                   seed = 1234, target_column = "DISCONT")
# model_lm_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "lm",
#                                  seed = 1234, target_column = "DISCONT")
# model_rf_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "rf",
#                                  seed = 1234, target_column = "DISCONT")
# 
# # (BAG, linear, rf) + 78 features + new discont value (generated with GuanRank)
# model_bag_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "bag",
#                                    seed = 1234, target_column = "DISCONT")
# model_lm_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "lm",
#                                   seed = 1234, target_column = "DISCONT")
# model_rf_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
#                                   seed = 1234, target_column = "DISCONT")

res2 <- data.frame(identity = rep(c("model_bag_4_death", "model_bag_4_discont", "model_bag_4_ensemble",
                                    "model_bag_78_death", "model_bag_78_discont", "model_bag_78_ensemble",
                                    
                                    "model_lm_4_death", "model_lm_4_discont", "model_lm_4_ensemble",
                                    "model_lm_78_death", "model_lm_78_discont", "model_lm_78_ensemble",
                                    
                                    "model_rf_4_death", "model_rf_4_discont", "model_rf_4_ensemble",
                                    "model_rf_78_death", "model_rf_78_discont", "model_rf_78_ensemble",
                                    
                                    "model_cox_4_death", "model_cox_4_discont", "model_cox_4_ensemble",
                                    "model_cox_78_death", "model_cox_78_discont", "model_cox_78_ensemble",
                                    
                                    "model_logit_4_death", "model_logit_4_discont", "model_logit_4_ensemble",
                                    "model_logit_78_death", "model_logit_78_discont", "model_logit_78_ensemble"),
                                  each = 2),
                   value = c(model_bag_4_death$roc$auc, model_bag_4_death$pr$auc.integral,
                             model_bag_4_discont$roc$auc, model_bag_4_discont$pr$auc.integral,
                             model_bag_4_ensemble$roc$auc, model_bag_4_ensemble$pr$auc.integral,
                             model_bag_78_death$roc$auc, model_bag_78_death$pr$auc.integral,
                             model_bag_78_discont$roc$auc, model_bag_78_discont$pr$auc.integral,
                             model_bag_78_ensemble$roc$auc, model_bag_78_ensemble$pr$auc.integral,
                             
                             model_lm_4_death$roc$auc, model_lm_4_death$pr$auc.integral,
                             model_lm_4_discont$roc$auc, model_lm_4_discont$pr$auc.integral,
                             model_lm_4_ensemble$roc$auc, model_lm_4_ensemble$pr$auc.integral,
                             model_lm_78_death$roc$auc, model_lm_78_death$pr$auc.integral,
                             model_lm_78_discont$roc$auc, model_lm_78_discont$pr$auc.integral,
                             model_lm_78_ensemble$roc$auc, model_lm_78_ensemble$pr$auc.integral,
                             
                             model_rf_4_death$roc$auc, model_rf_4_death$pr$auc.integral,
                             model_rf_4_discont$roc$auc, model_rf_4_discont$pr$auc.integral,
                             model_rf_4_ensemble$roc$auc, model_rf_4_ensemble$pr$auc.integral,
                             model_rf_78_death$roc$auc, model_rf_78_death$pr$auc.integral,
                             model_rf_78_discont$roc$auc, model_rf_78_discont$pr$auc.integral,
                             model_rf_78_ensemble$roc$auc, model_rf_78_ensemble$pr$auc.integral,
                             
                             model_cox_4_death$roc$auc, model_cox_4_death$pr$auc.integral,
                             model_cox_4_discont$roc$auc, model_cox_4_discont$pr$auc.integral,
                             model_cox_4_ensemble$roc$auc, model_cox_4_ensemble$pr$auc.integral,
                             model_cox_78_death$roc$auc, model_cox_78_death$pr$auc.integral,
                             model_cox_78_discont$roc$auc, model_cox_78_discont$pr$auc.integral,
                             model_cox_78_ensemble$roc$auc, model_cox_78_ensemble$pr$auc.integral,
                             
                             model_logit_4_death$roc$auc, model_logit_4_death$pr$auc.integral,
                             model_logit_4_discont$roc$auc, model_logit_4_discont$pr$auc.integral,
                             model_logit_4_ensemble$roc$auc, model_logit_4_discont$pr$auc.integral,
                             model_logit_78_death$roc$auc, model_logit_78_death$pr$auc.integral,
                             model_logit_78_discont$roc$auc, model_logit_78_discont$pr$auc.integral,
                             model_logit_78_ensemble$roc$auc, model_logit_78_ensemble$pr$auc.integral),
                   stringsAsFactors = F)

res2$model <- sapply(res2$identity, function(x) {
  if (grepl("bag", x)) return("BAG-CART")
  if (grepl("cox", x)) return("Cox")
  if (grepl("lm", x)) return("Linear")
  if (grepl("logit", x)) return("Logistic")
  if (grepl("rf", x)) return("RF")
})
res2$feature_num <- unlist(lapply(strsplit(res2$identity, "_"), "[", 3))
res2$gold_standard <- unlist(lapply(strsplit(res2$identity, "_"), "[", 4))
res2$model_feature <- paste0(res2$model, " (", res2$feature_num, ")")
res2$curve <- rep(c("AUC", "AUPRC"), 30)
p2_AUC <- ggplot(subset(res2, curve == "AUC"))
p2_AUC <- p2_AUC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("CEL + VEN --> ASC") + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_color_discrete("Predictor") + 
  scale_y_continuous(limits = c(0.4, 0.7))
p2_AUPRC <- ggplot(subset(res2, curve == "AUPRC"))
p2_AUPRC <- p2_AUPRC + geom_point(aes(model_feature, value, color = gold_standard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve") + ggtitle("CEL + VEN --> ASC") + 
  scale_color_discrete("Predictor") + scale_y_continuous(limits = c(0, 0.3)) + 
  geom_hline(yintercept = model_bag_4_death$pr$rand$auc.integral, linetype = "dashed")




# ---------------

grid_arrange_share_legend <- function(..., nrow, ncol, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(
    position,
    "bottom" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight)
    ),
    "right" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 2,
      widths = unit.c(unit(1, "npc") - lwidth, lwidth)
    )
  )
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}
grid_arrange_share_legend(p_AUC, p1_AUC, p2_AUC, nrow = 1, ncol = 3, position = "right")
grid_arrange_share_legend(p_AUPRC, p1_AUPRC, p2_AUPRC, nrow = 1, ncol = 3, position = "right")

# In Cohort cv, 25 times -------------- 
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
    
    # model_bag_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "bag",
    #                               seed = 1234, target_column = "DISCONT")
    # model_bag_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "bag",
    #                                 seed = 1234, target_column = "DISCONT")
    # model_bag_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "bag",
    #                                 seed = 1234, target_column = "DISCONT")
    model_cox_4_death <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 11:14)], test[, c(5, 11:14)],
                                  model = "cox", seed = 1234, target_column = "DISCONT")
    model_cox_4_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 11:14)], test[, c(5, 11:14)],
                                    model = "cox", seed = 1234, target_column = "DISCONT")
    model_cox_9_discont <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 7:15)], test[, c(5, 7:15)],
                                    model = "cox", seed = 1234, target_column = "DISCONT")
    model_cox_9_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 7:15)], test[, c(5, 7:15)],
                                    model = "cox", seed = 1234, target_column = "DISCONT")
    # model_cox_4_ensemble <- modelFitEnsemble(train[, c(3:6, 11:14)], test[, c(5, 11:14)],
    #                                 model = "cox", seed = 1234, target_column = "DISCONT")
    model_lm_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "lm",
                                 seed = 1234, target_column = "DISCONT")
    model_lm_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "lm",
                                   seed = 1234, target_column = "DISCONT")
    model_lm_9_death <- modelFit(DEATH ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "lm",
                                 seed = 1234, target_column = "DISCONT")
    model_lm_9_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "lm",
                                   seed = 1234, target_column = "DISCONT")
    # model_lm_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "lm",
    #                                seed = 1234, target_column = "DISCONT")
    model_logit_4_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "logit",
                                    seed = 1234, target_column = "DISCONT")
    model_logit_4_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "logit",
                                      seed = 1234, target_column = "DISCONT")
    model_logit_9_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:15)], test[, c(5, 7:15)], model = "logit",
                                    seed = 1234, target_column = "DISCONT")
    model_logit_9_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:15)], test[, c(5, 7:15)], model = "logit",
                                      seed = 1234, target_column = "DISCONT")
    # model_logit_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "logit",
    #                                   seed = 1234, target_column = "DISCONT")
    # model_rf_4_death <- modelFit(DEATH ~ ., train[, c(3, 11:14)], test[, c(5, 11:14)], model = "rf",
    #                              seed = 1234, target_column = "DISCONT")
    # model_rf_4_discont <- modelFit(DISCONT ~ ., train[, c(5, 11:14)], test[, c(5, 11:14)], model = "rf",
    #                                seed = 1234, target_column = "DISCONT")
    # model_rf_4_ensemble <- modelFitEnsemble(train[, c(3, 5, 11:14)], test[, c(5, 11:14)], model = "rf",
    #                                seed = 1234, target_column = "DISCONT")
    # model_bag_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "bag",
    #                                seed = 1234, target_column = "DISCONT")
    # model_bag_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "bag",
    #                                  seed = 1234, target_column = "DISCONT")
    # model_bag_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "bag",
    #                                  seed = 1234, target_column = "DISCONT")
    model_cox_78_death <- modelFit(Surv(DEATH_day, DEATH) ~ ., train[, c(3:4, 7:84)], test[, c(5, 7:84)],
                                   model = "cox", seed = 1234, target_column = "DISCONT")
    model_cox_78_discont <- modelFit(Surv(DISCONT_day, DISCONT) ~ ., train[, c(5:6, 7:84)], test[, c(5, 7:84)],
                                     model = "cox", seed = 1234, target_column = "DISCONT")
    # model_cox_78_ensemble <- modelFitEnsemble(train[, c(3:6, 7:84)], test[, c(5, 7:84)],
    #                                  model = "cox", seed = 1234, target_column = "DISCONT")
    model_lm_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "lm",
                                  seed = 1234, target_column = "DISCONT")
    model_lm_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "lm",
                                    seed = 1234, target_column = "DISCONT")
    # model_lm_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "lm",
    #                                 seed = 1234, target_column = "DISCONT")
    model_logit_78_death <- modelFit(as.factor(DEATH) ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "logit",
                                     seed = 1234, target_column = "DISCONT")
    model_logit_78_discont <- modelFit(as.factor(DISCONT) ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "logit",
                                       seed = 1234, target_column = "DISCONT")
    # model_logit_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "logit",
    #                                    seed = 1234, target_column = "DISCONT")
    # model_rf_78_death <- modelFit(DEATH ~ ., train[, c(3, 7:84)], test[, c(5, 7:84)], model = "rf",
    #                               seed = 1234, target_column = "DISCONT")
    # model_rf_78_discont <- modelFit(DISCONT ~ ., train[, c(5, 7:84)], test[, c(5, 7:84)], model = "rf",
    #                                 seed = 1234, target_column = "DISCONT")
    # model_rf_78_ensemble <- modelFitEnsemble(train[, c(3, 5, 7:84)], test[, c(5, 7:84)], model = "rf",
    #                                 seed = 1234, target_column = "DISCONT")
    # model_bag_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "bag",
    #                                   seed = 1234, target_column = "DISCONT")
    # model_lm_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "lm",
    #                                  seed = 1234, target_column = "DISCONT")
    # model_rf_4_discontGR <- modelFit(rankValue ~ ., train[, c(2, 11:14)], test[, c(5, 11:14)], model = "rf",
    #                                  seed = 1234, target_column = "DISCONT")
    # model_bag_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "bag",
    #                                    seed = 1234, target_column = "DISCONT")
    # model_lm_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "lm",
    #                                   seed = 1234, target_column = "DISCONT")
    # model_rf_78_discontGR <- modelFit(rankValue ~ ., train[, c(2, 7:84)], test[, c(5, 7:84)], model = "rf",
    #                                   seed = 1234, target_column = "DISCONT")
    
    res <- data.frame(identity = rep(c("model_lm_4_death", "model_lm_4_discont", 
                                       "model_lm_9_death", "model_lm_9_discont",
                                       "model_lm_78_death", "model_lm_78_discont",
                                       
                                       "model_cox_4_death", "model_cox_4_discont", 
                                       "model_cox_9_death", "model_cox_9_discont",
                                       "model_cox_78_death", "model_cox_78_discont", 
                                       
                                       "model_logit_4_death", "model_logit_4_discont", 
                                       "model_logit_9_death", "model_logit_9_discont",
                                       "model_logit_78_death", "model_logit_78_discont"),
                                     each = 2),
                      value = c(model_lm_4_death$roc$auc, model_lm_4_death$pr$auc.integral,
                                model_lm_4_discont$roc$auc, model_lm_4_discont$pr$auc.integral,
                                model_lm_9_death$roc$auc, model_lm_9_death$pr$auc.integral,
                                model_lm_9_discont$roc$auc, model_lm_9_discont$pr$auc.integral,
                                model_lm_78_death$roc$auc, model_lm_78_death$pr$auc.integral,
                                model_lm_78_discont$roc$auc, model_lm_78_discont$pr$auc.integral,
                                
                                model_cox_4_death$roc$auc, model_cox_4_death$pr$auc.integral,
                                model_cox_4_discont$roc$auc, model_cox_4_discont$pr$auc.integral,
                                model_cox_9_death$roc$auc, model_cox_9_death$pr$auc.integral,
                                model_cox_9_discont$roc$auc, model_cox_9_discont$pr$auc.integral,
                                model_cox_78_death$roc$auc, model_cox_78_death$pr$auc.integral,
                                model_cox_78_discont$roc$auc, model_cox_78_discont$pr$auc.integral,
                              
                                model_logit_4_death$roc$auc, model_logit_4_death$pr$auc.integral,
                                model_logit_4_discont$roc$auc, model_logit_4_discont$pr$auc.integral,
                                model_logit_9_death$roc$auc, model_logit_9_death$pr$auc.integral,
                                model_logit_9_discont$roc$auc, model_logit_9_discont$pr$auc.integral,
                                model_logit_78_death$roc$auc, model_logit_78_death$pr$auc.integral,
                                model_logit_78_discont$roc$auc, model_logit_78_discont$pr$auc.integral),
                      stringsAsFactors = F)
    
    res$model <- sapply(res$identity, function(x) {
      
      if (grepl("cox", x)) return("Cox")
      if (grepl("lm", x)) return("Linear")
      if (grepl("logit", x)) return("Logistic")
      
    })
    res$feature_num <- unlist(lapply(strsplit(res$identity, "_"), "[", 3))
    res$gold_standard <- unlist(lapply(strsplit(res$identity, "_"), "[", 4))
    res$model_feature <- paste0(res$model, " (", res$feature_num, ")")
    res$curve <- rep(c("AUC", "AUPRC"), 18)
    
    result <- rbind(result, res)
  }
  return(result)
}
ASCinCohortCV <- data.frame()
CELinCohortCV <- data.frame()
VENinCohortCV <- data.frame()
for (i in 1:5) {
  ASCinCohortCV <- rbind(ASCinCohortCV, kFoldCV(ASC, 5, i))
  CELinCohortCV <- rbind(CELinCohortCV, kFoldCV(CEL, 5, i))
  VENinCohortCV <- rbind(VENinCohortCV, kFoldCV(VEN, 5, i))
}
p_AUC <- ggplot(subset(ASCinCohortCV, curve == "AUC"))
p_AUC <- p_AUC + geom_boxplot(aes(model_feature, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve", title = "ASC") + scale_y_continuous(limits = c(0.2, 0.8)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed")
p_AUPRC <- ggplot(subset(ASCinCohortCV, curve == "AUPRC"))
p_AUPRC <- p_AUPRC + geom_boxplot(aes(model_feature, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve", title = "ASC") + scale_y_continuous(limits = c(0, 0.5))

p1_AUC <- ggplot(subset(CELinCohortCV, curve == "AUC"))
p1_AUC <- p1_AUC + geom_boxplot(aes(model_feature, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve", title = "CEL") + scale_y_continuous(limits = c(0.2, 0.8)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed")
p1_AUPRC <- ggplot(subset(CELinCohortCV, curve == "AUPRC"))
p1_AUPRC <- p1_AUPRC + geom_boxplot(aes(model_feature, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve", title = "CEL") + scale_y_continuous(limits = c(0, 0.5))

p2_AUC <- ggplot(subset(VENinCohortCV, curve == "AUC"))
p2_AUC <- p2_AUC + geom_boxplot(aes(model_feature, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve", title = "VEN") + scale_y_continuous(limits = c(0.2, 0.8)) +
  geom_hline(yintercept = 0.5, linetype = "dashed")
p2_AUPRC <- ggplot(subset(VENinCohortCV, curve == "AUPRC"))
p2_AUPRC <- p2_AUPRC + geom_boxplot(aes(model_feature, value, fill = gold_standard), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + scale_fill_discrete("Predictor") + 
  labs(x = "Model (Feature Number)", y = "Area Under Curve", title = "VEN") + scale_y_continuous(limits = c(0, 0.5))

grid_arrange_share_legend(p_AUC, p1_AUC, p2_AUC, nrow = 1, ncol = 3, position = "right")
grid_arrange_share_legend(p_AUPRC, p1_AUPRC, p2_AUPRC, nrow = 1, ncol = 3, position = "right")

cor(death_3_standard$V2, as.numeric(as.character(core[core$RPT %in% death_3_standard$V1, "DISCONT"])))
