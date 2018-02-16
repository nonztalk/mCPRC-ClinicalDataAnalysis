library(ipred)
library(PRROC)
library(randomForest)

## 1.1 Halabi model parameters
train <- rbind(habini_ASC_death, habini_CEL_death)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_VEN_death
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

## 1.1.1 choose four variables as orignial matlab file said
## 1.1.1.1 use linear regression model as the original
Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 7:10])
pred_Discont <- predict(Fit_Discont, test[, 7:10])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## 1.1.1.2 use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 7:10])
pred_Discont <- predict(Fit_Discont, test[, 7:10])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

# 1.1.1.3 use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 7:10])
pred_Discont <- predict(Fit_Discont, test[, 7:10])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

# ROC Curve
plot(roc_Death_lm, legend = F, lwd = 1, lty = 1, col = 2, auc.main = F,
     main = "ASC+CEL ROC curve\nALB+HB+PSA+ALP", xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(roc_Discont_lm, add = T, lwd = 1, lty = 1, col = 3)
plot(roc_Death_bag, legend = F, lwd = 1, lty = 2, col = 2, auc.main = F, add = T)
plot(roc_Discont_bag, add = T, lwd = 1, lty = 2, col = 3)
plot(roc_Death_rf, legend = F, lwd = 1, lty = 3, col = 2, auc.main = F, add = T)
plot(roc_Discont_rf, add = T, lwd = 1, lty = 3, col = 3)
abline(a = 0, b = 1)
legend("topleft", 
       c(paste0("DEATH, ", "Linear Regression ", "AUC:", round(roc_Death_lm$auc, 4)),
         paste0("DISCONT, ", "Linear Regression ", "AUC:", round(roc_Discont_lm$auc, 4)),
         paste0("DEATH, ", "Bag-Regression Tree ", "AUC:", round(roc_Death_bag$auc, 4)),
         paste0("DISCONT, ", "Bag-Regression Tree ", "AUC:", round(roc_Discont_bag$auc, 4)),
         paste0("DEATH, ", "Random Forest ", "AUC:", round(roc_Death_rf$auc, 4)),
         paste0("DISCONT, ", "Random Forest ", "AUC:", round(roc_Discont_rf$auc, 4))),
       col = c(2, 3, 2, 3, 2, 3), lty = c(1, 1, 2, 2, 3, 3), cex = 0.7, bty = "n", y.intersp = 0.7)

# PR Curve
plot(pr_Death_lm, legend = F, lwd = 1, lty = 1, col = 2, auc.main = F,
     main = "ASC+CEL PR curve\nALB+HB+PSA+ALP", xlab = "False Positive Rate", ylab = "True Positive Rate")
plot(pr_Discont_lm, add = T, lwd = 1, lty = 1, col = 3)
plot(pr_Death_bag, legend = F, lwd = 1, lty = 2, col = 2, auc.main = F, add = T)
plot(pr_Discont_bag, add = T, lwd = 1, lty = 2, col = 3)
plot(pr_Death_rf, legend = F, lwd = 1, lty = 3, col = 2, auc.main = F, add = T)
plot(pr_Discont_rf, add = T, lwd = 1, lty = 3, col = 3)
legend(x = 0.5, y = 1, 
       c(paste0("DEATH, ", "Linear Regression ", "AUPRC:", round(pr_Death_lm$auc.integral, 4)),
         paste0("DISCONT, ", "Linear Regression ", "AUPRC:", round(pr_Discont_lm$auc.integral, 4)),
         paste0("DEATH, ", "Bag-Regression Tree ", "AUPRC:", round(pr_Death_bag$auc.integral, 4)),
         paste0("DISCONT, ", "Bag-Regression Tree ", "AUPRC:", round(pr_Discont_bag$auc.integral, 4)),
         paste0("DEATH, ", "Random Forest ", "AUPRC:", round(pr_Death_rf$auc.integral, 4)),
         paste0("DISCONT, ", "Random Forest ", "AUPRC:", round(pr_Discont_rf$auc.integral, 4))),
       col = c(2, 3, 2, 3, 2, 3), lty = c(1, 1, 2, 2, 3, 3), cex = 0.7, bty = "n", y.intersp = 0.7)

# ## 1.1.2 use all variable in halabi parameters
# ## 1.1.2.1 use linear regression model
# Fit_Death <- lm(formula = DEATH ~ ANALGESICS+LDH+BONE+ECOG_C+ALB+HB+PSA+ALP+LYMPH_NODES, data = train)
# Fit_Discont <- lm(formula = DISCONT ~ ANALGESICS+LDH+BONE+ECOG_C+ALB+HB+PSA+ALP+LYMPH_NODES, data = train)
# pred_Death <- predict(Fit_Death, test[, c(-1, -2)])
# pred_Discont <- predict(Fit_Discont, test[, c(-1, -2)])
# roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## 1.1.2.2 use bootstrap aggregated classification and regression tree
# set.seed(1357)
# Fit_Death <- bagging(formula = DEATH ~ ANALGESICS+LDH+BONE+ECOG_C+ALB+HB+PSA+ALP+LYMPH_NODES, data = train)
# Fit_Discont <- bagging(formula = DISCONT ~ ANALGESICS+LDH+BONE+ECOG_C+ALB+HB+PSA+ALP+LYMPH_NODES, data = train)
# pred_Death <- predict(Fit_Death, test[, c(-1, -2)])
# pred_Discont <- predict(Fit_Discont, test[, c(-1, -2)])
# roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# # 1.1.2.3 use random forest 
# set.seed(2468)
# Fit_Death <- randomForest(formula = DEATH ~ ANALGESICS+LDH+BONE+ECOG_C+ALB+HB+PSA+ALP+LYMPH_NODES, data = train)
# Fit_Discont <- randomForest(formula = DISCONT ~ ANALGESICS+LDH+BONE+ECOG_C+ALB+HB+PSA+ALP+LYMPH_NODES, data = train)
# pred_Death <- predict(Fit_Death, test[, c(-1, -2)])
# pred_Discont <- predict(Fit_Discont, test[, c(-1, -2)])
# roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# plot(roc_Death_lm, legend = F, lwd = 1, lty = 1, col = 2, auc.main = F,
#      main = "ROC curve\nAll Variables", xlab = "False Positive Rate", ylab = "True Positive Rate")
# plot(roc_Discont_lm, add = T, lwd = 1, lty = 1, col = 3)
# plot(roc_Death_bag, legend = F, lwd = 1, lty = 2, col = 2, auc.main = F, add = T)
# plot(roc_Discont_bag, add = T, lwd = 1, lty = 2, col = 3)
# plot(roc_Death_rf, legend = F, lwd = 1, lty = 3, col = 2, auc.main = F, add = T)
# plot(roc_Discont_rf, add = T, lwd = 1, lty = 3, col = 3)
# abline(a = 0, b = 1)
# legend("topleft", 
#        c(paste0("DEATH, ", "Linear Regression ", "AUC:", round(roc_Death_lm$auc, 4)),
#          paste0("DISCONT, ", "Linear Regression ", "AUC:", round(roc_Discont_lm$auc, 4)),
#          paste0("DEATH, ", "Bag-Regression Tree ", "AUC:", round(roc_Death_bag$auc, 4)),
#          paste0("DISCONT, ", "Bag-Regression Tree ", "AUC:", round(roc_Discont_bag$auc, 4)),
#          paste0("DEATH, ", "Random Forest ", "AUC:", round(roc_Death_rf$auc, 4)),
#          paste0("DISCONT, ", "Random Forest ", "AUC:", round(roc_Discont_rf$auc, 4))),
#        col = c(2, 3, 2, 3, 2, 3), lty = c(1, 1, 2, 2, 3, 3), cex = 0.65, bty = "n", y.intersp = 0.7)
# 
# 
# plot(pr_Death_lm, legend = F, lwd = 1, lty = 1, col = 2, auc.main = F,
#      main = "PR curve\nAll Variables", xlab = "False Positive Rate", ylab = "True Positive Rate")
# plot(pr_Discont_lm, add = T, lwd = 1, lty = 1, col = 3)
# plot(pr_Death_bag, legend = F, lwd = 1, lty = 2, col = 2, auc.main = F, add = T)
# plot(pr_Discont_bag, add = T, lwd = 1, lty = 2, col = 3)
# plot(pr_Death_rf, legend = F, lwd = 1, lty = 3, col = 2, auc.main = F, add = T)
# plot(pr_Discont_rf, add = T, lwd = 1, lty = 3, col = 3)
# legend(x = 0.5, y = 1, 
#        c(paste0("DEATH, ", "Linear Regression ", "AUPRC:", round(pr_Death_lm$auc.integral, 4)),
#          paste0("DISCONT, ", "Linear Regression ", "AUPRC:", round(pr_Discont_lm$auc.integral, 4)),
#          paste0("DEATH, ", "Bag-Regression Tree ", "AUPRC:", round(pr_Death_bag$auc.integral, 4)),
#          paste0("DISCONT, ", "Bag-Regression Tree ", "AUPRC:", round(pr_Discont_bag$auc.integral, 4)),
#          paste0("DEATH, ", "Random Forest ", "AUPRC:", round(pr_Death_rf$auc.integral, 4)),
#          paste0("DISCONT, ", "Random Forest ", "AUPRC:", round(pr_Discont_rf$auc.integral, 4))),
#        col = c(2, 3, 2, 3, 2, 3), lty = c(1, 1, 2, 2, 3, 3), cex = 0.7, bty = "n", y.intersp = 0.7)
# 
# # Important variables same as 1.1.1, so I don't try another