library(MASS)
library(e1071)
setwd("~/Desktop/GuanLab/")

# ------------------------------------------ Plot 2 ----------------------------------------------

# CV InCohort ------------------

# 25 times 5-fold cv (waiting)
kFoldCV <- function(data, fold = 5, seed) {
  require(ipred)
  require(PRROC)
  require(randomForest)
  
  # generate groups within cohort
  set.seed(seed)
  data <- data[sample(nrow(data)), ]
  folds <- cut(seq(1, nrow(data)), breaks = fold, labels = F)
  
  AUC_DEATH_lm <- vector(length = 5)
  AUC_DISCONT_lm <- vector(length = 5)
  AUPRC_DEATH_lm <- vector(length = 5)
  AUPRC_DISCONT_lm <- vector(length = 5)
  
  AUC_DEATH_bag <- vector(length = 5)
  AUC_DISCONT_bag <- vector(length = 5)
  AUPRC_DEATH_bag <- vector(length = 5)
  AUPRC_DISCONT_bag <- vector(length = 5)
  
  AUC_DEATH_rf <- vector(length = 5)
  AUC_DISCONT_rf <- vector(length = 5)
  AUPRC_DEATH_rf <- vector(length = 5)
  AUPRC_DISCONT_rf <- vector(length = 5)
  
  # AUC_DEATH_cox <- vector(length = 5)
  # AUC_DISCONT_cox <- vector(length = 5)
  # AUPRC_DEATH_cox <- vector(length = 5)
  # AUPRC_DISCONT_cox <- vector(length = 5)
  
  # AUC_DEATH_lr <- vector(length = 5)
  # AUC_DISCONT_lr <- vector(length = 5)
  # AUPRC_DEATH_lr <- vector(length = 5)
  # AUPRC_DISCONT_lr <- vector(length = 5)
  
  for (i in 1:fold) {
    testIndexes <- which(folds == i, arr.ind = T)
    test <- data[testIndexes, ]
    train <- data[-testIndexes, ]
    
    # base learner train and test
    train$DISCONT <- as.numeric(as.character(train$DISCONT))
    train <- na.omit(train)
    test <- test[!grepl("\\.", test$DISCONT), ]
    test$DISCONT <- as.numeric(as.character(test$DISCONT))
    
    Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
    Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
    pred_Death <- predict(Fit_Death, test[, 8:11])
    pred_Discont <- predict(Fit_Discont, test[, 8:11])
    roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    
    AUC_DEATH_lm[i] <- roc_Death_lm$auc
    AUC_DISCONT_lm[i] <- roc_Discont_lm$auc
    AUPRC_DEATH_lm[i] <- pr_Death_lm$auc.integral
    AUPRC_DISCONT_lm[i] <- pr_Discont_lm$auc.integral
    
    set.seed(1234)
    Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
    Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
    pred_Death <- predict(Fit_Death, test[, 8:11])
    pred_Discont <- predict(Fit_Discont, test[, 8:11])
    roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    
    AUC_DEATH_bag[i] <- roc_Death_bag$auc
    AUC_DISCONT_bag[i] <- roc_Discont_bag$auc
    AUPRC_DEATH_bag[i] <- pr_Death_bag$auc.integral
    AUPRC_DISCONT_bag[i] <- pr_Discont_bag$auc.integral
    
    set.seed(2345)
    Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
    Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
    pred_Death <- predict(Fit_Death, test[, 8:11])
    pred_Discont <- predict(Fit_Discont, test[, 8:11])
    roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    
    AUC_DEATH_rf[i] <- roc_Death_rf$auc
    AUC_DISCONT_rf[i] <- roc_Discont_rf$auc
    AUPRC_DEATH_rf[i] <- pr_Death_rf$auc.integral
    AUPRC_DISCONT_rf[i] <- pr_Discont_rf$auc.integral
    
    # Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
    # Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
    # pred_Death <- predict(Fit_Death, test[, 8:11])
    # pred_Discont <- predict(Fit_Discont, test[, 8:11])
    # roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    # roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    # pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    # pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    # 
    # AUC_DEATH_cox[i] <- roc_Death_cox$auc
    # AUC_DISCONT_cox[i] <- roc_Discont_cox$auc
    # AUPRC_DEATH_cox[i] <- pr_Death_cox$auc.integral
    # AUPRC_DISCONT_cox[i] <- pr_Discont_cox$auc.integral
    
    # Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
    # Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
    # pred_Death <- predict(Fit_Death, test[, 8:11])
    # pred_Discont <- predict(Fit_Discont, test[, 8:11])
    # roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    # roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    # pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    # pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    
    # AUC_DEATH_lr[i] <- roc_Death_lr$auc
    # AUC_DISCONT_lr[i] <- roc_Discont_lr$auc
    # AUPRC_DEATH_lr[i] <- pr_Death_lr$auc.integral
    # AUPRC_DISCONT_lr[i] <- pr_Discont_lr$auc.integral
    
  }
  return(data.frame(AUC_DEATH_lm, AUC_DISCONT_lm, AUC_DEATH_bag, AUC_DISCONT_bag,
                    AUC_DEATH_rf, AUC_DISCONT_rf, AUPRC_DEATH_lm, AUPRC_DISCONT_lm,
                    AUPRC_DEATH_bag, AUPRC_DISCONT_bag, AUPRC_DEATH_rf, AUPRC_DISCONT_rf))
  # AUC_DEATH_cox, AUC_DISCONT_cox, AUC_DEATH_lr, AUC_DISCONT_lr,
  # AUPRC_DEATH_cox, AUPRC_DISCONT_cox, AUPRC_DEATH_lr, AUPRC_DISCONT_lr))
}
reshapeData <- function(data) {
  require(reshape2)
  data <- melt(data, value.name = "Value")
  data$CurveType <- sapply(data[, "variable"], function(x) ifelse(grepl("^AUC", x), "AUC", "AUPRC"))
  data$GoldStandard <- sapply(data[, "variable"], function(x) ifelse(grepl("_DEATH_", x), "Death", "Discontinuation"))
  data$Model <- sapply(data[, "variable"],
                       function(x) {
                         if (grepl("lm$", x)) return("Linear Regression")
                         if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                         if (grepl("rf$", x)) return("Random Forest")
                         # if (grepl("cox$", x)) return("Cox")
                         # if (grepl("lr$", x)) return("Logistic Regression")
                       })
  return(data)
}

CELinCohortCV <- data.frame()
VENinCohortCV <- data.frame()
ASCinCohortCV <- data.frame()

for (i in 1:25) {
  CELinCohortCV <- rbind(CELinCohortCV, reshapeData(kFoldCV(habini_CEL_death_3, 5, i)))
  VENinCohortCV <- rbind(VENinCohortCV, reshapeData(kFoldCV(habini_VEN_death_3, 5, i)))
  ASCinCohortCV <- rbind(ASCinCohortCV, reshapeData(kFoldCV(habini_ASC_death_3, 5, i)))
}

# 25 times 5-fold cv -- AUC
p <- ggplot(data = subset(CELinCohortCV, CurveType == "AUC"))
p <- p + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0.1, 0.9)) + geom_hline(yintercept = 0.5, linetype = "dashed")
p1 <- ggplot(data = subset(VENinCohortCV, CurveType == "AUC"))
p1 <- p1 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0.1, 0.9)) + geom_hline(yintercept = 0.5, linetype = "dashed")
p2 <- ggplot(data = subset(ASCinCohortCV, CurveType == "AUC"))
p2 <- p2 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0.1, 0.9)) + geom_hline(yintercept = 0.5, linetype = "dashed")

grid.arrange(p, p1, p2, nrow = 2, ncol = 2)

# 25 times 5-fold cv -- AUPRC
p <- ggplot(data = subset(CELinCohortCV, CurveType == "AUPRC"))
p <- p + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  scale_y_continuous(limits = c(0, 0.4))
p1 <- ggplot(data = subset(VENinCohortCV, CurveType == "AUPRC"))
p1 <- p1 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  scale_y_continuous(limits = c(0, 0.4))
p2 <- ggplot(data = subset(ASCinCohortCV, CurveType == "AUPRC"))
p2 <- p2 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  scale_y_continuous(limits = c(0, 0.4))

grid.arrange(p, p1, p2, nrow = 2, ncol = 2)


# ASC + CEL ------------------
train <- rbind(habini_ASC_death_3, habini_CEL_death_3)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_VEN_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use more classification algorithm
## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p <- ggplot(subset(data_for_plot, Curve == "AUC"))
p <- p + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- ASC + CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp <- pp + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- ASC + CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) + 
  scale_y_continuous(limits = c(0, 0.3))

# ASC + VEN ------------------ 
train <- rbind(habini_ASC_death_3, habini_VEN_death_3)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_CEL_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use more classification algorithm
## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p1 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p1 <- p1 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- ASC + VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp1 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp1 <- pp1 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- ASC + VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) + 
  scale_y_continuous(limits = c(0, 0.3))


# CEL + VEN ------------------ 
train <- rbind(habini_VEN_death_3, habini_CEL_death_3)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_ASC_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use more classification algorithm
## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                if (grepl("cox$", x)) return("Cox")
                                if (grepl("lr$", x)) return("Logistic Regression")
                              })
p2 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p2 <- p2 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- VEN + CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp2 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp2 <- pp2 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- VEN + CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) + 
  scale_y_continuous(limits = c(0, 0.3))


# ---------------
grid.arrange(p, p1, p2, nrow = 1, ncol = 3)
tiff("AUPRC_CrossCohort.tiff", width = 1618, height = 531)
grid.arrange(pp, pp1, pp2, nrow = 1, ncol = 3)
# ---------------

# In order to eliminate the bias between cohorts
# 1. try to train model in single cohort and test in the other two
# ASC --> CEL ------------------
train <- habini_ASC_death_3
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_CEL_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use more classification algorithm
## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p <- ggplot(subset(data_for_plot, Curve == "AUC"))
p <- p + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "ASC --> CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp <- pp + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "ASC --> CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0, 0.3))



# ASC --> VEN -----------------
train <- habini_ASC_death_3
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_VEN_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use more classification algorithm
## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p1 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p1 <- p1 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "ASC --> VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp1 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp1 <- pp1 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "ASC --> VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0, 0.3))

# VEN --> ASC ------------------ 
train <- habini_VEN_death_3
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_ASC_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

# ## use more classification algorithm
# ## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p2 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p2 <- p2 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "VEN --> ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp2 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp2 <- pp2 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "VEN --> ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0, 0.3))


# VEN --> CEL ------------------------
train <- habini_VEN_death_3
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_CEL_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

# ## use more classification algorithm
# ## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p3 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p3 <- p3 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "VEN --> CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp3 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp3 <- pp3 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "VEN --> CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0, 0.3))

# CEL --> ASC ----------------------
train <- habini_CEL_death_3
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_ASC_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

# ## use more classification algorithm
# ## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p4 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p4 <- p4 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "CEL --> ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp4 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp4 <- pp4 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "CEL --> ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0, 0.3))

# CEL --> VEN ------------------ 
train <- habini_CEL_death_3
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
test <- habini_VEN_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

# ## use more classification algorithm
# ## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                                     # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                                     # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
                                      # roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
                                      # roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p5 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p5 <- p5 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "CEL --> VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp5 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp5 <- pp5 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard), alpha = 0.8) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "CEL --> VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0, 0.3))







# ----------------
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
grid_arrange_share_legend(p, p1, p4, p5, p2, p3, nrow = 3, ncol = 2, position = "right")
grid_arrange_share_legend(pp, pp1, pp4, pp5, pp2, pp3, nrow = 3, ncol = 2, position = "right")

# ----------------
# ASC + CEL --> VEN ------------------
train <- rbind(habini_ASC_death_3, habini_CEL_death_3)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
train[, 8:11] <- lapply(train[, 8:11], scale)
test <- habini_VEN_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))
test[, 8:11] <- lapply(test[, 8:11], scale)

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use more classification algorithm
## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                            # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                            # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
# roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
# roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p <- ggplot(subset(data_for_plot, Curve == "AUC"))
p <- p + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- ASC + CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp <- pp + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- ASC + CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) + 
  scale_y_continuous(limits = c(0, 0.3))

# ASC + VEN --> CEL ------------------ 
train <- rbind(habini_ASC_death_3, habini_VEN_death_3)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
train[, 8:11] <- lapply(train[, 8:11], scale)
test <- habini_CEL_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))
test[, 8:11] <- lapply(test[, 8:11], scale)

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use more classification algorithm
## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# 
# ## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                            # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                            # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
# roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
# roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                # if (grepl("cox$", x)) return("Cox")
                                # if (grepl("lr$", x)) return("Logistic Regression")
                              })
p1 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p1 <- p1 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- ASC + VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp1 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp1 <- pp1 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- ASC + VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) + 
  scale_y_continuous(limits = c(0, 0.3))


# CEL + VEN --> ASC ------------------ 
train <- rbind(habini_VEN_death_3, habini_CEL_death_3)
train$DISCONT <- as.numeric(as.character(train$DISCONT))
train <- na.omit(train)
train[, 8:11] <- lapply(train[, 8:11], function(x) c(scale(x)))
test <- habini_ASC_death_3
test <- test[!grepl("\\.", test$DISCONT), ]
test$DISCONT <- as.numeric(as.character(test$DISCONT))
test[, 8] <- 0
test[, 9:11] <- lapply(test[, 9:11], function(x) c(scale(x)))

Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_lm <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_lm <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_lm <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_lm <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use bootstrap aggregated classification and regression tree
set.seed(1234)
Fit_Death <- bagging(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
Fit_Discont <- bagging(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_bag <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_bag <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_bag <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_bag <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use random forest 
set.seed(2345)
Fit_Death <- randomForest(formula = DEATH ~ ALB+HB+PSA+ALP, data = train )
Fit_Discont <- randomForest(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train )
pred_Death <- predict(Fit_Death, test[, 8:11])
pred_Discont <- predict(Fit_Discont, test[, 8:11])
roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## use more classification algorithm
## Cox model
# Fit_Death <- coxph(Surv(LKADT_P, DEATH) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# Fit_Discont <- coxph(Surv(LKADT_P, DISCONT) ~ ALB + HB + PSA + ALP, ties = "breslow", data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_cox <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_cox <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_cox <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_cox <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

## Logistic Regression
# Fit_Death <- glm(as.factor(DEATH) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# Fit_Discont <- glm(as.factor(DISCONT) ~ ALB + HB + PSA + ALP, family = binomial(link = "logit"), data = train)
# pred_Death <- predict(Fit_Death, test[, 8:11])
# pred_Discont <- predict(Fit_Discont, test[, 8:11])
# roc_Death_lr <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# roc_Discont_lr <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
# pr_Death_lr <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
# pr_Discont_lr <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)

data_for_plot <- data.frame(Type = c("roc_Death_lm", "roc_Discont_lm", "pr_Death_lm", "pr_Discont_lm", 
                                     "roc_Death_bag", "roc_Discont_bag", "pr_Death_bag", "pr_Discont_bag",
                                     "roc_Death_rf", "roc_Discont_rf", "pr_Death_rf", "pr_Discont_rf"),
                            # "roc_Death_cox", "roc_Discont_cox", "pr_Death_cox", "pr_Discont_cox",
                            # "roc_Death_lr", "roc_Discont_lr", "pr_Death_lr", "pr_Discont_lr"),
                            Value = c(roc_Death_lm$auc, roc_Discont_lm$auc, pr_Death_lm$auc.integral, pr_Discont_lm$auc.integral,
                                      roc_Death_bag$auc, roc_Discont_bag$auc, pr_Death_bag$auc.integral, pr_Discont_bag$auc.integral,
                                      roc_Death_rf$auc, roc_Discont_rf$auc, pr_Death_rf$auc.integral, pr_Discont_rf$auc.integral))
# roc_Death_cox$auc, roc_Discont_cox$auc, pr_Death_cox$auc.integral, pr_Discont_cox$auc.integral,
# roc_Death_lr$auc, roc_Discont_lr$auc, pr_Death_lr$auc.integral, pr_Discont_lr$auc.integral))
data_for_plot$Curve <- ifelse(grepl("^roc", data_for_plot$Type), "AUC", "AUPRC")
data_for_plot$GoldStandard <- ifelse(grepl("_Death_", data_for_plot$Type), "Death", "Discontinuation")
data_for_plot$Model <- sapply(data_for_plot$Type,
                              function(x) {
                                if (grepl("lm$", x)) return("Linear Regression")
                                if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                                if (grepl("rf$", x)) return("Random Forest")
                                if (grepl("cox$", x)) return("Cox")
                                if (grepl("lr$", x)) return("Logistic Regression")
                              })
p2 <- ggplot(subset(data_for_plot, Curve == "AUC"))
p2 <- p2 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- VEN + CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "Linear", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + scale_y_continuous(limits = c(0, 0.7))
pp2 <- ggplot(subset(data_for_plot, Curve == "AUPRC"))
pp2 <- pp2 + geom_point(aes(x = Model, y = Value, fill = GoldStandard, color = GoldStandard)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- VEN + CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF", "Cox" = "Cox", "Logistic Regression" = "Logistic")) + 
  scale_y_continuous(limits = c(0, 0.3))




grid.arrange(p, p1, p2, nrow = 1, ncol = 3)
grid.arrange(pp, pp1, pp2, nrow = 1, ncol = 3)
