library(ggplot2)
library(gridExtra)
kFoldCV <- function(data, fold = 5) {
  require(ipred)
  require(PRROC)
  require(randomForest)
  
  set.seed(1234)
  data <- data[sample(nrow(data)), ]
  folds <- cut(seq(1, nrow(data)), breaks = fold, labels = F)
  
  for (i in 1:fold) {
    testIndexes <- which(folds == i, arr.ind = T)
    test <- data[testIndexes, ]
    train <- data[-testIndexes, ]
    
    train$DISCONT <- as.numeric(as.character(train$DISCONT))
    train <- na.omit(train)
    test <- test[!grepl("\\.", test$DISCONT), ]
    test$DISCONT <- as.numeric(as.character(test$DISCONT))
    
    Fit_Death <- lm(formula = DEATH ~ ALB+HB+PSA+ALP, data = train)
    Fit_Discont <- lm(formula = DISCONT ~ ALB+HB+PSA+ALP, data = train)
    pred_Death <- predict(Fit_Death, test[, 7:10])
    pred_Discont <- predict(Fit_Discont, test[, 7:10])
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
    pred_Death <- predict(Fit_Death, test[, 7:10])
    pred_Discont <- predict(Fit_Discont, test[, 7:10])
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
    pred_Death <- predict(Fit_Death, test[, 7:10])
    pred_Discont <- predict(Fit_Discont, test[, 7:10])
    roc_Death_rf <- roc.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    roc_Discont_rf <- roc.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    pr_Death_rf <- pr.curve(scores.class0 = pred_Death, weights.class0 = test$DISCONT, curve = T)
    pr_Discont_rf <- pr.curve(scores.class0 = pred_Discont, weights.class0 = test$DISCONT, curve = T)
    
    AUC_DEATH_rf[i] <- roc_Death_rf$auc
    AUC_DISCONT_rf[i] <- roc_Discont_rf$auc
    AUPRC_DEATH_rf[i] <- pr_Death_rf$auc.integral
    AUPRC_DISCONT_rf[i] <- pr_Discont_rf$auc.integral
  }
  return(data.frame(AUC_DEATH_lm, AUC_DISCONT_lm, AUC_DEATH_bag, AUC_DISCONT_bag,
                    AUC_DEATH_rf, AUC_DISCONT_rf, AUPRC_DEATH_lm, AUPRC_DISCONT_lm,
                    AUPRC_DEATH_bag, AUPRC_DISCONT_bag, AUPRC_DEATH_rf, AUPRC_DISCONT_rf))
}

CELinCohortCV <- kFoldCV(habini_CEL_death, 5)
VENinCohortCV <- kFoldCV(habini_VEN_death, 5)
ASCinCohortCV <- kFoldCV(habini_ASC_death, 5)

CELinCohortCV_GR <- kFoldCV(habini_CEL_death_GR, 5)
VENinCohortCV_GR <- kFoldCV(habini_VEN_death_GR, 5)
ASCinCohortCV_GR <- kFoldCV(habini_ASC_death_GR, 5)

reshapeData <- function(data, GR = TRUE) {
  require(reshape2)
  data <- melt(data, value.name = "Value")
  data$CurveType <- sapply(data[, "variable"], function(x) ifelse(grepl("^AUC", x), "AUC", "AUPRC"))
  data$GoldStandard <- sapply(data[, "variable"], function(x) ifelse(grepl("_DEATH_", x), "Death", "Discontinuation"))
  data$Model <- sapply(data[, "variable"],
                  function(x) {
                    if (grepl("lm$", x)) return("Linear Regression")
                    if (grepl("bag$", x)) return("Bootstrap-Aggregated CART")
                    if (grepl("rf$", x)) return("Random Forest")
                  })
  data$GuanRank <- ifelse(GR, "Used", "Not Used")
  return(data)
}

CELinCohortCV <- reshapeData(CELinCohortCV, GR = F)
VENinCohortCV <- reshapeData(VENinCohortCV, GR = F)
ASCinCohortCV <- reshapeData(ASCinCohortCV, GR = F)

CELinCohortCV_GR <- reshapeData(CELinCohortCV_GR)
VENinCohortCV_GR <- reshapeData(VENinCohortCV_GR)
ASCinCohortCV_GR <- reshapeData(ASCinCohortCV_GR)

CELinCohortCV <- rbind(CELinCohortCV, CELinCohortCV_GR)
VENinCohortCV <- rbind(VENinCohortCV, VENinCohortCV_GR)
ASCinCohortCV <- rbind(ASCinCohortCV, ASCinCohortCV_GR)

p <- ggplot(data = subset(CELinCohortCV, CurveType == "AUC"))
p <- p + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0.4, 0.8))
p1 <- ggplot(data = subset(VENinCohortCV, CurveType == "AUC"))
p1 <- p1 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0.4, 0.8))
p2 <- ggplot(data = subset(ASCinCohortCV, CurveType == "AUC"))
p2 <- p2 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0.4, 0.8))
p3 <- ggplot(data = subset(rbind(CELinCohortCV, VENinCohortCV, ASCinCohortCV), CurveType == "AUC"))
p3 <- p3 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUC of different model -- Overall") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0.4, 0.8))
gridExtra::grid.arrange(p, p1, p2, p3, nrow = 2, ncol = 2)

p <- ggplot(data = subset(CELinCohortCV, CurveType == "AUPRC"))
p <- p + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0, 0.4))
p1 <- ggplot(data = subset(VENinCohortCV, CurveType == "AUPRC"))
p1 <- p1 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0, 0.4))
p2 <- ggplot(data = subset(ASCinCohortCV, CurveType == "AUPRC"))
p2 <- p2 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0, 0.4))
p3 <- ggplot(data = subset(rbind(CELinCohortCV, VENinCohortCV, ASCinCohortCV), CurveType == "AUPRC"))
p3 <- p3 + geom_boxplot(aes(x = Model, y = Value, fill = GoldStandard), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + labs(y = "Area Under Curve", title = "AUPRC of different model -- Overall") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0, 0.4))
gridExtra::grid.arrange(p, p1, p2, p3, nrow = 2, ncol = 2)

p <- ggplot(data = subset(CELinCohortCV, CurveType == "AUC" & GoldStandard == "Death"))
p <- p + geom_boxplot(aes(x = Model, y = Value, color = GuanRank), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 10)) + labs(y = "Area Under Curve", title = "AUC of GuanRank and different model -- CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0.4, 0.8))
p1 <- ggplot(data = subset(VENinCohortCV, CurveType == "AUC" & GoldStandard == "Death"))
p1 <- p1 + geom_boxplot(aes(x = Model, y = Value, color = GuanRank), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 10)) + labs(y = "Area Under Curve", title = "AUC of GuanRank and different model -- VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0.4, 0.8))
p2 <- ggplot(data = subset(ASCinCohortCV, CurveType == "AUC" & GoldStandard == "Death"))
p2 <- p2 + geom_boxplot(aes(x = Model, y = Value, color = GuanRank), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 10)) + labs(y = "Area Under Curve", title = "AUC of GuanRank and different model -- ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0.4, 0.8))
p3 <- ggplot(data = subset(rbind(CELinCohortCV, VENinCohortCV, ASCinCohortCV), CurveType == "AUC" & GoldStandard == "Death"))
p3 <- p3 + geom_boxplot(aes(x = Model, y = Value, color = GuanRank), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 10)) + labs(y = "Area Under Curve", title = "AUC of GuanRank and different model -- Overall") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0.4, 0.8))
gridExtra::grid.arrange(p, p1, p2, p3, nrow = 2, ncol = 2)

p <- ggplot(data = subset(CELinCohortCV, CurveType == "AUPRC" & GoldStandard == "Death"))
p <- p + geom_boxplot(aes(x = Model, y = Value, color = GuanRank), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 10)) + labs(y = "Area Under Curve", title = "AUPRC of GuanRank and different model -- CEL") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0, 0.4))
p1 <- ggplot(data = subset(VENinCohortCV, CurveType == "AUPRC" & GoldStandard == "Death"))
p1 <- p1 + geom_boxplot(aes(x = Model, y = Value, color = GuanRank), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 10)) + labs(y = "Area Under Curve", title = "AUPRC of GuanRank and different model -- VEN") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0, 0.4))
p2 <- ggplot(data = subset(ASCinCohortCV, CurveType == "AUPRC" & GoldStandard == "Death"))
p2 <- p2 + geom_boxplot(aes(x = Model, y = Value, color = GuanRank), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 10)) + labs(y = "Area Under Curve", title = "AUPRC of GuanRank and different model -- ASC") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) +
  scale_y_continuous(limits = c(0, 0.4))
p3 <- ggplot(data = subset(rbind(CELinCohortCV, VENinCohortCV, ASCinCohortCV), CurveType == "AUPRC" & GoldStandard == "Death"))
p3 <- p3 + geom_boxplot(aes(x = Model, y = Value, color = GuanRank), alpha = 0.5) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(size = 10)) + labs(y = "Area Under Curve", title = "AUPRC of GuanRank and different model -- Overall") + 
  scale_x_discrete("Model", labels = c("Bootstrap-Aggregated CART" = "BAG-CART", "Linear Regression" = "LR", "Random Forest" = "RF")) + 
  scale_y_continuous(limits = c(0, 0.4))
gridExtra::grid.arrange(p, p1, p2, p3, nrow = 2, ncol = 2)


