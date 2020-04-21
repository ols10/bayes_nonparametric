### Estimation based on real data
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(caret)
library(xgboost)
library(earth)
library(glm2)
library(glmnet)
library(SuperLearner)
library(mcmc)
library(ggthemes)
library(origami)
library(fastDummies)
library(future)
library(tidyverse)
library(mvtnorm)

options(mc.cores = 6)

set.seed(539)


## Simulate a dataset similar to CAESAR.
n <- 1700
X <- rmvnorm(n, mean = rep(0, 200))
colnames(X) <- paste0('X', 1:200)
L <- rbinom(n, 1, prob = plogis(1 + rowMeans(X)))
Y <- rbinom(n, 1, prob = plogis(1 - rowMeans(X)))
Y[L == 0] <- NA
dat <- data.frame(X, L, Y)

## Set up cross-validation function for each fold
cv_fun <- function(fold, data, lib){

    train_data <- training(data)
    valid_data <- data %>% validation() %>% select(-L, -Y) %>% data.matrix()

    y_train <- train_data$L
    x_train <- train_data %>% select(-L, -Y) %>% data.matrix()
    g_fit <- SuperLearner(Y = y_train, X = x_train, family = binomial(),
                          SL.library = lib, verbose = TRUE)

    y2_train <- train_data %>% filter(!is.na(Y)) %>% pull(Y)
    x2_train <- train_data %>% filter(!is.na(Y)) %>% select(-L, -Y) %>% data.matrix()
    m_fit <- SuperLearner(Y = y2_train, X = x2_train, family = binomial(),
                          SL.library = lib)

    g_pred <- as.vector(predict(g_fit, valid_data)$pred)
    m_pred <- as.vector(predict(m_fit, valid_data)$pred)

    return(list(pred = data.frame(g = g_pred, m = m_pred),
                fold = fold$validation_set))

}

## Set up learners
SL.caretRF <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'rf',  tuneLength = 20,
             trControl =  caret::trainControl(method = "cv", number = 5, search = 'random',
                                              verboseIter = TRUE), ...)
}
SL.caretXGB <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'xgbTree', tuneLength = 30,
             trControl =  caret::trainControl(method = "cv", number = 5, search = 'random',
                                              verboseIter = TRUE), ...)
}
screen.ttest50 <- function(Y, X, family, obsWeights, id, rank = 50, ...)
    screen.ttest(Y, X, family, obsWeights, id, rank = rank)

lib1  <- list(
    c('SL.glm', 'screen.ttest50'),
    c('SL.earth', 'screen.ttest50'),
    'SL.caretRF',
    'SL.caretXGB',
    'SL.glmnet')


## Cross-fit the initital estimates
folds <- make_folds(dat, fold_fun = function(n)folds_vfold(n, V = 10))

plan(multicore)
cv_results1 <- cross_validate(cv_fun = cv_fun, folds = folds, data = dat,
                              lib = lib1)
cv_res1 <- cv_results1$pred[order(cv_results1$fold), ]

## Compute prior distribution
beta_reparam <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

beta_res <- beta_reparam(0.3, 0.16)
pi <- function(theta) dbeta(theta, beta_res$alpha, beta_res$beta)

## Plot the prior distribution
pdf('prior.pdf')
curve(pi, xlab = '', ylab = '', )
dev.off()

## Compute the posterior distribution
dat$Y[is.na(dat$Y)] <- -9999
pp1 <- posterior(Y = dat$Y, L = dat$L, X, cv_res1$m, cv_res1$g, pi, 1e6)
theta_post1 <- pp1$post$theta


tmle_post1 <- pp1$tmle
pp <- pp1
tmle_post <- tmle_post1
theta_post <- theta_post1
est_val <-  mean(pp$post$theta)
sd_val <-  sd(pp$post$theta)
lower <- quantile(pp$post$theta, 0.025)
upper <- quantile(pp$post$theta, 0.975)
lower_99 <- quantile(pp$post$theta, 0.005)
upper_99 <- quantile(pp$post$theta, 0.995)

## Plot the posterior distribution
title_str <- paste("Posterior distribution for ", expression(theta), sep = " ")
subtitle_str <- str_c(95, "% ", "Credible interval: [",
                     round(lower, 3), ", ", round(upper, 3), "]",
                     sep = "")

subtitle_str <- str_c("with ", 95, "% ", "and ", 99, "% ", "credible intervals",sep = "")
data <- with(density(theta_post), data.frame(x, y))
breaks <- c(round(lower_99, 3), round(lower, 3), round(est_val, 3),
           round(upper, 3), round(upper_99, 3))

est_den <- ggplot(data, aes(x, y)) +
    geom_line(size = 1.5) +
    geom_area(data = subset(data, x > lower & x < upper), fill = 'grey', alpha = 0.7) +
    geom_area(data = subset(data, x > quantile(pp$post$theta, 0.005) &
                                  x < quantile(pp$post$theta, 0.995)), fill = 'grey', alpha = 0.5) +
    geom_vline(xintercept = est_val, color = "black", size = 1.5) +
    xlab("") + ylab("") +
    xlim(c(0, 1)) +
    scale_x_continuous(limits = c(min(breaks) * 0.9, max(breaks) * 1.1),
                       breaks = breaks, labels = as.character(breaks)) +
    theme_igray() +
    stat_function(fun = dnorm, args = list(mean = tmle_post1$estimate,
                                           sd = sd(tmle_post$infun) / sqrt(length(tmle_post$infun))),
                  linetype = 'dashed', size = 1.5) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.6),
          line = element_line(size = 4),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 30))

pdf('posterior.pdf')
est_den
dev.off()

## Plot the ROC curves

library(ROCR)

dfg <- data.frame(predictions = cv_res1$g, labels = dat$L)
predg <- prediction(dfg$predictions, dfg$labels)
perfg <- performance(predg, "tpr", "fpr")
aucg <- performance(predg, measure = "auc")

dfm <- data.frame(predictions = cv_res1$m[dat$L == 1], labels = dat$Y[dat$L == 1])
predm <- prediction(dfm$predictions, dfm$labels)
perfm <- performance(predm, "tpr", "fpr")
aucm <- performance(predm, measure = "auc")

pdf('rocs.pdf', width = 10, height = 5)
par(mfrow = c(1, 2), cex = 1.6)
plot(perfg, colorize = FALSE, main = 'Prob. source is known (g)',
     lwd = 3)
text(x = 0.6, y = 0.2,
     labels = paste0('AUC = ', round(aucg@y.values[[1]], 3)))
plot(perfm, colorize = FALSE, main = 'Prob. of cardiac sources (m)',
     lwd = 3, ylab = '')
text(x = 0.6, y = 0.2,
     labels = paste0('AUC = ', round(aucm@y.values[[1]], 3)))
dev.off()

