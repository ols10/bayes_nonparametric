rm(list = ls())

source("functions_mod.r")
source("utils.r")

## if (!require(mcmc)) install.packages('mcmc')
library(mcmc)
library(MCMCpack)
library(dplyr)
set.seed(539)


N <- c(400, 900, 1600, 2500, 4900) # sample size
p_dim <- c(4, 100)
reps <- 1000
seeds <- sample(1e7, reps)
g_val <- c("correct", "incorrect")
m_val <- c("correct", "incorrect")
prior_dist <- c("correct_small", "correct_large", "incorrect_small", "incorrect_large")

tasks <- expand.grid(g_val = g_val, m_val = m_val,
                    prior = prior_dist, seed = seeds, N = N, p_dim = p_dim,
                    stringsAsFactors = FALSE)

## true_theta <- t(sapply(p_dim, function(p_dim)true(1e6, p = p_dim)))
## save(true_theta, file = 'true.rds')
load('true.rds')


## function returns output for one set of parameters
sim <- function(i) {

    set.seed(tasks[i, 'seed'])
    n <- tasks[i, 'N']
    g_val <- tasks[i, 'g_val']
    m_val <- tasks[i, 'm_val']
    prior <- tasks[i, 'prior']
    p_dim <- tasks[i, 'p_dim']

    p_miss <- round(p_dim / 2, 0)

    data <- gendata(n, p = p_dim)

    L <- data$L
    X <- data.frame(data$X)
    Y <- data$Y

    if (g_val == "correct"){
        g <- predict(glm(L ~ ., family = binomial, data = X), type = 'response')
    } else {
        g <- predict(glm(L ~ ., family = binomial, data = X[, 1:p_miss]), type = 'response')
    }

    if (m_val == "correct") {
        m <- predict(glm(Y[L == 1] ~ ., family = binomial, data = X[L == 1, ]),
                    newdata = X, type = 'response')
    } else {
        m <- predict(glm(Y[L == 1] ~ ., family = binomial, data = X[L == 1, 1:p_miss]),
                    newdata = X[, 1:p_miss], type = 'response')
    }

    theta_init <- true_theta[true_theta[, 'p'] == p_dim, 'truth']

    if (prior == "correct_small") {
        beta_res <- beta_reparam(theta_init, 0.1*(1 - theta_init)*theta_init)
        pi <- function(theta) dbeta(theta, beta_res$alpha, beta_res$beta)
    } else if (prior == "correct_large") {
        beta_res <- beta_reparam(theta_init, 0.9*(1 - theta_init)*theta_init)
        pi <- function(theta) dbeta(theta, beta_res$alpha, beta_res$beta)
    } else if (prior == "incorrect_small"){
        beta_res <- beta_reparam(1 - theta_init, 0.1*(1 - theta_init)*theta_init)
        pi <- function(theta) dbeta(theta, beta_res$alpha, beta_res$beta)
    } else {
        beta_res <- beta_reparam(1 - theta_init, 0.9*(1 - theta_init)*theta_init)
        pi <- function(theta) dbeta(theta, beta_res$alpha, beta_res$beta)
    }

    pp <- posterior(data$Y, data$L, m, g, pi, 1e4)

    if(p_dim == 4 && prior == "correct_large") {

        if (g_val == "correct"){
            g_fit <- MCMClogit(L ~ ., data = X, b0 = 0, B0 = 0.01)
            g_mat <- as.matrix(g_fit) %*% t(as.matrix(cbind(rep(1, n), X)))
        } else {
            g_fit <- MCMClogit(L ~ ., data = X[, 1:p_miss], b0 = 0, B0 = 0.01)
            g_mat <- as.matrix(g_fit) %*% t(as.matrix(cbind(rep(1, n), X[, 1:p_miss])))
        }

        if (m_val == "correct") {
            m_fit <- MCMClogit(Y[L == 1] ~ ., data = X[L == 1, ], b0 = 0, B0 = 0.01)
            m_mat <- as.matrix(m_fit) %*% t(as.matrix(cbind(rep(1, n), X)))
        } else {
            m_fit <- MCMClogit(Y[L == 1] ~ ., data = X[L == 1, 1:p_miss], b0 = 0, B0 = 0.01)
            m_mat <- as.matrix(m_fit) %*% t(as.matrix(cbind(rep(1, n), X[, 1:p_miss])))
        }

        theta_mat <- (1 - plogis(g_mat)) * plogis(m_mat)
        theta_post <- rowMeans(theta_mat) / (1 - mean(L))

    } else {
        theta_post <- NULL
    }

    out <- list(posterior_est = pp, g_val = g_val, m_val = m_val, prior = prior, n = n, p = p_dim,
               param_post = theta_post)

    return(out)
}

funslave <- function(j){
    index <- (1:nrow(tasks))[(0:(nrow(tasks)-1)) %/% (nrow(tasks) / 8000) + 1 == j]
    out <- lapply(index, function(k){
        cat('doing task ', k, ' of ', nrow(tasks), '\n', file = 'progress.txt')
        out <- try(sim(k))
        if(inherits(out, 'try-error'))
            cat('error in task ', k, '\n', file = 'errors.txt')
        return(out)
    })

    return(out)

}
