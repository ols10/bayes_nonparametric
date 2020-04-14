## Ivan's functions without extra code

source("utils.r")

################################################
#   LFS parametric submodels and derivatives   #
################################################

pmod <- function(g, m, p){
  Dx <- (1 - g) / mean(1 - g) * (m - theta(g, m, p))
  return(function(eps) p * exp(eps * Dx) / sum(p * exp(eps * Dx)))
}

dpmod <- function(g, m, p){
  Dx <- (1 - g) / mean(1 - g) * (m - theta(g, m, p))
  pfun <- pmod(g, m, p)
  return(function(eps)
    pfun(eps) * (Dx - sum(Dx * exp(eps * Dx) * p) / sum(exp(eps * Dx) * p)))
}

mmod <- function(g, m, p){
  H <- (1 - g) / g * 1 / mean(1 - g)
  return(function(eps) plogis(qlogis(m) + eps * H))
}

dmmod <- function(g, m, p){
  H <- (1 - g) / g * 1 / mean(1 - g)
  mfun <- mmod(g, m, p)
  return(function(eps) H * mfun(eps) * (1 - mfun(eps)))
}

gmod <- function(g, m, p){
  M <- - (m - theta(g, m, p)) / mean(1 - g)
  return(function(eps) plogis(qlogis(g) + eps * M))
}

dgmod <- function(g, m, p){
  M <- - (m - theta(g, m, p)) / mean(1 - g)
  gfun <- gmod(g, m, p)
  return(function(eps) M * gfun(eps) * (1 - gfun(eps)))
}

#################################
#     parameter definition      #
#################################

theta <- function(g, m, p) sum((1 - g) * m * p) / sum((1 - g) * p)

#############################################
#     transformation eps -> theta(eps)      #
#############################################

thetaf <- function(g, m, p){
  Vectorize(function(eps) theta(gmod(g, m, p)(eps), mmod(g, m, p)(eps), pmod(g, m, p)(eps)))
}

##########################################
#       jacobian of theta(eps)           #
##########################################

dthetaf  <- function(g, m, p){
  r <- mean(g)
  dgmodfun <- dgmod(g, m, p)
  dmmodfun <- dmmod(g, m, p)
  dpmodfun <- dpmod(g, m, p)
  gmodfun <- gmod(g, m, p)
  mmodfun <- mmod(g, m, p)
  pmodfun <- pmod(g, m, p)
  return(function(eps) {
    sum(- dgmodfun(eps) * mmodfun(eps) * pmodfun(eps) +
          (1 - gmodfun(eps)) * dmmodfun(eps) * pmodfun(eps) +
          (1 - gmodfun(eps)) * mmodfun(eps) * dpmodfun(eps))
  })
}


##########################################
#       efficient influence fn           #
##########################################

infun <- function(Y, L, g, m, p) {
  H <- (1 - g) / g * 1 / sum((1 - g) * p)
  M <- - (m - theta(g, m, p)) / sum((1 - g) * p)
  return(H * L * (Y - m) - M * (1 - L))
}

#########################################
#            initial TMLE               #
#########################################

tmle <- function(Y, L, m, g) {

  n <- length(L)
  m <- bound01(m)
  g <- bound01(g)

  iter <- 1
  crit <- TRUE

  while(crit & iter < 20) {

    H <- (1 - g) / g * 1 / mean(1 - g)
    M <- - (m - theta(g, m, 1 / n)) / mean(1 - g)
    eps1 <- glm(Y ~ 0 + offset(qlogis(m)) + H, subset = L == 1, family = binomial)
    eps2 <- glm(L ~ 0 + offset(qlogis(g)) + M, family = binomial)
    m <- bound01(predict(eps1, newdata = data.frame(m = m, H = H), type = 'response'))
    g <- bound01(predict(eps2, newdata = data.frame(g = g, M = M), type = 'response'))
    crit <- abs(mean(infun(Y, L, g, m, 1/n))) > 1e-3 / n
    iter <- iter + 1

  }

  return(list(m = m, g = g, estimate = theta(g, m, 1/n),
              infun = infun(Y, L, g, m, 1/n)))
}

################################
#   posterior distribution     #
################################

posterior <- function(Y, L, m, g, pi, nsim = 1e4) {

  n <- length(Y)
  p <- 1 / n

  etmle <- tmle(Y, L, m, g)
  m <- etmle$m
  g <- etmle$g

  gmodfun <- gmod(g, m, p)
  mmodfun <- mmod(g, m, p)
  pmodfun <- pmod(g, m, p)
  thetafun <- thetaf(g, m, p)
  jacobian <- dthetaf(g, m, p)


  loglikelihood <- function(eps){
    loglik   <- log(((mmodfun(eps)^Y) * (1 - mmodfun(eps))^(1 - Y))^L *
                      (gmodfun(eps)^L) * (1 - gmodfun(eps))^(1 - L) * pmodfun(eps))
    return(-sum(loglik))
  }
  loglikelihood <- Vectorize(loglikelihood)

  pieps <- function(eps) abs(jacobian(eps)) * pi(thetafun(eps))
  pieps <- Vectorize(pieps)

  posteps <- function(eps) log(pieps(eps)) - loglikelihood(eps)

  library(mcmc)
  epsilon <- metrop(posteps, 0, nsim)$batch
  ## epsilon   <- mh(nsim, posteps, 0, .1)
  thetaval  <- thetafun(epsilon)

  output    <- list(tmle = etmle, post = data.frame(eps = epsilon, theta = thetaval))
  return(output)

}



