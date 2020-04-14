bound01 <- function(p) pmax(pmin(p, 1 - 1e-9), 1e-9)

rebeta <- function(n, a, b, shape1, shape2){
  return(a + rbeta(n, shape1, shape2) * (b - a))
}

debeta <- function(x, a, b, shape1, shape2){
  (x - a)^(shape1 - 1) * (b - x)^(shape2 - 1) /
      (beta(shape1, shape2) * (b - a)^(shape1 + shape2 - 1))
}

pebeta <- function(x, a, b, shape1, shape2){
  dens <- function(y) drbeta(y, a, b, shape1, shape2)
  return(integrate(dens, a, x)$value)
}

mfun <- function(X) {
    p <- ncol(X)
    if(p == 4) fac <- 2 else fac <- 10
    lin <- - rowSums(X) / fac
    return(plogis(- 1/2 - lin))
}

gfun <- function(X) {
    p <- ncol(X)
    if(p == 4) fac <- 2 else fac <- 10
    lin <- rowSums(X) / fac
    return(plogis(lin))
}

gendata <- function(n, p = 100) {

    X1 <- matrix(runif(n * p / 2, -1, 1), nrow = n)
    X2 <- matrix(rnorm(n * p / 2), nrow = n)

    X <- cbind(X1, X2)

    L <- rbinom(n, 1, gfun(X))
    Y <- rbinom(n, 1, mfun(X))
    Y[L == 0] <- 9999

    return(list(X = X, L = L, Y = Y))
}

true <- function(n, p = 100) {

    X1 <- matrix(runif(n * p / 2, -1, 1), nrow = n)
    X2 <- matrix(rnorm(n * p / 2), nrow = n)

    X <- cbind(X1, X2)

    truth <- mean((1 - gfun(X)) * mfun(X)) / mean(1 - gfun(X))
    effb  <- mean((1 - gfun(X)) * ((1 - gfun(X)) / gfun(X) * mfun(X) * (1 - mfun(X))
                  + (mfun(X) - truth)^2)) / (mean(1 - gfun(X)))^2

    return(c(p = p, truth = truth, effb = effb))

}


beta_reparam <- function(mu, var) {
    alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta = alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
}
