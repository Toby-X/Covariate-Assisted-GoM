install.packages(setdiff(
    c("doSNOW", "tictoc", "gtools", "sirt"),
    rownames(installed.packages())
))

source("../utils.r")
source("../cogom.r")
library(doSNOW)
library(gtools)
library(tictoc)
library(sirt)

## parallel computing parameters
numCores <- 32L
cl <- makeCluster(numCores)
registerDoSNOW(cl)

# parallel experiments
m <- 100

pi_gen <- function(pi, K) {
    pi <- rdirichlet(1, rep(1, K))
    return(pi)
}

GoM_simulate <- function(A, X, K, Pi, Theta) {
    ## StackSVD
    stack.est <- gom.cov(A, X, K, 1)
    idx.stack <- find_best_idx(stack.est$Pi, Pi)
    pi.err.stack <- mean(abs(stack.est$Pi[, idx.stack] - Pi))
    theta.err.stack <- mean(abs(stack.est$Theta[, idx.stack] - Theta))
    pi.err2.stack <- max(apply(stack.est$Pi[, idx.stack] - Pi, 1, norm, "2"))
    theta.err2.stack <- max(apply(stack.est$Theta[, idx.stack] - Theta, 1, norm, "2"))


    ## Covariate-assisted include the time of finding alpha
    alpha_seq <- seq(from = 0, to = 1, length = 20)
    alpha0 <- find_best_alpha(alpha_seq, A, X, K)
    cov.est <- gom.cov(A, X, K, alpha0)
    idx.cov <- find_best_idx(cov.est$Pi, Pi)
    pi.err.cov <- mean(abs(cov.est$Pi[, idx.cov] - Pi))
    theta.err.cov <- mean(abs(cov.est$Theta[, idx.cov] - Theta))
    pi.err2.cov <- max(apply(cov.est$Pi[, idx.cov] - Pi, 1, norm, "2"))
    theta.err2.cov <- max(apply(cov.est$Theta[, idx.cov] - Theta, 1, norm, "2"))

    res <- list(
        stack.pi = pi.err.stack, stack.theta = theta.err.stack,
        stack.pi2 = pi.err2.stack, stack.theta2 = theta.err2.stack,
        cov.pi = pi.err.cov, cov.theta = theta.err.cov,
        cov.pi2 = pi.err2.cov, cov.theta2 = theta.err2.cov,
        alpha.cov = alpha0
    )

    res
}

# N=500,K=3 Simulation
res_N500_s05 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "RSpectra")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N / 10
    W <- N / 20

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    X <- X_t + matrix(rnorm(N * W, 0, .5), nrow = N)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

res_N500_s1 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "RSpectra")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N / 10
    W <- N / 20

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    X <- X_t + matrix(rnorm(N * W, 0, 1), nrow = N)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

res_N500_s1.5 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "RSpectra")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N / 10
    W <- N / 20

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    X <- X_t + matrix(rnorm(N * W, 0, 1.5), nrow = N)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

res_N500_s2 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "RSpectra")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N / 10
    W <- N / 20

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    X <- X_t + matrix(rnorm(N * W, 0, 2), nrow = N)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

res_N500_JW1 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "RSpectra")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N / 10
    W <- N / 10

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    X <- X_t + matrix(rnorm(N * W, 0, 1), nrow = N)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

res_N500_JW5 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "RSpectra")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N / 10
    W <- ceiling(N / 50)

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    X <- X_t + matrix(rnorm(N * W, 0, 1), nrow = N)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

res_N500_JW10 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "RSpectra")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N / 10
    W <- ceiling(N / 100)

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    X <- X_t + matrix(rnorm(N * W, 0, 1), nrow = N)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}


save.image("GoM_com_stack.RData")
