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
    ## Covariate-assisted include the time of finding alpha
    alpha_seq <- seq(from = 0, to = 1, length = 20)
    tic()
    alpha0 <- find_best_alpha(alpha_seq, A, X, K)
    cov.est <- gom.cov.svd(A, X, K, alpha0) # there is no hollow, but hollowing, needs to change
    cov.time <- toc()
    cov.time <- cov.time$toc - cov.time$tic
    idx.cov <- find_best_idx(cov.est$Pi, Pi)
    pi.err.cov <- mean(abs(cov.est$Pi[, idx.cov] - Pi))
    theta.err.cov <- mean(abs(cov.est$Theta[, idx.cov] - Theta))
    pi.err2.cov <- max(apply(cov.est$Pi[, idx.cov] - Pi, 1, norm, "2"))
    theta.err2.cov <- max(apply(cov.est$Theta[, idx.cov] - Theta, 1, norm, "2"))

    ## Covariate-assisted include the time of finding alpha
    ## the tuning params are all set as follows
    # alpha_seq = seq(from=0,to=1,length=100)
    tic()
    alpha2 <- find_best_alpha(alpha_seq, A, X, K) # the code needs modification
    hetero.est <- gom.cov(A, X, K, alpha2)
    hetero.time <- toc()
    hetero.time <- hetero.time$toc - hetero.time$tic
    idx.hetero <- find_best_idx(hetero.est$Pi, Pi)
    pi.err.hetero <- mean(abs(hetero.est$Pi[, idx.hetero] - Pi))
    theta.err.hetero <- mean(abs(hetero.est$Theta[, idx.hetero] - Theta))
    pi.err2.hetero <- max(apply(hetero.est$Pi[, idx.hetero] - Pi, 1, norm, "2"))
    theta.err2.hetero <- max(apply(hetero.est$Theta[, idx.hetero] - Theta, 1, norm, "2"))

    res <- list(
        cov.pi = pi.err.cov, cov.theta = theta.err.cov, cov.time = cov.time,
        cov.pi2 = pi.err2.cov, cov.theta2 = theta.err2.cov,
        hetero.pi = pi.err.hetero, hetero.theta = theta.err.hetero,
        hetero.time = hetero.time,
        hetero.pi2 = pi.err2.hetero, hetero.theta2 = theta.err2.hetero,
        alpha.hetero = alpha2, alpha.cov = alpha0
    )
    # alpha is also returned to check if the range is good enough
    return(res)
}

# N=200,K=3 Simulation
hetero_N200_K3 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "mvtnorm")
) %dopar% {
    set.seed(i)
    N <- 200
    K <- 3
    J <- N * 20
    W <- N * 10
    beta <- 8

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    v <- runif(N)
    sigma2 <- N * v^beta / sum(v^beta)
    X <- X_t + t(rmvnorm(W, sigma = diag(sigma2)))
    X <- apply(X, 2, scale)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

# N=500,K=3 Simulation
hetero_N500_K3 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "mvtnorm")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N * 20
    W <- N * 10
    beta <- 8

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    v <- runif(N)
    sigma2 <- N * v^beta / sum(v^beta)
    X <- X_t + t(rmvnorm(W, sigma = diag(sigma2)))
    X <- apply(X, 2, scale)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

# N=1000,K=3 Simulation
hetero_N1000_K3 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc", "mvtnorm")
) %dopar% {
    set.seed(i)
    N <- 1000
    K <- 3
    J <- N * 20
    W <- N * 10
    beta <- 8

    # generate values
    Pi <- matrix(rep(0, N * K), nrow = N)
    Pi <- t(apply(Pi, 1, pi_gen, K = K))
    Theta <- matrix(runif(J * K), ncol = K)
    M <- matrix(rnorm(W * K, 0, 1), nrow = W)
    Pi[1:K, ] <- diag(rep(1, K))

    A_t <- Pi %*% t(Theta)
    X_t <- Pi %*% t(M)

    A <- matrix(rbinom(N * J, 1, A_t), nrow = N)
    v <- runif(N)
    sigma2 <- N * v^beta / sum(v^beta)
    X <- X_t + t(rmvnorm(W, sigma = diag(sigma2)))
    X <- apply(X, 2, scale)

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

# N=200,K=3 Simulation
missing_N200_K3 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc")
) %dopar% {
    set.seed(i)
    N <- 200
    K <- 3
    J <- N * 20
    W <- N * 10
    # sample rate
    p <- 0.2

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
    X <- apply(X, 2, scale)

    A.mis <- matrix(rbinom(N * J, 1, p), nrow = N)
    X.mis <- matrix(rbinom(N * W, 1, p), nrow = N)
    A[A.mis == 0] <- 0
    X[X.mis == 0] <- 0

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

# N=500,K=3 Simulation
missing_N500_K3 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc")
) %dopar% {
    set.seed(i)
    N <- 500
    K <- 3
    J <- N * 20
    W <- N * 10
    # sample rate
    p <- 0.2

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
    X <- apply(X, 2, scale)

    A.mis <- matrix(rbinom(N * J, 1, p), nrow = N)
    X.mis <- matrix(rbinom(N * W, 1, p), nrow = N)
    A[A.mis == 0] <- 0
    X[X.mis == 0] <- 0

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}

# N=1000,K=3 Simulation
missing_N1000_K3 <- foreach(
    i = 1:m, .combine = rbind,
    .packages = c("sirt", "gtools", "tictoc")
) %dopar% {
    set.seed(i)
    N <- 1000
    K <- 3
    J <- N * 20
    W <- N * 10
    # sample rate
    p <- 0.2

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
    X <- apply(X, 2, scale)

    A.mis <- matrix(rbinom(N * J, 1, p), nrow = N)
    X.mis <- matrix(rbinom(N * W, 1, p), nrow = N)
    A[A.mis == 0] <- 0
    X[X.mis == 0] <- 0

    res <- GoM_simulate(A, X, K, Pi, Theta)
    return(res)
}


save.image("GoM_hetero.RData")
