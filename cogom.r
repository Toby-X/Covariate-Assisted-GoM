source("utils.R")

.vertex_hunting <- function(U) {
  K <- ncol(U)
  S <- c()
  Y <- U
  IK <- diag(rep(1, K))
  for (k in 1:K) {
    l <- apply(Y, 1, norm, "2")
    S <- c(S, which.max(l))
    u <- Y[S[k], ] / norm(Y[S[k], ], "2")
    u <- data.matrix(u)
    Y <- Y %*% (IK - u %*% t(u))
  }
  
  S
}

gom.svd <- function(A, K, eps=1e-3) {
  svd.A <- RSpectra::svds(A, K)
  S <- .vertex_hunting(svd.A$u)
  US <- svd.A$u[S, ]
  Pi <- svd.A$u %*% solve(US)
  Pi <- t(apply(Pi, 1, .proj))
  Theta <- svd.A$v %*% (svd.A$d * t(svd.A$u)) %*% (Pi %*% solve(t(Pi) %*% Pi))
  Theta <- pmin(Theta, 1 - eps)
  Theta <- pmax(Theta, eps)

  list(Pi = Pi, Theta = Theta, S = S)
}

gom.cov.svd <- function(A, X, K, alpha, eps=1e-3) {
  L <- cbind(A, sqrt(alpha) * X)
  svd_L <- RSpectra::svds(L, K)

  S <- .vertex_hunting(svd_L$u)
  US <- svd_L$u[S, ]
  Pi <- svd_L$u %*% solve(US)
  Pi <- t(apply(Pi, 1, .proj))

  Theta_M <- svd_L$v %*% (svd_L$d * t(svd_L$u)) %*% (Pi %*% solve(t(Pi) %*% Pi))
  Theta <- Theta_M[, 1:ncol(A)]
  M <- Theta_M[, (ncol(A) + 1):ncol(L)] / sqrt(alpha)
  Theta <- pmin(Theta, 1 - eps)
  Theta <- pmax(Theta, eps)

  list(Pi = Pi, Theta = Theta, M = M, S = S)
}

gom.cov <- function(A, X, K, alpha, hiter=10, eps = 1e-3) {
  G_A <- HeteroPCA(A %*% t(A), K, eps = eps)
  G_X <- HeteroPCA(X %*% t(X), K, eps = eps)
  G <- G_A + alpha * G_X

  svd_G <- RSpectra::svds(G, K)
  S <- .vertex_hunting(svd_G$u)
  US <- svd_G$u[S, ]
  Pi <- svd_G$u %*% solve(US)
  Pi <- t(apply(Pi, 1, .proj))
  
  res_svd <- gom.cov.svd(A, K, eps)
  loss_mae <- t(apply(Pi, 2, function(x) colMeans(abs(x - res_svd$Pi))))
  od <- clue::solve_LSAP(loss_mae)

  list(Pi = Pi, Theta = res_svd$Theta[, od], M = res_svd$M[,od], S = res_svd$S)
}

find_best_alpha <- function(alpha_seq, A, X, K, eps = 1e-3, nfold = 5) {
  err <- sapply(alpha_seq, .recovery.cv, A, X, K, eps, nfold)
  model.loess <- stats::loess(err ~ alpha_seq)
  idx.best <- which.min(model.loess$fitted)
  
  alpha_seq[idx.best]
}

.recovery.cv <- function(alpha, A, X, K, eps = 1e-3, nfold = 5) {
  ## create K-fold missing matrix
  n <- nrow(A)
  idx.A <- rep(1:nfold, length.out = nrow(A) * ncol(A))
  idx.X <- rep(1:nfold, length.out = nrow(X) * ncol(X))
  idx.A <- sample(idx.A)
  idx.X <- sample(idx.X)
  err <- rep(0, nfold)

  for (kk in 1:nfold) {
    A.tmp <- A
    X.tmp <- X
    bool_mat.A <- matrix(idx.A == kk, nrow = n)
    bool_mat.X <- matrix(idx.X == kk, nrow = n)
    A.tmp[as.logical(bool_mat.A)] <- 0
    X.tmp[as.logical(bool_mat.X)] <- 0
    A.tmp <- A.tmp / sum(idx.A!=kk) * length(idx.A)
    X.tmp <- X.tmp / sum(idx.X!=kk) * length(idx.X)
    cov.tmp <- gom.cov(A.tmp, X.tmp, K, alpha, eps)
    A.est <- cov.tmp$Pi %*% t(cov.tmp$Theta)
    err[kk] <- mean(abs(A[bool_mat.A] - A.est[bool_mat.A]))
  }

  mean(err)
}
