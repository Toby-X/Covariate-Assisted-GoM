.proj <- function(y, a = 1) {
    v <- y
    rho <- (sum(y) - a) / length(y)
    vt_len <- length(v)
    repeat{
        v <- v[which(v > rho)]
        v_len <- length(v)
        rho <- (sum(v) - a) / v_len
        if (v_len == vt_len) {
            break
        } else {
            vt_len <- v_len
        }
    }
    tau <- rho
    # K <- length(v)
    y <- pmax(y - tau, 0)
    return(y)
}

Prune <- function(U, r, q, e) {
    l <- apply(U, 1, norm, "2")
    quan_l <- quantile(l, 1 - q)
    P0 <- which(l >= quan_l)
    U0 <- U[P0, ]
    x <- rowMeans(FNN::knnx.dist(U, U0, r))
    quan_x <- quantile(x, 1 - e)
    P <- which(x > quan_x)
    P <- P0[P]
    return(P)
}

HeteroPCA <- function(Sigma, K, max.iter = 20, eps = 1e-8) {
    Sigma0 <- Sigma - diag(diag(Sigma))
    Sigma <- Sigma0
    dev <- 10
    iter <- 0
    while (dev > eps && iter < max.iter) {
        evd.tmp <- RSpectra::svds(Sigma, K)
        Sigma.tmp <- evd.tmp$u %*% (evd.tmp$d * t(evd.tmp$v))
        dev <- abs(mean(diag(Sigma) - diag(Sigma.tmp)))
        Sigma <- Sigma0 + diag(diag(Sigma.tmp))
        iter <- iter + 1
    }
    evd <- RSpectra::svds(Sigma, K)
    return(list(U = evd$u, iter = iter))
}
