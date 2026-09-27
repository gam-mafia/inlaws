library(inlaws)
jet <- getFromNamespace(".zinb_logjet", "inlaws")
spec <- getFromNamespace(".zinb_jet_spec", "inlaws")
near <- function(x, y, tol = 2e-5) {
  stopifnot(
    length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
    max(abs(x - y) / pmax(1, abs(x), abs(y))) < tol
  )
}
differentiate <- function(f, x, j, h = 1e-4) {
  e <- numeric(length(x))
  e[j] <- h
  (f(x - 2 * e) - 8 * f(x - e) + 8 * f(x + e) - f(x + 2 * e)) / (12 * h)
}
grid <- expand.grid(y = c(0, 1, 7, 51, 1000), u = c(-18, -4, 1, 4), b = c(-30, -2, 1, 30), v = c(-3, 1, 8, 23))
eta <- as.matrix(grid[, c("u", "b", "v")])
y <- grid$y
s <- spec(4)
packed <- lapply(1:4, function(r) {
  ii <- which(rowSums(s$a) == r)
  s$a[ii[order(-s$a[ii, 1], -s$a[ii, 2])], , drop = FALSE]
})
for (h in c(FALSE, TRUE)) {
  d <- jet(y, eta, h, 4)
  for (j in 1:3) {
    ep <- em <- eta
    ep[, j] <- ep[, j] + 1e-4
    em[, j] <- em[, j] - 1e-4
    dp <- jet(y, ep, h, 3)
    dm <- jet(y, em, h, 3)
    near((dp$l0 - dm$l0) / 2e-4, d$l1[, j])
    for (r in 1:3) {
      for (k in seq_len(nrow(packed[[r]]))) {
        a <- packed[[r]][k, ]
        a[j] <- a[j] + 1
        ix <- which(apply(packed[[r + 1]], 1, function(b) all(a == b)))
        near((dp[[r + 1]][, k] - dm[[r + 1]][, k]) / 2e-4, d[[r + 2]][, ix])
      }
    }
  }
  if (h) {
    for (r in 2:4) {
      mixed <- packed[[r]][, 2] > 0 & packed[[r]][, 2] < r
      stopifnot(all(d[[r + 1]][, mixed] == 0))
    }
    # Tiny-mean positive counts concentrate at one; curvature must not vanish
    # by subtracting two rounded unit scores.
    z <- jet(1, matrix(c(-35, 0, 1), 1), TRUE, 4)
    stopifnot(z$l1[1, 1] < 0, z$l2[1, 1] < 0)
    near(z$l1[1, 1] / exp(-35), -(1 + exp(-1)) / 2, 1e-6)
  }
}

# Independent coefficient and smoothing-direction contractions, including
# overlapping coefficient indices and weighted sandwich scores.
set.seed(130)
for (h in c(FALSE, TRUE)) {
  for (shared in c(FALSE, TRUE)) {
    fam <- if (h) zanb() else zinb()
    n <- 12
    X <- cbind(1, seq(-1, 1, length.out = n), 1, runif(n), 1, rnorm(n))
    attr(X, "lpi") <- list(1:2, if (shared) c(2, 3, 4) else 3:4, 5:6)
    beta <- c(1, .2, .5, -.1, .7, .1)
    y <- rep(c(0, 1, 4, 10), 3)
    wt <- c(0, .4, rep(1, 5), rep(2, 5))
    off <- replicate(3, runif(n, -.2, .2), simplify = FALSE)
    ll <- function(b, deriv = 1, ...) fam$ll(y, X, b, wt, fam, offset = off, deriv = deriv, ...)
    first <- ll(beta)
    for (j in seq_along(beta)) {
      near(first$lb[j], differentiate(function(b) ll(b, 0)$l, beta, j))
      near(first$lbb[, j], differentiate(function(b) ll(b)$lb, beta, j))
    }
    directions <- matrix(rnorm(12), 6, 2)
    second <- matrix(rnorm(18), 6, 3)
    V <- crossprod(matrix(rnorm(36), 6)) + diag(6)
    R <- chol(solve(V), pivot = TRUE)
    full <- ll(beta, 4, d1b = directions, d2b = second, fh = R, D = rep(1, 6))
    contracted <- ll(beta, 2, d1b = directions, fh = V)
    for (j in 1:2) {
      fd <- differentiate(function(a) ll(beta + a * directions[, j])$lbb, 0, 1)
      near(full$d1H[[j]], fd)
      near(contracted$d1H[j], sum(V * fd))
    }
    kk <- 0
    for (j in 1:2) {
      for (k in j:2) {
        kk <- kk + 1
        fd <- differentiate(function(a) ll(beta + a * directions[, k], 3, d1b = directions)$d1H[[j]], 0, 1)
        fd <- fd + differentiate(function(a) ll(beta + a * second[, kk])$lbb, 0, 1)
        near(full$trHid2H[kk], sum(V * fd))
      }
    }
    scores <- t(vapply(seq_len(n), function(i) {
      w <- numeric(n)
      w[i] <- wt[i]
      fam$ll(y, X, beta, w, fam, offset = off, deriv = 1)$lb
    }, numeric(6)))
    near(fam$sandwich(y, X, beta, wt, fam, off), crossprod(scores))
    jj <- attr(X, "lpi")
    eta <- sapply(1:3, function(j) drop(X[, jj[[j]], drop = FALSE] %*% beta[jj[[j]]]) + off[[j]])
    near(first$l, fam$ll(y, X, beta, wt, fam, eta = eta)$l)
    eta[1, ] <- Inf
    near(first$l, fam$ll(y, X, beta, wt, fam, eta = eta)$l)
  }
}

# Continuity through both truncation-series switch points, through order four.
for (mu in c(.01, expm1(.1))) {
  e <- matrix(c(log(mu), .4, 0), 1)
  ep <- em <- e
  ep[, 1] <- ep[, 1] + 1e-8
  em[, 1] <- em[, 1] - 1e-8
  dp <- jet(2, ep, TRUE, 4)
  dm <- jet(2, em, TRUE, 4)
  for (k in seq_along(dp)) near(dp[[k]], dm[[k]], 1e-6)
}
