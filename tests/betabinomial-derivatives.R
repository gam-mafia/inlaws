library(inlaws)

near <- function(x, y, tol = 2e-6) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
            max(abs(x - y)) <= tol * max(1, abs(x), abs(y)))
}
differentiate <- function(f, x, j, h = 1e-3) {
  e <- numeric(length(x)); e[j] <- h
  (f(x - 2 * e) - 8 * f(x - e) + 8 * f(x + e) - f(x + 2 * e)) / (12 * h)
}
dd <- getFromNamespace(".bb_derivatives", "inlaws")
rising <- getFromNamespace(".bb_rising", "inlaws")

# An independent finite-product likelihood is stable even near the binomial
# limit, unlike direct subtraction of lbeta values.
logp <- function(y, m, eta) {
  mu <- plogis(eta[1]); phi <- exp(eta[2])
  prodlog <- function(n, a) sum(log1p((seq_len(n) - 1) / a))
  dbinom(y, m, mu, log = TRUE) + prodlog(y, mu * phi) +
    prodlog(m - y, (1 - mu) * phi) - prodlog(m, phi)
}
for (m in c(1, 20, 64, 65, 500)) for (phi in c(1e-6, .01, 1, 20, 2000, 1e12)) {
  for (mu in c(.0001, .2, .5, .99)) for (y in unique(c(0, round(m * .3), m))) {
    eta <- c(qlogis(mu), log(phi))
    a <- dd(y, m, matrix(eta, 1L))
    near(a$l0, logp(y, m, eta), 1e-8)
    for (r in 1:4) for (j in 1:2) {
      f <- function(e) {
        if (r == 1L) return(logp(y, m, e))
        dd(y, m, matrix(e, 1L), r - 1L)[[paste0("l", r - 1L)]]
      }
      near(differentiate(f, eta, j), a[[paste0("l", r)]][, seq_len(r) + j - 1L], 3e-5)
    }
  }
}
# Check both numerical switching boundaries and their derivative continuity.
for (m in c(65, 1000)) for (a in c(8, m / .05)) {
  eps <- 1e-8
  near(rising(log(a) - eps, m), rising(log(a) + eps, m), 1e-7)
}
for (m in c(2, 20, 100)) for (mu in c(.1, .5, .9)) for (phi in c(.01, 3, 1e10)) {
  y <- 0:m; eta <- cbind(rep(qlogis(mu), m + 1), rep(log(phi), m + 1))
  p <- exp(dd(y, rep(m, m + 1), eta, 0L)$l0)
  near(sum(p), 1, 1e-8)
  near(sum(y * p), m * mu, 1e-8)
  near(sum((y - m * mu)^2 * p), m * mu * (1 - mu) * (1 + (m - 1) / (1 + phi)), 1e-8)
  if (phi == 1e10) near(p, dbinom(y, m, mu), 1e-7)
}
for (m in 0:1) {
  a <- dd(m, m, matrix(c(-2, -20), 1L))
  for (r in 1:4) stopifnot(all(a[[paste0("l", r)]][, -1L] == 0))
}
stopifnot(all(is.finite(unlist(dd(c(0, 20, 8), rep(20, 3),
  cbind(c(-700, 700, 30), c(-700, -700, 700)))))))

# Coefficient contractions, weights, offsets, and shared predictor coefficients.
set.seed(725)
for (shared in c(FALSE, TRUE)) {
  nr <- 14L
  X <- cbind(1, runif(nr), 1, rnorm(nr))
  attr(X, "lpi") <- if (shared) list(c(1, 2), c(2, 3, 4)) else list(1:2, 3:4)
  y <- sample(0:8, nr, TRUE); m <- sample(10:90, nr, TRUE)
  wt <- c(0, runif(nr - 1, .5, 2)); off <- list(runif(nr, -.2, .2), runif(nr, -.2, .2))
  family <- betabinomial(); family$trials <- m
  beta <- c(-.2, .3, 1, .2)
  ll <- family$ll(y, X, beta, wt, family, off, deriv = 1)
  for (j in 1:4) {
    near(differentiate(function(b) family$ll(y, X, b, wt, family, off)$l, beta, j), ll$lb[j])
    near(differentiate(function(b) family$ll(y, X, b, wt, family, off, deriv = 1)$lb,
                       beta, j), ll$lbb[, j])
  }
  d1b <- matrix(rnorm(8), 4, 2); d2b <- matrix(rnorm(12), 4, 3)
  V <- crossprod(matrix(rnorm(16), 4)) + diag(4)
  D <- runif(4, .5, 2); Hi <- D * t(D * V)
  ll3 <- family$ll(y, X, beta, wt, family, off, deriv = 3, d1b = d1b)
  ll4 <- family$ll(y, X, beta, wt, family, off, deriv = 4, d1b = d1b,
                  d2b = d2b, fh = chol(solve(V), pivot = TRUE), D = D)
  ll4e <- family$ll(y, X, beta, wt, family, off, deriv = 4, d1b = d1b,
                   d2b = d2b, fh = eigen(solve(V)), D = D)
  for (j in 1:2) {
    near(differentiate(function(z) family$ll(y, X, beta + d1b %*% z, wt, family,
      off, deriv = 1)$lbb, c(0, 0), j), ll3$d1H[[j]])
  }
  kk <- 0L
  for (j in 1:2) for (k in j:2) {
    kk <- kk + 1L
    fd <- differentiate(function(z) family$ll(y, X, beta + d1b %*% z, wt,
      family, off, deriv = 3, d1b = d1b)$d1H[[j]], c(0, 0), k)
    fd <- fd + differentiate(function(z) family$ll(y, X, beta + z * d2b[, kk],
      wt, family, off, deriv = 1)$lbb, 0, 1)
    near(ll4$trHid2H[kk], sum(Hi * fd), 1e-5)
    near(ll4e$trHid2H[kk], sum(Hi * fd), 1e-5)
  }
  stopifnot(all(is.finite(unlist(ll4))))
  # Inactive rows remain harmless even when their predictors overflow.
  eta <- matrix(0, nr, 2); eta[1, ] <- Inf
  stopifnot(family$ll(y, X, beta, wt, family, eta = eta, deriv = 1)$l0[1] == 0)
}

# Precision derivatives retain their leading term at the binomial limit.
for (phi in c(1e8, 1e12, 1e16)) {
  a <- dd(0, 20, matrix(c(qlogis(.3), log(phi)), 1))
  expected <- choose(20, 2) * .3 / .7
  near(a$l1[1, 2] * phi, -expected, 1e-6)
  near(a$l2[1, 3] * phi, expected, 1e-6)
}
