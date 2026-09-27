library(inlaws)

near <- function(x, y, tol = 2e-6) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
            max(abs(x - y)) <= tol * max(1, abs(x), abs(y)))
}
differentiate <- function(f, x, j, h = 1e-4 * max(abs(x[j]), 0.01)) {
  e <- numeric(length(x)); e[j] <- h
  (f(x - 2 * e) - 8 * f(x - e) + 8 * f(x + e) - f(x + 2 * e)) / (12 * h)
}
natural <- getFromNamespace(".nbls_derivatives", "inlaws")

# Check each derivative against a derivative of the preceding order. Scale
# both coordinates multiplicatively so tiny theta derivatives remain tested.
grid <- expand.grid(y = c(0, 1, 7, 51, 500), mu = c(0.03, 3, 80),
                    theta = c(0.02, 2, 150, 1e5, 1e10))
for (i in seq_len(nrow(grid))) {
  y <- grid$y[i]; mu <- grid$mu[i]; theta <- grid$theta[i]
  dd <- natural(y, mu, theta)
  for (r in 1:4) for (q in 0:r) {
    p <- r - q
    coordinate <- if (p > 0) 1L else 2L
    pp <- p - (coordinate == 1L); qq <- q - (coordinate == 2L)
    f <- function(z) {
      if (r == 1L) {
        # An independent finite-product likelihood avoids numerical noise in
        # dnbinom's size derivative very close to the Poisson limit.
        m <- mu * z[1]; t <- theta * z[2]
        return(sum(log1p((seq_len(y) - 1) / t)) - lgamma(y + 1) +
                 y * log(m) - (t + y) * log1p(m / t))
      }
      natural(y, mu * z[1], theta * z[2], r - 1L)[[r - 1L]][1L, qq + 1L] *
        mu^pp * theta^qq
    }
    fd <- differentiate(f, c(1, 1), coordinate)
    near(fd, dd[[r]][1L, q + 1L] * mu^p * theta^q, tol = 2e-5)
  }
}

# Check link transformations and coefficient-level contractions, with offsets,
# fractional weights, and a zero-weight row. Include shared coefficients.
set.seed(130)
for (shared in c(FALSE, TRUE)) for (link in c("log", "identity", "sqrt")) {
  n <- 12L
  X <- cbind(1, seq(-1, 1, length.out = n), rnorm(n), 1, runif(n))
  attr(X, "lpi") <- list(1:3, if (shared) 3:5 else 4:5)
  beta <- c(2, .1, -.1, 1, .2)
  y <- rep(c(0, 2, 4, 10), 3)
  w <- c(0, .4, rep(1, 5), rep(2, 5))
  off <- list(runif(n, -.2, .2), runif(n, -.4, .4))
  fam <- nbls(list(link, "log"))
  ll <- function(b, deriv = 1L, ...) fam$ll(y, X, b, w, fam,
                                          offset = off, deriv = deriv, ...)
  first <- ll(beta)
  for (j in seq_along(beta)) {
    near(first$lb[j], differentiate(function(b) ll(b, 0)$l, beta, j))
    near(first$lbb[, j], differentiate(function(b) ll(b)$lb, beta, j))
  }
  directions <- matrix(rnorm(10), 5, 2)
  second <- matrix(rnorm(15), 5, 3)
  V <- crossprod(matrix(rnorm(25), 5)) + diag(5)
  R <- chol(solve(V), pivot = TRUE)
  full <- ll(beta, 4, d1b = directions, d2b = second, fh = R, D = rep(1, 5))
  contracted <- ll(beta, 2, d1b = directions, fh = V)
  for (j in 1:2) {
    fd <- differentiate(function(a) ll(beta + a * directions[, j])$lbb, 0, 1, h = 1e-4)
    near(full$d1H[[j]], fd)
    near(contracted$d1H[j], sum(V * fd))
  }
  kk <- 0L
  for (j in 1:2) for (k in j:2) {
    kk <- kk + 1L
    # Differentiate the third derivative contraction and add the change in the
    # coefficient direction induced by second derivatives of beta(rho).
    fd <- differentiate(function(a) ll(beta + a * directions[, k], 3,
      d1b = directions)$d1H[[j]], 0, 1, h = 1e-4)
    fd <- fd + differentiate(function(a) ll(beta + a * second[, kk])$lbb,
                              0, 1, h = 1e-4)
    near(full$trHid2H[kk], sum(V * fd), tol = 1e-5)
  }
  # Sandwich is the cross-product of weighted individual coefficient scores.
  scores <- t(vapply(seq_len(n), function(i) {
    wi <- numeric(n); wi[i] <- w[i]
    fam$ll(y, X, beta, wi, fam, offset = off, deriv = 1)$lb
  }, numeric(5)))
  near(fam$sandwich(y, X, beta, w, fam, off), crossprod(scores))
  # Omitting a zero-weight row must give exactly the same likelihood.
  X1 <- X[-1, , drop = FALSE]; attr(X1, "lpi") <- attr(X, "lpi")
  omitted <- fam$ll(y[-1], X1, beta, w[-1], fam,
                    offset = lapply(off, `[`, -1), deriv = 1)
  near(first$l, omitted$l)
  near(first$lb, omitted$lb)
  near(first$lbb, omitted$lbb)
  # Explicit eta and coefficient construction agree, including offsets.
  jj <- attr(X, "lpi")
  eta <- sapply(1:2, function(j) drop(X[, jj[[j]], drop = FALSE] %*%
                                      beta[jj[[j]]]) + off[[j]])
  near(first$l, fam$ll(y, X, beta, w, fam, eta = eta)$l)
  eta[1, ] <- c(Inf, -Inf)
  near(first$l, fam$ll(y, X, beta, w, fam, eta = eta)$l)
}

# Fixed-size mean derivatives agree with nb()'s deviance derivatives.
y <- c(0, 1, 3, 20); mu <- c(.3, 1, 5, 12); theta <- 2.3
nbd <- mgcv::nb(theta = theta)$Dd(y, mu, log(theta), rep(1, 4), level = 2)
dd <- natural(y, mu, rep(theta, 4))
for (k in 1:4) near(dd[[k]][, 1], -nbd[[c("Dmu", "Dmu2", "Dmu3", "Dmu4")[k]]] / 2)

# Near-Poisson predictor derivatives retain their small, nonzero leading term.
# For y=0, mu=2, log L = -theta log(1 + 2/theta).
f <- nbls(); X <- matrix(c(1, 0, 0, 1), 1, 4)
attr(X, "lpi") <- list(1:2, 3:4)
for (theta in c(1e6, 1e10, 1e14)) {
  z <- f$ll(0, X, c(log(2), 0, 0, log(theta)), 1, f, deriv = 1)
  stopifnot(abs(z$lb[4] * theta + 2) < 1e-4,
            abs(z$lbb[4, 4] * theta - 2) < 1e-4)
}
