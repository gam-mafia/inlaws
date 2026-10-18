library("inlaws")

near <- function(x, y, tol = 2e-6) {
  stopifnot(
    length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
    max(abs(x - y)) <= tol * max(1, abs(x), abs(y))
  )
}
differentiate <- function(f, x, j, h = 1e-4) {
  e <- numeric(length(x))
  e[j] <- h
  (f(x - 2 * e) - 8 * f(x - e) + 8 * f(x + e) - f(x + 2 * e)) / (12 * h)
}
natural <- getFromNamespace(".betals_derivatives", "inlaws")

# Every predictor derivative is checked against a derivative of the preceding
# order. Include small shapes, near-boundary means, and the Stirling transition.
grid <- expand.grid(
  y = c(1e-7, .1, .5, .9, 1 - 1e-7),
  eta = c(-12, -2, 0, 2, 12),
  phi = c(1e-4, .1, 2, 16, 32, 100, 1e5)
)
eta <- cbind(grid$eta, log(grid$phi))
dd <- natural(grid$y, eta, 4L)
for (r in 1:4) {
  for (q in 0:r) {
    p <- r - q
    coordinate <- if (p) 1L else 2L
    qq <- q - (coordinate == 2L)
    f <- function(z) {
      candidate <- eta
      candidate[, coordinate] <- candidate[, coordinate] + z
      prev <- natural(grid$y, candidate, r - 1L)
      if (r == 1L) prev$l0 else prev[[r]][, qq + 1L]
    }
    fd <- differentiate(f, 0, 1)
    # Rowwise tolerances prevent large derivatives from hiding failures elsewhere.
    stopifnot(all(abs(fd - dd[[r + 1L]][, q + 1L]) <=
      2e-5 * pmax(1, abs(fd), abs(dd[[r + 1L]][, q + 1L]))))
  }
}

# Independent mean score/Hessian and precision derivatives from betar().
y <- c(.02, .2, .6, .98)
mu <- c(.01, .3, .5, .9)
phi <- c(.1, 2, 10, 100)
z <- natural(y, cbind(qlogis(mu), log(phi)))
ref <- mgcv::betar()$Dd(y, mu, log(phi), rep(1, 4), level = 2)
m1 <- mu * (1 - mu)
m2 <- m1 * (1 - 2 * mu)
near(z$l1[, 1], -ref$Dmu * m1 / 2)
near(z$l2[, 1], -(ref$Dmu2 * m1^2 + ref$Dmu * m2) / 2)
near(z$l1[, 2], -ref$Dth / 2)
near(z$l2[, 3], -ref$Dth2 / 2)
near(z$l2[, 2], -ref$Dmuth * m1 / 2)

# At the symmetric mean, precision curvature remains small and nonzero even
# when the leading gamma terms would cancel to floating-point zero.
for (phi in c(1e6, 1e10, 1e14)) {
  z <- natural(.5, matrix(c(0, log(phi)), 1))
  near(z$l1[, 2], .5 + .25 / phi, 1e-12)
  near(z$l2[, 3] * phi, -.25, 1e-10)
  near(z$l3[, 4] * phi, .25, 1e-10)
  near(z$l4[, 5] * phi, -.25, 1e-10)
  near(c(z$l1[, 1], z$l3[, 1]), c(0, 0), 1e-12)
}

# Check link transformations and coefficient-level contractions, with offsets,
# fractional weights, and a zero-weight row. Use separate predictor coefficients.
set.seed(130)
{
  n <- 12L
  X <- cbind(1, seq(-1, 1, length.out = n), rnorm(n), 1, runif(n))
  attr(X, "lpi") <- list(1:3, 4:5)
  beta <- c(.2, .1, -.1, log(10), .2)
  y <- rep(c(.05, .2, .5, .9), 3)
  w <- c(0, .4, rep(1, 5), rep(2, 5))
  off <- list(runif(n, -.2, .2), runif(n, -.4, .4))
  fam <- betals()
  ll <- function(b, deriv = 1L, ...) {
    fam$ll(y, X, b, w, fam,
      offset = off, deriv = deriv, ...
    )
  }
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
  for (j in 1:2) {
    for (k in j:2) {
      kk <- kk + 1L
      # Differentiate the third derivative contraction and add the change in the
      # coefficient direction induced by second derivatives of beta(rho).
      fd <- differentiate(function(a) {
        ll(beta + a * directions[, k], 3,
          d1b = directions
        )$d1H[[j]]
      }, 0, 1, h = 1e-4)
      fd <- fd + differentiate(function(a) ll(beta + a * second[, kk])$lbb,
        0, 1,
        h = 1e-4
      )
      near(full$trHid2H[kk], sum(V * fd), tol = 1e-5)
    }
  }
  # Sandwich is the cross-product of weighted individual coefficient scores.
  scores <- t(vapply(seq_len(n), function(i) {
    wi <- numeric(n)
    wi[i] <- w[i]
    fam$ll(y, X, beta, wi, fam, offset = off, deriv = 1)$lb
  }, numeric(5)))
  near(fam$sandwich(y, X, beta, w, fam, off), crossprod(scores))
  # Omitting a zero-weight row must give exactly the same likelihood.
  X1 <- X[-1, , drop = FALSE]
  attr(X1, "lpi") <- attr(X, "lpi")
  omitted <- fam$ll(y[-1], X1, beta, w[-1], fam,
    offset = lapply(off, `[`, -1), deriv = 1
  )
  near(first$l, omitted$l)
  near(first$lb, omitted$lb)
  near(first$lbb, omitted$lbb)
  # Explicit eta and coefficient construction agree, including offsets.
  jj <- attr(X, "lpi")
  eta <- sapply(1:2, function(j) {
    drop(X[, jj[[j]], drop = FALSE] %*%
      beta[jj[[j]]]) + off[[j]]
  })
  near(first$l, fam$ll(y, X, beta, w, fam, eta = eta)$l)
  eta[1, ] <- c(Inf, -Inf)
  near(first$l, fam$ll(y, X, beta, w, fam, eta = eta)$l)
}
