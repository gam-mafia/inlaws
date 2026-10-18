library("inlaws")
evaluate <- getFromNamespace(".cgen_gamma_eval", "inlaws")
saturated <- getFromNamespace(".cgen_gamma_saturated", "inlaws")
near <- function(a, b, tol = 1e-5) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))

# Every censoring type, Q=0 and both signs, differentiated through order four.
grid <- expand.grid(sigma = c(.2, .7, 2), Q = c(-.2, 0, 1e-7, .1, 1, 3), kind = 1:4)
eta <- cbind(log(2), log(grid$sigma), log1p(grid$sigma * grid$Q))
y <- rep(c(.3, 1.3, 4), length.out = nrow(eta))
attr(y, "censor") <- cbind(y, -Inf, Inf, y * 1.5)[cbind(seq_along(y), grid$kind)]
r <- evaluate(y, eta, 4L)
stopifnot(all(is.finite(r$ell)))
near(r$ell[, 1], evaluate(y, eta)$ell[, 1], 1e-10)
for (k in 1:3) {
  ep <- em <- eta; ep[, k] <- ep[, k] + 1e-5; em[, k] <- em[, k] - 1e-5
  rp <- evaluate(y, ep, 4L); rm <- evaluate(y, em, 4L)
  for (i in which(rowSums(r$powers) > 0 & r$powers[, k] > 0)) {
    p <- r$powers[i, ]; p[k] <- p[k] - 1L
    ix <- which(rowSums(abs(sweep(r$powers, 2, p))) == 0)
    a <- r$ell[, i] * r$powers[i, k]; b <- (rp$ell[, ix] - rm$ell[, ix]) / 2e-5
    stopifnot(max(abs(a - b) / (1 + abs(a))) < 5e-5)
  }
}
shape <- which(rowSums(abs(sweep(r$powers, 2, c(0, 0, 1)))) == 0)
stopifnot(any(abs(r$ell[grid$Q == 0 & grid$kind != 1, shape]) > .1))

# Independent gamma transformation / log-normal probabilities and densities.
for (i in which(abs(grid$Q) > 1e-6 | grid$Q == 0)) {
  s <- grid$sigma[i]; Q <- grid$Q[i]; k <- Q^-2
  a <- if (Q == 0) log(2) - s^2 / 2 else
    log(2) - lgamma(k + s/Q) + lgamma(k) + (s/Q) * log(k)
  cdf <- function(x, lower = TRUE) {
    if (Q == 0) plnorm(x, a, s, lower.tail = lower, log.p = TRUE) else
      pgamma(exp(log(k) + Q * (log(x) - a)/s), k,
        lower.tail = if (Q > 0) lower else !lower, log.p = TRUE)
  }
  ref <- switch(grid$kind[i],
    if (Q == 0) dlnorm(y[i], a, s, log = TRUE) else {
      logg <- log(k) + Q * (log(y[i]) - a)/s
      dgamma(exp(logg), k, log = TRUE) + logg + log(abs(Q)) - log(s) - log(y[i])
    }, cdf(y[i]), cdf(y[i], FALSE), {
      big <- cdf(y[i] * 1.5); small <- cdf(y[i])
      if (big > -log(2)) { big <- cdf(y[i], FALSE); small <- cdf(y[i] * 1.5, FALSE) }
      big + log(-expm1(small - big))
    })
  near(r$ell[i, 1], ref, 2e-9)
}

# Far tails remain finite on the log scale, including probabilities < exp(-800).
f <- cgen_gammalss()
for (z in c(-40, -8, 8, 40)) for (lower in c(FALSE, TRUE)) {
  value <- exp(z - .5); yy <- structure(value, censor = if (lower) -Inf else Inf)
  rr <- evaluate(yy, matrix(c(0, 0, 0), 1), 4L)
  near(rr$ell[1, 1], pnorm(z, lower.tail = lower, log.p = TRUE), 1e-10)
  stopifnot(all(is.finite(rr$ell)))
}
# Narrow intervals whose probabilities disappear in ordinary CDF subtraction.
for (bound in c(.5, 2, exp(600))) {
  upper <- bound * (1 + 1e-12)
  yy <- structure(bound, censor = upper)
  ee <- matrix(c(log(bound) + .5, 0, 0), 1)
  rr <- evaluate(yy, ee, 4L)
  ref <- dnorm(0, log = TRUE) + log(log1p((upper - bound)/bound))
  near(rr$ell[1, 1], ref, 1e-10)
  stopifnot(all(is.finite(rr$ell)))
}
# Zero lower bounds, full-support intervals and exact-vector equivalence.
mu <- cbind(c(1, 2, 3), .7, c(.5, 1, 2))
yy <- cbind(c(0, 0, 1), c(2, Inf, 1))
rr <- evaluate(yy, log(mu), 4L)
near(rr$ell[1, 1], f$cdf(2, mu[1, , drop = FALSE], logp = TRUE))
stopifnot(all(rr$ell[2, ] == 0))
g <- gen_gammalss()
near(evaluate(1:3, log(mu), 4L)$ell,
  getFromNamespace(".gen_gamma_eval", "inlaws")(1:3, log(mu), 4L)$ell)

# Saturated interval means agree with independent scalar maximization.
for (Q in c(-.3, 0, .7, 2)) {
  s <- .6; mu <- matrix(c(2, s, 1 + s * Q), 1)
  yy <- structure(.8, censor = 3)
  sat <- saturated(yy, mu)$loglik
  opt <- optimize(function(u) evaluate(yy, cbind(u, log(s), log(1 + s * Q)))$ell[1, 1],
    c(-5, 5), maximum = TRUE, tol = 1e-10)
  near(sat, opt$objective, 1e-9)
}

# Super-exponential tails narrower than one ulp of the log boundary.
for (Q in c(-2, 2)) {
  ee <- matrix(c(log(2), log(.2), log1p(.2*Q)), 1)
  a <- getFromNamespace(".gen_gamma_core", "inlaws")(ee)$location[1, 1]
  bound <- exp(a + .2*35/Q)
  yy <- structure(bound, censor = if (Q < 0) -Inf else Inf)
  rr <- evaluate(yy, ee, 4L)
  stopifnot(all(is.finite(rr$ell)))
  near(rr$ell[1, 1], f$cdf(bound, exp(ee), logp = TRUE, lower.tail = Q < 0), 1e-10)
  for (k in 1:3) {
    ep <- em <- ee; ep[, k] <- ep[, k]+1e-5; em[, k] <- em[, k]-1e-5
    rp <- evaluate(yy, ep, 4L); rm <- evaluate(yy, em, 4L)
    for (i in which(rowSums(rr$powers) > 0 & rr$powers[, k] > 0)) {
      p <- rr$powers[i, ]; p[k] <- p[k]-1L
      ix <- which(rowSums(abs(sweep(rr$powers, 2, p))) == 0)
      a <- rr$ell[, i]*rr$powers[i, k]; b <- (rp$ell[, ix]-rm$ell[, ix])/2e-5
      stopifnot(max(abs(a-b)/(1+abs(a))) < 1e-4)
    }
  }
}

# Out-of-range optimizer trials are rejected without distribution errors.
bad <- evaluate(structure(1, censor = -Inf), matrix(c(-1000, 0, 0), 1))
stopifnot(all(is.nan(bad$ell)))

# A unit-probability tail remains usable when its opposite log probability
# is itself beyond double-precision range.
ee <- matrix(c(log(2), log(.03), log(.01)), 1)
rr <- evaluate(structure(.01, censor = Inf), ee, 4L)
stopifnot(all(rr$ell == 0))
