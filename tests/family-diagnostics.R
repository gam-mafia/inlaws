library("inlaws")

near <- function(x, y, tol = 1e-7) stopifnot(isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE)))
diagnostic <- function(f, mu, y = NULL, eta = NULL) f$gratia_diagnostics(mu, rep(1, NROW(mu)), 1, y, eta)

# Constructor factories use the same mutable parameter state as native hooks.
cases <- list(nbls = list(nbls(), cbind(c(.5, 2, 4), 2)),
  cmp = list(cmp(1.2), c(.5, 2, 4)), cmpls = list(cmpls(), cbind(c(.5, 2, 4), 1.2)),
  gpoisson = list(gpoisson(2), c(.5, 2, 4)), pig = list(pig(.5), c(.5, 2, 4)),
  pigls = list(pigls(), cbind(c(.5, 2, 4), .5)),
  zinb = list(zinb(), cbind(c(.5, 2, 4), .7, 2)), zanb = list(zanb(), cbind(c(.5, 2, 4), .7, 2)))
for (z in cases) {
  f <- z[[1]]; mu <- z[[2]]
  stopifnot(is.function(f$gratia_diagnostics))
  d <- diagnostic(f, mu)
  near(d$cdf_left(c(0, .5, 2)), f$cdf(c(-1, 0, 1), mu))
  near(exp(d$cdf_left(c(0, .5, 2), TRUE, FALSE)), 1 - f$cdf(c(-1, 0, 1), mu))
  clone <- unserialize(serialize(f, NULL))
  near(diagnostic(clone, mu)$variance, d$variance)
  if (is.function(clone$putTheta)) {
    clone$putTheta(clone$getTheta() + .2)
    near(diagnostic(clone, mu)$variance, clone$variance(mu))
    stopifnot(!isTRUE(all.equal(diagnostic(clone, mu)$variance, d$variance)))
  }
}

# CMP moments against a separate truncated power series and root solve.
mu <- c(.5, 2, 4); nu <- 1.2; k <- 0:300
reference <- lapply(mu, function(m) {
  probabilities <- function(a) { z <- a * k - nu * lgamma(k + 1); p <- exp(z - max(z)); p / sum(p) }
  a <- uniroot(function(a) sum(k * probabilities(a)) - m, c(-10, 10), tol = 1e-12)$root
  p <- probabilities(a)
  c(mean = sum(k * p), variance = sum((k - m)^2 * p))
})
for (f in list(cmp(nu), cmpls())) {
  d <- diagnostic(f, if (inherits(f, "general.family")) cbind(mu, nu) else mu)
  near(d$mean, vapply(reference, `[[`, numeric(1), "mean"))
  near(d$variance, vapply(reference, `[[`, numeric(1), "variance"))
}

# Beta-binomial margins against a separate PMF, including zero trial totals.
f <- betabinomial(); mu <- cbind(c(.2, .5, .8), c(3, 4, 5)); trials <- c(0, 7, 12)
y <- cbind(c(0, 3, 12), trials - c(0, 3, 12))
d <- diagnostic(f, mu, y)
for (i in 1:3) {
  k <- 0:trials[i]; a <- mu[i, 1] * mu[i, 2]; b <- (1 - mu[i, 1]) * mu[i, 2]
  pmf <- choose(trials[i], k) * beta(k + a, trials[i] - k + b) / beta(a, b)
  q <- y[i, 1]
  near(d$cdf(y[, 1])[i], sum(pmf[k <= q]))
  near(d$cdf_left(y[, 1])[i], sum(pmf[k < q]))
  near(d$mean[i], sum(k * pmf))
  near(d$variance[i], sum((k - sum(k * pmf))^2 * pmf))
}
stopifnot(all(d$qf(c(0, .5, 1)) <= trials))
near(d$cdf(c(-1, 2.5, Inf)) + d$cdf(c(-1, 2.5, Inf), lower_tail = FALSE), rep(1, 3))

# Unknown observation processes retain intervals without guessing thresholds.
y <- structure(c(.4, 1, 2), censor = c(-Inf, 1, Inf))
for (z in list(list(cgamma(.6), rep(2, 3)), list(cgammals(), cbind(rep(2, 3), .6)),
              list(clognormal(.5), rep(2, 3)), list(clognormalls(), cbind(rep(2, 3), .5)))) {
  d <- diagnostic(z[[1]], z[[2]], y)
  near(d$bounds(y)$lower, c(-Inf, 1, 2)); near(d$bounds(y)$upper, c(.4, 1, Inf))
  stopifnot(!d$reference_available, !diagnostic(z[[1]], z[[2]], 1:3)$reference_available)
}
f <- cgamma(.6, observe = censoring_process(.4, 3))
f <- unserialize(serialize(f, NULL))
observed <- diagnostic(f, rep(2, 3), y)$observe(c(.1, 1, 4))
near(observed, c(.4, 1, 3)); near(attr(observed, "censor"), c(-Inf, 1, Inf))
stopifnot(inherits(try(censoring_process(c(.1, .2)), silent = TRUE), "try-error"))

# Full endpoint moments from independently integrated beta components.
f <- ordbeta(); mu <- c(.2, .5, .8); th <- f$getTheta(TRUE)
p0 <- plogis(th[2] - qlogis(mu)); p1 <- plogis(qlogis(mu) - th[3]); interior <- 1 - p0 - p1
first <- second <- numeric(3)
for (i in 1:3) {
  density <- function(x) dbeta(x, mu[i] * th[1], (1 - mu[i]) * th[1])
  first[i] <- p1[i] + interior[i] * integrate(function(x) x * density(x), 0, 1)$value
  second[i] <- p1[i] + interior[i] * integrate(function(x) x^2 * density(x), 0, 1)$value
}
d <- diagnostic(f, mu)
near(d$mean, first); near(d$variance, second - first^2)
near(d$cdf_left(c(0, .5, 1)), c(0, f$cdf(.5, mu)[2], 1 - p1[3]))

# Heavy-tail moment existence is separate from CDF availability.
d <- diagnostic(gpd(), cbind(2, c(.2, .6, 1.2)))
stopifnot(!d$mean_available, !d$variance_available, is.na(d$mean[3]), is.infinite(d$variance[2]))
near(d$variance[1], 4 / (.8^2 * .6))
for (f in list(log_logistic(.6), log_logisticls())) {
  d <- diagnostic(f, if (inherits(f, "general.family")) cbind(1:3, .6) else 1:3)
  stopifnot(d$continuous, !d$variance_available, all(is.infinite(d$variance)))
}

# Multivariate means, beta margins, and joint draws on the correct response scale.
eta <- rbind(c(.5, -.7, log(4)), c(-.2, .3, log(3)))
p <- exp(cbind(0, eta[, 1:2])); p <- p / rowSums(p); phi <- exp(eta[, 3])
y <- rbind(c(2, 3, 5), c(4, 0, 3))
for (f in list(dirichlet(2), dirmult(2))) {
  is_count <- f$family == "dirmult"
  d <- diagnostic(f, eta, if (is_count) y else y / rowSums(y), eta)
  total <- if (is_count) rowSums(y) else 1
  near(d$mean, p * total)
  near(d$variance, p * (1 - p) * if (is_count) total * (total + phi) / (1 + phi) else 1 / (1 + phi))
  stopifnot(!d$independent, d$joint_available, is.null(d$qf))
  draw <- d$rd(); near(rowSums(draw), rep_len(total, 2))
  if (!is_count) {
    near(d$cdf(y / rowSums(y)), pbeta(y / rowSums(y), p * phi, (1 - p) * phi))
    stopifnot(inherits(try(d$diagnostic_residuals(y, "deviance"), silent = TRUE), "try-error"))
  }
}

f <- cumulative_link(R = 4, link = "probit")
eta <- c(-1, 0, 1); y <- c(1, 2, 4)
d <- diagnostic(f, rep(99, 3), y, eta)
prob <- sapply(1:4, function(k) f$cdf(k, eta) - f$cdf(k - 1, eta))
near(d$mean, as.vector(prob %*% (1:4)))
near(d$variance, as.vector(prob %*% (1:4)^2) - d$mean^2)
near(d$diagnostic_residuals(y, "response"), y - d$mean)
near(d$diagnostic_residuals(y, "pearson"), (y - d$mean) / sqrt(d$variance))

# A single observation must keep the matrix shape of a joint response draw.
d <- diagnostic(dirmult(2), eta = matrix(c(.2, -.4, log(3)), 1),
  mu = matrix(c(.2, -.4, log(3)), 1), y = matrix(c(0, 1, 4), 1))
stopifnot(identical(dim(d$rd()), c(1L, 3L)), identical(dim(d$cdf(1)), c(1L, 3L)))

# Sample covariance checks the joint generator against the analytic covariance.
set.seed(735)
for (f in list(dirichlet(2), dirmult(2))) {
  pars <- matrix(c(.2, -.4, log(3)), 1)
  d <- diagnostic(f, pars, matrix(c(2, 3, 5), 1), pars)
  draws <- t(replicate(5000, as.numeric(d$rd())))
  p <- exp(c(0, .2, -.4)); p <- p / sum(p)
  factor <- if (f$family == "dirmult") 10 * (10 + 3) / (1 + 3) else 1 / (1 + 3)
  expected <- (diag(p) - tcrossprod(p)) * factor
  near(colMeans(draws), as.numeric(d$mean), .06)
  near(cov(draws), expected, .08)
}
