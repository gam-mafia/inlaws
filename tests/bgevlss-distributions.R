library("inlaws")

near <- function(a, b, tol = 2e-12) stopifnot(all(abs(a - b) <= tol*(1 + abs(b))))
v <- function(p) -log(-log(p))
f <- bgevlss()
p <- c(.0001, .02, .05, .0500001, .1, .1999999, .2, .25, .5, .75, .99, 1 - 1e-10)
mu <- cbind(2, 1.4, rep(.15, length(p)))
q <- f$qf(p, mu)
near(f$cdf(q, mu), p)
near(f$cdf(q, mu, lower.tail = FALSE), 1 - p)
near(f$qf(log(p), mu, log.p = TRUE), q)
near(f$qf(log1p(-p), mu, lower.tail = FALSE, log.p = TRUE), q)
near(q[p == .5], 2)
near(q[p == .75] - q[p == .25], 1.4)
stopifnot(identical(f$qf(c(0, 1), mu[1:2, ]), c(-Inf, Inf)),
          identical(f$qf(c(0, 1), mu[1:2, ], lower.tail = FALSE), c(Inf, -Inf)),
          identical(f$cdf(c(-Inf, Inf), mu[1:2, ], log_p = TRUE), c(-Inf, 0)),
          all(is.na(f$cdf(c(NA, 0), rbind(mu[1, ], c(0, NA, 0))))),
          is.nan(suppressWarnings(f$qf(-.1, mu[1, , drop = FALSE]))),
          inherits(try(bgevlss(pa = .2, pb = .2), silent = TRUE), "try-error"),
          inherits(try(bgevlss(pb = .3), silent = TRUE), "try-error"))

# The extreme upper tail includes values beyond lower-CDF rounding and
# survival underflow. A very small spread also tests right-tail scaling.
lp <- c(-1e-14, -.01, -20, -40, -200, -800, -1200)
mu_tail <- cbind(0, 1e-50, rep(.15, length(lp)))
x <- f$qf(lp, mu_tail, log.p = TRUE, lower.tail = FALSE)
stopifnot(all(is.finite(x)))
near(f$cdf(x, mu_tail, log_p = TRUE, lower.tail = FALSE), lp)
near(f$cdf(f$qf(-1200, mu[1, , drop = FALSE], log.p = TRUE),
           mu[1, , drop = FALSE], log_p = TRUE), -1200)
stopifnot(is.finite(f$cdf(1e308, matrix(c(0, 1e-50, .15), 1),
                          log_p = TRUE, lower.tail = FALSE)))

# Reference agreement at default and custom probabilities, including points
# in both tails and the blend. evgam is optional outside these tests.
for (args in list(list(), list(pa = .04, pb = .16, alpha = .7, beta = .6))) {
  f0 <- do.call(bgevlss, args)
  pars <- cbind(rep(2, length(p)), 1.4, rep(c(.05, .18, .4), length.out = length(p)))
  q0 <- f0$qf(p, pars)
  near(f0$cdf(q0, pars), p)
  a <- if (length(args)) args$alpha else .5
  b <- if (length(args)) args$beta else .5
  near(f0$qf(a, pars), pars[, 1])
  near(f0$qf(1-b/2, pars) - f0$qf(b/2, pars), pars[, 2])
  if (requireNamespace("evgam", quietly = TRUE)) {
    params <- c(list(location = pars[, 1], scale = pars[, 2], shape = pars[, 3]), args)
    near(do.call(evgam::pbgev, c(list(q = q0), params)), p)
    idx <- p >= .02 & p <= .99
    ref <- do.call(evgam::dbgev, c(list(x = q0[idx], log = TRUE),
      lapply(params, function(z) if (length(z) > 1L) z[idx] else z)))
    eta <- cbind(pars[idx, 1L], log(pars[idx, 2:3, drop = FALSE]))
    near(f0$ll(q0[idx], NULL, NULL, 1, f0, eta = eta)$l0, ref, 1e-11)
  }
}

# At zero shape the response distribution is exactly Gumbel.
gumbel <- mu
gumbel[, 3L] <- 0
sigma <- 1.4/(v(.75) - v(.25))
near(f$qf(p, gumbel), 2 + sigma*(v(p) - v(.5)))
near(f$cdf(q, gumbel), exp(-exp(-(q - 2)/sigma - v(.5))))
dm <- f$gratia_diagnostics(gumbel[1, , drop = FALSE])
near(dm$mean, 2 + sigma*(-digamma(1) - v(.5)))
near(dm$variance, pi^2/6*sigma^2)
tiny <- gumbel
tiny[, 3L] <- 1e-12
near(f$qf(p, tiny), f$qf(p, gumbel), 2e-11)

subnormal <- gumbel[1, , drop = FALSE]
subnormal[, 3L] <- .Machine$double.xmin*.01
subnormal_moments <- f$gratia_diagnostics(subnormal)
near(subnormal_moments$mean, dm$mean)
near(subnormal_moments$variance, dm$variance)

# Check the moments independently by integration of the density, split
# at the two joins. Use nonzero location and nonunit spread.
par <- mu[1, , drop = FALSE]
diagnostics <- f$gratia_diagnostics(par)
cuts <- c(-Inf, f$qf(.05, par), f$qf(.2, par), Inf)
eta <- c(2, log(1.4), log(.15))
for (k in 0:2) {
  center <- if (k == 2L) diagnostics$mean else 0
  integral <- sum(vapply(1:3, function(i) integrate(function(x) {
    e <- matrix(rep(eta, each = length(x)), ncol = 3L)
    (x - center)^k*exp(f$ll(x, NULL, NULL, 1, f, eta = e)$l0)
  }, cuts[i], cuts[i+1L], rel.tol = 1e-8)$value, numeric(1)))
  near(integral, c(1, diagnostics$mean, diagnostics$variance)[k+1L], 2e-8)
}
moments <- f$gratia_diagnostics(cbind(0, 1, c(.49, .5 - 1e-10, .5, .99, 1 - 1e-10, 1, 1.2)))
stopifnot(all(is.finite(moments$mean[1:5])), all(is.na(moments$mean[6:7])),
          all(is.finite(moments$variance[1:2])), all(is.infinite(moments$variance[3:7])),
          moments$variance[2] > 1e8, moments$mean[5] > 1e8,
          !moments$mean_available, !moments$variance_available,
          all(is.finite(moments$cdf(0))), all(is.finite(moments$qf(.9))),
          all(is.na(moments$diagnostic_residuals(rep(0, 7), "pearson")[3:7])))
set.seed(481)
par <- par[rep(1, 12000), , drop = FALSE]
d <- f$gratia_diagnostics(par)
r <- d$rd()
u <- d$cdf(r)
stopifnot(abs(mean(r) - diagnostics$mean) < .04,
          abs(var(r) - diagnostics$variance) < .17,
          max(abs(sort(u) - ppoints(length(u)))) < .02)
