library("inlaws")

near <- function(x, y, tol = 1e-8) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE)))
}
mu <- c(.3, 1.2, 4, 5)
k <- c(.4, 1, 2.6, 5)
q <- c(.001, .5, 3, 8)
lambda <- mu / gamma(1 + 1/k)
p <- c(1e-14, .03, .8, 1 - 1e-10)
f <- weibulls()
for (tail in c(FALSE, TRUE)) for (logp in c(FALSE, TRUE)) {
  pp <- if (logp) log(p) else p
  near(f$cdf(q, cbind(mu, k), logp = logp, lower.tail = tail),
       pweibull(q, k, lambda, lower.tail = tail, log.p = logp))
  near(f$qf(pp, cbind(mu, k), log.p = logp, lower.tail = tail),
       qweibull(pp, k, lambda, lower.tail = tail, log.p = logp))
  for (j in seq_along(k)) {
    g <- weibull(k[j])
    near(g$cdf(q, mu, logp = logp, lower.tail = tail),
         f$cdf(q, cbind(mu, k[j]), logp = logp, lower.tail = tail))
    near(g$qf(pp, mu, log.p = logp, lower.tail = tail),
         f$qf(pp, cbind(mu, k[j]), log.p = logp, lower.tail = tail))
  }
}
set.seed(17)
draw <- f$rd(cbind(mu, k), wt = 8, scale = 2)
set.seed(17)
near(draw, qweibull(runif(length(mu)), k, lambda))
diag <- f$gratia_diagnostics(cbind(mu, k), 1, 1, q, NULL)
stopifnot(diag$continuous)
near(diag$mean, mu)
near(diag$variance, lambda^2 * (gamma(1 + 2/k) - gamma(1 + 1/k)^2))

g <- weibull(1)
near(g$cdf(q, mu), pexp(q, rate = 1/mu))
near(g$qf(p, mu), qexp(p, rate = 1/mu))
near(g$variance(mu), mu^2)
near(g$aic(q, mu, wt = 1, dev = 0), -2 * sum(dexp(q, rate = 1/mu, log = TRUE)))
near(g$dev.resids(q, mu, 1), 2 * (q/mu - 1 - log(q/mu)))
near(g$cdf(c(-Inf, 0, Inf, NA_real_), mu), c(0, 0, 1, NA_real_))
near(g$qf(c(0, 1), 1), c(0, Inf))

# Work below the range of ordinary probabilities, and when gamma(1 + 1/k)
# overflows. A large mean and shape can have representable variance even
# though squared means overflow and squared coefficients of variation underflow.
near(g$cdf(g$qf(-1000, 1e300, log.p = TRUE), 1e300, logp = TRUE), -1000)
near(weibull(.001)$cdf(1, 1, logp = TRUE, lower.tail = FALSE),
     -exp(.001 * lgamma(1001)))
near(weibull(1e250)$variance(1e250), pi^2/6)
stopifnot(identical(weibull(1e-200)$qf(c(0, .9, 1), 1), c(0, 0, Inf)))

g <- weibull(-2)
stopifnot(g$n.theta == 1L)
near(g$getTheta(TRUE), 2)
g$putTheta(log(3))
diag <- g$gratia_diagnostics(mu, 1, 1, q, NULL)
near(diag$cdf(q), weibull(3)$cdf(q, mu))
near(diag$variance, g$variance(mu))
