library(inlaws)

near <- function(x, y, tol = 1e-10) {
  stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = tol)))
}

# The helpers use the family's stored dispersion, not mgcv's external scale
# (which is one). Replication weights must not alter the latent distribution.
mu <- c(.2, 1, 4)
p <- c(.01, .5, .99)
f <- cgamma(theta = .4)
for (phi in c(.4, 1.7)) {
  f$putTheta(log(phi))
  q <- qgamma(p, shape = 1 / phi, scale = mu * phi)
  near(f$qf(p, mu), q)
  near(f$qf(.5, mu, wt = c(0, 2, 8), scale = 99),
       qgamma(.5, shape = 1 / phi, scale = mu * phi))
  near(f$cdf(q, mu), p)
  near(f$cdf(q, mu, wt = c(0, 2, 8), scale = 99, logp = TRUE), log(p))
  near(f$cdf(q, mu, lower.tail = FALSE), 1 - p)
  logp <- c(-1000, -20, -1)
  upper_q <- f$qf(logp, mu, lower.tail = FALSE, log.p = TRUE)
  near(f$cdf(upper_q, mu, lower.tail = FALSE, logp = TRUE), logp)
  set.seed(521)
  expected <- rgamma(length(mu), shape = 1 / phi, scale = mu * phi)
  set.seed(521)
  near(f$rd(mu, wt = c(0, 2, 8), scale = 99), expected)
}
near(f$cdf(c(-1, 0, Inf), mu), c(0, 0, 1))
near(f$qf(c(0, 1), 1), c(0, Inf))
stopifnot(length(f$rd(numeric())) == 0L,
          is.na(f$cdf(NA_real_, 1)), is.na(f$qf(NA_real_, 1)))
# The usual fixing helpers must leave the supplied components in place.
stopifnot(identical(mgcv::fix.family.rd(f)$rd, f$rd),
          identical(mgcv::fix.family.qf(f)$qf, f$qf))

# Distribution callbacks retain the fitted censored model's dispersion.
set.seed(76)
dat <- data.frame(x = runif(150))
y <- rgamma(150, shape = 2, scale = exp(.2 + dat$x) / 2)
dat$yc <- cbind(y, y)
censored <- y < .5
dat$yc[censored, 1] <- .5
dat$yc[censored, 2] <- -Inf
fit <- mgcv::gam(yc ~ s(x, k = 6), data = dat, family = cgamma(), method = "REML")
phi <- fit$family$getTheta(TRUE)
mu <- fitted(fit)
near(fit$family$qf(.5, mu, fit$prior.weights, fit$sig2),
     qgamma(.5, 1 / phi, scale = mu * phi))
set.seed(901)
expected <- replicate(3, rgamma(length(mu), 1 / phi, scale = mu * phi))
set.seed(901)
near(replicate(3, fit$family$rd(mu)), expected)
set.seed(902)
expected <- replicate(3, qgamma(runif(length(mu)), 1 / phi, scale = mu * phi))
set.seed(902)
near(replicate(3, fit$family$qf(runif(length(mu)), mu)), expected)
probs <- c(.05, .5, .95)
near(sapply(probs, function(p) fit$family$qf(p, mu)), sapply(probs, function(p) qgamma(p, 1 / phi, scale = mu * phi)))
nd <- data.frame(x = c(.1, .5, .9))
new_mu <- predict(fit, nd, type = "response")
near(sapply(probs, function(p) fit$family$qf(p, new_mu)),
     sapply(probs, function(p) qgamma(p, 1 / phi, scale = new_mu * phi)))
# Fitted family state must survive serialization without reverting to the
# starting dispersion or depending on an external model object.
restored <- unserialize(serialize(fit$family, NULL))
near(restored$cdf(y, mu), pgamma(y, 1 / phi, scale = mu * phi))

# A censored PIT randomizes over its probability interval, rather than treating
# the supplied censoring limit as an exact observation. Its inverse stays in
# the corresponding response interval (left, right, finite, and exact).
lo <- c(0, 3, .7, 1)
hi <- c(.5, Inf, 1.4, 1)
m <- rep(1, 4)
prob_lo <- fit$family$cdf(lo, m)
prob_hi <- fit$family$cdf(hi, m)
u <- prob_lo + c(.2, .4, .6, .8) * (prob_hi - prob_lo)
y_draw <- fit$family$qf(u, m)
# CDF/quantile inversion can round to either side of an exact observation.
exact <- lo == hi
stopifnot(all(y_draw[!exact] >= lo[!exact]),
          all(y_draw[!exact] <= hi[!exact]),
          all(is.finite(qnorm(u))))
near(y_draw[exact], lo[exact])
