library(inlaws)

near <- function(a, b, tol = 1e-12) {
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
}

fam <- clognormal(theta = 0.7)
mu <- c(0.5, 1, 3)
p <- c(0.01, 0.4, 0.99)
q <- c(0.2, 1, 5)

# Both tail directions and both probability scales, including vector means
# and the fitted family's updated sigma rather than its initial value.
for (sigma in c(0.7, 1.2)) {
  fam$putTheta(log(sigma))
  meanlog <- log(mu) - sigma^2/2
  for (lower in c(TRUE, FALSE)) {
    for (logged in c(FALSE, TRUE)) {
      prob <- if (logged) log(p) else p
      quant <- fam$qf(prob, mu, lower.tail = lower, log.p = logged)
      near(quant, qlnorm(prob, meanlog, sigma, lower.tail = lower, log.p = logged))
      near(fam$cdf(q, mu, logp = logged, lower.tail = lower),
           plnorm(q, meanlog, sigma, lower.tail = lower, log.p = logged))
      near(fam$cdf(quant, mu, logp = logged, lower.tail = lower), prob)
    }
  }
  # Defaults and the fifth positional CDF argument retain the established
  # CDF interface. Case weights and external scale do not affect the law.
  near(fam$qf(p, mu), qlnorm(p, meanlog, sigma))
  near(fam$cdf(q, mu), plnorm(q, meanlog, sigma))
  near(fam$cdf(q, mu, 1, 1, TRUE), plnorm(q, meanlog, sigma, log.p = TRUE))
  near(fam$qf(p, mu, wt = c(0, 2, 5), scale = 4), fam$qf(p, mu))
  near(fam$cdf(q, mu, wt = c(0, 2, 5), scale = 4), fam$cdf(q, mu))
}

# Unit log-SD and mu = exp(1/2) give a standard normal log response.
fam <- clognormal(theta = 1)
mu <- exp(0.5)
stopifnot(1 - fam$cdf(exp(10), mu) == 0)
upper <- fam$cdf(exp(10), mu, lower.tail = FALSE)
stopifnot(upper > 0)
near(upper / pnorm(10, lower.tail = FALSE), 1)

stopifnot(fam$cdf(exp(-40), mu) == 0)
log_lower <- fam$cdf(exp(-40), mu, logp = TRUE)
stopifnot(is.finite(log_lower))
near(log_lower, pnorm(-40, log.p = TRUE))
log_upper <- fam$cdf(exp(40), mu, lower.tail = FALSE, logp = TRUE)
near(log_upper, log_lower)

# Avoid rounding a near-one input to one or underflowing a tiny input to zero.
stopifnot(is.infinite(fam$qf(1 - 1e-20, mu)))
upper_quantile <- fam$qf(1e-20, mu, lower.tail = FALSE)
stopifnot(is.finite(upper_quantile))
near(fam$cdf(upper_quantile, mu, lower.tail = FALSE) / 1e-20, 1)
for (lower in c(TRUE, FALSE)) {
  quant <- fam$qf(-1000, mu, lower.tail = lower, log.p = TRUE)
  stopifnot(is.finite(quant), quant > 0)
  near(fam$cdf(quant, mu, lower.tail = lower, logp = TRUE), -1000)
  endpoint <- if (lower) c(0, Inf) else c(Inf, 0)
  near(fam$qf(c(0, 1), mu, lower.tail = lower), endpoint)
  near(fam$qf(c(-Inf, 0), mu, lower.tail = lower, log.p = TRUE), endpoint)
  probability <- if (lower) c(0, 1) else c(1, 0)
  near(fam$cdf(c(0, Inf), mu, lower.tail = lower), probability)
  near(fam$cdf(c(0, Inf), mu, lower.tail = lower, logp = TRUE), log(probability))
}

# mgcv must preserve the extended quantile callback interface.
near(mgcv::fix.family.qf(fam)$qf(1e-20, mu, lower.tail = FALSE), upper_quantile)
