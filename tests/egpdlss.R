library("inlaws")

near <- function(x, y, tolerance = 1e-7) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)))
}
expect_warning <- function(expr) {
  warned <- FALSE
  out <- withCallingHandlers(expr, warning = function(w) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  })
  stopifnot(warned)
  out
}
expect_error <- function(expr, pattern) {
  out <- tryCatch(expr, error = identity)
  stopifnot(inherits(out, "error"), grepl(pattern, conditionMessage(out)))
}

f <- egpdlss()
g <- gpdlss()
stopifnot(identical(gpd, gpdlss), identical(g$family, "gpd"), f$nlp == 3L)
expect_error(egpdlss(c("log", "log", "log")), "requires links")
expect_error(egpdlss(c("log", "identity")), "requires links")

# Compare the distribution and log likelihood with the original family
# across negative, near-zero and positive shape values.
mu <- cbind(c(.5, 1, 2, 3, 4, 5), c(-.4, -1e-10, 0, 1e-10, .2, 1), 1)
eta <- cbind(log(mu[, 1]), log(mu[, 2] + .5), 0)
for (p in c(0, 1e-12, .01, .4, .95, 1 - 1e-12, 1)) {
  q <- f$qf(p, mu)
  near(q, g$qf(p, mu[, 1:2]))
  near(f$cdf(q, mu), g$cdf(q, mu[, 1:2]))
  near(f$cdf(q, mu, logp = TRUE), g$cdf(q, mu[, 1:2], logp = TRUE))
}
y <- f$qf(.7, mu)
wt <- c(0, 1, 2, 1, .3, 3)
dg <- g$ll(y, NULL, NULL, wt, g, eta = eta[, 1:2], deriv = 4, ncv = TRUE)
de <- f$ll(y, NULL, NULL, wt, f, eta = eta, deriv = 4, ncv = TRUE)
near(de$l0, dg$l0)
for (order in 1:4) {
  # Select the entries with no derivative in the bulk predictor.
  pure <- 1L + c(0L, cumsum(seq_len(order)))
  near(de[[paste0("l", order)]][, pure], dg[[paste0("l", order)]])
}

# An independent log density, and finite differences of packed derivatives
# in the three actual predictors through fourth order for both tail links.
logdensity <- function(y, sigma, xi, kappa) {
  t <- y / sigma
  h <- t
  nz <- xi != 0
  h[nz] <- log1p(xi[nz] * t[nz]) / xi[nz]
  log(kappa) - log(sigma) - (1 + xi) * h + (kappa - 1) * log(-expm1(-h))
}
packed <- function(d) {
  if (!d) return(matrix(integer(), 1, 0))
  p <- as.matrix(expand.grid(rep(list(1:3), d)))
  p <- p[apply(p, 1, function(x) all(diff(x) >= 0)), , drop = FALSE]
  p[do.call(base::order, as.data.frame(p)), , drop = FALSE]
}
fd <- function(f, theta, j, h = 2e-5) {
  step <- numeric(3); step[j] <- h
  (f(theta - 2 * step) - 8 * f(theta - step) +
     8 * f(theta + step) - f(theta + 2 * step)) / (12 * h)
}
for (link in c("logshift", "identity")) {
  f <- egpdlss(c("log", link, "log"))
  for (shape in c(-.45, -1e-8, 0, .4)) for (kappa in c(.45, 1, 2.4)) {
    y <- c(.0001, .2, .9, 2.4)
    theta <- c(log(1.2), f$linfo[[2]]$linkfun(shape), log(kappa))
    l <- function(theta, d) {
      pars <- matrix(rep(theta, each = length(y)), ncol = 3L)
      code <- c(0L, 1L, 1L, 3L, 4L)[d + 1L]
      f$ll(y, NULL, NULL, NULL, f, eta = pars, deriv = code, ncv = TRUE)[[paste0("l", d)]]
    }
    near(l(theta, 0L), logdensity(y, 1.2, rep(shape, length(y)), kappa))
    for (d in 1:4) {
      actual <- l(theta, d)
      previous <- packed(d - 1L)
      target <- packed(d)
      for (j in 1:3) {
        cols <- apply(cbind(previous, j), 1, function(v) {
          which(apply(target, 1, function(row) all(row == sort(v))))
        })
        near(fd(function(p) l(p, d - 1L), theta, j), actual[, cols], 8e-5)
      }
    }
  }
}

f <- egpdlss()
mu <- cbind(c(.6, 2, 4, 1), c(-.3, 0, .2, .4), c(.6, 1.6, 3, 1.1))
for (p in c(.001, .15, .6, .95)) {
  q <- f$qf(p, mu)
  near(f$cdf(q, mu), rep(p, nrow(mu)))
  near(g$cdf(q, mu[, 1:2])^mu[, 3], rep(p, nrow(mu)))
  near(f$qf(log1p(-p), mu, log.p = TRUE, lower.tail = FALSE), q)
}
near(f$qf(0, mu), rep(0, nrow(mu)))
near(f$qf(1, mu), c(2, Inf, Inf, Inf))
near(f$cdf(c(-Inf, 0, Inf, Inf), mu), c(0, 0, 1, 1))
stopifnot(is.nan(expect_warning(f$qf(-1, mu))[1]),
          is.nan(expect_warning(f$cdf(1, cbind(1, 0, 0)))),
          is.na(f$qf(NA_real_, cbind(1, 0, 1))),
          is.na(f$cdf(1, cbind(1, 0, NA_real_))))
set.seed(144); u <- runif(nrow(mu))
set.seed(144); near(f$rd(mu, wt = 2, scale = 3), f$qf(u, mu))

# Large kappa preserves a representable log upper tail after the GPD
# survival underflows. Tiny lower quantiles retain the scale factor.
m <- cbind(1, 0, 1e300)
near(f$cdf(1000, m, logp = TRUE, lower.tail = FALSE), log(1e300) - 1000)
near(f$qf(log(1e300) - 1000, m, log.p = TRUE, lower.tail = FALSE), 1000)
m <- cbind(1e100, 0, 1)
near(log(f$qf(-900, m, log.p = TRUE)), log(1e100) - 900)
near(f$cdf(f$qf(-900, m, log.p = TRUE), m, logp = TRUE), -900)
m <- cbind(1, 0, .1)
near(f$cdf(1e-300, m, logp = TRUE), .1 * log(1e-300))
m <- cbind(2, -.4, 2)
x <- 5 - 1e-7
near(f$cdf(x, m, logp = TRUE, lower.tail = FALSE), log(2) + 2.5 * log((5 - x)/5))
stopifnot(f$cdf(5, m) == 1, f$cdf(5 + 1e-7, m, lower.tail = FALSE) == 0)
et <- matrix(c(log(2), log(.1), log(2)), nrow = 1)
stopifnot(is.finite(f$ll(x, NULL, NULL, 1, f, eta = et)$l),
          f$ll(5 + 1e-7, NULL, NULL, 1, f, eta = et)$l == -Inf,
          f$ll(5 + 1e-7, NULL, NULL, 0, f, eta = et)$l == 0)

# Check moments by numerical integration over the exponential variate.
mu <- cbind(2, c(-.4, -1e-8, 0, 1e-8, .2, .51, 1.1), 1.7)
di <- f$gratia_diagnostics(mu, 1, 1, NULL, NULL)
stopifnot(!di$mean_available, !di$variance_available, is.na(di$mean[7]),
          all(is.infinite(di$variance[6:7])), di$continuous)
near(di$qf(di$cdf(1, log_p = TRUE, lower_tail = FALSE), log_p = TRUE, lower_tail = FALSE),
     rep(1, nrow(mu)))
for (i in 1:5) {
  xi <- mu[i, 2]
  # G(y) = 1 - exp(-h); derivative of G(y)^kappa with respect to h.
  y <- function(h) if (xi == 0) 2 * h else 2 * expm1(xi*h)/xi
  reference <- function(k) integrate(function(h) {
    y(h)^k * 1.7 * exp(-h) * (-expm1(-h))^.7
  }, 0, 100, rel.tol = 1e-9)$value
  m1 <- reference(1)
  near(di$mean[i], m1)
  near(di$variance[i], reference(2) - m1^2)
}
o <- list(y = rep(1, nrow(mu)), fitted.values = mu, prior.weights = c(rep(1, 6), 0))
near(f$residuals(o, "response"), 1 - di$mean)
rp <- expect_warning(f$residuals(o, "pearson"))
stopifnot(is.na(rp[6]), rp[7] == 0)
expect_error(f$residuals(o, "deviance"), "unavailable")
expect_error(gam(list(y ~ 1, ~ 1, ~ 1), data = data.frame(y = c(0, 1, 2)),
                 family = f), "strictly positive")

# Small kappa limits of the moments follow from the derivative of log Beta.
k <- 1e-20
d <- f$gratia_diagnostics(cbind(2, .2, k), 1, 1, NULL, NULL)
first <- digamma(1) - digamma(.8)
second <- digamma(1) - digamma(.6)
near(d$mean / k, 2 * first / .2)
near(d$variance / k, 4 * (second - 2*first) / .2^2)

# An endpoint smaller than double range must not replace the lower bound.
stopifnot(f$cdf(0, cbind(1e-300, -1e100, .8)) == 0,
          g$cdf(0, cbind(1e-300, -1e100)) == 0)
