library(inlaws)
near <- function(x, y, tol = 1e-7) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}
expect_error <- function(expr, pattern) {
  err <- tryCatch(expr, error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err), fixed = TRUE))
}
expect_warning <- function(expr, pattern) {
  seen <- FALSE
  ans <- withCallingHandlers(expr, warning = function(w) {
    if (grepl(pattern, conditionMessage(w), fixed = TRUE)) seen <<- TRUE
    invokeRestart("muffleWarning")
  })
  stopifnot(seen)
  ans
}
D <- getFromNamespace(".gpd_derivatives", "inlaws")

# Five-point numerical derivatives use a step chosen independently of the
# implementation. Differentiating each derivative order checks ALL mixed
# entries through order four, including the exponential limit and series seam.
fd <- function(f, x, j, h = 1e-4) {
  e <- numeric(length(x)); e[j] <- h
  (f(x - 2 * e) - 8 * f(x - e) + 8 * f(x + e) - f(x + 2 * e)) / (12 * h)
}
for (xi in c(-0.49, -0.4, -1e-8, 0, 1e-8, 0.2, 0.5, 1.5)) {
  y <- c(0, 0.2, 1, 2)
  theta <- c(0.3, xi)
  for (ord in 1:4) {
    f <- function(p) as.numeric(D(y, rep(p[1], 4), rep(p[2], 4), ord - 1)[[ord]])
    actual <- D(y, rep(theta[1], 4), rep(theta[2], 4), ord)[[ord + 1L]]
    near(matrix(fd(f, theta, 1), nrow = 4), actual[, 1:ord], 2e-6)
    near(matrix(fd(f, theta, 2), nrow = 4), actual[, 2:(ord + 1L)], 2e-6)
  }
}
for (v in c(-0.50001, -0.5, -0.49999, 0.49999, 0.5, 0.50001, -0.95)) {
  theta <- c(0, v)
  for (ord in 1:4) {
    f <- function(p) as.numeric(D(1, p[1], p[2], ord - 1)[[ord]])
    actual <- D(1, theta[1], theta[2], ord)[[ord + 1L]]
    near(fd(f, theta, 1, 1e-5), actual[, 1:ord], 1e-5)
    near(fd(f, theta, 2, 1e-5), actual[, 2:(ord + 1L)], 1e-5)
  }
}

# Independent likelihood reference, also used below for optim-based MLEs.
logdensity <- function(y, sigma, xi) {
  if (any(sigma <= 0 | 1 + xi * y / sigma <= 0)) return(rep(-Inf, length(y)))
  if (length(xi) == 1 && xi == 0) return(-log(sigma) - y / sigma)
  -log(sigma) - (1 + 1 / xi) * log1p(xi * y / sigma)
}
y <- c(0.1, 0.3, 1, 2)
for (xi in c(-.4, 0, .2, 1.5))
  near(D(y, rep(.3, 4), rep(xi, 4), 0)$l0, logdensity(y, exp(.3), xi))

# Link transformation and coefficient assembly, including weighted observations
# and offsets. Check transformed third/fourth derivatives separately as well.
X <- cbind(1, c(-.5, -.2, .1, .6), 1, c(.2, -.1, .3, -.3))
attr(X, "lpi") <- list(1:2, 3:4)
wt <- c(0, .5, 2, 1)
off <- list(rep(.15, 4), rep(-.1, 4))
for (link in c("identity", "logshift")) {
  fam <- gpd(c("log", link))
  b <- c(.3, .1, fam$linfo[[2]]$linkfun(.1), .1)
  ll <- function(b, deriv = 0) fam$ll(y, X, b, wt, fam, off, deriv = deriv)
  d <- ll(b, 1)
  for (j in 1:4) {
    near(fd(function(b) ll(b)$l, b, j), d$lb[j], 1e-6)
    near(fd(function(b) ll(b, 1)$lb, b, j), d$lbb[, j], 1e-6)
  }
  transform <- function(p) {
    xi <- fam$linfo[[2]]$linkinv(p[2]); sh <- fam$linfo[[2]]
    ds <- D(1.1, p[1], xi, 4)
    mgcv::gamlss.etamu(ds$l1, ds$l2, ds$l3, ds$l4,
      cbind(1, sh$mu.eta(p[2])), cbind(0, sh$d2link(xi)),
      cbind(0, sh$d3link(xi)), cbind(0, sh$d4link(xi)),
      fam$tri$i2, fam$tri$i3, fam$tri$i4, 3)
  }
  for (shape in c(-.4, 0, .3)) {
    p <- c(.2, fam$linfo[[2]]$linkfun(shape))
    ds <- transform(p)
    for (ord in 2:4) {
      near(fd(function(p) as.numeric(transform(p)[[ord - 1]]), p, 1),
           ds[[ord]][, 1:ord], 1e-6)
      near(fd(function(p) as.numeric(transform(p)[[ord - 1]]), p, 2),
           ds[[ord]][, 2:(ord + 1L)], 1e-6)
    }
  }
}

fam <- gpd()
stopifnot(inherits(fam, "general.family"), fam$nlp == 2L,
          is.function(fam$cdf), is.function(fam$qf), is.function(fam$rd))
for (bad in list("log", c("identity", "identity"), c("log", "log"),
                 list("log", NA_character_), list(c("log", "log"), "identity")))
  expect_error(gpd(bad), "gpd requires links")
mu <- cbind(c(.5, 1, 2, 5, 1, 1), c(-.4, 0, 1e-12, .2, 1, 2))
for (p in c(0, 1e-12, .01, .5, .99, 1 - 1e-12, 1)) {
  q <- fam$qf(p, mu)
  near(fam$cdf(q, mu), rep(p, nrow(mu)), 1e-9)
  near(exp(fam$cdf(q, mu, logp = TRUE)), rep(p, nrow(mu)), 1e-9)
}
near(fam$qf(c(.2, .7), cbind(c(2, 3), 0)), qexp(c(.2, .7), rate = 1 / c(2, 3)))
near(fam$cdf(c(.2, .7), cbind(c(2, 3), 0)), pexp(c(.2, .7), rate = 1 / c(2, 3)))
near(fam$qf(1, matrix(c(2, -.4), 1)), 5)
near(fam$cdf(c(-Inf, -1, 0, 5, 6, Inf), cbind(rep(2, 6), -.4)), c(0, 0, 0, 1, 1, 1))
near(fam$cdf(1e308, cbind(c(1e-200, 1e-200), c(0, .5))), c(1, 1))
stopifnot(is.na(fam$qf(NA_real_, matrix(c(1, 0), 1))),
          is.na(fam$cdf(1, matrix(c(NA_real_, 0), 1))))
stopifnot(is.nan(expect_warning(fam$qf(2, matrix(c(1, 0), 1)), "NaNs")),
          is.nan(expect_warning(fam$cdf(1, matrix(c(-1, 0), 1)), "NaNs")))
expect_error(fam$qf(c(.2, .3, .4), mu), "one value per row")
expect_error(fam$cdf(1, c(1, 0)), "matrix")
set.seed(321); u <- runif(nrow(mu))
set.seed(321); near(fam$rd(mu), fam$qf(u, mu))
near(fam$qf(.5, mu, wt = 20, scale = 10), fam$qf(.5, mu))

# Ordinary unpenalized MLE fits compared with independently optimized likelihood.
set.seed(534)
for (shape in c(-.3, 0, .3)) {
  d <- data.frame(y = fam$rd(cbind(rep(2, 600), shape)))
  opt <- optim(c(log(2), .05), function(p) {
    ll <- logdensity(d$y, exp(p[1]), p[2])
    if (any(!is.finite(ll))) return(1e100)
    -sum(ll)
  }, method = "Nelder-Mead", control = list(reltol = 1e-12, maxit = 3000))
  for (link in c("identity", "logshift")) {
    b <- gam(list(y ~ 1, ~ 1), data = d, family = gpd(c("log", link)), method = "REML")
    near(b$fitted.values[1, ], c(exp(opt$par[1]), opt$par[2]), 2e-4)
    near(as.numeric(logLik(b)), -opt$value, 1e-7)
    stopifnot(all(is.finite(residuals(b))), is.na(b$null.deviance))
  }
}

# Smooth scale AND shape, including negative, near-zero and positive shapes.
set.seed(142)
d <- data.frame(x = runif(700))
truth <- cbind(exp(.3 + .5 * sin(6 * d$x)), -.2 + .4 * d$x)
d$y <- fam$rd(truth)
for (link in c("identity", "logshift")) {
  b <- gam(list(y ~ s(x, k = 5), ~ s(x, k = 4)), data = d,
           family = gpd(c("log", link)), method = "REML")
  stopifnot(all(is.finite(coef(b))), all(1 + b$fitted.values[, 2] * d$y / b$fitted.values[, 1] > 0),
            mean(abs(log(b$fitted.values[, 1] / truth[, 1]))) < .3)
  pr <- predict(b, type = "response", se.fit = TRUE)
  lp <- predict(b, type = "link", se.fit = TRUE)
  near(pr$fit, b$fitted.values)
  near(pr$fit[, 1], exp(lp$fit[, 1]))
  near(pr$fit[, 2], b$family$linfo[[2]]$linkinv(lp$fit[, 2]))
  near(pr$se.fit[, 1], exp(lp$fit[, 1]) * lp$se.fit[, 1])
  near(pr$se.fit[, 2], b$family$linfo[[2]]$mu.eta(lp$fit[, 2]) * lp$se.fit[, 2])
  nd <- data.frame(x = c(.1, .4, .8))
  pm <- predict(b, nd, type = "response")
  for (p in c(.1, .9)) near(fam$cdf(fam$qf(p, pm), pm), rep(p, nrow(pm)))
  set.seed(4)
  sim <- fam$rd(pm)
  set.seed(4)
  near(sim, fam$qf(runif(nrow(pm)), pm))
  stopifnot(length(sim) == 3L, all(sim >= 0))
  pit <- fam$cdf(d$y, b$fitted.values)
  stopifnot(all(is.finite(qnorm(pit))), abs(mean(pit) - .5) < .05)
}

# Integer likelihood weights agree with replicated observations. A zero-weight
# observation outside the fitted support cannot affect the result.
set.seed(832)
d <- data.frame(y = fam$rd(cbind(rep(1, 150), -.25)), w = rep(1:3, 50))
b <- gam(list(y ~ 1, ~ 1), data = d, weights = w, family = fam, method = "REML")
br <- gam(list(y ~ 1, ~ 1), data = d[rep(seq_len(nrow(d)), d$w), ], family = fam, method = "REML")
near(coef(b), coef(br), 1e-6)
dz <- rbind(d, data.frame(y = 1e10, w = 0))
bz <- gam(list(y ~ 1, ~ 1), data = dz, weights = w, family = fam, method = "REML")
near(coef(b), coef(bz), 1e-6)
stopifnot(tail(residuals(bz), 1) == 0)
bs <- gam(list(y ~ 1, ~ 1), data = d, weights = w, family = fam,
          start = coef(b), method = "REML")
near(coef(b), coef(bs), 1e-6)
d$o1 <- .2; d$o2 <- -.1
bo <- gam(list(y ~ 1 + offset(o1), ~ 1 + offset(o2)), data = d,
          weights = w, family = fam, method = "REML")
near(coef(bo) + c(.2, -.1), coef(b), 1e-6)
near(predict(bo, type = "response"), b$fitted.values, 1e-6)
pm <- predict(bo, d[1:4, ], type = "response")
near(fam$cdf(fam$qf(.5, pm), pm), rep(.5, 4))

for (yy in list(c(-1, 1, 2), c(Inf, 1, 2), c(0, 0, 0)))
  expect_error(gam(list(y ~ 1, ~ 1), data = data.frame(y = yy), family = fam), "gpd requires")
expect_error(gam(list(y ~ 1, ~ 1), data = d, family = fam, weights = rep(0, nrow(d))), "weights")
expect_error(gam(list(y ~ 1, ~ 1), data = d, family = fam, weights = rep(-1, nrow(d))), "negative")
expect_error(gam(list(y ~ 1, ~ 1), data = d, family = fam, start = c(-10, -10)), "starting coefficients")
dna <- rbind(d, d[1, ]); dna$y[nrow(dna)] <- NA
bn <- gam(list(y ~ 1, ~ 1), data = dna, family = fam, weights = w, method = "REML", na.action = na.omit)
near(coef(bn), coef(b), 1e-6)

# Residual conventions and warnings for nonexistent moments.
o <- list(y = c(0, 1, 2, 3), fitted.values = cbind(rep(2, 4), 0), prior.weights = c(1, 1, 2, 0))
r <- fam$residuals(o)
stopifnot(r[1] == -Inf, r[3] == 0, r[4] == 0)
near(r[2], -sqrt(2 * (log(2) + .5 - 1)))
near(fam$residuals(o, "response"), o$y - 2)
near(fam$residuals(o, "pearson"), (o$y - 2) * sqrt(o$prior.weights) / 2)
o$fitted.values[, 2] <- c(.6, 1, -1.1, 2)
r <- expect_warning(fam$residuals(o, "pearson"), "moments")
stopifnot(is.na(r[1]), is.na(r[2]), r[4] == 0)
r <- expect_warning(fam$residuals(o, "response"), "moments")
stopifnot(is.na(r[2]), is.na(r[4]))
r <- expect_warning(fam$residuals(o), "shape <= -1")
stopifnot(is.na(r[3]))

# Extreme but representable quantiles/CDFs must not overflow intermediate
# products. These references evaluate the defining equations on the log scale.
q <- fam$qf(1 - exp(-20), matrix(c(1e-200, 40), 1))
stopifnot(is.finite(q), q > 0)
near(fam$cdf(q, matrix(c(1e-200, 40), 1)), 1 - exp(-20))
near(fam$cdf(1e308, matrix(c(1, 10000), 1)),
     -expm1(-(log(10000) + log(1e308)) / 10000))
near(log(fam$qf(.9, matrix(c(1, -1e308), 1))), -log(1e308))
near(log(fam$cdf(1e-200, matrix(c(1, 1e-308), 1))), -200 * log(10))
cat("gpd derivatives, fits, distribution helpers and diagnostics: passed\n")
