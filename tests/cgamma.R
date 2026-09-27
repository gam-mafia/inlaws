library(inlaws)

near <- function(x, y, tol = 2e-6) {
  err <- max(abs(x - y) / (1 + abs(y)))
  if (!is.finite(err) || err > tol) stop("relative derivative/error check: ", err)
}
expect_error <- function(expr, text) {
  err <- tryCatch(expr, error = identity)
  stopifnot(inherits(err, "error"), grepl(text, conditionMessage(err), fixed = TRUE))
}
# Richardson differentiation avoids a dependency on a differentiation package.
first <- function(f, x, h = 1e-4) {
  d1 <- (f(x + h) - f(x - h)) / (2 * h)
  d2 <- (f(x + h / 2) - f(x - h / 2)) / h
  (4 * d2 - d1) / 3
}
second <- function(f, x, h = 0.002) {
  d1 <- (f(x + h) - 2 * f(x) + f(x - h)) / h^2
  d2 <- (f(x + h / 2) - 2 * f(x) + f(x - h / 2)) / (h / 2)^2
  (4 * d2 - d1) / 3
}

# All censoring cases, including support endpoints, uninformative observations,
# a narrow interval, and a zero-weight row.
y <- cbind(c(.8, .5, 3, .7, 0, 0, 1, 2),
           c(.8, -Inf, Inf, 1.4, .2, Inf, 1 + 1e-9, -Inf))
mu <- c(.9, 1.2, 2, 1, 1, 1, .8, 2)
w <- c(1, .3, 2, 1, 1, 1, 1, 0)
fam <- cgamma()
for (phi in c(.1, .7, 3)) {
  th <- log(phi)
  d <- fam$Dd(y, mu, th, w, 2)
  near(d$Dmu, first(function(e) fam$dev.resids(y, mu + e, w, th), 0))
  near(d$Dth, first(function(t) fam$dev.resids(y, mu, w, t), th))
  for (pair in list(c("Dmu2", "Dmu"), c("Dmu3", "Dmu2"), c("Dmu4", "Dmu3"))) {
    near(d[[pair[1]]], first(function(e) fam$Dd(y, mu + e, th, w, 2)[[pair[2]]], 0))
  }
  for (pair in list(c("Dth2", "Dth"), c("Dmuth", "Dmu"), c("Dmuth2", "Dmuth"),
                   c("Dmu2th", "Dmu2"), c("Dmu2th2", "Dmu2th"), c("Dmu3th", "Dmu3"))) {
    near(d[[pair[1]]], first(function(t) fam$Dd(y, mu, t, w, 2)[[pair[2]]], th))
  }
  sat <- fam$ls(y, w, th, 1)
  near(sat$lsth1, first(function(t) fam$ls(y, w, t, 1)$ls, th))
  near(sat$lsth2, first(function(t) fam$ls(y, w, t, 1)$lsth1, th))
  near(sum(fam$dev.resids(y, mu, w, th)) / 2 - sat$ls,
       fam$aic(y, mu, th, w, 0) / 2)
  stopifnot(all(vapply(d, function(z) z[8] == 0, logical(1))))
}

# Independent derivatives of base-R gamma tail probabilities, including tails
# whose ordinary (non-log) probabilities underflow to zero.
tail_derivs <- getFromNamespace(".cg_tail", "inlaws")
for (k in c(.03, .2, 1, 2, 10, 100, 1000)) {
  for (ratio in c(.001, .1, .9, 1, 1.1, 2, 10, 100)) {
    for (lower in c(TRUE, FALSE)) {
      theta <- -log(k)
      got <- tail_derivs(k * ratio, k, lower)
      f <- function(t) pgamma(exp(-t) * ratio, exp(-t), lower.tail = lower, log.p = TRUE)
      near(got[1, 1], f(theta), 1e-12)
      near(got[1, 2], first(f, theta), 1e-6)
      near(got[1, 3], second(f, theta), 3e-5)
    }
  }
}

# An independent integral check of shape derivatives: differentiate the density
# under the integral, instead of differentiating incomplete-gamma recurrences.
loglik <- getFromNamespace(".cg_loglik", "inlaws")
for (k in c(.7, 2, 20)) for (bounds in list(c(0, .2), c(.3, 2), c(4, Inf))) {
  m <- 1.3
  lo <- bounds[1]; hi <- bounds[2]
  response <- if (is.infinite(hi)) cbind(lo, Inf) else cbind(lo, hi)
  z <- loglik(response, m, -log(k))
  score <- function(y) -k * (log(k) + 1 - digamma(k) + log(y / m) - y / m)
  density <- function(y) exp(dgamma(y, k, scale = m / k, log = TRUE) - z[1, 1])
  s1 <- integrate(function(y) density(y) * score(y), lo, hi,
                  rel.tol = 1e-8, subdivisions = 1000)$value
  s2 <- integrate(function(y) density(y) *
                    (score(y)^2 - score(y) + k - k^2 * trigamma(k)), lo, hi,
                  rel.tol = 1e-8, subdivisions = 1000)$value - s1^2
  near(z[1, 2:3], c(s1, s2), 2e-6)
}

# Far upper tails need a different mean-derivative evaluation: cancellations
# between probability ratios otherwise corrupt third and fourth derivatives.
for (k in c(.03, 1, 2, 100, 1000)) {
  x <- max(101, 10 * k)
  for (upper in c(Inf, 1.001 * x / k, 2 * x / k)) {
    yy <- cbind(x / k, upper)
    th <- -log(k)
    d <- fam$Dd(yy, 1, th, 1, 2)
    for (pair in list(c("Dmu2", "Dmu"), c("Dmu3", "Dmu2"), c("Dmu4", "Dmu3"))) {
      near(d[[pair[1]]], first(function(m) fam$Dd(yy, m, th, 1, 2)[[pair[2]]], 1))
    }
    for (pair in list(c("Dmuth", "Dmu"), c("Dmuth2", "Dmuth"),
                     c("Dmu2th", "Dmu2"), c("Dmu2th2", "Dmu2th"), c("Dmu3th", "Dmu3"))) {
      near(d[[pair[1]]], first(function(t) fam$Dd(yy, 1, t, 1, 2)[[pair[2]]], th))
    }
  }
}
# At shape two, Q(2,x) = exp(-x) * (1+x), giving a fully independent
# closed-form derivative oracle even when ordinary probabilities underflow.
d <- fam$Dd(cbind(50000, Inf), 1, -log(2), 1, 2)
for (r in 1:4) {
  reference <- -2 * ((-1)^(r + 1) * factorial(r) * 1e5 +
                     (-1)^(r - 1) * factorial(r - 1) * ((1 + 1e5)^(-r) - 1))
  near(d[[c("Dmu", "Dmu2", "Dmu3", "Dmu4")[r]]], reference, 1e-10)
}

# A narrow interval approaches density times width; its derivatives approach
# those of an exact observation. Check all mean and mixed derivatives.
theta <- log(.4)
yn <- cbind(1, 1 + 1e-10)
ye <- mean(yn)
dn <- fam$Dd(yn, .8, theta, 1, 2)
de <- fam$Dd(ye, .8, theta, 1, 2)
for (nm in setdiff(names(de), c("EDmu2", "EDmu2th"))) near(dn[[nm]], de[[nm]], 1e-6)
near(fam$aic(yn, .8, theta, 1, 0),
     -2 * (dgamma(ye, 2.5, scale = .8 / 2.5, log = TRUE) + log(diff(as.numeric(yn)))), 1e-10)

# Exact likelihood and deviance agree with the standard gamma formulas.
phi <- .4
ye <- c(.2, .8, 2, 5)
me <- c(.5, .9, 1, 3)
we <- c(.5, 1, 2, 0)
near(fam$aic(ye, me, log(phi), we, 0), -2 * sum(we * dgamma(ye, 1 / phi, scale = me * phi, log = TRUE)))
near(fam$dev.resids(ye, me, we, log(phi)), Gamma()$dev.resids(ye, me, we) / phi)

# Fixed-dispersion smooth fits must reproduce ordinary gamma fits.
set.seed(402)
dat <- data.frame(x = runif(250))
dat$y <- rgamma(250, shape = 3, scale = exp(.5 + sin(2 * pi * dat$x)) / 3)
standard <- mgcv::gam(y ~ s(x, k = 8), data = dat, family = Gamma(link = "log"),
                      scale = 1 / 3, method = "REML")
fixed <- mgcv::gam(y ~ s(x, k = 8), data = dat, family = cgamma(theta = 1 / 3), method = "REML")
near(coef(fixed), coef(standard), 2e-5)
# Extended-family deviance includes inverse dispersion, so the smoothing
# parameter is expressed in different units.
near(fixed$sp / 3, standard$sp, 2e-5)
near(predict(fixed, type = "response"), predict(standard, type = "response"), 2e-5)
near(fixed$edf, standard$edf, 2e-5)
near(vcov(fixed), vcov(standard), 2e-5)
near(fixed$family$getTheta(TRUE), 1 / 3)
stopifnot(fixed$scale == 1)

# With no censoring and estimated dispersion, recover the ordinary Gamma
# model under both REML and ML. Check both accepted response encodings.
# Gamma stores its likelihood dispersion in reml.scale; sig2/scale instead
# uses the separate Fletcher estimator by default. Compare like with like.
dat$exact <- cbind(dat$y, dat$y)
prediction_data <- data.frame(x = seq(0, 1, length.out = 101))
for (method in c("REML", "ML")) {
  ordinary <- mgcv::gam(y ~ s(x, k = 8), data = dat, family = Gamma(link = "log"), method = method)
  for (response in c("y", "exact")) {
    candidate <- mgcv::gam(stats::reformulate("s(x, k = 8)", response), data = dat,
                           family = cgamma(), method = method)
    phi <- candidate$family$getTheta(TRUE)
    near(phi, ordinary$reml.scale, 2e-5)
    near(coef(candidate), coef(ordinary), 2e-5)
    near(fitted(candidate), fitted(ordinary), 2e-5)
    near(candidate$edf, ordinary$edf, 2e-5)
    near(candidate$sp * phi, ordinary$sp, 2e-5)
    near(candidate$gcv.ubre, ordinary$gcv.ubre, 2e-7)
    near(candidate$deviance * phi, ordinary$deviance, 2e-5)
    near(vcov(candidate), vcov(ordinary) * ordinary$reml.scale / ordinary$sig2, 2e-5)
    p0 <- predict(ordinary, prediction_data, type = "response", se.fit = TRUE)
    p1 <- predict(candidate, prediction_data, type = "response", se.fit = TRUE)
    near(p1$fit, p0$fit, 2e-5)
    near(p1$se.fit, p0$se.fit * sqrt(ordinary$reml.scale / ordinary$sig2), 2e-5)
  }
}

# Independent MLE for a regression with mixed censoring, fitted using a plain
# base-R likelihood with density/CDF/survival terms, rather than family helpers.
dat$yc <- cbind(dat$y, dat$y)
il <- dat$y < .6
ir <- dat$y > 4
ii <- seq_len(nrow(dat)) %% 4 == 0 & !il & !ir
dat$yc[il, 1] <- .6; dat$yc[il, 2] <- -Inf
dat$yc[ir, 1] <- 4; dat$yc[ir, 2] <- Inf
dat$yc[ii, 1] <- floor(dat$y[ii] / .2) * .2
dat$yc[ii, 2] <- dat$yc[ii, 1] + .2
exact <- !(il | ir | ii)
nll <- function(par) {
  m <- exp(par[1] + par[2] * dat$x)
  k <- exp(-par[3])
  ll <- numeric(nrow(dat))
  ll[exact] <- dgamma(dat$y[exact], k, scale = m[exact] / k, log = TRUE)
  ll[il] <- pgamma(.6, k, scale = m[il] / k, log.p = TRUE)
  ll[ir] <- pgamma(4, k, scale = m[ir] / k, lower.tail = FALSE, log.p = TRUE)
  ll[ii] <- log(pgamma(dat$yc[ii, 2], k, scale = m[ii] / k) -
                  pgamma(dat$yc[ii, 1], k, scale = m[ii] / k))
  -sum(ll)
}
fit <- mgcv::gam(yc ~ x, data = dat, family = cgamma(), method = "ML")
opt <- optim(c(0, 0, 0), nll, method = "BFGS", control = list(reltol = 1e-12))
near(c(coef(fit), fit$family$getTheta()), opt$par, 1e-5)
near(fit$family$aic(dat$yc, fitted(fit), wt = rep(1, nrow(dat)), dev = 0) / 2, opt$value)
# Distinct family instances must not share dispersion state.
f1 <- cgamma(); f2 <- cgamma(); f1$putTheta(log(.2))
near(f2$getTheta(TRUE), 1)

# Smooth fits, offsets, predictions and matrix/attribute subsetting.
fit <- mgcv::gam(yc ~ s(x, k = 8) + offset(rep(.2, nrow(dat))), data = dat,
                 family = cgamma(), method = "REML")
stopifnot(fit$converged, fit$family$getTheta(TRUE) > 0)
pr <- predict(fit, type = "response", se.fit = TRUE)
near(pr$fit, exp(predict(fit, type = "link")))
stopifnot(all(is.finite(pr$se.fit)), all(pr$se.fit >= 0))
# The null deviance must optimise the censored likelihood, not use the mean
# of the supplied censoring limits as an intercept estimate.
null <- optimize(function(eta) sum(fit$family$dev.resids(dat$yc,
                  rep(exp(eta + .2), nrow(dat)), rep(1, nrow(dat)))), c(-5, 5))
near(fit$null.deviance, null$objective)
# bam's two fitting paths use the working-curvature fields differently from
# gam. Check both estimated-dispersion paths on the mixed-censoring response.
for (discrete in c(FALSE, TRUE)) {
  b <- mgcv::bam(yc ~ s(x, k = 8), data = dat, family = cgamma(),
                 method = "fREML", discrete = discrete)
  stopifnot(all(is.finite(coef(b))), all(is.finite(fitted(b))),
            is.finite(b$family$getTheta(TRUE)), b$family$getTheta(TRUE) > 0)
}
fitted_y <- fam$subsety(dat$yc, 1)
stopifnot(identical(dim(fitted_y), c(1L, 2L)))
ya <- dat$yc[, 1]; attr(ya, "censor") <- dat$yc[, 2]
stopifnot(identical(attr(fam$subsety(ya, c(2, 4)), "censor"), dat$yc[c(2, 4), 2]))

expect_error(cgamma(link = "identity"), "only the log link")
expect_error(cgamma(theta = NA_real_), "finite numeric dispersion")
expect_error(cgamma(theta = c(1, 2)), "finite numeric dispersion")
expect_error(fam$response(c(0, 1)), "positive exact")
expect_error(fam$response(cbind(1, .5)), "ordered interval")
expect_error(fam$response(cbind(0, -Inf)), "positive exact")
expect_error(fam$response(matrix(1, 2, 3)), "two-column")
expect_error(mgcv::gam(y ~ x, data = dat, weights = rep(0, nrow(dat)), family = cgamma()), "at least one positive")
