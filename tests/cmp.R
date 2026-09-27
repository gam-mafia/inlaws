library(inlaws)
near_cmp <- function(x, y, tol = 1e-7) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
            all(abs(x - y) <= tol * (1 + abs(y))))
}
error_cmp <- function(expr, pattern) {
  err <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err)))
}
moments <- getFromNamespace(".cmp_moments", "inlaws")
control <- getFromNamespace(".cmp_control", "inlaws")()

# Independent sums, with support extended past the adaptive algorithm's support.
grid <- expand.grid(mu = c(1e-6, .01, .2, 1, 2.5, 10, 100),
                    nu = c(.1, .3, 1, 3, 10, 30))
for (i in seq_len(nrow(grid))) {
  m <- grid$mu[i]; nu <- grid$nu[i]
  z <- moments(m, nu)[[1]]
  j <- 0:(max(z$y) + 100)
  lp <- j * z$a - nu * lgamma(j + 1)
  w <- exp(lp - max(lp)); p <- w / sum(w)
  stopifnot(abs(z$mean - m) <= 1e-9 * m)
  near_cmp(z$logZ, max(lp) + log(sum(w)))
  near_cmp(z$variance, sum(p * (j - m)^2))
  near_cmp(z$p, p[match(z$y, j)])
}
# Helpers accept observation-specific dispersion for the future location-scale family.
near_cmp(vapply(moments(c(.1, 2, 10), c(.2, 1, 5)), `[[`, numeric(1), "mean"), c(.1, 2, 10))

# Analytic derivatives, checked by independent differencing of neighbouring
# likelihood values and lower-order analytic derivatives (not fourth differences).
f <- cmp(control = list(sum_tol = 1e-12, mean_tol = 1e-12))
y <- c(0, 1, 3, 8); mu <- c(.2, .8, 2.5, 10); w <- c(.5, 1, 2, 1)
for (nu in c(.3, 1, 3)) {
  th <- log(nu); h <- 1e-5
  d <- f$Dd(y, mu, th, w, 2)
  plus <- f$Dd(y, mu + h, th, w, 2); minus <- f$Dd(y, mu - h, th, w, 2)
  tp <- f$Dd(y, mu, th + h, w, 2); tm <- f$Dd(y, mu, th - h, w, 2)
  checks <- list(
    Dmu = (f$dev.resids(y, mu+h, w, th)-f$dev.resids(y, mu-h, w, th))/(2*h),
    Dmu2 = (plus$Dmu-minus$Dmu)/(2*h),
    Dmu3 = (plus$Dmu2-minus$Dmu2)/(2*h),
    Dmu4 = (plus$Dmu3-minus$Dmu3)/(2*h),
    Dth = (f$dev.resids(y, mu, w, th+h)-f$dev.resids(y, mu, w, th-h))/(2*h),
    Dth2 = (tp$Dth-tm$Dth)/(2*h), Dmuth = (tp$Dmu-tm$Dmu)/(2*h),
    Dmu2th = (tp$Dmu2-tm$Dmu2)/(2*h), Dmuth2 = (tp$Dmuth-tm$Dmuth)/(2*h),
    Dmu2th2 = (tp$Dmu2th-tm$Dmu2th)/(2*h), Dmu3th = (tp$Dmu3-tm$Dmu3)/(2*h))
  for (nm in names(checks)) near_cmp(d[[nm]], checks[[nm]], 2e-6)
  s <- f$ls(y, w, th, 1); sp <- f$ls(y, w, th+h, 1); sm <- f$ls(y, w, th-h, 1)
  near_cmp(s$lsth1, (sp$ls-sm$ls)/(2*h), 2e-6)
  near_cmp(s$lsth2, (sp$lsth1-sm$lsth1)/(2*h), 2e-6)
  near_cmp(sum(s$LSTH1), s$lsth1)
  near_cmp(sum(f$dev.resids(y, mu, w, th))/2 - s$ls,
           f$aic(y, mu, theta = th, wt = w, dev = NULL)/2)
  # Expected curvature and its mixed derivative, averaging over an independent PMF.
  z <- moments(2, nu)[[1]]; j <- z$y
  dj <- f$Dd(j, rep(2, length(j)), th, 1, 2)
  near_cmp(sum(z$p*dj$Dmu), 0)
  near_cmp(sum(z$p*dj$Dmu2), dj$EDmu2[1])
  near_cmp(sum(z$p*dj$Dmu3), dj$EDmu3[1])
  near_cmp(sum(z$p*dj$Dmu2th), dj$EDmu2th[1])
}

# Poisson identity including offset, non-unit weights, smoothing and diagnostics.
set.seed(240925)
dat <- data.frame(x = runif(90), exposure = runif(90, .5, 2))
dat$mu <- exp(.3 + sin(6 * dat$x)) * dat$exposure
dat$y <- rpois(90, dat$mu)
dat$w <- rep(c(0, .5, 1, 2, 1), length.out = 90)
form <- y ~ s(x, k = 6) + offset(log(exposure))
bp <- gam(form, data = dat, weights = w, family = poisson(), method = "REML")
bc <- gam(form, data = dat, weights = w, family = cmp(1), method = "REML")
near_cmp(coef(bc), coef(bp), 2e-5)
near_cmp(bc$sp, bp$sp, 2e-5)
near_cmp(fitted(bc), fitted(bp), 2e-5)
near_cmp(as.numeric(logLik(bc)), as.numeric(logLik(bp)), 2e-6)
near_cmp(bc$deviance, bp$deviance, 2e-6)
near_cmp(bc$null.deviance, bp$null.deviance, 2e-6)
near_cmp(AIC(bc), AIC(bp), 2e-6)
near_cmp(residuals(bc, type = "deviance"), residuals(bp, type = "deviance"), 2e-5)
near_cmp(predict(bc, newdata = dat[1:8, ], type = "response"), fitted(bc)[1:8])
stopifnot(all(is.finite(predict(bc, se.fit = TRUE)$se.fit)))
pp <- cmp(1)
near_cmp(pp$variance(c(.01, 1, 10)), c(.01, 1, 10))
stopifnot(length(pp$variance(numeric())) == 0L,
          length(pp$qf(numeric(), 1)) == 0L,
          pp$aic(2, 1, wt = 0, dev = 0) == 0)
near_cmp(pp$cdf(c(-1, 0, 3, Inf), c(.1, .2, 2, 10)), ppois(c(-1, 0, 3, Inf), c(.1, .2, 2, 10)))
near_cmp(pp$cdf(3, 2, logp = TRUE), ppois(3, 2, log.p = TRUE))
near_cmp(pp$qf(c(.01, .5, .9), c(.1, 2, 10)), qpois(c(.01, .5, .9), c(.1, 2, 10)))
stopifnot(pp$qf(0, 2) == 0, is.infinite(pp$qf(1, 2)), pp$qf(1, 0) == 0,
          pp$cdf(0, 0) == 1, pp$cdf(-1, 0) == 0, is.na(pp$cdf(NA, 1)))
near_cmp(pp$dev.resids(c(0, 1000), c(.2, Inf), c(1, 0)), c(.4, 0))
stopifnot(all(vapply(pp$Dd(1000, Inf, 0, 0, 2), function(x) x == 0, logical(1))))
set.seed(2)
ss <- replicate(2, bc$family$rd(fitted(bc)))
stopifnot(identical(dim(ss), c(90L, 2L)), all(ss >= 0), all(ss == floor(ss)))

# Estimated dispersion: independent mean inversion + direct-sum likelihood.
# Intercept-only mean MLE equals the sample mean; profile nu independently.
set.seed(320)
for (nu in c(.5, 3)) {
  yy <- cmp(nu)$rd(rep(3, 160))
  fit <- gam(yy ~ 1, family = cmp(), method = "REML")
  reference <- function(t, restricted = TRUE) {
    jj <- 0:500; nu0 <- exp(t)
    eval <- function(a) {
      lp <- jj*a - nu0*lgamma(jj+1); ww <- exp(lp-max(lp))
      list(mean = sum(jj*ww)/sum(ww), logZ = max(lp)+log(sum(ww)),
           variance = sum(ww*(jj-mean(yy))^2)/sum(ww))
    }
    a <- uniroot(function(a) eval(a)$mean - mean(yy), c(-20, 20), tol=1e-11)$root
    z <- eval(a)
    -sum(yy*a - nu0*lgamma(yy+1) - z$logZ) +
      if (restricted) .5*log(length(yy)*mean(yy)^2/z$variance) else 0
  }
  ref <- optimize(reference, c(log(.1), log(10)), tol=1e-8)
  near_cmp(fit$family$getTheta(), ref$minimum, 2e-5)
  near_cmp(as.numeric(logLik(fit)), -reference(ref$minimum, FALSE), 2e-6)
  near_cmp(fitted(fit), rep(mean(yy), length(yy)))
}

# Smooth-model recovery and independent family state when reusing an object.
set.seed(81)
dat$y <- cmp(2)$rd(dat$mu)
fam <- cmp()
b1 <- gam(form, data = dat, family = fam, method = "REML")
old <- b1$family$getTheta(TRUE)
stopifnot(b1$converged, old > .8, old < 5,
          sqrt(mean((log(fitted(b1)) - log(dat$mu))^2)) < .6)
dat$y <- cmp(.5)$rd(dat$mu)
b2 <- gam(form, data = dat, family = fam, method = "REML")
stopifnot(b2$converged, b2$family$getTheta(TRUE) < 1,
          identical(old, b1$family$getTheta(TRUE)), fam$getTheta(TRUE) == 1)
stopifnot(cmp(-2)$getTheta(TRUE) == 2, cmp(-2)$n.theta == 1L,
          cmp(2)$n.theta == 0L, cmp(0)$n.theta == 1L)

error_cmp(cmp(link = "identity"), "log link")
error_cmp(cmp(theta = NA_real_), "finite")
error_cmp(cmp(control = list(typo = 1)), "control")
error_cmp(cmp(control = list(max_terms = 1.5)), "control")
error_cmp(cmp(control = list(sum_tol = 0)), "control")
error_cmp(cmp(control = list(max_terms = 5))$variance(20), "max_terms")
error_cmp(cmp(2, control = list(max_iter = 1))$variance(2), "max_iter")
error_cmp(cmp(1)$qf(1.2, 2), "probabilities")
error_cmp(gam(c(0, 1, 2, 3, 4) ~ 1, family = cmp(1), method = "REML", scale = 2), "scale is fixed")
error_cmp(gam(c(0, .5, 1) ~ 1, family = cmp(1), method = "REML"), "integers")
error_cmp(gam(c(-1, 1, 2) ~ 1, family = cmp(1), method = "REML"), "integers")
error_cmp(gam(c(0, 1, 2) ~ 1, weights = c(1, -1, 1), family = cmp(1), method = "REML"), "negative|nonnegative")
error_cmp(gam(rep(0, 5) ~ 1, family = cmp(1), method = "REML"), "positive count")
