library(inlaws)

near <- function(x, y, tol = 1e-7) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol,
    scale = max(1, mean(abs(y[is.finite(y)]))))))
}
error <- function(expr, pattern) {
  e <- tryCatch(expr, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e), fixed = TRUE))
}
ll <- getFromNamespace(".gp_logpmf", "inlaws")
saturated <- getFromNamespace(".gp_saturated", "inlaws")

# Constructor and invalid inputs.
stopifnot(inherits(gp(), "extended.family"), gp()$n.theta == 1L,
          gp(2)$n.theta == 0L, length(gp(1)$getTheta()) == 0L)
near(gp()$getTheta(TRUE), 2)
near(gp(0)$getTheta(TRUE), 2)
near(gp(-3)$getTheta(TRUE), 3)
for (bad in list(NA_real_, Inf, -Inf, .5, -.5, -1, c(1, 2), "2", numeric()))
  error(gp(bad), "theta must")
error(gp(link = "identity"), "log link")

# Independent direct PMF, normalization, and exact GP-1 moments.
for (phi in c(1, 1 + 1e-7, 2, 10, 50)) {
  f <- gp(phi)
  for (m in c(.05, 1, 20, 500)) {
    y <- 0:100000
    p <- exp(ll(y, rep(m, length(y)), log(phi - 1)))
    near(sum(p), 1, 1e-10)
    near(sum(y * p), m, 1e-9)
    near(sum((y - m)^2 * p), m * phi, 1e-8)
    j <- 0:25
    a <- m / sqrt(phi); b <- 1 - 1 / sqrt(phi)
    direct <- a * (a + b*j)^(j - 1) * exp(-a - b*j) / factorial(j)
    near(exp(ll(j, rep(m, length(j)), log(phi - 1))), direct, 1e-10)
    near(f$cdf(j, m), cumsum(p)[j + 1], 1e-10)
    probs <- c(.001, .1, .5, .9, .999)
    qs <- f$qf(probs, m)
    stopifnot(all(f$cdf(qs, m) >= probs - 1e-13),
              all(f$cdf(qs - 1, m) < probs + 1e-13))
    near(f$cdf(j, m, logp = TRUE), log(f$cdf(j, m)), 1e-10)
  }
}
f <- gp(2)
near(f$cdf(c(-Inf, -1, 0, 1, Inf), 0), c(0, 0, 1, 1, 1))
near(f$qf(c(0, .1, 1), 0), c(0, 0, 0))
near(f$qf(c(0, 1), 2), c(0, Inf))
stopifnot(is.na(f$cdf(NA_real_, 1)), is.na(f$qf(NA_real_, 1)))
stopifnot(is.nan(suppressWarnings(f$qf(-.1, 2))))
near(f$cdf(c(.2, 1.9, 2.1), 2), f$cdf(0:2, 2))
error(f$cdf(1, -1), "nonnegative")
error(f$qf(1:3 / 4, c(1, 2)), "equal lengths")
error(f$cdf(1e6, 1), "million terms")
stopifnot(identical(f$qf(numeric(), 1), numeric()))
# Underflow in an ordinary probability should not erase a finite log CDF.
near(f$cdf(0, 2000, logp = TRUE), -2000 / sqrt(2))

# Check all fitted derivatives by differencing lower-order analytic quantities.
# The first derivative checks instead difference the independent likelihood.
y <- c(0, 1, 2, 7, 40)
mu <- c(.2, .7, 4, 9, 32)
w <- c(.5, 2, 1, 0, 1.5)
for (phi in c(1 + 1e-5, 1.5, 5, 30)) {
  t <- log(phi - 1); fam <- gp(phi); h <- 1e-5
  d <- fam$Dd(y, mu, t, w, 2)
  mp <- fam$Dd(y, mu + h, t, w, 2)
  mm <- fam$Dd(y, mu - h, t, w, 2)
  tp <- fam$Dd(y, mu, t + h, w, 2)
  tm <- fam$Dd(y, mu, t - h, w, 2)
  near(d$Dmu, -2*w*(ll(y, mu+h, t)-ll(y, mu-h, t))/(2*h), 2e-6)
  near(d$Dth, (fam$dev.resids(y, mu, w, t+h)-fam$dev.resids(y, mu, w, t-h))/(2*h), 2e-6)
  for (pair in list(c("Dmu2", "Dmu"), c("Dmu3", "Dmu2"), c("Dmu4", "Dmu3")))
    near(d[[pair[1]]], (mp[[pair[2]]]-mm[[pair[2]]])/(2*h), 2e-6)
  for (pair in list(c("Dth2", "Dth"), c("Dmuth", "Dmu"),
                   c("Dmu2th", "Dmu2"), c("Dmu3th", "Dmu3"),
                   c("Dmuth2", "Dmuth"), c("Dmu2th2", "Dmu2th")))
    near(d[[pair[1]]], (tp[[pair[2]]]-tm[[pair[2]]])/(2*h), 2e-6)
  sat <- saturated(y, t)
  sp <- saturated(y, t+h); sm <- saturated(y, t-h)
  near(sat[, 2], (sp[, 1]-sm[, 1])/(2*h), 2e-6)
  near(sat[, 3], (sp[, 2]-sm[, 2])/(2*h), 2e-6)
  for (i in which(y > 0)) {
    opt <- optimize(function(eta) ll(y[i], exp(eta), t), c(-20, 20), maximum = TRUE)
    near(sat[i, 1], opt$objective, 1e-7)
  }
  ls <- fam$ls(y, w, t, 1)
  near(ls$ls, sum(w*sat[, 1]))
  near(ls$lsth1, sum(w*sat[, 2]))
  near(ls$lsth2, sum(w*sat[, 3]))
  # Zero-weight observations must not produce 0*Inf or 0*NaN.
  badmu <- mu; badmu[w == 0] <- Inf
  near(fam$dev.resids(y, badmu, w), fam$dev.resids(y, mu, w))
  stopifnot(all(vapply(fam$Dd(y, badmu, t, w, 2), function(x) all(is.finite(x)), logical(1))))
}

# Simulation uses independent branching rather than CDF inversion.
set.seed(461)
for (phi in c(1, 1.01, 2, 10)) {
  z <- gp(phi)$rd(rep(4, 60000))
  stopifnot(all(z >= 0 & z == floor(z)), abs(mean(z)-4) < .09,
            abs(var(z)/(4*phi)-1) < .06)
}

# Fitting and the exact Poisson boundary.
set.seed(12)
dat <- data.frame(x = runif(400), exposure = runif(400, .5, 2))
dat$y <- gp(2)$rd(exp(.4 + sin(6*dat$x)) * dat$exposure)
form <- y ~ s(x, k = 6) + offset(log(exposure))
precise <- gam.control(epsilon = 1e-10, newton = list(conv.tol = 1e-10))
for (method in c("REML", "ML")) {
  p <- gam(form, data = dat, family = gp(1), method = method, control = precise)
  q <- gam(form, data = dat, family = poisson(), method = method, control = precise)
  near(coef(p), coef(q), 1e-5)
  near(logLik(p), logLik(q), 1e-7)
  near(deviance(p), deviance(q), 1e-7)
  near(p$null.deviance, q$null.deviance, 1e-6)
  near(p$sp, q$sp, 1e-5)
  near(p$Vp, q$Vp, 1e-5)
  fit <- gam(form, data = dat, family = gp(), method = method)
  stopifnot(fit$converged, fit$family$getTheta(TRUE) > 1,
            all(is.finite(predict(fit, se.fit = TRUE)$se.fit)), fit$sig2 == 1)
  for (type in c("response", "pearson", "deviance", "working"))
    stopifnot(all(is.finite(residuals(fit, type = type))))
  near(as.numeric(logLik(fit)), sum(ll(dat$y, fitted(fit), fit$family$getTheta())))
}

# Low-count unpenalized REML also exercises negative observed curvature.
low <- gam(y ~ x + offset(log(exposure)), data = dat, family = gp(), method = "REML")
stopifnot(low$converged)

# An independent unpenalized likelihood optimizer tests dispersion estimation.
# Higher counts avoid mgcv's zero-column ML SVD bug (see gp documentation).
set.seed(52)
dat <- data.frame(x = runif(400), exposure = runif(400, .8, 1.2))
dat$y <- gp(2)$rd(exp(4 + .4*dat$x) * dat$exposure)
X <- model.matrix(~x, dat); off <- log(dat$exposure)
objective <- function(par) -sum(ll(dat$y, exp(drop(X %*% par[1:2])+off), par[3]))
fit <- gam(y ~ x + offset(log(exposure)), data = dat, family = gp(), method = "ML")
opt <- optim(c(coef(fit), fit$family$getTheta()), objective, method = "BFGS",
             control = list(reltol = 1e-12))
near(c(coef(fit), fit$family$getTheta()), opt$par, 1e-4)
near(-as.numeric(logLik(fit)), opt$value, 1e-8)
# Fixed dispersion against its own independent optimizer.
fixed <- gam(y ~ x + offset(log(exposure)), data = dat, family = gp(2), method = "ML")
optf <- optim(coef(fixed), function(b) objective(c(b, 0)), method = "BFGS")
near(coef(fixed), optf$par, 1e-5)

# Integer weights agree with replication; fractional weights agree with ML.
dat$w <- rep(0:3, length.out = nrow(dat))
weighted <- gam(y ~ x, data = dat, weights = w, family = gp(), method = "ML")
replicated <- gam(y ~ x, data = dat[rep(seq_len(nrow(dat)), dat$w), ], family = gp(), method = "ML")
near(coef(weighted), coef(replicated), 1e-5)
near(weighted$family$getTheta(TRUE), replicated$family$getTheta(TRUE), 1e-5)
dat$w <- dat$w / 2
fractional <- gam(y ~ x, data = dat, weights = w, family = gp(), method = "ML")
fw <- function(par) -sum(dat$w * ll(dat$y, exp(drop(X %*% par[1:2])), par[3]))
optw <- optim(c(coef(fractional), fractional$family$getTheta()), fw, method = "BFGS")
near(c(coef(fractional), fractional$family$getTheta()), optw$par, 1e-4)

# Reusing a template must not mutate it or a previous fit.
template <- gp(-3)
f1 <- gam(y ~ x, data = dat, family = template, method = "REML")
saved <- f1$family$getTheta()
f2 <- gam(y ~ 1, data = dat, family = template, method = "REML")
near(template$getTheta(TRUE), 3)
near(f1$family$getTheta(), saved)
near(f1$family$getTheta(TRUE), 1+exp(saved))
stopifnot(all(is.finite(f1$family$rd(fitted(f1)))),
          all(is.finite(f1$family$qf(.5, fitted(f1)))))
nointercept <- gam(y ~ x - 1 + offset(log(exposure)), data = dat, family = gp(2), method = "REML")
near(nointercept$null.deviance, sum(nointercept$family$dev.resids(dat$y, dat$exposure, rep(1, nrow(dat)))))
for (bad in list(c(-1, 2), c(.5, 1), c(Inf, 1)))
  error(gam(y ~ 1, data = data.frame(y = bad), family = gp(), method = "REML"), "integer responses")
error(gam(y ~ 1, data = data.frame(y = c(0, 0)), family = gp(), method = "REML"), "positive count")
error(gam(y ~ x, data = dat, family = gp(), method = "REML", scale = 2), "scale is fixed")
# Estimated dispersion can approach the Poisson boundary without switching.
set.seed(44)
pd <- data.frame(y = rpois(500, 4))
edge <- gam(y ~ 1, data = pd, family = gp(), method = "REML")
stopifnot(edge$converged, is.finite(edge$family$getTheta()), edge$family$getTheta(TRUE) >= 1)
cat("Generalized Poisson distribution, derivatives and fitting: passed\n")
