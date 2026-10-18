library("inlaws")
near <- function(a, b, tol = 1e-5) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
fails <- function(expr) stopifnot(inherits(tryCatch(expr, error = identity), "error"))
loglik <- function(y, eta, wt = 1) {
  mu <- plogis(eta[, 1])
  phi <- exp(eta[, 2])
  wt * dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE)
}
set.seed(29)
n <- 240L
dat <- data.frame(x = runif(n), z = runif(n), off1 = runif(n, -.1, .1), off2 = runif(n, -.1, .1))
pars <- cbind(plogis(.4 + sin(6 * dat$x) + dat$off1), exp(2 + .8 * dat$z + dat$off2))
fam <- betals()
dat$y <- fam$rd(pars)
forms <- list(y ~ s(x, k = 6) + offset(off1), ~ s(z, k = 5) + offset(off2))
nei <- list(a = seq_len(n), ma = seq.int(4L, n, 4L), d = as.vector(apply(matrix(seq_len(n), 4L), 2L, rev)), md = seq.int(4L, n, 4L))
check <- function(f, method) {
  stopifnot(
    identical(f$method, method), identical(f$outer.info$conv, "full convergence"),
    all(is.finite(coef(f))), all(is.finite(f$gcv.ubre)), identical(dim(fitted(f)), c(n, 2L)),
    all(fitted(f)[, 1] > 0 & fitted(f)[, 1] < 1), all(fitted(f)[, 2] > 0)
  )
  p <- predict(f, type = "response", se.fit = TRUE)
  l <- predict(f, se.fit = TRUE)
  near(p$fit, fitted(f))
  near(p$fit[, 1], plogis(l$fit[, 1]))
  near(p$fit[, 2], exp(l$fit[, 2]))
  near(p$se.fit[, 1], p$fit[, 1] * (1 - p$fit[, 1]) * l$se.fit[, 1])
  near(p$se.fit[, 2], p$fit[, 2] * l$se.fit[, 2])
  near(logLik(f), sum(loglik(dat$y, l$fit)), 1e-7)
  stopifnot(all(is.finite(residuals(f))), is.finite(summary(f)$dev.expl))
}
reml <- gam(forms, data = dat, family = fam, method = "REML")
check(reml, "REML")
bfgs <- gam(forms, data = dat, family = betals(), method = "REML", optimizer = c("outer", "bfgs"))
check(bfgs, "REML")
near(coef(reml), coef(bfgs), 1e-3)
loo <- gam(forms, data = dat, family = betals(), method = "NCV")
check(loo, "NCV")
grouped <- gam(forms, data = dat, family = betals(), method = "NCV", nei = nei)
check(grouped, "NCV")
# NCV offset/subset handling and derivatives, independently recomputed.
G <- gam(forms, data = dat, family = betals(), fit = FALSE)
X <- predict(reml, type = "lpmatrix")
beta <- coef(reml)
yy <- dat$y
ll <- fam$ll(yy, X, beta, rep(1, n), fam, offset = G$offset, deriv = 1, ncv = TRUE)
ll$gamma <- 1
cvfam <- fam
cvfam$qapprox <- FALSE
H <- -ll$lbb + diag(1, length(beta))
cv <- cvfam$ncv(X, yy, rep(1, n), nei, beta, cvfam, ll, R = chol(H), offset = G$offset)
cveta <- attr(cv$NCV, "eta.cv")
near(cv$NCV, -sum(loglik(yy[nei$d], cveta)), 1e-7)
llg <- ll
llg$gamma <- 1.3
cvg <- cvfam$ncv(X, yy, rep(1, n), nei, beta, cvfam, llg, R = chol(H), offset = G$offset)
near(cvg$NCV, -1.3 * sum(loglik(yy[nei$d], cveta)) + .3 * sum(loglik(yy[nei$d], predict(reml)[nei$d, ])), 1e-7)
direction <- matrix(seq(-.05, .05, length.out = length(beta)), ncol = 1)
ll3 <- fam$ll(yy, X, beta, rep(1, n), fam, offset = G$offset, deriv = 3, d1b = direction, ncv = TRUE)
ll3$gamma <- 1
cv3 <- cvfam$ncv(X, yy, rep(1, n), nei, beta, cvfam, ll3, R = chol(H), offset = G$offset, dH = ll3$d1H, db = direction, deriv = TRUE)
score <- function(step) {
  b <- beta + step * direction[, 1]
  lr <- fam$ll(yy, X, b, rep(1, n), fam, offset = G$offset, deriv = 1, ncv = TRUE)
  lr$gamma <- 1
  cvfam$ncv(X, yy, rep(1, n), nei, b, cvfam, lr, R = chol(-lr$lbb + diag(1, length(beta))), offset = G$offset)$NCV
}
near(cv3$NCV1, (score(1e-4) - score(-1e-4)) / 2e-4, 1e-5)
# Sandwich must be the cross-product of the observation scores.
sm <- X
jj <- attr(X, "lpi")
for (k in 1:2) sm[, jj[[k]]] <- X[, jj[[k]], drop = FALSE] * ll$l1[, k]
near(fam$sandwich(yy, X, beta, rep(1, n), fam, G$offset), crossprod(sm), 1e-7)
# Replication weights and a direct unpenalized maximum likelihood reference.
dat$w <- rep(c(0, 1, 2), length.out = n)
a <- gam(list(y ~ x, ~z), data = dat, weights = w, family = betals(), method = "REML")
b <- gam(list(y ~ x, ~z), data = dat[rep(seq_len(n), dat$w), ], family = betals(), method = "REML")
near(coef(a), coef(b))
near(logLik(a), logLik(b))
opt <- optim(coef(a) + .01, function(b) -sum(loglik(yy, cbind(b[1] + b[2] * dat$x, b[3] + b[4] * dat$z), dat$w)), method = "BFGS", control = list(reltol = 1e-12))
near(coef(a), opt$par, 1e-4)
near(logLik(a), -opt$value, 1e-6)
# Constant precision uses the same likelihood as betar().
cst <- gam(list(y ~ x, ~1), data = dat, family = betals(), method = "REML")
mu <- fitted(cst)
f <- mgcv::betar(theta = mu[1, 2])
near(f$aic(yy, mu[, 1], wt = rep(1, n), dev = 0), -2 * as.numeric(logLik(cst)))
# Residuals and diagnostics agree with beta moments.
v <- mu[, 1] * (1 - mu[, 1]) / (1 + mu[, 2])
near(residuals(cst, type = "pearson"), (yy - mu[, 1]) / sqrt(v))
near(residuals(cst, type = "response"), yy - mu[, 1])
near(sum(residuals(cst)^2), deviance(cst))
d <- cst$family$gratia_diagnostics(mu, rep(1, n), 1, yy, predict(cst))
near(d$mean, mu[, 1])
near(d$variance, v)
stopifnot(d$continuous)
dat$x[1] <- NA
missing <- gam(list(y ~ x, ~z), data = dat, family = betals(), method = "REML", na.action = na.exclude)
stopifnot(nobs(missing) == n - 1L, length(residuals(missing)) == n, is.na(residuals(missing)[1]))
copy <- unserialize(serialize(reml, NULL))
nd <- dat[2:8, c("x", "z", "off1", "off2")]
near(predict(copy, nd, type = "response"), predict(reml, nd, type = "response"))
for (bad_weight in list(rep(-1, n), rep(0, n))) fails(gam(forms, data = dat, weights = bad_weight, family = betals()))
fails(gam(list(y ~ 1, ~1, 1 + 2 ~ z - 1), data = dat, family = betals()))
fails(bam(forms, data = dat, family = betals()))
# Invalid data and links are rejected explicitly, including endpoints.
for (bad in list(c(0, .5), c(.5, 1), c(-.1, .5), c(.5, Inf), matrix(.5, 2, 2))) {
  fails(gam(list(y ~ 1, ~1), data = list(y = bad), family = betals()))
}
for (bad in list(list("log", "log"), list("logit"), list("logit", NA), list(1, "log"))) fails(betals(bad))
fails(gam(forms, data = dat, family = betals(), method = "QNCV"))
fails(gam(y ~ x, data = dat, family = betals()))
# Distribution functions agree with base R, including log upper tails.
par <- cbind(c(.02, .3, .8), c(.2, 8, 100))
a <- par[, 1] * par[, 2]
b <- (1 - par[, 1]) * par[, 2]
for (lower in c(TRUE, FALSE)) {
  for (logp in c(TRUE, FALSE)) {
    p <- c(.001, .3, .999)
    if (logp) p <- log(p)
    near(fam$qf(p, par, wt = 99, scale = 2, lower.tail = lower, log.p = logp), qbeta(p, a, b, lower.tail = lower, log.p = logp))
    q <- c(.01, .4, .95)
    near(fam$cdf(q, par, wt = 99, scale = 2, lower.tail = lower, logp = logp), pbeta(q, a, b, lower.tail = lower, log.p = logp))
  }
}
near(fam$cdf(c(-Inf, 0, Inf), par), c(0, 0, 1))
set.seed(17)
near(fam$rd(par), {
  set.seed(17)
  rbeta(3, a, b)
})
for (bad in list(c(.5, 3), matrix(c(0, 3), 1), matrix(c(.5, -1), 1))) fails(fam$rd(bad))
set.seed(19)
draws <- fam$rd(cbind(rep(.3, 30000), rep(8, 30000)))
near(mean(draws), .3, .006)
near(var(draws), .3 * .7 / 9, .04)
# Saturation holds precision fixed and generally differs from setting mu=y.
getdev <- getFromNamespace(".betals_deviance", "inlaws")
y <- c(.02, .2, .5, .8, .98)
mu <- cbind(rep(.4, 5), c(.1, 2, 10, 100, 1e4))
dev <- getdev(y, mu, rep(1, 5))
reference <- vapply(seq_along(y), function(i) {
  phi <- mu[i, 2]
  sat <- optimize(function(a) dbeta(y[i], plogis(a) * phi, plogis(-a) * phi, log = TRUE), c(-25, 25), maximum = TRUE, tol = 1e-10)$objective
  2 * (sat - dbeta(y[i], mu[i, 1] * phi, (1 - mu[i, 1]) * phi, log = TRUE))
}, numeric(1))
near(dev, reference, 1e-7)
