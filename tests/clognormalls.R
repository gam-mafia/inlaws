library(inlaws)

near <- function(a, b, tol = 1e-6) {
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
}
fails <- function(expr) stopifnot(inherits(tryCatch(expr, error = identity), "error"))
evaluate <- getFromNamespace(".clnormal_eval", "inlaws")

# Independent observed-data log likelihood, used below to check the original
# response measure, case weights, and NCV's censor-preserving row selection.
loglik <- function(y, eta, wt = 1) {
  if (!is.matrix(y)) y <- cbind(y, y)
  location <- eta[, 1] - exp(2*eta[, 2])/2
  sigma <- exp(eta[, 2])
  exact <- y[, 1] == y[, 2]
  left <- y[, 2] == -Inf; right <- y[, 2] == Inf
  interval <- !(exact | left | right)
  ans <- numeric(nrow(y))
  ans[exact] <- dlnorm(y[exact, 1], location[exact], sigma[exact], log = TRUE)
  ans[left] <- plnorm(y[left, 1], location[left], sigma[left], log.p = TRUE)
  ans[right] <- plnorm(y[right, 1], location[right], sigma[right], lower.tail = FALSE, log.p = TRUE)
  ans[interval] <- log(plnorm(y[interval, 2], location[interval], sigma[interval]) -
                       plnorm(y[interval, 1], location[interval], sigma[interval]))
  ans * wt
}

# Check every mixed derivative through total order four. In particular the
# new third/fourth scale derivatives were not needed by clognormal().
y <- cbind(c(.7, 1.2, 2, .8, 1.4, exp(10), exp(-10)),
           c(.7, -Inf, Inf, 1.3, 1.40001, Inf, -Inf))
a <- c(.1, -.2, .3, .15, -.3, 0, 0)
t <- log(c(.6, .7, .8, .9, 1.1, .6, .6))
w <- c(1, 2, .5, 3, 1, .7, 1)
r <- evaluate(y, a, t, 4L, predictor = TRUE)
take <- function(r, a, b) r$ell[, which(r$powers[, 1] == a & r$powers[, 2] == b)] * factorial(a)*factorial(b)
near(r$ell[, 1], loglik(y, cbind(a, t)), 1e-6)
h <- 2e-5
ap <- evaluate(y, a+h, t, 4L, predictor = TRUE)
am <- evaluate(y, a-h, t, 4L, predictor = TRUE)
tp <- evaluate(y, a, t+h, 4L, predictor = TRUE)
tm <- evaluate(y, a, t-h, 4L, predictor = TRUE)
for (k in 1:4) for (b in 0:k) {
  aa <- k-b
  if (aa > 0) near(take(r, aa, b), (take(ap, aa-1, b)-take(am, aa-1, b))/(2*h), 1e-4)
  if (b > 0) near(take(r, aa, b), (take(tp, aa, b-1)-take(tm, aa, b-1))/(2*h), 1e-4)
}

fam <- clognormalls()
# Coefficient gradients/Hessians and the direct-eta NCV path, with offsets.
X <- cbind(1, seq(-1, 1, length.out = nrow(y)), 1)
attr(X, "lpi") <- list(1:2, 3L)
beta <- c(.1, -.2, log(.8))
off <- list(rep(.2, nrow(y)), seq(-.1, .1, length.out = nrow(y)))
eta <- cbind(drop(X[,1:2] %*% beta[1:2])+off[[1]], beta[3]+off[[2]])
ll <- fam$ll(y, X, beta, w, fam, offset = off, deriv = 1, ncv = TRUE)
near(ll$l, sum(loglik(y, eta, w)), 1e-6)
near(ll$l, fam$ll(y, X, beta, w, fam, offset = off, eta = eta)$l)
for (k in seq_along(beta)) {
  bp <- bm <- beta; bp[k] <- bp[k]+h; bm[k] <- bm[k]-h
  lp <- fam$ll(y, X, bp, w, fam, offset = off, deriv = 1)
  lm <- fam$ll(y, X, bm, w, fam, offset = off, deriv = 1)
  near(ll$lb[k], (lp$l-lm$l)/(2*h), 1e-5)
  near(ll$lbb[,k], (lp$lb-lm$lb)/(2*h), 1e-5)
}

set.seed(17)
n <- 240L
dat <- data.frame(x = runif(n), z = runif(n), off1 = runif(n, -.1, .1),
                  off2 = runif(n, -.05, .05))
true_mu <- exp(.4 + sin(6*dat$x) + dat$off1)
true_sigma <- exp(-1 + .6*dat$z + dat$off2)
raw <- rlnorm(n, log(true_mu)-true_sigma^2/2, true_sigma)
dat$y <- raw
# Constant scale/no smooths: compare against ordinary log-Gaussian likelihood
# optimization, and against clognormal's ML family-parameter estimate.
constant <- gam(list(y ~ x + offset(off1), ~ 1), data = dat,
                family = clognormalls(), method = "REML")
extended <- gam(y ~ x + offset(off1), data = dat,
                family = clognormal(), method = "ML")
near(fitted(constant)[,1], fitted(extended), 1e-4)
near(fitted(constant)[,2], rep(extended$family$getTheta(TRUE), n), 1e-4)
near(as.numeric(logLik(constant)), as.numeric(logLik(extended)), 1e-5)

# All censoring types, with varying scale and offsets in both predictors.
yc <- cbind(raw, raw)
kind <- seq_len(n) %% 4L
left <- kind == 1L & raw < 1
right <- kind == 2L & raw > 2
interval <- kind == 3L
stopifnot(any(left), any(right), any(interval))
yc[left, ] <- cbind(rep(1, sum(left)), rep(-Inf, sum(left)))
yc[right, ] <- cbind(rep(2, sum(right)), rep(Inf, sum(right)))
yc[interval, ] <- cbind(raw[interval]*.9, raw[interval]*1.1)
dat$y <- yc
forms <- list(y ~ s(x, k = 6) + offset(off1), ~ s(z, k = 5) + offset(off2))
nei <- list(a = seq_len(n), ma = seq.int(4L, n, 4L),
            d = as.vector(apply(matrix(seq_len(n), 4L), 2L, rev)), md = seq.int(4L, n, 4L))
check_fit <- function(fit, method) {
  stopifnot(identical(fit$method, method),
            identical(fit$outer.info$conv, "full convergence"),
            all(is.finite(coef(fit))), all(is.finite(fit$gcv.ubre)),
            identical(dim(fitted(fit)), c(n, 2L)), all(fitted(fit) > 0))
  pr <- predict(fit, type = "response", se.fit = TRUE)
  lp <- predict(fit, se.fit = TRUE)
  near(pr$fit, fitted(fit))
  near(pr$fit, exp(lp$fit))
  near(pr$se.fit, pr$fit*lp$se.fit)
  near(as.numeric(logLik(fit)), sum(loglik(yc, lp$fit)), 1e-6)
  stopifnot(all(is.finite(residuals(fit))),
    all(is.na(residuals(fit, type = "response")[yc[,1] != yc[,2]])))
  stopifnot(is.finite(summary(fit)$dev.expl))
}
reml <- gam(forms, data = dat, family = clognormalls(), method = "REML")
check_fit(reml, "REML")
loo <- gam(forms, data = dat, family = clognormalls(), method = "NCV")
check_fit(loo, "NCV")
grouped <- gam(forms, data = dat, family = clognormalls(), method = "NCV", nei = nei)
check_fit(grouped, "NCV")
# Different censor patterns must survive gamlss.ncv's reordered response
# subsets. Re-evaluate the reported CV predictors with an independent likelihood.
G <- gam(forms, data = dat, family = clognormalls(), fit = FALSE)
X <- predict(reml, type = "lpmatrix")
beta <- coef(reml)
ll <- fam$ll(yc, X, beta, rep(1, n), fam, offset = G$offset, deriv = 1, ncv = TRUE)
ll$gamma <- 1
cvfam <- fam; cvfam$qapprox <- FALSE
H <- -ll$lbb + diag(1, length(beta))
cv <- cvfam$ncv(X, yc, rep(1, n), nei, beta, cvfam, ll, R = chol(H), offset = G$offset)
cveta <- attr(cv$NCV, "eta.cv")
near(as.numeric(cv$NCV), -sum(loglik(yc[nei$d, , drop = FALSE], cveta)), 1e-6)
# gamma's full-data score adjustment must include offsets as well.
ll_gamma <- ll; ll_gamma$gamma <- 1.3
cv_gamma <- cvfam$ncv(X, yc, rep(1, n), nei, beta, cvfam, ll_gamma,
                      R = chol(H), offset = G$offset)
near(as.numeric(cv_gamma$NCV),
     -1.3*sum(loglik(yc[nei$d, , drop = FALSE], cveta)) +
       .3*sum(loglik(yc[nei$d, , drop = FALSE], predict(reml)[nei$d, ])), 1e-6)

# Differentiate the NCV score along a coefficient direction, including the
# corresponding change in the Hessian used for neighbourhood updates.
direction <- matrix(seq(-.05, .05, length.out = length(beta)), ncol = 1)
ll3 <- fam$ll(yc, X, beta, rep(1, n), fam, offset = G$offset,
              deriv = 3, d1b = direction, ncv = TRUE)
ll3$gamma <- 1
cv3 <- cvfam$ncv(X, yc, rep(1, n), nei, beta, cvfam, ll3, R = chol(H),
                 offset = G$offset, dH = ll3$d1H, db = direction, deriv = TRUE)
score <- function(step) {
  b <- beta + step*direction[, 1]
  lr <- fam$ll(yc, X, b, rep(1, n), fam, offset = G$offset, deriv = 1, ncv = TRUE)
  lr$gamma <- 1
  cvfam$ncv(X, yc, rep(1, n), nei, b, cvfam, lr,
    R = chol(-lr$lbb + diag(1, length(beta))), offset = G$offset)$NCV
}
near(cv3$NCV1, (score(1e-4) - score(-1e-4))/2e-4, 1e-5)

# Frequency weights are likelihood weights, including zero-weight rows.
dat$w <- rep(c(0, 1, 2), length.out = n)
weighted <- gam(list(y ~ x, ~ z), data = dat, weights = w,
                family = clognormalls(), method = "REML")
repdat <- dat[rep(seq_len(n), dat$w), ]
replicated <- gam(list(y ~ x, ~ z), data = repdat, family = clognormalls(), method = "REML")
near(coef(weighted), coef(replicated), 1e-5)
near(as.numeric(logLik(weighted)), as.numeric(logLik(replicated)), 1e-5)
# Independently maximize the censored likelihood with both parameters varying.
objective <- function(b) -sum(loglik(yc,
  cbind(b[1] + b[2]*dat$x, b[3] + b[4]*dat$z), dat$w))
opt <- optim(as.numeric(coef(weighted)) + .01, objective, method = "BFGS",
             control = list(reltol = 1e-12))
near(coef(weighted), opt$par, 1e-4)
near(as.numeric(logLik(weighted)), -opt$value, 1e-6)
# Matrix responses keep their censoring state after missing-value removal,
# subsetting and serialization, including predictions without a response.
dat$x[1] <- NA
missing <- gam(list(y ~ x, ~ z), data = dat, family = clognormalls(),
               method = "REML", na.action = na.exclude)
stopifnot(nobs(missing) == n-1L, length(residuals(missing)) == n, is.na(residuals(missing)[1]))
subset_fit <- gam(list(y ~ x, ~ z), data = dat, family = clognormalls(),
                  method = "REML", subset = seq_len(n) > 4L)
stopifnot(nobs(subset_fit) == n-4L)
copy <- unserialize(serialize(reml, NULL))
newdata <- dat[2:8, c("x", "z", "off1", "off2")]
near(predict(copy, newdata, type = "response"), predict(reml, newdata, type = "response"))

# Distribution components use row-specific means and sigma, and retain the
# upper-tail/log-probability interface of clognormal().
mu <- fitted(reml)
for (lower in c(TRUE, FALSE)) {
  quant <- fam$qf(-20, mu, lower.tail = lower, log.p = TRUE)
  near(fam$cdf(quant, mu, lower.tail = lower, logp = TRUE), rep(-20, n))
}
near(fam$qf(.5, mu), exp(log(mu[,1])-mu[,2]^2/2))
set.seed(10)
actual <- reml$family$rd(mu)
set.seed(10)
near(actual, rlnorm(n, log(mu[, 1]) - mu[, 2]^2/2, mu[, 2]))

fails(clognormalls(link = "log"))
fails(clognormalls(link = list("log", "identity")))
fails(gam(list(y ~ x, ~ z), data = dat, weights = rep(-1, n), family = clognormalls()))
fails(gam(list(y ~ x, ~ z), data = dat, weights = rep(0, n), family = clognormalls()))
fails(gam(list(y ~ 1, ~ 1, 1+2 ~ x-1), data = dat, family = clognormalls()))
fails(bam(forms, data = dat, family = clognormalls()))
