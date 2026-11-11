library("inlaws")

near <- function(x, y, tol = 1e-7) {
  stopifnot(all(abs(x - y) <= tol * pmax(1, abs(x), abs(y))))
}
differentiate <- function(f, x, h = 1e-5) (f(x + h) - f(x - h)) / (2 * h)
eval_weibull <- getFromNamespace(".weibull_eval", "inlaws")

# Column order within degree d is (d,0), (d-1,1), ..., (0,d).
y <- c(0.8, 1.8, 2.1, 1, 1)
eta <- c(0.3, 0.4, 0.9, 0.4, 0.002)
theta <- log(c(.4, 1, 3, .00009, 102))
r <- eval_weibull(y, eta, theta, 4L)
ok <- c(1:3, 5L)
near(r[ok, 1L], dweibull(y[ok], exp(theta[ok]),
  exp(eta[ok] - lgamma(1 + exp(-theta[ok]))), log = TRUE))
for (degree in 1:4) {
  old <- (degree - 1L) * degree / 2L + seq_len(degree)
  current <- degree * (degree + 1L) / 2L + seq_len(degree + 1L)
  near(r[, current[-length(current)]], differentiate(
    function(x) eval_weibull(y, x, theta, degree - 1L)[, old, drop = FALSE], eta, h = 1e-6),
    tol = 8e-5)
  near(r[, current[-1L]], differentiate(
    function(x) eval_weibull(y, eta, x, degree - 1L)[, old, drop = FALSE], theta),
    tol = 5e-5)
}

# The extended-family interface differentiates in mu and log(k).
f <- weibull()
y <- c(1, 1.4, 3.8, 2.1)
mu <- c(1.4, 1.2, 1.8, 3)
w <- c(0, 2, .5, 1)
th <- log(1.6)
d <- function(mu, th) f$dev.resids(y, mu, w, th)
dd <- function(mu, th) f$Dd(y, mu, th, w, 2L)
z <- dd(mu, th)
near(z$Dmu, differentiate(function(m) d(m, th), mu))
near(z$Dth, differentiate(function(t) d(mu, t), th))
for (pair in list(c("Dmu", "Dmu2"), c("Dth", "Dmuth"), c("Dmu2", "Dmu3"),
                 c("Dmuth", "Dmu2th"), c("Dmu3", "Dmu4"), c("Dmuth2", "Dmu2th2"),
                 c("Dmu2th", "Dmu3th"), c("EDmu2", "EDmu3")))
  near(z[[pair[2L]]], differentiate(function(m) dd(m, th)[[pair[1L]]], mu))
for (pair in list(c("Dth", "Dth2"), c("Dmu", "Dmuth"), c("Dmu2", "Dmu2th"),
                 c("Dmuth", "Dmuth2"), c("Dmu2th", "Dmu2th2"),
                 c("Dmu3", "Dmu3th"), c("EDmu2", "EDmu2th")))
  near(z[[pair[2L]]], differentiate(function(t) dd(mu, t)[[pair[1L]]], th))
stopifnot(all(vapply(z, function(a) a[1L] == 0, logical(1))))
sat <- function(t) f$ls(y, w, t, 1)
near(sat(th)$lsth1, differentiate(function(t) sat(t)$ls, th))
near(sat(th)$lsth2, differentiate(function(t) sat(t)$lsth1, th))
near(sat(th)$ls - sum(d(mu, th)) / 2,
     sum(w * dweibull(y, exp(th), mu / gamma(1 + exp(-th)), log = TRUE)))

# At extreme predictors the conventional scale may underflow while the
# likelihood stays finite. A far-left hazard can also underflow while its
# fourth mean derivative remains representable.
k <- c(1e-200, 1e100)
eta <- c(600, 8e-98)
r <- eval_weibull(c(1, 1), eta, log(k), 4L)
stopifnot(all(is.finite(r)))
near(r[1, 1] * k[1], -exp(-1))
near(r[2, 11], -exp(-800 - 0.5772156649015329 + 4 * log(k[2])), tol = 1e-11)
# Tiny nonzero deviance near saturation survives expm1(t) - t cancellation.
near(weibull(1)$dev.resids(1 + 1e-12, 1, 1) / (1e-12)^2, 1, tol = 2e-4)

# The extended-family fourth mean derivative can be representable even
# when mu^4 overflows. Rescaling y and mu together leaves t invariant.
mu <- 1e90
th <- log(1e30)
ordinary <- f$Dd(1, 1, th, 1, 2L)$Dmu4
scaled <- f$Dd(mu, mu, th, 1, 2L)$Dmu4
near(scaled / (ordinary / mu / mu / mu / mu), 1)
