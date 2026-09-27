library(inlaws)
near <- function(x, y, tol = 2e-5) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
            max(abs(x - y)) <= tol * max(1, abs(x), abs(y)))
}
diff5 <- function(f, x, j, h = 1e-4) {
  e <- numeric(length(x)); e[j] <- h
  (f(x - 2*e) - 8*f(x-e) + 8*f(x+e) - f(x+2*e)) / (12*h)
}
jet <- getFromNamespace(".pig_lljet", "inlaws")
pmf <- getFromNamespace(".pig_logpmf", "inlaws")
saturated <- getFromNamespace(".pig_saturated", "inlaws")
idx <- getFromNamespace(".pig_index", "inlaws")
col <- function(p, q) (p+q)*(p+q+1)/2 + q+1
# Differentiate each preceding analytic derivative, in relative mean units.
grid <- expand.grid(y = c(0, 1, 5, 40), mu = c(.03, 2, 30),
                    t = log(c(1e-8, .03, 1, 20)))
for (i in seq_len(nrow(grid))) {
  y <- grid$y[i]; mu <- grid$mu[i]; t <- grid$t[i]
  z <- jet(y, mu, t)
  near(z[, 1], pmf(y, mu, t), 1e-10)
  for (k in 2:15) {
    p <- idx[k, 1]; q <- idx[k, 2]
    j <- if (p > 0) 1L else 2L
    pp <- p - (j == 1); qq <- q - (j == 2)
    f <- function(x) jet(y, mu*x[1], t+x[2])[, col(pp, qq)] *
      factorial(pp)*factorial(qq)*mu^pp
    near(diff5(f, c(1, 0), j), z[, k]*factorial(p)*factorial(q)*mu^p)
  }
}
# Saturated envelope derivatives include movement of the optimizing mean.
for (t in log(c(.001, .3, 4))) {
  s <- saturated(c(0, 1, 7, 30), t)
  near(s[, 2], diff5(function(x) saturated(c(0, 1, 7, 30), x)[, 1], t, 1))
  near(s[, 3], diff5(function(x) saturated(c(0, 1, 7, 30), x)[, 2], t, 1))
}
# Full extended-family callback contracts, including fractional/zero weights.
y <- c(0, 1, 4, 12); mu <- c(.2, 2, 3, 8); w <- c(.5, 0, 2, 1)
f <- pig(.7); t <- f$getTheta()
dd <- f$Dd(y, mu, t, w, 2)
near(dd$Dmu, diff5(function(z) f$dev.resids(y, mu+z, w, t), 0, 1))
for (nm in c("Dmu", "Dmu2", "Dmuth", "Dmu2th")) {
  target <- c(Dmu="Dmu2", Dmu2="Dmu3", Dmuth="Dmu2th", Dmu2th="Dmu3th")[[nm]]
  near(dd[[target]], diff5(function(z) f$Dd(y, mu+z, t, w, 2)[[nm]], 0, 1))
}
near(dd$Dmu4, diff5(function(z) f$Dd(y, mu+z, t, w, 2)$Dmu3, 0, 1))
for (nm in c("Dth", "Dmu", "Dmuth", "Dmu2", "Dmu2th")) {
  target <- c(Dth="Dth2", Dmu="Dmuth", Dmuth="Dmuth2", Dmu2="Dmu2th", Dmu2th="Dmu2th2")[[nm]]
  near(dd[[target]], diff5(function(z) f$Dd(y, mu, z, w, 2)[[nm]], t, 1))
}
near(dd$Dth, diff5(function(z) f$dev.resids(y, mu, w, z), t, 1))
ls <- f$ls(y, w, t, 1)
near(ls$lsth1, diff5(function(z) f$ls(y, w, z, 1)$ls, t, 1))
near(ls$lsth2, diff5(function(z) f$ls(y, w, z, 1)$lsth1, t, 1))
near(sum(ls$LSTH1), ls$lsth1)
near(ls$ls - sum(f$dev.resids(y, mu, w, t))/2, sum(w*pmf(y, mu, t)))
near(dd$EDmu3, diff5(function(z) f$Dd(y, mu+z, t, w, 1)$EDmu2, 0, 1))
near(dd$EDmu2th, diff5(function(z) f$Dd(y, mu, z, w, 1)$EDmu2, t, 1))
stopifnot(all(dd$EDmu2[w > 0] > 0), all(vapply(dd, function(x) x[2] == 0, logical(1))))
# Zero-weight invalid trial means do not contaminate any callback.
mu[2] <- Inf
stopifnot(all(is.finite(unlist(f$Dd(y, mu, t, w, 2)))),
          all(is.finite(f$dev.resids(y, mu, w))))
# Leading nonzero dispersion derivatives survive close to Poisson.
for (th in c(1e-6, 1e-10, 1e-14)) {
  z <- jet(0, 2, log(th))
  near(z[, 3]/th, 2, 1e-5)
  near(2*z[, 6]/th, 2, 1e-5)
}
