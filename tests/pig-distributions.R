library(inlaws)
near <- function(x, y, tol = 1e-7) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
            max(abs(x-y)) <= tol*max(1,abs(x),abs(y)))
}
error <- function(expr) stopifnot(inherits(tryCatch(expr, error=identity), "error"))
pmf <- getFromNamespace(".pig_logpmf", "inlaws")
# Independent integration in log mixing-variable coordinates.
for (mu in c(.2, 3, 15)) for (th in c(.05, .5, 3)) for (y in c(0, 1, 6, 30)) {
  integrand <- function(z) {
    x <- exp(z)
    exp(dpois(y, mu*x, log=TRUE) - .5*log(2*pi*th) - .5*z -
          (x-1)^2/(2*th*x))
  }
  p <- integrate(integrand, -30, 20, rel.tol=1e-10, abs.tol=1e-14)$value
  near(exp(pmf(y, mu, log(th))), p, 1e-10)
}
# A scaled Bessel calculation is independent of the probability recurrence.
for (mu in c(.1, 4, 80)) for (th in c(.02, .8, 10)) {
  y <- c(0:10, 50, 150); s <- sqrt(1+2*mu*th)
  lp <- .5*log(2/(pi*th)) + y*log(mu) - (y-.5)*log(s) -
    lgamma(y+1) - 2*mu/(1+s) + log(besselK(s/th, y-.5, expon.scaled=TRUE))
  keep <- is.finite(lp)
  near(pmf(y[keep], mu, log(th)), lp[keep])
}
for (th in c(.01, .3, 2)) {
  y <- 0:3000; p <- exp(pmf(y, 3, log(th)))
  near(sum(p), 1); near(sum(y*p), 3); near(sum((y-3)^2*p), 3+th*9)
}
for (mu in c(.01, 1, 50)) near(pmf(0:80, mu, log(1e-10)), dpois(0:80, mu, log=TRUE))
f <- pig(.6)
y <- 0:100; p <- exp(pmf(y, 3, log(.6)))
near(f$cdf(y, 3), cumsum(p)); near(f$cdf(y+.7, 3), cumsum(p))
near(f$cdf(y, 3, logp=TRUE), log(cumsum(p)))
probs <- c(.001, .1, .5, .9, .999)
q <- f$qf(probs, 3)
stopifnot(all(f$cdf(q,3) >= probs), all(f$cdf(q-1,3) < probs),
          identical(f$qf(c(0,1),3), c(0,Inf)),
          all(f$qf(c(0,.5,1),0) == 0),
          identical(f$cdf(c(-Inf,-1,Inf),3), c(0,0,1)),
          identical(f$cdf(c(-1,0,1),0), c(0,1,1)),
          is.na(f$cdf(NA_real_,3)), is.na(f$qf(.5,NA_real_)))
near(f$cdf(4,c(1,2,3)), vapply(1:3,function(m) f$cdf(4,m),numeric(1)))
error(f$cdf(1:2,1:3)); error(f$qf(.5,-1)); error(f$cdf(1e6,3))
set.seed(91); x <- f$rd(rep(4,150000))
near(mean(x), 4, .015); near(var(x), 4+.6*16, .025)
set.seed(23); x <- f$rd(rep(2,20))
set.seed(23); near(x,f$rd(rep(2,20),wt=100,scale=3))
stopifnot(all(f$rd(rep(0,20)) == 0), length(f$rd(numeric())) == 0)

error(pig(1e308)$cdf(1,1e308))
error(pmf(0,1e308,log(1e308)))
