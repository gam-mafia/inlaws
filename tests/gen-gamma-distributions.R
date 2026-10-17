library("inlaws")
near <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(x) stopifnot(inherits(tryCatch(x,error=identity),"error"))
f <- gen_gammalss(); p <- c(0,1e-12,.001,.1,.5,.9,.999,1-1e-12,1)
for (Q in c(-1.5,-.5,-1e-8,0,1e-8,1e-5,1e-4,.5,1,4)) {
  mu <- cbind(rep(2,length(p)),.5,1+.5*Q)
  for (tail in c(FALSE,TRUE)) for (logp in c(FALSE,TRUE)) {
    prob <- if(logp) log(p) else p
    q <- f$qf(prob,mu,lower.tail=tail,log.p=logp)
    near(f$cdf(q,mu,lower.tail=tail),p,1e-9)
    lp <- f$cdf(q,mu,lower.tail=tail,logp=TRUE)
    near(lp[is.finite(log(p))],log(p)[is.finite(log(p))],1e-7)
    stopifnot(q[if(tail) 1L else length(p)] == 0,
              q[if(tail) length(p) else 1L] == Inf)
  }
  near(f$cdf(c(-Inf,-1,0,rep(Inf,length(p)-3)),mu),c(0,0,0,rep(1,length(p)-3)))
}
# A tiny gamma quantile can underflow although Y remains representable.
mu <- cbind(rep(1,3),1,101)
q <- f$qf(c(.1,.5,.9),mu)
stopifnot(all(q>0)); near(f$cdf(q,mu),c(.1,.5,.9),1e-10)
q <- c(.001,.2,1,3,20); n <- length(q); m <- 2; s <- .5
# Named submodels: lognormal, gamma, Weibull and exponential.
mu <- cbind(rep(m,n),s,1)
near(f$cdf(q,mu),plnorm(q,log(m)-s^2/2,s))
near(f$qf(.9,mu),qlnorm(rep(.9,n),log(m)-s^2/2,s))
mu[,3] <- 1+s^2
near(f$cdf(q,mu),pgamma(q,shape=1/s^2,scale=m*s^2))
mu[,3] <- 1+s
near(f$cdf(q,mu),pweibull(q,shape=1/s,scale=m/gamma(1+s)))
mu <- cbind(rep(m,n),1,2)
near(f$cdf(q,mu),pexp(q,rate=1/m))
# Moments and PIT provide independent checks on joint callback conventions.
variance <- getFromNamespace(".gen_gamma_variance","inlaws")
near(variance(cbind(m,s,1)),m^2*expm1(s^2))
near(variance(cbind(m,s,1+s^2)),m^2*s^2)
set.seed(62)
for (Q in c(-.3,0,.5,1,4)) {
  pars <- cbind(rep(m,20000),.4,1+.4*Q); y <- f$rd(pars)
  stopifnot(length(y)==nrow(pars),all(is.finite(y)),all(y>0))
  near(mean(y),m,.03); near(var(y),variance(pars[1,,drop=FALSE]),.07)
  near(mean(f$cdf(y,pars)),.5,.01)
}
d <- f$gratia_diagnostics(cbind(c(1,2,3),.5,c(.3,.5,1)),1,1,NULL,NULL)
near(d$mean,1:3)
stopifnot(d$continuous,!d$variance_available,all(is.infinite(d$variance[1:2])))
near(d$variance[3],9*expm1(.25))
clone <- unserialize(serialize(f,NULL)); near(clone$cdf(q,mu),f$cdf(q,mu))
mu[1,1] <- NA; stopifnot(is.na(f$cdf(q,mu)[1]),is.na(f$qf(.5,mu)[1]),is.na(f$rd(mu)[1]))
for (bad in list(c(1,1,1),matrix(1,2,2),cbind(1,0,1),cbind(1,1,-1),cbind(Inf,1,1))) fails(f$cdf(1,bad))
fails(f$qf(c(.1,.2),matrix(1,3,3)))
fails(gen_gammalss(list("log","identity","log")))
stopifnot(all(is.nan(suppressWarnings(f$qf(-.1,cbind(1,.5,1))))))
for(r in c(.8,1,1.5)) stopifnot(is.nan(suppressWarnings(f$qf(-.1,cbind(1,.5,r)))))
