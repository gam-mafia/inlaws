library("inlaws")
near <- function(a,b,tol=1e-7) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch(expr,error=identity),"error"))
f <- scatlss()
mu <- cbind(c(-2,0,3),c(.2,1,4),c(3.1,8,1e8)); q <- c(-3,.2,8)
for(tail in c(FALSE,TRUE)) for(logp in c(FALSE,TRUE)) {
  near(f$cdf(q,mu,lower.tail=tail,logp=logp),pt((q-mu[,1])/mu[,2],mu[,3],lower.tail=tail,log.p=logp))
  p <- c(.01,.4,.999); if(logp) p <- log(p)
  near(f$qf(p,mu,lower.tail=tail,log.p=logp),mu[,1]+mu[,2]*qt(p,mu[,3],lower.tail=tail,log.p=logp))
  near(f$cdf(f$qf(p,mu,lower.tail=tail,log.p=logp),mu,lower.tail=tail,logp=logp),p)
}
near(f$qf(c(0,.5,1),mu),c(-Inf,0,Inf))
near(f$cdf(c(-Inf,NA,Inf),mu),c(0,NA,1))
near(f$cdf(0,mu),f$cdf(rep(0,3),mu))
p <- mu; p[2,2] <- NA
stopifnot(is.na(f$cdf(0,p)[2]),is.na(f$qf(.5,p)[2]))
# Simulation uses the requested Student-t law, ignoring likelihood weights.
set.seed(88); a <- f$rd(mu)
set.seed(88); near(a,mu[,1]+mu[,2]*rt(3,mu[,3]))
set.seed(88); near(a,f$rd(mu,wt=100,scale=10))
set.seed(97); many <- matrix(rep(c(1.2,.7,9),each=100000),ncol=3)
y <- f$rd(many)
near(mean(y),1.2,.01); near(var(y),.7^2*9/7,.02)
d <- f$gratia_diagnostics(mu,rep(1,3),1,q,NULL)
near(d$mean,mu[,1]); near(d$variance,mu[,2]^2*mu[,3]/(mu[,3]-2))
stopifnot(isTRUE(d$continuous))
# The diagnostics factory binds moments; probability/sampler hooks remain native.
clone <- unserialize(serialize(f,NULL))
near(clone$cdf(q,mu),f$cdf(q,mu)); near(clone$gratia_diagnostics(mu,1,1,q,NULL)$variance,d$variance)
for(bound in list(2,-1,Inf,NA_real_,numeric(),c(3,4),"3")) fails(scatlss(bound))
for(pars in list(cbind(0,0,4),cbind(Inf,1,4),cbind(0,1,2),matrix(1,2,2))) fails(f$cdf(0,pars))
fails(f$cdf(1:2,mu)); fails(f$qf(1:2,mu))
# Deviances retain close-to-location residuals, even when z^2/nu underflows.
for(nu in c(4,1e8,1e308)) {
  z <- c(-1e-9,1e-9,1e-8,1e-6)
  # Independent local-curvature expansion; its relative remainder is
  # bounded by z^2/(2*nu), below 1.3e-13 on this grid.
  expected <- (1+1/nu)*z^2
  actual <- f$dev.resids(z,cbind(rep(0,length(z)),1,nu),1)
  near(actual/expected,rep(1,length(z)),1e-12)
}
# Huge but finite residuals need no representable squared residual.
near(f$dev.resids(1e160,matrix(c(0,1,4),1),1),
     2*(dt(0,4,log=TRUE)-dt(1e160,4,log=TRUE)),1e-12)
# A finite shifted-log component must not produce infinite reported df.
big <- scatlss(1e308); X <- matrix(1,1,3)
attr(X,"lpi") <- list(1L,2L,3L)
trial <- big$ll(0,X,c(0,0,log(1e308)),1,big,deriv=0)
stopifnot(identical(trial$l,-Inf))
