library("inlaws")
near <- function(a,b,tol=1e-7) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch(expr,error=identity),"error"))
moments <- getFromNamespace(".hp_moments","inlaws")
mass <- getFromNamespace(".hp_logmass","inlaws")
f <- hyper_poissonls()

# Independent wide-support normalizer and uniroot inversion, without family numerics.
reference <- function(mu,sigma) {
  k <- 0:4000
  probs <- function(a) {
    lw <- k*a+lgamma(sigma)-lgamma(sigma+k)
    w <- exp(lw-max(lw)); w/sum(w)
  }
  a <- uniroot(function(a) sum(k*probs(a))-mu,
               c(-50,log((mu+1)*(sigma+1))+5),tol=1e-13)$root
  p <- probs(a)
  list(a=a,p=p,mean=sum(k*p),variance=sum((k-mu)^2*p))
}
grid <- expand.grid(mu=c(.01,.5,2,20,100),sigma=c(.01,.3,1,3,100))
z <- moments(grid$mu,grid$sigma,f$control)
for(i in seq_len(nrow(grid))) {
  m <- grid$mu[i]; s <- grid$sigma[i]; ref <- reference(m,s); v <- z[[i]]
  near(v$a,ref$a,1e-9); near(v$mean,m,1e-9); near(v$variance,ref$variance,1e-8)
  near(v$variance,exp(v$a)+(exp(v$a)-s+1)*m-m^2,1e-7)
  yy <- 0:200
  pm <- exp(vapply(yy,function(y) mass(y,v),numeric(1)))
  near(pm,ref$p[yy+1L],1e-9)
  use <- which(pm[-1]>1e-280 & pm[-length(pm)]>1e-280)
  near((pm[-1]/pm[-length(pm)])[use],(exp(v$a)/(s+head(yy,-1L)))[use],1e-8)
  pars <- cbind(rep(m,7),s); q <- c(-1,0,.5,2,10,100,Inf)
  expected <- vapply(q,function(q) sum(ref$p[(0:4000)<=q]),numeric(1))
  near(f$cdf(q,pars),expected,1e-9)
  near(f$cdf(q,pars)+f$cdf(q,pars,lower.tail=FALSE),rep(1,7),1e-12)
}
# Exact Poisson reductions include tiny tails beyond the normalizer support.
for(m in c(1e-20,1e-16,1e-10)) {
  zz<-moments(m,1,f$control)[[1]]
  near(mass(0,zz)/(-m),1,1e-12)
  near(f$dev.resids(0,matrix(c(m,1),1),1)/(2*m),1,1e-12)
}
# A large dispersion with small mean needs a small geometric-scale window.
large<-moments(2,1e12,f$control)[[1]]
near(large$mean,2,1e-9);near(large$variance,6,1e-8)
stopifnot(length(large$y)<1000L)
pars <- cbind(c(.01,2,100),1)
for(lower in c(TRUE,FALSE)) for(logp in c(TRUE,FALSE)) {
  near(f$cdf(c(5,100,300),pars,lower.tail=lower,logp=logp),
       ppois(c(5,100,300),pars[,1],lower.tail=lower,log.p=logp),1e-9)
  p <- c(.01,.4,.99999); if(logp) p <- log(p)
  near(f$qf(p,pars,lower.tail=lower,log.p=logp),
       qpois(p,pars[,1],lower.tail=lower,log.p=logp))
}
near(f$qf(c(-1000,-100,-10),pars,lower.tail=FALSE,log.p=TRUE),
     qpois(c(-1000,-100,-10),pars[,1],lower.tail=FALSE,log.p=TRUE))
# At CDF steps a strict approximate comparison can advance the quantile by one.
# Poisson callbacks must agree with native routines, also at tiny log tails.
for(m in c(.01,.5,2,20,100)) for(lower in c(TRUE,FALSE)) for(logp in c(TRUE,FALSE)) {
  k<-0:200;pa<-cbind(rep(m,length(k)),1)
  p<-ppois(k,m,lower.tail=lower,log.p=logp)
  near(f$cdf(k,pa,lower.tail=lower,logp=logp),p,1e-13)
  near(f$qf(p,pa,lower.tail=lower,log.p=logp),
       qpois(p,m,lower.tail=lower,log.p=logp),1e-13)
}
for(lower in c(TRUE,FALSE)) {
  k<-0:100;pa<-cbind(rep(20,length(k)),20)
  p<-f$cdf(k,pa,lower.tail=lower,logp=TRUE)
  use<-p>log(.001)&p<log(.99)
  near(f$qf(p[use],pa[use,,drop=FALSE],lower.tail=lower,log.p=TRUE),k[use])
}
for(s in c(.1,1,10)) {
  p <- c(.0001,.1,.5,.9,.9999); pars <- cbind(rep(3,5),s)
  q <- f$qf(p,pars)
  stopifnot(all(f$cdf(q,pars)>=p-1e-10),all(f$cdf(q-1,pars)<p+1e-10))
}
pars <- cbind(c(0,2,3),c(.1,1,3))
near(f$qf(c(1,0,1),pars),c(0,0,Inf))
near(f$cdf(c(-Inf,NA,Inf),pars),c(0,NA,1))
near(f$cdf(0,pars),f$cdf(rep(0,3),pars))
missing <- pars; missing[2,2] <- NA
stopifnot(is.na(f$cdf(1,missing)[2]),is.na(f$qf(.5,missing)[2]))
stopifnot(length(f$cdf(numeric(),pars))==0L,length(f$qf(.5,matrix(numeric(),0,2)))==0L)
for(p in list(cbind(-1,1),cbind(1,0),cbind(Inf,1),matrix(1,2,3))) fails(f$cdf(0,p))
fails(f$cdf(1:2,pars)); fails(f$qf(2,pars)); fails(f$qf(.1,pars,log.p=TRUE))
fails(hyper_poissonls(link=list("identity","log")))
for(ctl in list(list(bad=1),list(sum_tol=1),list(max_terms=1.2),list(mean_tol=NA)))
  fails(hyper_poissonls(control=ctl))
fails(moments(2,3,hyper_poissonls(control=list(max_iter=1))$control))
fails(moments(100,1,hyper_poissonls(control=list(max_terms=2))$control))

# Native diagnostics bind the same parameters and preserve discrete left limits.
d <- f$gratia_diagnostics(pars,c(0,.5,2),1,c(0,1,5),NULL)
near(d$mean,pars[,1])
near(d$variance,c(0,2,moments(3,3,f$control)[[1]]$variance))
stopifnot(identical(d$continuous,FALSE),is.function(d$rd),is.function(d$qf))
for(lower in c(TRUE,FALSE)) for(logp in c(TRUE,FALSE)) {
  near(d$cdf_left(c(0,.5,2),logp,lower),f$cdf(c(-1,0,1),pars,logp=logp,lower.tail=lower))
  near(d$cdf(c(0,1,5),logp,lower),f$cdf(c(0,1,5),pars,logp=logp,lower.tail=lower))
}
clone <- unserialize(serialize(f,NULL))
near(clone$gratia_diagnostics(pars,1,1,NULL,NULL)$variance,d$variance)
set.seed(23); draw <- f$rd(pars)
set.seed(23); near(f$rd(pars,wt=100,scale=3),draw)
set.seed(84)
for(s in c(.2,3)) {
  many <- cbind(rep(3,4000),s); y <- f$rd(many)
  v <- moments(3,s,f$control)[[1]]$variance
  stopifnot(abs(mean(y)-3)<6*sqrt(v/length(y)),abs(var(y)-v)<.15*v)
}
