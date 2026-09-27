library(inlaws)
near<-function(a,b,tol=1e-8)stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails<-function(expr)stopifnot(inherits(tryCatch(expr,error=identity),"error"))
for (s in c(.001,.1,.3,.49,.5,.8,.99)) {
  mu<-c(.5,1,3)
  loc<-log(mu*sinpi(s)/(pi*s))
  f<-log_logistic(s); g<-log_logisticls(); pars<-cbind(mu,s)
  for (lt in c(TRUE,FALSE)) for(lp in c(TRUE,FALSE)) {
    p<-c(.01,.4,.99);if(lp)p<-log(p)
    q<-f$qf(p,mu,lower.tail=lt,log.p=lp)
    near(q,exp(qlogis(p,loc,s,lower.tail=lt,log.p=lp)))
    near(f$cdf(q,mu,lower.tail=lt,logp=lp),p,1e-7)
    near(g$qf(p,pars,lower.tail=lt,log.p=lp),q)
    near(g$cdf(q,pars,lower.tail=lt,logp=lp),p,1e-7)
  }
  near(f$cdf(c(-1,0,Inf),mu),c(0,0,1))
  near(f$qf(c(0,1),1),c(0,Inf))
  near(f$cdf(c(0,Inf),1,logp=TRUE),c(-Inf,0))
  # Deep log tails without cancellation from subtracting a CDF from one.
  for (lt in c(TRUE,FALSE)) {
    lp<-c(-50,-100,-500)
    q<-f$qf(lp,mu,lower.tail=lt,log.p=TRUE)
    near(f$cdf(q,mu,lower.tail=lt,logp=TRUE),lp,1e-7)
  }
  set.seed(17); a<-f$rd(mu)
  set.seed(17);near(a,exp(rlogis(length(mu),loc,s)))
  set.seed(17);near(a,g$rd(pars))
  near(f$qf(.7,mu,wt=100,scale=42),f$qf(.7,mu))
  near(g$cdf(mu,pars,wt=100,scale=42),g$cdf(mu,pars))
  if(s<.5)near(f$variance(mu),mu^2*(tan(pi*s)/(pi*s)-1)) else
    stopifnot(all(is.infinite(f$variance(mu))))
}
# Integrate in log-response coordinates, with independent logistic density.
for (s in c(.1,.3,.7)) {
  mu<-2;loc<-log(mu*sinpi(s)/(pi*s))
  near(integrate(function(z)dlogis(z,loc,s),-100,100,rel.tol=1e-10)$value,1)
  near(integrate(function(z)exp(z)*dlogis(z,loc,s),-100,100,rel.tol=1e-10)$value,mu)
  if(s<.5) {
    moment<-integrate(function(z)exp(2*z)*dlogis(z,loc,s),-100,100,rel.tol=1e-10)$value
    near(moment-mu^2,log_logistic(s)$variance(mu))
  }
}
# Tiny finite variance must not be rounded to zero by tan(x)/x - 1.
near(log_logistic(1e-9)$variance(1)/(pi*1e-9)^2,1/3)
for (f in list(log_logistic(),log_logistic(0),log_logistic(-.4))) {
  stopifnot(f$n.theta==1L)
  near(f$getTheta(TRUE),if(identical(f$ini.theta,qlogis(.4))) .4 else .25)
}
f<-log_logistic(); f$putTheta(qlogis(.6)); near(f$getTheta(TRUE),.6)
near(f$qf(.7,2),log_logistic(.6)$qf(.7,2))
fails(log_logisticls()$rd(c(1,.2)))
fails(log_logisticls()$rd(cbind(1,1)))
fails(f$rd(-1)); fails(f$putTheta(Inf))
