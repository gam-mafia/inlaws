library(inlaws)
near<-function(x,y,tol=1e-8) stopifnot(length(x)==length(y),all(is.finite(x)),all(is.finite(y)),
  max(abs(x-y)/pmax(1,abs(x),abs(y)))<tol)
error<-function(expr) stopifnot(inherits(tryCatch(expr,error=identity),"error"))
f<-pigls();pars<-cbind(c(.1,1,3,10),c(.02,.5,1,3));y<-c(0,1,4,15)
pmf<-getFromNamespace(".pig_logpmf","inlaws");sat<-getFromNamespace(".pig_saturated","inlaws")
near(pmf(y,pars[,1],log(pars[,2])),vapply(1:4,function(i)pmf(y[i],pars[i,1],log(pars[i,2])),numeric(1)))
near(sat(y,log(pars[,2])),do.call(rbind,lapply(1:4,function(i)sat(y[i],log(pars[i,2])))))
for(logp in c(FALSE,TRUE)) near(f$cdf(y,pars,logp=logp),vapply(1:4,function(i)
  pig(pars[i,2])$cdf(y[i],pars[i,1],logp=logp),numeric(1)))
for(p in list(.7,c(.1,.2,.8,.95))) {
  q<-f$qf(p,pars)
  near(q,vapply(1:4,function(i)pig(pars[i,2])$qf(rep_len(p,4)[i],pars[i,1]),numeric(1)))
  stopifnot(all(f$cdf(q,pars)>=p),all(f$cdf(q-1,pars)<p))
}
near(f$cdf(3,pars),vapply(1:4,function(i)pig(pars[i,2])$cdf(3,pars[i,1]),numeric(1)))
near(f$qf(c(.1,.5,.9),pars[2,,drop=FALSE]),pig(pars[2,2])$qf(c(.1,.5,.9),pars[2,1]))
stopifnot(f$qf(1,pars[1,,drop=FALSE])==Inf,f$qf(0,pars[1,,drop=FALSE])==0,
  f$qf(1,cbind(0,1))==0,f$cdf(0,cbind(0,1))==1,
  is.na(f$cdf(1,cbind(2,NA_real_))),is.na(f$rd(cbind(NA_real_,1))))
set.seed(7);a<-f$rd(pars)
set.seed(7);b<-getFromNamespace(".pig_random","inlaws")(pars[,1],log(pars[,2]))
near(a,b)
# Distinct row-specific dispersions reproduce their analytic moments.
set.seed(137);large<-pars[rep(1:4,each=60000),];draws<-f$rd(large)
for(i in 1:4) {
  z<-draws[seq_len(60000)+(i-1)*60000]
  near(mean(z),pars[i,1],.03)
  near(var(z),pars[i,1]+pars[i,2]*pars[i,1]^2,.065)
}
error(f$cdf(1,c(2,1)));error(f$qf(.5,cbind(1,-1)))
error(f$rd(cbind(-1,1)));error(f$cdf(1:3,pars))
stopifnot(length(f$rd(matrix(numeric(),0,2)))==0,
          length(f$cdf(1,matrix(numeric(),0,2)))==0)
