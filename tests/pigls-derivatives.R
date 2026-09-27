library(inlaws)
near <- function(x,y,tol=2e-5) stopifnot(length(x)==length(y),all(is.finite(x)),all(is.finite(y)),
  max(abs(x-y)/pmax(1,abs(x),abs(y))) < tol)
differentiate <- function(f,x,j,h=1e-4) {
  e<-numeric(length(x));e[j]<-h
  (f(x-2*e)-8*f(x-e)+8*f(x+e)-f(x+2*e))/(12*h)
}
jet<-getFromNamespace(".pig_lljet","inlaws")
pmf<-getFromNamespace(".pig_logpmf","inlaws")
idx<-getFromNamespace(".pig_index","inlaws")
column<-function(p,q) (p+q)*(p+q+1)/2+q+1
# All mixed predictor derivatives through order four, with row-specific theta.
grid<-expand.grid(y=c(0,1,7,40),mu=c(.03,2,30),theta=c(1e-8,.05,1,10))
z<-jet(grid$y,grid$mu,log(grid$theta),predictor=TRUE)
near(z[,1],pmf(grid$y,grid$mu,log(grid$theta)),1e-10)
for(k in 2:15) {
  p<-idx[k,1];q<-idx[k,2];j<-if(p>0) 1 else 2
  pp<-p-(j==1);qq<-q-(j==2)
  fn<-function(x) jet(grid$y,grid$mu*exp(x[1]),log(grid$theta)+x[2],predictor=TRUE)[,column(pp,qq)]*
    factorial(pp)*factorial(qq)
  near(differentiate(fn,c(0,0),j),z[,k]*factorial(p)*factorial(q))
}
# Coefficient-level contractions, offsets, fractional weights and zero weights.
set.seed(42); n<-14L
X<-cbind(1,seq(-1,1,length.out=n),rnorm(n),1,runif(n))
attr(X,"lpi")<-list(1:3,4:5)
beta<-c(.7,.1,-.1,-.2,.3); y<-rep(c(0,1,4,8),length.out=n)
w<-c(0,.4,rep(1,6),rep(2,6));off<-list(runif(n,-.2,.2),runif(n,-.3,.3))
fam<-pigls()
ll<-function(b,deriv=1,...) fam$ll(y,X,b,w,fam,offset=off,deriv=deriv,...)
base<-ll(beta)
for(j in seq_along(beta)) {
  near(base$lb[j],differentiate(function(b) ll(b,0)$l,beta,j))
  near(base$lbb[,j],differentiate(function(b) ll(b)$lb,beta,j))
}
directions<-matrix(rnorm(10),5,2);second<-matrix(rnorm(15),5,3)
V<-crossprod(matrix(rnorm(25),5))+diag(5);R<-chol(solve(V),pivot=TRUE)
full<-ll(beta,4,d1b=directions,d2b=second,fh=R,D=rep(1,5))
contracted<-ll(beta,2,d1b=directions,fh=V)
for(j in 1:2) {
  fd<-differentiate(function(a) ll(beta+a*directions[,j])$lbb,0,1)
  near(full$d1H[[j]],fd);near(contracted$d1H[j],sum(V*fd))
}
kk<-0L
for(j in 1:2)for(k in j:2) {
  kk<-kk+1L
  fd<-differentiate(function(a) ll(beta+a*directions[,k],3,d1b=directions)$d1H[[j]],0,1)+
    differentiate(function(a) ll(beta+a*second[,kk])$lbb,0,1)
  near(full$trHid2H[kk],sum(V*fd))
}
scores<-t(vapply(seq_len(n),function(i) {
  wi<-numeric(n);wi[i]<-w[i];fam$ll(y,X,beta,wi,fam,offset=off,deriv=1)$lb
},numeric(5)))
near(fam$sandwich(y,X,beta,w,fam,off),crossprod(scores))
eta<-cbind(drop(X[,1:3]%*%beta[1:3])+off[[1]],drop(X[,4:5]%*%beta[4:5])+off[[2]])
near(base$l,fam$ll(y,X,beta,w,fam,eta=eta)$l)
eta[1,]<-c(Inf,-Inf)
near(base$l,fam$ll(y,X,beta,w,fam,eta=eta)$l)
Xs<-X[-1,,drop=FALSE];attr(Xs,"lpi")<-attr(X,"lpi")
omit<-fam$ll(y[-1],Xs,beta,w[-1],fam,offset=lapply(off,`[`,-1),deriv=1)
near(base$l,omit$l);near(base$lb,omit$lb);near(base$lbb,omit$lbb)
eta[2,]<-c(1000,-1000)
stopifnot(identical(fam$ll(y,X,beta,w,fam,eta=eta)$l,-Inf))
# Near-Poisson dispersion derivatives retain the small nonzero leading term.
for(th in c(1e-6,1e-10,1e-14)) {
  z<-jet(0,2,log(th),predictor=TRUE)
  near(z[,3]/th,2);near(2*z[,6]/th,2)
}
