library("inlaws")
near <- function(a,b,tol=1e-6) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
jet <- getFromNamespace(".hp_lljet","inlaws")
moments <- getFromNamespace(".hp_moments","inlaws")
index <- getFromNamespace(".hp_index","inlaws")(4L)
f <- hyper_poissonls()
evaluate <- function(y,eta) {
  p <- exp(eta)
  jet(y,p[,1],p[,2],moments(p[,1],p[,2],f$control))
}
grid <- expand.grid(mu=c(.05,.5,2,20),sigma=c(.03,.5,1,3,30),y=c(0,1,7))
eta <- cbind(log(grid$mu),log(grid$sigma)); r <- evaluate(grid$y,eta)
# Difference all lower analytic partials, including every mixed fourth partial.
for(j in 1:2) {
  ep <- em <- eta; ep[,j]<-ep[,j]+1e-5; em[,j]<-em[,j]-1e-5
  rp <- evaluate(grid$y,ep); rm <- evaluate(grid$y,em)
  for(k in which(index[,j]>0)) {
    p <- index[k,]; p[j]<-p[j]-1L
    ix <- which(rowSums(abs(sweep(index,2,p)))==0)
    near(r[,k]*index[k,j],(rp[,ix]-rm[,ix])/2e-5,2e-4)
  }
}
# Independent likelihood, normalization and inversion for tensor finite differences.
reference <- function(y,eta) {
  mu<-exp(eta[1]); s<-exp(eta[2]); k<-0:1500
  prob <- function(a) { v<-a*k+lgamma(s)-lgamma(s+k); w<-exp(v-max(v)); w/sum(w) }
  a<-uniroot(function(a) sum(k*prob(a))-mu,c(-40,log((mu+1)*(s+1))+4),tol=1e-14)$root
  log(prob(a)[y+1L])
}
fd <- function(y,eta,p,h) {
  stencils <- lapply(p,function(d) {
    if(!d) return(list(x=0,w=1))
    x <- -4:4; target<-numeric(9); target[d+1L]<-factorial(d)
    list(x=x,w=solve(outer(0:8,x,function(k,x) x^k),target))
  })
  ix<-expand.grid(seq_along(stencils[[1]]$x),seq_along(stencils[[2]]$x))
  sum(vapply(seq_len(nrow(ix)),function(i) {
    a<-stencils[[1]]; b<-stencils[[2]]; j<-ix[i,1]; k<-ix[i,2]
    a$w[j]*b$w[k]*reference(y,eta+h*c(a$x[j],b$x[k]))/h^sum(p)
  },numeric(1)))
}
for(s in c(.2,1,5)) {
  e<-log(c(2,s)); a<-evaluate(3,matrix(e,1))
  for(k in 2:nrow(index)) {
    p<-index[k,]; ref<-fd(3,e,p,.08)
    near(a[,k]*prod(factorial(p)),ref,2e-3)
  }
}
# Coefficient derivatives and sandwich accumulate shared terms and weights.
set.seed(53); n<-15L; X<-cbind(1,runif(n),1,runif(n),runif(n))
attr(X,"lpi")<-list(c(1,2,5),c(3,4,5)); beta<-c(1,.2,.3,.2,.1)
y<-sample(0:8,n,TRUE); wt<-rep(c(0,.5,2),5)
off<-list(seq(-.1,.1,length.out=n),rep(.2,n))
l<-f$ll(y,X,beta,wt,f,off,deriv=1,ncv=TRUE)
score<-matrix(0,n,ncol(X)); jj<-attr(X,"lpi")
for(j in 1:2) score[,jj[[j]]]<-score[,jj[[j]],drop=FALSE]+X[,jj[[j]],drop=FALSE]*l$l1[,j]
near(l$lb,colSums(score)); near(f$sandwich(y,X,beta,wt,f,off),crossprod(score))
for(k in seq_along(beta)) {
  bp<-bm<-beta; bp[k]<-bp[k]+1e-5; bm[k]<-bm[k]-1e-5
  lp<-f$ll(y,X,bp,wt,f,off,deriv=1); lm<-f$ll(y,X,bm,wt,f,off,deriv=1)
  near(l$lb[k],(lp$l-lm$l)/2e-5); near(l$lbb[,k],(lp$lb-lm$lb)/2e-5)
}
et<-matrix(0,n,2); et[1,]<-1000; y[1]<-1e300
z<-f$ll(y,X,beta,wt,f,eta=et,deriv=4,ncv=TRUE)
stopifnot(z$l0[1]==0,all(z$l1[1,]==0),all(z$l4[1,]==0))

# Large observation-to-mode gaps exercise the log-Gamma/polygamma path.
# Tight inversion avoids amplifying root tolerance when differencing tiny
# dispersion derivatives of very unlikely observations.
large<-expand.grid(mu=c(1001,2000),sigma=c(.1,1,30,1e4),y=c(0,7,3000))
ee<-cbind(log(large$mu),log(large$sigma))
ctl<-hyper_poissonls(control=list(mean_tol=1e-13,sum_tol=1e-13))$control
large_jet<-function(e) {
  p<-exp(e)
  jet(large$y,p[,1],p[,2],moments(p[,1],p[,2],ctl))
}
a<-large_jet(ee)
near(a[large$sigma==1,1],dpois(large$y[large$sigma==1],large$mu[large$sigma==1],log=TRUE),1e-10)
for(j in 1:2) {
  ep<-em<-ee;ep[,j]<-ep[,j]+1e-4;em[,j]<-em[,j]-1e-4
  ap<-large_jet(ep);am<-large_jet(em)
  for(k in which(index[,j]>0)) {
    p<-index[k,];p[j]<-p[j]-1L
    ix<-which(rowSums(abs(sweep(index,2,p)))==0)
    expected<-(ap[,ix]-am[,ix])/2e-4
    stopifnot(max(abs(a[,k]*index[k,j]-expected)/pmax(1,abs(expected)))<4e-4)
  }
}
