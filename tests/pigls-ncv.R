library(inlaws)
near<-function(x,y,tol=2e-5) stopifnot(length(x)==length(y),all(is.finite(x)),all(is.finite(y)),
  max(abs(x-y)/pmax(1,abs(x),abs(y))) < tol)
pmf<-getFromNamespace(".pig_logpmf","inlaws")
set.seed(762);n<-450L
d<-data.frame(x=runif(n),z=runif(n),a=runif(n,-.2,.2),b=runif(n,-.2,.2))
d$y<-pigls()$rd(cbind(exp(.7+.6*sin(6*d$x)+d$a),exp(-.3+1.5*(d$z-.5)+d$b)))
form<-list(y~s(x,k=5)+offset(a),~s(z,k=4)+offset(b))
loo<-list(a=seq_len(n),ma=seq_len(n),d=seq_len(n),md=seq_len(n))
grouped<-list(a=seq_len(n),ma=seq.int(3,n,3),
  d=as.vector(apply(matrix(seq_len(n),3),2,rev)),md=seq.int(3,n,3))
check<-function(f)stopifnot(identical(f$method,"NCV"),
  identical(f$outer.info$conv,"full convergence"),all(is.finite(coef(f))),
  all(is.finite(f$Vp)),all(is.finite(f$sp)),all(f$sp>0),all(is.finite(f$gcv.ubre)))
a<-gam(form,data=d,family=pigls(),method="NCV");check(a)
b<-gam(form,data=d,family=pigls(),method="NCV",nei=loo);check(b)
near(coef(a),coef(b));near(a$sp,b$sp);near(a$gcv.ubre,b$gcv.ubre)
c<-gam(form,data=d,family=pigls(),method="NCV",nei=grouped,gamma=1.2);check(c)
# Actual smoothing-parameter derivatives, with repeated prediction rows,
# different deletion/prediction sets, and offsets in both predictors.
nei<-list(a=c(4,5,6,10,11,12,4,5,6,19,20,21),ma=c(3,6,9,12),
          d=c(5,11,12,5,20),md=c(1,3,4,5),jackknife=FALSE)
ns<-asNamespace("mgcv");setup<-get("Sl.setup",ns)
repara<-get("Sl.initial.repara",ns);fit5<-get("gam.fit5",ns)
fam<-pigls();fam$qapprox<-FALSE
G<-gam(form,data=d,family=fam,fit=FALSE)
Sl<-setup(G);X<-repara(Sl,G$X,both.sides=FALSE)
rho<-rep(log(2),length(G$sp));start<-NULL
call<-function(rho,method="REML",deriv=1) fit5(X,G$y,rho,Sl,
  weights=G$w,offset=G$offset,family=fam,deriv=deriv,scoreType=method,Mp=-1,
  nei=nei,gamma=1.2,start=start,control=gam.control(epsilon=1e-10))
base<-call(rho,deriv=2);start<-base$coefficients
for(j in seq_along(rho)) {
  step<-rep(0,length(rho));step[j]<-.005
  fmm<-call(rho-2*step);fm<-call(rho-step);fp<-call(rho+step);fpp<-call(rho+2*step)
  near(base$REML1[j],(fmm$REML-8*fm$REML+8*fp$REML-fpp$REML)/.06,1e-4)
  near(base$REML2[,j],(fmm$REML1-8*fm$REML1+8*fp$REML1-fpp$REML1)/.06,3e-4)
}
f<-call(rho,"NCV")
for(j in seq_along(rho)) {
  rp<-rm<-rho;rp[j]<-rp[j]+1e-4;rm[j]<-rm[j]-1e-4
  near(f$NCV1[j],(call(rp,"NCV")$NCV-call(rm,"NCV")$NCV)/2e-4,3e-4)
}
# Reconstruct the reported NCV score from independent probability evaluations.
ecv<-attr(f$NCV,"eta.cv");eta<-f$linear.predictors;ix<-nei$d
near(as.numeric(f$NCV),-1.2*sum(pmf(d$y[ix],exp(ecv[,1]),ecv[,2]))+
  .2*sum(pmf(d$y[ix],exp(eta[ix,1]),eta[ix,2])))
# Exact refits with the same basis and penalties check the CV predictors.
exact<-ecv*0;prevA<-prevD<-0
for(k in seq_along(nei$ma)) {
  drop<-nei$a[(prevA+1):nei$ma[k]];pred<-nei$d[(prevD+1):nei$md[k]]
  w<-G$w;w[drop]<-0
  ff<-fit5(X,G$y,rho,Sl,weights=w,offset=G$offset,family=fam,deriv=0,
    Mp=-1,gamma=1.2,start=start,control=gam.control(epsilon=1e-10))
  exact[(prevD+1):nei$md[k],]<-ff$linear.predictors[pred,]
  prevA<-nei$ma[k];prevD<-nei$md[k]
}
near(ecv,exact,.015)
