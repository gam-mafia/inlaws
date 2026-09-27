library(inlaws)
near<-function(a,b,tol=1e-5)stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails<-function(expr)stopifnot(inherits(tryCatch(expr,error=identity),"error"))
loglik<-function(y,eta,wt=1) {
  mu<-exp(eta[,1]);s<-plogis(eta[,2])
  wt*(dlogis(log(y),log(mu*sinpi(s)/(pi*s)),s,log=TRUE)-log(y))
}
set.seed(29);n<-240L
dat<-data.frame(x=runif(n),z=runif(n),off1=runif(n,-.1,.1),off2=runif(n,-.1,.1))
pars<-cbind(exp(.4+sin(6*dat$x)+dat$off1),plogis(-1.6+.8*dat$z+dat$off2))
fam<-log_logisticls(); dat$y<-fam$rd(pars)
forms<-list(y~s(x,k=6)+offset(off1),~s(z,k=5)+offset(off2))
nei<-list(a=seq_len(n),ma=seq.int(4L,n,4L),d=as.vector(apply(matrix(seq_len(n),4L),2L,rev)),md=seq.int(4L,n,4L))
check<-function(f,method) {
  stopifnot(identical(f$method,method),identical(f$outer.info$conv,"full convergence"),
    all(is.finite(coef(f))),all(is.finite(f$gcv.ubre)),identical(dim(fitted(f)),c(n,2L)),
    all(fitted(f)[,1]>0),all(fitted(f)[,2]>0 & fitted(f)[,2]<1))
  p<-predict(f,type="response",se.fit=TRUE);l<-predict(f,se.fit=TRUE)
  near(p$fit,fitted(f));near(p$fit[,1],exp(l$fit[,1]));near(p$fit[,2],plogis(l$fit[,2]))
  near(p$se.fit[,1],p$fit[,1]*l$se.fit[,1]);near(p$se.fit[,2],p$fit[,2]*(1-p$fit[,2])*l$se.fit[,2])
  near(logLik(f),sum(loglik(dat$y,l$fit)),1e-7)
  stopifnot(all(is.finite(residuals(f))),is.finite(summary(f)$dev.expl))
}
reml<-gam(forms,data=dat,family=fam,method="REML");check(reml,"REML")
bfgs<-gam(forms,data=dat,family=log_logisticls(),method="REML",optimizer=c("outer","bfgs"));check(bfgs,"REML")
near(coef(reml),coef(bfgs),1e-3)
loo<-gam(forms,data=dat,family=log_logisticls(),method="NCV");check(loo,"NCV")
grouped<-gam(forms,data=dat,family=log_logisticls(),method="NCV",nei=nei);check(grouped,"NCV")
# NCV offset/subset handling and derivatives, independently recomputed.
G<-gam(forms,data=dat,family=log_logisticls(),fit=FALSE)
X<-predict(reml,type="lpmatrix"); beta<-coef(reml); yy<-dat$y
ll<-fam$ll(yy,X,beta,rep(1,n),fam,offset=G$offset,deriv=1,ncv=TRUE);ll$gamma<-1
cvfam<-fam;cvfam$qapprox<-FALSE
H<--ll$lbb+diag(1,length(beta))
cv<-cvfam$ncv(X,yy,rep(1,n),nei,beta,cvfam,ll,R=chol(H),offset=G$offset)
cveta<-attr(cv$NCV,"eta.cv")
near(cv$NCV,-sum(loglik(yy[nei$d],cveta)),1e-7)
llg<-ll;llg$gamma<-1.3
cvg<-cvfam$ncv(X,yy,rep(1,n),nei,beta,cvfam,llg,R=chol(H),offset=G$offset)
near(cvg$NCV,-1.3*sum(loglik(yy[nei$d],cveta))+.3*sum(loglik(yy[nei$d],predict(reml)[nei$d,])),1e-7)
direction<-matrix(seq(-.05,.05,length.out=length(beta)),ncol=1)
ll3<-fam$ll(yy,X,beta,rep(1,n),fam,offset=G$offset,deriv=3,d1b=direction,ncv=TRUE);ll3$gamma<-1
cv3<-cvfam$ncv(X,yy,rep(1,n),nei,beta,cvfam,ll3,R=chol(H),offset=G$offset,dH=ll3$d1H,db=direction,deriv=TRUE)
score<-function(step) {
  b<-beta+step*direction[,1]
  lr<-fam$ll(yy,X,b,rep(1,n),fam,offset=G$offset,deriv=1,ncv=TRUE);lr$gamma<-1
  cvfam$ncv(X,yy,rep(1,n),nei,b,cvfam,lr,R=chol(-lr$lbb+diag(1,length(beta))),offset=G$offset)$NCV
}
near(cv3$NCV1,(score(1e-4)-score(-1e-4))/2e-4,1e-5)
# Sandwich must be the cross-product of the observation scores.
sm<-X;jj<-attr(X,"lpi")
for(k in 1:2)sm[,jj[[k]]]<-X[,jj[[k]],drop=FALSE]*ll$l1[,k]
near(fam$sandwich(yy,X,beta,rep(1,n),fam,G$offset),crossprod(sm),1e-7)
# Replication weights and a direct unpenalized maximum likelihood reference.
dat$w<-rep(c(0,1,2),length.out=n)
a<-gam(list(y~x,~z),data=dat,weights=w,family=log_logisticls(),method="REML")
b<-gam(list(y~x,~z),data=dat[rep(seq_len(n),dat$w),],family=log_logisticls(),method="REML")
near(coef(a),coef(b));near(logLik(a),logLik(b))
opt<-optim(coef(a)+.01,function(b)-sum(loglik(yy,cbind(b[1]+b[2]*dat$x,b[3]+b[4]*dat$z),dat$w)),method="BFGS",control=list(reltol=1e-12))
near(coef(a),opt$par,1e-4);near(logLik(a),-opt$value,1e-6)
# Constant scale uses the same likelihood as the extended family.
cst<-gam(list(y~x,~1),data=dat,family=log_logisticls(),method="REML")
mu<-fitted(cst)
f<-log_logistic(mu[1,2])
near(f$aic(yy,mu[,1],wt=rep(1,n),dev=0),-2*as.numeric(logLik(cst)))
# Diagnostics, infinite variances, omission and prediction after serialization.
s<-mu[,2];v<-mu[,1]^2*(tan(pi*s)/(pi*s)-1)
near(residuals(cst,type="pearson"),(yy-mu[,1])/sqrt(v))
copy<-cst;copy$fitted.values[,2]<-.7
stopifnot(all(is.na(residuals(copy,type="pearson"))),all(is.finite(residuals(copy))))
dat$x[1]<-NA
missing<-gam(list(y~x,~z),data=dat,family=log_logisticls(),method="REML",na.action=na.exclude)
stopifnot(nobs(missing)==n-1L,length(residuals(missing))==n,is.na(residuals(missing)[1]))
copy<-unserialize(serialize(reml,NULL));nd<-dat[2:8,c("x","z","off1","off2")]
near(predict(copy,nd,type="response"),predict(reml,nd,type="response"))
for(bad_weight in list(rep(-1,n),rep(0,n)))fails(gam(forms,data=dat,weights=bad_weight,family=log_logisticls()))
fails(gam(list(y~1,~1,1+2~z-1),data=dat,family=log_logisticls()))
fails(bam(forms,data=dat,family=log_logisticls()))
