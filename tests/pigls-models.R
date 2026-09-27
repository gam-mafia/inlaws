library(inlaws)
near<-function(x,y,tol=1e-5) stopifnot(length(x)==length(y),all(is.finite(x)),all(is.finite(y)),
  max(abs(x-y)/pmax(1,abs(x),abs(y))) < tol)
error<-function(expr,pattern) {
  e<-tryCatch(expr,error=identity);stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e)))
}
pmf<-getFromNamespace(".pig_logpmf","inlaws")
sat<-getFromNamespace(".pig_saturated","inlaws")
set.seed(721);n<-800L
d<-data.frame(x=runif(n),z=runif(n),a=runif(n,-.15,.15),b=runif(n,-.15,.15))
mu<-exp(.8+.7*sin(2*pi*d$x)+d$a);theta<-exp(-.5+1.2*(d$z-.5)+d$b)
d$y<-pigls()$rd(cbind(mu,theta))
form<-list(y~s(x,k=6)+offset(a),~s(z,k=4)+offset(b))
check<-function(f,method) {
  stopifnot(identical(f$method,method),identical(f$outer.info$conv,"full convergence"),
    all(is.finite(coef(f))),all(is.finite(f$Vp)),all(is.finite(f$sp)),all(f$sp>0),
    identical(colnames(fitted(f)),c("mu","theta")),all(fitted(f)>0))
  pr<-predict(f,type="response",se.fit=TRUE); lp<-predict(f,type="link",se.fit=TRUE)
  near(pr$fit,fitted(f));near(pr$fit,exp(lp$fit));near(pr$se.fit,pr$fit*lp$se.fit)
  near(as.numeric(logLik(f)),sum(pmf(d$y,pr$fit[,1],lp$fit[,2])))
  near(AIC(f),-2*as.numeric(logLik(f))+2*attr(logLik(f),"df"))
  near(f$deviance,sum(residuals(f,"deviance")^2))
}
f<-gam(form,data=d,family=pigls(),method="REML");check(f,"REML")
g<-gam(form,data=d,family=pigls(),method="REML",optimizer=c("outer","bfgs"));check(g,"REML")
near(f$gcv.ubre,g$gcv.ubre,1e-5);near(log(fitted(f)),log(fitted(g)),.002)
stopifnot(cor(log(fitted(f)[,1]),log(mu))>.9,cor(log(fitted(f)[,2]),log(theta))>.75)
# Covariance-based prediction checks on new data, with both offsets.
nd<-d[1:12,c("x","z","a","b")];X<-predict(f,nd,type="lpmatrix")
pr<-predict(f,nd,type="response",se.fit=TRUE);lp<-predict(f,nd,type="link",se.fit=TRUE)
for(j in 1:2) {
  jj<-attr(X,"lpi")[[j]];xx<-X[,jj,drop=FALSE]
  se<-sqrt(rowSums((xx%*%f$Vp[jj,jj,drop=FALSE])*xx))
  near(lp$se.fit[,j],se);near(pr$se.fit[,j],se*pr$fit[,j])
}
near(residuals(f,"response"),d$y-fitted(f)[,1])
near(residuals(f,"pearson"),(d$y-fitted(f)[,1])/sqrt(fitted(f)[,1]+fitted(f)[,2]*fitted(f)[,1]^2))
# Null deviance retains each fitted theta and mean offset.
t<-log(fitted(f)[,2]);satll<-sum(sat(d$y,t,derivatives=FALSE)[,1])
opt<-optimize(function(a)-sum(pmf(d$y,exp(a+d$a),t)),c(-5,5),tol=1e-9)
near(f$null.deviance,2*(satll+opt$objective))
# Supplied starts and model serialization retain parameters and predictions.
again<-gam(form,data=d,family=pigls(),method="REML",start=coef(f),sp=f$sp)
near(fitted(again),fitted(f),1e-4)
copy<-unserialize(serialize(f,NULL));near(predict(copy,nd,type="response"),pr$fit)
# Fully parametric likelihood weighting and joint likelihood optimization.
d$w<-rep(0:3,length.out=n)
pform<-list(y~x+offset(a),~z+offset(b))
a<-gam(pform,data=d,weights=w,family=pigls(),method="REML")
b<-gam(pform,data=d[rep(seq_len(n),d$w),],family=pigls(),method="REML")
c<-gam(pform,data=d[d$w>0,],weights=w,family=pigls(),method="REML")
near(coef(a),coef(b));near(coef(a),coef(c))
near(as.numeric(logLik(a)),as.numeric(logLik(b)))
objective<-function(z)-sum(d$w*pmf(d$y,exp(z[1]+z[2]*d$x+d$a),z[3]+z[4]*d$z+d$b))
mle<-optim(coef(a)+.01,objective,method="BFGS",control=list(reltol=1e-12))
near(coef(a),mle$par,2e-4);near(as.numeric(logLik(a)),-mle$value)
# Constant dispersion remains an estimated coefficient, not a fixed argument.
constant<-gam(list(y~x,~1),data=d,family=pigls(),method="REML")
stopifnot(diff(range(fitted(constant)[,2]))==0)
# Without a mean intercept, null deviance retains only the mean offset.
no_intercept<-gam(list(y~0+x+offset(a),~z+offset(b)),data=d,family=pigls(),method="REML")
near(no_intercept$null.deviance,sum(no_intercept$family$dev.resids(d$y,
  cbind(exp(d$a),fitted(no_intercept)[,2]),rep(1,n))))
# Repeated covariates with distinct coefficient blocks are supported.
same_covariate<-gam(list(y~x,~x),data=d,family=pigls(),method="REML")
stopifnot(all(is.finite(coef(same_covariate))))
Xs<-predict(same_covariate,type="lpmatrix")
near(same_covariate$family$ll(d$y,Xs,coef(same_covariate),rep(1,n),
  same_covariate$family,deriv=1)$lb,rep(0,4),1e-4)
# Missing observations and subsetting use the ordinary gam mechanisms.
d$x[1]<-NA
missing<-gam(pform,data=d,family=pigls(),method="REML",na.action=na.exclude)
stopifnot(nobs(missing)==n-1L,length(residuals(missing))==n,is.na(residuals(missing)[1]))
subset<-gam(pform,data=d,family=pigls(),method="REML",subset=seq_len(n)>3)
stopifnot(nobs(subset)==n-3L)
for(link in list("log",list("log","identity"),list(NA_character_,"log"))) error(pigls(link),"two.*log")
for(y in list(c(-1,1),c(.5,1),c(Inf,1),c(0,0)))
  error(gam(list(y~1,~1),data=data.frame(y=y),family=pigls()),"integer|all-zero")
for(w in list(c(-1,1),c(Inf,1),c(0,0)))
  error(gam(list(y~1,~1),data=data.frame(y=c(1,2)),weights=w,family=pigls()),"weights")
error(gam(y~z,data=d,family=pigls()),"two formulae")
error(gam(list(y~1,~1,1+2~z-1),data=d,family=pigls()),"shared")
error(bam(form,data=d,family=pigls()),"general families|requires gam")
error(gam(form,data=d,family=pigls(),method="QNCV"),"QNCV")
