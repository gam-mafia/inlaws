library(inlaws)
near <- function(a,b,tol=1e-5) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch(expr,error=identity),"error"))
loglik <- function(y,mu,s) dlogis(log(y),log(mu*sinpi(s)/(pi*s)),s,log=TRUE)-log(y)
set.seed(27)
n <- 200L
d <- data.frame(x=runif(n),z=runif(n),o=runif(n,-.1,.1),os=runif(n,-.1,.1))
d$y <- log_logistic(.25)$rd(exp(.4+sin(6*d$x)+d$o))
form <- y~s(x,k=6)+offset(o)
check <- function(f,method,fixed=FALSE) {
  message("Checking ",deparse(f$call)[1],"; ",method,"; fixed=",fixed)
  stopifnot(identical(f$method,method),isTRUE(f$converged),all(is.finite(coef(f))),
            all(is.finite(f$gcv.ubre)),f$family$n.theta==as.integer(!fixed))
  sc<-f$family$getTheta(TRUE)
  stopifnot(sc>0,sc<1)
  near(as.numeric(logLik(f)),sum(f$prior.weights*loglik(f$y,fitted(f),sc)),1e-6)
  near(f$deviance,sum(f$family$dev.resids(f$y,fitted(f),f$prior.weights)),1e-6)
  pr<-predict(f,type="response",se.fit=TRUE); lp<-predict(f,se.fit=TRUE)
  near(pr$fit,exp(lp$fit)); near(pr$se.fit,pr$fit*lp$se.fit)
  stopifnot(all(is.finite(residuals(f))),all(is.finite(residuals(f,type="response"))))
}
# Match nb(): callers need no external scale argument. Keep scale=1 only
# for absolute logLik reference checks: mgcv 1.9-4 mislabels default bam's
# fixed scale in its outer reporting code, without changing the fit or AIC.
check_default <- function(default, explicit) {
  stopifnot(isTRUE(default$converged),identical(default$method,explicit$method))
  near(coef(default),coef(explicit),1e-8)
  near(default$sp,explicit$sp,1e-8)
  near(default$family$getTheta(),explicit$family$getTheta(),1e-8)
  near(default$Vp,explicit$Vp,1e-8)
  near(default$gcv.ubre,explicit$gcv.ubre,1e-8)
  near(default$deviance,explicit$deviance,1e-8)
  near(AIC(default),AIC(explicit),1e-8)
  near(predict(default,type="response"),fitted(explicit),1e-8)
  near(default$sig2,1)
}
for (fixed in c(FALSE,TRUE)) {
  theta <- if(fixed) .25 else NULL
  for (method in c("REML","ML","NCV")) {
    f<-gam(form,data=d,family=log_logistic(theta),method=method)
    check(f,method,fixed)
    if(method=="NCV" && !fixed) loo<-f
  }
  for (method in c("REML","ML","fREML")) {
    f<-bam(form,data=d,family=log_logistic(theta),method=method,scale=1,nthreads=1)
    check(f,method,fixed)
    default<-bam(form,data=d,family=log_logistic(theta),method=method,nthreads=1)
    check_default(default,f)
    if(method=="fREML") ordinary<-f
  }
  f<-bam(form,data=d,family=log_logistic(theta),method="fREML",discrete=TRUE,scale=1,nthreads=1)
  default<-bam(form,data=d,family=log_logistic(theta),method="fREML",discrete=TRUE,nthreads=1)
  check_default(default,f)
  check(f,"fREML",fixed); near(fitted(f),fitted(ordinary),.005)
  near(f$family$getTheta(TRUE),ordinary$family$getTheta(TRUE),.005)
  if(packageVersion("mgcv")>="1.9.4") {
    f<-bam(form,data=d,family=log_logistic(theta),method="NCV",discrete=TRUE,scale=1,nthreads=1)
    check(f,"NCV",fixed)
    default<-bam(form,data=d,family=log_logistic(theta),method="NCV",discrete=TRUE,nthreads=1)
    check_default(default,f)
  }
}
nei<-list(a=seq_len(n),ma=seq_len(n),d=seq_len(n),md=seq_len(n))
f<-gam(form,data=d,family=log_logistic(),method="NCV",nei=nei)
near(coef(f),coef(loo));near(f$sp,loo$sp)
grouped<-list(a=seq_len(n),ma=seq.int(4L,n,4L),d=as.vector(apply(matrix(seq_len(n),4L),2L,rev)),md=seq.int(4L,n,4L))
for (theta in list(NULL,.25)) {
  f<-gam(form,data=d,family=log_logistic(theta),method="NCV",nei=grouped)
  check(f,"NCV",!is.null(theta))
  if(packageVersion("mgcv")>="1.9.4") {
    f<-bam(form,data=d,family=log_logistic(theta),method="NCV",discrete=TRUE,nei=grouped,scale=1,nthreads=1)
    check(f,"NCV",!is.null(theta))
    default<-bam(form,data=d,family=log_logistic(theta),method="NCV",discrete=TRUE,nei=grouped,nthreads=1)
    check_default(default,f)
  }
}
# Reusing a constructor must not mutate previous models or the template.
for (engine in c("gam","bam","discrete")) {
  template<-log_logistic()
  fit <- function(dat) if(engine=="gam") gam(form,data=dat,family=template,method="REML") else
    bam(form,data=dat,family=template,method="fREML",discrete=engine=="discrete",nthreads=1)
  a<-fit(d); old<-a$family$getTheta(); oldll<-logLik(a)
  dd<-d; dd$y<-d$y^1.5; b<-fit(dd)
  near(a$family$getTheta(),old);near(logLik(a),oldll)
  near(template$getTheta(TRUE),.25)
}
# Integer replication, zero-row omission and non-integer likelihood weights.
d$w<-rep(c(0,1,2,3),length.out=n)
a<-gam(y~x+offset(o),data=d,weights=w,family=log_logistic(),method="REML")
b<-gam(y~x+offset(o),data=d[rep(seq_len(n),d$w),],family=log_logistic(),method="REML")
c<-gam(y~x+offset(o),data=d[d$w>0,],weights=w,family=log_logistic(),method="REML")
near(coef(a),coef(b),1e-4);near(coef(a),coef(c),1e-4)
d$w<-rep(c(0,.5,1,1.5),length.out=n)
f<-gam(form,data=d,weights=w,family=log_logistic(),method="REML");check(f,"REML")
# Unpenalized ML with fixed scale has the same coefficients as direct MLE.
f<-gam(y~x+offset(o),data=d,family=log_logistic(.25),method="ML")
opt<-optim(coef(f),function(b)-sum(loglik(d$y,exp(b[1]+b[2]*d$x+d$o),.25)),method="BFGS")
near(coef(f),opt$par,1e-5)
# Infinite response variance still supports likelihood fitting.
f<-gam(form,data=d,family=log_logistic(.7),method="REML")
stopifnot(all(is.na(residuals(f,type="pearson"))),all(is.finite(residuals(f))))
for (x in list(NA_real_,Inf,c(.1,.2),"a",-1,1)) fails(log_logistic(x))
fails(log_logistic(link="identity"));fails(log_logisticls(list("log","log")))
for(y in list(c(0,1),c(-1,1),c(Inf,1),matrix(1,2,2)))
  fails(gam(y~1,family=log_logistic()))
fails(gam(form,data=d,family=log_logistic(),method="REML",scale=2))
