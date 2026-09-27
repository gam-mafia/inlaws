library(inlaws)
near <- function(x,y,tol=1e-5) {
  stopifnot(length(x)==length(y),all(is.finite(x)),all(is.finite(y)),
            max(abs(x-y)) <= tol*max(1,abs(x),abs(y)))
}
error <- function(expr, pattern) {
  e <- tryCatch(expr,error=identity)
  stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e)))
}
pmf <- getFromNamespace(".pig_logpmf","inlaws")
set.seed(912)
n <- 180L; d <- data.frame(x=runif(n), o=runif(n,-.2,.2))
d$mu <- exp(.3 + sin(2*pi*d$x) + d$o)
d$y <- pig(.6)$rd(d$mu)
form <- y ~ s(x,k=6) + offset(o)
check <- function(f, method, fixed=FALSE) {
  message("Checking ", deparse(f$call)[1], "; ", method, "; fixed=", fixed)
  stopifnot(identical(f$method,method),isTRUE(f$converged),
    all(is.finite(coef(f))),all(is.finite(f$Vp)),
    all(is.finite(f$sp)),all(f$sp>0),all(is.finite(f$gcv.ubre)),
    is.finite(f$family$getTheta(TRUE)), f$family$getTheta(TRUE)>0,
    f$family$n.theta == if (fixed) 0L else 1L)
  if (fixed) near(f$family$getTheta(TRUE),.6)
  pr <- predict(f,type="response",se.fit=TRUE)
  near(pr$fit,fitted(f)); near(pr$fit,exp(predict(f,type="link")))
  stopifnot(all(is.finite(pr$se.fit)),all(pr$se.fit>=0))
  near(f$deviance,sum(residuals(f,"deviance")^2))
  near(as.numeric(logLik(f)),sum(pmf(d$y,fitted(f),f$family$getTheta())),1e-6)
  near(AIC(f),-2*as.numeric(logLik(f))+2*attr(logLik(f),"df"))
}
# Default bam calls must fit the same model as explicit scale=1, just as
# nb() does. Keep the latter for absolute logLik checks above: mgcv 1.9-4
# mislabels default bam's fixed scale, affecting logLik but not the fit/AIC.
check_default <- function(default, explicit) {
  stopifnot(isTRUE(default$converged), identical(default$method, explicit$method))
  near(coef(default), coef(explicit), 1e-8)
  near(default$sp, explicit$sp, 1e-8)
  near(default$family$getTheta(TRUE), explicit$family$getTheta(TRUE), 1e-8)
  near(default$Vp, explicit$Vp, 1e-8)
  near(default$gcv.ubre, explicit$gcv.ubre, 1e-8)
  near(default$deviance, explicit$deviance, 1e-8)
  near(AIC(default), AIC(explicit), 1e-8)
  near(predict(default, type="response"), fitted(explicit), 1e-8)
  near(default$sig2, 1)
}
for (fixed in c(FALSE,TRUE)) {
  theta <- if (fixed) .6 else NULL
  for (method in c("REML","ML","NCV")) {
    f <- gam(form,data=d,family=pig(theta),method=method)
    check(f,method,fixed)
    if (method=="NCV" && !fixed) loo <- f
    if (method=="REML" && !fixed) reml <- f
  }
  for (method in c("REML","ML","fREML")) {
    f <- bam(form,data=d,family=pig(theta),method=method,scale=1,nthreads=1)
    check(f,method,fixed)
    default <- bam(form,data=d,family=pig(theta),method=method,nthreads=1)
    check_default(default,f)
    if (method=="fREML") ordinary <- f
  }
  discrete <- bam(form,data=d,family=pig(theta),method="fREML",discrete=TRUE,scale=1,nthreads=1)
  check(discrete,"fREML",fixed)
  default <- bam(form,data=d,family=pig(theta),method="fREML",discrete=TRUE,nthreads=1)
  check_default(default,discrete)
  near(fitted(discrete),fitted(ordinary),.005)
  near(discrete$family$getTheta(TRUE),ordinary$family$getTheta(TRUE),.005)
}
# Explicit singleton neighbourhoods match default leave-one-out NCV.
nei <- list(a=seq_len(n),ma=seq_len(n),d=seq_len(n),md=seq_len(n))
f <- gam(form,data=d,family=pig(),method="NCV",nei=nei)
check(f,"NCV"); near(coef(f),coef(loo)); near(f$sp,loo$sp)
near(f$gcv.ubre,loo$gcv.ubre); near(f$family$getTheta(),loo$family$getTheta())
if (packageVersion("mgcv") >= "1.9.4") {
  implicit <- bam(form,data=d,family=pig(),method="NCV",discrete=TRUE,scale=1,nthreads=1)
  explicit <- bam(form,data=d,family=pig(),method="NCV",nei=nei,discrete=TRUE,scale=1,nthreads=1)
  check(implicit,"NCV"); check(explicit,"NCV")
  near(coef(implicit),coef(explicit)); near(implicit$sp,explicit$sp)
  near(implicit$gcv.ubre,explicit$gcv.ubre)
}
nei$ma <- nei$md <- seq.int(3,n,3)
for (fixed in c(FALSE,TRUE)) {
  theta <- if (fixed) .6 else NULL
  f <- gam(form,data=d,family=pig(theta),method="NCV",nei=nei)
  check(f,"NCV",fixed)
  if (packageVersion("mgcv") >= "1.9.4") {
    f <- bam(form,data=d,family=pig(theta),method="NCV",nei=nei,
             discrete=TRUE,scale=1,nthreads=1)
    check(f,"NCV",fixed)
    default <- bam(form,data=d,family=pig(theta),method="NCV",nei=nei,
                   discrete=TRUE,nthreads=1)
    check_default(default,f)
  }
}
# Reusing a template does not change previously fitted dispersions.
template <- pig(-.4)
a <- gam(form,data=d,family=template,method="REML")
theta <- a$family$getTheta(); predictions <- fitted(a)
b <- gam(y~1,data=d,family=template,method="ML")
near(a$family$getTheta(),theta); near(fitted(a),predictions)
near(template$getTheta(TRUE),.4)
# Both bam paths must also isolate templates and previous fitted objects.
for (discrete in c(FALSE, TRUE)) {
  template <- pig(-.4)
  first <- bam(form,data=d,family=template,method="fREML",discrete=discrete,scale=1,nthreads=1)
  th <- first$family$getTheta()
  variance <- first$family$variance(fitted(first))
  second <- bam(y~s(x,k=4),data=d,family=template,method="fREML",discrete=discrete,scale=1,nthreads=1)
  near(first$family$getTheta(),th)
  near(first$family$variance(fitted(first)),variance)
  near(template$getTheta(TRUE),.4)
  near(as.numeric(logLik(first)),sum(pmf(d$y,fitted(first),th)))
}
# Intercept-only fit agrees with independent direct likelihood optimization.
mle <- optim(c(coef(b),b$family$getTheta())+.05,function(z)
  -sum(pmf(d$y,exp(z[1]),z[2])),method="BFGS",control=list(reltol=1e-11))
near(c(coef(b),b$family$getTheta()),mle$par,2e-4)
# Likelihood weighting equals replication, including offsets and zero weights.
d$w <- rep(0:3,length.out=n)
a <- gam(y~x+offset(o),data=d,weights=w,family=pig(),method="REML")
b <- gam(y~x+offset(o),data=d[rep(seq_len(n),d$w),],family=pig(),method="REML")
c <- gam(y~x+offset(o),data=d[d$w>0,],weights=w,family=pig(),method="REML")
near(coef(a),coef(b)); near(coef(a),coef(c))
near(a$family$getTheta(),b$family$getTheta()); near(a$family$getTheta(),c$family$getTheta())
near(residuals(a,"response"),d$y-fitted(a))
near(residuals(a,"pearson"),sqrt(d$w)*(d$y-fitted(a))/sqrt(a$family$variance(fitted(a))))
# Null deviance retains offset and the fitted global dispersion.
f <- reml
opt <- optimize(function(z) sum(f$family$dev.resids(d$y,exp(z+d$o),rep(1,n))),c(-10,10),tol=1e-9)
near(f$null.deviance,opt$objective)
pr <- predict(f,newdata=d[1:12,],type="response",se.fit=TRUE)
lp <- predict(f,newdata=d[1:12,],type="link",se.fit=TRUE)
near(pr$se.fit,pr$fit*lp$se.fit)
for (x in list(NA_real_,Inf,c(1,2),"a")) error(pig(x),"theta")
error(pig(link="identity"),"log link")
for (yy in list(c(-1,1),c(.5,1),c(Inf,1),c(0,0)))
  error(gam(y~1,data=data.frame(y=yy),family=pig()),"responses|positive count")
for (w in list(c(-1,1),c(Inf,1),c(0,0)))
  error(gam(y~1,data=data.frame(y=c(1,2)),weights=w,family=pig()),"weights|positive count")
error(gam(form,data=d,family=pig(),method="REML",scale=2),"scale is fixed")

# Exercise genuine binning (fewer discretization points than unique x),
# and ordinary bam block processing, at a moderate data size.
set.seed(27)
large <- data.frame(x=runif(1500))
large$y <- pig(.6)$rd(exp(.2 + sin(2*pi*large$x)))
ordinary <- bam(y~s(x,k=6),data=large,family=pig(.6),method="fREML",
                scale=1,chunk.size=200,nthreads=1)
binned <- bam(y~s(x,k=6),data=large,family=pig(.6),method="fREML",
              discrete=50,scale=1,nthreads=1)
stopifnot(isTRUE(ordinary$converged),isTRUE(binned$converged))
near(predict(ordinary,large,type="response"),predict(binned,large,type="response"),.02)
near(as.numeric(logLik(binned)),sum(pmf(large$y,fitted(binned),log(.6))))
