library("inlaws")
near <- function(a,b,tol=1e-5) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch(expr,error=identity),"error"))
set.seed(46); n <- 600L
d <- data.frame(x=runif(n),z=runif(n),a=runif(n,-.1,.1),b=runif(n,-.1,.1),c=runif(n,-.1,.1))
f <- scatlss(); d$y <- f$rd(cbind(sin(6*d$x)+d$a,exp(-.5+.2*d$z+d$b),3+exp(.3+d$x+d$c)))
form <- list(y~s(x,k=6)+offset(a),~z+offset(b),~x+offset(c))
complete_data <- d
check <- function(fit,method) {
  stopifnot(identical(fit$method,method),identical(fit$outer.info$conv,"full convergence"),
    all(is.finite(coef(fit))),is.finite(fit$gcv.ubre),
    identical(dim(fitted(fit)),c(n,3L)),all(fitted(fit)[,2]>0),all(fitted(fit)[,3]>=3),
    identical(colnames(fitted(fit)),c("mu","sigma","nu")))
  stopifnot(all(is.finite(fit$Vp)), all(is.finite(fit$Ve)))
  if (method %in% c("NCV", "QNCV"))
    stopifnot(identical(fit$optimizer, c("outer", "bfgs")))
  p <- predict(fit,type="response",se.fit=TRUE); l <- predict(fit,type="link",se.fit=TRUE)
  near(p$fit,fitted(fit)); near(p$fit,cbind(l$fit[,1],exp(l$fit[,2]),3+exp(l$fit[,3])))
  near(p$se.fit,l$se.fit*cbind(1,exp(l$fit[,2]),exp(l$fit[,3])))
  near(logLik(fit),sum(fit$prior.weights*(dt((fit$y-p$fit[,1])/p$fit[,2],p$fit[,3],log=TRUE)-log(p$fit[,2]))),1e-8)
  near(AIC(fit),-2*as.numeric(logLik(fit))+2*attr(logLik(fit),"df"))
  dist <- fit$family$gratia_diagnostics(p$fit,fit$prior.weights,1,fit$y,l$fit)
  near(residuals(fit,type="response"),fit$y-dist$mean)
  near(residuals(fit,type="pearson"),(fit$y-dist$mean)*sqrt(fit$prior.weights/dist$variance))
  stopifnot(all(is.finite(residuals(fit))),is.finite(summary(fit)$dev.expl))
}
a <- gam(form,data=d,family=f,method="REML"); check(a,"REML")
b <- gam(form,data=d,family=f,method="REML",optimizer=c("outer","bfgs")); check(b,"REML")
near(coef(a),coef(b),.003)
for(method in c("NCV","QNCV")) check(gam(form,data=d,family=f,method=method),method)
nei <- list(a=seq_len(n),ma=seq.int(3L,n,3L),
  d=as.vector(apply(matrix(seq_len(n),3L),2L,rev)),md=seq.int(3L,n,3L),jackknife=FALSE)
for(method in c("NCV","QNCV")) check(gam(form,data=d,family=f,method=method,nei=nei),method)
# NCV diagnostics requested during REML also bypass the unrelated extended guard.
r <- gam(form,data=d,family=f,method="REML",nei=nei); check(r,"REML")
stopifnot(is.finite(r$NCV))
for(opt in c("newton","bfgs")) {
  g <- gam(list(y~s(x,k=6),~1,~1),data=d,family=f,method="REML",optimizer=c("outer",opt))
  check(g,"REML"); stopifnot(length(unique(fitted(g)[,3]))==1L)
}
smooth <- gam(list(y~s(x,k=5),~s(z,k=4),~s(x,k=4)),data=d,family=f,method="REML")
check(smooth,"REML")
d$w <- rep(c(0,1,2),length.out=n)
forms <- list(y~x+offset(a),~z+offset(b),~1+offset(c))
u <- gam(forms,data=d,weights=w,family=f,method="REML")
v <- gam(forms,data=d[rep(seq_len(n),d$w),],family=f,method="REML")
near(coef(u),coef(v)); near(logLik(u),logLik(v))
objective <- function(beta) {
  mu <- beta[1]+beta[2]*d$x+d$a
  sig <- exp(beta[3]+beta[4]*d$z+d$b); df <- 3+exp(beta[5]+d$c)
  -sum(d$w*(dt((d$y-mu)/sig,df,log=TRUE)-log(sig)))
}
opt <- optim(coef(u)+.01,objective,method="BFGS",control=list(reltol=1e-12))
near(coef(u),opt$par,5e-4); near(logLik(u),-opt$value,1e-8)
# Entirely constant model and custom lower bound, independent likelihood.
c0 <- gam(list(y~1,~1,~1),data=d,family=scatlss(2.2),method="REML")
stopifnot(all(fitted(c0)[,3]>2.2))
near(logLik(c0),sum(dt((d$y-fitted(c0)[,1])/fitted(c0)[,2],fitted(c0)[,3],log=TRUE)-log(fitted(c0)[,2])),1e-10)
# Missing rows in responses, predictors, offsets and weights remain aligned
# after subset selection. Compare with the explicitly cleaned model frame.
d$x[4] <- NA; d$y[7] <- NA; d$c[8] <- NA; d$w[9] <- NA
subset <- seq_len(n) %% 3 != 0
miss <- gam(forms,data=d,weights=w,subset=subset,na.action=na.exclude,family=f,method="REML")
keep <- subset & complete.cases(d)
clean <- gam(forms,data=d[keep,],weights=w,family=f,method="REML")
near(coef(miss),coef(clean)); near(logLik(miss),logLik(clean))
stopifnot(length(miss$y)==sum(keep),nobs(miss)==sum(keep & d$w>0),length(residuals(miss))==sum(subset),sum(is.na(residuals(miss)))==sum(subset & !complete.cases(d)))
copy <- unserialize(serialize(a,NULL)); nd <- d[20:25,]
near(predict(copy,nd,type="response"),predict(a,nd,type="response"))
for(bad_weights in list(rep(-1,n),rep(0,n),rep(Inf,n))) fails(gam(forms,data=d,weights=bad_weights,family=f))
for(y in list(rep(Inf,n),matrix(1,n,2))) { bad <- d; bad$y <- y; fails(gam(forms,data=bad,family=f)) }
fails(gam(list(y~1,~1),data=d,family=f))
fails(gam(list(y~1,~1,~1,1+3~z-1),data=d,family=f))
for(discrete in c(FALSE,TRUE)) fails(suppressWarnings(bam(form,data=d,family=f,method="fREML",discrete=discrete)))
# Default leave-one-out neighbourhoods are constructed after row omission.
missing_form <- list(y~s(x,k=6)+offset(a),~z+offset(b),~1+offset(c))
for(method in c("NCV","QNCV")) {
  omitted <- gam(missing_form,data=d,weights=w,subset=subset,na.action=na.exclude,family=f,method=method)
  explicit <- gam(missing_form,data=d[keep,],weights=w,family=f,method=method)
  stopifnot(omitted$method==method,explicit$method==method,
    identical(omitted$outer.info$conv,"full convergence"))
  near(coef(omitted),coef(explicit),1e-6); near(logLik(omitted),logLik(explicit),1e-8)
}
# mgcv 1.9-4 nanei() corrupts cumulative boundaries for partially omitted
# explicit groups. Reject all affected request routes, preserving REML alone.
if(packageVersion("mgcv")=="1.9.4") {
  nsub <- sum(subset)
  groups <- list(a=seq_len(nsub),ma=pmin(seq.int(3L,nsub+2L,3L),nsub),
    d=seq_len(nsub),md=pmin(seq.int(3L,nsub+2L,3L),nsub),jackknife=FALSE)
  for(method in c("NCV","QNCV","REML")) {
    err <- tryCatch(gam(missing_form,data=d,weights=w,subset=subset,na.action=na.exclude,
      family=f,method=method,nei=groups),error=identity)
    stopifnot(inherits(err,"error"),grepl("mgcv 1.9-4",conditionMessage(err),fixed=TRUE),
      grepl("then construct nei",conditionMessage(err),fixed=TRUE))
  }
}
# Rebuilding explicit groups on complete cases restores both criteria.
nclean <- sum(keep)
groups <- list(a=seq_len(nclean),ma=pmin(seq.int(3L,nclean+2L,3L),nclean),
  d=seq_len(nclean),md=pmin(seq.int(3L,nclean+2L,3L),nclean),jackknife=FALSE)
for(method in c("NCV","QNCV")) {
  complete <- gam(missing_form,data=d[keep,],weights=w,family=f,method=method,nei=groups)
  stopifnot(complete$method==method,identical(complete$outer.info$conv,"full convergence"))
}
# An entirely parametric model has no smoothness parameters. Its NCV/QNCV
# scores work with jackknife disabled; 1.9-4's empty-penalty jackknife fails.
nd <- d[keep,]; nr <- seq_len(nrow(nd))
plain_nei <- list(a=nr,ma=nr,d=nr,md=nr,jackknife=FALSE)
for(method in c("NCV","QNCV")) {
  plain <- gam(forms,data=nd,weights=w,family=f,method=method,nei=plain_nei)
  stopifnot(plain$method==method,is.finite(plain$gcv.ubre))
  if(packageVersion("mgcv")=="1.9.4") {
    err <- tryCatch(gam(forms,data=nd,weights=w,family=f,method=method),error=identity)
    stopifnot(inherits(err,"error"),grepl("empty-penalty",conditionMessage(err),fixed=TRUE))
  }
}
# Shared preprocessing allows correctly adjusted whole-group omissions,
# preserving the original optimizer hook and all three predictor offsets.
grouped <- function(n) list(a=seq_len(n),ma=seq.int(3L,n,3L),
  d=as.vector(apply(matrix(seq_len(n),3L),2L,rev)),
  md=seq.int(3L,n,3L),jackknife=FALSE)
whole <- complete_data; whole$x[1:3] <- NA
clean_whole <- whole[complete.cases(whole),]
for(method in c("NCV","QNCV","REML")) {
  fit <- gam(form,data=whole,family=f,method=method,nei=grouped(nrow(whole)))
  rebuilt <- gam(form,data=clean_whole,family=f,method=method,nei=grouped(nrow(clean_whole)))
  for(model in list(fit,rebuilt)) stopifnot(model$method==method,
    identical(model$optimizer,c("outer","bfgs")),
    identical(model$outer.info$conv,"full convergence"),
    all(is.finite(model$Vp)),is.finite(model$gcv.ubre),is.finite(model$NCV))
  near(coef(fit),coef(rebuilt)); near(fit$NCV,rebuilt$NCV)
}
