library("inlaws")
near <- function(a,b,tol=1e-5) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch(expr,error=identity),"error"))
eval_ll <- getFromNamespace(".gen_gamma_eval","inlaws")
set.seed(23); n <- 450L
d <- data.frame(x=runif(n),z=runif(n),a=runif(n,-.1,.1),b=runif(n,-.1,.1),c=runif(n,-.05,.05))
f <- gen_gammalss()
d$y <- f$rd(cbind(exp(.5+.5*sin(6*d$x)+d$a),exp(-.7+.3*d$z+d$b),exp(.2+.1*d$x+d$c)))
form <- list(y~s(x,k=5)+offset(a),~s(z,k=4)+offset(b),~x+offset(c))
check <- function(fit,method) {
  stopifnot(identical(fit$method,method),identical(fit$outer.info$conv,"full convergence"),
    all(is.finite(coef(fit))),is.finite(fit$gcv.ubre),
    identical(dim(fitted(fit)),c(n,3L)),all(fitted(fit)>0),
    identical(colnames(fitted(fit)),c("mu","sigma","r")))
  p <- predict(fit,type="response",se.fit=TRUE); l <- predict(fit,type="link",se.fit=TRUE)
  near(p$fit,fitted(fit)); near(p$fit,exp(l$fit)); near(p$se.fit,p$fit*l$se.fit)
  near(logLik(fit),sum(eval_ll(d$y,l$fit)$ell[,1]),1e-8)
  near(AIC(fit),-2*as.numeric(logLik(fit))+2*attr(logLik(fit),"df"))
  stopifnot(all(is.finite(residuals(fit))),is.finite(summary(fit)$dev.expl))
  dist <- fit$family$gratia_diagnostics(p$fit,fit$prior.weights,1,fit$y,l$fit)
  near(dist$mean,p$fit[,1]); near(residuals(fit,type="response"),d$y-dist$mean)
  near(residuals(fit,type="pearson"),(d$y-dist$mean)/sqrt(dist$variance))
}
a <- gam(form,data=d,family=f,method="REML"); check(a,"REML")
b <- gam(form,data=d,family=f,method="REML",optimizer=c("outer","bfgs")); check(b,"REML")
near(coef(a),coef(b),.001)
for(method in c("NCV","QNCV")) check(gam(form,data=d,family=f,method=method),method)
# Optimized grouped criteria, with prediction rows reordered within groups.
nei <- list(a=seq_len(n),ma=seq.int(3L,n,3L),
  d=as.vector(apply(matrix(seq_len(n),3L),2L,rev)),md=seq.int(3L,n,3L),jackknife=FALSE)
for(method in c("NCV","QNCV")) check(gam(form,data=d,family=f,method=method,nei=nei),method)
# Scale and transformed shape can both have smooth predictors.
smooth <- gam(list(y~s(x,k=5),~s(z,k=4),~s(x,k=4)),data=d,family=f,method="REML")
stopifnot(smooth$method=="REML",all(is.finite(coef(smooth))))
# Replication weights and independent unpenalized likelihood maximization.
d$w <- rep(c(0,1,2),length.out=n)
forms <- list(y~x,~z,~1)
u <- gam(forms,data=d,weights=w,family=f,method="REML")
v <- gam(forms,data=d[rep(seq_len(n),d$w),],family=f,method="REML")
near(coef(u),coef(v)); near(logLik(u),logLik(v))
objective <- function(b) {
  eta <- cbind(b[1]+b[2]*d$x,b[3]+b[4]*d$z,rep(b[5],n))
  -sum(d$w*eval_ll(d$y,eta)$ell[,1])
}
opt <- optim(coef(u)+.01,objective,method="BFGS",control=list(reltol=1e-12))
near(coef(u),opt$par,2e-4); near(logLik(u),-opt$value,1e-7)
# Negative shape preserves a mean but need not preserve a variance.
copy <- a; copy$fitted.values[,3] <- .4
stopifnot(all(is.na(residuals(copy,type="pearson"))),all(is.finite(residuals(copy,type="deviance"))))
# Missing rows, predictor offsets and serialization stay aligned.
d$x[1] <- NA
missing <- gam(form,data=d,family=f,method="REML",na.action=na.exclude)
stopifnot(nobs(missing)==n-1L,length(residuals(missing))==n,is.na(residuals(missing)[1]))
copy <- unserialize(serialize(a,NULL)); nd <- d[2:8,c("x","z","a","b","c")]
near(predict(copy,nd,type="response"),predict(a,nd,type="response"))
for(bad_weights in list(rep(-1,n),rep(0,n))) fails(gam(forms,data=d,weights=bad_weights,family=f))
for(y in list(rep(0,n),rep(-1,n),rep(Inf,n))) {
  bad <- d; bad$y <- y; fails(gam(forms,data=bad,family=f))
}
fails(gam(list(y~1,~1),data=d,family=f))
fails(gam(list(y~1,~1,~1,1+3~z-1),data=d,family=f))
# bam routes must fail explicitly, rather than silently fitting another model.
for(discrete in c(FALSE,TRUE)) for(method in c("fREML","NCV"))
  fails(suppressWarnings(mgcv::bam(form,data=d,family=f,method=method,discrete=discrete)))

# Both nuisance predictors may be constants with either REML optimizer.
constant_form <- list(y ~ s(x, k = 5), ~ 1, ~ 1)
constant <- gam(constant_form, data = d, family = f, method = "REML")
constant_bfgs <- gam(constant_form, data = d, family = f, method = "REML",
  optimizer = c("outer", "bfgs"))
near(coef(constant), coef(constant_bfgs), .001)
