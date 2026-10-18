library("inlaws")
near <- function(x, y, tol = 1e-7) stopifnot(length(x) == length(y), max(abs(as.numeric(x)-as.numeric(y))) < tol)
fails <- function(expr, pattern) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
derivatives <- getFromNamespace(".bi_derivatives", "inlaws")
prediction <- getFromNamespace(".bi_prediction", "inlaws")

# Independent probability/density integration and differentiation for every mode.
for (mode in c("both", "zero", "one")) {
  f <- zoibetals(mode); ep <- f$endpoints; q <- f$nlp
  e <- c(.3, log(9), -.7, .4)[seq_len(q)]
  prob <- setNames(c(0, 1, 0), c("p0", "pI", "p1"))
  prob[ep] <- exp(e[-(1:2)]); prob <- prob/sum(prob)
  m <- plogis(e[1]); ph <- exp(e[2])
  den <- function(x) prob["pI"]*dbeta(x, m*ph, (1-m)*ph)
  near(prob["p0"]+prob["p1"]+integrate(den, 0, 1)$value, 1)
  mn <- prob["p1"]+integrate(function(x) x*den(x), 0, 1)$value
  vr <- prob["p1"]+integrate(function(x) x*x*den(x), 0, 1)$value-mn^2
  et <- matrix(e, 1)
  d <- f$gratia_diagnostics(et, 1, 1, NULL, NULL)
  near(d$mean, mn); near(d$variance, vr)
  for (y in c(0, .17, .83, 1)) {
    if ((y == 0 && mode == "one") || (y == 1 && mode == "zero")) next
    z <- derivatives(y, et, 1.7, ep, f$tri$i2)
    expected <- if (y == 0) log(prob["p0"]) else if (y == 1) log(prob["p1"]) else log(den(y))
    near(z$l, 1.7*expected)
    step <- 1e-5
    for (j in seq_len(q)) {
      plus <- minus <- et; plus[j] <- plus[j]+step; minus[j] <- minus[j]-step
      zp <- derivatives(y, plus, 1.7, ep, f$tri$i2)
      zm <- derivatives(y, minus, 1.7, ep, f$tri$i2)
      near(z$l1[j], (zp$l-zm$l)/(2*step), 2e-6)
      near(z$l2[f$tri$i2[j, ]], (zp$l1-zm$l1)/(2*step), 2e-6)
      np <- prediction(plus, ep)$fit; nm <- prediction(minus, ep)$fit
      near(prediction(et, ep)$gradient[, , j], (np-nm)/(2*step), 2e-6)
    }
  }
  # Endpoint jumps and upper/lower tails, with NA propagation.
  near(d$cdf(0), prob["p0"]); near(d$cdf_left(0), 0)
  near(d$cdf(1), 1); near(d$cdf_left(1), 1-prob["p1"])
  for (at in c(-Inf, -.1, 0, .3, 1, 1.1, Inf)) {
    near(d$cdf(at)+d$cdf(at, lower_tail=FALSE), 1)
    near(d$cdf_left(at)+d$cdf_left(at, lower_tail=FALSE), 1)
  }
  stopifnot(is.na(d$cdf(NA_real_)), is.na(d$cdf_left(NA_real_)), is.na(d$qf(NA_real_)))
  near(d$qf(c(0)), 0); near(d$qf(1), 1)
  for (u in c(.05, .3, .6, .95)) {
    at <- d$qf(u)
    stopifnot(d$cdf_left(at) <= u+1e-10, d$cdf(at) >= u-1e-10)
    near(d$qf(log(u), log_p=TRUE), at)
    near(d$qf(log1p(-u), log_p=TRUE, lower_tail=FALSE), at)
  }
  # Degenerate observations do not define a finite model, but their likelihood
  # derivatives at finite parameter values are still well defined.
  z <- derivatives(c(if (mode == "one") .3 else 0, .2, 1), rbind(e, e, rep(Inf,q)), c(1,1,0), ep, f$tri$i2)
  stopifnot(is.finite(z$l), all(z$l1[3, ] == 0), all(z$l2[3, ] == 0))
  set.seed(271)
  en <- matrix(rep(e, each=30000), 30000, q)
  draw <- f$rd(en)
  near(mean(draw), mn, .008); near(var(draw), vr, .003)
  near(mean(draw == 0), prob["p0"], .008); near(mean(draw == 1), prob["p1"], .008)
  ds <- f$gratia_diagnostics(en, 1, 1, draw, en)
  pit <- runif(length(draw), ds$cdf_left(draw), ds$cdf(draw))
  near(mean(pit), .5, .008); near(var(pit), 1/12, .003)
  stopifnot(max(abs(sort(pit)-ppoints(length(pit)))) < .015)
  copy <- unserialize(serialize(f, NULL))
  near(copy$cdf(.4, et), f$cdf(.4, et))
}
# Direct upper-tail evaluation must retain small endpoint masses.
f <- zoibetals(); et <- matrix(c(0, log(8), -40, -40), 1)
d <- f$gratia_diagnostics(et, 1, 1, NULL, NULL)
near(d$cdf_left(1, log_p=TRUE, lower_tail=FALSE), -40, 1e-12)

set.seed(104)
n <- 650
dat <- data.frame(x=runif(n,-1,1), z=runif(n,-1,1))
fits <- list()
for (mode in c("both", "zero", "one")) {
  f <- zoibetals(mode); q <- f$nlp; ep <- f$endpoints
  eta <- cbind(.2+dat$x, 2+.8*dat$z, -.6+.7*dat$z, -.9-.5*dat$z)[,seq_len(q)]
  dat$y <- f$rd(eta)
  forms <- c(list(y~x, ~z), rep(list(~z), q-2L))
  b <- gam(forms, family=f, data=dat, method="REML")
  fits[[mode]] <- b
  stopifnot(b$method == "REML", all(is.finite(coef(b))))
  # Independent beta maximum likelihood on the interior (no REML nuisance
  # parameter estimation or separate smoothing selections in this comparison).
  mid <- dat$y > 0 & dat$y < 1
  A <- model.matrix(~x, dat[mid, ]); B <- model.matrix(~z, dat[mid, ])
  beta_ll <- function(a) {
    mu <- plogis(A%*%a[1:2]); phi <- exp(B%*%a[3:4])
    -sum(dbeta(dat$y[mid], mu*phi, (1-mu)*phi, log=TRUE))
  }
  opt <- optim(coef(b)[1:4]+.02, beta_ll, method="BFGS",
               control=list(reltol=1e-12, maxit=1000))
  near(coef(b)[1:4], opt$par, 3e-5)
  near(beta_ll(coef(b)[1:4]), opt$value, 1e-9)
  if (q == 3L) {
    dat$member <- dat$y == if (mode == "zero") 0 else 1
    independent <- glm(member~z, data=dat, family=binomial())
    near(tail(coef(b),2), coef(independent), 1e-6)
    membership_ll <- as.numeric(logLik(independent))
  } else {
    Z <- model.matrix(~z, dat)
    membership <- function(a) {
      pp <- exp(cbind(Z%*%a[1:2], 0, Z%*%a[3:4])); pp <- pp/rowSums(pp)
      cl <- ifelse(dat$y == 0, 1, ifelse(dat$y == 1, 3, 2))
      -sum(log(pp[cbind(seq_len(n),cl)]))
    }
    mem <- optim(tail(coef(b),4)+.02, membership, method="BFGS", control=list(reltol=1e-12))
    near(tail(coef(b),4), mem$par, 3e-5)
    membership_ll <- -mem$value
  }
  near(as.numeric(logLik(b)), -opt$value+membership_ll, 1e-8)
  pr <- predict_zoibetals(b,type="all",se.fit=TRUE)
  stopifnot(identical(colnames(pr$fit), c("mu","phi","p0","pI","p1","mean","variance")),
            all(is.finite(pr$se.fit)), is.na(b$deviance), is.na(b$null.deviance))
  near(rowSums(pr$fit[,3:5]),rep(1,n))
  near(predict_zoibetals(b), pr$fit[,6])
  near(predict_zoibetals(b,type="variance"), pr$fit[,7])
  near(residuals(b,type="response"), dat$y-pr$fit[,6])
  near(residuals(b,type="pearson"), (dat$y-pr$fit[,6])/sqrt(pr$fit[,7]))
  if (mode != "both") {
    disabled <- if (mode == "zero") "p1" else "p0"
    stopifnot(all(pr$fit[,disabled] == 0), all(pr$se.fit[,disabled] == 0))
  }
  nd <- dat[1:4,]
  near(predict_zoibetals(b,nd,type="all",block.size=1),pr$fit[1:4,])
  # Full covariance delta method, checked independently by coefficient differences.
  J <- matrix(0,4,length(coef(b))); h <- 1e-5
  for (j in seq_along(coef(b))) {
    bp <- bm <- b; bp$coefficients[j] <- coef(b)[j]+h; bm$coefficients[j] <- coef(b)[j]-h
    J[,j] <- (predict_zoibetals(bp,nd)-predict_zoibetals(bm,nd))/(2*h)
  }
  near(pr$se.fit[1:4,6], sqrt(rowSums((J%*%vcov(b))*J)), 2e-6)
  # Likelihood weights equal replication; fractional weights are also supported.
  dat$w <- rep(c(0,1,2,3),length.out=n)
  w <- gam(forms,family=f,data=dat,weights=w,method="REML")
  r <- gam(forms,family=f,data=dat[rep(seq_len(n),dat$w),],method="REML")
  near(coef(w),coef(r),2e-6); near(logLik(w),logLik(r),2e-7)
  frac <- gam(forms,family=f,data=dat,weights=dat$w/2,method="REML")
  near(coef(w),coef(frac),2e-6)
  # Offsets in every predictor: algebraically equivalent to coefficient shifts.
  dat$o1 <- .2*dat$x; dat$o2 <- -.3*dat$z; dat$o3 <- .4*dat$z; dat$o4 <- -.2*dat$z
  of <- c(list(y~x+offset(o1),~z+offset(o2)),list(~z+offset(o3)))
  if (q == 4L) of <- c(of,list(~z+offset(o4)))
  bo <- gam(of,family=f,data=dat,method="REML")
  near(predict_zoibetals(bo,type="all"),pr$fit,2e-5)
  dm <- dat; dm$y[2] <- NA; dm$z[3] <- NA
  missing <- gam(forms,family=f,data=dm,na.action=na.exclude,method="REML")
  pm <- predict_zoibetals(missing,type="all")
  stopifnot(nrow(missing$model) == n-2, nrow(pm) == n, all(is.na(pm[2:3,])))
  # Selection is assessed separately from the unpenalized factorization above.
  smooth <- c(list(y~s(x,k=5),~s(z,k=5)),rep(list(~s(z,k=5)),q-2L))
  bs <- gam(smooth,family=f,data=dat,method="REML",optimizer="efs")
  stopifnot(all(is.finite(bs$sp)), all(bs$sp>0), bs$outer.info$conv == "full convergence",
            all(is.finite(predict_zoibetals(bs,se.fit=TRUE)$se.fit)))
  # All-fixed smoothing parameters and shared coefficients use the same Hessian.
  G <- gam(smooth,family=f,data=dat,method="REML",fit=FALSE)
  fixed <- gam(G=G,sp=rep(.5,q),method="REML")
  stopifnot(all(is.finite(coef(fixed))))
  fails(gam(G=G,sp=c(.5,rep(-1,q-1)),method="REML"),"all smoothing parameters")
}
b <- fits$both; dat <- b$model
shared <- gam(list(y~1,~z,~z,~z,1+3~x-1),data=dat,family=zoibetals(),method="REML")
stopifnot(all(is.finite(coef(shared))),all(is.finite(vcov(shared,sandwich=TRUE))))
stopifnot(all(is.finite(predict_zoibetals(shared,se.fit=TRUE)$se.fit)))
fails(residuals(b),"deviance residuals are undefined")
fails(zoibetals("invalid"),"arg")
fails(predict_zoibetals(lm(1:5~1)),"gam fitted")
fails(b$family$rd(matrix(0,2,7)),"link-scale")
fails(b$family$ll(b$y,model.matrix(b),coef(b),b$prior.weights,b$family,deriv=2),"EFS REML")
fails(b$family$ll(b$y,model.matrix(b),coef(b),b$prior.weights,b$family,ncv=TRUE),"EFS REML")
for (mode in c("both","zero","one")) {
  f <- zoibetals(mode)
  v <- f$validate
  yy <- c(0,.2,.4,.7,1)
  if (mode == "zero") yy <- yy[-5]
  if (mode == "one") yy <- yy[-1]
  v(yy,rep(1,length(yy)))
  fails(v(c(yy,Inf),rep(1,length(yy)+1)),"finite numeric")
  fails(v(c(yy,-.1),rep(1,length(yy)+1)),"responses in")
  fails(v(yy,rep(0,length(yy))),"positive total")
  fails(v(yy,rep(-1,length(yy))),"nonnegative")
  fails(v(c(0,.2,.2,1),rep(1,4)),if (mode == "both") "distinct" else "excluded")
  fails(v(c(.2,.4,.7),rep(1,3)),"active endpoint")
  for (end in f$endpoints) {
    w <- rep(1,length(yy)); w[yy == if (end == "p0") 0 else 1] <- 0
    fails(v(yy,w),"active endpoint")
  }
  if (mode != "both") fails(v(c(yy,if(mode=="zero") 1 else 0),c(rep(1,length(yy)),0)),"excluded")
}
# Sparse, but present, endpoints are retained in the requested model.
ds <- data.frame(y=c(0,1,seq(.1,.9,length.out=98)))
sparse <- gam(list(y~1,~1,~1,~1),data=ds,family=zoibetals(),method="REML")
near(predict_zoibetals(sparse,type="probabilities")[1,],c(.01,.98,.01),1e-6)
# Smooth EFS offsets must agree between fitting likelihood and prediction.
dat <- fits$both$model
for (j in 1:4) dat[[paste0("off",j)]] <- sin(dat$x*j)/5
bs <- gam(list(y~s(x,k=5)+offset(off1),~z+offset(off2),
               ~s(z,k=5)+offset(off3),~z+offset(off4)),
          data=dat,family=zoibetals(),method="REML")
near(predict(bs,type="link"),bs$linear.predictors,1e-8)
reference <- prediction(predict(bs,type="link"),c("p0","p1"))$fit
ii <- dat$y>0 & dat$y<1
ll <- log(reference[,"p0"]); ll[dat$y==1] <- log(reference[dat$y==1,"p1"])
ll[ii] <- log(reference[ii,"pI"])+dbeta(dat$y[ii],
  reference[ii,"mu"]*reference[ii,"phi"],(1-reference[ii,"mu"])*reference[ii,"phi"],log=TRUE)
near(sum(ll),as.numeric(logLik(bs)),1e-8)
fails(gam(list(y~s(x,k=5,id=1)+s(z,k=5,id=1),~1,~1,~1),
          data=dat,family=zoibetals(),method="REML"),"all smoothing parameters")
# QQ simulation uses the link-scale callback and complete Pearson discrepancy.
pdf(file=tempfile(fileext=".pdf"))
qq.gam(fits$both,type="pearson",rep=5)
invisible(dev.off())
# Native factory contracts are mandatory; integration checks use the newer
# extensible gratia interface only when that interface is installed.
if (requireNamespace("gratia",quietly=TRUE) &&
    exists("response_distribution",envir=asNamespace("gratia"),inherits=FALSE)) {
  for (b in fits) {
    pit <- gratia::quantile_residuals(b,type="pit",seed=91)
    qr <- gratia::quantile_residuals(b,type="quantile",seed=91)
    stopifnot(length(pit)==length(b$y),all(is.finite(qr)),all(pit>=0 & pit<=1))
    near(qnorm(pit),qr,1e-10)
    dr <- gratia::dharma_residuals(b,n=20,seed=92)
    stopifnot(length(dr)==length(b$y),all(dr>=0 & dr<=1))
  }
}
# A single parameter row can evaluate an entire CDF/quantile grid.
for (mode in c("both","zero","one")) {
  f <- zoibetals(mode); et <- matrix(c(.4,log(7),-.8,.3)[seq_len(f$nlp)],1)
  d <- f$gratia_diagnostics(et,1,1,NULL,NULL)
  grid <- c(-Inf,0,.2,.7,1,Inf)
  near(d$cdf(grid),vapply(grid,d$cdf,numeric(1)))
  near(d$cdf_left(grid),vapply(grid,d$cdf_left,numeric(1)))
  u <- c(0,.1,.3,.7,.99,1)
  near(d$qf(u),vapply(u,d$qf,numeric(1)))
  stopifnot(length(d$cdf(numeric()))==0,length(d$qf(numeric()))==0)
  # Additional derivative regimes: small shapes, strong membership logits,
  # concentrated beta responses, and predictors away from a zero intercept.
  for (ee in list(c(-3,log(.6),2,-4),c(2,log(80),-3,1))) {
    e <- matrix(ee[seq_len(f$nlp)],1); h <- 1e-5
    z <- derivatives(.31,e,.6,f$endpoints,f$tri$i2)
    for (j in seq_len(f$nlp)) {
      plus <- minus <- e; plus[j] <- plus[j]+h; minus[j] <- minus[j]-h
      zp <- derivatives(.31,plus,.6,f$endpoints,f$tri$i2)
      zm <- derivatives(.31,minus,.6,f$endpoints,f$tri$i2)
      near(z$l1[j],(zp$l-zm$l)/(2*h),2e-6)
      near(z$l2[f$tri$i2[j,]],(zp$l1-zm$l1)/(2*h),2e-6)
    }
  }
}
f <- zoibetals(); et <- matrix(c(0,log(8),-40,-40),1)
d <- f$gratia_diagnostics(et,1,1,NULL,NULL)
# The log CDF left of one retains the tiny atom even when 1-p1 rounds to one.
stopifnot(d$cdf_left(1,log_p=TRUE)<0,
          abs(d$cdf_left(1,log_p=TRUE)/(-exp(-40))-1)<1e-12)
