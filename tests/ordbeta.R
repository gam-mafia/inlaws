library("inlaws")

near <- function(a, b, tol = 1e-7) {
  stopifnot(length(a) == length(b), all(is.finite(a)), all(is.finite(b)),
            max(abs(a-b)/(1+abs(a)+abs(b))) < tol)
}
near_difference <- function(a, plus, minus, h, tol = 3e-6) {
  fd <- (plus-minus)/(2*h)
  # Subtraction of large nearly equal derivatives limits attainable FD
  # accuracy, especially for cutpoint effects on beta-dominated curvature.
  rounding <- 64*.Machine$double.eps*(abs(plus)+abs(minus))/h
  stopifnot(all(is.finite(a)),all(is.finite(fd)),
            all(abs(a-fd) < tol*(1+abs(a)+abs(fd))+rounding))
}
expect_error <- function(expr, pattern) {
  e <- tryCatch(expr, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
make_family <- getFromNamespace(".ordbeta_family", "inlaws")
components <- getFromNamespace(".ordbeta_components", "inlaws")

# Independent response transformations: do not use the package's component
# or slope functions, so prediction and numerical Jacobian checks cannot share
# the implementation's probability or derivative mistakes.
response_values <- function(eta, par) {
  mu <- plogis(eta)
  p0 <- plogis(par[2]-eta)
  p1 <- plogis(eta-par[3])
  pI <- plogis(eta-par[2])-plogis(eta-par[3])
  cbind(mu = mu, p0 = p0, pI = pI, p1 = p1, mean = p1+pI*mu)
}
prediction_contract <- function(object, nd) {
  X <- predict(object, nd, type = "lpmatrix")
  eta <- drop(X %*% coef(object))+nd$off
  par <- object$family$getTheta(TRUE)
  expected <- response_values(eta, par)
  actual <- predict(object, nd, type = "response", se.fit = TRUE)
  stopifnot(identical(dim(actual$fit), c(nrow(nd), 5L)),
            identical(dim(actual$se.fit), dim(actual$fit)))
  near(actual$fit, expected)
  near(predict(object, nd, type = "link"), eta)
  near(predict(object, nd, type = "link", se.fit = TRUE)$se.fit,
       sqrt(rowSums((X %*% object$Vp)*X)))
  # Independently perturb coefficients, holding nuisance parameters fixed.
  h <- 1e-5
  J <- array(0, c(nrow(nd), 5L, ncol(X)))
  for (j in seq_len(ncol(X))) {
    bp <- bm <- coef(object); bp[j] <- bp[j]+h; bm[j] <- bm[j]-h
    J[, , j] <- (response_values(drop(X %*% bp)+nd$off, par)-
                  response_values(drop(X %*% bm)+nd$off, par))/(2*h)
  }
  for (k in 1:5) {
    jac <- matrix(J[, k, ], nrow(nd), ncol(X))
    near(actual$se.fit[, k], sqrt(rowSums((jac %*% object$Vp)*jac)), 1e-7)
  }
  for (size in c(1L, 5L, 1000L)) {
    blocked <- predict(object, nd, type = "response", se.fit = TRUE, block.size = size)
    near(blocked$fit, actual$fit); near(blocked$se.fit, actual$se.fit)
    near(predict(object, nd, type = "response", block.size = size), actual$fit)
    for (type in c("response", "conditional", "probabilities")) {
      columns <- switch(type, response = 5L, conditional = 1L, probabilities = 2:4)
      helper <- predict_ordbeta(object, nd, type = type, se.fit = TRUE, block.size = size)
      near(helper$fit, actual$fit[, columns]); near(helper$se.fit, actual$se.fit[, columns])
      near(predict_ordbeta(object, nd, type = type, block.size = size), helper$fit)
      if (type == "probabilities") {
        stopifnot(identical(colnames(helper$fit), c("p0", "pI", "p1")),
                  identical(colnames(helper$se.fit), colnames(helper$fit)))
      } else stopifnot(is.null(dim(helper$fit)), is.null(dim(helper$se.fit)))
    }
  }
  # Direct callback contract: supplied eta takes precedence, and X, beta, off,
  # and Vb are used without reaching into the fitted object.
  cb <- object$family$predict
  direct <- cb(object$family, se = TRUE, X = X, beta = coef(object),
               off = nd$off, Vb = object$Vp)
  near(direct$fit, actual$fit); near(direct$se.fit, actual$se.fit)
  stopifnot(identical(colnames(direct$fit), colnames(expected)),
            identical(colnames(direct$se.fit), colnames(expected)))
  supplied <- cb(object$family, se = TRUE, eta = eta, X = X,
                 beta = coef(object)+1, off = 100, Vb = 4*object$Vp)
  near(supplied$fit, actual$fit); near(supplied$se.fit, 2*actual$se.fit)
  near(cb(object$family, eta = eta)$fit, actual$fit)
  near(cb(object$family, X = X, beta = coef(object))$fit,
       response_values(drop(X %*% coef(object)), par))
  near(cb(object$family, X = X, beta = coef(object), off = .3)$fit,
       response_values(drop(X %*% coef(object))+.3, par))
  for (column in c(if (!is.finite(par[2])) 2L, if (!is.finite(par[3])) 4L)) {
    stopifnot(all(actual$fit[, column] == 0), all(actual$se.fit[, column] == 0))
  }
  if (object$family$mode == "beta") {
    stopifnot(all(actual$fit[, 3] == 1), all(actual$se.fit[, 3] == 0))
    near(actual$fit[, 1], actual$fit[, 5])
    near(actual$se.fit[, 1], actual$se.fit[, 5])
  }
}

# Check all entries of the derivative contract, including packed cross terms.
for (mode in c("full", "lower", "upper", "beta")) {
  theta <- switch(mode, full=c(log(9), -.8, log(2)), lower=c(log(9), -.8),
                  upper=c(log(9), 1.2), beta=log(9))
  f <- make_family(mode, theta)
  y <- switch(mode, full=c(0, .02, .4, .95, 1), lower=c(0, .02, .4, .95, .99),
              upper=c(.001, .02, .4, .95, 1), beta=c(.001, .02, .4, .95, .99))
  w <- c(.5, 2, 0, 1, 3)
  for (mu in list(c(.07, .2, .5, .8, .92), c(.001, .02, .5, .97, .999))) {
    h <- pmin(mu, 1-mu)*1e-4
    d <- f$Dd(y, mu, theta, w, 2)
    plus <- f$Dd(y, mu+h, theta, w, 2)
    minus <- f$Dd(y, mu-h, theta, w, 2)
    near(d$Dmu, (f$dev.resids(y, mu+h, w)-f$dev.resids(y, mu-h, w))/(2*h), 2e-6)
    for (pair in list(c("Dmu2", "Dmu"), c("Dmu3", "Dmu2"), c("Dmu4", "Dmu3"),
                     c("Dmuth", "Dth"), c("Dmu2th", "Dmuth"), c("Dmu3th", "Dmu2th"),
                     c("Dmuth2", "Dth2"), c("Dmu2th2", "Dmuth2"), c("EDmu3", "EDmu2")))
      near(d[[pair[1]]], (plus[[pair[2]]]-minus[[pair[2]]])/(2*h), 3e-6)
    nt <- length(theta)
    packed <- which(upper.tri(matrix(0, nt, nt), diag=TRUE), arr.ind=TRUE)
    packed <- packed[order(packed[,1], packed[,2]), , drop=FALSE]
    # A slightly larger theta step avoids subtraction loss when a large
    # mu derivative has only a small dependence on precision near a bound.
    h <- 1e-4
    for (j in seq_along(theta)) {
      tp <- tm <- theta; tp[j] <- tp[j]+h; tm[j] <- tm[j]-h
      plus <- f$Dd(y, mu, tp, w, 2); minus <- f$Dd(y, mu, tm, w, 2)
      near(d$Dth[,j], (f$dev.resids(y, mu, w, tp)-f$dev.resids(y, mu, w, tm))/(2*h), 2e-6)
      for (pair in list(c("Dmuth", "Dmu"), c("Dmu2th", "Dmu2"),
                       c("Dmu3th", "Dmu3"), c("EDmu2th", "EDmu2")))
        near_difference(d[[pair[1]]][,j], plus[[pair[2]]], minus[[pair[2]]], h)
      for (k in which(packed[,2] == j)) {
        i <- packed[k,1]
        for (pair in list(c("Dth2", "Dth"), c("Dmuth2", "Dmuth"), c("Dmu2th2", "Dmu2th")))
          near_difference(d[[pair[1]]][,k], plus[[pair[2]]][,i], minus[[pair[2]]][,i], h)
      }
    }
    stopifnot(all(d$EDmu2 >= 0))
    for (level in 0:1) {
      low <- f$Dd(y, mu, theta, w, level)
      for (field in names(low)) near(low[[field]], d[[field]])
    }
    one <- f$Dd(y[1], mu[1], theta, w[1], 2)
    stopifnot(identical(dim(one$Dth), c(1L, nt)))
    for (field in names(one)) near(one[[field]], if(is.matrix(d[[field]])) d[[field]][1,,drop=FALSE] else d[[field]][1])
  }
  # The observed curvature averages to the analytic Fisher information.
  m <- .37; p <- components(qlogis(m), f$getTheta(TRUE)); phi <- exp(theta[1])
  interior <- integrate(function(y) f$Dd(y, rep(m,length(y)), theta, rep(1,length(y)), 0)$Dmu2 *
                          dbeta(y, m*phi, (1-m)*phi), 0, 1, rel.tol=1e-6)$value
  obs <- interior*p[1,2]
  if (p[1,1] > 0) obs <- obs+p[1,1]*f$Dd(0,m,theta,1,0)$Dmu2
  if (p[1,3] > 0) obs <- obs+p[1,3]*f$Dd(1,m,theta,1,0)$Dmu2
  near(obs, f$Dd(.4,m,theta,1,0)$EDmu2, 1e-5)
}

# Distribution methods: atoms, continuous quantiles, moments and log tails.
set.seed(420)
for (mode in c("full", "lower", "upper", "beta")) {
  theta <- switch(mode, full=c(log(12), -1, log(2)), lower=c(log(12), -1),
                  upper=c(log(12), 1), beta=log(12))
  f <- make_family(mode,theta); mu <- .4; p <- components(qlogis(mu),f$getTheta(TRUE))
  near(sum(p),1)
  near(f$cdf(0,mu),p[1,1]); near(f$cdf(1,mu),1)
  near(f$cdf(-1,mu),0); near(f$cdf(2,mu),1)
  near(f$cdf(1-1e-10,mu),1-p[1,3],1e-6)
  probs <- p[1,1]+p[1,2]*c(.1,.5,.9)
  near(f$cdf(f$qf(probs,mu),mu),probs)
  near(f$qf(c(0,1),mu),c(0,1))
  if(p[1,1]>0) near(f$qf(p[1,1]/2,mu),0)
  if(p[1,3]>0) near(f$qf(1-p[1,3]/2,mu),1)
  draws <- f$rd(rep(mu,40000))
  near(mean(draws),p[1,3]+p[1,2]*mu,.008)
  near(var(draws),f$variance(mu),.008)
  near(mean(draws==0),p[1,1],.008); near(mean(draws==1),p[1,3],.008)
  near(exp(f$cdf(.3,mu,logp=TRUE)), f$cdf(.3,mu))
}
f <- make_family("upper",c(log(10),1))
stopifnot(is.finite(f$cdf(1e-100,.8,logp=TRUE)))

# Tight gaps and far-tail cutpoints exercise stable probability arithmetic.
for (cuts in list(c(-20,20),c(5,5.0001),c(-30,-29))) {
  theta <- c(log(4),cuts[1],log(diff(cuts)))
  f <- make_family("full",theta)
  m <- c(.03,.5,.97); y <- c(0,.4,1)
  dd <- f$Dd(y,m,theta,rep(1,3),2)
  stopifnot(all(vapply(dd,function(x) all(is.finite(x)),logical(1))))
  near(rowSums(components(qlogis(m),f$getTheta(TRUE))),rep(1,3))
  h <- 1e-4
  for(j in seq_along(theta)) {
    tp <- tm <- theta; tp[j] <- tp[j]+h; tm[j] <- tm[j]-h
    near_difference(dd$Dth[,j],f$dev.resids(y,m,1,tp),f$dev.resids(y,m,1,tm),h)
  }
}

# Independent simulation with substantial observations in every component.
set.seed(75)
n <- 350
d <- data.frame(x=runif(n,-1,1), off=runif(n,-.2,.2))
eta <- .2+1.1*d$x+d$off
mu <- plogis(eta); u <- runif(n)
d$y <- rbeta(n,mu*12,(1-mu)*12)
d$y[u < plogis(-1-eta)] <- 0
d$y[u > 1-plogis(eta-1.3)] <- 1
template <- ordbeta()
fit <- gam(y~x+offset(off),data=d,family=template,method="ML")
stopifnot(fit$converged,fit$scale==1, fit$family$n.theta==3)
stopifnot(is.null(attr(fit$family,"G")))
saved <- fit$family$getTheta(TRUE)
reuse <- gam(y~1,data=d,family=fit$family,method="ML")
near(fit$family$getTheta(TRUE),saved)
near(template$getTheta(),c(log(10),-1,log(2)))

# Independent objective, including every likelihood constant, for MLE check.
objective <- function(par, data) {
  eta <- par[1]+par[2]*data$x+data$off
  phi <- exp(par[3]); lo <- par[4]; hi <- lo+exp(par[5]); mu <- plogis(eta)
  y <- data$y; z <- y==0; o <- y==1; i <- !z & !o
  ll <- numeric(length(y))
  ll[z] <- plogis(lo-eta[z],log.p=TRUE); ll[o] <- plogis(eta[o]-hi,log.p=TRUE)
  ll[i] <- log(plogis(eta[i]-lo)-plogis(eta[i]-hi))+
    dbeta(y[i],mu[i]*phi,(1-mu[i])*phi,log=TRUE)
  -sum(ll)
}
par <- c(coef(fit),fit$family$getTheta())
ref <- optim(c(0,0,log(5),-.5,log(1)),objective,data=d,method="BFGS",control=list(reltol=1e-11),hessian=TRUE)
stopifnot(ref$convergence==0)
near(par,ref$par,2e-4)
near(as.numeric(logLik(fit)),-ref$value,1e-7)

# Nonuniform frequency weights are equivalent to replicated observations.
d$freq <- rep(1:3,length.out=nrow(d))
weighted <- gam(y~x+offset(off),data=d,weights=freq,family=ordbeta(),method="ML")
replicated <- gam(y~x+offset(off),data=d[rep(seq_len(nrow(d)),d$freq),],family=ordbeta(),method="ML")
near(coef(weighted),coef(replicated),1e-6)
near(as.numeric(logLik(weighted)),as.numeric(logLik(replicated)),1e-7)
near(sum(residuals(weighted,type="deviance")^2),deviance(weighted),1e-6)

# Reduced cases, zero weights, absent interiors, and family isolation.
for (mode in c("lower","upper","beta")) {
  sub <- d[switch(mode,lower=d$y<1,upper=d$y>0,beta=d$y>0 & d$y<1),]
  g <- suppressMessages(gam(y~x+offset(off),data=sub,family=template,method="ML"))
  stopifnot(g$converged,g$family$mode==mode,g$family$n.theta==if(mode=="beta") 1 else 2)
  natural <- g$family$getTheta(TRUE)
  stopifnot(natural[2]<natural[3],natural[1]>0)
  if(mode %in% c("upper","beta")) stopifnot(natural[2]==-Inf)
  if(mode %in% c("lower","beta")) stopifnot(natural[3]==Inf)
  near(fit$family$getTheta(TRUE),saved)
  prediction_contract(g, d[1:12, c("x", "off")])
  if(mode=="beta") {
    reference <- gam(y~x+offset(off),data=sub,family=mgcv::betar(),method="ML")
    near(coef(g),coef(reference),1e-5)
    near(g$family$getTheta(TRUE)[1],reference$family$getTheta(TRUE),1e-5)
    near(as.numeric(logLik(g)),as.numeric(logLik(reference)),1e-7)
  }
  d$w <- as.numeric(switch(mode,lower=d$y<1,upper=d$y>0,beta=d$y>0 & d$y<1))
  gw <- suppressMessages(gam(y~x+offset(off),data=d,weights=w,family=template,method="ML"))
  near(coef(gw),coef(g),1e-5)
  stopifnot(all(is.finite(gw$family$dev.resids(d$y,gw$fitted.values,d$w))))
}
expect_error(gam(y~x,data=d[d$y %in% c(0,1),],family=ordbeta(),method="REML"),"precision is unidentified")
expect_error(ordbeta(link="identity"),"logit")
expect_error(ordbeta(start=list(phi=-1)),"positive")
expect_error(ordbeta(start=list(cutpoints=c(1,0))),"ordered")
expect_error(ordbeta(start=list(other=1)),"named list")
expect_error(gam(y~x,data=d,family=ordbeta(start=list(cutpoints=c(-Inf,1))),method="ML"),"active starting")
bad <- d; bad$y[1] <- -1
expect_error(gam(y~x,data=bad,family=ordbeta(),method="ML"),"responses")
bad$y[1] <- Inf
expect_error(gam(y~x,data=bad,family=ordbeta(),method="ML"),"responses")
gstart <- gam(y~x+offset(off),data=d,family=ordbeta(start=list(phi=8,cutpoints=c(-1,1))),method="ML")
near(coef(gstart),coef(fit),1e-5)

# Smooths, prediction transformations/SEs, offsets and missing values.
set.seed(142)
smooth_data <- data.frame(x=runif(400))
truth <- sin(2*pi*smooth_data$x)
mu <- plogis(truth); u <- runif(400)
smooth_data$y <- rbeta(400,mu*15,(1-mu)*15)
smooth_data$y[u<plogis(-1-truth)] <- 0
smooth_data$y[u>1-plogis(truth-1)] <- 1
for (method in c("ML","REML")) {
  sm <- gam(y~s(x,k=7),data=smooth_data,family=ordbeta(),method=method)
  stopifnot(sm$converged,sqrt(mean((predict(sm)-truth)^2))<.3)
}
uc <- predict_ordbeta(sm,se.fit=TRUE,unconditional=TRUE)
uc_ref <- predict(sm,type="response",se.fit=TRUE,unconditional=TRUE)
near(uc$fit,uc_ref$fit[,5]); near(uc$se.fit,uc_ref$se.fit[,5])
# The callback must use the covariance supplied by predict.gam(), including Vc.
X <- predict(sm,type="lpmatrix")
eta <- drop(X %*% coef(sm)); h <- 1e-5
slopes <- (response_values(eta+h,sm$family$getTheta(TRUE))-
             response_values(eta-h,sm$family$getTheta(TRUE)))/(2*h)
stopifnot(!is.null(sm$Vc))
near(uc_ref$se.fit,abs(slopes)*sqrt(rowSums((X %*% sm$Vc)*X)),1e-7)
excluded <- predict(sm,type="response",exclude='s(x)',se.fit=TRUE)
near(excluded$fit,response_values(rep(coef(sm)[1],nrow(smooth_data)),sm$family$getTheta(TRUE)))
near(predict_ordbeta(sm,exclude='s(x)'),excluded$fit[,5])
near(predict_ordbeta(sm,exclude='s(x)',se.fit=TRUE)$se.fit,excluded$se.fit[,5])
nd <- d[1:12,c("x","off")]
prediction_contract(fit,nd)
prediction_contract(fit,nd[1,,drop=FALSE])
nd_zero <- nd; nd_zero$off <- 0
prediction_contract(fit,nd_zero)
near(predict(fit,type="response"),response_values(fit$linear.predictors,fit$family$getTheta(TRUE)))
lp <- predict(fit,nd,se.fit=TRUE)
lp <- lapply(lp, as.numeric)
p <- predict_ordbeta(fit,nd,type="probabilities",se.fit=TRUE)
mu <- predict_ordbeta(fit,nd,type="conditional",se.fit=TRUE)
pr <- predict_ordbeta(fit,nd,se.fit=TRUE)
stopifnot(is.null(dim(mu$fit)),is.null(dim(pr$fit)),identical(dim(p$fit),c(nrow(nd),3L)))
near(rowSums(p$fit),rep(1,nrow(nd)))
near(pr$fit,p$fit[,3]+p$fit[,2]*mu$fit)
near(pr$fit,predict(fit,nd,type="response")[,5])
near(pr$se.fit,predict(fit,nd,type="response",se.fit=TRUE)$se.fit[,5])
near(mu$fit,plogis(lp$fit))
near(mu$se.fit,plogis(lp$fit)*plogis(-lp$fit)*lp$se.fit)
h <- 1e-5
pp <- components(lp$fit+h,fit$family$getTheta(TRUE)); pm <- components(lp$fit-h,fit$family$getTheta(TRUE))
near(p$se.fit,abs((pp-pm)/(2*h))*lp$se.fit,1e-7)
near(pr$se.fit,abs((pp[,3]+pp[,2]*plogis(lp$fit+h)-pm[,3]-pm[,2]*plogis(lp$fit-h))/(2*h))*lp$se.fit,1e-7)
near(residuals(fit,type="response"),d$y-predict(fit,type="response")[,5])
near(sum(residuals(fit,type="deviance")^2),deviance(fit),1e-6)
near(fitted(fit),plogis(predict(fit)))
nd$x[2] <- NA
stopifnot(is.na(predict_ordbeta(fit,nd)[2]))
for (size in c(1L,5L,1000L)) {
  missing_pred <- predict(fit,nd,type="response",se.fit=TRUE,block.size=size)
  stopifnot(all(is.na(missing_pred$fit[2,])),all(is.na(missing_pred$se.fit[2,])),
            identical(rownames(missing_pred$fit),rownames(nd)))
  keep <- seq_len(nrow(nd)) != 2L
  omitted <- predict(fit,nd,type="response",se.fit=TRUE,na.action=na.omit,block.size=size)
  near(omitted$fit,missing_pred$fit[keep,]); near(omitted$se.fit,missing_pred$se.fit[keep,])
  for (type in c("response","conditional","probabilities")) {
    cols <- switch(type,response=5L,conditional=1L,probabilities=2:4)
    helper <- predict_ordbeta(fit,nd,type=type,se.fit=TRUE,block.size=size)
    stopifnot(identical(unname(is.na(helper$fit)),unname(is.na(missing_pred$fit[,cols]))),
              identical(unname(is.na(helper$se.fit)),unname(is.na(missing_pred$se.fit[,cols]))))
    near(helper$fit[!is.na(helper$fit)],missing_pred$fit[,cols][!is.na(helper$fit)])
    near(helper$se.fit[!is.na(helper$se.fit)],missing_pred$se.fit[,cols][!is.na(helper$se.fit)])
    helper_omit <- predict_ordbeta(fit,nd,type=type,se.fit=TRUE,na.action=na.omit,block.size=size)
    near(helper_omit$fit,omitted$fit[,cols]); near(helper_omit$se.fit,omitted$se.fit[,cols])
  }
}
missing <- d; missing$y[1] <- NA
gm <- gam(y~x+offset(off),data=missing,family=ordbeta(),method="ML",na.action=na.exclude)
stopifnot(length(predict_ordbeta(gm))==nrow(d),is.na(predict_ordbeta(gm)[1]))
for (action in list(na.exclude,na.omit)) {
  gm <- gam(y~x+offset(off),data=missing,family=ordbeta(),method="ML",na.action=action)
  pred <- predict(gm,type="response",se.fit=TRUE)
  for (type in c("response","conditional","probabilities")) {
    cols <- switch(type,response=5L,conditional=1L,probabilities=2:4)
    helper <- predict_ordbeta(gm,type=type,se.fit=TRUE)
    stopifnot(identical(unname(is.na(helper$fit)),unname(is.na(pred$fit[,cols]))),
              identical(unname(is.na(helper$se.fit)),unname(is.na(pred$se.fit[,cols]))))
    near(helper$fit[!is.na(helper$fit)],pred$fit[,cols][!is.na(pred$fit[,cols])])
    near(helper$se.fit[!is.na(helper$se.fit)],pred$se.fit[,cols][!is.na(pred$se.fit[,cols])])
  }
}
expect_error(predict_ordbeta(lm(1:5~1)),"gam fitted")

# Distribution callbacks use the interior beta mean, not the response mean.
for (new in c(FALSE,TRUE)) {
  mu <- if (new) plogis(predict(fit, d[1:15, ], type = "link")) else fitted(fit)
  set.seed(19)
  draws <- replicate(3, fit$family$rd(mu))
  quant <- sapply(c(.1, .5, .9), function(p) fit$family$qf(p, mu))
  stopifnot(all(is.finite(draws)),all(draws>=0 & draws<=1),ncol(draws)==3,
            all(is.finite(quant)),all(quant>=0 & quant<=1),ncol(quant)==3)
}
