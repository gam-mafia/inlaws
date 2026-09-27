library(inlaws)
near <- function(a, b, tol = 1e-6) {
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
}
fails <- function(expr) stopifnot(inherits(tryCatch(expr, error = identity), "error"))
f <- clognormal(theta = 0.7)
y <- cbind(c(0.5, 1.1, 2, 0.7, 1.4, exp(8), exp(-9)),
           c(0.5, -Inf, Inf, 1.3, 1.40001, Inf, -Inf))
mu <- c(1.2, 1.4, 0.8, 1.1, 0.9, 1, 1)
w <- c(1, 2, 0.5, 3, 1, 0.7, 1)
th <- log(0.7)
z <- f$prepare.response(y)
ll <- c(dlnorm(y[1,1], log(mu[1])-.7^2/2, .7, log=TRUE),
        plnorm(y[2,1], log(mu[2])-.7^2/2, .7, log.p=TRUE),
        plnorm(y[3,1], log(mu[3])-.7^2/2, .7, lower.tail=FALSE, log.p=TRUE),
        log(plnorm(y[4,2], log(mu[4])-.7^2/2, .7)-plnorm(y[4,1], log(mu[4])-.7^2/2, .7)),
        log(plnorm(y[5,2], log(mu[5])-.7^2/2, .7)-plnorm(y[5,1], log(mu[5])-.7^2/2, .7)),
        plnorm(y[6,1], -.7^2/2, .7, lower.tail=FALSE, log.p=TRUE),
        plnorm(y[7,1], -.7^2/2, .7, log.p=TRUE))
near(f$aic(z, mu, wt=w), -2*sum(w*ll))
near(sum(f$dev.resids(z, mu, w))/2 - f$ls(z,w,th,1)$ls, -sum(w*ll))
# All first through fourth derivative slots, independently differenced along
# both axes. These include mixed censoring, narrow intervals and tail events.
d <- f$Dd(z,mu,th,w,2)
h <- 2e-5
mp <- f$Dd(z,mu+h,th,w,2); mm <- f$Dd(z,mu-h,th,w,2)
tp <- f$Dd(z,mu,th+h,w,2); tm <- f$Dd(z,mu,th-h,w,2)
near(d$Dmu, (f$dev.resids(z,mu+h,w)-f$dev.resids(z,mu-h,w))/(2*h), 1e-5)
near(d$Dth, (f$dev.resids(z,mu,w,th+h)-f$dev.resids(z,mu,w,th-h))/(2*h), 1e-5)
for (pair in list(c("Dmu2","Dmu"),c("Dmu3","Dmu2"),c("Dmu4","Dmu3"),
                  c("Dmuth","Dth"),c("Dmu2th","Dmuth"),c("Dmu3th","Dmu2th"),
                  c("Dmuth2","Dth2"),c("Dmu2th2","Dmuth2"))) {
  near(d[[pair[1]]], (mp[[pair[2]]]-mm[[pair[2]]])/(2*h), 2e-4)
}
for (pair in list(c("Dth2","Dth"),c("Dmuth","Dmu"),c("Dmu2th","Dmu2"),
                  c("Dmu3th","Dmu3"),c("Dmuth2","Dmuth"),c("Dmu2th2","Dmu2th"))) {
  near(d[[pair[1]]], (tp[[pair[2]]]-tm[[pair[2]]])/(2*h), 2e-4)
}
s0 <- f$ls(z,w,th,1); sp <- f$ls(z,w,th+h,1); sm <- f$ls(z,w,th-h,1)
near(s0$lsth1,(sp$ls-sm$ls)/(2*h)); near(s0$lsth2,(sp$lsth1-sm$lsth1)/(2*h))
# No shared theta state; serialization preserves fitted family state.
g <- clognormal(); g$putTheta(log(2)); near(f$getTheta(TRUE),.7)
f2 <- unserialize(serialize(f,NULL)); near(f2$getTheta(TRUE),.7)
# Validation and row-safe censoring attributes.
for (bad in list(0,-1,Inf,cbind(1,0),cbind(1,.5),matrix(1,2,3))) fails(f$prepare.response(bad))
for (bad in list(0,NA,Inf,c(1,2))) fails(clognormal(theta=bad))
fails(clognormal(link="identity"))
near(attr(f$subsety(z,c(4,2)),"censor"),y[c(4,2),2])
stopifnot(is.matrix(f$subsety(y,1)))
near(f$qf(.5,mu,w,1),exp(log(mu)-.7^2/2))
near(f$cdf(f$qf(.3,mu,w,1),mu,w,1),rep(.3,length(mu)))
near(f$variance(mu),mu^2*expm1(.7^2))
set.seed(921)
sim <- f$rd(rep(2,50000),1,1); stopifnot(abs(mean(sim)-2)<.03)

set.seed(92)
n <- 220
dat <- data.frame(x=runif(n), off=runif(n,-.1,.1))
raw <- rlnorm(n, .2+sin(6*dat$x)+dat$off, .5)
dat$y <- raw
# Uncensored fits agree with the usual Gaussian log-response fit.
b <- gam(y~s(x,k=7)+offset(off),data=dat,family=clognormal(),method="REML")
ref <- gam(log(y)~s(x,k=7)+offset(off),data=dat,family=gaussian(),method="REML")
near(b$family$getTheta(TRUE),sqrt(ref$sig2),1e-4)
near(predict(b),predict(ref)+ref$sig2/2,1e-4)
near(fitted(b),predict(b,type="response"))
near(predict(b,type="response",se.fit=TRUE)$se.fit,
     exp(predict(b))*predict(b,se.fit=TRUE)$se.fit)
near(as.numeric(logLik(b)),sum(dlnorm(raw,log(fitted(b))-b$family$getTheta(TRUE)^2/2,
                                   b$family$getTheta(TRUE),log=TRUE)),1e-6)
# Exact observations agree with direct likelihood optimization, including sigma.
mle <- gam(y~x,data=dat,family=clognormal(),method="ML")
op <- optim(c(.2,0,log(.5)),function(p) -sum(dlnorm(raw,p[1]+p[2]*dat$x-exp(2*p[3])/2,
                  exp(p[3]),log=TRUE)),method="BFGS",control=list(reltol=1e-12))
near(coef(mle),op$par[1:2],1e-4); near(mle$family$getTheta(),op$par[3],1e-4)
# Mixed censoring with observation-specific limits; a cnorm comparison is a
# test oracle only. Transform finite entries without taking log(-Inf).
yc <- cbind(raw,raw)
id <- seq_len(n) %% 4
il <- id==1 & raw<1; yc[il,] <- cbind(rep(1,sum(il)),rep(-Inf,sum(il)))
ir <- id==2 & raw>2; yc[ir,] <- cbind(rep(2,sum(ir)),rep(Inf,sum(ir)))
ii <- id==3; yc[ii,] <- cbind(raw[ii]*.9,raw[ii]*1.1)
dat$y <- yc
logged <- yc; logged[is.finite(logged)] <- log(logged[is.finite(logged)])
dat$logged <- logged
bc <- gam(y~s(x,k=7)+offset(off),data=dat,family=clognormal(),method="REML")
br <- gam(logged~s(x,k=7)+offset(off),data=dat,family=cnorm(),method="REML")
near(bc$family$getTheta(TRUE),br$family$getTheta(TRUE),2e-3)
near(predict(bc),predict(br)+br$family$getTheta(TRUE)^2/2,2e-3)
stopifnot(all(is.finite(residuals(bc))), all(is.na(residuals(bc,type="response")[yc[,1]!=yc[,2]])))
# Fixed sigma, frequency weights, zero weights, missing rows and subset.
dat$w <- rep(c(0,1,2),length.out=n)
bw <- gam(y~x,data=dat,weights=w,family=clognormal(.5),method="ML")
repdat <- dat[rep(seq_len(n),dat$w),]
brep <- gam(y~x,data=repdat,family=clognormal(.5),method="ML")
near(coef(bw),coef(brep),1e-5); near(as.numeric(logLik(bw)),as.numeric(logLik(brep)),1e-5)
fails(gam(y~x,data=dat,weights=rep(-1,n),family=clognormal()))
dat$x[1] <- NA
bn <- gam(y~x,data=dat,subset=seq_len(n)>5,family=clognormal(.5),method="ML")
stopifnot(nobs(bn)==n-5)
bna <- gam(y~x,data=dat,family=clognormal(.5),method="ML",na.action=na.exclude)
stopifnot(length(residuals(bna))==n,is.na(residuals(bna)[1]))
# Both bam paths must retain censoring data when processing row blocks.
dat$x[1] <- .5
for (discrete in c(FALSE,TRUE)) {
  bb <- bam(y~s(x,k=7)+offset(off),data=dat,family=clognormal(),
            method="fREML",discrete=discrete,chunk.size=55)
  stopifnot(all(is.finite(coef(bb))),is.finite(bb$family$getTheta(TRUE)))
  near(predict(bb,type="response"),fitted(bb),1e-6)
  near(predict(bb),predict(bc),.04)
}
# Extreme and almost coincident finite bounds remain finite even when taking
# logs would round both endpoints to the same floating-point number.
extreme <- cbind(c(exp(10), exp(-10), 1e100),
                 c(exp(10.01), exp(-9.99), 1e100*(1+1e-14)))
ed <- f$Dd(extreme,rep(1,3),th,rep(1,3),2)
stopifnot(all(is.finite(unlist(ed))), all(is.finite(f$dev.resids(extreme,rep(1,3),rep(1,3)))))
# A tiny interval has log mass log(width) + log density to first order.
lo <- extreme[3,1]; hi <- extreme[3,2]
near(f$aic(extreme[3,,drop=FALSE],1,wt=1),
     -2*(log(hi-lo)+dlnorm(lo,-.7^2/2,.7,log=TRUE)),1e-8)
# An estimated sigma also respects frequency weights.
bwe <- gam(y~x,data=dat,weights=w,family=clognormal(-.4),method="ML")
repdat <- dat[rep(seq_len(n),dat$w),]
brepe <- gam(y~x,data=repdat,family=clognormal(),method="ML")
near(coef(bwe),coef(brepe),1e-5)
near(bwe$family$getTheta(),brepe$family$getTheta(),1e-5)
# High-order tail derivatives need stable log-CDF evaluation as well as a
# finite likelihood. Check beyond the range where Mills-ratio subtraction is
# numerically safe, for both one-sided and finite-interval observations.
yt <- cbind(c(exp(100),exp(-100),exp(100)),c(Inf,-Inf,exp(100.1)))
dt <- f$Dd(yt,rep(1,3),th,rep(1,3),2)
pt <- f$Dd(yt,rep(1,3),th+h,rep(1,3),2)
mt <- f$Dd(yt,rep(1,3),th-h,rep(1,3),2)
near(dt$Dmu2th2,(pt$Dmu2th-mt$Dmu2th)/(2*h),1e-6)
near(f$aic(yt[1,,drop=FALSE],1,wt=1),
     -2*plnorm(yt[1,1],-.7^2/2,.7,lower.tail=FALSE,log.p=TRUE))
