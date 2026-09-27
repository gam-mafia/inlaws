library(inlaws)

near <- function(x,y,tol=1e-7) {
  stopifnot(isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=tol)))
}
expect_error <- function(expr,text) {
  e <- tryCatch(expr,error=identity)
  stopifnot(inherits(e,"error"),grepl(text,conditionMessage(e),fixed=TRUE))
}
dderiv <- getFromNamespace(".dirmult_deriv","inlaws")
prob <- getFromNamespace(".dirmult_prob","inlaws")
rising <- getFromNamespace(".dirmult_rising","inlaws")

## Independent density on moderate inputs (including zeros, totals of one,
## non-unit weights and the >64-count branch of the rising-factorial helper).
reference <- function(y,eta) {
  p <- exp(cbind(0,eta[,-ncol(eta),drop=FALSE]))
  p <- p/rowSums(p); phi <- exp(eta[,ncol(eta)])
  alpha <- p*phi; N <- rowSums(y)
  lgamma(N+1)-rowSums(lgamma(y+1))+lgamma(phi)-lgamma(N+phi)+
    rowSums(lgamma(alpha+y)-lgamma(alpha))
}
set.seed(812)
y <- matrix(sample(0:100,36,replace=TRUE),12,3)
y[1,] <- c(0,1,0)
eta <- matrix(rnorm(36),12,3)
w <- seq(0,2,length.out=12)
z <- dderiv(y,eta,w)
near(z$l0,reference(y,eta),1e-10)
near(z$l,sum(w*reference(y,eta)),1e-10)

## Central differences of an independent likelihood and of analytic scores.
i2 <- dirmult(2)$tri$i2
for(j in 1:3) {
  ep <- em <- eta; ep[,j] <- ep[,j]+1e-5; em[,j] <- em[,j]-1e-5
  near(z$l1[,j],w*(reference(y,ep)-reference(y,em))/2e-5,2e-6)
  a <- dderiv(y,ep,w); b <- dderiv(y,em,w)
  for(k in 1:3) near(z$l2[,i2[j,k]],(a$l1[,k]-b$l1[,k])/2e-5,2e-6)
}
## Check the two-category beta-binomial identity.
y2 <- cbind(0:20,20:0)
e2 <- cbind(seq(-2,2,length.out=21),rep(log(3),21))
p2 <- plogis(e2[,1]); a <- 3*(1-p2); b <- 3*p2
near(dderiv(y2,e2)$l0,lchoose(20,y2[,1])+
       lbeta(y2[,1]+a,y2[,2]+b)-lbeta(a,b),1e-10)
## Totals of one have exactly zero concentration score and Hessian.
one <- diag(3); eo <- cbind(c(-2,1,4),c(0,3,-1),c(-10,0,30))
z1 <- dderiv(one,eo)
near(z1$l0,log(diag(prob(eo)$p)),1e-10)
near(z1$l1[,3],rep(0,3),1e-12)
near(z1$l2[,i2[3,3]],rep(0,3),1e-12)

## Stable recurrence reference covers both branches of the asymptotic switch,
## tiny concentrations, large counts and effectively multinomial data.
for(m in c(0,1,64,65,1000)) for(la in c(-700,-10,0,10,log(max(1,m)*1e3)+c(-.01,.01),50)) {
  k <- if(m>1) seq_len(m-1) else numeric(0)
  zz <- log(k)-la
  rr <- sum(pmax(zz,0)+log1p(exp(-abs(zz))))
  q <- c(if(m>0) 1 else numeric(0),plogis(-zz))
  h <- rising(la,m)
  near(h$r,rr,2e-7); near(h$A,sum(q),2e-7); near(h$B,-sum(q*q),2e-7)
}
em <- eta; em[,3] <- 50
pm <- prob(em)$p
near(dderiv(y,em)$l0,lgamma(rowSums(y)+1)-rowSums(lgamma(y+1))+rowSums(y*log(pm)),1e-10)
extreme <- matrix(c(-700,700,-20,700,-700,20),2,3,byrow=TRUE)
stopifnot(all(is.finite(unlist(dderiv(matrix(c(0,3,7,4,0,8),2,3,byrow=TRUE),extreme)))))

## Simulated count compositions, with identifiable covariate-dependent phi.
set.seed(731)
n <- 500L
d <- data.frame(x=runif(n),z=runif(n),o=rnorm(n,sd=.2),w=sample(1:3,n,TRUE))
truth <- cbind(.5*sin(2*pi*d$x),.7*d$x-.2,1+1.5*d$z)
p <- prob(truth)$p
N <- sample(20:40,n,TRUE)
d$y <- t(vapply(seq_len(n),function(i) {
  as.numeric(rmultinom(1,N[i],rgamma(3,shape=exp(truth[i,3])*p[i,])))
},numeric(3)))
colnames(d$y) <- c("a","b","c")

## Intercept-only fit against independent, unpenalized likelihood optimization.
b0 <- gam(list(y~1,~1,~1),family=dirmult(2),data=d,weights=w,method="REML")
opt <- optim(c(0,0,1),function(b) -sum(d$w*reference(d$y,
             matrix(b,n,3,byrow=TRUE))),method="L-BFGS-B",lower=c(-8,-8,-5),
             upper=c(8,8,8),control=list(factr=1e4,pgtol=1e-6))
near(coef(b0),opt$par,2e-5)
near(as.numeric(logLik(b0)),-opt$value,1e-9)
near(AIC(b0),-2*as.numeric(logLik(b0))+2*sum(b0$edf),1e-10)
near(b0$null.deviance,b0$deviance,1e-7)
## Integer weights agree with replicated rows.
dr <- d[rep(seq_len(n),d$w),]
br <- gam(list(y~1,~1,~1),family=dirmult(2),data=dr,method="REML")
near(coef(br),coef(b0),1e-6)

b <- gam(list(y~s(x,k=6)+offset(o),~s(x,k=6)+offset(.3*o),~z+offset(.2*o)),
         family=dirmult(2),data=d,weights=w,method="REML",optimizer="efs")
stopifnot(all(is.finite(coef(b))),all(is.finite(vcov(b))),
          b$outer.info$conv=="full convergence")
e <- predict(b,type="link")
near(b$l,sum(d$w*reference(d$y,e)),1e-9)
near(b$linear.predictors,e,1e-10)
stopifnot(abs(coef(b)["z.2"]-1.5)<.5)
pr <- predict(b,type="response",se.fit=TRUE)
stopifnot(identical(dim(pr$fit),c(n,4L)),all(is.finite(pr$se.fit)),all(pr$se.fit>=0))
near(rowSums(pr$fit[,1:3]),rep(1,n),1e-12)
near(pr$fit[,4],exp(e[,3]),1e-12)
nd <- d[1:7,c("x","z","o")]
near(predict(b,newdata=nd,type="response"),pr$fit[1:7,],1e-10)
near(predict(b,newdata=d[,c("x","z","o")],type="response",block.size=37),pr$fit,1e-10)
## Numerical coefficient Jacobian checks response SEs, including covariance.
X <- predict(b,newdata=nd,type="lpmatrix")
f <- b$family
base <- f$predict(f,X=X,beta=coef(b),off=list(nd$o,.3*nd$o,.2*nd$o),se=TRUE,Vb=vcov(b))
for(k in 1:4) {
  J <- matrix(0,nrow(nd),length(coef(b)))
  for(j in seq_along(coef(b))) {
    bp <- bm <- coef(b); bp[j] <- bp[j]+1e-5; bm[j] <- bm[j]-1e-5
    J[,j] <- (f$predict(f,X=X,beta=bp,off=list(nd$o,.3*nd$o,.2*nd$o))$fit[,k]-
              f$predict(f,X=X,beta=bm,off=list(nd$o,.3*nd$o,.2*nd$o))$fit[,k])/2e-5
  }
  near(base$se.fit[,k],sqrt(rowSums((J%*%vcov(b))*J)),1e-7)
}
near(residuals(b,"response"),d$y-N*pr$fit[,1:3],1e-10)
near(sum(residuals(b)^2),b$deviance,1e-10)
stopifnot(all(is.finite(residuals(b,"pearson"))))

## Fixed smoothing parameters, predictor-shared coefficients, zero weights,
## and na.exclude (an entire count row is dropped).
Gfix <- gam(list(y~s(x,k=5),~s(x,k=5),~1),family=dirmult(2),data=d,
            method="REML",fit=FALSE)
bfix <- gam(G=Gfix,method="REML",sp=c(.5,2))
near(bfix$full.sp,c(.5,2),1e-12)
bs <- gam(list(y~1,~1,~1,1+2~s(x,k=5)-1),family=dirmult(2),data=d,method="REML")
stopifnot(all(is.finite(predict(bs,type="response",se.fit=TRUE)$se.fit)))
dm <- d; dm$y[4,2] <- NA; dm$x[8] <- NA
bm <- gam(list(y~x,~x,~z),family=dirmult(2),data=dm,method="REML",na.action=na.exclude)
stopifnot(nrow(bm$y)==n-2L,nrow(predict(bm,type="response"))==n,
          all(is.na(predict(bm,type="response")[c(4,8),])))
dw <- d; dw$w <- rep(1,n); dw$w[1:10] <- 0
bz <- gam(list(y~1,~1,~1),family=dirmult(2),data=dw,weights=w,method="REML")
bsub <- gam(list(y~1,~1,~1),family=dirmult(2),data=d[-(1:10),],method="REML")
near(coef(bz),coef(bsub),1e-6)

## Validation and explicit unsupported modes.
validate <- dirmult(2)$validate
expect_error(dirmult(0),"positive integer")
expect_error(dirmult(1.5),"positive integer")
expect_error(validate(y[,1:2],rep(1,12)),"K+1")
ybad <- y; ybad[1,1] <- .1
expect_error(validate(ybad,rep(1,12)),"integers")
ybad[1,1] <- -1
expect_error(validate(ybad,rep(1,12)),"integers")
ybad[1,1] <- Inf
expect_error(validate(ybad,rep(1,12)),"finite")
ybad <- y; ybad[1,] <- 0
expect_error(validate(ybad,rep(1,12)),"totals")
ybad <- y; ybad[,1] <- 0
expect_error(validate(ybad,rep(1,12)),"categories")
expect_error(validate(y,rep(-1,12)),"weights")
expect_error(validate(y,rep(0,12)),"weights")
expect_error(validate(diag(3),rep(1,3)),"not identifiable")
expect_error(gam(G=Gfix,sp=c(.5,-1),method="REML"),"all smoothing parameters")
expect_error(bam(list(y~x,~x,~z),family=dirmult(2),data=d),"general families not supported")

## A smooth concentration predictor, and a category permutation with identical
## intercept models, protect the parameterization and predictor ordering.
bphi <- gam(list(y~s(x,k=5),~s(x,k=5),~s(z,k=5)),
            family=dirmult(2),data=d,method="REML")
stopifnot(all(is.finite(coef(bphi))),cor(predict(bphi)[,3],truth[,3])>.8)
dp <- d; dp$y <- dp$y[,c(3,1,2)]
bp <- gam(list(y~1,~1,~1),family=dirmult(2),data=dp,weights=w,method="REML")
near(predict(bp,type="response")[,c(2,3,1,4)],predict(b0,type="response"),2e-6)

## More than three categories: verify packed cross-Hessian indexing.
set.seed(9)
y6 <- matrix(sample(0:15,48,TRUE),8,6)
e6 <- matrix(rnorm(48),8,6)
z6 <- dderiv(y6,e6); i6 <- dirmult(5)$tri$i2
for(j in 1:6) {
  ep <- em <- e6; ep[,j] <- ep[,j]+1e-5; em[,j] <- em[,j]-1e-5
  for(k in 1:6) near(z6$l2[,i6[j,k]],
                    (dderiv(y6,ep)$l1[,k]-dderiv(y6,em)$l1[,k])/2e-5,2e-6)
}
