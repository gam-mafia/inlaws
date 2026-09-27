library(inlaws)
near <- function(x,y,tol=2e-5) stopifnot(length(x)==length(y), all(is.finite(x)),
  all(is.finite(y)),all(abs(x-y)<=tol*(1+abs(y))))
error <- function(expr,pattern) {
  e <- tryCatch(force(expr),error=identity)
  stopifnot(inherits(e,'error'),grepl(pattern,conditionMessage(e)))
}
f <- cmpls()
m <- c(.01,.3,2,10)
near(f$cdf(2,cbind(m,1)),ppois(2,m))
near(f$qf(.8,cbind(m,1)),qpois(.8,m))
near(f$cdf(2,cbind(m,1),logp=TRUE),ppois(2,m,log.p=TRUE))
stopifnot(f$qf(1,cbind(0,2))==0,is.infinite(f$qf(1,cbind(2,2))),
  f$qf(0,cbind(2,2))==0,is.na(f$cdf(NA,cbind(2,2))),
  is.na(f$cdf(2,cbind(NA,2))),length(f$rd(matrix(numeric(),0,2)))==0)
for(nu in c(.3,1,4)) {
  near(f$cdf(3,cbind(m,nu)),cmp(nu)$cdf(3,m))
  near(f$qf(.7,cbind(m,nu)),cmp(nu)$qf(.7,m))
}
set.seed(937)
n <- 100; d <- data.frame(x=runif(n),z=runif(n),a=runif(n,-.1,.1),b=runif(n,-.1,.1))
d$y <- f$rd(cbind(exp(1+.5*d$x+d$a),exp(.2+.5*d$z+d$b)))
d$w <- rep(c(0,.5,1,2),length.out=n)
# Unpenalized fit agrees with direct likelihood optimization.
a <- gam(list(y~x+offset(a),~z+offset(b)),data=d,weights=w,family=f,method='REML')
X <- predict(a,type='lpmatrix')
obj <- function(beta) -f$ll(d$y,X,beta,d$w,f,offset=list(d$a,d$b))$l
ref <- optim(coef(a),obj,method='BFGS',control=list(reltol=1e-10))
near(coef(a),ref$par,1e-4)
near(as.numeric(logLik(a)),-ref$value)
near(predict(a,type='response'),fitted(a))
near(predict(a,type='response'),exp(predict(a,type='link')))
stopifnot(identical(colnames(fitted(a)),c('mu','nu')),
 all(is.finite(predict(a,se.fit=TRUE,type='response')$se.fit)),
 all(is.finite(residuals(a))),all(is.finite(residuals(a,type='pearson'))),
 is.finite(a$null.deviance),is.finite(AIC(a)))
error(cmpls(link='log'),'two')
error(cmpls(control=list(typo=1)),'control')
error(f$qf(1.1,cbind(1,1)),'probabilities')
error(f$cdf(1,cbind(1,0)),'positive')
error(gam(list(y~1,~1),data=data.frame(y=rep(0,10)),family=f,method='REML'),'all-zero')
error(gam(list(y~1,~1),data=data.frame(y=c(0,.5,2,3)),family=f,method='REML'),'integer')

# REML outer Newton and BFGS agree for a smooth mean predictor.
set.seed(20)
d <- data.frame(x=runif(70))
d$y <- cmpls()$rd(cbind(exp(1+d$x),exp(.4+d$x)))
form <- list(y~s(x,k=4),~x)
a <- gam(form,data=d,family=cmpls(),method="REML")
b <- gam(form,data=d,family=cmpls(),method="REML",optimizer=c("outer","bfgs"))
near(coef(a),coef(b),.002)
# A nearly linear smooth can put the optimum at infinite smoothing;
# compare the fitted model and criterion rather than diverging sp values.
near(fitted(a),fitted(b),.002)
near(a$gcv.ubre,b$gcv.ubre,1e-5)
