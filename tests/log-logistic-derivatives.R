library(inlaws)
near <- function(a, b, tol = 1e-5) {
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
}
evaluate <- getFromNamespace(".llog_eval", "inlaws")
take <- function(r, a, b, field = "ell")
  r[[field]][, which(r$powers[, 1] == a & r$powers[, 2] == b)]*factorial(a)*factorial(b)
# Independent density on the original response scale.
loglik <- function(y, mu, s) dlogis(log(y), log(mu*sinpi(s)/(pi*s)), s, log = TRUE)-log(y)
y <- c(.1, 1, 4, exp(-50), exp(50), .8, 1.1, 2)
mu <- c(.3, 1.5, 3, 2, 1, .7, 1.2, 3)
t <- qlogis(c(.15, .3, .7, .2, .8, .001, .99, .5))
h <- 2e-6
for (predictor in c(FALSE, TRUE)) {
  u <- if (predictor) log(mu) else mu
  r <- evaluate(y, u, t, 4L, predictor)
  near(r$ell[, 1], loglik(y, mu, plogis(t)), 1e-10)
  up <- evaluate(y, u+h, t, 4L, predictor); um <- evaluate(y, u-h, t, 4L, predictor)
  tp <- evaluate(y, u, t+h, 4L, predictor); tm <- evaluate(y, u, t-h, 4L, predictor)
  for (k in 1:4) for (b in 0:k) {
    a <- k-b
    if (a) near(take(r,a,b), (take(up,a-1,b)-take(um,a-1,b))/(2*h), 5e-5)
    if (b) near(take(r,a,b), (take(tp,a,b-1)-take(tm,a,b-1))/(2*h), 5e-5)
  }
}
# Boundary stability, including density tails and very narrow distributions.
for (tt in c(-25, -8, 8, 25)) {
  r <- evaluate(c(.5, 1, 2), rep(0,3), tt, 4L, TRUE)
  stopifnot(all(is.finite(r$ell)), all(is.finite(r$sat)))
  rp <- evaluate(c(.5, 1, 2), rep(0,3), tt+1e-5, 4L, TRUE)
  rm <- evaluate(c(.5, 1, 2), rep(0,3), tt-1e-5, 4L, TRUE)
  for (k in 1:4) near(take(r,0,k),(take(rp,0,k-1)-take(rm,0,k-1))/2e-5,1e-4)
}
# Extended-family fields, saturated scale derivatives and expected curvature.
f <- log_logistic(.3); w <- c(0, .5, 2, 1, 3, 1, 1, 1)
r <- evaluate(y, mu, f$getTheta(), 4L)
d <- f$Dd(y, mu, f$getTheta(), w, 2)
map <- list(Dmu=c(1,0), Dmu2=c(2,0), Dmu3=c(3,0), Dmu4=c(4,0),
            Dth=c(0,1), Dth2=c(0,2), Dmuth=c(1,1), Dmuth2=c(1,2),
            Dmu2th=c(2,1), Dmu2th2=c(2,2), Dmu3th=c(3,1))
for (nm in names(map)) {
  a<-map[[nm]][1]; b<-map[[nm]][2]
  near(d[[nm]], 2*w*(take(r,a,b,"sat")-take(r,a,b)))
}
near(d$EDmu2,2*w/(3*mu^2*.3^2))
near(d$EDmu3,-2*d$EDmu2/mu)
near(d$EDmu2th,-2*(1-.3)*d$EDmu2)
sat <- f$ls(y,w,f$getTheta(),1)
near(sat$ls,sum(w*r$sat[,1])); near(sat$lsth1,sum(w*take(r,0,1,"sat")))
near(sat$lsth2,sum(w*take(r,0,2,"sat")))
# General-family coefficient contractions, offsets and direct-eta evaluation.
f <- log_logisticls()
X <- cbind(1, seq(-1,1,length.out=length(y)), 1)
attr(X,"lpi") <- list(1:2,3L)
beta <- c(.2,-.1,-1)
off <- list(rep(.1,length(y)),seq(-.2,.2,length.out=length(y)))
eta <- cbind(drop(X[,1:2]%*%beta[1:2])+off[[1]], beta[3]+off[[2]])
r <- f$ll(y,X,beta,w,f,offset=off,deriv=1,ncv=TRUE)
near(r$l,sum(w*loglik(y,exp(eta[,1]),plogis(eta[,2]))))
near(r$l,f$ll(y,X,beta,w,f,offset=off,eta=eta)$l)
for (k in seq_along(beta)) {
  bp<-bm<-beta; bp[k]<-bp[k]+h; bm[k]<-bm[k]-h
  rp<-f$ll(y,X,bp,w,f,offset=off,deriv=1); rm<-f$ll(y,X,bm,w,f,offset=off,deriv=1)
  near(r$lb[k],(rp$l-rm$l)/(2*h)); near(r$lbb[,k],(rp$lb-rm$lb)/(2*h))
}
