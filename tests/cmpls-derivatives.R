library(inlaws)
near <- function(x, y, tol = 3e-5) stopifnot(length(x) == length(y),
  all(is.finite(x)), all(is.finite(y)), all(abs(x-y) <= tol * (1+abs(y))))
f <- cmpls(control = list(sum_tol = 1e-12, mean_tol = 1e-12))
y <- c(0, 1, 4, 12); w <- c(.5, 2, 1, 0)
X <- cbind(1, c(-1, .2, .6, 1), 1, c(.5, -.3, 1, 2))
attr(X, 'lpi') <- list(1:2, 3:4)
b <- c(.8, .3, -.2, .4)
eta <- cbind(drop(X[,1:2] %*% b[1:2]), drop(X[,3:4] %*% b[3:4]))
eval <- function(e, deriv=4) f$ll(y, X, b, w, f, eta=e, deriv=deriv, ncv=TRUE)
z <- eval(eta); h <- 1e-5
# Differentiate each packed derivative through order four in both predictors.
for (j in 1:2) {
  ep <- em <- eta; ep[,j] <- ep[,j]+h; em[,j] <- em[,j]-h
  pp <- eval(ep); mm <- eval(em)
  near(z$l1[,j], (pp$l0-mm$l0)/(2*h))
  for (k in 1:3) for (q in 0:k) {
    col <- q+1L+(j==2L)
    near(z[[paste0('l',k+1)]][,col],
      (pp[[paste0('l',k)]][,q+1]-mm[[paste0('l',k)]][,q+1])/(2*h))
  }
}
ll <- function(b) f$ll(y,X,b,w,f,deriv=1)
z <- ll(b)
for (j in 1:4) {
  bp <- bm <- b; bp[j] <- bp[j]+h; bm[j] <- bm[j]-h
  pp <- ll(bp); mm <- ll(bm)
  near(z$lb[j], (pp$l-mm$l)/(2*h))
  near(z$lbb[,j], (pp$lb-mm$lb)/(2*h))
}
# Independent PMF summation with independent mean inversion.
for (nu in c(.3,1,3)) {
  mu <- c(.2, 1, 4, 10)
  e <- cbind(log(mu),log(nu)); obs <- eval(e)
  for (i in 1:3) {
    k <- 0:500
    logz <- function(a) {v <- k*a-nu*lgamma(k+1); max(v)+log(sum(exp(v-max(v))))}
    meanfun <- function(a) {v <- k*a-nu*lgamma(k+1); p <- exp(v-max(v)); sum(k*p)/sum(p)}
    a <- uniroot(function(a) meanfun(a)-mu[i], c(-20,20),tol=1e-12)$root
    near(obs$l0[i],w[i]*(y[i]*a-nu*lgamma(y[i]+1)-logz(a)))
  }
}
# Zero weights mask even overflowing predictors.
e <- eta; e[4,] <- Inf
zz <- eval(e)
near(zz$l0,z <- eval(eta)$l0)
stopifnot(all(zz$l4[4,] == 0))
