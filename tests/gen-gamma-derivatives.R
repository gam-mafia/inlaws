library("inlaws")
evaluate <- getFromNamespace(".gen_gamma_eval", "inlaws")
near <- function(a, b, tol = 1e-5) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
# All mixed partials through order four, including exactly Q=0 and both sides
# of the Stirling and power-series transitions. Difference lower derivatives,
# not four nested numerical differences of log densities.
grid <- expand.grid(sigma = c(.2, .7, 2), Q = c(-.2, 0, 1e-7, .1, 1, 3))
grid <- rbind(grid, data.frame(sigma = 1, Q = c(1/sqrt(20) + c(-1,1)*1e-4, .25+c(-1,1)*1e-4)))
eta <- cbind(log(2), log(grid$sigma), log1p(grid$sigma*grid$Q))
y <- rep(c(.3, 1.3, 4), length.out = nrow(eta))
r <- evaluate(y, eta, 4L)
stopifnot(all(is.finite(r$ell)))
for (k in 1:3) {
  ep <- em <- eta; ep[, k] <- ep[, k]+1e-5; em[, k] <- em[, k]-1e-5
  rp <- evaluate(y, ep, 4L); rm <- evaluate(y, em, 4L)
  for (i in which(rowSums(r$powers)>0 & r$powers[, k]>0)) {
    p <- r$powers[i, ]; p[k] <- p[k]-1L
    ix <- which(rowSums(abs(sweep(r$powers, 2, p))) == 0)
    near(r$ell[, i]*r$powers[i, k], (rp$ell[, ix]-rm$ell[, ix])/2e-5, 3e-4)
  }
}
# Independent gamma density plus change-of-variables Jacobian away from Q=0.
for (i in which(abs(grid$Q) > 1e-6)) {
  s <- grid$sigma[i]; Q <- grid$Q[i]; k <- Q^-2
  a <- log(2)-lgamma(k+s/Q)+lgamma(k)+(s/Q)*log(k)
  logg <- log(k)+Q*(log(y[i])-a)/s
  ref <- dgamma(exp(logg), k, log = TRUE)+logg+log(abs(Q))-log(s)-log(y[i])
  near(r$ell[i, 1], ref, 1e-10)
}
ix <- which(grid$Q == 0)
near(r$ell[ix, 1], dlnorm(y[ix], log(2)-grid$sigma[ix]^2/2, grid$sigma[ix], log=TRUE), 1e-12)
# Shape information must survive exactly at the log-normal limit.
shape <- which(rowSums(abs(sweep(r$powers, 2, c(0,0,1)))) == 0)
stopifnot(any(abs(r$ell[ix, shape]) > .1))
# Packed coefficient derivatives and sandwich against independent differences.
f <- gen_gammalss(); set.seed(52); n <- 12L
X <- cbind(1, runif(n), 1, runif(n), 1, runif(n)); attr(X,"lpi") <- list(1:2,3:4,5:6)
beta <- c(.3,.2,-.7,.1,.15,.1); y <- exp(rnorm(n)); w <- rep(c(0,1,2),4)
l <- f$ll(y,X,beta,w,f,deriv=1,ncv=TRUE)
score <- X
for(k in 1:3) score[,attr(X,"lpi")[[k]]] <- X[,attr(X,"lpi")[[k]]]*l$l1[,k]
near(l$lb, colSums(score)); near(f$sandwich(y,X,beta,w,f),crossprod(score),1e-12)
for(k in seq_along(beta)) {
  bp <- bm <- beta; bp[k] <- bp[k]+1e-5; bm[k] <- bm[k]-1e-5
  lp <- f$ll(y,X,bp,w,f,deriv=1); lm <- f$ll(y,X,bm,w,f,deriv=1)
  near(l$lb[k],(lp$l-lm$l)/2e-5)
  near(l$lbb[,k],(lp$lb-lm$lb)/2e-5)
}
# Unused extreme rows must contribute exactly zero at every derivative order.
y[1] <- .Machine$double.xmax
et <- matrix(0,n,3); et[1,] <- 1000
z <- f$ll(y,X,beta,w,f,eta=et,deriv=4,ncv=TRUE)
stopifnot(z$l0[1]==0,all(z$l1[1,]==0),all(z$l4[1,]==0))
