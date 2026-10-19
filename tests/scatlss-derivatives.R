library("inlaws")
near <- function(a, b, tol = 1e-6) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
evaluate <- getFromNamespace(".scatlss_eval", "inlaws")
# Density references include the scale Jacobian and both sides of each series
# transition. Finite nu corrections persist even at enormous degrees of freedom.
grid <- expand.grid(sigma = c(.03, .7, 3), nu = c(3.01, 6, 40-1e-4, 40+1e-4, 1e3, 1e8, 1e100), z = c(0, .5-1e-5, .5+1e-5, 3, 30))
transition <- expand.grid(sigma=c(.03,3),nu=c(3.01,8,40-1e-4,40+1e-4,1e4),
                          edge=c(-1e-5,1e-5))
transition$z <- sqrt(transition$nu)*(.5+transition$edge)
grid <- rbind(grid,transition[,c("sigma","nu","z")])
eta <- cbind(.2, log(grid$sigma), log(grid$nu-3))
y <- .2 + grid$sigma*grid$z
r <- evaluate(y, eta, 4L)
near(r$ell[, 1], dt(grid$z, grid$nu, log = TRUE)-log(grid$sigma), 1e-11)
stopifnot(all(is.finite(r$ell)))
shape <- which(rowSums(abs(sweep(r$powers,2,c(0,0,1))))==0)
stopifnot(all(r$ell[grid$nu==1e100 & grid$z==0,shape] != 0))
# Difference every lower analytic partial to check all mixed fourth derivatives.
for (k in 1:3) {
  h <- rep(1e-5, nrow(eta)); if (k == 1L) h <- h*grid$sigma
  ep <- em <- eta; ep[, k] <- ep[, k]+h; em[, k] <- em[, k]-h
  rp <- evaluate(y, ep, 4L); rm <- evaluate(y, em, 4L)
  for (i in which(rowSums(r$powers)>0 & r$powers[, k]>0)) {
    p <- r$powers[i, ]; p[k] <- p[k]-1L
    ix <- which(rowSums(abs(sweep(r$powers, 2, p))) == 0)
    multiplier <- grid$sigma^r$powers[i, 1]
    near(multiplier*r$ell[, i]*r$powers[i, k], multiplier*(rp$ell[, ix]-rm$ell[, ix])/(2*h), 3e-4)
  }
}
# Independent tensor-product finite differences of stats::dt(), without
# calling any analytic derivative. Richardson extrapolation reduces bias.
fd <- function(y, eta, powers, h) {
  stencils <- lapply(powers, function(d) {
    if (d == 0) return(list(nodes=0, weights=1))
    nodes <- seq.int(-max(2L,d), max(2L,d))
    a <- outer(0:(length(nodes)-1L), nodes, function(k,x) x^k)
    target <- numeric(length(nodes)); target[d+1L] <- factorial(d)
    list(nodes=nodes, weights=solve(a,target))
  })
  ix <- expand.grid(lapply(stencils, function(s) seq_along(s$nodes)))
  ee <- matrix(eta,nrow(ix),3L,byrow=TRUE); w <- rep(1,nrow(ix))
  step <- c(exp(eta[2]),1,1)*h
  for(k in 1:3) {
    ee[,k] <- ee[,k]+step[k]*stencils[[k]]$nodes[ix[,k]]
    w <- w*stencils[[k]]$weights[ix[,k]]/step[k]^powers[k]
  }
  sum(w*(dt((y-ee[,1])/exp(ee[,2]),3+exp(ee[,3]),log=TRUE)-ee[,2]))
}
for (nu in c(3.1, 8, 39.99, 40.01, 1e4)) {
  e <- c(.2,log(.4),log(nu-3)); y0 <- .2+.4*2.3
  rr <- evaluate(y0,matrix(e,1),4L)
  for(i in 2:nrow(rr$powers)) {
    p <- rr$powers[i,]; a <- fd(y0,e,p,.04); b <- fd(y0,e,p,.02)
    ref <- (4*b-a)/3
    near(.4^p[1]*rr$ell[1,i]*prod(factorial(p)), .4^p[1]*ref, 8e-4)
  }
}
# Coefficient gradient/Hessian and sandwich with offsets and zero weights.
f <- scatlss(); set.seed(52); n <- 18L
X <- cbind(1,runif(n),1,runif(n),1,runif(n)); attr(X,"lpi") <- list(1:2,3:4,5:6)
beta <- c(.3,.2,-.7,.1,1,.1); y <- rt(n,5); w <- rep(c(0,1,2),6)
off <- list(seq(-.1,.1,length.out=n),rep(.1,n),rep(-.2,n))
l <- f$ll(y,X,beta,w,f,off,deriv=1,ncv=TRUE)
score <- X
for(k in 1:3) score[,attr(X,"lpi")[[k]]] <- X[,attr(X,"lpi")[[k]]]*l$l1[,k]
near(l$lb,colSums(score)); near(f$sandwich(y,X,beta,w,f,off),crossprod(score),1e-12)
for(k in seq_along(beta)) {
  bp <- bm <- beta; bp[k] <- bp[k]+1e-5; bm[k] <- bm[k]-1e-5
  lp <- f$ll(y,X,bp,w,f,off,deriv=1); lm <- f$ll(y,X,bm,w,f,off,deriv=1)
  near(l$lb[k],(lp$l-lm$l)/2e-5); near(l$lbb[,k],(lp$lb-lm$lb)/2e-5)
}
y[1] <- .Machine$double.xmax; et <- matrix(0,n,3); et[1,] <- 1000
z <- f$ll(y,X,beta,w,f,eta=et,deriv=4,ncv=TRUE)
stopifnot(z$l0[1]==0,all(z$l1[1,]==0),all(z$l4[1,]==0))
# Fixed-parameter scat comparisons do not equate different REML treatments.
nu <- 7; sig <- .8; yy <- seq(-4,4,length.out=n); mu <- rep(.2,n)
s <- mgcv::scat(theta=c(nu,sig))
eta <- cbind(mu,log(sig),log(nu-3))
near(sum(evaluate(yy,eta)$ell[,1]),-.5*s$aic(yy,mu,wt=rep(1,n),dev=0),1e-12)
near(f$dev.resids(yy,cbind(mu,sig,nu),w),s$dev.resids(yy,mu,w),1e-12)
# Custom lower bound is present in every link derivative.
near(evaluate(yy,eta,min.df=2.2)$ell[,1], dt((yy-mu)/sig,2.2+exp(eta[,3]),log=TRUE)-log(sig),1e-12)
