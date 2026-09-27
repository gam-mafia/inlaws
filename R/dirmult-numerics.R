## Stable log rising factorial minus m*log(a), and its first two
## derivatives with respect to a, scaled by a and a^2 respectively.
## Counts m are non-negative integers. Working with log(a) also avoids
## underflow of small category concentrations.
.dirmult_rising <- function(la,m) {
  r <- A <- B <- numeric(length(m))
  small <- m<=64
  if (any(small)) {
    ii <- which(small & m>0)
    if (length(ii)) {
      A[ii] <- 1; B[ii] <- -1
      for (k in seq_len(max(m[ii])-1)) {
        jj <- ii[m[ii]>k]
        z <- log(k)-la[jj]
        r[jj] <- r[jj] + pmax(z,0) + log1p(exp(-abs(z)))
        q <- stats::plogis(-z)
        A[jj] <- A[jj]+q; B[jj] <- B[jj]-q*q
      }
    }
  }
  ## Series in (m-1)/a, with a uniformly small remainder. In particular
  ## this does not subtract two almost equal lgamma/digamma values.
  ii <- which(!small & la>log(m)+log(1e3))
  if (length(ii)) {
    t <- m[ii]; ia <- exp(-la[ii]); u <- (t-1)*ia
    s1 <- t*u/2
    s2 <- t*u*(2*t-1)*ia/6
    s3 <- t*t*u*u*ia/4
    s4 <- t*u*(2*t-1)*ia*(3*t*t-3*t-1)*ia*ia/30
    r[ii] <- s1-s2/2+s3/3-s4/4
    A[ii] <- t-s1+s2-s3+s4
    B[ii] <- -t+2*s1-3*s2+4*s3-5*s4
  }
  ii <- which(!small & !(la>log(m)+log(1e3)))
  if (length(ii)) {
    a <- exp(la[ii]); t <- m[ii]
    tiny <- a<1
    jj <- ii[tiny]; aa <- a[tiny]; tt <- t[tiny]
    r[jj] <- lgamma(aa+tt)-lgamma(1+aa)-(tt-1)*la[jj]
    A[jj] <- 1+aa*(digamma(aa+tt)-digamma(1+aa))
    B[jj] <- -1+aa^2*(trigamma(aa+tt)-trigamma(1+aa))
    jj <- ii[!tiny]; aa <- a[!tiny]; tt <- t[!tiny]
    r[jj] <- lgamma(aa+tt)-lgamma(aa)-tt*la[jj]
    A[jj] <- aa*(digamma(aa+tt)-digamma(aa))
    B[jj] <- aa^2*(trigamma(aa+tt)-trigamma(aa))
  }
  list(r=r,A=A,B=B)
}

.dirmult_prob <- function(eta) {
  z <- cbind(0,eta[,-ncol(eta),drop=FALSE])
  z <- z-apply(z,1,max)
  lp <- z-log(rowSums(exp(z)))
  list(p=exp(lp),lp=lp,logphi=eta[,ncol(eta)])
}

## Observation derivatives on the linear predictor scale. If J maps
## predictors to log(alpha), H = J'diag(A+B)J + sum(A_j H(log(alpha_j))).
.dirmult_deriv <- function(y,eta,wt=rep(1,nrow(y)),deriv=TRUE) {
  n <- nrow(y); C <- ncol(y); K <- C-1L
  pr <- .dirmult_prob(eta); p <- pr$p
  total <- rowSums(y)
  a0 <- .dirmult_rising(pr$logphi,total)
  a <- .dirmult_rising(as.vector(pr$lp+pr$logphi),as.vector(y))
  A <- matrix(a$A,n,C); D <- matrix(a$A+a$B,n,C)
  l0 <- lgamma(total+1)-rowSums(lgamma(y+1)) +
        rowSums(y*pr$lp+matrix(a$r,n,C))-a0$r
  ret <- list(l=sum(wt*l0),l0=l0)
  if (!deriv) return(ret)
  sa <- rowSums(A); sd <- rowSums(D)
  l1 <- cbind(A[,-1,drop=FALSE]-p[,-1,drop=FALSE]*sa,sa-a0$A)
  l2 <- matrix(0,n,C*(C+1)/2)
  h <- 0L
  for (j in seq_len(C)) for (k in j:C) {
    h <- h+1L
    if (k<C) {
      pj <- p[,j+1]; pk <- p[,k+1]
      l2[,h] <- pj*pk*(sd+sa)-pj*D[,k+1]-pk*D[,j+1]
      if (j==k) l2[,h] <- l2[,h]+D[,j+1]-sa*pj
    } else if (j<C) l2[,h] <- D[,j+1]-p[,j+1]*sd else
      l2[,h] <- sd-a0$A-a0$B
  }
  ret$l1 <- l1*wt; ret$l2 <- l2*wt
  ret
}


## Intercept-only reference with the same offsets and row weights.
.dirmult_null <- function(y,wt,offset,eta) {
  C <- ncol(y); n <- nrow(y)
  off <- matrix(0,n,C)
  for (j in seq_len(C)) if (length(offset)>=j && !is.null(offset[[j]]))
    off[,j] <- offset[[j]]
  p <- colSums(y*wt)+0.5; p <- p/sum(p)
  b <- c(log(p[-1]/p[1]),stats::weighted.mean(eta[,C],wt))
  fn <- function(b) -.dirmult_deriv(y,sweep(off,2,b,"+"),wt,FALSE)$l
  gr <- function(b) -colSums(.dirmult_deriv(y,sweep(off,2,b,"+"),wt)$l1)
  fit <- stats::optim(b,fn,gr,method="BFGS",control=list(maxit=200,reltol=1e-10))
  if (fit$convergence!=0) {
    warning("dirmult intercept-only reference did not converge")
    return(NA_real_)
  }
  2*fit$value
}
