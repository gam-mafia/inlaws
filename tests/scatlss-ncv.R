library("inlaws")
near <- function(a,b,tol=1e-5) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
set.seed(23); n <- 450L
d <- data.frame(x=runif(n),z=runif(n),a=runif(n,-.1,.1),b=runif(n,-.1,.1),c=runif(n,-.05,.05))
f <- scatlss()
d$y <- f$rd(cbind(sin(6*d$x)+d$a,exp(-.5+.2*d$z+d$b),3+exp(.3+d$x+d$c)))
form <- list(y~s(x,k=5)+offset(a),~s(z,k=4)+offset(b),~x+offset(c))
ns <- asNamespace("mgcv"); G <- gam(form,data=d,family=f,fit=FALSE)
Sl <- get("Sl.setup",ns)(G); X <- get("Sl.initial.repara",ns)(Sl,G$X,both.sides=FALSE)
rho <- rep(log(2),length(G$sp)); start <- NULL
# Different numbers of deletion/prediction rows, including reordered/repeated
# predictions and repeated deletion groups, plus offsets in all three predictors.
nei <- list(a=c(4,5,6,10,11,12,4,5,6,19,20,21),ma=c(3,6,9,12),
            d=c(5,11,12,5,20),md=c(1,3,4,5),jackknife=FALSE)
call <- function(rho,method="REML",deriv=1,q=FALSE,neighbours=nei) {
  ff <- f; ff$qapprox <- q
  get("gam.fit5",ns)(X,G$y,rho,Sl,weights=G$w,offset=G$offset,family=ff,
    deriv=deriv,scoreType=method,Mp=-1,nei=neighbours,gamma=1.2,start=start,
    control=gam.control(epsilon=1e-12))
}
base <- call(rho,deriv=2); start <- base$coefficients
for(j in seq_along(rho)) {
  step <- numeric(length(rho)); step[j] <- .01
  mm <- call(rho-2*step); m <- call(rho-step); p <- call(rho+step); pp <- call(rho+2*step)
  near(base$REML1[j],(mm$REML-8*m$REML+8*p$REML-pp$REML)/.12,5e-5)
  near(base$REML2[,j],(mm$REML1-8*m$REML1+8*p$REML1-pp$REML1)/.12,2e-4)
}
for(q in c(FALSE,TRUE)) {
  fit <- call(rho,"NCV",q=q)
  for(j in seq_along(rho)) {
    rp <- rm <- rho; rp[j] <- rp[j]+1e-4; rm[j] <- rm[j]-1e-4
    near(fit$NCV1[j],(call(rp,"NCV",q=q)$NCV-call(rm,"NCV",q=q)$NCV)/2e-4,1e-5)
  }
  # Independently assemble each signed deletion Hessian and solve directly;
  # this guards against the mgcv 1.9-4 Cholesky update discrepancy.
  ll <- f$ll(G$y,X,start,G$w,f,G$offset,deriv=1,ncv=TRUE)
  eta <- getFromNamespace(".scatlss_eta","inlaws")(X,start,G$offset)
  # The penalized negative Hessian is penalty minus log-likelihood Hessian.
  H <- base$St - base$lbb
  direct <- eta[nei$d,,drop=FALSE]*0; pa <- pd <- 0L; jj <- attr(X,"lpi")
  for(i in seq_along(nei$ma)) {
    a <- nei$a[seq.int(pa+1L,nei$ma[i])]; id <- seq.int(pd+1L,nei$md[i]); pred <- nei$d[id]
    hd <- H; g <- numeric(ncol(X)); col <- 0L
    for(j in 1:3) {
      xx <- X[a,jj[[j]],drop=FALSE]; g[jj[[j]]] <- drop(crossprod(xx,ll$l1[a,j]))
      for(k in j:3) {
        col <- col+1L
        block <- crossprod(xx,X[a,jj[[k]],drop=FALSE]*ll$l2[a,col])
        hd[jj[[j]],jj[[k]]] <- hd[jj[[j]],jj[[k]],drop=FALSE]+block
        if(j!=k) hd[jj[[k]],jj[[j]]] <- hd[jj[[k]],jj[[j]],drop=FALSE]+t(block)
      }
    }
    delta <- solve(hd,g)
    for(j in 1:3) direct[id,j] <- eta[pred,j]-drop(X[pred,jj[[j]],drop=FALSE]%*%delta[jj[[j]]])
    pa <- nei$ma[i]; pd <- nei$md[i]
  }
  near(attr(fit$NCV,"eta.cv"),direct,1e-8)
  if(!q) {
    lp <- getFromNamespace(".scatlss_eval","inlaws")(G$y[nei$d],direct)$ell[,1]
    near(fit$NCV,-1.2*sum(lp)+.2*sum(ll$l0[nei$d]),1e-8)
    # Actual deleted-data refits use the same basis and penalties.
    pa <- pd <- 0L
    for(i in seq_along(nei$ma)) {
      a <- nei$a[seq.int(pa+1L,nei$ma[i])]; id <- seq.int(pd+1L,nei$md[i]); w <- G$w; w[a] <- 0
      ref <- get("gam.fit5",ns)(X,G$y,rho,Sl,weights=w,offset=G$offset,
        family=f,deriv=0,Mp=-1,start=start,control=gam.control(epsilon=1e-12))
      stopifnot(max(abs(direct[id,,drop=FALSE]-ref$linear.predictors[nei$d[id],,drop=FALSE])) < .1)
      pa <- nei$ma[i]; pd <- nei$md[i]
    }
  }
}
# Negative deriv requests return coefficient perturbations for jackknife
# uncertainty. They must describe the same deletions as the scoring route.
ll <- f$ll(G$y,X,start,G$w,f,G$offset,deriv=1,ncv=TRUE); ll$gamma <- 1
ff <- f; ff$qapprox <- FALSE
jk <- ff$ncv(X,G$y,G$w,nei,start,ff,ll,R=chol(base$St-base$lbb),offset=G$offset,deriv=-1)
near(attr(jk$NCV,"eta.cv"),direct,1e-8)
stopifnot(identical(dim(attr(jk$NCV,"deta.cv")),c(length(nei$ma),ncol(X))))

# An independent reference computes every deleted Hessian explicitly and
# evaluates Student-t densities directly; no native deletion update is called.
reference <- function(rho, q, neighbours = nei, gamma = 1.2) {
  b <- call(rho, "REML")
  be <- b$coefficients; jj <- attr(X,"lpi"); p <- ncol(X)
  eta <- getFromNamespace(".scatlss_eta","inlaws")(X,be,G$offset)
  ll <- f$ll(G$y,X,be,G$w,f,G$offset,deriv=1,ncv=TRUE)
  H <- b$St-b$lbb; nd <- length(neighbours$d)
  cv <- matrix(0,nd,3); changes <- matrix(0,length(neighbours$ma),p)
  pa <- pd <- 0L
  for(i in seq_along(neighbours$ma)) {
    a <- neighbours$a[seq.int(pa+1L,neighbours$ma[i])]
    id <- seq.int(pd+1L,neighbours$md[i]); pred <- neighbours$d[id]
    # Observation-level design maps all three predictors simultaneously.
    hd <- H; g <- numeric(p)
    for(row in a) {
      Z <- matrix(0,3,p)
      for(k in 1:3) Z[k,jj[[k]]] <- X[row,jj[[k]]]
      g <- g + drop(crossprod(Z,ll$l1[row,]))
      hd <- hd + crossprod(Z,matrix(ll$l2[row,f$tri$i2],3,3) %*% Z)
    }
    changes[i,] <- solve(hd,g)
    for(k in 1:3) cv[id,k] <- eta[pred,k]-drop(X[pred,jj[[k]],drop=FALSE]%*%changes[i,jj[[k]]])
    pa <- neighbours$ma[i]; pd <- neighbours$md[i]
  }
  ix <- neighbours$d
  logdensity <- function(eta) dt((G$y[ix]-eta[,1])/exp(eta[,2]),3+exp(eta[,3]),log=TRUE)-eta[,2]
  base <- G$w[ix]*logdensity(eta[ix,,drop=FALSE])
  if(!q) rows <- -gamma*G$w[ix]*logdensity(cv) - (1-gamma)*base else {
    delta <- cv-eta[ix,,drop=FALSE]
    rows <- -base
    for(i in seq_len(nd)) rows[i] <- rows[i] - gamma*(sum(ll$l1[ix[i],]*delta[i,]) +
      .5*drop(crossprod(delta[i,],matrix(ll$l2[ix[i],f$tri$i2],3,3)%*%delta[i,])))
  }
  list(score=sum(rows),rows=rows,eta=cv,changes=changes)
}
ll <- f$ll(G$y,X,start,G$w,f,G$offset,deriv=1,ncv=TRUE)
# Tail observations have convex location log likelihood, and signed mixed
# curvature must be retained rather than forced positive for a downdate.
stopifnot(any(ll$l2[,1]>0),any(ll$l2[,2]<0),any(ll$l2[,2]>0))
for(q in c(FALSE,TRUE)) {
  fit <- call(rho,"NCV",q=q); ref <- reference(rho,q)
  near(fit$NCV,ref$score,1e-8)
  near(attr(fit$NCV,"eta.cv"),ref$eta,1e-8)
  row_deriv <- matrix(0,length(nei$d),length(rho))
  for(k in seq_along(rho)) {
    rp <- rm <- rho; rp[k] <- rp[k]+1e-4; rm[k] <- rm[k]-1e-4
    plus <- reference(rp,q); minus <- reference(rm,q)
    near(fit$NCV1[k],(plus$score-minus$score)/2e-4,2e-5)
    row_deriv[,k] <- (plus$rows-minus$rows)/2e-4
  }
  near(fit$Vg,crossprod(row_deriv),2e-5)
}
near(attr(jk$NCV,"deta.cv"),reference(rho,FALSE)$changes,1e-8)
# Default leave-one-out structure uses the same signed reference and gamma adjustment.
loo <- list(a=seq_len(n),ma=seq_len(n),d=seq_len(n),md=seq_len(n),jackknife=FALSE)
for(q in c(FALSE,TRUE)) {
  fit <- call(rho,"NCV",q=q,neighbours=loo)
  near(fit$NCV,reference(rho,q,loo)$score,1e-7)
}
# Nonuniform, zero and fractional weights in grouped deletion solves.
G$w <- rep(c(0,.5,2,1),length.out=n)
base <- call(rho,deriv=2); start <- base$coefficients
for(q in c(FALSE,TRUE)) {
  fit <- call(rho,"NCV",q=q); ref <- reference(rho,q)
  near(fit$NCV,ref$score,1e-7)
  for(k in seq_along(rho)) {
    rp <- rm <- rho; rp[k] <- rp[k]+1e-4; rm[k] <- rm[k]-1e-4
    near(fit$NCV1[k],(reference(rp,q)$score-reference(rm,q)$score)/2e-4,3e-5)
  }
}
