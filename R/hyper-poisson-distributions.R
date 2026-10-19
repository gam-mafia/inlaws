.hp_parameters <- function(mu) {
  if (!is.matrix(mu) || !is.numeric(mu) || ncol(mu) != 2L)
    stop("hyper-Poisson parameters must be a numeric matrix with mu and sigma columns")
  if (any(!is.na(mu[,1]) & (!is.finite(mu[,1]) | mu[,1] < 0)) ||
      any(!is.na(mu[,2]) & (!is.finite(mu[,2]) | mu[,2] <= 0)))
    stop("hyper-Poisson requires finite nonnegative means and positive dispersions")
}

.hp_log1mexp <- function(x) {
  ifelse(x < -log(2), log1p(-exp(x)), log(-expm1(x)))
}

# Sum a tail starting on the requested side of the mode, on a relative
# log scale. Its monotone recurrence supplies a geometric remainder bound.
.hp_tail <- function(k, z, up, control) {
  at <- k; lw <- 0; top <- 0; sw <- 1
  for (i in seq_len(control$max_terms)) {
    if (!up && at == 0) return(min(0,.hp_logmass(k,z)+top+log(sw)))
    lr <- if (up) z$a-log(z$sigma+at) else log(z$sigma+at-1)-z$a
    # Relative error of the remaining sum, including its first omitted term.
    if (lr < 0 && lw+lr-log(-expm1(lr)) <= log(control$sum_tol)+top+log(sw))
      return(min(0,.hp_logmass(k,z)+top+log(sw)))
    lw <- lw+lr; at <- at + if (up) 1 else -1
    if (lw > top) { sw <- sw*exp(top-lw)+1; top <- lw } else sw <- sw+exp(lw-top)
  }
  stop("hyper-Poisson tail summation failed: max_terms exhausted",call.=FALSE)
}

.hp_cdf_one <- function(q,z,lower.tail,control) {
  if (q < 0) return(if (lower.tail) -Inf else 0)
  if (q == Inf || z$mean == 0) return(if (lower.tail) 0 else -Inf)
  q <- floor(q)
  # Sum the side away from the mode and complement on the log scale.
  low <- q < z$mode
  lp <- if (low) .hp_tail(q,z,FALSE,control) else .hp_tail(q+1,z,TRUE,control)
  if (low == lower.tail) lp else .hp_log1mexp(lp)
}

.hp_quantile_one <- function(lp,z,lower.tail,control) {
  if (z$mean == 0) return(0)
  if (lp == -Inf) return(if (lower.tail) 0 else Inf)
  if (lp == 0) return(if (lower.tail) Inf else 0)
  # Use the smaller target tail, preserving upper log probabilities near zero.
  small_lower <- if (lower.tail) lp <= -log(2) else lp > -log(2)
  target <- if (small_lower == lower.tail) lp else .hp_log1mexp(lp)
  # A stable complement can still differ by an ulp when a probability lies
  # exactly on a discrete CDF step. Keep that rounding from advancing a count.
  fuzz <- 8*.Machine$double.eps*max(1,abs(target))
  reached <- function(k) {
    v <- .hp_cdf_one(k,z,small_lower,control)
    if (small_lower) v >= target-fuzz else v <= target+fuzz
  }
  if (reached(0)) return(0)
  lo <- 0; hi <- max(1,ceiling(z$mean))
  for (iter in seq_len(control$max_iter)) {
    if (reached(hi)) break
    lo <- hi; hi <- 2*hi+1
    if (hi > 2^50 || iter == control$max_iter)
      stop("hyper-Poisson quantile bracketing failed",call.=FALSE)
  }
  while (hi-lo > 1) {
    k <- floor((hi+lo)/2)
    if (reached(k)) hi <- k else lo <- k
  }
  hi
}

.hp_probability <- function(q,mu,control,quantile=FALSE,lower.tail=TRUE,log.p=FALSE) {
  .hp_parameters(mu)
  if (!is.numeric(q) || !is.null(dim(q))) stop("hyper-Poisson probability arguments must be numeric vectors")
  for (flag in list(lower.tail,log.p))
    if (!is.logical(flag) || length(flag)!=1L || is.na(flag)) stop("tail and log flags must be TRUE or FALSE")
  n <- nrow(mu)
  if (!length(q) || !n) return(numeric())
  if (!length(q) %in% c(1L,n)) stop("probability arguments must be scalar or one per parameter row")
  q <- rep_len(q,n)
  if (quantile && any(!is.na(q) & if (log.p) q>0 else q<0 | q>1))
    stop("hyper-Poisson probabilities must lie in [0, 1] (or have nonpositive log probabilities)")
  out <- rep(NA_real_,n)
  keep <- which(!is.na(q) & !is.na(mu[,1]) & !is.na(mu[,2]))
  # The Poisson submodel has native probability functions, including their
  # discrete-boundary conventions. This shortcut is confined to callbacks:
  # the likelihood still retains every dispersion derivative at sigma = 1.
  poisson <- keep[mu[keep,2]==1]
  if (length(poisson)) {
    out[poisson] <- if (quantile)
      stats::qpois(q[poisson],mu[poisson,1],lower.tail=lower.tail,log.p=log.p) else
      stats::ppois(q[poisson],mu[poisson,1],lower.tail=lower.tail,log.p=log.p)
  }
  keep <- keep[mu[keep,2]!=1]
  z <- .hp_moments(mu[keep,1],mu[keep,2],control)
  for (j in seq_along(keep)) {
    i <- keep[j]
    if (quantile) out[i] <- .hp_quantile_one(if (log.p) q[i] else log(q[i]),z[[j]],lower.tail,control)
    else {
      lp <- .hp_cdf_one(q[i],z[[j]],lower.tail,control)
      out[i] <- if (log.p) lp else exp(lp)
    }
  }
  out
}

.hp_diagnostics <- function(family,mu,wt,scale,y,eta) {
  .hp_parameters(mu)
  variance <- rep(NA_real_,nrow(mu)); keep <- which(!is.na(mu[,1]) & !is.na(mu[,2]))
  variance[keep] <- vapply(.hp_moments(mu[keep,1],mu[keep,2],family$control),`[[`,numeric(1),"variance")
  cdf <- function(q,log_p=FALSE,lower_tail=TRUE)
    family$cdf(q,mu,wt,scale,logp=log_p,lower.tail=lower_tail)
  list(mean=mu[,1],variance=variance,continuous=FALSE,
    cdf=cdf,
    cdf_left=function(q,log_p=FALSE,lower_tail=TRUE) cdf(ceiling(q)-1,log_p,lower_tail),
    qf=function(p,log_p=FALSE,lower_tail=TRUE)
      family$qf(p,mu,wt,scale,lower.tail=lower_tail,log.p=log_p),
    rd=function() family$rd(mu,wt,scale),
    diagnostic_residuals=function(y,type) {
      r <- y-mu[,1]
      if (type=="response") return(r)
      if (type=="pearson") return(ifelse(wt==0,0,r*sqrt(wt/variance)))
      if (type=="deviance") return(sign(r)*sqrt(family$dev.resids(y,mu,wt)))
      stop("unsupported hyper-Poisson diagnostic residual type")
    })
}
