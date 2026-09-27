# Retain the cgamma response contract while giving general.family a vector
# (mgcv counts length(y) before evaluating initialize).
.cgammals_response <- function(y) {
  r <- .cg_response(y)
  out <- r$y
  attr(out, "censor") <- r$censor
  out
}

# Bivariate Taylor coefficients in (log(mu), log(phi)), through total order
# four. The polynomial arithmetic is shared with the lognormal evaluator;
# gamma-specific probability recurrences below never numerically difference.
.cgammals_eval <- function(y, eta1, eta2, order = 0L) {
  r <- .cg_response(y)
  n <- length(r$y)
  eta1 <- rep_len(eta1, n); eta2 <- rep_len(eta2, n)
  j <- .clnormal_jets(order, order)
  C <- j$constant; M <- j$multiply
  inv <- function(x) j$compose(x, lapply(0:order, function(k) (-1)^k/x[, 1]^(k+1)))
  u <- C(eta1); t <- C(eta2)
  if (order) {
    u[, which(j$powers[, 1] == 1 & j$powers[, 2] == 0)] <- 1
    t[, which(j$powers[, 1] == 0 & j$powers[, 2] == 1)] <- 1
  }
  K <- j$exp(-t)
  lgammajet <- function(k) {
    coeff <- list(lgamma(k[, 1]))
    if (order) for (a in seq_len(order))
      coeff[[a+1L]] <- psigamma(k[, 1], deriv = a-1L)/factorial(a)
    j$compose(k, coeff)
  }
  # log(exp(a)-exp(b)); b may represent a zero probability.
  logdiff <- function(a, b) {
    out <- a
    ii <- is.finite(b[, 1])
    if (any(ii)) {
      delta <- b[ii, , drop = FALSE] - a[ii, , drop = FALSE]
      z <- -j$exp(delta)
      z[, 1] <- -expm1(delta[, 1])
      out[ii, ] <- a[ii, , drop = FALSE] + j$log(z)
    }
    out
  }
  # Calculate whichever incomplete-gamma expansion converges naturally, then
  # complement if needed. Values are replaced with R's accurate log tails.
  tail <- function(bound, U, T, k, lower) {
    x0 <- exp(log(bound) - U[, 1] - T[, 1])
    value <- stats::pgamma(x0, k[, 1], lower.tail = lower, log.p = TRUE)
    out <- C(value)
    if (!order) return(out)
    for (series in c(TRUE, FALSE)) {
      ii <- which(x0 > 0 & is.finite(x0) & ((x0 < k[, 1]+1) == series))
      if (!length(ii)) next
      logx <- C(log(bound[ii])) - U[ii, , drop = FALSE] - T[ii, , drop = FALSE]
      X <- j$exp(logx); kk <- k[ii, , drop = FALSE]
      one <- C(rep(1, length(ii)))
      pref <- M(kk, logx) - X - lgammajet(kk)
      converged <- FALSE
      if (series) {
        term <- total <- one
        for (a in seq_len(10000L)) {
          term <- M(term, M(X, inv(kk + a*one)))
          total <- total + term
          if (all(is.finite(term)) && all(abs(term) < 2e-14*(1+abs(total)))) {
            converged <- TRUE
            break
          }
          if (any(!is.finite(total))) break
        }
        ans <- pref - j$log(kk) + j$log(total)
      } else {
        b <- X + one - kk
        cc <- C(rep(1e300, length(ii)))
        d <- inv(b); h <- d
        for (a in seq_len(10000L)) {
          aa <- a*(kk - a*one)
          b <- b + 2*one
          d <- inv(b + M(aa, d))
          cc <- b + M(aa, inv(cc))
          delta <- M(cc, d)
          h <- M(h, delta)
          if (all(is.finite(delta)) && all(abs(delta-one) < 2e-14)) {
            converged <- TRUE
            break
          }
          if (any(!is.finite(h))) break
        }
        ans <- pref + j$log(h)
      }
      if (!converged || any(!is.finite(ans)))
        stop("cgammals: incomplete-gamma derivatives failed to converge")
      ans[, 1] <- stats::pgamma(x0[ii], kk[, 1], lower.tail = series, log.p = TRUE)
      if (series != lower) ans <- logdiff(0*one, ans)
      out[ii, ] <- ans
    }
    out[, 1] <- value
    out
  }
  density <- function(y, U, T, k) {
    logx <- C(log(y)) - U - T
    ans <- M(k, logx) - j$exp(logx) - lgammajet(k) - C(log(y))
    ans[, 1] <- stats::dgamma(y, k[, 1], scale = exp(U[, 1]+T[, 1]), log = TRUE)
    ans
  }
  ell <- C(numeric(n))
  ii <- which(r$exact)
  if (length(ii)) ell[ii, ] <- density(r$y[ii], u[ii, , drop = FALSE],
                                      t[ii, , drop = FALSE], K[ii, , drop = FALSE])
  narrow <- !r$exact & .cg_narrow(r$lo, r$hi, exp(eta1), K[, 1])
  regular <- !r$exact & !narrow
  lower <- stats::pgamma(exp(log(r$lo)-eta1-eta2), K[, 1]) < 0.5
  for (lt in c(TRUE, FALSE)) {
    ii <- which(regular & lower == lt)
    if (!length(ii)) next
    U <- u[ii, , drop = FALSE]; T <- t[ii, , drop = FALSE]; k <- K[ii, , drop = FALSE]
    a <- tail(if (lt) r$hi[ii] else r$lo[ii], U, T, k, lt)
    b <- tail(if (lt) r$lo[ii] else r$hi[ii], U, T, k, lt)
    ell[ii, ] <- logdiff(a, b)
  }
  ii <- which(narrow)
  if (length(ii)) {
    rule <- .gauss_legendre(16L)
    U <- u[ii, , drop = FALSE]; T <- t[ii, , drop = FALSE]; k <- K[ii, , drop = FALSE]
    width <- r$hi[ii]-r$lo[ii]
    centre <- density(r$lo[ii]+width/2, U, T, k)
    total <- C(numeric(length(ii)))
    for (a in seq_along(rule$nodes)) {
      yy <- r$lo[ii] + width*(rule$nodes[a]+1)/2
      total <- total + rule$weights[a]*j$exp(density(yy, U, T, k)-centre)
    }
    ell[ii, ] <- j$log(total) + centre + C(log(width/2))
  }
  list(ell = ell, powers = j$powers)
}
