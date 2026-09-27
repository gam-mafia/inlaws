# Pure-R CMP reference numerics. Taylor coefficients include factorial divisors.
.cmp_control <- function(control = list()) {
  defaults <- list(sum_tol = 1e-10, mean_tol = 1e-10,
                   max_terms = 100000L, max_iter = 100L)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))) ||
       anyDuplicated(names(control)) || any(!names(control) %in% names(defaults)))))
    stop("invalid or unknown cmp control", call. = FALSE)
  defaults[names(control)] <- control
  for (nm in names(defaults)) {
    x <- defaults[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 ||
        (nm %in% c("sum_tol", "mean_tol") && x >= 1) ||
        (nm %in% c("max_terms", "max_iter") && (x != floor(x) || x > .Machine$integer.max)))
      stop("invalid cmp control: ", nm, call. = FALSE)
  }
  defaults
}

.cmp_index <- function(order) {
  do.call(rbind, lapply(0:order, function(d) cbind(d:0, 0:d)))
}
.cmp_i4 <- .cmp_index(4L)
.cmp_i5 <- .cmp_index(5L)

# Relative weights around a mode. Bound both omitted tails for every absolute
# mixed moment through order five about (mode, log(mode!)). See inst/maths/cmp.
.cmp_series <- function(a, nu, control) {
  if (!is.finite(a) || !is.finite(nu) || nu <= 0)
    stop("CMP requires finite log(lambda) and positive finite nu", call. = FALSE)
  mode <- floor(exp(a / nu))
  if (!is.finite(mode) || mode > 2^50)
    stop("CMP mode exceeds numerical range", call. = FALSE)
  lm <- lgamma(mode + 1)
  radius <- min(control$max_terms, max(8, ceiling(8 * sqrt((mode + 1) / nu))))
  lo <- max(0, mode - floor((radius - 1) / 2))
  hi <- lo + radius - 1
  repeat {
    if (hi - lo + 1 > control$max_terms)
      stop("CMP summation failed: max_terms exhausted", call. = FALSE)
    j <- seq.int(lo, hi)
    lf <- lgamma(j + 1)
    lw <- (j - mode) * a - nu * (lf - lm)
    w <- exp(lw)
    sw <- sum(w)
    dx <- abs(j - mode); dt <- abs(lf - lm)
    ok <- TRUE
    for (h in seq_len(nrow(.cmp_i5))) {
      p <- .cmp_i5[h, 1]; q <- .cmp_i5[h, 2]
      sm <- sum(w * dx^p * dt^q)
      # Upper tail starts at k. log(j) <= log(k) * j/k, k >= 2.
      k <- hi + 1; d <- k - mode
      lr <- a - nu * log(k + 1) + (p + q) * log1p(1 / d) + q * log1p(1 / k)
      ub <- Inf
      if (k >= 3 && lr < 0) {
        logb <- (k - mode) * a - nu * (lgamma(k + 1) - lm) +
          (p + q) * log(d) + q * log(log(k)) - log(-expm1(lr))
        ub <- exp(logb)
      }
      lb <- 0
      if (lo > 0) {
        k <- lo - 1; d <- mode - k
        lr <- if (k == 0) -Inf else -a + nu * log(k) + (p + q) * log1p(1 / d)
        lb <- Inf
        if (lr < 0) {
          factor <- if (q == 0) 0 else q * log(log(max(mode, 1)))
          lb <- exp((k - mode) * a - nu * (lgamma(k + 1) - lm) +
                      (p + q) * log(d) + factor - log(-expm1(lr)))
        }
      }
      # An identically zero absolute moment needs no relative accuracy.
      target <- control$sum_tol * max(sm, .Machine$double.xmin)
      if (!is.finite(ub + lb) || ub + lb > target) { ok <- FALSE; break }
    }
    if (ok) break
    extra <- max(8, ceiling((hi - lo + 1) / 2))
    lo <- max(0, lo - extra); hi <- hi + extra
  }
  pr <- w / sw
  mean <- sum(pr * j)
  list(a = a, nu = nu, y = j, lf = lf, p = pr,
       logZ = mode * a - nu * lm + log(sw), mean = mean,
       variance = sum(pr * (j - mean)^2))
}

.cmp_mean_one <- function(mu, nu, control) {
  if (!is.finite(mu) || mu < 0 || !is.finite(nu) || nu <= 0)
    stop("CMP requires a nonnegative finite mean and positive finite nu", call. = FALSE)
  if (mu == 0) return(list(a = -Inf, nu = nu, y = 0, lf = 0, p = 1,
                          logZ = 0, mean = 0, variance = 0))
  a <- if (mu < 1) log(mu) else nu * log(mu + max(0, (nu - 1) / (2 * nu)))
  if (!is.finite(a)) a <- log(mu)
  lower <- -Inf; upper <- Inf
  for (iter in seq_len(control$max_iter)) {
    z <- .cmp_series(a, nu, control)
    err <- z$mean - mu
    if (abs(err) <= control$mean_tol * mu) return(z)
    if (err < 0) lower <- a else upper <- a
    step <- err / z$variance
    next_a <- a - step
    if (!is.finite(next_a) || next_a <= lower || next_a >= upper || abs(step) > 10) {
      next_a <- if (is.finite(lower) && is.finite(upper)) (lower + upper) / 2 else
        a + if (err < 0) 1 else -1
    }
    a <- next_a
  }
  stop("CMP mean inversion failed: max_iter exhausted", call. = FALSE)
}

# Strict scalar recycling, also used by the distribution callbacks.
.cmp_recycle <- function(...) {
  x <- list(...)
  if (any(lengths(x) == 0L)) return(lapply(x, function(z) z[integer()]))
  n <- max(lengths(x))
  if (any(!lengths(x) %in% c(1L, n))) stop("incompatible CMP argument lengths", call. = FALSE)
  lapply(x, rep_len, length.out = n)
}
.cmp_moments <- function(mu, nu, control = .cmp_control()) {
  v <- .cmp_recycle(mu, nu)
  Map(function(m, n) .cmp_mean_one(m, n, control), v[[1]], v[[2]])
}

# Polynomial multiplication on a total-degree triangle, vectorised over rows.
.cmp_mul <- function(x, y, index = .cmp_i4) {
  out <- x * 0
  for (i in seq_len(nrow(index))) {
    if (!any(x[, i] != 0)) next
    for (j in seq_len(nrow(index))) {
      ij <- index[i, ] + index[j, ]
      if (sum(ij) > max(rowSums(index)) || !any(y[, j] != 0)) next
      k <- (sum(ij) * (sum(ij) + 1)) / 2 + ij[2] + 1
      out[, k] <- out[, k] + x[, i] * y[, j]
    }
  }
  out
}

# Taylor coefficients of A(a + da, nu + dn). Centering protects cumulants
# against cancellation; coefficients are the log of the moment generating jet.
.cmp_A <- function(z) {
  n <- length(z); m <- matrix(0, n, nrow(.cmp_i5))
  means <- vapply(z, `[[`, numeric(1), "mean")
  ets <- vapply(z, function(s) sum(s$p * s$lf), numeric(1))
  for (r in seq_len(n)) {
    x <- z[[r]]$y - means[r]; t <- -(z[[r]]$lf - ets[r])
    for (k in 4:ncol(m)) {
      i <- .cmp_i5[k, 1]; j <- .cmp_i5[k, 2]
      m[r, k] <- sum(z[[r]]$p * x^i * t^j) / (factorial(i) * factorial(j))
    }
  }
  # m has no constant or linear term, so log(1+m) = m - m^2/2 to order 5.
  ans <- m - .cmp_mul(m, m, .cmp_i5) / 2
  ans[, 1] <- vapply(z, `[[`, numeric(1), "logZ")
  ans[, 2] <- means; ans[, 3] <- -ets
  ans
}
.cmp_compose <- function(coef, x, y) {
  one <- x * 0; one[, 1] <- 1
  xp <- yp <- vector("list", 5); xp[[1]] <- yp[[1]] <- one
  for (i in 2:5) {
    xp[[i]] <- .cmp_mul(xp[[i - 1]], x)
    yp[[i]] <- .cmp_mul(yp[[i - 1]], y)
  }
  out <- one * 0
  for (k in seq_len(ncol(coef))) {
    i <- .cmp_i4[k, 1]; j <- .cmp_i4[k, 2]
    out <- out + coef[, k] * .cmp_mul(xp[[i + 1]], yp[[j + 1]])
  }
  out
}

# Likelihood jet in (mu, log(nu)), or (log(mu), log(nu)) for predictors,
# obtained by formal implicit inversion.
.cmp_lljet <- function(y, mu, nu, z, predictor = FALSE) {
  A <- .cmp_A(z); n <- length(z)
  dnu <- dmu <- da <- matrix(0, n, 15)
  for (k in 1:4) dnu[, k * (k + 1) / 2 + k + 1] <- nu / factorial(k)
  if (predictor) {
    for (k in 1:4) dmu[, k * (k + 1) / 2 + 1] <- mu / factorial(k)
  } else dmu[, 2] <- 1
  Aa <- A[, 1:15, drop = FALSE] * 0
  for (k in 1:15) {
    ij <- .cmp_i4[k, ] + c(1, 0)
    col <- sum(ij) * (sum(ij) + 1) / 2 + ij[2] + 1
    Aa[, k] <- A[, col] * ij[1]
  }
  Aa[, 1] <- 0
  variance <- vapply(z, `[[`, numeric(1), "variance")
  if (any(!is.finite(variance) | variance <= 0))
    stop("CMP variance is numerically zero or nonfinite", call. = FALSE)
  for (i in 1:4) da <- da + (dmu - .cmp_compose(Aa, da, dnu)) / variance
  Ac <- A[, 1:15, drop = FALSE]; Ac[, 1] <- 0
  out <- y * da - lgamma(y + 1) * dnu - .cmp_compose(Ac, da, dnu)
  out[, 1] <- y * vapply(z, `[[`, numeric(1), "a") - nu * lgamma(y + 1) - A[, 1]
  if (any(!is.finite(out)) || any(!is.finite(da)))
    stop("CMP likelihood derivatives exceed numerical range", call. = FALSE)
  # Return da for exact expected observed derivatives (likelihood affine in y).
  list(ll = out, da = da)
}
