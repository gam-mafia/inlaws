# GP-1 likelihood and distribution calculations, independent of gam().
# t = log(phi - 1); t = -Inf denotes the exact Poisson boundary.
.gpoisson_parameters <- function(t) {
  lp <- pmax(t, 0) + log1p(exp(-abs(t)))
  c <- exp(-lp / 2)
  h <- stats::plogis(t)
  list(c = c, b = -expm1(-lp / 2), cp = -c * h / 2,
       cpp = c * (3 * h^2 / 4 - h / 2), h = h, lp = lp)
}

.gpoisson_logpmf <- function(y, mu, t) {
  p <- .gpoisson_parameters(t)
  if (t == -Inf) return(stats::dpois(y, mu, log = TRUE))
  z <- p$c * mu + p$b * y
  ans <- stats::dpois(y, z, log = TRUE) + log(mu) - p$lp / 2 - log(z)
  ans[y == 0] <- -p$c * mu[y == 0]
  ans[mu == 0 & y > 0] <- -Inf
  ans
}

# Ordinary analytic derivatives of log p, with respect to mu and t.
.gpoisson_lld <- function(y, mu, t, level = 2L) {
  p <- .gpoisson_parameters(t)
  z <- p$c * mu + p$b * y
  u <- p$c / z
  ans <- list(m = 1 / mu + (y - 1) * u - p$c,
              mm = -1 / mu^2 - (y - 1) * u^2)
  if (level > 0) {
    v <- p$cp * (mu - y) / z
    up <- (p$cp / z) * (y / z)
    ans$t <- -p$h / 2 + (y - 1) * v - p$cp * (mu - y)
    ans$mt <- (y - 1) * up - p$cp
    ans$mmm <- 2 / mu^3 + 2 * (y - 1) * u^3
    ans$mmt <- -2 * (y - 1) * u * up
  }
  if (level > 1) {
    upp <- (p$cpp / z) * (y / z) - 2 * up * v
    ans$tt <- -p$h * (1 - p$h) / 2 +
      (y - 1) * (p$cpp * (mu - y) / z - v^2) - p$cpp * (mu - y)
    ans$mtt <- (y - 1) * upp - p$cpp
    ans$mmtt <- -2 * (y - 1) * (up^2 + u * upp)
    ans$mmmt <- 6 * (y - 1) * u^2 * up
    ans$mmmm <- -6 / mu^4 - 6 * (y - 1) * u^4
  }
  # At zero the log likelihood is exactly -c*mu. Avoid cancellation.
  zero <- y == 0
  for (nm in names(ans)) ans[[nm]][zero] <- switch(nm,
    m = -p$c, t = -p$cp * mu[zero], mt = -p$cp,
    tt = -p$cpp * mu[zero], mtt = -p$cpp, 0)
  ans
}

.gpoisson_saturated <- function(y, t) {
  out <- matrix(0, length(y), 3L)
  pos <- y > 0
  if (any(pos)) {
    p <- .gpoisson_parameters(t)
    # Solve a^2 - y*(1-b)*a - b*y = 0; mu = a/(1-b).
    v <- y[pos] * p$c
    a <- (v + sqrt(v^2 + 4 * p$b * y[pos])) / 2
    mu <- a / p$c
    d <- .gpoisson_lld(y[pos], mu, t)
    out[pos, 1] <- .gpoisson_logpmf(y[pos], mu, t)
    out[pos, 2] <- d$t
    # Envelope theorem, including movement of the saturated mean.
    out[pos, 3] <- d$tt - d$mt^2 / d$mm
  }
  out
}

.gpoisson_recycle <- function(x, mu) {
  if (!is.numeric(x) || !is.numeric(mu)) stop("gpoisson distribution arguments must be numeric")
  if (!length(x) || !length(mu)) return(list(numeric(), numeric()))
  n <- max(length(x), length(mu))
  if (!all(c(length(x), length(mu)) %in% c(1L, n)))
    stop("gpoisson distribution arguments must have equal lengths or be scalars")
  if (any(!is.na(mu) & (!is.finite(mu) | mu < 0)))
    stop("gpoisson distribution means must be finite and nonnegative")
  list(rep_len(x, n), rep_len(mu, n))
}

.gpoisson_logadd <- function(a, b) {
  hi <- pmax(a, b)
  ans <- hi + log1p(exp(pmin(a, b) - hi))
  ans[hi == -Inf] <- -Inf
  ans
}

# Sum the actual probabilities, never a renormalized finite approximation.
# The finite CDF needs no tail approximation; quantiles stop at their crossing.
.gpoisson_probability <- function(x, mu, t, quantile = FALSE) {
  v <- .gpoisson_recycle(x, mu); x <- v[[1]]; mu <- v[[2]]
  if (quantile && any(!is.na(x) & (x < 0 | x > 1)))
    warning("NaNs produced")
  if (t == -Inf) {
    if (quantile) return(stats::qpois(x, mu))
    return(stats::ppois(x, mu, log.p = TRUE))
  }
  vapply(seq_along(x), function(i) {
    xi <- x[i]; mi <- mu[i]
    if (is.na(xi) || is.na(mi)) return(NA_real_)
    if (quantile) {
      if (xi < 0 || xi > 1) return(NaN)
      if (xi == 0 || mi == 0) return(0)
      if (xi == 1) return(Inf)
      target <- log(xi); last <- 1000000L - 1L
    } else {
      if (xi < 0) return(-Inf)
      if (mi == 0 || xi == Inf) return(0)
      last <- floor(xi)
      if (last >= 1000000L) stop("gpoisson CDF summation exceeds one million terms")
    }
    total <- -Inf
    for (first in seq.int(0, last, by = 512L)) {
      j <- seq.int(first, min(last, first + 511L))
      lp <- .gpoisson_logpmf(j, rep(mi, length(j)), t)
      peak <- max(lp)
      cumulative <- if (peak == -Inf) rep(-Inf, length(j)) else
        peak + log(cumsum(exp(lp - peak)))
      cumulative <- .gpoisson_logadd(total, cumulative)
      if (quantile && any(cumulative >= target))
        return(j[which(cumulative >= target)[1L]])
      total <- cumulative[length(j)]
    }
    if (quantile) stop("gpoisson quantile summation failed after one million terms")
    if (!is.finite(total) || total > 1e-10)
      stop("gpoisson CDF summation lost numerical accuracy")
    min(0, total)
  }, numeric(1))
}

# Poisson immigration followed by Poisson offspring, summed to extinction.
.gpoisson_random <- function(mu, t) {
  mu <- .gpoisson_recycle(mu, mu)[[2]]
  if (t == -Inf) return(stats::rpois(length(mu), mu))
  p <- .gpoisson_parameters(t)
  if (p$b >= 1) stop("gpoisson dispersion exceeds simulation precision")
  total <- generation <- stats::rpois(length(mu), mu * p$c)
  for (i in seq_len(1000000L)) {
    active <- which(!is.na(generation) & generation > 0)
    if (!length(active)) return(total)
    generation[active] <- stats::rpois(length(active), p$b * generation[active])
    total[active] <- total[active] + generation[active]
    if (any(!is.finite(total[active]) | total[active] > 2^53))
      stop("gpoisson simulation exceeds exact integer range")
  }
  stop("gpoisson simulation failed to reach extinction after one million generations")
}
