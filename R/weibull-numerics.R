.weibull_response <- function(y) {
  if (!is.numeric(y) || !is.null(dim(y)) || !is.null(attr(y, "censor")) ||
      any(!is.finite(y) | y <= 0))
    stop("Weibull families require a finite, positive, uncensored numeric response vector",
         call. = FALSE)
  y
}

.weibull_check_weights <- function(w) {
  if (any(!is.finite(w) | w < 0) || !is.finite(sum(w)) || !any(w > 0))
    stop("Weibull families require finite nonnegative weights with finite positive total",
         call. = FALSE)
}

# Coefficients of log Gamma(1 + u) about u = 0: -Euler and (-1)^j zeta(j)/j.
.weibull_lgamma_series <- c(-0.5772156649015329, 0.8224670334241132,
  -0.4006856343865314, 0.2705808084277845, -0.2073855510286740,
  0.1695571769974082, -0.1440498967688461, 0.1255096695247430,
  -0.1113342658695647, 0.1000994575127818)

# H(theta) = k * log Gamma(1 + 1/k), and derivatives in theta = log(k).
# Use the expansion at 1 to preserve corrections when 1 + 1/k rounds to 1.
# At small k, multiply Stirling's expansion through by k before evaluating.
.weibull_shift <- function(theta, order = 0L) {
  k <- exp(theta)
  h <- matrix(NA_real_, length(theta), order + 1L)
  large <- which(theta > log(100))
  small <- which(theta < -log(10000))
  mid <- which(abs(theta) < Inf & theta <= log(100) & theta >= -log(10000))
  u <- 1 / k[mid]
  d <- list(lgamma(1 + u))
  if (order >= 1L) d[[2L]] <- -u * digamma(1 + u)
  if (order >= 2L) {
    a <- u^2 * trigamma(1 + u)
    d[[3L]] <- -d[[2L]] + a
  }
  if (order >= 3L) {
    b <- u^3 * psigamma(1 + u, 2)
    d[[4L]] <- d[[2L]] - 3 * a - b
  }
  if (order >= 4L)
    d[[5L]] <- -d[[2L]] + 7 * a + 6 * b + u^4 * psigamma(1 + u, 3)
  powers <- seq_along(.weibull_lgamma_series) - 1L
  for (j in 0:order) {
    h[mid, j + 1L] <- k[mid] * Reduce(`+`, Map(`*`, as.list(choose(j, 0:j)), d[seq_len(j + 1L)]))
    terms <- .weibull_lgamma_series * (-powers)^j
    h[large, j + 1L] <- drop(exp(-outer(theta[large], powers)) %*% terms)
    ks <- k[small]
    h[small, j + 1L] <- switch(j + 1L, -theta[small] - 1, -1, 0, 0, 0) +
      ks * ((log(2 * pi) - theta[small] - j) / 2) +
      2^j * ks^2 / 12 - 4^j * ks^4 / 360 + 6^j * ks^6 / 1260
  }
  h
}

.weibull_log_ratio <- function(y, mu) {
  result <- log(y) - log(mu)
  close <- which(abs(result) < 0.5)
  result[close] <- log1p((y[close] - mu[close]) / mu[close])
  result
}

.weibull_deviance_core <- function(t) {
  dev <- expm1(t) - t
  dev[is.infinite(t)] <- Inf
  small <- which(abs(t) < 1e-4)
  x <- t[small]
  dev[small] <- x^2 * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2
  dev
}

# Compute the raw mixed partials, in the order (a,b) with a+b increasing
# and a decreasing within each order. Write t = k(log(y)-log(mu)) + H.
# Partial a in log(mu) of exp(t) is (-k)^a exp(t). Its theta derivatives
# are exp(t)(-k)^a times the exponential Bell polynomials in t' + a, t'', ...
.weibull_eval <- function(y, eta, theta, order = 0L, mu = NULL, deviance = FALSE) {
  theta <- rep_len(theta, length(y))
  k <- exp(theta)
  shift <- .weibull_shift(theta, order)
  delta <- if (is.null(mu)) log(y) - eta else .weibull_log_ratio(y, mu)
  v <- k * delta
  t <- v + shift[, 1L]
  sat <- theta - log(y) - 1
  out <- matrix(0, length(y), (order + 1L) * (order + 2L) / 2L)
  out[, 1L] <- if (deviance) .weibull_deviance_core(t) else sat - .weibull_deviance_core(t)
  if (order == 0L) return(out)
  td <- shift[, -1L, drop = FALSE] + v
  # Scale polynomial arguments before forming powers. This keeps the
  # exponential product representable when exp(t) itself underflows.
  for (a in 0:order) {
    d <- td
    d[, 1L] <- d[, 1L] + a
    norm <- pmax(1, abs(d[, 1L]))
    for (j in seq_len(order - a)) norm <- pmax(norm, abs(d[, j]))
    z <- d / norm
    bell <- list(rep(1, length(y)))
    for (b in 0:(order - a)) {
      if (a + b == 0L) next
      if (b > 0L) {
        terms <- lapply(seq_len(b), function(j)
          choose(b - 1L, j - 1L) * (z[, j] / norm^(j - 1L)) * bell[[b - j + 1L]])
        bell[[b + 1L]] <- Reduce(`+`, terms)
      }
      pol <- bell[[b + 1L]]
      log_factor <- t + a * theta + b * log(norm)
      exponential <- (-1)^a * exp(log_factor) * pol
      bad <- which(!is.finite(exponential))
      exponential[bad] <- (-1)^a * sign(pol[bad]) *
        exp(log_factor[bad] + log(abs(pol[bad])))
      exponential[which(pol == 0 & is.finite(t))] <- 0
      linear <- if (a == 0L) td[, b] else if (a == 1L) -k else 0
      derivative <- linear - exponential
      # Retain exp(t) - 1 near saturation in the first derivatives.
      if (a + b == 1L) derivative <- -linear * expm1(t)
      if (deviance) derivative <- -derivative
      else if (a == 0L && b == 1L) derivative <- derivative + 1
      degree <- a + b
      out[, degree * (degree + 1L) / 2L + b + 1L] <- derivative
    }
  }
  out
}

.weibull_parameters <- function(mu, k = NULL) {
  if (is.null(k)) {
    if (!is.matrix(mu) || ncol(mu) != 2L)
      stop("mu must be a two-column matrix containing mean and shape", call. = FALSE)
    k <- mu[, 2L]
    mu <- mu[, 1L]
  }
  if (!is.numeric(mu) || !is.numeric(k) ||
      any(!is.finite(mu) | mu <= 0) || any(!is.finite(k) | k <= 0))
    stop("Weibull mean and shape must be positive and finite", call. = FALSE)
  list(mu = mu, k = rep_len(k, length(mu)))
}

.weibull_probability <- function(x, mu, k = NULL, quantile = FALSE,
                                 lower.tail = TRUE, log.p = FALSE) {
  par <- .weibull_parameters(mu, k)
  if (!length(x) || !length(par$mu)) return(numeric())
  n <- max(length(x), length(par$mu))
  x <- rep_len(x, n)
  mu <- rep_len(par$mu, n)
  k <- rep_len(par$k, n)
  shift <- .weibull_shift(log(k))[, 1L]
  if (quantile) {
    log_hazard <- log(stats::qexp(x, lower.tail = lower.tail, log.p = log.p))
    if (log.p && lower.tail) {
      tiny <- which(x < -36)
      log_hazard[tiny] <- x[tiny]
    }
    q <- exp(log(mu) + (log_hazard - shift) / k)
    return(q)
  }
  t <- k * .weibull_log_ratio(pmax(0, x), mu) + shift
  if (!lower.tail) {
    log_survival <- -exp(t)
    return(if (log.p) log_survival else exp(log_survival))
  }
  prob <- -expm1(-exp(t))
  if (!log.p) return(prob)
  log_prob <- log(prob)
  tiny <- which(t < -36)
  log_prob[tiny] <- t[tiny]
  log_prob
}

.weibull_random <- function(mu, k = NULL) {
  .weibull_probability(stats::runif(NROW(mu)), mu, k, quantile = TRUE)
}

.weibull_log_variance <- function(mu, k) {
  k <- rep_len(k, length(mu))
  u <- 1 / k
  ratio <- numeric(length(mu))
  mid <- which(u > .005 & u < 10000)
  ratio[mid] <- lgamma(1 + 2 * u[mid]) - 2 * lgamma(1 + u[mid])
  large <- which(u >= 10000)
  # Stirling difference, with leading terms combined before subtraction.
  ratio[large] <- 2 * log(2) * u[large] - .5 * (log(pi) - log(k[large])) -
    k[large] / 8 + k[large]^3 / 192
  ratio[is.infinite(u)] <- Inf
  log_cv2 <- ratio
  ordinary <- which(ratio < 36)
  log_cv2[ordinary] <- log(expm1(ratio[ordinary]))
  small <- which(u <= .005)
  j <- 2:length(.weibull_lgamma_series)
  correction <- drop(outer(u[small], j - 2L, `^`) %*%
    (.weibull_lgamma_series[j] * (2^j - 2)))
  # Evaluate expm1(u^2 * correction)/u^2 by its continuous limit at zero.
  z <- u[small]^2 * correction
  relative <- rep(1, length(z))
  nonzero <- z != 0
  relative[nonzero] <- expm1(z[nonzero]) / z[nonzero]
  log_cv2[small] <- log(correction * relative) - 2 * log(k[small])
  2 * log(mu) + log_cv2
}

.weibull_deviance <- function(y, mu, k, wt) {
  wt <- rep_len(wt, length(y))
  keep <- wt > 0
  out <- numeric(length(y))
  k <- rep_len(k, length(y))[keep]
  t <- k * .weibull_log_ratio(y[keep], mu[keep]) + .weibull_shift(log(k))[, 1L]
  out[keep] <- 2 * (wt[keep] * .weibull_deviance_core(t))
  out
}

.weibull_residuals <- function(y, mu, k, wt, type) {
  wt <- rep_len(wt, length(y))
  if (type == "response") return(y - mu)
  k <- rep_len(k, length(y))
  out <- numeric(length(y))
  keep <- which(wt > 0)
  y <- y[keep]; mu <- mu[keep]; k <- k[keep]; w <- wt[keep]
  if (type == "deviance") {
    t <- k * .weibull_log_ratio(y, mu) + .weibull_shift(log(k))[, 1L]
    out[keep] <- sign(t) * sqrt(2 * (w * .weibull_deviance_core(t)))
  } else {
    diff <- y - mu
    out[keep] <- sign(diff) * exp(log(abs(diff)) +
      (log(w) - .weibull_log_variance(mu, k)) / 2)
  }
  out
}

.weibull_null_deviance <- function(y, k, w, offset, intercept) {
  keep <- w > 0
  k <- rep_len(k, length(y))[keep]
  offset <- if (length(offset)) rep_len(offset, length(y))[keep] else numeric(sum(keep))
  y <- y[keep]
  w <- w[keep]
  t0 <- k * (log(y) - offset) + .weibull_shift(log(k))[, 1L]
  b <- 0
  if (intercept) {
    limits <- range(t0 / k)
    if (diff(limits) == 0) return(0)
    log_wk <- log(w) + log(k)
    # The score is monotone. Compare its positive and negative terms on
    # the log scale using the known roots at the individual saturated means.
    log_total <- function(x) {
      m <- max(x)
      if (is.infinite(m)) return(m)
      m + log(sum(exp(x - m)))
    }
    score <- function(b) log_total(log_wk + (t0 - k * b)) - log_total(log_wk)
    b <- stats::uniroot(score, limits, f.lower = 1, f.upper = -1, tol = 1e-10)$root
  }
  sum(2 * (w * .weibull_deviance_core(t0 - k * b)))
}

.weibull_log_link <- function() {
  link <- stats::make.link("log")
  link$linkinv <- link$mu.eta <- exp
  link
}
