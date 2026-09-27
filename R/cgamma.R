# Small second-order forward derivatives with respect to log dispersion.
# Columns contain the value, first derivative and second derivative (not
# Taylor coefficients). These helpers are private to the censored gamma family.
.cg_jet <- function(x, d1 = 0, d2 = 0) {
  cbind(x, rep_len(d1, length(x)), rep_len(d2, length(x)))
}

.cg_mul <- function(a, b) {
  cbind(a[, 1] * b[, 1],
        a[, 2] * b[, 1] + a[, 1] * b[, 2],
        a[, 3] * b[, 1] + 2 * a[, 2] * b[, 2] + a[, 1] * b[, 3])
}

.cg_inv <- function(a) {
  u <- a[, 2] / a[, 1]
  cbind(1 / a[, 1], -u / a[, 1],
        (2 * u^2 - a[, 3] / a[, 1]) / a[, 1])
}

.cg_log <- function(a) {
  u <- a[, 2] / a[, 1]
  cbind(log(a[, 1]), u, a[, 3] / a[, 1] - u^2)
}

.cg_exp <- function(a) {
  v <- exp(a[, 1])
  cbind(v, v * a[, 2], v * (a[, 3] + a[, 2]^2))
}

# log(exp(a) - exp(b)), where a >= b. Used on the smaller tail.
.cg_logdiff <- function(a, b) {
  v <- a[, 1] + log(-expm1(b[, 1] - a[, 1]))
  r <- exp(b[, 1] - v)
  d <- a[, 2] - b[, 2]
  cbind(v, a[, 2] + r * d,
        a[, 3] + r * (a[, 3] - b[, 3]) - r * (1 + r) * d^2)
}

# Log of x times the unit-rate gamma density, with x = k * bound / mu
# and k = exp(-theta). The derivatives keep bound and mu fixed.
.cg_boundary <- function(x, k) {
  a <- log(x) - digamma(k)
  s <- x - k * (a + 1)
  cbind(stats::dgamma(x, shape = k, log = TRUE) + log(x), s,
        -s + k - k^2 * trigamma(k))
}

# Differentiate the convergent incomplete-gamma series (lower tail) or
# continued fraction (upper tail). Values themselves come from R's pgamma.
# No finite differences of the likelihood are used in model fitting.
.cg_tail <- function(x, k, lower.tail = TRUE, deriv = TRUE) {
  k <- rep_len(k, length(x))
  value <- stats::pgamma(x, shape = k, lower.tail = lower.tail, log.p = TRUE)
  out <- .cg_jet(value)
  finite <- which(x > 0 & is.finite(x))
  if (!deriv || !length(finite)) return(out)
  for (series in c(TRUE, FALSE)) {
    ii <- finite[(x[finite] < k[finite] + 1) == series]
    if (!length(ii)) next
    xx <- x[ii]
    kk <- k[ii]
    X <- .cg_jet(xx, -xx, xx)
    K <- .cg_jet(kk, -kk, kk)
    one <- .cg_jet(rep(1, length(ii)))
    pref <- .cg_boundary(xx, kk)
    converged <- FALSE
    if (series) {
      term <- total <- one
      for (j in seq_len(10000L)) {
        term <- .cg_mul(term, .cg_mul(X, .cg_inv(K + j * one)))
        total <- total + term
        if (all(abs(term) < 2e-14 * (1 + abs(total)))) {
          converged <- TRUE
          break
        }
      }
      ans <- pref - .cg_log(K) + .cg_log(total)
    } else {
      b <- X + one - K
      c <- .cg_jet(rep(1e300, length(ii)))
      d <- .cg_inv(b)
      h <- d
      for (j in seq_len(10000L)) {
        a <- j * (K - j * one)
        b <- b + 2 * one
        d <- .cg_inv(b + .cg_mul(a, d))
        c <- b + .cg_mul(a, .cg_inv(c))
        delta <- .cg_mul(c, d)
        h <- .cg_mul(h, delta)
        if (all(abs(delta - one) < 2e-14)) {
          converged <- TRUE
          break
        }
      }
      ans <- pref + .cg_log(h)
    }
    if (!converged || any(!is.finite(ans))) {
      stop("cgamma: incomplete-gamma derivatives failed to converge")
    }
    # Use pgamma's accurate value, especially before complementing a tail.
    ans[, 1] <- stats::pgamma(xx, kk, lower.tail = series, log.p = TRUE)
    if (series != lower.tail) ans <- .cg_logdiff(0 * one, ans)
    out[ii, ] <- ans
  }
  out[, 1] <- value
  out
}

.cg_response <- function(y) {
  if (is.matrix(y)) {
    if (ncol(y) != 2L) stop("cgamma requires a vector or a two-column response")
    censor <- y[, 2]
    y <- y[, 1]
  } else {
    censor <- attr(y, "censor")
    if (is.null(censor)) censor <- y
  }
  if (!is.numeric(y) || !is.numeric(censor) || length(censor) != length(y) ||
      anyNA(y) || anyNA(censor) || any(!is.finite(y))) {
    stop("cgamma requires finite response values and non-missing censoring bounds")
  }
  exact <- censor == y
  left <- censor == -Inf
  right <- censor == Inf
  interval <- !(exact | left | right)
  if (any(y[exact | left] <= 0) || any(y[right | interval] < 0) ||
      any(censor[interval] <= y[interval])) {
    stop(paste("cgamma requires positive exact observations and left thresholds,",
               "non-negative lower bounds, and ordered interval bounds"))
  }
  lo <- hi <- y
  lo[left] <- 0
  hi[right] <- Inf
  hi[interval] <- censor[interval]
  list(y = y, censor = censor, exact = exact, lo = lo, hi = hi)
}

.cg_narrow <- function(lo, hi, mu, k) {
  is.finite(hi) & lo > 0 & (hi - lo) < 1e-4 * lo & k * (hi - lo) / mu < 1
}

# Narrow intervals suffer cancellation even with log-CDF subtraction, notably
# in derivatives of boundary differences. Integrate the density and its
# analytic derivatives together using a short Gauss-Legendre rule instead.
.cg_quadrature <- function(lo, hi, mu, k) {
  rule <- .gauss_legendre(16)
  y <- lo + (hi - lo) * (rule$nodes + 1) / 2
  x <- k * y / mu
  X <- .cg_jet(x, -x, x)
  K <- .cg_jet(rep(k, length(x)), -k, k)
  one <- .cg_jet(rep(1, length(x)))
  lf <- .cg_boundary(x, k)
  lf[, 1] <- lf[, 1] - log(y)
  centre <- max(lf[, 1])
  lf[, 1] <- lf[, 1] - centre
  f <- rule$weights * .cg_exp(lf)
  total <- matrix(colSums(f), 1L)
  logp <- .cg_log(total)
  logp[, 1] <- logp[, 1] + centre + log((hi - lo) / 2)
  p <- list(X - K)
  p[[2]] <- .cg_mul(X - K - one, p[[1]]) - X
  p[[3]] <- .cg_mul(X - K - 2 * one, p[[2]]) -
    .cg_mul(X, 2 * X - 2 * K - 2 * one)
  p[[4]] <- .cg_mul(X - K - 3 * one, p[[3]]) -
    .cg_mul(X, 3 * .cg_mul(X, X) - 6 * .cg_mul(K + 2 * one, X) +
              3 * .cg_mul(K + one, K + 2 * one))
  ratios <- lapply(seq_len(4L), function(j) {
    .cg_mul(matrix(colSums(.cg_mul(f, p[[j]])), 1L), .cg_inv(total)) / mu^j
  })
  list(logp = logp, ratios = ratios)
}

.cg_probability <- function(lo, hi, mu, k, deriv = TRUE) {
  x0 <- k * lo / mu
  x1 <- k * hi / mu
  # Choose the tail in which subtracting two probabilities loses least precision.
  lower <- stats::pgamma(x0, k) < 0.5
  out <- matrix(0, length(lo), 3L)
  narrow <- .cg_narrow(lo, hi, mu, k)
  for (lt in c(TRUE, FALSE)) {
    ii <- which(lower == lt & !narrow)
    if (!length(ii)) next
    a <- .cg_tail(if (lt) x1[ii] else x0[ii], k, lt, deriv)
    b <- .cg_tail(if (lt) x0[ii] else x1[ii], k, lt, deriv)
    out[ii, ] <- .cg_logdiff(a, b)
  }
  for (i in which(narrow)) out[i, ] <- .cg_quadrature(lo[i], hi[i], mu[i], k)$logp
  out
}

.cg_loglik <- function(y, mu, theta, deriv = TRUE) {
  r <- .cg_response(y)
  mu <- rep_len(mu, length(r$y))
  k <- exp(-theta)
  out <- matrix(0, length(mu), 3L)
  ii <- which(r$exact)
  if (length(ii)) {
    t <- r$y[ii] / mu[ii]
    s <- -k * (log(k) + 1 - digamma(k) + log(t) - t)
    out[ii, ] <- cbind(stats::dgamma(r$y[ii], k, scale = mu[ii] / k, log = TRUE),
                       s, -s + k - k^2 * trigamma(k))
  }
  ii <- which(!r$exact)
  if (length(ii)) out[ii, ] <- .cg_probability(r$lo[ii], r$hi[ii], mu[ii], k, deriv)
  out
}

# Saturated means for finite intervals are logarithmic means of their bounds:
# equality of the two boundary densities gives mu = (u-l)/log(u/l).
# For one-sided intervals the supremum probability is one.
.cg_saturated <- function(y, theta, deriv = TRUE) {
  r <- .cg_response(y)
  out <- matrix(0, length(r$y), 3L)
  ii <- which(r$exact)
  if (length(ii)) out[ii, ] <- .cg_loglik(r$y[ii], r$y[ii], theta, deriv)
  ii <- which(!r$exact & r$lo > 0 & is.finite(r$hi))
  if (length(ii)) {
    mu <- (r$hi[ii] - r$lo[ii]) / log1p((r$hi[ii] - r$lo[ii]) / r$lo[ii])
    out[ii, ] <- .cg_probability(r$lo[ii], r$hi[ii], mu, exp(-theta), deriv)
  }
  out
}

.cg_log_mean <- function(r) {
  a <- r[[1]]; b <- r[[2]]; c <- r[[3]]
  aa <- .cg_mul(a, a)
  list(a, b - aa, c - 3 * .cg_mul(a, b) + 2 * .cg_mul(aa, a),
       r[[4]] - 4 * .cg_mul(a, c) - 3 * .cg_mul(b, b) +
         12 * .cg_mul(aa, b) - 6 * .cg_mul(aa, aa))
}

# In the far upper tail, fourth derivatives formed from probability ratios
# cancel terms of order x^4 to recover a result of order x. Instead differentiate
# Q = G/x * S, S = 1 + (k-1)/x + (k-1)(k-2)/x^2 + ... .
# Used only for x > max(100, 5*k), where the required terms decrease rapidly.
# The expansion is asymptotic: stop on small derivative terms, not after a
# fixed number of terms, and fail explicitly if that tolerance is not reached.
.cg_farupper <- function(bound, mu, k) {
  x <- k * bound / mu
  X <- .cg_jet(x, -x, x)
  K <- .cg_jet(k, -k, k)
  one <- .cg_jet(1)
  term <- one
  sums <- c(list(one), rep(list(.cg_jet(0)), 4L))
  converged <- FALSE
  for (j in seq_len(1000L)) {
    term <- .cg_mul(term, .cg_mul(K - j * one, .cg_inv(X)))
    small <- TRUE
    for (r in 0:4) {
      factor <- if (r == 0) 1 else if (j < r) 0 else prod(j - 0:(r - 1)) / mu^r
      increment <- factor * term
      sums[[r + 1L]] <- sums[[r + 1L]] + increment
      small <- small && all(abs(increment) < 2e-14 * (1 + abs(sums[[r + 1L]])))
    }
    if (!all(is.finite(term))) break
    if (j >= 4L && isTRUE(small)) {
      converged <- TRUE
      break
    }
  }
  if (!converged) stop("cgamma: upper-tail derivative expansion failed to converge")
  inv <- .cg_inv(sums[[1]])
  ratios <- lapply(sums[-1], function(s) .cg_mul(s, inv))
  out <- .cg_log_mean(ratios)
  for (r in 1:4) {
    out[[r]] <- out[[r]] + (-1)^r * factorial(r - 1) *
      (K - one - r * X) / mu^r
  }
  list(logp = .cg_boundary(x, k) - .cg_log(X) + .cg_log(sums[[1]]), mean = out)
}

.cg_farinterval <- function(lo, hi, mu, k) {
  lower <- .cg_farupper(lo, mu, k)
  if (!is.finite(hi)) return(lower$mean)
  upper <- .cg_farupper(hi, mu, k)
  delta <- Map(`-`, upper$mean, lower$mean)
  e <- .cg_exp(upper$logp - lower$logp)
  denominator <- .cg_inv(.cg_jet(1) - e)
  a <- delta[[1]]; b <- delta[[2]]; c <- delta[[3]]
  aa <- .cg_mul(a, a)
  bell <- list(a, b + aa, c + 3 * .cg_mul(a, b) + .cg_mul(aa, a),
               delta[[4]] + 4 * .cg_mul(a, c) + 3 * .cg_mul(b, b) +
                 6 * .cg_mul(aa, b) + .cg_mul(aa, aa))
  ratios <- lapply(bell, function(b) -.cg_mul(.cg_mul(e, b), denominator))
  Map(`+`, lower$mean, .cg_log_mean(ratios))
}

.cg_mean_derivs <- function(y, mu, theta, level) {
  r <- .cg_response(y)
  k <- exp(-theta)
  n <- length(mu)
  order <- level + 2L
  ans <- lapply(seq_len(order), function(j) matrix(0, n, 3L))
  ii <- which(r$exact)
  if (length(ii)) for (j in seq_len(order)) {
    v <- (-1)^j * factorial(j - 1) * k * (1 - j * r$y[ii] / mu[ii]) / mu[ii]^j
    ans[[j]][ii, ] <- cbind(v, -v, v)
  }
  narrow <- .cg_narrow(r$lo, r$hi, mu, k)
  far <- !r$exact & !narrow & k * r$lo / mu > max(100, 5 * k)
  for (i in which(far)) {
    dd <- .cg_farinterval(r$lo[i], r$hi[i], mu[i], k)
    for (j in seq_len(order)) ans[[j]][i, ] <- dd[[j]]
  }
  ii <- which(!r$exact & !far)
  if (length(ii)) {
    logp <- .cg_probability(r$lo[ii], r$hi[ii], mu[ii], k, level > 0)
    ratios <- lapply(seq_len(order), function(j) matrix(0, length(ii), 3L))
    for (side in c("lo", "hi")) {
      bound <- r[[side]][ii]
      jj <- which(bound > 0 & is.finite(bound))
      if (!length(jj)) next
      x <- k * bound[jj] / mu[ii[jj]]
      X <- .cg_jet(x, -x, x)
      K <- .cg_jet(rep(k, length(jj)), -k, k)
      one <- .cg_jet(rep(1, length(jj)))
      g <- .cg_exp(.cg_boundary(x, k) - logp[jj, , drop = FALSE])
      # Polynomials multiplying the boundary density in d^j P / d mu^j.
      b <- list(one, X - K - one)
      if (order >= 3L) b[[3]] <- .cg_mul(X - K - 2 * one, b[[2]]) - X
      if (order >= 4L) {
        b[[4]] <- .cg_mul(X - K - 3 * one, b[[3]]) -
          .cg_mul(X, 2 * X - 2 * K - 4 * one)
      }
      for (j in seq_len(order)) {
        ratios[[j]][jj, ] <- ratios[[j]][jj, ] +
          (if (side == "lo") 1 else -1) * .cg_mul(g, b[[j]]) / mu[ii[jj]]^j
      }
    }
    for (j in which(.cg_narrow(r$lo[ii], r$hi[ii], mu[ii], k))) {
      q <- .cg_quadrature(r$lo[ii[j]], r$hi[ii[j]], mu[ii[j]], k)
      for (a in seq_len(order)) ratios[[a]][j, ] <- q$ratios[[a]]
    }
    a <- ratios[[1]]
    b <- ratios[[2]]
    aa <- .cg_mul(a, a)
    ans[[1]][ii, ] <- a
    ans[[2]][ii, ] <- b - aa
    if (order >= 3L) {
      c <- ratios[[3]]
      ans[[3]][ii, ] <- c - 3 * .cg_mul(a, b) + 2 * .cg_mul(aa, a)
    }
    if (order >= 4L) {
      ans[[4]][ii, ] <- ratios[[4]] - 4 * .cg_mul(a, c) -
        3 * .cg_mul(b, b) + 12 * .cg_mul(aa, b) - 6 * .cg_mul(aa, aa)
    }
  }
  ans
}

#' Censored gamma distribution for generalized additive models
#'
#' A gamma family with a log link and a single global dispersion parameter,
#' supporting exact, left-censored, right-censored and interval-censored data.
#'
#' @param theta Dispersion parameter, the reciprocal of gamma shape. `NULL` or
#'   zero estimates dispersion starting at one; a positive value fixes it;
#'   a negative value estimates it starting at `abs(theta)`.
#' @param link Only `"log"` is supported.
#'
#' @details
#' The conditional mean is `mu = exp(eta)` and variance is `theta * mu^2`.
#' The gamma shape is `1 / theta` and scale is `mu * theta`. Dispersion is
#' estimated as a family parameter on the log scale. The additional mgcv scale
#' parameter is fixed at one; `family$getTheta(TRUE)` returns dispersion.
#'
#' An ordinary vector response specifies exact observations. For censoring use
#' a two-column matrix, following [mgcv::cnorm]: equal entries denote an exact
#' observation; a second entry of `-Inf` or `Inf` denotes left or right censoring
#' at the first entry; otherwise the entries are ordered lower and upper bounds.
#' Exact observations and left thresholds must be positive. Lower bounds may
#' be zero. Censoring is assumed non-informative conditional on the covariates.
#'
#' Prior weights multiply log-likelihood contributions (replication weights);
#' they do not change individual gamma dispersions. Zero weights are allowed.
#'
#' Response predictions are ordinary plug-in estimates of the conditional
#' arithmetic mean, `exp(eta)`, with delta-method standard errors. They do not
#' include an epsilon correction for nonlinear functions of uncertain model
#' coefficients or integrate over new random effects. Unlike a lognormal model,
#' no residual-variance adjustment is needed to obtain the conditional mean.
#'
#' The family includes `rd(mu, wt = 1, scale = 1)`,
#' `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)`, and
#' `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)`.
#' These evaluate the latent, uncensored gamma distribution using the supplied
#' response-scale mean and the family's current dispersion. `rd` draws one
#' response per element of `mu`. `wt` and `scale` are accepted for compatibility
#' with mgcv: replication weights do not alter the response distribution, and
#' dispersion comes from `getTheta(TRUE)`, not the external scale argument.
#' `qf` and `cdf` support upper tails and log probabilities; the CDF's `logp`
#' spelling follows the mgcv CDF convention.
#'
#' These functions do not apply censoring, draw model parameters, or condition
#' on an observation's censoring event. For a randomized quantile residual of
#' an observation censored to `[l, u]`, draw uniformly between `cdf(l, mu)` and
#' `cdf(u, mu)`, then apply [stats::qnorm()]. Use bounds zero or infinity for
#' one-sided censoring. For exact observations use `qnorm(cdf(y, mu))`.
#'
#' @return An object of class `extended.family` and `family` for [mgcv::gam]
#'   or [mgcv::bam].
#' @seealso [mgcv::cnorm], [mgcv::cpois], [stats::Gamma]
#' @export
#' @examples
#' set.seed(12)
#' dat <- data.frame(x = runif(200))
#' mu <- exp(0.5 + sin(2 * pi * dat$x))
#' y <- rgamma(200, shape = 3, scale = mu / 3)
#' dat$yc <- cbind(y, y)
#' censored <- y < 0.4
#' dat$yc[censored, 1] <- 0.4
#' dat$yc[censored, 2] <- -Inf
#' fit <- mgcv::gam(yc ~ s(x, k = 8), data = dat,
#'                  family = cgamma(), method = "REML")
#' fit$family$getTheta(TRUE) # estimated dispersion; shape is its reciprocal
#' predict(fit, type = "response")
cgamma <- function(theta = NULL, link = "log") {
  if (!identical(link, "log")) stop("cgamma supports only the log link")
  if (!is.null(theta) && (!is.numeric(theta) || length(theta) != 1L ||
                         !is.finite(theta))) {
    stop("theta must be NULL or a single finite numeric dispersion")
  }
  estimate <- is.null(theta) || theta <= 0
  initial <- if (is.null(theta) || theta == 0) 0 else log(abs(theta))
  state <- new.env(parent = emptyenv())
  state$theta <- initial
  linkobj <- stats::make.link("log")
  getTheta <- function(trans = FALSE) if (trans) exp(state$theta) else state$theta
  putTheta <- function(theta) state$theta <- theta
  rd <- function(mu, wt = 1, scale = 1) {
    phi <- getTheta(TRUE)
    stats::rgamma(length(mu), shape = 1 / phi, scale = mu * phi)
  }
  qf <- function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    phi <- getTheta(TRUE)
    stats::qgamma(p, shape = 1 / phi, scale = mu * phi,
                  lower.tail = lower.tail, log.p = log.p)
  }
  cdf <- function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE) {
    phi <- getTheta(TRUE)
    stats::pgamma(q, shape = 1 / phi, scale = mu * phi,
                  lower.tail = lower.tail, log.p = logp)
  }
  dev.resids <- function(y, mu, wt, theta = NULL) {
    if (is.null(theta)) theta <- getTheta()
    d <- 2 * wt * (.cg_saturated(y, theta, FALSE)[, 1] -
                    .cg_loglik(y, mu, theta, FALSE)[, 1])
    d[wt == 0] <- 0
    pmax(d, 0)
  }
  Dd <- function(y, mu, theta, wt, level = 0) {
    wt <- rep_len(wt, length(mu))
    d <- .cg_mean_derivs(y, mu, theta, level)
    d <- lapply(d, function(z) {
      z <- -2 * wt * z
      z[wt == 0, ] <- 0
      z
    })
    curvature <- d[[2]][, 1] + d[[1]][, 1] / mu
    # Censored rows use positive observed curvature on the log-link scale.
    # Exact rows use Gamma's expected information, also used for mgcv's EDF.
    # bam uses EDmu2 without the score term in the link conversion.
    out <- list(Dmu = d[[1]][, 1], Dmu2 = d[[2]][, 1],
                EDmu2 = pmax(curvature, 0))
    exact <- .cg_response(y)$exact
    out$EDmu2[exact] <- 2 * wt[exact] * exp(-theta) / mu[exact]^2
    if (level > 0) {
      dt <- 2 * wt * (.cg_saturated(y, theta) - .cg_loglik(y, mu, theta))
      dt[wt == 0, ] <- 0
      out$Dth <- dt[, 2]
      out$Dmuth <- d[[1]][, 2]
      out$Dmu3 <- d[[3]][, 1]
      out$Dmu2th <- d[[2]][, 2]
      out$EDmu2th <- ifelse(curvature > 0, d[[2]][, 2] + d[[1]][, 2] / mu, 0)
      out$EDmu2th[exact] <- -out$EDmu2[exact]
    }
    if (level > 1) {
      out$Dmu4 <- d[[4]][, 1]
      out$Dth2 <- dt[, 3]
      out$Dmuth2 <- d[[1]][, 3]
      out$Dmu2th2 <- d[[2]][, 3]
      out$Dmu3th <- d[[3]][, 2]
    }
    out
  }
  ls <- function(y, w, theta, scale) {
    z <- w * .cg_saturated(y, theta)
    z[w == 0, ] <- 0
    list(ls = sum(z[, 1]), lsth1 = sum(z[, 2]),
         LSTH1 = z[, 2, drop = FALSE], lsth2 = sum(z[, 3]))
  }
  aic <- function(y, mu, theta = NULL, wt, dev) {
    if (is.null(theta)) theta <- getTheta()
    z <- wt * .cg_loglik(y, mu, theta, FALSE)[, 1]
    z[wt == 0] <- 0
    -2 * sum(z)
  }
  initialize <- expression({
    rr <- family$response(y)
    y <- rr$y
    attr(y, "censor") <- rr$censor
    if (any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
      stop("cgamma requires finite non-negative weights with at least one positive weight")
    }
    mustart <- y
    left <- rr$censor == -Inf
    mustart[left] <- y[left] / 2
    interval <- !rr$exact & is.finite(rr$hi)
    mustart[interval] <- (rr$lo[interval] + rr$hi[interval]) / 2
    positive <- mustart[mustart > 0]
    fallback <- if (length(positive)) stats::median(positive) else 1
    mustart <- pmax(mustart, fallback * 0.01)
  })
  subsety <- function(y, ind) {
    if (is.matrix(y)) return(y[ind, , drop = FALSE])
    censor <- attr(y, "censor")
    y <- y[ind]
    if (!is.null(censor)) attr(y, "censor") <- censor[ind]
    y
  }
  postproc <- function(family, y, prior.weights, fitted, linear.predictors,
                       offset, intercept) {
    find_null_deviance <- utils::getFromNamespace("find.null.dev", "mgcv")
    list(null.deviance = find_null_deviance(family, y, eta = linear.predictors,
                                           offset = offset, weights = prior.weights),
         family = paste0("Censored Gamma(dispersion = ",
                         signif(family$getTheta(TRUE), 4), ")"))
  }
  structure(list(family = "Censored Gamma", link = "log",
                 linkfun = linkobj$linkfun, linkinv = linkobj$linkinv,
                 mu.eta = linkobj$mu.eta, valideta = linkobj$valideta,
                 validmu = function(mu) all(is.finite(mu) & mu > 0),
                 variance = function(mu) mu^2 * getTheta(TRUE),
                 dev.resids = dev.resids, Dd = Dd, ls = ls, aic = aic,
                 initialize = initialize, response = .cg_response,
                 subsety = subsety, postproc = postproc,
                 getTheta = getTheta, putTheta = putTheta,
                 rd = rd, qf = qf, cdf = cdf,
                 n.theta = as.integer(estimate), ini.theta = initial,
                 scale = 1, no.r.sq = TRUE),
            class = c("extended.family", "family"))
}
