# Log sinc(pi*x) and its derivatives. The small-x series avoids subtracting
# singular terms in the derivatives; upper-bound evaluation uses 1-s instead.
.llog_sinc <- function(x, order = 0L) {
  ans <- vector("list", order + 1L)
  small <- x < 0.05
  ans[[1]] <- log(sinpi(x)/(pi*x))
  if (order >= 1L) ans[[2]] <- pi*cospi(x)/sinpi(x) - 1/x
  if (order >= 2L) ans[[3]] <- -pi^2/sinpi(x)^2 + 1/x^2
  if (order >= 3L) ans[[4]] <- 2*pi^3*cospi(x)/sinpi(x)^3 - 2/x^3
  if (order >= 4L) ans[[5]] <- -2*pi^4*(1 + 2*cospi(x)^2)/sinpi(x)^4 + 6/x^4
  # log(sin(a)/a) = -sum zeta(2k)/(k*pi^(2k)) * a^(2k).
  coefficients <- c(-pi^2/6, -pi^4/180, -pi^6/2835,
                    -pi^8/37800, -pi^10/467775, -691*pi^12/3831077250)
  if (any(small)) for (d in 0:order) {
    v <- numeric(sum(small))
    for (k in seq_along(coefficients)) if (2*k >= d)
      v <- v + coefficients[k]*factorial(2*k)/factorial(2*k-d)*x[small]^(2*k-d)
    ans[[d+1L]][small] <- v
  }
  ans
}

.llog_correction <- function(s) {
  upper <- s > 0.5
  ans <- .llog_sinc(pmin(s, 1-s))[[1]]
  ans[upper] <- ans[upper] + log1p(-s[upper]) - log(s[upper])
  ans
}

.llog_response <- function(y) {
  if (!is.numeric(y) || !is.null(dim(y)) || any(!is.finite(y) | y <= 0))
    stop("log-logistic families require a finite positive numeric response vector")
  y
}

.llog_variance <- function(mu, s) {
  s <- rep_len(s, length(mu))
  v <- rep(Inf, length(mu))
  ok <- s < 0.5
  a <- pi*s[ok]
  f <- tan(a)/a - 1
  small <- a < 0.01
  f[small] <- a[small]^2*(1/3 + a[small]^2*(2/15 + a[small]^2*(17/315 + a[small]^2*62/2835)))
  v[ok] <- mu[ok]^2*f
  v
}

# Shared analytic Taylor evaluator, in (mu, logit(s)) or (log(mu), logit(s)).
# Coefficients are derivatives divided by factorials, not finite differences.
.llog_eval <- function(y, mu, theta, order = 0L, predictor = FALSE) {
  n <- length(y)
  mu <- rep_len(mu, n); theta <- rep_len(theta, n)
  j <- .clnormal_jets(order, order)
  C <- j$constant; M <- j$multiply
  u <- C(mu); t <- C(theta)
  if (order) {
    u[, which(j$powers[, 1] == 1 & j$powers[, 2] == 0)] <- 1
    t[, which(j$powers[, 1] == 0 & j$powers[, 2] == 1)] <- 1
  }
  softplus <- function(x) {
    v <- x[, 1]; p <- stats::plogis(v); q <- stats::plogis(-v)
    ds <- list(pmax(v, 0) + log1p(exp(-abs(v))), p, p*q,
               p*q*(q-p), p*q*(1-6*p*q))
    j$compose(x, lapply(0:order, function(k) ds[[k+1L]]/factorial(k)))
  }
  logS <- -softplus(-t)
  S <- j$exp(logS); U <- j$exp(-softplus(t))
  # sin(pi*s) = sin(pi*(1-s)); use the small argument at either boundary.
  upper <- theta > 0
  V <- S; V[upper, ] <- U[upper, ]
  ds <- .llog_sinc(V[, 1], order)
  correction <- j$compose(V, lapply(0:order, function(k) ds[[k+1L]]/factorial(k)))
  correction[upper, ] <- correction[upper, , drop = FALSE] - t[upper, , drop = FALSE]
  loc <- (if (predictor) u else j$log(u)) + correction
  z <- M(C(log(y))-loc, j$exp(-logS))
  p <- stats::plogis(z[, 1]); q <- stats::plogis(-z[, 1])
  ds <- list(stats::dlogis(z[, 1], log = TRUE), q-p, -2*p*q,
             -2*p*q*(q-p), -2*p*q*(1-6*p*q))
  density <- j$compose(z, lapply(0:order, function(k) ds[[k+1L]]/factorial(k)))
  list(ell = density - logS - C(log(y)),
       sat = -logS - C(log(4) + log(y)), powers = j$powers)
}

.llog_parameters <- function(mu, s = NULL) {
  if (is.null(s)) {
    if (!is.matrix(mu) || ncol(mu) != 2L)
      stop("mu must be a two-column matrix containing mean and scale")
    s <- mu[, 2]; mu <- mu[, 1]
  }
  if (!is.numeric(mu) || !is.numeric(s) || any(!is.finite(mu) | mu <= 0) ||
      any(!is.finite(s) | s <= 0 | s >= 1))
    stop("mean must be positive and finite and scale must lie strictly between zero and one")
  list(location = log(mu) + .llog_correction(s), s = s, n = length(mu))
}
.llog_random <- function(mu, s = NULL) {
  p <- .llog_parameters(mu, s)
  exp(stats::rlogis(p$n, p$location, p$s))
}
.llog_quantile <- function(p, mu, s = NULL, lower.tail = TRUE, log.p = FALSE) {
  par <- .llog_parameters(mu, s)
  exp(stats::qlogis(p, par$location, par$s, lower.tail = lower.tail, log.p = log.p))
}
.llog_cdf <- function(q, mu, s = NULL, logp = FALSE, lower.tail = TRUE) {
  par <- .llog_parameters(mu, s)
  # Negative arguments are outside the support, just like zero.
  stats::plogis(log(pmax(q, 0)), par$location, par$s, log.p = logp, lower.tail = lower.tail)
}
