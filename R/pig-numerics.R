# Taylor coefficients in (mu, log(theta)); coefficients include factorials.
.pig_index <- do.call(rbind, lapply(0:4, function(k) cbind(k:0, 0:k)))
.pig_products <- lapply(seq_len(15L), function(k) {
  ij <- .pig_index[k, ]
  a <- which(.pig_index[, 1] <= ij[1] & .pig_index[, 2] <= ij[2])
  b <- vapply(a, function(j) {
    z <- ij - .pig_index[j, ]; sum(z) * (sum(z) + 1) / 2 + z[2] + 1
  }, numeric(1))
  cbind(a, b)
})
.pig_mul <- function(a, b) {
  out <- a * 0
  for (k in seq_len(ncol(a))) {
    p <- .pig_products[[k]]
    out[, k] <- rowSums(a[, p[, 1], drop = FALSE] * b[, p[, 2], drop = FALSE])
  }
  out
}
.pig_const <- function(x, n = length(x)) {
  a <- matrix(0, n, 15L); a[, 1] <- x; a
}
.pig_unary <- function(a, power = NULL) {
  # Expansion around the constant in relative coordinates avoids large
  # intermediate powers in reciprocals and square roots.
  u <- a / a[, 1]; u[, 1] <- 0
  term <- .pig_const(1, nrow(a)); out <- term
  if (is.null(power)) out[, 1] <- log(a[, 1])
  coefficient <- 1
  for (k in 1:4) {
    term <- .pig_mul(term, u)
    coefficient <- if (is.null(power)) (-1)^(k + 1) / k else
      coefficient * (power - k + 1) / k
    out <- out + coefficient * term
  }
  if (!is.null(power)) out <- out * a[, 1]^power
  out
}

# Half-integer Bessel recurrence. With s=sqrt(1+2*mu*theta), a=theta/s,
# h_0=1, h_k=(2*k-1)*a+1/h_{k-1},
# log p_y = -2*mu/(1+s) + y*log(mu/s) - log(y!) + sum_{k=1}^{y-1} log h_k.
# This is finite, positive arithmetic, including as theta approaches zero.
.pig_logpmf <- function(y, mu, t) {
  n <- max(length(y), length(mu), length(t))
  y <- rep_len(y, n); mu <- rep_len(mu, n); t <- rep_len(t, n)
  theta <- exp(t); s <- sqrt(1 + 2 * mu * theta)
  if (any(!is.finite(s))) stop("PIG parameters exceed numerical range")
  out <- -2 * mu / (1 + s)
  pos <- y > 0 & mu > 0
  out[pos] <- out[pos] + y[pos] * (log(mu[pos]) - log(s[pos])) - lgamma(y[pos] + 1)
  out[y > 0 & mu == 0] <- -Inf
  h <- rep(1, n); a <- theta / s
  if (any(pos) && max(y[pos]) > 1) for (k in seq_len(max(y[pos]) - 1)) {
    ii <- which(pos & y > k)
    h[ii] <- (2 * k - 1) * a[ii] + 1 / h[ii]
    out[ii] <- out[ii] + log(h[ii])
  }
  out
}

.pig_lljet <- function(y, mu, t, predictor = FALSE) {
  n <- length(mu)
  m <- .pig_const(mu)
  if (predictor) {
    for (k in 1:4) m[, k * (k + 1) / 2 + 1] <- mu / factorial(k)
  } else m[, 2] <- 1
  th <- .pig_const(exp(t), n)
  for (k in 1:4) th[, k * (k + 1) / 2 + k + 1] <- exp(t) / factorial(k)
  one <- .pig_const(1, n)
  s <- .pig_unary(one + 2 * .pig_mul(m, th), .5)
  a <- .pig_mul(th, .pig_unary(s, -1))
  out <- -2 * .pig_mul(m, .pig_unary(one + s, -1)) +
    y * (.pig_unary(m) - .pig_unary(s))
  out[, 1] <- out[, 1] - lgamma(y + 1)
  h <- one
  if (length(y) && max(y) > 1) for (k in seq_len(max(y) - 1)) {
    ii <- which(y > k)
    hh <- (2 * k - 1) * a[ii, , drop = FALSE] +
      .pig_unary(h[ii, , drop = FALSE], -1)
    h[ii, ] <- hh
    out[ii, ] <- out[ii, , drop = FALSE] + .pig_unary(hh)
  }
  if (any(!is.finite(out))) stop("PIG likelihood derivatives exceed numerical range")
  out
}

.pig_lld <- function(y, mu, t, level = 2L) {
  z <- .pig_lljet(y, mu, t)
  cols <- c(m = 2, mm = 4, t = 3, mt = 5, mmm = 7, mmt = 8,
            mmmm = 11, tt = 6, mtt = 9, mmtt = 13, mmmt = 12)
  lapply(cols, function(k) z[, k] * prod(factorial(.pig_index[k, ])))
}

.pig_score <- function(y, mu, t) {
  theta <- exp(t); s <- sqrt(1 + 2 * mu * theta)
  a <- theta / s; da <- -a * theta / s^2
  out <- -1 / s + y * (1 / mu - theta / s^2)
  h <- 1; dh <- 0
  if (y > 1) for (k in seq_len(y - 1)) {
    dh <- (2 * k - 1) * da - dh / h^2
    h <- (2 * k - 1) * a + 1 / h
    out <- out + dh / h
  }
  out
}
.pig_saturated <- function(y, t, derivatives = TRUE) {
  scalar <- length(t) == 1L
  u <- if (scalar) unique(y) else y
  ts <- rep_len(t, length(u))
  ans <- matrix(0, length(u), 3L)
  ii <- which(u > 0); mus <- numeric(length(ii))
  for (j in seq_along(ii)) {
    yi <- u[ii[j]]; ti <- ts[ii[j]]
    score <- function(logmu) .pig_score(yi, exp(logmu), ti) * exp(logmu)
    lo <- log(yi) - 1; hi <- log(yi) + 1
    while (score(lo) < 0) lo <- lo - 2
    while (score(hi) > 0) {
      hi <- hi + 2
      if (hi > 350) stop("PIG saturated mean exceeds numerical range")
    }
    mus[j] <- exp(stats::uniroot(score, c(lo, hi), tol = 1e-10)$root)
  }
  if (length(ii)) {
    ans[ii, 1] <- .pig_logpmf(u[ii], mus, ts[ii])
    if (derivatives) {
      dd <- .pig_lld(u[ii], mus, ts[ii])
      ans[ii, 2:3] <- cbind(dd$t, dd$tt - dd$mt^2 / dd$mm)
    }
  }
  if (scalar) ans[match(y, u), , drop = FALSE] else ans
}

.pig_recycle <- function(x, mu) {
  if (!is.numeric(x) || !is.numeric(mu)) stop("PIG distribution arguments must be numeric")
  if (!length(x) || !length(mu)) return(list(numeric(), numeric()))
  n <- max(length(x), length(mu))
  if (!all(c(length(x), length(mu)) %in% c(1L, n)))
    stop("PIG distribution arguments must have equal lengths or be scalars")
  if (any(!is.na(mu) & (!is.finite(mu) | mu < 0)))
    stop("PIG distribution means must be finite and nonnegative")
  list(rep_len(x, n), rep_len(mu, n))
}
.pig_probability <- function(x, mu, t, quantile = FALSE) {
  v <- .pig_recycle(x, mu); x <- v[[1]]; mu <- v[[2]]
  if (!length(x)) return(numeric())
  if (!length(t) || !length(t) %in% c(1L, length(x)))
    stop("PIG dispersion must be scalar or match the distribution arguments")
  t <- rep_len(t, length(x))
  if (quantile && any(!is.na(x) & (x < 0 | x > 1))) warning("NaNs produced")
  vapply(seq_along(x), function(i) {
    xi <- x[i]; mi <- mu[i]; ti <- t[i]
    if (is.na(xi) || is.na(mi) || is.na(ti)) return(NA_real_)
    if (quantile) {
      if (xi < 0 || xi > 1) return(NaN)
      if (xi == 0 || mi == 0) return(0)
      if (xi == 1) return(Inf)
      target <- log(xi); last <- 999999L
    } else {
      if (xi < 0) return(-Inf)
      if (mi == 0 || xi == Inf) return(0)
      last <- floor(xi)
      if (last >= 1000000L) stop("PIG CDF summation exceeds one million terms")
    }
    s <- sqrt(1 + 2 * mi * exp(ti))
    if (!is.finite(s)) stop("PIG parameters exceed numerical range")
    a <- exp(ti) / s
    lp <- total <- -2 * mi / (1 + s); h <- 1
    if (quantile && total >= target) return(0)
    if (last > 0) for (k in seq_len(last)) {
      if (k > 1) h <- (2 * k - 3) * a + 1 / h
      lp <- lp + log(mi) - log(s) - log(k) + log(h)
      peak <- max(total, lp)
      total <- peak + log1p(exp(min(total, lp) - peak))
      if (quantile && total >= target) return(as.numeric(k))
    }
    if (quantile) stop("PIG quantile summation failed after one million terms")
    if (!is.finite(total) || total > 1e-8) stop("PIG CDF summation lost numerical accuracy")
    min(0, total)
  }, numeric(1))
}
.pig_random <- function(mu, t) {
  mu <- .pig_recycle(mu, mu)[[2]]
  # Michael--Schucany--Haas inverse Gaussian sampler, with a rationalized
  # small root to avoid subtraction of nearly equal numbers.
  v <- exp(t) * stats::rnorm(length(mu))^2 / 2
  if (any(!is.finite(v))) stop("PIG simulation exceeds numerical range")
  z <- 1 / (1 + v + sqrt(v) * sqrt(2 + v))
  if (any(!is.finite(z) | z <= 0)) stop("PIG simulation exceeds numerical range")
  flip <- stats::runif(length(mu)) > 1 / (1 + z)
  z[flip] <- 1 / z[flip]
  rate <- mu * z
  if (any(!is.na(rate) & (!is.finite(rate) | rate > 2^53)))
    stop("PIG simulation exceeds exact integer range")
  stats::rpois(length(mu), rate)
}
