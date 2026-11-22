## Let F = log((a)_m) and la = log(a). Return r = F - m*la, the scaled
## natural-shape derivatives A = a*F_a, B = a^2*F_aa, and log-shape
## derivatives D2, D3 and D4. A is also the first log-shape derivative of F.
## dA = A - m is evaluated separately for the multinomial limit.
.dirmult_rising <- function(la, m, order = 2L) {
  r <- A <- dA <- B <- D2 <- D3 <- D4 <- numeric(length(m))

  ii <- which(m > 0 & m <= 64)
  if (length(ii)) {
    if (order >= 1L) A[ii] <- 1
    if (order >= 2L) B[ii] <- -1
    for (k in seq_len(max(m[ii]) - 1L)) {
      jj <- ii[m[ii] > k]
      z <- log(k) - la[jj]
      v <- exp(-abs(z))
      r[jj] <- r[jj] + pmax(z, 0) + log1p(v)
      if (order >= 1L) {
        q <- ifelse(z > 0, v, 1) / (1 + v)
        s <- ifelse(z > 0, 1, v) / (1 + v)
        A[jj] <- A[jj] + q
        dA[jj] <- dA[jj] - s
      }
      if (order >= 2L) {
        curvature <- q * s
        B[jj] <- B[jj] - q^2
        D2[jj] <- D2[jj] + curvature
      }
      if (order >= 3L) D3[jj] <- D3[jj] + curvature * (s - q)
      if (order >= 4L) D4[jj] <- D4[jj] + curvature * (1 - 6 * curvature)
    }
  }

  log_m <- log(m)
  log_ratio <- log_m - la
  series <- log_ratio < log(0.05)
  ii <- which(m > 64 & series)
  if (length(ii)) {
    # Expand log(1 + k/a) in k/a. The normalized power sums avoid large
    # powers of the count. Multiply the terms on the log scale so that
    # m * (m/a)^j can remain representable even if (m/a)^j underflows.
    h <- 1 / m[ii]
    averages <- matrix(1, length(ii), 17L)
    terms <- matrix(0, length(ii), order + 1L)
    for (j in 1:16) {
      value <- 1
      for (k in 0:(j - 1L))
        value <- value - choose(j + 1L, k) * averages[, k + 1L] * h^(j - k)
      averages[, j + 1L] <- value / (j + 1L)
      term <- (-1)^(j + 1L) *
        exp(log_m[ii] + log(averages[, j + 1L]) + j * log_ratio[ii]) / j
      for (k in 0:order)
        terms[, k + 1L] <- terms[, k + 1L] + (-j)^k * term
    }
    r[ii] <- terms[, 1L]
    if (order >= 1L) {
      dA[ii] <- terms[, 2L]
      A[ii] <- m[ii] + dA[ii]
    }
    if (order >= 2L) {
      D2[ii] <- terms[, 3L]
      B[ii] <- D2[ii] - A[ii]
    }
    if (order >= 3L) D3[ii] <- terms[, 4L]
    if (order >= 4L) D4[ii] <- terms[, 5L]
  }

  ii <- which(m > 64 & !series & la < log(32))
  if (length(ii)) {
    a <- exp(la[ii])
    # The k = 0 term is exactly la. Shifting it out avoids singular
    # polygamma evaluations even when exp(la) underflows.
    upper <- a + m[ii]
    lower <- a + 1
    log_gamma <- lgamma(upper)
    overflow <- is.infinite(log_gamma)
    if (any(!overflow)) {
      jj <- ii[!overflow]
      r[jj] <- log_gamma[!overflow] - lgamma(lower[!overflow]) - (m[jj] - 1) * la[jj]
    }
    if (any(overflow)) {
      # Combine the leading terms before multiplying by the count.
      # The reciprocal correction to Stirling's formula is negligible here.
      jj <- ii[overflow]
      log_upper <- log(upper[overflow])
      r[jj] <- (m[jj] - 1) * (log_upper - la[jj] - 1) +
        (lower[overflow] - 0.5) * log_upper - lower[overflow] +
        log(2*pi) / 2 - lgamma(lower[overflow])
    }
    if (order >= 1L) {
      scaled <- lapply(seq_len(order), function(k) {
        delta <- psigamma(upper, deriv = k - 1L) - psigamma(lower, deriv = k - 1L)
        sign(delta) * exp(k * la[ii] + log(abs(delta)))
      })
      A[ii] <- 1 + scaled[[1L]]
      dA[ii] <- A[ii] - m[ii]
    }
    if (order >= 2L) {
      B[ii] <- -1 + scaled[[2L]]
      D2[ii] <- scaled[[1L]] + scaled[[2L]]
    }
    if (order >= 3L) D3[ii] <- scaled[[1L]] + 3 * scaled[[2L]] + scaled[[3L]]
    if (order >= 4L) D4[ii] <- scaled[[1L]] + 7 * scaled[[2L]] +
      6 * scaled[[3L]] + scaled[[4L]]
  }

  ii <- which(m > 64 & !series & la >= log(32))
  if (length(ii)) {
    # Difference of Stirling expansions for log Gamma(a+m) and log Gamma(a).
    # Express the result and derivatives in ratios. Neither a+m nor powers
    # of a need to fit in a double.
    ratio <- exp(log_ratio[ii])
    L <- log1p(ratio)
    s <- -expm1(-L)
    t <- exp(-L)
    power <- seq.int(1L, 13L, by = 2L)
    coefficient <- c(1/12, -1/360, 1/1260, -1/1680, 1/1188, -691/360360, 1/156)
    inverse_powers <- outer(exp(-la[ii]), power, `^`)
    correction <- matrix(0, length(ii), order + 1L)
    for (k in 0:order) {
      difference <- expm1(-outer(L, power + k))
      correction[, k + 1L] <- drop((inverse_powers * difference) %*% coefficient)
      coefficient <- -coefficient * (power + k)
    }
    r[ii] <- m[ii] * (L / ratio + L - 1) - L / 2 + correction[, 1L]
    if (order >= 1L) {
      A[ii] <- m[ii] * (L / ratio) + s / 2 + correction[, 2L]
      dA[ii] <- m[ii] * (L / ratio - 1) + s / 2 + correction[, 2L]
    }
    if (order >= 2L) {
      B[ii] <- -m[ii] * (s / ratio) - s * (2 - s) / 2 + correction[, 3L]
      D2[ii] <- m[ii] * ((L - s) / ratio) - s * t / 2 + rowSums(correction[, 2:3, drop = FALSE])
    }
    if (order >= 3L) D3[ii] <- m[ii] * ((L - s - s^2) / ratio) + s * t * (1 - 2 * s) / 2 +
      drop(correction[, 2:4, drop = FALSE] %*% c(1, 3, 1))
    if (order >= 4L) D4[ii] <- m[ii] * ((L - s - 2 * s^3) / ratio) - s * t * (1 - 6 * s * t) / 2 +
      drop(correction[, 2:5, drop = FALSE] %*% c(1, 7, 6, 1))
  }

  out <- list(r = r)
  if (order >= 1L) out$A <- A
  if (order >= 2L) {
    out$B <- B
    out$D2 <- D2
  }
  if (order >= 3L) out$D3 <- D3
  if (order >= 4L) out$D4 <- D4
  if (order >= 1L) out$dA <- dA
  out
}
