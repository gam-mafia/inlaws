# Taylor coefficients are derivatives divided by multi-factorials. Cache only
# convolution tables, never observation-dependent values.
.gen_gamma_jets <- local({
  cache <- new.env(parent = emptyenv())
  function(order) {
    key <- as.character(order)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
    powers <- as.matrix(expand.grid(a = 0:order, b = 0:order, c = 0:order))
    powers <- powers[rowSums(powers) <= order, , drop = FALSE]
    nc <- nrow(powers)
    pairs <- lapply(seq_len(nc), function(k) which(
      outer(powers[, 1], powers[, 1], `+`) == powers[k, 1] &
      outer(powers[, 2], powers[, 2], `+`) == powers[k, 2] &
      outer(powers[, 3], powers[, 3], `+`) == powers[k, 3], arr.ind = TRUE))
    constant <- function(x) cbind(x, matrix(0, length(x), nc - 1L))
    multiply <- function(x, y) {
      out <- matrix(0, nrow(x), nc)
      for (k in seq_len(nc)) {
        z <- pairs[[k]]
        out[, k] <- rowSums(x[, z[, 1], drop = FALSE] * y[, z[, 2], drop = FALSE])
      }
      out
    }
    compose <- function(x, ds) {
      dx <- x; dx[, 1] <- 0
      out <- constant(ds[[1]]); p <- constant(rep(1, nrow(x)))
      if (order) for (k in seq_len(order)) {
        p <- multiply(p, dx); out <- out + p * ds[[k + 1L]]
      }
      out
    }
    exponential <- function(x) compose(x, lapply(0:order, function(k) exp(x[, 1])/factorial(k)))
    expminusone <- function(x) {
      out <- exponential(x); out[, 1] <- expm1(x[, 1]); out
    }
    logarithm <- function(x) compose(x, c(list(log(x[, 1])),
      lapply(seq_len(order), function(k) (-1)^(k - 1)/(k*x[, 1]^k))))
    inverse <- function(x) compose(x, lapply(0:order, function(k) (-1)^k/x[, 1]^(k + 1L)))
    lg <- function(x) compose(x, c(list(lgamma(x[, 1])),
      lapply(seq_len(order), function(k) psigamma(x[, 1], k - 1L)/factorial(k))))
    polynomial <- function(x, coefficients) {
      ds <- lapply(0:order, function(d) {
        ans <- numeric(nrow(x))
        for (k in rev(seq.int(d, length(coefficients) - 1L)))
          ans <- ans*x[, 1] + coefficients[k + 1L]*choose(k, d)
        ans
      })
      compose(x, ds)
    }
    out <- list(powers = powers, constant = constant, multiply = multiply,
      compose = compose, exp = exponential, expm1 = expminusone, log = logarithm,
      inverse = inverse, lgamma = lg, polynomial = polynomial)
    assign(key, out, cache); out
  }
})

# Bernoulli(2j) / {2j (2j-1)} in the Stirling correction for log Gamma.
.gen_gamma_stirling <- c(1/12, -1/360, 1/1260, -1/1680, 1/1188,
                         -691/360360, 1/156, -3617/122400)

# Stable mean correction, analytic at Q=0. The Stirling difference requires
# BOTH gamma arguments >= 20; otherwise evaluate the gamma functions directly.
.gen_gamma_core <- function(eta, order = 0L) {
  j <- .gen_gamma_jets(order); C <- j$constant; M <- j$multiply
  n <- nrow(eta)
  u <- C(eta[, 1]); t <- C(eta[, 2]); v <- C(eta[, 3])
  if (order) for (k in 1:3) {
    p <- integer(3); p[k] <- 1L
    ix <- which(rowSums(abs(sweep(j$powers, 2, p))) == 0L)
    if (k == 1L) u[, ix] <- 1 else if (k == 2L) t[, ix] <- 1 else v[, ix] <- 1
  }
  one <- C(rep(1, n))
  h <- j$expm1(v); Q <- M(h, j$exp(-t)); Q2 <- M(Q, Q)
  k0 <- 1/Q[, 1]^2
  correction <- normalizer <- C(numeric(n))
  asym <- k0 >= 20 & k0*exp(eta[, 3]) >= 20
  if (any(asym)) {
    ix <- which(asym); hh <- h[ix, , drop = FALSE]; vv <- v[ix, , drop = FALSE]
    # F(h) = {(1+h) log(1+h)-h}/h^2, F(0) = 1/2.
    f <- C(numeric(length(ix))); small <- abs(hh[, 1]) <= .25
    if (any(small)) f[small, ] <- j$polynomial(hh[small, , drop = FALSE],
      (-1)^(0:60)/((1:61)*(2:62)))
    if (any(!small)) {
      a <- hh[!small, , drop = FALSE]; b <- vv[!small, , drop = FALSE]
      f[!small, ] <- M(M(one[ix[!small], , drop = FALSE] + a, b) - a,
                       M(j$inverse(a), j$inverse(a)))
    }
    cc <- M(j$exp(2*t[ix, , drop = FALSE]), f) - .5*vv
    qq <- Q2[ix, , drop = FALSE]; q4 <- M(qq, qq); power <- qq
    aa <- C(rep(-.5*log(2*pi), length(ix)))
    for (s in seq_along(.gen_gamma_stirling)) {
      b <- .gen_gamma_stirling[s]
      cc <- cc + b*M(power, j$expm1(-(2*s - 1)*vv))
      aa <- aa - b*power; power <- M(power, q4)
    }
    correction[ix, ] <- cc; normalizer[ix, ] <- aa
  }
  if (any(!asym)) {
    ix <- which(!asym); kk <- j$inverse(Q2[ix, , drop = FALSE])
    vv <- v[ix, , drop = FALSE]; hh <- h[ix, , drop = FALSE]
    kr <- M(kk, j$exp(vv))
    correction[ix, ] <- j$lgamma(kr) - j$lgamma(kk) - M(M(kk, hh), j$log(kk))
  }
  # Normalize using Q alone, including when the second gamma argument is small.
  regular <- k0 < 20; ix <- which(regular)
  if (length(ix)) {
    kk <- j$inverse(Q2[ix, , drop = FALSE]); q <- Q[ix, , drop = FALSE]
    q[q[, 1] < 0, ] <- -q[q[, 1] < 0, , drop = FALSE]
    normalizer[ix, ] <- j$log(q) - j$lgamma(kk) + M(kk, j$log(kk)) - kk
  }
  ix <- which(!asym & !regular)
  if (length(ix)) {
    qq <- Q2[ix, , drop = FALSE]; q4 <- M(qq, qq); power <- qq
    aa <- C(rep(-.5*log(2*pi), length(ix)))
    for (b in .gen_gamma_stirling) { aa <- aa - b*power; power <- M(power, q4) }
    normalizer[ix, ] <- aa
  }
  list(location = u - correction, sigma = j$exp(t), Q = Q,
       normalizer = normalizer, correction = correction, t = t, j = j)
}

.gen_gamma_log_density <- function(log_y, r) {
  j <- r$j; M <- j$multiply; C <- j$constant
  z <- M(C(log_y) - r$location, j$exp(-r$t)); x <- M(r$Q, z)
  # E(x) = {exp(x)-1-x}/x^2 is analytic at zero.
  e <- C(numeric(length(log_y))); small <- abs(x[, 1]) <= 1
  if (any(small)) e[small, ] <- j$polynomial(x[small, , drop = FALSE], 1/gamma(3:35))
  if (any(!small)) {
    xx <- x[!small, , drop = FALSE]; inv <- j$inverse(xx)
    e[!small, ] <- M(j$expm1(xx) - xx, M(inv, inv))
  }
  r$normalizer - r$t - M(M(z, z), e)
}

.gen_gamma_eval <- function(y, eta, order = 0L) {
  r <- .gen_gamma_core(eta, order); C <- r$j$constant
  list(ell = .gen_gamma_log_density(log(y), r) - C(log(y)),
       sat = r$normalizer - r$t - C(log(y)), powers = r$j$powers,
       location = r$location[, 1], Q = r$Q[, 1])
}

.gen_gamma_parameters <- function(mu) {
  if (!is.matrix(mu) || !is.numeric(mu) || ncol(mu) != 3L)
    stop("gen_gammalss parameters must be a numeric matrix with mean, sigma and r columns")
  if (any(!is.na(mu) & (!is.finite(mu) | mu <= 0)))
    stop("gen_gammalss requires positive finite mean, sigma and r")
  mu
}

.gen_gamma_distribution <- function(mu) {
  .gen_gamma_parameters(mu)
  out <- list(location = rep(NA_real_, nrow(mu)), sigma = mu[, 2],
              Q = rep(NA_real_, nrow(mu)), correction = rep(NA_real_, nrow(mu)))
  ok <- stats::complete.cases(mu)
  if (any(ok)) {
    r <- .gen_gamma_core(log(mu[ok, , drop = FALSE]))
    out$location[ok] <- r$location[, 1]; out$Q[ok] <- r$Q[, 1]
    out$correction[ok] <- r$correction[, 1]
  }
  out
}

.gen_gamma_probability <- function(x, mu, quantile = FALSE, lower.tail = TRUE, log.p = FALSE) {
  a <- .gen_gamma_distribution(mu); n <- nrow(mu)
  if (!is.numeric(x) || !length(x) %in% c(1L, n))
    stop("probabilities or quantiles must be numeric, scalar or one per parameter row")
  x <- rep_len(x, n); out <- rep(NA_real_, n)
  ok <- stats::complete.cases(mu) & !is.na(x)
  normal <- ok & a$Q == 0
  # Second-order expansion of the log-gamma law avoids gamma-argument
  # rounding at huge shape. It retains the signed shape correction.
  tiny <- ok & abs(a$Q) <= 1e-5 & a$Q != 0
  if (any(tiny)) {
    ii <- which(tiny); Q <- a$Q[ii]
    if (quantile) {
      zz <- stats::qnorm(x[ii], lower.tail = lower.tail, log.p = log.p)
      w <- zz; use <- is.finite(zz) & abs(Q)*(1+abs(zz))^3 <= .01
      tiny[ii[is.finite(zz) & !use]] <- FALSE
      w[use] <- zz[use] - Q[use]*(zz[use]^2+2)/6 + Q[use]^2*(zz[use]^3+5*zz[use])/36
      out[ii] <- exp(a$location[ii] + a$sigma[ii]*w)
    } else {
      zz <- (log(pmax(0, x[ii])) - a$location[ii])/a$sigma[ii]
      lp <- stats::pnorm(zz, lower.tail = lower.tail, log.p = TRUE)
      use <- is.finite(zz) & abs(Q)*(1+abs(zz))^3 <= .01
      tiny[ii[is.finite(zz) & !use]] <- FALSE
      d <- Q[use]*(zz[use]^2+2)/6 - Q[use]^2*(zz[use]^5+2*zz[use]^3+6*zz[use])/72
      if (!lower.tail) d <- -d
      lp[use] <- lp[use] + log1p(d*exp(stats::dnorm(zz[use], log = TRUE)-lp[use]))
      out[ii] <- if (log.p) lp else exp(lp)
    }
  }
  if (any(normal)) out[normal] <- if (quantile)
    stats::qlnorm(x[normal], a$location[normal], a$sigma[normal], lower.tail, log.p) else
    stats::plnorm(pmax(0, x[normal]), a$location[normal], a$sigma[normal], lower.tail, log.p)
  for (positive in c(TRUE, FALSE)) {
    ii <- which(ok & !normal & !tiny & (a$Q > 0) == positive)
    if (!length(ii)) next
    Q <- a$Q[ii]; k <- 1/Q^2
    tail <- if (positive) lower.tail else !lower.tail
    if (quantile) {
      g <- stats::qgamma(x[ii], shape = k, lower.tail = tail, log.p = log.p)
      w <- (log(g)-log(k))/Q
      under <- !is.na(g) & g == 0 & is.finite(x[ii])
      if (any(under)) {
        lp <- if (log.p) x[ii][under] else log(x[ii][under])
        if (!tail) lp <- if (log.p) log(-expm1(x[ii][under])) else log1p(-x[ii][under])
        # When the gamma quantile underflows, its leading small-x CDF
        # gives an accurate log quantile; the transformed Y may be finite.
        logg <- (lp + lgamma(k[under]+1))/k[under]
        w[under] <- (logg-log(k[under]))/Q[under]
      }
      close <- is.finite(g) & abs(g/k-1) < .5
      w[close] <- log1p((g[close]-k[close])/k[close])/Q[close]
      out[ii] <- exp(a$location[ii] + a$sigma[ii]*w)
    } else {
      z <- (log(pmax(0, x[ii])) - a$location[ii])/a$sigma[ii]
      logg <- log(k) + Q*z
      out[ii] <- stats::pgamma(exp(logg), shape = k, lower.tail = tail, log.p = log.p)
      under <- logg < log(.Machine$double.xmin)
      if (any(under)) {
        lp <- pmin(0, k[under]*logg[under] - lgamma(k[under]+1))
        if (!tail) lp <- log(-expm1(lp))
        out[ii[under]] <- if (log.p) lp else exp(lp)
      }
    }
  }
  out
}

.gen_gamma_variance <- function(mu) {
  a <- .gen_gamma_distribution(mu); out <- rep(NA_real_, nrow(mu))
  ok <- stats::complete.cases(mu); finite <- ok & mu[, 3] > .5
  out[ok & !finite] <- Inf
  if (any(finite)) {
    ix <- which(finite)
    eta <- cbind(0, log(2*mu[ix, 2]), log(2*mu[ix, 3]-1))
    c2 <- .gen_gamma_core(eta)$correction[, 1]
    out[ix] <- mu[ix, 1]^2*expm1(c2 - 2*a$correction[ix])
  }
  out
}
