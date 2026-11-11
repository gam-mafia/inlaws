# v(p) is the standard Gumbel quantile. Work relative to the location
# quantile so that translating the response does not change the blend width.
.bgev_geometry <- function(t, xi, spec, j) {
  v <- spec$v
  C <- j$constant
  n <- nrow(t)
  log_neg_exprel <- function(x) {
    ans <- C(numeric(n))
    small <- x[, 1L] <= .5
    if (any(small)) {
      z <- x[small, , drop = FALSE]
      ans[small, ] <- j$log(j$polynomial(-z, 1/gamma(2:21)))
    }
    if (any(!small)) {
      z <- x[!small, , drop = FALSE]
      ans[!small, ] <- j$log(-j$expm1(-z)) - j$log(z)
    }
    ans
  }
  log_base <- t - C(rep(log(v[["hi"]] - v[["lo"]]), n)) -
    log_neg_exprel((v[["hi"]] - v[["lo"]])*xi)
  log_increment <- function(lo, hi) {
    d <- hi - lo
    log_base + C(rep(log(d), n)) + (hi - v[["hi"]])*xi + log_neg_exprel(d*xi)
  }
  relative <- function(at) {
    if (at == v[["alpha"]]) return(C(numeric(n)))
    ends <- sort(c(at, v[["alpha"]]))
    sign(at - v[["alpha"]])*j$exp(log_increment(ends[1L], ends[2L]))
  }
  log_width <- log_increment(v[["a"]], v[["b"]])
  list(a = relative(v[["a"]]), b = relative(v[["b"]]),
       log_width = log_width, log_sigma_b = log_base + (v[["b"]] - v[["hi"]])*xi)
}

# log1p(x)/x with its analytic continuation at zero. Factoring the
# likelihood this way retains derivatives in the Gumbel limit.
.bgev_log1prel <- function(x, j) {
  out <- j$constant(numeric(nrow(x)))
  small <- abs(x[, 1L]) <= .125
  if (any(small)) out[small, ] <- j$polynomial(
    x[small, , drop = FALSE], (-1)^(0:24)/(1:25))
  if (any(!small)) {
    z <- x[!small, , drop = FALSE]
    log_z <- .bgev_log1p(z, j)
    # Divide by the base value first to avoid inverse powers of large x.
    base <- z[, 1L]
    out[!small, ] <- j$multiply(log_z/base, j$inverse(z/base))
  }
  out
}

.bgev_log1p <- function(x, j) {
  out <- j$log(j$constant(rep(1, nrow(x))) + x)
  out[, 1L] <- log1p(x[, 1L])
  out
}

# Inside the blend, GEV and Gumbel cumulative hazards are exp(log_f)
# and exp(log_g). The GEV log hazard is convex and lies below the
# Gumbel secant, so the derivative of the beta weight adds density.
.bgev_blend <- function(w, xi, spec, j, density = FALSE) {
  n <- nrow(w)
  C <- j$constant
  M <- j$multiply
  one <- C(rep(1, n))
  d <- spec$d
  x <- d*xi
  small <- which(x[, 1L] <= .5)
  large <- which(x[, 1L] > .5)
  decrement <- relative_gradient <- C(numeric(n))
  if (length(small)) {
    z <- x[small, , drop = FALSE]
    w0 <- w[small, , drop = FALSE]
    e <- j$polynomial(z, 1/gamma(2:21))
    u <- M(w0, j$expm1(z))
    decrement[small, ] <- d*M(M(w0, e), .bgev_log1prel(u, j))
    if (density) relative_gradient[small, ] <- M(e, j$inverse(one[small, , drop = FALSE] + u))
  }
  if (length(large)) {
    z <- x[large, , drop = FALSE]
    w0 <- w[large, , drop = FALSE]
    tail <- one[large, , drop = FALSE] - w0
    denominator <- w0 + M(tail, j$exp(-z))
    inverse_x <- j$inverse(z)
    decrement[large, ] <- d*M(z + j$log(denominator), inverse_x)
    if (density) relative_gradient[large, ] <- M(-j$expm1(-z), M(inverse_x, j$inverse(denominator)))
  }
  log_g <- C(rep(-spec$v[["a"]], n)) - d*w
  log_f <- C(rep(-spec$v[["a"]], n)) - decrement
  delta <- j$expm1(log_f - log_g)

  # Evaluate the beta CDF polynomial from the nearest end, where it
  # vanishes to fifth order. Its derivative is 630*w^4*(1-w)^4.
  left <- w[, 1L] <= .5
  edge <- w
  edge[!left, ] <- (one - w)[!left, , drop = FALSE]
  edge2 <- M(edge, edge)
  edge4 <- M(edge2, edge2)
  edge5 <- M(edge4, edge)
  small_weight <- M(edge5, j$polynomial(edge, c(126, -420, 540, -315, 70)))
  weight <- complement <- small_weight
  weight[!left, ] <- (one - small_weight)[!left, , drop = FALSE]
  complement[left, ] <- (one - small_weight)[left, , drop = FALSE]
  log_rate <- log_g
  if (any(left)) log_rate[left, ] <- log_g[left, , drop = FALSE] +
    .bgev_log1p(M(weight, delta)[left, , drop = FALSE], j)
  if (any(!left)) log_rate[!left, ] <- log_g[!left, , drop = FALSE] +
    j$log((complement + M(weight, j$exp(log_f - log_g)))[!left, , drop = FALSE])
  if (!density) return(log_rate)

  w_product <- M(w, one - w)
  wp2 <- M(w_product, w_product)
  derivative <- 630*M(wp2, wp2)
  relative_f <- M(j$exp(log_f - log_g), relative_gradient)
  gradient <- d*(complement + M(weight, relative_f)) - M(derivative, delta)
  list(log_rate = log_rate, log_gradient = log_g + j$log(gradient))
}

.bgev_eval <- function(y, eta, order, spec) {
  j <- .gen_gamma_jets(order)
  C <- j$constant
  M <- j$multiply
  predictors <- lapply(1:3, function(k) {
    x <- C(eta[, k])
    if (order) x[, which(j$powers[, k] == 1L & rowSums(j$powers) == 1L)] <- 1
    x
  })
  xi <- j$exp(predictors[[3L]])
  g <- .bgev_geometry(predictors[[2L]], xi, spec, j)
  centered <- C(y) - predictors[[1L]]
  left <- which(centered[, 1L] <= g$a[, 1L])
  right <- which(centered[, 1L] >= g$b[, 1L] & centered[, 1L] > g$a[, 1L])
  middle <- setdiff(seq_along(y), c(left, right))
  ell <- C(rep(-Inf, length(y)))
  if (length(left)) {
    w <- M((centered - g$a)[left, , drop = FALSE], j$exp(-g$log_width[left, , drop = FALSE]))
    log_rate <- C(rep(-spec$v[["a"]], length(left))) - spec$d*w
    ell[left, ] <- log_rate - j$exp(log_rate) +
      C(rep(log(spec$d), length(left))) - g$log_width[left, , drop = FALSE]
  }
  if (length(right)) {
    log_sigma <- g$log_sigma_b[right, , drop = FALSE]
    z <- M((centered - g$b)[right, , drop = FALSE], j$exp(-log_sigma))
    x <- M(xi[right, , drop = FALSE], z)
    log_rate <- C(rep(-spec$v[["b"]], length(right))) - M(z, .bgev_log1prel(x, j))
    ell[right, ] <- log_rate - j$exp(log_rate) - log_sigma - .bgev_log1p(x, j)
  }
  if (length(middle)) {
    log_width <- g$log_width[middle, , drop = FALSE]
    w <- M((centered - g$a)[middle, , drop = FALSE], j$exp(-log_width))
    blend <- .bgev_blend(w, xi[middle, , drop = FALSE], spec, j, density = TRUE)
    ell[middle, ] <- blend$log_gradient - j$exp(blend$log_rate) - log_width
  }
  list(ell = ell, powers = j$powers)
}

.bgev_parameters <- function(mu) {
  if (!is.matrix(mu) || !is.numeric(mu) || ncol(mu) != 3L)
    stop("bgevlss parameters must be a numeric matrix with location, spread and shape columns")
  if (any(!is.finite(mu[!is.na(mu)])) ||
      any(mu[, 2L] <= 0, na.rm = TRUE) || any(mu[, 3L] < 0, na.rm = TRUE))
    stop("bgevlss requires finite location, positive spread and nonnegative shape")
}

.bgev_blend_rate <- function(w, xi, spec) {
  # The scalar/vector path also serves the bracketed inverse and moments.
  j <- .gen_gamma_jets(0L)
  .bgev_blend(j$constant(w), j$constant(xi), spec, j)[, 1L]
}

.bgev_probability <- function(x, mu, spec, quantile = FALSE, lower.tail = TRUE, log.p = FALSE) {
  .bgev_parameters(mu)
  n <- nrow(mu)
  if (!is.numeric(x) || !length(x) %in% c(1L, n))
    stop("probabilities or quantiles must be numeric, scalar or one per parameter row")
  x <- rep_len(x, n)
  out <- rep(NA_real_, n)
  ok <- stats::complete.cases(mu) & !is.na(x)
  if (!any(ok)) return(out)
  x <- x[ok]
  loc <- mu[ok, 1L]
  xi <- mu[ok, 3L]
  j <- .gen_gamma_jets(0L)
  g <- lapply(.bgev_geometry(j$constant(log(mu[ok, 2L])), j$constant(xi), spec, j), drop)
  width <- exp(g$log_width)
  va <- spec$v[["a"]]
  vb <- spec$v[["b"]]
  d <- spec$d
  if (quantile) {
    invalid <- if (log.p) x > 0 else x < 0 | x > 1
    if (any(invalid)) warning("NaNs produced")
    x[invalid] <- NaN
    lp <- if (log.p) x else log(x)
    log_rate <- log(-lp)
    if (!lower.tail) {
      # log(-log(1-p)) = log(p) + O(p) for tiny survival probabilities.
      log_rate <- lp
      usual <- which(lp >= -36 & lp <= log(.5))
      high <- which(lp > log(.5))
      log_rate[usual] <- log(-log1p(-exp(lp[usual])))
      log_rate[high] <- log(-log(-expm1(lp[high])))
    }
    ans <- log_rate
    ans[log_rate == Inf & !is.na(log_rate)] <- -Inf
    ans[log_rate == -Inf & !is.na(log_rate)] <- Inf
    left <- which(is.finite(log_rate) & log_rate >= -va)
    right <- which(is.finite(log_rate) & log_rate <= -vb)
    middle <- which(is.finite(log_rate) & log_rate < -va & log_rate > -vb)
    ans[left] <- loc[left] + g$a[left] - width[left]*(log_rate[left] + va)/d
    z <- -log_rate[right] - vb
    xz <- xi[right]*z
    # Evaluate the increment in log space when exp(xz) alone would overflow.
    log_exprel <- numeric(length(xz))
    moderate <- xz > 0 & xz <= .5
    large <- xz > .5 & is.finite(xz)
    log_exprel[is.infinite(xz)] <- Inf
    log_exprel[moderate] <- log(expm1(xz[moderate])/xz[moderate])
    log_exprel[large] <- xz[large] + log(-expm1(-xz[large])) - log(xz[large])
    ans[right] <- loc[right] + g$b[right] +
      exp(g$log_sigma_b[right] + log(z) + log_exprel)
    if (length(middle)) {
      # Solve in the standardized blend, independently of response units.
      lo <- numeric(length(middle))
      hi <- rep(1, length(middle))
      for (iteration in seq_len(50L)) {
        mid <- (lo + hi)/2
        above <- .bgev_blend_rate(mid, xi[middle], spec) > log_rate[middle]
        lo[above] <- mid[above]
        hi[!above] <- mid[!above]
      }
      ans[middle] <- loc[middle] + g$a[middle] + width[middle]*(lo + hi)/2
    }
  } else {
    centered <- x - loc
    log_rate <- ifelse(is.infinite(x), -x, NA_real_)
    left <- which(is.finite(x) & centered <= g$a)
    right <- which(is.finite(x) & centered >= g$b & centered > g$a)
    middle <- which(is.finite(x) & centered > g$a & centered < g$b)
    log_rate[left] <- -va - d*((centered[left] - g$a[left])/width[left])
    z <- centered[right] - g$b[right]
    log_z <- log(z) - g$log_sigma_b[right]
    log_xz <- log(xi[right]) + log_z
    decrement <- exp(log_z)
    moderate <- which(log_xz > -36 & log_xz < 30)
    large <- which(log_xz >= 30)
    decrement[moderate] <- log1p(exp(log_xz[moderate]))/xi[right[moderate]]
    decrement[large] <- (log_xz[large] + log1p(exp(-log_xz[large])))/xi[right[large]]
    log_rate[right] <- -vb - decrement
    log_rate[middle] <- .bgev_blend_rate(
      (centered[middle] - g$a[middle])/width[middle], xi[middle], spec)
    ans <- -exp(log_rate)
    if (!lower.tail) {
      usual <- which(log_rate >= -36)
      tiny <- which(log_rate < -36)
      ans[usual] <- log(-expm1(ans[usual]))
      ans[tiny] <- log_rate[tiny]
    }
    if (!log.p) ans <- exp(ans)
  }
  out[ok] <- ans
  out
}

.bgev_moments <- function(mu, spec) {
  .bgev_parameters(mu)
  n <- nrow(mu)
  m <- v <- rep(NA_real_, n)
  valid <- stats::complete.cases(mu)
  v[valid & mu[, 3L] >= .5] <- Inf
  rows <- which(valid & mu[, 3L] < 1)
  if (!length(rows)) return(list(mean = m, variance = v))

  # Use (Y - b)/width to avoid cancellation from the location and scale.
  # On the left use t = -log(G(Y)) = ta + u, an exponential variable.
  ta <- exp(-spec$v[["a"]])
  tb <- exp(-spec$v[["b"]])
  d <- spec$d
  left <- vapply(1:2, function(k) stats::integrate(function(u)
    (-1 - log1p(u/ta)/d)^k*exp(-u - ta), 0, Inf, rel.tol = 1e-10)$value,
    numeric(1))
  rule <- .gauss_legendre(32L)
  w <- (1 + rule$nodes)/2
  j <- .gen_gamma_jets(0L)
  shapes <- unique(mu[rows, 3L])
  std <- vapply(shapes, function(xi) {
    blend <- .bgev_blend(j$constant(w), j$constant(rep(xi, length(w))), spec, j,
                         density = TRUE)
    mass <- exp(blend$log_gradient[, 1L] - exp(blend$log_rate[, 1L]))*rule$weights/2
    sigma_ratio <- if (xi*d < .Machine$double.eps) 1/d else xi/(-expm1(-xi*d))
    raw <- c(NA_real_, Inf)
    for (k in seq_len(if (xi < .5) 2L else 1L)) {
      # On the right s = log(tb / t). Integrate the leading tail term
      # exactly, leaving an integrand with extra exponential decay even
      # as xi approaches a moment-existence boundary.
      leading <- if (k == 1L) 1/(1 - xi) else 2/((1 - xi)*(1 - 2*xi))
      correction <- stats::integrate(function(s) {
        inc <- s
        regular <- xi*s >= .Machine$double.eps
        inc[regular] <- -expm1(-xi*s[regular])/xi
        inc^k * exp(-(1 - k*xi)*s) * expm1(-tb*exp(-s))
      }, 0, Inf, rel.tol = 1e-10)$value
      raw[k] <- left[k] + sum(mass*(w - 1)^k) +
        tb*sigma_ratio^k*(leading + correction)
    }
    raw
  }, numeric(2))
  raw <- std[, match(mu[rows, 3L], shapes), drop = FALSE]
  g <- .bgev_geometry(j$constant(log(mu[rows, 2L])), j$constant(mu[rows, 3L]), spec, j)
  width <- exp(g$log_width[, 1L])
  m[rows] <- mu[rows, 1L] + (g$b[, 1L] + width*raw[1L, ])
  finite_var <- mu[rows, 3L] < .5
  v[rows[finite_var]] <- exp(2*g$log_width[finite_var, 1L] +
    log(pmax(0, raw[2L, finite_var] - raw[1L, finite_var]^2)))
  list(mean = m, variance = v)
}
