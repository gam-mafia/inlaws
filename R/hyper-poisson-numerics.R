# Hyper-Poisson numerics. Taylor coefficients include factorial divisors.
.hp_control <- function(control = list()) {
  defaults <- list(sum_tol = 1e-10, mean_tol = 1e-10,
                   max_terms = 100000L, max_iter = 100L)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) ||
       any(!nzchar(names(control))) || anyDuplicated(names(control)) ||
       any(!names(control) %in% names(defaults)))))
    stop("invalid or unknown hyper-Poisson control", call. = FALSE)
  defaults[names(control)] <- control
  for (nm in names(defaults)) {
    x <- defaults[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 ||
        (nm %in% c("sum_tol", "mean_tol") && x >= 1) ||
        (nm %in% c("max_terms", "max_iter") &&
         (x != floor(x) || x > .Machine$integer.max)))
      stop("invalid hyper-Poisson control: ", nm, call. = FALSE)
  }
  defaults
}

.hp_index <- function(order) {
  do.call(rbind, lapply(0:order, function(d) cbind(d:0, 0:d)))
}

# The rth log(sigma) derivatives of -log(sigma + k), divided by r!.
# Ratios avoid cancellation of polygamma differences at large sigma.
.hp_step <- function(k, sigma, order) {
  u <- sigma / (sigma + k)
  v <- k / (sigma + k)
  b <- matrix(0, length(k), order)
  if (order >= 1L) b[, 1] <- -u
  if (order >= 2L) b[, 2] <- -u * v / 2
  if (order >= 3L) b[, 3] <- -u * v * (1 - 2*u) / 6
  if (order >= 4L) b[, 4] <- -u * v * (1 - 6*u + 6*u^2) / 24
  if (order >= 5L) b[, 5] <- -u * v * (1 - 14*u + 36*u^2 - 24*u^3) / 120
  b
}

# b_y - b_base, where b_y = log Gamma(sigma) - log Gamma(sigma+y).
.hp_bdiff <- function(y, base, sigma, order = 0L) {
  n <- abs(y - base)
  if (n == 0) return(numeric(order + 1L))
  lo <- min(y, base); hi <- max(y, base)
  sign <- if (y > base) 1 else -1
  if (n <= 1000) {
    k <- seq.int(lo, hi - 1)
    return(sign * c(-sum(log(sigma + k)),
                    if (order) colSums(.hp_step(k, sigma, order))))
  }
  x <- sigma + lo
  # Stable log-Gamma ratio using Stirling corrections when x is large.
  correction <- function(x) (1/12 + (-1/360 + (1/1260 - 1/(1680*x^2))/x^2)/x^2)/x
  lr <- if (x >= 10) n*log(x) + (x+n-.5)*log1p(n/x) - n +
    correction(x+n) - correction(x) else lgamma(x+n) - lgamma(x)
  out <- c(-lr, numeric(order))
  if (order) {
    # Stirling numbers of the second kind convert sigma derivatives to log(sigma).
    S <- matrix(0, 5, 5); S[1, 1] <- 1
    for (r in 2:5) for (k in 1:r)
      S[r, k] <- k*S[r-1, k] + if (k > 1) S[r-1, k-1] else 0
    for (r in seq_len(order)) for (k in seq_len(r)) out[r+1L] <- out[r+1L] +
      S[r,k] * sigma^k * (psigamma(x, k-1L) - psigamma(x+n, k-1L)) / factorial(r)
  }
  sign * out
}

.hp_series <- function(a, sigma, control) {
  lambda <- exp(a)
  mode <- max(0, floor(lambda - sigma) + 1)
  if (!is.finite(lambda) || !is.finite(mode) || mode > 2^50)
    stop("hyper-Poisson mode exceeds numerical range", call. = FALSE)
  # Below sigma the successive ratios are bounded by a geometric law.
  # Its spread avoids a huge initial window when sigma is large but mu small.
  ratio <- lambda/sigma
  spread <- if (ratio < 1) sqrt(ratio)/(1-ratio) else sqrt(lambda+1)
  width <- min(control$max_terms, max(16, ceiling(12*spread)))
  lo <- max(0, mode - floor((width - 1)/2)); hi <- lo + width - 1
  repeat {
    if (hi - lo + 1 > control$max_terms)
      stop("hyper-Poisson summation failed: max_terms exhausted", call. = FALSE)
    j <- seq.int(lo, hi); at <- mode - lo + 1L
    # Relative log weights from the recurrence; no large log-Gamma subtraction.
    lw <- c(0, cumsum(a - log(sigma + utils::head(j, -1L))))
    lw <- lw - lw[at]
    w <- exp(lw); d <- abs(j - mode)
    ok <- TRUE
    for (r in 0:5) {
      du <- hi + 1 - mode
      lr <- a - log(sigma + hi + 1) + r*log1p(1/du)
      ub <- if (lr < 0) exp(utils::tail(lw, 1) + a - log(sigma + hi) +
        r*log(du) - log(-expm1(lr))) else Inf
      lb <- 0
      if (lo > 0) {
        dl <- mode - lo + 1
        lr <- if (lo == 1) -Inf else log(sigma + lo - 2) - a + r*log1p(1/dl)
        lb <- if (lr < 0) exp(lw[1] + log(sigma + lo - 1) - a +
          r*log(dl) - log(-expm1(lr))) else Inf
      }
      if (!is.finite(ub + lb) || ub + lb > control$sum_tol*sum(w*(1+d)^r)) {
        ok <- FALSE; break
      }
    }
    if (ok) break
    extra <- max(8, ceiling(length(j)/2))
    lo <- max(0, lo-extra); hi <- hi+extra
  }
  logsw <- log1p(sum(w[-at]))
  p <- w/exp(logsw); mean <- sum(p*j)
  list(a = a, sigma = sigma, mode = mode, y = j, p = p,
       logsw = logsw, mean = mean, variance = sum(p*(j-mean)^2))
}

.hp_mean_one <- function(mu, sigma, control) {
  if (!is.finite(mu) || mu < 0 || !is.finite(sigma) || sigma <= 0)
    stop("hyper-Poisson requires finite nonnegative means and positive dispersions", call. = FALSE)
  if (mu == 0) return(list(a = -Inf, sigma = sigma, mode = 0, y = 0,
                          p = 1, logsw = 0, mean = 0, variance = 0))
  # Small means have lambda ~ mu*sigma; large sigma approaches geometric.
  a <- log(mu) + log(mu + sigma) - log1p(mu)
  lower <- -Inf; upper <- Inf
  for (iter in seq_len(control$max_iter)) {
    z <- .hp_series(a, sigma, control)
    err <- z$mean - mu
    if (abs(err) <= control$mean_tol*mu) return(z)
    if (err < 0) lower <- a else upper <- a
    step <- err/z$variance; next_a <- a-step
    if (!is.finite(next_a) || next_a <= lower || next_a >= upper || abs(step) > 5) {
      next_a <- if (is.finite(lower) && is.finite(upper)) (lower+upper)/2 else
        a + if (err < 0) 1 else -1
    }
    a <- next_a
  }
  stop("hyper-Poisson mean inversion failed: max_iter exhausted", call. = FALSE)
}

.hp_moments <- function(mu, sigma, control) {
  if (length(mu) != length(sigma)) stop("incompatible hyper-Poisson parameter lengths")
  keys <- paste(sprintf("%.17g", mu), sprintf("%.17g", sigma))
  ii <- match(keys, unique(keys)); first <- which(!duplicated(keys))
  z <- lapply(first, function(i) .hp_mean_one(mu[i], sigma[i], control))
  z[ii]
}

.hp_logmass <- function(y, z) {
  if (is.na(y)) return(NA_real_)
  if (!is.finite(y) || y < 0 || y != floor(y)) return(-Inf)
  if (z$mean == 0) return(if (y == 0) 0 else -Inf)
  (y-z$mode)*z$a + .hp_bdiff(y, z$mode, z$sigma) - z$logsw
}

# Cached convolution plans for bivariate Taylor polynomials.
.hp_spec <- local({
  cache <- vector("list", 6L)
  function(order) {
    if (!is.null(cache[[order+1L]])) return(cache[[order+1L]])
    index <- .hp_index(order)
    pairs <- lapply(seq_len(nrow(index)), function(k) {
      i <- which(index[,1] <= index[k,1] & index[,2] <= index[k,2])
      ij <- matrix(index[k,], length(i), 2, byrow = TRUE) - index[i,,drop=FALSE]
      cbind(i, rowSums(ij)*(rowSums(ij)+1)/2 + ij[,2]+1L)
    })
    cache[[order+1L]] <<- list(index = index, pairs = pairs)
    cache[[order+1L]]
  }
})
.hp_mul <- function(x, y, order = 4L) {
  pairs <- .hp_spec(order)$pairs
  out <- x*0
  for (k in seq_along(pairs)) {
    p <- pairs[[k]]
    out[,k] <- rowSums(x[,p[,1],drop=FALSE]*y[,p[,2],drop=FALSE])
  }
  out
}
.hp_compose <- function(coef, x, t) {
  one <- x*0; one[,1] <- 1
  xp <- tp <- vector("list", 5); xp[[1]] <- tp[[1]] <- one
  for (k in 2:5) {
    xp[[k]] <- .hp_mul(xp[[k-1L]], x)
    tp[[k]] <- .hp_mul(tp[[k-1L]], t)
  }
  out <- one*0; index <- .hp_spec(4L)$index
  for (k in seq_len(ncol(coef))) out <- out +
    coef[,k]*.hp_mul(xp[[index[k,1]+1L]], tp[[index[k,2]+1L]])
  out
}

# Normalizer variation G(da,dt), after removing mean*da and the
# mode's log-Gamma jet. Centering makes the remaining MGF start at degree two.
.hp_normalizer <- function(z) {
  index <- .hp_spec(5L)$index
  m <- matrix(0, length(z), nrow(index)); means <- matrix(0, length(z), 5)
  for (i in seq_along(z)) {
    s <- z[[i]]; j <- s$y; at <- s$mode-j[1]+1L
    b <- rbind(numeric(5), apply(.hp_step(utils::head(j,-1L), s$sigma, 5L), 2, cumsum))
    b <- sweep(b, 2, b[at,], "-")
    means[i,] <- colSums(b*s$p)
    b <- sweep(b, 2, means[i,], "-")
    e <- matrix(0, length(j), 6L); e[,1] <- 1
    for (q in 1:5) for (k in 1:q) e[,q+1L] <- e[,q+1L] + k*b[,k]*e[,q-k+1L]/q
    x <- j-s$mean
    for (k in 4:nrow(index)) {
      p <- index[k,1]; q <- index[k,2]
      m[i,k] <- sum(s$p*x^p*e[,q+1L])/factorial(p)
    }
  }
  out <- m-.hp_mul(m,m,5L)/2
  for (q in 1:5) out[,q*(q+1)/2+q+1] <- out[,q*(q+1)/2+q+1]+means[,q]
  out
}

.hp_lljet <- function(y, mu, sigma, z, order = 4L) {
  n <- length(y)
  out <- matrix(0,n,15L)
  out[,1] <- vapply(seq_len(n), function(i) .hp_logmass(y[i],z[[i]]), numeric(1))
  if (!order) return(out[,1,drop=FALSE])
  G <- .hp_normalizer(z)
  Ga <- G[,1:15,drop=FALSE]*0
  index <- .hp_spec(4L)$index
  for (k in 2:15) {
    ij <- index[k,]+c(1,0)
    col <- sum(ij)*(sum(ij)+1)/2+ij[2]+1L
    Ga[,k] <- G[,col]*ij[1]
  }
  dm <- dt <- da <- matrix(0,n,15L); dt[,3] <- 1
  for (k in 1:4) dm[,k*(k+1)/2+1] <- mu/factorial(k)
  variance <- vapply(z, `[[`, numeric(1), "variance")
  if (any(!is.finite(variance) | variance <= 0))
    stop("hyper-Poisson variance is numerically zero or nonfinite", call. = FALSE)
  for (k in 1:4) da <- da+(dm-.hp_compose(Ga,da,dt))/variance
  bc <- t(vapply(seq_len(n), function(i)
    .hp_bdiff(y[i],z[[i]]$mode,sigma[i],4L)[-1L], numeric(4)))
  ans <- (y-mu)*da-.hp_compose(G[,1:15,drop=FALSE],da,dt)
  for (k in 1:4) ans[,k*(k+1)/2+k+1] <- ans[,k*(k+1)/2+k+1]+bc[,k]
  ans[,1] <- out[,1]
  if (any(!is.finite(ans))) stop("hyper-Poisson derivatives exceed numerical range", call. = FALSE)
  ans[,seq_len((order+1L)*(order+2L)/2),drop=FALSE]
}
