# Student-t log density as analytic Taylor coefficients. The existing jet
# algebra is distribution-independent and stores only convolution tables.
.scatlss_eval <- function(y, eta, order = 0L, min.df = 3) {
  j <- .gen_gamma_jets(order); C <- j$constant; M <- j$multiply
  n <- nrow(eta)
  predictors <- lapply(1:3, function(k) {
    x <- C(eta[, k])
    if (order) {
      p <- integer(3); p[k] <- 1L
      x[, which(rowSums(abs(sweep(j$powers, 2, p))) == 0)] <- 1
    }
    x
  })
  u <- predictors[[1]]; t <- predictors[[2]]; v <- predictors[[3]]
  one <- C(rep(1, n))
  log1p_jet <- function(x) {
    z <- j$log(C(rep(1, nrow(x))) + x)
    z[, 1] <- log1p(x[, 1]); z
  }
  # Reciprocal df avoids overflowing powers of nu in fourth derivatives.
  h <- C(numeric(n)); positive <- eta[, 3] > 0
  if (any(positive)) {
    e <- j$exp(-v[positive, , drop = FALSE])
    h[positive, ] <- M(e, j$inverse(C(rep(1, sum(positive))) + min.df*e))
  }
  if (any(!positive)) h[!positive, ] <- j$inverse(
    C(rep(min.df, sum(!positive))) + j$exp(v[!positive, , drop = FALSE]))
  normalizer <- C(numeric(n)); asym <- h[, 1] <= 1/40
  if (any(asym)) {
    # log Gamma((nu+1)/2) - log Gamma(nu/2) - log(pi*nu)/2.
    # Stirling's ratio series through nu^-15 retains finite-df information.
    coefficients <- c(-.5*log(2*pi), -1/4, 0, 1/24, 0, -1/20, 0,
      17/112, 0, -31/36, 0, 691/88, 0, -5461/52, 0, 929569/480)
    normalizer[asym, ] <- j$polynomial(h[asym, , drop = FALSE], coefficients)
  }
  if (any(!asym)) {
    df <- j$inverse(h[!asym, , drop = FALSE])
    normalizer[!asym, ] <- j$lgamma((df + C(rep(1, nrow(df))))/2) -
      j$lgamma(df/2) - .5*j$log(df) - C(rep(.5*log(pi), nrow(df)))
  }
  z <- M(C(y) - u, j$exp(-t)); z2 <- M(z, z); a <- M(z2, h)
  # log1p(a)/a is analytic at zero, including extremely large finite nu.
  f <- C(numeric(n)); small <- a[, 1] <= .25
  if (any(small)) f[small, ] <- j$polynomial(a[small, , drop = FALSE], (-1)^(0:40)/(1:41))
  if (any(!small)) {
    aa <- a[!small, , drop = FALSE]
    f[!small, ] <- M(log1p_jet(aa), j$inverse(aa))
  }
  sat <- normalizer - t
  list(ell = sat - .5*M(M(one + h, z2), f), sat = sat, powers = j$powers)
}

.scatlss_parameters <- function(mu) {
  if (!is.matrix(mu) || !is.numeric(mu) || ncol(mu) != 3L)
    stop("scatlss parameters must be a numeric matrix with mu, sigma and nu columns")
  if (any(!is.na(mu[, 1]) & !is.finite(mu[, 1])) ||
      any(!is.na(mu[, 2]) & (!is.finite(mu[, 2]) | mu[, 2] <= 0)) ||
      any(!is.na(mu[, 3]) & (!is.finite(mu[, 3]) | mu[, 3] <= 2)))
    stop("scatlss requires finite location, positive finite scale and finite df > 2")
  mu
}

.scatlss_probability <- function(x, mu, quantile = FALSE, lower.tail = TRUE, log.p = FALSE) {
  .scatlss_parameters(mu)
  n <- nrow(mu)
  if (!is.numeric(x) || !length(x) %in% c(1L, n))
    stop("probabilities or quantiles must be numeric, scalar or one per parameter row")
  x <- rep_len(x, n)
  if (quantile) mu[, 1] + mu[, 2]*stats::qt(x, mu[, 3], lower.tail = lower.tail, log.p = log.p) else
    stats::pt((x-mu[, 1])/mu[, 2], mu[, 3], lower.tail = lower.tail, log.p = log.p)
}

.scatlss_variance <- function(mu) mu[, 2]^2 * (1 + 2/(mu[, 3]-2))
