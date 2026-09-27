#' Generalized Pareto models for threshold excesses
#'
#' A two-predictor family for [mgcv::gam()], modelling the scale and shape of
#' nonnegative excesses above a fixed threshold.
#'
#' @param link Two link names, in scale, shape order. Scale must use `"log"`.
#'   Shape uses `"logshift"` (the default) or `"identity"`.
#' @details
#' Supply a list of two formulae: the response and scale predictor first, then
#' the shape predictor. Select observations above a threshold `u` and supply
#' `excess = y - u`. Location is fixed at zero on this excess scale. Thresholds
#' may vary between observations but are treated as known. This family does not
#' estimate exceedance frequencies or threshold uncertainty; return periods
#' require additional information about exceedance frequency.
#'
#' The density is \eqn{\sigma^{-1}(1+\xi z/\sigma)^{-1/\xi-1}}, for
#' \eqn{z\geq 0}, \eqn{\sigma>0} and \eqn{1+\xi z/\sigma>0}. At
#' \eqn{\xi=0} this is the exponential density with scale \eqn{\sigma}.
#' The shape link `"logshift"` is \eqn{\log(\xi+0.5)}, restricting shape to
#' \eqn{\xi>-0.5}. The identity link removes this restriction, allowing
#' nonregular likelihood behaviour and potentially unstable estimation.
#' Both links still require every positively weighted observation to lie in
#' the distribution's support. Prior weights multiply log likelihoods.
#'
#' Fitted values and `predict(..., type = "response")` are matrices with
#' columns for scale and shape (in that order), both on their natural scales. They are not
#' conditional means. Link predictions are log scale and linked shape.
#' Only `gam()` with ordinary model matrices is supported, not `bam()`.
#'
#' The family includes `rd(mu, wt = NULL, scale = NULL)`,
#' `qf(p, mu, wt = NULL, scale = NULL)` and
#' `cdf(q, mu, wt = NULL, scale = NULL, logp = FALSE)`. Here `mu` is a two-column
#' matrix of natural scale and shape parameters. `p` and `q` can be scalars or
#' have one value per row of `mu`; `rd` returns one draw per row. `wt` and
#' `scale` are compatibility arguments and do not change the distribution.
#' Quantiles and simulations are excesses: add the appropriate threshold to
#' recover original measurement units. Invalid probabilities or parameters
#' produce `NaN` with a warning; missing inputs propagate to the output.
#'
#' Response residuals require shape less than one; Pearson residuals require
#' shape less than one half. Undefined moments produce `NA` with a warning.
#' Deviance residuals compare with a saturated scale equal to the observed
#' excess, holding fitted shape fixed, and use the sign of `excess - scale`.
#' They are undefined for shape at or below minus one, and are negative infinity
#' at zero excess (the saturated likelihood is unbounded). Zero-weight
#' deviance/Pearson residuals are zero. Null deviance is not supplied.
#' Probability-integral-transform residuals can instead be calculated with
#' `qnorm(family$cdf(excess, fitted(fit)))`.
#'
#' @return An mgcv `general.family` object.
#' @export
#' @examples
#' set.seed(12)
#' d <- data.frame(time = runif(300))
#' fam <- gpd()
#' pars <- cbind(exp(0.5 * sin(2 * pi * d$time)), 0.1)
#' d$y <- 10 + fam$rd(pars)
#' pot <- subset(d, y > 10)
#' pot$excess <- pot$y - 10
#' fit <- mgcv::gam(list(excess ~ s(time, k = 6), ~ 1),
#'                  family = gpd(), data = pot, method = "REML")
#' predict(fit, type = "response")[1:5, ]
#' fit$family$rd(fitted(fit))[1:5]
#' fit$family$qf(0.95, fitted(fit))[1:5]
#' pit <- fit$family$cdf(pot$excess, fitted(fit))
#' stats::qqnorm(stats::qnorm(pit))
#' # Unrestricted shape link:
#' gpd(link = list("log", "identity"))
gpd <- function(link = list("log", "logshift")) {
  if (!(is.list(link) || is.character(link)) || length(link) != 2L ||
      any(lengths(link) != 1L) ||
      !identical(link[[1L]], "log") ||
      !is.character(link[[2L]]) || is.na(link[[2L]]) ||
      !link[[2L]] %in% c("logshift", "identity")) {
    stop('gpd requires links "log" and either "logshift" or "identity"')
  }
  shape <- stats::make.link("identity")
  shape$d2link <- shape$d3link <- shape$d4link <- function(mu) rep(0, length(mu))
  if (link[[2L]] == "logshift") {
    shape$linkfun <- function(mu) log(mu + 0.5)
    shape$linkinv <- function(eta) exp(eta) - 0.5
    shape$mu.eta <- function(eta) exp(eta)
    shape$valideta <- function(eta) all(is.finite(eta))
    shape$name <- "logshift"
    shape$d2link <- function(mu) -1 / (mu + 0.5)^2
    shape$d3link <- function(mu) 2 / (mu + 0.5)^3
    shape$d4link <- function(mu) -6 / (mu + 0.5)^4
  }
  scale <- stats::make.link("log")
  scale$d2link <- function(mu) -1 / mu^2
  scale$d3link <- function(mu) 2 / mu^3
  scale$d4link <- function(mu) -6 / mu^4
  structure(list(
    family = "gpd", link = unlist(link), nlp = 2L,
    linfo = list(scale, shape), tri = mgcv::trind.generator(2),
    ll = .gpd_ll, initialize.coef = .gpd_start,
    initialize = expression({
      n <- rep(1, nobs)
      start <- family$initialize.coef(y, x, weights, offset, start, family)
    }),
    postproc = expression({
      object$null.deviance <- NA_real_
      colnames(object$fitted.values) <- c("scale", "shape")
    }),
    predict = .gpd_predict, residuals = .gpd_residuals,
    rd = .gpd_rd, qf = .gpd_qf, cdf = .gpd_cdf,
    d2link = 1, d3link = 1, d4link = 1, ls = 1,
    available.derivs = 2, discrete.ok = FALSE
  ), class = c("general.family", "extended.family", "family"))
}

# H = log(1 + xi * t) / xi, with H(t, 0) = t.
# For pure shape derivatives, the power series in v = xi*t avoids subtractive
# cancellation. Differentiating the series b times gives coefficient
# (-1)^(j+b) (j+b)! / (j! (j+b+1)) for t^(b+1) v^j.
# Eighty terms at |v| <= 1/2 suffice through order four in double precision.
.gpd_H <- function(t, xi, order = 0L) {
  v <- xi * t
  out <- matrix(0, length(t), order + 1L)
  small <- abs(v) <= 0.5
  for (b in 0:order) {
    if (any(small)) {
      z <- v[small]
      s <- numeric(length(z))
      for (j in 80:0) {
        fac <- if (b == 0L) 1 else prod(j + seq_len(b))
        s <- s * z + (-1)^(j + b) * fac / (j + b + 1)
      }
      out[small, b + 1L] <- t[small]^(b + 1L) * s
    }
    if (any(!small)) {
      k <- !small
      out[k, b + 1L] <- if (b == 0L) log1p(v[k]) / xi[k] else
        ((-1)^(b - 1L) * factorial(b - 1L) * (t[k] / (1 + v[k]))^b -
           b * out[k, b]) / xi[k]
    }
  }
  out
}

# Analytic mixed derivatives in rho = log(sigma) and xi, in mgcv's packed order.
# For a >= 1, rho derivatives of H are rational functions. Shape derivatives
# of those rational functions give every mixed derivative through order four.
.gpd_derivatives <- function(y, rho, xi, order = 4L) {
  t <- y * exp(-rho)
  v <- xi * t
  w <- 1 + v
  h0 <- .gpd_H(t, xi, order)
  H <- function(a, b) {
    if (a == 0L) return(h0[, b + 1L])
    if (a == 1L) return((-1)^(b + 1L) * factorial(b) * (t / w)^(b + 1L))
    if (a == 2L) return((-1)^b * factorial(b + 1L) * (t / w)^(b + 1L) / w)
    if (a == 3L && b == 0L) return(-t * (1 - v) / w^3)
    if (a == 3L && b == 1L) return(t^2 * (4 - 2 * v) / w^4)
    t * (1 - 4 * v + v^2) / w^4
  }
  ans <- list(l0 = -rho - (1 + xi) * h0[, 1L])
  for (d in seq_len(order)) {
    ans[[paste0("l", d)]] <- vapply(0:d, function(b) {
      a <- d - b
      -(1 + xi) * H(a, b) - (if (b > 0L) b * H(a, b - 1L) else 0) -
        as.numeric(a == 1L && b == 0L)
    }, numeric(length(y)))
    # vapply simplifies to a vector when there is only one observation.
    dim(ans[[paste0("l", d)]]) <- c(length(y), d + 1L)
  }
  ans
}

.gpd_eta <- function(X, coef, offset) {
  if (!is.matrix(X)) stop("gpd requires an ordinary model matrix")
  jj <- attr(X, "lpi")
  eta <- matrix(0, nrow(X), 2L)
  for (j in 1:2) {
    eta[, j] <- X[, jj[[j]], drop = FALSE] %*% coef[jj[[j]]]
    if (length(offset) >= j && !is.null(offset[[j]])) eta[, j] <- eta[, j] + offset[[j]]
  }
  eta
}

.gpd_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                    d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                    D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.null(eta)) eta <- .gpd_eta(X, coef, offset)
  if (is.null(wt)) wt <- rep(1, length(y))
  rho <- eta[, 1L]
  xi <- family$linfo[[2L]]$linkinv(eta[, 2L])
  active <- wt > 0
  # Zero weights must contribute exactly zero even if parameters put these
  # observations outside support. Use benign values before any arithmetic.
  yy <- y; yy[!active] <- 0
  rho[!active] <- xi[!active] <- 0
  sigma <- exp(rho)
  valid <- is.finite(rho) & is.finite(sigma) & sigma > 0 & is.finite(xi) &
    is.finite(yy / sigma) & (1 + xi * (yy / sigma) > 0)
  if (family$link[2L] == "logshift") valid <- valid & (xi > -0.5)
  if (!all(valid)) {
    l0 <- rep(0, length(y)); l0[active & !valid] <- -Inf
    return(list(l = -Inf, l0 = l0))
  }
  ord <- if (deriv == 0L) 0L else if (deriv == 1L) 2L else if (deriv < 4L) 3L else 4L
  ds <- .gpd_derivatives(yy, rho, xi, ord)
  l0 <- ds$l0 * wt
  if (!deriv) return(list(l = sum(l0), l0 = l0))
  sh <- family$linfo[[2L]]
  es <- eta[, 2L]; es[!active] <- sh$linkfun(0)
  # The first likelihood argument is already log scale, so its link here is
  # identity. Public linfo still maps log scale to the natural parameter sigma.
  ig1 <- cbind(1, sh$mu.eta(es))
  g2 <- cbind(0, sh$d2link(xi))
  g3 <- if (ord >= 3L) cbind(0, sh$d3link(xi)) else 0
  g4 <- if (ord >= 4L) cbind(0, sh$d4link(xi)) else 0
  tri <- family$tri
  de <- mgcv::gamlss.etamu(ds$l1 * wt, ds$l2 * wt,
                          if (ord >= 3L) ds$l3 * wt else 0,
                          if (ord >= 4L) ds$l4 * wt else 0,
                          ig1, g2, g3, g4, tri$i2, tri$i3, tri$i4, deriv - 1L)
  ret <- mgcv::gamlss.gH(X, attr(X, "lpi"), de$l1, de$l2, tri$i2,
                        l3 = de$l3, i3 = tri$i3, l4 = de$l4, i4 = tri$i4,
                        d1b = d1b, d2b = d2b, deriv = deriv - 1L,
                        fh = fh, D = D, sandwich = sandwich)
  if (ncv) { ret$l1 <- de$l1; ret$l2 <- de$l2; ret$l3 <- de$l3 }
  ret$l <- sum(l0)
  ret$l0 <- l0
  ret
}

.gpd_start <- function(y, X, wt, offset, start, family) {
  if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) || any(y < 0))
    stop("gpd requires finite nonnegative excesses")
  if (is.null(wt)) wt <- rep(1, length(y))
  if (length(wt) != length(y) || any(!is.finite(wt)) || any(wt < 0) || !any(wt > 0))
    stop("gpd requires finite nonnegative weights with at least one positive weight")
  if (!any(y[wt > 0] > 0)) stop("gpd requires at least one positive excess with positive weight")
  if (!is.matrix(X)) stop("gpd requires an ordinary model matrix")
  if (!is.null(start)) {
    if (length(start) != ncol(X) || any(!is.finite(start)) ||
        !is.finite(family$ll(y, X, start, wt, family, offset)$l))
      stop("gpd starting coefficients must be finite and put the response within support")
    return(start)
  }
  jj <- attr(X, "lpi")
  # Project constant, interior parameter targets onto each predictor, allowing
  # offsets and no-intercept designs. Search several scales/shapes and retain
  # the best finite likelihood; never alter the observations to obtain a start.
  sw <- sqrt(wt / max(wt))
  project <- function(target, j) {
    if (length(offset) >= j && !is.null(offset[[j]])) target <- target - offset[[j]]
    b <- qr.coef(qr(X[, jj[[j]], drop = FALSE] * sw), target * sw)
    b[is.na(b)] <- 0
    b
  }
  base <- sum(y * (wt / sum(wt)))
  best <- -Inf; ans <- NULL
  for (shape in c(0, 0.1, 0.5)) for (mult in c(1, 2, 10, 100)) {
    b <- numeric(ncol(X))
    b[jj[[1L]]] <- project(rep(log(base) + log(mult), length(y)), 1L)
    b[jj[[2L]]] <- project(rep(family$linfo[[2L]]$linkfun(shape), length(y)), 2L)
    value <- family$ll(y, X, b, wt, family, offset)$l
    if (is.finite(value) && value > best) { ans <- b; best <- value }
  }
  if (is.null(ans)) stop("Could not find support-valid gpd starting coefficients; supply start")
  ans
}

.gpd_predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                         beta = NULL, off = NULL, Vb = NULL) {
  if (is.null(eta)) eta <- .gpd_eta(X, beta, off) else se <- FALSE
  fit <- cbind(scale = exp(eta[, 1L]), shape = family$linfo[[2L]]$linkinv(eta[, 2L]))
  out <- list(fit = fit)
  if (se) {
    jj <- attr(X, "lpi")
    s <- fit
    for (j in 1:2) {
      xj <- X[, jj[[j]], drop = FALSE]
      v <- rowSums((xj %*% Vb[jj[[j]], jj[[j]], drop = FALSE]) * xj)
      s[, j] <- abs(family$linfo[[j]]$mu.eta(eta[, j])) * sqrt(pmax(0, v))
    }
    out$se.fit <- s
  }
  out
}

.gpd_args <- function(x, mu) {
  if (!is.matrix(mu) || !is.numeric(mu) || ncol(mu) != 2L)
    stop("mu must be a numeric matrix with scale and shape columns")
  if (!is.numeric(x) || !(length(x) %in% c(1L, nrow(mu))))
    stop("Input must be numeric and scalar or have one value per row of mu")
  x <- rep_len(x, nrow(mu))
  missing <- is.na(x) | is.na(mu[, 1L]) | is.na(mu[, 2L])
  bad <- !missing & (!is.finite(mu[, 1L]) | mu[, 1L] <= 0 | !is.finite(mu[, 2L]))
  list(x = x, sigma = mu[, 1L], xi = mu[, 2L], missing = missing, bad = bad)
}

.gpd_qf <- function(p, mu, wt = NULL, scale = NULL) {
  a <- .gpd_args(p, mu)
  p <- a$x; xi <- a$xi; sigma <- a$sigma
  bad <- a$bad | (!a$missing & (p < 0 | p > 1))
  if (any(bad)) warning("NaNs produced by invalid gpd probabilities or parameters")
  out <- rep(NA_real_, length(p)); out[bad] <- NaN
  good <- !a$missing & !bad
  out[good & p == 0] <- 0
  end <- good & p == 1
  out[end] <- ifelse(xi[end] < 0, -sigma[end] / xi[end], Inf)
  k <- good & p > 0 & p < 1
  h <- -log1p(-p[k]); v <- xi[k] * h
  ratio <- rep(1, length(v)); nz <- v != 0 & abs(v) <= 0.5
  ratio[nz] <- expm1(v[nz]) / v[nz]
  out[k] <- sigma[k] * h * ratio
  # Evaluate large-magnitude shape products on the log scale: intermediate
  # exp(v) or q/sigma can overflow even when the final result is representable.
  positive <- k; positive[k] <- v > 0.5
  vp <- v[v > 0.5]
  out[positive] <- exp(log(sigma[positive]) - log(xi[positive]) +
                        vp + log1p(-exp(-vp)))
  negative <- k; negative[k] <- v < -0.5
  out[negative] <- exp(log(sigma[negative]) - log(-xi[negative]) +
                        log(-expm1(v[v < -0.5])))
  out
}

.gpd_cdf <- function(q, mu, wt = NULL, scale = NULL, logp = FALSE) {
  a <- .gpd_args(q, mu)
  q <- a$x; xi <- a$xi; sigma <- a$sigma
  if (any(a$bad)) warning("NaNs produced by invalid gpd parameters")
  out <- rep(NA_real_, length(q)); out[a$bad] <- NaN
  good <- !a$missing & !a$bad
  out[good & q <= 0] <- if (logp) -Inf else 0
  end <- good & (q == Inf | (xi < 0 & q >= -sigma / xi))
  out[end] <- if (logp) 0 else 1
  k <- good & q > 0 & !end
  t <- q[k] / sigma[k]
  h <- rep(Inf, length(t))
  finite <- is.finite(t)
  h[finite] <- .gpd_H(t[finite], xi[k][finite])[, 1L]
  positive <- xi[k] > 0 & (!is.finite(t) | !is.finite(xi[k] * t))
  if (any(positive)) {
    logv <- log(xi[k][positive]) + log(q[k][positive]) - log(sigma[k][positive])
    h[positive] <- (pmax(logv, 0) + log1p(exp(-abs(logv)))) / xi[k][positive]
  }
  if (logp) out[k] <- ifelse(h > log(2), log1p(-exp(-h)), log(-expm1(-h))) else
    out[k] <- -expm1(-h)
  out
}

.gpd_rd <- function(mu, wt = NULL, scale = NULL) {
  .gpd_qf(stats::runif(nrow(mu)), mu, wt, scale)
}

.gpd_residuals <- function(object, type = c("deviance", "pearson", "response")) {
  type <- match.arg(type)
  y <- object$y; sigma <- object$fitted.values[, 1L]; xi <- object$fitted.values[, 2L]
  wt <- object$prior.weights
  if (is.null(wt)) wt <- rep(1, length(y))
  out <- rep(NA_real_, length(y))
  if (type == "deviance") {
    ok <- xi > -1 & wt > 0
    pos <- ok & y > 0
    sat <- -log(y[pos]) - (1 + xi[pos]) * .gpd_H(rep(1, sum(pos)), xi[pos])[, 1L]
    fit <- -log(sigma[pos]) - (1 + xi[pos]) * .gpd_H(y[pos] / sigma[pos], xi[pos])[, 1L]
    out[pos] <- sign(y[pos] - sigma[pos]) * sqrt(pmax(0, 2 * wt[pos] * (sat - fit)))
    out[ok & y == 0] <- -Inf
    if (any(xi <= -1 & wt > 0)) warning("gpd deviance residuals undefined for shape <= -1")
  } else {
    ok <- xi < (if (type == "response") 1 else 0.5)
    out[ok] <- y[ok] - sigma[ok] / (1 - xi[ok])
    if (type == "pearson") out[ok] <- out[ok] * sqrt(wt[ok]) *
      (1 - xi[ok]) * sqrt(1 - 2 * xi[ok]) / sigma[ok]
    if (any(!ok & wt > 0)) warning("gpd residuals undefined because required moments do not exist")
  }
  if (type != "response") out[wt == 0] <- 0
  out
}
