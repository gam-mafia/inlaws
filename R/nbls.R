#' Negative binomial location and size family
#'
#' Model the mean and size of a negative binomial distribution with two
#' additive predictors. The variance is `mu + mu^2 / theta`: larger `theta`
#' means less overdispersion.
#'
#' @param link Two links, for the mean and size respectively. The mean accepts
#'   `"log"`, `"identity"`, or `"sqrt"`; size requires `"log"`.
#' @return A `general.family` object for [mgcv::gam()].
#' @md
#' @details Supply a list of two formulae, with the response in the first.
#'   The second predictor models `log(theta)`; use `~ 1` for constant size.
#'   Both predictors can contain offsets. Fitted values and response predictions
#'   have two columns, in the order `mu`, `theta`. Link predictions contain
#'   `g(mu)`, `log(theta)`. Response prediction standard errors use the delta
#'   method. Identity and square-root mean links require positive predictors.
#'
#'   REML fitting supports outer Newton, outer BFGS, and extended Fellner--Schall
#'   (`optimizer = "efs"`), using analytic derivatives through fourth order.
#'   NCV, `bam()`, and `gamm()` are not supported. Observations must be
#'   nonnegative integer counts, with finite nonnegative prior weights.
#'   An all-zero response among positive-weight observations is not identifiable
#'   and is rejected. Size can tend to infinity for Poisson or underdispersed
#'   data; it is not artificially bounded, and finite size estimates and their
#'   standard errors need not be identifiable in this case.
#'
#'   Deviance residuals compare the fitted mean to the saturated mean, keeping
#'   each fitted size fixed. Null deviance optimizes a constant mean predictor
#'   with the mean offset and fitted sizes retained: it is a conditional
#'   mean-fit diagnostic, not a measure of fit of both predictors together.
#'   An intercept-only size predictor need not give the same REML estimates as
#'   [mgcv::nb()], because size is a regression coefficient here rather than an
#'   extra family parameter in the marginal likelihood.
#'
#'   The family supplies random generation (`rd`), quantile (`qf`), and CDF
#'   (`cdf`) callbacks. Prior weights
#'   weight the likelihood; they do not change the distribution simulated for
#'   an individual observation. Mathematical details are in
#'   `inst/maths/nbls.md` in the package source.
#' @export
#' @examples
#' set.seed(42)
#' dat <- data.frame(x = runif(500), z = runif(500))
#' dat$y <- rnbinom(500, mu = exp(1 + sin(2 * pi * dat$x)),
#'                 size = exp(1 + 2 * dat$z))
#' fit <- mgcv::gam(list(y ~ s(x, k = 8), ~ s(z, k = 6)),
#'                  family = nbls(), method = "REML", data = dat)
#' head(predict(fit, type = "response"))
#' head(fit$family$rd(fitted(fit)))
nbls <- function(link = list("log", "log")) {
  if (length(link) != 2L ||
      !all(vapply(link, function(x) is.character(x) && length(x) == 1L &&
                    !is.na(x), logical(1))) ||
      !link[[1L]] %in% c("log", "identity", "sqrt") || link[[2L]] != "log") {
    stop('nbls requires two links: "log", "identity", or "sqrt" for mu, and "log" for theta')
  }
  linfo <- lapply(link, function(lnk) {
    z <- stats::make.link(lnk)
    # Unlike make.link("log"), do not floor exp(eta): its derivatives must
    # describe exactly the same function as its inverse link.
    if (lnk == "log") z$linkinv <- z$mu.eta <- function(eta) exp(eta)
    f <- mgcv::fix.family.link(structure(list(link = lnk, canonical = "none",
      linkfun = z$linkfun, mu.eta = z$mu.eta), class = "family"))
    z[c("d2link", "d3link", "d4link")] <- f[c("d2link", "d3link", "d4link")]
    z
  })
  initialize <- expression({
    n <- rep.int(1, nobs)
    start <- family$init(y, x, E, weights, offset, start, family)
  })
  postproc <- expression({
    colnames(object$fitted.values) <- c("mu", "theta")
    object$null.deviance <- object$family$null.deviance(object, G$offset)
  })
  structure(list(family = "nbls", link = unlist(link), nlp = 2L,
    linfo = linfo, tri = mgcv::trind.generator(2), ll = .nbls_ll,
    # mgcv 1.9-3 passes G$offxset (a typo) to its general-family EFS fitter.
    # Its supported preinitialize hook can supply this alias without changing
    # mgcv or capturing training offsets in the likelihood closure.
    preinitialize = function(G) list(offxset = G$offset),
    initialize = initialize, init = .nbls_initialize, postproc = postproc,
    null.deviance = .nbls_null_deviance, residuals = .nbls_residuals,
    sandwich = function(y, X, coef, wt, family, offset = NULL) {
      family$ll(y, X, coef, wt, family, offset = offset,
                deriv = 1, sandwich = TRUE)$lbb
    },
    rd = function(mu, wt = NULL, scale = NULL) {
      stats::rnbinom(nrow(mu), mu = mu[, 1L], size = mu[, 2L])
    },
    qf = function(p, mu, wt = NULL, scale = NULL) {
      stats::qnbinom(p, mu = mu[, 1L], size = mu[, 2L])
    },
    cdf = function(q, mu, wt = NULL, scale = NULL, logp = FALSE) {
      stats::pnbinom(q, mu = mu[, 1L], size = mu[, 2L], log.p = logp)
    },
    d2link = 1, d3link = 1, d4link = 1, ls = 1,
    available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
}

# Pure theta derivatives. In the Poisson limit, evaluate a convergent series
# instead of subtracting almost equal polygamma and rational terms.
.nbls_theta_derivatives <- function(y, mu, theta, order) {
  ans <- matrix(0, length(y), order)
  series <- theta > 100 * pmax(1, y, mu)
  ii <- which(!series)
  if (length(ii)) {
    yy <- y[ii]; mm <- mu[ii]; tt <- theta[ii]; z <- mm + tt
    for (q in seq_len(order)) {
      delta <- psigamma(tt + yy, deriv = q - 1L) - psigamma(tt, deriv = q - 1L)
      # Exact recurrence for small counts also covers y = 0 without cancellation.
      small <- yy <= 50
      delta[small] <- 0
      if (any(small)) for (k in 0:49) {
        use <- small & yy > k
        delta[use] <- delta[use] + (-1)^(q - 1L) * factorial(q - 1L) /
          (tt[use] + k)^q
      }
      if (q == 1L) {
        ans[ii, q] <- delta - log1p(mm / tt) + (mm - yy) / z
      } else {
        ans[ii, q] <- delta + (-1)^q * factorial(q - 2L) *
          (1 / tt^(q - 1L) - 1 / z^(q - 1L)) +
          (-1)^q * factorial(q - 1L) * (yy - mm) / z^q
      }
    }
  }
  ii <- which(series)
  if (length(ii)) {
    yy <- y[ii]; mm <- mu[ii]; tt <- theta[ii]
    # powers[, r + 1] = sum_{k=0}^{y-1} (k/theta)^r, obtained from
    # telescoping (k+1)^(r+1) - k^(r+1), without large integer powers.
    powers <- matrix(0, length(ii), 13L)
    powers[, 1L] <- yy
    for (r in 1:12) {
      sr <- yy * (yy / tt)^r
      for (j in 0:(r - 1L)) sr <- sr - choose(r + 1L, j) *
        powers[, j + 1L] / tt^(r - j)
      powers[, r + 1L] <- sr / (r + 1L)
      cr <- (-1)^(r + 1L) * ((powers[, r + 1L] - yy * (mm / tt)^r) / r +
                               mm * (mm / tt)^r / (r + 1L))
      for (q in seq_len(order)) ans[ii, q] <- ans[ii, q] +
        (-1)^q * prod(r:(r + q - 1L)) * cr / tt^q
    }
  }
  ans
}

# Natural-parameter derivatives, packed mu...mu, mu...theta, ..., theta...theta.
.nbls_derivatives <- function(y, mu, theta, order = 4L) {
  pure <- .nbls_theta_derivatives(y, mu, theta, order)
  z <- mu + theta
  out <- vector("list", order)
  for (r in seq_len(order)) {
    a <- matrix(0, length(y), r + 1L)
    for (q in 0:r) {
      p <- r - q
      if (p == 0L) {
        a[, q + 1L] <- pure[, q]
      } else if (q == 0L) {
        a[, 1L] <- if (p == 1L) (theta / z) * (y - mu) / mu else
          (-1)^(p - 1L) * factorial(p - 1L) *
          (y / mu^p - (y + theta) / z^p)
      } else {
        a[, q + 1L] <- (-1)^(p + q) * factorial(p - 1L) *
          (prod(p:(p + q - 1L)) * (y - mu) / z^(p + q) +
             prod((p - 1L):(p + q - 2L)) / z^(p + q - 1L))
      }
    }
    out[[r]] <- a
  }
  names(out) <- paste0("l", seq_len(order))
  out
}

.nbls_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                      d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                      D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("nbls does not support discrete model matrices")
  if (ncv) stop("nbls does not support NCV")
  jj <- attr(X, "lpi")
  if (length(jj) != 2L) stop("nbls requires two linear predictors")
  n <- length(y)
  if (is.null(wt)) wt <- rep.int(1, n)
  if (is.null(eta)) {
    eta <- matrix(0, n, 2L)
    for (j in 1:2) {
      eta[, j] <- drop(X[, jj[[j]], drop = FALSE] %*% coef[jj[[j]]])
      if (length(offset) >= j && !is.null(offset[[j]])) eta[, j] <- eta[, j] + offset[[j]]
    }
  }
  active <- wt > 0
  # Inactive rows have no likelihood or derivative contribution, even if their
  # predictors overflow. Substitute harmless parameters before transforming.
  eta[!active, ] <- 1
  mu <- family$linfo[[1L]]$linkinv(eta[, 1L])
  theta <- family$linfo[[2L]]$linkinv(eta[, 2L])
  valid <- is.finite(mu) & mu > 0 & is.finite(theta) & theta > 0 &
    is.finite(eta[, 1L]) & is.finite(eta[, 2L])
  if (family$link[[1L]] != "log") valid <- valid & eta[, 1L] > 0
  if (!all(valid)) return(list(l = -Inf, l0 = rep(-Inf, n),
    lb = rep(0, ncol(X)), lbb = matrix(0, ncol(X), ncol(X))))
  l0 <- numeric(n)
  l0[active] <- wt[active] * stats::dnbinom(y[active], mu = mu[active],
                                         size = theta[active], log = TRUE)
  ret <- list()
  if (deriv > 0L) {
    order <- if (deriv == 1L) 2L else if (deriv < 4L) 3L else 4L
    yy <- y; yy[!active] <- 0
    dd <- .nbls_derivatives(yy, mu, theta, order)
    for (j in seq_along(dd)) dd[[j]] <- dd[[j]] * wt
    ig1 <- cbind(family$linfo[[1L]]$mu.eta(eta[, 1L]), theta)
    g <- lapply(2:order, function(k) cbind(
      family$linfo[[1L]][[paste0("d", k, "link")]](mu),
      family$linfo[[2L]][[paste0("d", k, "link")]](theta)))
    tri <- family$tri
    de <- mgcv::gamlss.etamu(dd$l1, dd$l2, dd$l3, dd$l4,
      ig1, g[[1L]], if (order >= 3L) g[[2L]] else NULL,
      if (order == 4L) g[[3L]] else NULL,
      tri$i2, tri$i3, tri$i4, deriv = deriv - 1L)
    shared <- deriv == 4L && anyDuplicated(unlist(jj)) > 0L
    ret <- mgcv::gamlss.gH(X, jj, de$l1, de$l2, tri$i2,
      l3 = de$l3, i3 = tri$i3, l4 = de$l4, i4 = tri$i4,
      d1b = d1b, d2b = d2b, deriv = if (shared) 2L else deriv - 1L, fh = fh, D = D,
      sandwich = sandwich)
    # mgcv 1.9-3/1.9-4's fourth-order trace accumulator overwrites overlapping
    # predictor columns. Accumulate these contractions explicitly for shared
    # coefficients; its gradient, Hessian, and third-order paths are additive.
    if (shared) ret$trHid2H <- .nbls_shared_trace(X, jj, de, tri, d1b, d2b, fh, D)
  }
  ret$l <- sum(l0)
  ret$l0 <- l0
  ret
}

.nbls_shared_trace <- function(X, jj, de, tri, d1b, d2b, fh, D) {
  if (is.null(D)) D <- rep(1, ncol(X))
  if (is.list(fh)) {
    ev <- fh$values
    ev <- ifelse(ev > 0, 1 / ev, 0)
    V <- (fh$vectors %*% (ev * t(fh$vectors))) * outer(D, D)
  } else {
    pivot <- attr(fh, "pivot")
    if (is.null(pivot)) pivot <- seq_len(ncol(X))
    V <- chol2inv(fh)[order(pivot), order(pivot), drop = FALSE] * outer(D, D)
  }
  x <- lapply(jj, function(j) X[, j, drop = FALSE])
  a <- lapply(seq_along(jj), function(j) x[[j]] %*% d1b[jj[[j]], , drop = FALSE])
  b <- lapply(seq_along(jj), function(j) x[[j]] %*% d2b[jj[[j]], , drop = FALSE])
  out <- numeric(ncol(d2b))
  for (i in seq_along(jj)) for (j in i:length(jj)) {
    h <- rowSums((x[[i]] %*% V[jj[[i]], jj[[j]], drop = FALSE]) * x[[j]])
    kk <- 0L
    for (k in seq_len(ncol(d1b))) for (l in k:ncol(d1b)) {
      kk <- kk + 1L
      v <- numeric(nrow(X))
      for (q in seq_along(jj)) {
        v <- v + de$l3[, tri$i3[i, j, q]] * b[[q]][, kk]
        for (s in seq_along(jj)) v <- v + de$l4[, tri$i4[i, j, q, s]] *
          a[[q]][, k] * a[[s]][, l]
      }
      out[kk] <- out[kk] + (if (i == j) 1 else 2) * sum(h * v)
    }
  }
  out
}

.nbls_initialize <- function(y, x, E, wt, offset, start, family) {
  if (is.null(wt)) wt <- rep.int(1, length(y))
  if (any(!is.finite(y) | y < 0 | y != floor(y)))
    stop("nbls requires finite nonnegative integer counts")
  if (any(!is.finite(wt) | wt < 0)) stop("nbls requires finite nonnegative weights")
  keep <- wt > 0
  if (!any(keep)) stop("nbls requires positive-weight observations")
  if (!any(y[keep] > 0)) stop("nbls cannot fit an all-zero positive-weight response")
  jj <- attr(x, "lpi")
  if (length(jj) != 2L) stop("nbls requires two linear predictors")
  if (!is.null(start)) return(start)
  off <- lapply(1:2, function(j) {
    if (length(offset) >= j && !is.null(offset[[j]])) rep_len(offset[[j]], length(y)) else
      numeric(length(y))
  })
  # Fit both blocks together so coefficients shared by predictors are handled
  # once. Penalties only regularize these starting regressions.
  design <- lapply(jj, function(ind) {
    a <- matrix(0, sum(keep), ncol(x))
    a[, ind] <- x[keep, ind, drop = FALSE]
    a * sqrt(wt[keep])
  })
  ee <- E
  if (length(ee) && sum(ee^2) > 0) ee <- ee *
    (0.01 * sqrt(sum(x[keep, , drop = FALSE]^2) / sum(ee^2)))
  solve_start <- function(target) {
    a <- do.call(rbind, c(design, list(ee)))
    b <- c(unlist(lapply(1:2, function(j)
      (target[[j]][keep] - off[[j]][keep]) * sqrt(wt[keep]))), rep(0, nrow(ee)))
    ans <- qr.coef(qr(a), b)
    ans[!is.finite(ans)] <- 0
    ans
  }
  target <- list(family$linfo[[1L]]$linkfun(pmax(y, 1/6)), rep(0, length(y)))
  start <- solve_start(target)
  eta <- drop(x[, jj[[1L]], drop = FALSE] %*% start[jj[[1L]]]) + off[[1L]]
  mu <- family$linfo[[1L]]$linkinv(eta)
  num <- sum(wt[keep] * mu[keep]^2)
  den <- sum(wt[keep] * ((y[keep] - mu[keep])^2 - mu[keep]))
  theta <- num / den
  if (!is.finite(theta) || theta <= 0) theta <- 1
  target[[2L]][] <- log(theta)
  start <- solve_start(target)
  # A log link always has an interior start. For constrained mean links, find
  # an interior point by fitting a constant positive predictor-scale target,
  # increasing it if offsets require this. If this fails, request a feasible
  # user-supplied start rather than changing the fitted link or mean.
  if (family$link[[1L]] != "log") {
    for (attempt in 0:12) {
      eta <- drop(x[, jj[[1L]], drop = FALSE] %*% start[jj[[1L]]]) + off[[1L]]
      if (all(is.finite(eta[keep]) & eta[keep] > 0)) break
      target[[1L]][] <- 2^attempt * max(1, mean(y[keep]), abs(off[[1L]][keep]))
      start <- solve_start(target)
    }
    eta <- drop(x[, jj[[1L]], drop = FALSE] %*% start[jj[[1L]]]) + off[[1L]]
    if (any(!is.finite(eta[keep]) | eta[keep] <= 0))
      stop("could not initialize a positive mean predictor; supply start or use a log link")
  }
  start
}

.nbls_deviance <- function(y, mu, theta, wt) {
  ans <- numeric(length(y))
  ii <- wt > 0
  ans[ii] <- pmax(0, 2 * wt[ii] *
    (stats::dnbinom(y[ii], mu = y[ii], size = theta[ii], log = TRUE) -
       stats::dnbinom(y[ii], mu = mu[ii], size = theta[ii], log = TRUE)))
  ans
}

.nbls_residuals <- function(object, type = c("deviance", "pearson", "response")) {
  type <- match.arg(type)
  mu <- object$fitted.values[, 1L]; theta <- object$fitted.values[, 2L]
  r <- object$y - mu
  if (type == "deviance") r <- sign(r) * sqrt(.nbls_deviance(object$y, mu, theta,
                                                            object$prior.weights))
  if (type == "pearson") r <- sqrt(object$prior.weights) * r / sqrt(mu + mu^2 / theta)
  r
}

.nbls_null_deviance <- function(object, offset = NULL) {
  ii <- object$prior.weights > 0
  y <- object$y[ii]; wt <- object$prior.weights[ii]
  theta <- object$fitted.values[ii, 2L]
  off <- if (length(offset) && !is.null(offset[[1L]]))
    rep_len(offset[[1L]], length(ii))[ii] else numeric(sum(ii))
  loglink <- object$family$link[[1L]] == "log"
  lower <- -min(off) + sqrt(.Machine$double.eps)
  inverse <- object$family$linfo[[1L]]$linkinv
  objective <- function(a) {
    eta <- (if (loglink) a else lower + exp(a)) + off
    mu <- inverse(eta)
    if (any(!is.finite(mu) | mu <= 0)) return(.Machine$double.xmax / 100)
    -sum(wt * stats::dnbinom(y, mu = mu, size = theta, log = TRUE))
  }
  initial <- object$family$linfo[[1L]]$linkfun(stats::weighted.mean(y, wt)) -
    stats::weighted.mean(off, wt)
  if (!loglink) initial <- log(max(initial - lower, 1e-3))
  opt <- stats::optim(initial, objective, method = "BFGS")
  eta <- (if (loglink) opt$par else lower + exp(opt$par)) + off
  sum(.nbls_deviance(y, inverse(eta), theta, wt))
}
