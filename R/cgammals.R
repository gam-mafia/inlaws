#' Censored gamma mean and dispersion family
#'
#' Model the conditional arithmetic mean and dispersion of a gamma response
#' with separate additive predictors, allowing non-informative censoring.
#'
#' @param link Two links, for the mean and dispersion respectively. Both must
#'   be `"log"`.
#' @return A `general.family` object for [mgcv::gam()].
#' @details Supply two formulae, with the response in the first and a one-sided
#'   dispersion formula second. Both support offsets. The model has mean
#'   `mu = exp(eta1)`, dispersion `phi = exp(eta2)`, shape `1/phi`, scale
#'   `mu*phi`, and variance `phi*mu^2`. There is no artificial dispersion bound;
#'   very small dispersion can produce singular fits.
#'
#'   Use a positive vector for exact responses or the two-column censoring
#'   matrix described in [cgamma()]. All censoring types may be mixed; zero
#'   lower bounds are supported. Missing rows are handled by `gam()`'s
#'   `na.action`. Finite nonnegative case weights multiply likelihood terms,
#'   including zero weights, and do not change the response distribution.
#'
#'   REML and NCV smoothing selection are supported, including neighbourhoods
#'   supplied through `nei`. `bam()`, `gamm()`, ML smoothing selection, QNCV,
#'   and shared coefficients between predictors are outside the supported
#'   scope. The same covariates and smooths may appear in both formulae.
#'
#'   Fitted values and response predictions have columns `mu`, `phi`; link
#'   predictions are their natural logarithms. Response standard errors use
#'   the delta method. Predictions do not integrate coefficient uncertainty
#'   or random effects. Unlike [mgcv::gammals()], response predictions return
#'   dispersion itself rather than log dispersion, and both links are ordinary
#'   log links.
#'
#'   The distribution callbacks take a two-column response-scale matrix `mu`:
#'   `rd(mu, wt = 1, scale = 1)`,
#'   `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)`, and
#'   `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)`.
#'   They describe uncensored responses and ignore `wt` and external `scale`.
#'
#'   Use `~ 1` for constant dispersion. REML estimates may differ from
#'   [cgamma()] because dispersion is a regression coefficient here rather
#'   than a family parameter. Log likelihood and AIC include density constants.
#'   Deviance compares with the saturated mean at each fitted dispersion;
#'   null deviance fits a constant mean with fitted dispersions and mean
#'   offsets retained. Response and Pearson residuals are `NA` for censored
#'   observations. Mathematical details are in `inst/maths/cgammals.md`.
#' @export
#' @examples
#' set.seed(12)
#' dat <- data.frame(x = runif(200), z = runif(200))
#' mu <- exp(0.5 + sin(2*pi*dat$x))
#' phi <- exp(-1 + dat$z)
#' y <- rgamma(200, shape = 1/phi, scale = mu*phi)
#' dat$y <- cbind(pmax(y, 0.4), ifelse(y < 0.4, -Inf, y))
#' fit <- mgcv::gam(list(y ~ s(x, k = 7), ~ s(z, k = 5)), data = dat,
#'                  family = cgammals(), method = "REML")
#' head(predict(fit, type = "response"))
cgammals <- function(link = list("log", "log")) {
  if (length(link) != 2L || !all(vapply(link, function(x)
      is.character(x) && length(x) == 1L && !is.na(x) && x == "log", logical(1))))
    stop('cgammals requires two "log" links')
  linfo <- lapply(link, function(x) {
    lk <- stats::make.link(x)
    # Use exact exponentials, consistently with the predictor derivatives.
    lk$linkinv <- lk$mu.eta <- function(eta) exp(eta)
    lk
  })
  initialize <- expression({
    y <- family$prepare.response(y)
    n <- rep.int(1, nobs)
    start <- family$init(y, x, E, weights, offset, start)
  })
  postproc <- expression({
    colnames(object$fitted.values) <- c("mu", "phi")
    object$null.deviance <- object$family$null.deviance(object, G$offset)
  })
  structure(list(family = "cgammals", link = unlist(link), nlp = 2L,
    linfo = linfo, tri = mgcv::trind.generator(2), ll = .cgammals_ll,
    ncv = .cgammals_ncv, initialize = initialize, init = .cgammals_init,
    prepare.response = .cgammals_response,
    # General-family fitters count length(y) before initialize is evaluated.
    preinitialize = function(G) list(y = .cgammals_response(G$y), offxset = G$offset),
    postproc = postproc, residuals = .cgammals_residuals,
    null.deviance = .cgammals_null_deviance,
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    rd = function(mu, wt = 1, scale = 1)
      stats::rgamma(nrow(mu), shape = 1/mu[, 2], scale = mu[, 1]*mu[, 2]),
    qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
      stats::qgamma(p, shape = 1/mu[, 2], scale = mu[, 1]*mu[, 2],
                    lower.tail = lower.tail, log.p = log.p),
    cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
      stats::pgamma(q, shape = 1/mu[, 2], scale = mu[, 1]*mu[, 2],
                    lower.tail = lower.tail, log.p = logp),
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
}

.cgammals_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                           d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                           D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("cgammals does not support discrete model matrices")
  jj <- attr(X, "lpi")
  if (is.null(eta)) {
    eta <- matrix(0, nrow(X), 2L)
    for (k in 1:2) {
      eta[, k] <- drop(X[, jj[[k]], drop = FALSE] %*% coef[jj[[k]]])
      if (length(offset) >= k && !is.null(offset[[k]]))
        eta[, k] <- eta[, k] + offset[[k]]
    }
  }
  # Direct eta already includes offsets in gamlss.ncv; do not add them twice.
  wt <- rep_len(if (is.null(wt)) 1 else wt, nrow(eta))
  eta[wt == 0, ] <- 0
  valid <- is.finite(eta[, 1]) & is.finite(eta[, 2]) &
    is.finite(exp(eta[, 1])) & exp(eta[, 1]) > 0 &
    is.finite(exp(eta[, 2])) & exp(eta[, 2]) > 0 &
    is.finite(exp(-eta[, 2])) & exp(-eta[, 2]) > 0
  if (!all(valid)) return(list(l = -Inf, l0 = rep(-Inf, nrow(eta))))
  order <- if (deriv == 0) 0L else if (deriv == 1) 2L else if (deriv < 4) 3L else 4L
  # Do not evaluate irrelevant extreme observations at zero-weight rows.
  yy <- .cgammals_response(y)
  bounds <- attr(yy, "censor")
  yy[wt == 0] <- bounds[wt == 0] <- 1
  attr(yy, "censor") <- bounds
  r <- .cgammals_eval(yy, eta[, 1], eta[, 2], order)
  ell <- r$ell * wt
  ell[wt == 0, ] <- 0
  ret <- list()
  if (deriv > 0) {
    packed <- lapply(seq_len(order), function(k) {
      out <- matrix(0, nrow(eta), k + 1L)
      for (b in 0:k) {
        a <- k - b
        out[, b + 1L] <- ell[, which(r$powers[, 1] == a & r$powers[, 2] == b)] *
          factorial(a)*factorial(b)
      }
      out
    })
    l3 <- if (order >= 3) packed[[3]] else 0
    l4 <- if (order >= 4) packed[[4]] else 0
    # NCV supplies eta for row subsets whose X may lack lpi. It only needs
    # observation derivatives, not coefficient-space contractions there.
    if (!is.null(jj)) ret <- mgcv::gamlss.gH(X, jj, packed[[1]], packed[[2]], family$tri$i2,
      l3 = l3, i3 = family$tri$i3, l4 = l4, i4 = family$tri$i4,
      d1b = d1b, d2b = d2b, deriv = deriv - 1L, fh = fh, D = D,
      sandwich = sandwich)
    if (ncv) {
      ret$l1 <- packed[[1]]; ret$l2 <- packed[[2]]; ret$l3 <- l3
      if (order >= 4) ret$l4 <- l4
    }
  }
  ret$l0 <- ell[, 1]; ret$l <- sum(ret$l0)
  ret
}

.cgammals_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                            Hi = NULL, R = NULL, offset = NULL, dH = NULL,
                            db = NULL, deriv = FALSE, nt = 1) {
  # gamlss.ncv subsets y with y[nei$d], dropping attributes (or flattening
  # matrices). Pass row indices and decode inside its likelihood callback.
  response <- .cgammals_response(y)
  bounds <- attr(response, "censor")
  likelihood <- family$ll
  fit_offset <- lapply(1:2, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]])) rep_len(offset[[k]], length(response)) else
      numeric(length(response))
  })
  family$ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                        d1b = 0, ..., eta = NULL) {
    yy <- response[y]
    attr(yy, "censor") <- bounds[y]
    off <- lapply(fit_offset, function(x) x[y])
    if (!is.null(eta)) for (k in 1:2) eta[, k] <- eta[, k] + off[[k]]
    likelihood(yy, X, coef, wt, family, offset = off, deriv = deriv,
               d1b = d1b, ..., eta = eta)
  }
  # Keep the helper's predictors offset-free and add offsets only in the
  # callback. Its gamma adjustment evaluates the full-data eta without
  # adding offsets, whereas its held-out eta normally includes them.
  ret <- utils::getFromNamespace("gamlss.ncv", "mgcv")(
    X, seq_along(response), wt, nei, beta, family, llf, H = H, Hi = Hi,
    R = R, offset = NULL, dH = dH, db = db, deriv = deriv, nt = nt)
  eta <- attr(ret$NCV, "eta.cv")
  if (!is.null(eta)) {
    for (k in 1:2) eta[, k] <- eta[, k] + fit_offset[[k]][nei$d]
    attr(ret$NCV, "eta.cv") <- eta
  }
  ret
}

.cgammals_init <- function(y, x, E, wt, offset, start) {
  if (is.list(x)) stop("cgammals requires gam(), not bam()")
  jj <- attr(x, "lpi")
  if (length(jj) != 2L) stop("cgammals requires two formulae")
  if (anyDuplicated(unlist(jj))) stop("cgammals does not support shared predictor coefficients")
  if (any(!is.finite(wt) | wt < 0) || !any(wt > 0))
    stop("cgammals requires finite non-negative weights with positive total")
  if (!is.null(start)) return(start)
  keep <- wt > 0
  off <- lapply(1:2, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]])) rep_len(offset[[k]], length(y)) else
      numeric(length(y))
  })
  # Censor bounds give rough, finite pseudo-responses only for initialization.
  rr <- .cg_response(y)
  pseudo <- rr$y
  pseudo[rr$censor == -Inf] <- rr$y[rr$censor == -Inf]/2
  ii <- !rr$exact & is.finite(rr$hi)
  pseudo[ii] <- rr$lo[ii]/2 + rr$hi[ii]/2
  positive <- pseudo[keep & pseudo > 0]
  fallback <- if (length(positive)) stats::median(positive) else 1
  pseudo <- pmax(pseudo, fallback * 0.01)
  target <- log(pseudo)
  solve_start <- function(k, target) {
    xx <- x[keep, jj[[k]], drop = FALSE] * sqrt(wt[keep])
    ee <- E[, jj[[k]], drop = FALSE]
    if (length(ee) && sum(ee^2) > 0)
      ee <- ee * (0.01 * sqrt(sum(xx^2)/sum(ee^2)))
    ans <- qr.coef(qr(rbind(xx, ee)),
      c((target[keep] - off[[k]][keep])*sqrt(wt[keep]), rep(0, nrow(ee))))
    ans[!is.finite(ans)] <- 0
    ans
  }
  start <- numeric(ncol(x))
  start[jj[[1]]] <- solve_start(1, target)
  loc <- drop(x[, jj[[1]], drop = FALSE] %*% start[jj[[1]]]) + off[[1]]
  phi <- stats::weighted.mean((pseudo[keep]/exp(loc[keep]) - 1)^2, wt[keep])
  phi <- max(0.05, min(phi, 5)) # Only an initialization safeguard.
  start[jj[[2]]] <- solve_start(2, rep(log(phi), length(y)))
  start
}

.cgammals_deviance <- function(y, mu, wt) {
  wt <- rep_len(wt, nrow(mu))
  keep <- wt > 0
  d <- numeric(nrow(mu))
  if (!any(keep)) return(d)
  r <- .cg_response(y)
  yy <- cbind(r$y[keep], r$censor[keep])
  actual <- .cgammals_eval(yy, log(mu[keep, 1]), log(mu[keep, 2]))$ell[, 1]
  rr <- .cg_response(yy)
  sat <- numeric(nrow(yy))
  ii <- rr$exact | (rr$lo > 0 & is.finite(rr$hi))
  if (any(ii)) {
    target <- .cgammals_saturated_mean(rr)
    sat[ii] <- .cgammals_eval(yy[ii, , drop = FALSE], log(target[ii]),
                              log(mu[keep, 2][ii]))$ell[, 1]
  }
  d[keep] <- pmax(0, 2*wt[keep]*(sat - actual))
  d
}

.cgammals_saturated_mean <- function(r) {
  target <- r$y
  ii <- !r$exact & r$lo > 0 & is.finite(r$hi)
  width <- r$hi[ii] - r$lo[ii]
  target[ii] <- width/log1p(width/r$lo[ii])
  target[!r$exact & r$lo == 0] <- 0
  target[is.infinite(r$hi)] <- Inf
  target
}

.cgammals_residuals <- function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
  type <- match.arg(type)
  r <- .cg_response(object$y)
  mu <- object$fitted.values[, 1]; phi <- object$fitted.values[, 2]
  if (type == "deviance") {
    sgn <- sign(.cgammals_saturated_mean(r) - mu)
    sgn[r$censor == -Inf] <- -1
    sgn[r$censor == Inf] <- 1
    return(sgn * sqrt(.cgammals_deviance(object$y, object$fitted.values, object$prior.weights)))
  }
  ans <- r$y - mu
  if (type != "response") ans <- ans*sqrt(object$prior.weights/(mu^2*phi))
  ans[!r$exact] <- NA_real_
  ans
}

.cgammals_null_deviance <- function(object, offset) {
  y <- .cgammals_response(object$y)
  w <- object$prior.weights
  mu <- object$fitted.values
  off <- if (length(offset) && !is.null(offset[[1]])) rep_len(offset[[1]], length(y)) else
    numeric(length(y))
  objective <- function(a) {
    candidate <- cbind(exp(a + off), mu[, 2])
    if (any(!is.finite(candidate) | candidate <= 0)) return(.Machine$double.xmax/100)
    sum(.cgammals_deviance(y, candidate, w))
  }
  initial <- stats::weighted.mean(log(mu[, 1]) - off, w)
  opt <- stats::optim(initial, objective, method = "BFGS")
  opt$value
}
