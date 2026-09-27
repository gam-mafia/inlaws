#' Log-logistic arithmetic mean and varying scale family
#'
#' Model the arithmetic mean and log-response logistic scale using two
#' additive predictors and the distribution described in [log_logistic()].
#' @param link The two links, `list("log", "logit")`, for mean and scale.
#'   Other links are not supported.
#' @details
#' Supply two formulae: the first includes the positive response and models
#' `log(mu)`; the second is one-sided and models `qlogis(s)`. Scale satisfies
#' `0 < s < 1`, preserving a finite arithmetic mean. Response variance is
#' infinite for `s >= 0.5`; these rows have `NA` Pearson residuals.
#' Responses are finite positive numeric vectors, without censoring.
#' Weights are likelihood replication weights, including zero weights.
#'
#' Supports `gam()` REML and NCV, including grouped neighbourhoods and
#' offsets in both predictors. `bam()`, `gamm()` and ML selection are not
#' supported. Separate coefficients are required for the two predictors;
#' shared coefficients via additional formulae are not supported.
#'
#' Fitted values and response predictions have columns `mu` and `s`.
#' Link predictions are `log(mu)` and `qlogis(s)`. Standard errors on the
#' response scale use the delta method. Predictions are plug-in conditional
#' values, without averaging over coefficient uncertainty or random effects.
#' `~ 1` specifies a constant scale predictor, but REML estimates need not
#' match [log_logistic()], which treats scale as a family parameter.
#'
#' The callbacks `rd(mu, wt = 1, scale = 1)`,
#' `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)` and
#' `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)` take a
#' two-column response-scale matrix as `mu`. They ignore `wt` and external
#' `scale`. `rd` draws one response per matrix row.
#'
#' Deviance maximizes each observation's likelihood over mean with scale
#' held fixed; its residual sign compares the observation to the median.
#' Null deviance fits a constant mean predictor retaining mean offsets and
#' fitted scales, providing a conditional mean-fit diagnostic. Log likelihood
#' and AIC include the original-response Jacobian and density constants.
#' @return A `general.family` object for [mgcv::gam()].
#' @seealso [log_logistic()]
#' @export
#' @examples
#' set.seed(22)
#' dat <- data.frame(x = runif(250), z = runif(250))
#' pars <- cbind(exp(sin(6*dat$x)), plogis(-1.5 + dat$z))
#' dat$y <- log_logisticls()$rd(pars)
#' fit <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(z, k = 5)), data = dat,
#'                  family = log_logisticls(), method = "REML")
#' head(predict(fit, type = "response"))
log_logisticls <- function(link = list("log", "logit")) {
  if (!identical(unname(unlist(link)), c("log", "logit")))
    stop('log_logisticls requires links "log" and "logit"')
  linfo <- lapply(link, stats::make.link)
  linfo[[1]]$linkinv <- linfo[[1]]$mu.eta <- function(eta) exp(eta)
  linfo[[2]]$linkinv <- stats::plogis
  linfo[[2]]$mu.eta <- function(eta) stats::plogis(eta)*stats::plogis(-eta)
  initialize <- expression({
    y <- family$prepare.response(y)
    n <- rep.int(1, nobs)
    start <- family$init(y, x, E, weights, offset, start)
  })
  postproc <- expression({
    colnames(object$fitted.values) <- c("mu", "s")
    object$null.deviance <- object$family$null.deviance(object, G$offset)
  })
  structure(list(family = "log_logisticls", link = unlist(link), nlp = 2L,
    linfo = linfo, tri = mgcv::trind.generator(2), ll = .llogls_ll,
    ncv = .llogls_ncv, initialize = initialize, init = .llogls_init,
    prepare.response = .llog_response,
    preinitialize = function(G) list(y = .llog_response(G$y), offxset = G$offset),
    postproc = postproc, residuals = .llogls_residuals,
    null.deviance = .llogls_null_deviance,
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    rd = function(mu, wt = 1, scale = 1) .llog_random(mu),
    qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
      .llog_quantile(p, mu, lower.tail = lower.tail, log.p = log.p),
    cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
      .llog_cdf(q, mu, logp = logp, lower.tail = lower.tail),
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
}

.llogls_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                           d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                           D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("log_logisticls does not support discrete model matrices")
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
    stats::plogis(eta[, 2]) > 0 & stats::plogis(eta[, 2]) < 1
  if (!all(valid)) return(list(l = -Inf, l0 = rep(-Inf, nrow(eta))))
  order <- if (deriv == 0) 0L else if (deriv == 1) 2L else if (deriv < 4) 3L else 4L
  y[wt == 0] <- 1
  r <- .llog_eval(y, eta[, 1], eta[, 2], order, predictor = TRUE)
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

.llogls_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                            Hi = NULL, R = NULL, offset = NULL, dH = NULL,
                            db = NULL, deriv = FALSE, nt = 1) {
  # Pass row indices to preserve correct offset selection for NCV subsets.
  response <- .llog_response(y)
  likelihood <- family$ll
  fit_offset <- lapply(1:2, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]])) rep_len(offset[[k]], length(response)) else
      numeric(length(response))
  })
  family$ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                        d1b = 0, ..., eta = NULL) {
    yy <- response[y]
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

.llogls_init <- function(y, x, E, wt, offset, start) {
  if (is.list(x)) stop("log_logisticls requires gam(), not bam()")
  jj <- attr(x, "lpi")
  if (length(jj) != 2L) stop("log_logisticls requires two formulae")
  if (anyDuplicated(unlist(jj))) stop("log_logisticls does not support shared predictor coefficients")
  if (any(!is.finite(wt) | wt < 0) || !any(wt > 0))
    stop("log_logisticls requires finite non-negative weights with positive total")
  if (!is.null(start)) return(start)
  keep <- wt > 0
  off <- lapply(1:2, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]])) rep_len(offset[[k]], length(y)) else
      numeric(length(y))
  })
  target <- log(y)
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
  sigma <- sqrt(stats::weighted.mean((target[keep] - loc[keep])^2, wt[keep]))
  sigma <- max(0.05, min(sigma*sqrt(3)/pi, 0.8))
  start[jj[[2]]] <- solve_start(2, rep(stats::qlogis(sigma), length(y)))
  sd <- stats::plogis(drop(x[, jj[[2]], drop = FALSE] %*% start[jj[[2]]]) + off[[2]])
  start[jj[[1]]] <- solve_start(1, target - .llog_correction(sd))
  start
}

.llogls_deviance <- function(y, mu, wt) {
  wt <- rep_len(wt, length(y))
  y[wt == 0] <- 1
  mu[wt == 0, 1] <- 1
  mu[wt == 0, 2] <- 0.25
  r <- .llog_eval(y, log(mu[, 1]), stats::qlogis(mu[, 2]), predictor = TRUE)
  d <- 2*wt*(r$sat[, 1] - r$ell[, 1])
  d[wt == 0] <- 0
  pmax(0, d)
}

.llogls_residuals <- function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values[, 1]; s <- object$fitted.values[, 2]
  if (type == "deviance") {
    location <- log(mu) + .llog_correction(s)
    return(sign(log(y)-location)*sqrt(.llogls_deviance(y, object$fitted.values, object$prior.weights)))
  }
  ans <- y - mu
  if (type != "response") {
    v <- .llog_variance(mu, s)
    ans <- ans*sqrt(object$prior.weights/v)
    ans[!is.finite(v)] <- NA_real_
  }
  ans
}

.llogls_null_deviance <- function(object, offset) {
  y <- .llog_response(object$y)
  w <- object$prior.weights
  mu <- object$fitted.values
  off <- if (length(offset) && !is.null(offset[[1]])) rep_len(offset[[1]], length(y)) else
    numeric(length(y))
  objective <- function(a) {
    candidate <- cbind(exp(a + off), mu[, 2])
    if (any(!is.finite(candidate) | candidate <= 0)) return(.Machine$double.xmax/100)
    sum(.llogls_deviance(y, candidate, w))
  }
  initial <- stats::weighted.mean(log(mu[, 1]) - off, w)
  opt <- stats::optim(initial, objective, method = "BFGS")
  opt$value
}
