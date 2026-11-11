#' Weibull arithmetic mean and varying shape family
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Model the arithmetic mean and shape of a positive Weibull response with
#' two additive predictors.
#' @param link Two links, both `"log"`.
#' @return A `general.family` object for [mgcv::gam()] and [mgcv::scasm()].
#' @details Supply two formulae. The first includes the response and models
#' `log(mu)`; the second models `log(k)`. A second formula of `~ 1` estimates
#' constant shape. The conventional Weibull scale is
#' `lambda = mu / Gamma(1 + 1/k)`. See [weibull()] for the distribution and
#' conditional interpretation of the mean. Both predictors accept offsets.
#' Weights multiply the exact uncensored likelihood. Zero weights are allowed.
#' The two predictors require separate coefficient blocks.
#'
#' In `gam()`, REML supports outer Newton, outer BFGS, and extended
#' Fellner--Schall (EFS). EFS requires the smoothing parameters to be all
#' free and unlinked, or all fixed. NCV and QNCV use outer BFGS, with
#' leave-one-out or explicit `nei` neighbourhoods. The criterion uses
#' one-step deletion predictions. QNCV also approximates held-out likelihood
#' contributions quadratically. Use explicit `nei` with `jackknife = FALSE`
#' for fully unpenalized NCV/QNCV models on mgcv 1.9-4.
#' The general-family `nlm` and `optim` branches of mgcv 1.9-4 and 1.9-5
#' require a null-model objective scale that mgcv leaves unset.
#'
#' In `scasm()`, REML uses EFS. mgcv 1.9-4 and 1.9-5 expose no
#' `method`, `optimizer` or `nei` arguments in `scasm()`, so NCV and QNCV
#' are unavailable through that fitter. If active constraints remove all
#' free directions of a smooth, the upstream EFS penalty trace can fail
#' on an empty projected penalty. `bam()` and `gamm()` are outside the
#' supported scope.
#'
#' Fitted values have columns `mu` and `k`. Response predictions return
#' both parameters in that order. Link predictions contain their logarithms.
#' Response standard errors use the delta method and describe uncertainty
#' in the estimated parameters.
#' The `rd`, `cdf`, and `qf` callbacks accept a two-column response-scale
#' parameter matrix as `mu`, and otherwise use the arguments of [weibull()].
#' Native gratia diagnostics provide the same response distribution and moments.
#'
#' Deviance maximizes each observation's likelihood over the mean with shape
#' held fixed. The saturated conventional scale equals the observation, so the
#' saturated arithmetic mean is `y * Gamma(1 + 1/k)`. Deviance residuals
#' are signed relative to the fitted conventional scale.
#' Null deviance fits a constant log mean, retaining mean offsets and fitted
#' shapes. With no mean intercept it uses the offsets alone.
#' @seealso [weibull()], [family_diagnostics]
#' @template ncv-neighbourhoods
#' @export
#' @examples
#' set.seed(78)
#' d <- data.frame(x = runif(220), z = runif(220))
#' f <- weibulls()
#' d$y <- f$rd(cbind(exp(sin(6 * d$x)), exp(0.4 + d$z)))
#' b <- mgcv::gam(list(y ~ s(x, k = 6), ~ z), data = d, family = f,
#'                method = "REML")
#' head(predict(b, type = "response", se.fit = TRUE)$fit)
#' if ("scasm" %in% getNamespaceExports("mgcv")) {
#'   d$y <- f$rd(cbind(exp(1.5 * d$x), exp(0.7 + 0.5 * d$z)))
#'   constrained <- mgcv::scasm(list(y ~ s(x, bs = "sc", xt = "m+", k = 6), ~ 1),
#'                               data = d, family = f)
#' }
weibulls <- function(link = list("log", "log")) {
  if (!identical(unname(unlist(link)), c("log", "log")))
    stop('weibulls requires two "log" links', call. = FALSE)
  log_link <- .weibull_log_link()
  family <- structure(list(family = "weibulls", link = c("log", "log"),
    nlp = 2L, linfo = list(log_link, log_link), tri = mgcv::trind.generator(2),
    ll = .weibulls_ll, ncv = .inlaws_two_predictor_ncv,
    loglik = .weibull_eval,
    check.nei = .inlaws_check_ncv_nei,
    presetup = expression({
      # scasm has no method, optimizer or nei formals.
      if (exists("method", inherits = FALSE)) {
        family$check.nei(nei, attr(mf, "na.action"))
        if (missing(method)) method <- "REML"
        if (!method %in% c("REML", "NCV", "QNCV"))
          stop("weibulls supports REML, NCV and QNCV", call. = FALSE)
        if (method %in% c("NCV", "QNCV") || !is.null(nei))
          optimizer <- c("outer", "bfgs")
      }
    }),
    preinitialize = .weibulls_preinitialize, init = .weibulls_init,
    initialize = expression({
      if (exists("penalized", inherits = FALSE) && !penalized &&
          scoreType == "NCV" && nei$jackknife > 2 &&
          utils::packageVersion("mgcv") == "1.9.4")
        stop("weibulls NCV/QNCV for a fully unpenalized model on mgcv 1.9-4 requires explicit nei with jackknife = FALSE",
             call. = FALSE)
      n <- rep.int(1, nobs)
      start <- family$init(y, x, E, weights, offset, start)
    }),
    postproc = expression({
      colnames(object$fitted.values) <- c("mu", "k")
      active <- object$prior.weights > 0
      eta <- object$linear.predictors[active, , drop = FALSE]
      ell <- object$family$loglik(object$y[active], eta[, 1L], eta[, 2L])[, 1L]
      object$aic <- 2 * (sum(object$edf) - sum(object$prior.weights[active] * ell))
      object$deviance <- sum(object$family$dev.resids(object$y,
        object$fitted.values, object$prior.weights))
      object$null.deviance <- object$family$null.deviance(object$y,
        object$fitted.values[, 2L], object$prior.weights, G$offset[[1L]],
        attr(G$pterms[[1L]], "intercept") == 1L)
    }),
    predict = .weibulls_predict, residuals = .weibulls_residuals,
    dev.resids = function(y, mu, wt) .weibull_deviance(y, mu[, 1L], mu[, 2L], wt),
    null.deviance = .weibull_null_deviance,
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1L, sandwich = TRUE)$lbb,
    rd = function(mu, wt = 1, scale = 1) .weibull_random(mu),
    cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
      .weibull_probability(q, mu, lower.tail = lower.tail, log.p = logp),
    qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
      .weibull_probability(p, mu, quantile = TRUE, lower.tail = lower.tail, log.p = log.p),
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
  .inlaws_diagnostic_family(family, "weibulls")
}

.weibulls_preinitialize <- function(G) {
  .weibull_response(G$y)
  .weibull_check_weights(G$w)
  jj <- attr(G$X, "lpi")
  if (length(jj) != 2L) stop("weibulls requires two formulae", call. = FALSE)
  if (anyDuplicated(unlist(jj)))
    stop("weibulls requires a separate coefficient block in each predictor", call. = FALSE)
  opt <- get0("optimizer", envir = parent.frame(), inherits = FALSE)
  method <- get0("method", envir = parent.frame(), inherits = FALSE)
  if (length(method) && !method %in% c("REML", "NCV", "GCV.Cp"))
    stop("weibulls supports REML, NCV and QNCV", call. = FALSE)
  scale <- get0("scale", envir = parent.frame(), inherits = FALSE)
  if ((!is.null(scale) && scale != 1) || G$family$scale != 1)
    stop("Weibull likelihood scale is fixed at one", call. = FALSE)
  if (length(opt) && opt[1L] == "efs" && length(G$sp) &&
      ((!is.null(G$L) && !isTRUE(all.equal(unname(G$L), diag(length(G$sp))))) ||
       any(G$lsp0 != 0)))
    stop("weibulls EFS requires all smoothing parameters free and unlinked, or all fixed",
         call. = FALSE)
  # mgcv 1.9-3's EFS path reads this spelling.
  list(offxset = G$offset)
}

.weibulls_eta <- function(X, coef, offset) {
  jj <- attr(X, "lpi")
  eta <- matrix(0, nrow(X), 2L)
  for (j in 1:2) {
    eta[, j] <- drop(X[, jj[[j]], drop = FALSE] %*% coef[jj[[j]]])
    if (length(offset) >= j && length(offset[[j]])) eta[, j] <- eta[, j] + offset[[j]]
  }
  eta
}

.weibulls_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                         d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                         D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("weibulls requires ordinary model matrices", call. = FALSE)
  supplied <- !is.null(eta)
  jj <- attr(X, "lpi")
  if (!supplied) eta <- .weibulls_eta(X, coef, offset)
  w <- rep_len(if (is.null(wt)) 1 else wt, length(y))
  active <- w > 0
  y[!active] <- 1
  eta[!active, ] <- 0
  pars <- exp(eta)
  bad <- function() list(l = -Inf, l0 = ifelse(active, -Inf, 0))
  if (any(!is.finite(pars) | pars <= 0)) return(bad())
  order <- if (!deriv) 0L else if (deriv <= 1L) 2L else if (deriv < 4L) 3L else 4L
  ell <- .weibull_eval(y, eta[, 1L], eta[, 2L], order) * w
  if (any(!is.finite(ell))) return(bad())
  out <- list(l = sum(ell[, 1L]), l0 = ell[, 1L])
  if (deriv) {
    l1 <- ell[, 2:3, drop = FALSE]
    l2 <- ell[, 4:6, drop = FALSE]
    l3 <- if (order >= 3L) ell[, 7:10, drop = FALSE] else 0
    l4 <- if (order >= 4L) ell[, 11:15, drop = FALSE] else 0
    if (!supplied) {
      # The trace-only contraction in mgcv drops dimensions for an
      # intercept-only predictor. Use its full Hessian derivative there.
      trace <- deriv == 2L && any(lengths(jj) == 1L)
      out <- c(out, mgcv::gamlss.gH(X, jj, l1, l2, family$tri$i2,
        l3 = l3, i3 = family$tri$i3, l4 = l4, i4 = family$tri$i4,
        d1b = d1b, d2b = d2b, deriv = if (trace) 2L else deriv - 1L,
        fh = fh, D = D, sandwich = sandwich))
      if (trace) out$d1H <- vapply(out$d1H, function(z) sum(z * fh), numeric(1))
    }
    if (ncv) out <- c(out, list(l1 = l1, l2 = l2, l3 = l3, l4 = l4))
  }
  out
}

.weibulls_init <- function(y, x, E, wt, offset, start) {
  if (!is.null(start)) return(start)
  jj <- attr(x, "lpi")
  if (is.null(wt)) wt <- rep(1, length(y))
  keep <- wt > 0
  off <- lapply(1:2, function(j) if (length(offset) >= j && length(offset[[j]]))
    rep_len(offset[[j]], length(y)) else numeric(length(y)))
  initial <- function(j, target) {
    xx <- x[keep, jj[[j]], drop = FALSE] * sqrt(wt[keep])
    penalty <- E[, jj[[j]], drop = FALSE]
    if (length(penalty) && sum(penalty^2) > 0)
      penalty <- penalty * (0.01 * sqrt(sum(xx^2) / sum(penalty^2)))
    b <- qr.coef(qr(rbind(xx, penalty)), c(
      (target[keep] - off[[j]][keep]) * sqrt(wt[keep]), rep(0, nrow(penalty))))
    b[!is.finite(b)] <- 0
    b
  }
  log_y <- log(y)
  beta <- numeric(ncol(x))
  beta[jj[[1L]]] <- initial(1L, log_y)
  eta <- .weibulls_eta(x, beta, off)
  sd_log <- sqrt(stats::weighted.mean((log_y[keep] - eta[keep, 1L])^2, wt[keep]))
  shape <- max(0.1, min(10, pi / (sqrt(6) * sd_log)))
  beta[jj[[2L]]] <- initial(2L, rep(log(shape), length(y)))
  theta <- .weibulls_eta(x, beta, off)[, 2L]
  # E[log(Y)] = log(mu) - (Euler + H)/k.
  adjustment <- (0.5772156649015329 + .weibull_shift(theta)[, 1L]) / exp(theta)
  beta[jj[[1L]]] <- initial(1L, log_y + adjustment)
  beta
}

.weibulls_predict <- function(family, se = FALSE, eta = NULL, y = NULL,
                              X = NULL, beta = NULL, off = NULL, Vb = NULL) {
  if (is.null(eta)) eta <- .weibulls_eta(X, beta, off) else se <- FALSE
  fit <- exp(eta)
  colnames(fit) <- c("mu", "k")
  if (!se) return(list(fit = fit))
  se_fit <- fit
  for (j in 1:2) {
    ix <- attr(X, "lpi")[[j]]
    xx <- X[, ix, drop = FALSE]
    se_fit[, j] <- fit[, j] * sqrt(pmax(0, rowSums((xx %*% Vb[ix, ix, drop = FALSE]) * xx)))
  }
  list(fit = fit, se.fit = se_fit)
}

.weibulls_residuals <- function(object, type = c("deviance", "pearson", "response",
                                               "scaled.pearson"), ...) {
  mu <- object$fitted.values
  .weibull_residuals(object$y, mu[, 1L], mu[, 2L], object$prior.weights, match.arg(type))
}
