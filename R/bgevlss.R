#' Blended GEV location, spread and shape family
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Fit the Castro-Camilo, Huser and Rue blended generalized extreme value
#' (bGEV) distribution for block maxima using three additive predictors.
#'
#' @param pa,pb Probabilities of the lower and upper blend boundaries.
#'   Require `0 < pa < pb <= min(alpha, beta / 2)`.
#' @param alpha Probability defining the location quantile, in `(0, 1)`.
#' @param beta Probability outside the central interval defining the spread,
#'   in `(0, 1)`.
#' @return An mgcv `general.family` with three predictors, distribution
#'   callbacks and native gratia diagnostics.
#'
#' @details Supply a list of three formulae with the response in the first.
#'   The identity link models the `alpha` quantile. A log link models the
#'   spread, the difference between the `1 - beta / 2` and `beta / 2`
#'   quantiles. A second log link models the positive tail shape `xi`. Use
#'   `~ 1` to estimate a constant parameter. Use separate coefficient blocks.
#'   All three predictors can include offsets. Responses must be finite
#'   numeric vectors. Finite nonnegative case weights multiply log likelihood
#'   contributions, with at least one positive weight.
#'
#'   The blend boundaries are the `pa` and `pb` quantiles of the positive-shape
#'   GEV. A Gumbel distribution matches its CDF at both points. The response
#'   CDF uses the Gumbel to the left, the GEV to the right, and interpolates
#'   the log CDFs with a Beta(5, 5) CDF weight between the boundaries. The
#'   density is the derivative of this CDF. Matching at both ends gives
#'   continuous log-density derivatives through order four. Analytic
#'   predictor derivatives retain the Gumbel limit as shape tends to zero.
#'   A shape near zero can be weakly identified. Extreme predictors or a
#'   blend width below the resolution of the response can cause numerical
#'   failure. The shape has no upper bound.
#'
#'   In [mgcv::gam()], REML uses outer Newton, outer BFGS or extended
#'   Fellner--Schall (`optimizer = "efs"`). EFS uses `gamma = 1` internally
#'   and requires all smoothing parameters free and unlinked, or all fixed.
#'   mgcv passes an unset objective scale to the outer `nlm` and `optim`
#'   optimizers for general families, so those optimizers are unavailable.
#'   Use outer Newton or BFGS when EFS is unsuitable.
#'   NCV and QNCV use outer BFGS with default leave-one-out or explicit
#'   `nei` neighbourhoods. They use a one-step approximation to deletion,
#'   and QNCV evaluates a quadratic approximation to the held-out log
#'   likelihood. Supplying `nei` with REML also selects outer BFGS.
#'   With mgcv 1.9-3 and 1.9-4, upstream jackknife reparameterization fails
#'   for fully unpenalized NCV/QNCV. Supply `nei` with `jackknife = FALSE`
#'   to evaluate the criterion without computing jackknife covariance.
#'
#'   Shape-constrained models in `mgcv::scasm()` were verified on mgcv 1.9-4 and 1.9-5
#'   with REML/EFS and a monotone smooth with binding inequalities. Its
#'   fitter hardwires REML/EFS and provides no NCV, QNCV, Newton or BFGS
#'   smoothing-selection route. Omit `method`, `optimizer` and `nei` from
#'   direct `scasm()` calls. In mgcv 1.9-5 the upstream fitter can fail when
#'   all coefficient directions of a smooth are constrained, or when no
#'   smooth is supplied. `scasm()` is unavailable in mgcv 1.9-3.
#'   `bam()`, `gamm()` and other criteria are outside the supported scope.
#'
#'   Fitted values and response predictions contain `location`, `spread`,
#'   and `shape`, in that order. mgcv can drop response prediction column
#'   names. Link predictions return location, log spread and log shape.
#'   Response prediction standard errors use the delta method and describe
#'   parameter estimation uncertainty.
#'
#'   Use the three-column response-scale parameter matrix `mu` with
#'   `rd(mu, wt = 1, scale = 1)`,
#'   `cdf(q, mu, wt = 1, scale = 1, log_p = FALSE, lower.tail = TRUE)` and
#'   `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)`.
#'   Supply a scalar or one value per row for `q` or `p`. Missing inputs
#'   propagate. These callbacks ignore weights and external scale and accept
#'   shape zero for the limiting Gumbel law. Quantiles are evaluated in closed
#'   form in both tails and by bracketed inversion in the blend.
#'
#'   Native `gratia_diagnostics()` provides probabilities, quantiles, simulation
#'   and the complete response moments. Mean exists for `xi < 1` and is
#'   reported as `NA` otherwise. Variance exists for `xi < 0.5` and is
#'   reported as `Inf` otherwise. Response and Pearson residuals require
#'   their corresponding moments. Deviance and null deviance are undefined.
#'   Use PIT or quantile residuals when moments are unavailable.
#' @references
#'   Castro-Camilo, D., Huser, R. and Rue, H. (2022). Practical strategies for
#'   generalized extreme value-based regression models for extremes.
#'   *Environmetrics*, 33(6), e2742.
#'   [DOI 10.1002/env.2742](https://doi.org/10.1002/env.2742)
#' @seealso [family_diagnostics], [gpdlss()], [mgcv::gam()]
#' @template ncv-neighbourhoods
#' @export
#' @examples
#' set.seed(405)
#' n <- 240
#' d <- data.frame(x = runif(n), z = runif(n), o = runif(n, -.12, .12))
#' f <- bgevlss()
#' d$y <- f$rd(cbind(4 * pmin(d$x, .45) + d$o,
#'                   exp(d$o + .2 * d$z), exp(-1.3 + d$o)))
#' b <- mgcv::gam(list(y ~ s(x, k = 6) + offset(o), ~ z + offset(o),
#'                     ~ offset(o)), data = d, family = f, method = "REML")
#' head(predict(b, type = "response"))
#' head(b$family$qf(.99, fitted(b))) # conditional 0.99 quantile
#' if ("scasm" %in% getNamespaceExports("mgcv")) {
#'   bc <- mgcv::scasm(list(y ~ s(x, bs = "sc", xt = "m+", k = 7) + offset(o),
#'                          ~ z + offset(o), ~ offset(o)), data = d, family = f)
#'   head(fitted(bc))
#' }
bgevlss <- function(pa = 0.05, pb = 0.20, alpha = 0.5, beta = 0.5) {
  scalar_probability <- function(p) is.numeric(p) && length(p) == 1L &&
    is.finite(p) && p > 0 && p < 1
  if (!all(vapply(list(pa, pb, alpha, beta), scalar_probability, logical(1))) ||
      pa >= pb || beta/2 == 0 || 1 - beta/2 == 1)
    stop("bgevlss requires 0 < pa < pb < 1 and alpha and beta in (0, 1)")
  if (pb > min(alpha, beta/2))
    stop("bgevlss requires pb <= min(alpha, beta/2) so the predictors retain their quantile interpretation")
  v <- stats::setNames(-log(-log(c(pa, pb, alpha, beta/2, 1 - beta/2))),
                       c("a", "b", "alpha", "lo", "hi"))
  if (v[["a"]] >= v[["b"]] || v[["lo"]] >= v[["hi"]])
    stop("bgevlss requires distinct blend and spread boundaries at machine precision")
  spec <- list(v = v, d = v[["b"]] - v[["a"]])
  log_link <- stats::make.link("log")
  log_link$linkinv <- log_link$mu.eta <- exp

  cdf <- function(q, mu, wt = 1, scale = 1, log_p = FALSE, lower.tail = TRUE)
    .bgev_probability(q, mu, spec, lower.tail = lower.tail, log.p = log_p)
  qf <- function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
    .bgev_probability(p, mu, spec, quantile = TRUE, lower.tail = lower.tail, log.p = log.p)
  rd <- function(mu, wt = 1, scale = 1) qf(stats::runif(NROW(mu)), mu)

  family <- structure(list(family = "bgevlss", link = c("identity", "log", "log"),
    nlp = 3L, linfo = list(stats::make.link("identity"), log_link, log_link),
    tri = mgcv::trind.generator(3), spec = spec,
    ll = .bgev_ll, ncv = .inlaws_ncv3,
    prepare.response = .bgev_response, init = .bgev_init,
    presetup = expression({
      if (exists("method", inherits = FALSE)) {
        if (missing(method)) method <- "REML"
        if (!method %in% c("REML", "NCV", "QNCV"))
          stop("bgevlss supports REML, NCV and QNCV in gam()")
        if (method %in% c("NCV", "QNCV") || !is.null(nei)) {
          if (identical(optimizer[1L], "efs"))
            stop("bgevlss NCV/QNCV in gam() requires BFGS; EFS uses REML without nei")
          optimizer <- c("outer", "bfgs")
        }
      } else if (any(c("method", "optimizer", "nei") %in% names(cl))) {
        stop("scasm() fits bgevlss with REML/EFS; omit method, optimizer and nei")
      }
    }),
    preinitialize = .bgev_preinitialize,
    initialize = expression({
      y <- family$prepare.response(y)
      if (exists("penalized", inherits = FALSE) && !penalized &&
          scoreType == "NCV" && nei$jackknife > 2 &&
          utils::packageVersion("mgcv") < "1.9.5")
        stop("bgevlss NCV/QNCV jackknife covariance for a fully unpenalized model is unavailable with mgcv 1.9-3 and 1.9-4; supply explicit nei with jackknife = FALSE")
      n <- rep.int(1, nobs)
      start <- family$init(y, x, E, weights, offset, start, family$spec)
    }),
    postproc = expression({
      colnames(object$fitted.values) <- c("location", "spread", "shape")
      object$deviance <- object$null.deviance <- NA_real_
    }),
    predict = .bgev_predict,
    residuals = function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
      type <- match.arg(type)
      .bgev_residuals(object$y, type, object$prior.weights,
                      .bgev_moments(object$fitted.values, spec))
    },
    rd = rd, cdf = cdf, qf = qf,
    gratia_diagnostics = function(mu, wt = 1, scale = 1, y = NULL, eta = NULL) {
      out <- .bgev_moments(mu, spec)
      out$continuous <- TRUE
      out$mean_available <- all(is.finite(out$mean))
      out$variance_available <- all(is.finite(out$variance))
      out$moment_reason <- "Blended GEV means require shape < 1 and variances require shape < 0.5."
      out$cdf <- function(q, log_p = FALSE, lower_tail = TRUE)
        cdf(q, mu, wt, scale, log_p = log_p, lower.tail = lower_tail)
      out$qf <- function(p, log_p = FALSE, lower_tail = TRUE)
        qf(p, mu, wt, scale, log.p = log_p, lower.tail = lower_tail)
      out$rd <- function() rd(mu, wt, scale)
      out$diagnostic_residuals <- function(y, type) .bgev_residuals(y, type, wt, out)
      out
    },
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))

  # scasm has a general-family presetup hook, but no nei argument.
  .inlaws_ncv_family(family)
}

.bgev_response <- function(y) {
  if (!is.numeric(y) || !is.null(dim(y)) || any(!is.finite(y)))
    stop("bgevlss requires a finite numeric response vector")
  y
}

.bgev_preinitialize <- function(G) {
  opt <- get0("optimizer", envir = parent.frame(), inherits = FALSE)
  method <- get0("method", envir = parent.frame(), inherits = FALSE)
  if (length(method) && !method %in% c("REML", "NCV", "QNCV"))
    stop("bgevlss supports REML, NCV and QNCV in gam()")
  if (length(opt) >= 2L && identical(opt[1L], "outer") && length(G$sp) &&
      opt[2L] %in% c("nlm", "optim"))
    stop("mgcv does not supply the required objective scale to nlm or optim for general families; use Newton, BFGS or EFS")
  if (identical(opt[1L], "efs") && length(G$sp) &&
      ((!is.null(G$L) && !isTRUE(all.equal(unname(G$L),
          diag(length(G$sp)), check.attributes = FALSE))) || any(G$lsp0 != 0)))
    stop("bgevlss EFS in gam() requires all smoothing parameters free and unlinked, or all fixed")
  list(y = .bgev_response(G$y), offxset = G$offset)
}

.bgev_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                     d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                     D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  supplied <- !is.null(eta)
  jj <- attr(X, "lpi")
  if (!supplied) eta <- .inlaws_gamlss3_eta(X, coef, offset, family$family)
  wt <- rep_len(if (is.null(wt)) 1 else wt, nrow(eta))
  active <- wt > 0
  eta[!active, ] <- 0
  y[!active] <- 0
  bad <- function() list(l = -Inf, l0 = ifelse(active, -Inf, 0))
  pars <- exp(eta[, 2:3, drop = FALSE])
  if (any(!is.finite(eta)) || any(!is.finite(pars) | pars <= 0)) return(bad())
  order <- if (!deriv) 0L else if (deriv == 1L) 2L else if (deriv < 4L) 3L else 4L
  ans <- suppressWarnings(.bgev_eval(y, eta, order, family$spec))
  ell <- ans$ell * wt
  ell[!active, ] <- 0
  if (any(!is.finite(ell))) return(bad())
  .inlaws_gamlss_jet_ll(ell, ans$powers, X, jj, family, supplied, deriv, d1b, d2b,
                       fh, D, ncv, sandwich)
}

.bgev_init <- function(y, X, E, wt, offset, start, spec) {
  jj <- attr(X, "lpi")
  if (length(jj) != 3L || anyDuplicated(unlist(jj)))
    stop("bgevlss requires three formulae with separate coefficient blocks")
  if (any(!is.finite(wt) | wt < 0) || !is.finite(sum(wt)) || !any(wt > 0))
    stop("bgevlss requires finite nonnegative weights with positive total")
  if (!is.null(start)) {
    if (length(start) != ncol(X) || any(!is.finite(start)))
      stop("invalid bgevlss starting coefficients")
    return(start)
  }
  keep <- wt > 0
  y <- y[keep]
  w <- sqrt(wt[keep])
  off <- lapply(1:3, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]]))
      rep_len(offset[[k]], length(keep))[keep]
    else numeric(length(y))
  })
  solve_start <- function(k, target) {
    x <- X[keep, jj[[k]], drop = FALSE]*w
    e <- E[, jj[[k]], drop = FALSE]
    if (length(e) && sum(e^2) > 0) e <- e*.1*sqrt(sum(x^2)/sum(e^2))
    ans <- qr.coef(qr(rbind(x, e)), c((target - off[[k]])*w, rep(0, nrow(e))))
    ans[is.na(ans)] <- 0
    ans
  }
  start <- numeric(ncol(X))
  start[jj[[1L]]] <- solve_start(1L, y)
  err <- y - drop(X[keep, jj[[1L]], drop = FALSE] %*% start[jj[[1L]]]) - off[[1L]]
  v <- spec$v
  # A central residual spread is less sensitive to the largest maxima.
  err_q <- stats::quantile(err, exp(-exp(-v[c("alpha", "lo", "hi")])), names = FALSE)
  spread <- max(err_q[3L] - err_q[2L], sqrt(.Machine$double.eps))
  start[jj[[1L]]] <- solve_start(1L, y + err_q[1L] - stats::weighted.mean(err, w^2))
  start[jj[[2L]]] <- solve_start(2L, log(spread))
  start[jj[[3L]]] <- solve_start(3L, log(.15))
  start
}

.bgev_predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                          beta = NULL, off = NULL, Vb = NULL) {
  if (is.null(eta)) eta <- .inlaws_gamlss3_eta(X, beta, off, family$family) else se <- FALSE
  fit <- cbind(eta[, 1L], exp(eta[, 2L]), exp(eta[, 3L]))
  colnames(fit) <- c("location", "spread", "shape")
  out <- list(fit = fit)
  if (se) {
    jj <- attr(X, "lpi")
    s <- fit
    multiplier <- cbind(1, fit[, 2:3, drop = FALSE])
    for (k in 1:3) {
      x <- X[, jj[[k]], drop = FALSE]
      s[, k] <- multiplier[, k]*sqrt(pmax(0,
        rowSums((x %*% Vb[jj[[k]], jj[[k]], drop = FALSE])*x)))
    }
    out$se.fit <- s
  }
  out
}

.bgev_residuals <- function(y, type, wt, moments) {
  if (type == "deviance")
    stop("bgevlss deviance residuals are undefined; use PIT, quantile, response or Pearson residuals")
  res <- y - moments$mean
  if (type != "response") {
    wt <- rep_len(wt, length(y))
    res <- ifelse(is.finite(moments$variance), sqrt(wt)*res/sqrt(moments$variance), NA_real_)
    res[wt == 0] <- 0
  }
  res
}
