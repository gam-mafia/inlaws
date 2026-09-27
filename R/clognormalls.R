#' Censored lognormal mean and scale family
#'
#' Model the conditional arithmetic mean and log-response standard deviation
#' with separate additive predictors, using an independent censored lognormal
#' likelihood shared with [clognormal()].
#'
#' @param link Two links, for the arithmetic mean and log-response standard
#'   deviation respectively. Both must be `"log"`.
#' @return A `general.family` object for [mgcv::gam()].
#' @details Supply a list of two formulae: the first specifies the response
#'   and mean predictor; the second is one-sided and specifies the scale
#'   predictor. Both can contain offsets. The model is
#'   `log(Y_i) ~ N(log(mu_i) - sigma_i^2/2, sigma_i^2)`, with
#'   `log(mu_i) = eta_1i` and `log(sigma_i) = eta_2i`. Thus `mu_i` is the
#'   conditional arithmetic mean of the uncensored response, and its variance
#'   is `mu_i^2 * expm1(sigma_i^2)`. `sigma_i` is the standard deviation of
#'   `log(Y_i)`, not of `Y_i`. No lower bound is imposed on sigma; use a modest
#'   scale model, since scale approaching zero can produce singular fits.
#'
#'   The response is supplied on its original positive scale. Use a vector
#'   for exact observations or the two-column censoring matrix described in
#'   [clognormal()]. Exact, left-, right-, and interval-censored observations
#'   can be mixed. Case weights multiply log-likelihood contributions and do
#'   not change the conditional distribution. Missing rows are handled by
#'   `gam()`'s `na.action`.
#'
#'   Smoothing parameter selection supports `method = "REML"` and
#'   `method = "NCV"`, including neighbourhoods supplied through `nei`.
#'   `bam()`, `gamm()`, and ML smoothing selection are not supported. Separate
#'   coefficients are required for the two predictors; shared coefficients
#'   specified with an additional numeric-left-hand-side formula are not
#'   supported. Ordinary covariates or smooths may appear in both formulae.
#'
#'   Fitted values and response predictions have two columns: arithmetic mean
#'   `mu`, then log-response standard deviation `sigma`. Link predictions give
#'   their natural logarithms. Response standard errors use the delta method.
#'   Predictions do not average over coefficient uncertainty or subject random
#'   effects. The distribution components accept the two-column response-scale
#'   matrix as `mu`: `rd(mu, wt = 1, scale = 1)`,
#'   `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)`, and
#'   `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)`.
#'   They describe uncensored responses and ignore `wt` and external `scale`.
#'
#'   Use `~ 1` for a constant scale predictor. This specifies the same response
#'   distribution as [clognormal()], but REML estimates can differ: a scale
#'   intercept is integrated as an unpenalized coefficient here, whereas
#'   `clognormal()` treats sigma as a family parameter. AIC and log likelihood
#'   include original-scale density constants and Jacobians.
#'
#'   Deviance residuals compare the fitted likelihood with the maximum over
#'   each observation's log location, keeping its fitted sigma fixed. Null
#'   deviance fits a constant mean predictor with fitted sigmas and the mean
#'   offset retained; it is a conditional mean-fit diagnostic. Response and
#'   Pearson residuals are `NA` for censored observations.
#' @export
#' @examples
#' set.seed(17)
#' dat <- data.frame(x = runif(250), z = runif(250))
#' mu <- exp(0.5 + sin(2*pi*dat$x))
#' sigma <- exp(-1 + dat$z)
#' y <- rlnorm(250, log(mu) - sigma^2/2, sigma)
#' dat$y <- cbind(pmax(y, 0.4), ifelse(y < 0.4, -Inf, y))
#' fit <- mgcv::gam(list(y ~ s(x, k = 7), ~ s(z, k = 5)), data = dat,
#'                  family = clognormalls(), method = "REML")
#' head(predict(fit, type = "response"))
clognormalls <- function(link = list("log", "log")) {
  if (length(link) != 2L || !all(vapply(link, function(x)
      is.character(x) && length(x) == 1L && !is.na(x) && x == "log", logical(1))))
    stop('clognormalls requires two "log" links')
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
    colnames(object$fitted.values) <- c("mu", "sigma")
    object$null.deviance <- object$family$null.deviance(object, G$offset)
  })
  structure(list(family = "clognormalls", link = unlist(link), nlp = 2L,
    linfo = linfo, tri = mgcv::trind.generator(2), ll = .clnormalls_ll,
    ncv = .clnormalls_ncv, initialize = initialize, init = .clnormalls_init,
    prepare.response = .clnormal_response,
    # General-family fitters count length(y) before initialize is evaluated.
    preinitialize = function(G) list(y = .clnormal_response(G$y), offxset = G$offset),
    postproc = postproc, residuals = .clnormalls_residuals,
    null.deviance = .clnormalls_null_deviance,
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    rd = function(mu, wt = 1, scale = 1)
      stats::rlnorm(nrow(mu), log(mu[, 1]) - mu[, 2]^2/2, mu[, 2]),
    qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
      stats::qlnorm(p, log(mu[, 1]) - mu[, 2]^2/2, mu[, 2],
                    lower.tail = lower.tail, log.p = log.p),
    cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
      stats::plnorm(q, log(mu[, 1]) - mu[, 2]^2/2, mu[, 2],
                    lower.tail = lower.tail, log.p = logp),
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
}

.clnormalls_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                           d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                           D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("clognormalls does not support discrete model matrices")
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
    is.finite(exp(2*eta[, 2])) & exp(2*eta[, 2]) > 0
  if (!all(valid)) return(list(l = -Inf, l0 = rep(-Inf, nrow(eta))))
  order <- if (deriv == 0) 0L else if (deriv == 1) 2L else if (deriv < 4) 3L else 4L
  r <- .clnormal_eval(y, eta[, 1], eta[, 2], order, predictor = TRUE)
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

.clnormalls_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                            Hi = NULL, R = NULL, offset = NULL, dH = NULL,
                            db = NULL, deriv = FALSE, nt = 1) {
  # gamlss.ncv subsets y with y[nei$d], dropping attributes (or flattening
  # matrices). Pass row indices and decode inside its likelihood callback.
  response <- .clnormal_response(y)
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

.clnormalls_init <- function(y, x, E, wt, offset, start) {
  if (is.list(x)) stop("clognormalls requires gam(), not bam()")
  jj <- attr(x, "lpi")
  if (length(jj) != 2L) stop("clognormalls requires two formulae")
  if (anyDuplicated(unlist(jj))) stop("clognormalls does not support shared predictor coefficients")
  if (any(!is.finite(wt) | wt < 0) || !any(wt > 0))
    stop("clognormalls requires finite non-negative weights with positive total")
  if (!is.null(start)) return(start)
  keep <- wt > 0
  off <- lapply(1:2, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]])) rep_len(offset[[k]], length(y)) else
      numeric(length(y))
  })
  # Censor bounds give rough, finite pseudo-responses only for initialization.
  target <- log(as.numeric(y)); bounds <- attr(y, "censor")
  target[bounds == -Inf] <- target[bounds == -Inf] - 0.5
  target[bounds == Inf] <- target[bounds == Inf] + 0.5
  ii <- is.finite(bounds) & bounds != y
  target[ii] <- (target[ii] + log(bounds[ii]))/2
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
  sigma <- max(0.1, min(sigma, 2))
  start[jj[[2]]] <- solve_start(2, rep(log(sigma), length(y)))
  sd <- exp(drop(x[, jj[[2]], drop = FALSE] %*% start[jj[[2]]]) + off[[2]])
  start[jj[[1]]] <- solve_start(1, target + sd^2/2)
  start
}

.clnormalls_deviance <- function(y, mu, wt) {
  r <- .clnormal_eval(y, log(mu[, 1]), log(mu[, 2]), predictor = TRUE)
  d <- 2*wt*(r$sat[, 1] - r$ell[, 1])
  d[wt == 0] <- 0
  pmax(0, d)
}

.clnormalls_residuals <- function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
  type <- match.arg(type)
  y <- .clnormal_response(object$y); bounds <- attr(y, "censor")
  mu <- object$fitted.values[, 1]; sigma <- object$fitted.values[, 2]
  if (type == "deviance") {
    location <- log(mu) - sigma^2/2
    target <- log(as.numeric(y))
    ii <- is.finite(bounds) & bounds != y
    target[ii] <- (target[ii] + log(bounds[ii]))/2
    sign <- sign(target - location)
    sign[bounds == -Inf] <- -1; sign[bounds == Inf] <- 1
    return(sign * sqrt(.clnormalls_deviance(y, object$fitted.values, object$prior.weights)))
  }
  ans <- as.numeric(y) - mu
  if (type != "response") ans <- ans*sqrt(object$prior.weights/(mu^2*expm1(sigma^2)))
  ans[bounds != y] <- NA_real_
  ans
}

.clnormalls_null_deviance <- function(object, offset) {
  y <- .clnormal_response(object$y)
  w <- object$prior.weights
  mu <- object$fitted.values
  off <- if (length(offset) && !is.null(offset[[1]])) rep_len(offset[[1]], length(y)) else
    numeric(length(y))
  objective <- function(a) {
    candidate <- cbind(exp(a + off), mu[, 2])
    if (any(!is.finite(candidate) | candidate <= 0)) return(.Machine$double.xmax/100)
    sum(.clnormalls_deviance(y, candidate, w))
  }
  initial <- stats::weighted.mean(log(mu[, 1]) - off, w)
  opt <- stats::optim(initial, objective, method = "BFGS")
  opt$value
}
