#' Ordered beta regression for bounded continuous responses
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' An extended family for [mgcv::gam()] and [mgcv::bam()] with exact masses
#' at zero and one and a beta density in the interior.
#'
#' @param link Only `"logit"` is supported. The link applies to the mean of
#'   the interior beta component, not to the complete response mean.
#' @param start Optional named list with `phi` (positive finite precision)
#'   and/or `cutpoints` (two ordered numeric cutpoints). These initialize
#'   estimation; they do not fix parameters. Active cutpoints must be finite.
#' @return An `extended.family` object for `gam()` and `bam()`.
#' @details
#' For predictor `eta`, let `mu = plogis(eta)`, with precision `phi` and
#' ordered cutpoints `c1 < c2`. The endpoint probabilities are
#' `p0 = plogis(c1 - eta)` and `p1 = plogis(eta - c2)`. The interior mass
#' is `pI = plogis(eta - c1) - plogis(eta - c2)`; conditional on being
#' interior, the response has beta shapes `mu * phi` and `(1 - mu) * phi`.
#'
#' Responses must be numeric and in `[0,1]`; they are never clipped or
#' automatically rescaled. After missing-value handling, positive-weight
#' observations determine the model: missing zeros fix `c1 = -Inf`, missing
#' ones fix `c2 = Inf`, and missing both endpoints gives ordinary beta
#' regression. With no interior observations precision is unidentified and
#' fitting fails. A reduction is reported in a message. Inactive starting
#' cutpoints are ignored with a message. Prior weights multiply the log
#' likelihood; they are not binomial denominators or beta precision weights.
#'
#' All identifiable nuisance parameters are estimated. The additional GLM
#' scale is fixed at one. `fit$family$getTheta(TRUE)` returns named `phi`,
#' `cut_lower`, and `cut_upper`, including infinite inactive cutpoints.
#' Sparse components, separation, or almost constant interior responses can
#' still produce weak identification or convergence problems.
#'
#' `predict(fit, type = "response")` returns a five-column matrix, in the
#' order `mu`, `p0`, `pI`, `p1`, `mean`, following [zoibetals()] conventions:
#' `mu` is the conditional interior beta mean; `p0`, `pI`, and `p1` are the
#' probabilities of zero, an interior response, and one; `mean = p1 + pI * mu`
#' is the complete response mean. This replaces the previous vector of complete
#' means: select column 5 or use [predict_ordbeta()] with `type = "response"`
#' for that vector. mgcv may discard response column names; the helper restores
#' them before selecting components. With `se.fit = TRUE`, `fit` and `se.fit`
#' have the same five columns. Link predictions remain a vector of `eta`.
#'
#' `fit$fitted.values` and the `mu` arguments of the family distribution
#' methods retain the interior mean. Prediction standard errors use the full
#' coefficient covariance and the delta method, treating fitted precision and
#' active cutpoints as fixed. Disabled endpoint probabilities and their standard
#' errors are exactly zero. `unconditional = TRUE` requests mgcv's smoothing
#' parameter uncertainty correction to the coefficient covariance, when
#' available; it does not add direct uncertainty in precision or cutpoints.
#' These are estimation standard errors, not predictive standard deviations
#' or full joint uncertainty intervals.
#'
#' Response and Pearson residuals use the complete mean and variance.
#' Deviance residuals use the saturated likelihood with nuisance parameters
#' held fixed. The null deviance optimizes a common intercept (when present),
#' retaining offsets and fitted nuisance parameters. The `rd`, `cdf`, and
#' `qf` family callbacks provide response draws, CDFs, and quantiles. Simulation
#' conditions on fitted parameters. The CDF is right-continuous; randomized
#' quantile residuals require randomization across its jumps at the endpoints.
#'
#' `gam()` supports REML, ML, NCV and QNCV. `bam()` supports REML, ML and
#' fREML, and discrete fitting supports fREML and NCV. Ordinary `bam()`
#' stabilizes fitting steps using positive working curvature equal to the larger
#' of the expected and observed information on the predictor scale. Discrete
#' NCV uses the working linear model. Use `scale = 1` in `bam()` so mgcv records
#' the scale as fixed in the fitted model. Fixed nuisance parameters and separate
#' predictors for precision or cutpoints are unavailable.
#' @references
#' Kubinec, R. (2023). Ordered Beta Regression: A Parsimonious, Well-Fitting
#' Model for Continuous Data with Lower and Upper Bounds. Political Analysis,
#' 31, 519--536. doi:10.1017/pan.2022.20.
#' @export
#' @examples
#' set.seed(4)
#' d <- data.frame(x = runif(200))
#' eta <- sin(2 * pi * d$x)
#' mu <- plogis(eta)
#' u <- runif(nrow(d))
#' d$y <- rbeta(nrow(d), mu * 10, (1 - mu) * 10)
#' d$y[u < plogis(-1 - eta)] <- 0
#' d$y[u > 1 - plogis(eta - 1)] <- 1
#' fit <- mgcv::gam(y ~ s(x, k = 6), data = d,
#'                  family = ordbeta(), method = "REML")
#' head(predict(fit, type = "response")) # mu, p0, pI, p1, mean
#' head(predict_ordbeta(fit)) # complete response mean as a vector
#' head(predict_ordbeta(fit, type = "probabilities"))
#' fit$family$getTheta(TRUE)
ordbeta <- function(link = "logit", start = NULL) {
  if (!identical(link, "logit")) stop("ordbeta supports only the logit link")
  if (!is.null(start)) {
    if (!is.list(start) || is.null(names(start)) || anyDuplicated(names(start)) ||
        any(!names(start) %in% c("phi", "cutpoints")))
      stop("start must be a named list containing phi and/or cutpoints")
    if (!is.null(start$phi) && (!is.numeric(start$phi) || length(start$phi) != 1L ||
        !is.finite(start$phi) || start$phi <= 0)) stop("starting phi must be positive and finite")
    if (!is.null(start$cutpoints) && (!is.numeric(start$cutpoints) || length(start$cutpoints) != 2L ||
        anyNA(start$cutpoints) || start$cutpoints[1] >= start$cutpoints[2]))
      stop("starting cutpoints must be two ordered numbers")
  }
  .ordbeta_family("full", c(log(10), -1, log(2)), start)
}

.ordbeta_family <- function(mode, initial, start = NULL) {
  # A configuration belongs to one fitted model. The bam wrapper and fitter
  # both need to see the same final mode and estimated nuisance parameters.
  theta <- initial
  needs_setup <- FALSE
  stabilize_bam <- FALSE
  getTheta <- function(trans = FALSE) {
    if (trans) .ordbeta_natural(theta, mode) else theta
  }
  putTheta <- function(theta) {
    assign("theta", unname(theta), envir = environment(getTheta))
  }
  nt <- length(initial)
  bases <- lapply(0:2, function(level) .ordbeta_basis(nt, level))
  fbasis <- .ordbeta_basis(nt, first = TRUE)
  dev.resids <- function(y, mu, wt, theta = NULL) {
    if (is.null(theta)) theta <- getTheta()
    wt <- rep_len(wt, length(y)); out <- numeric(length(y)); k <- wt > 0
    out[k] <- -2*wt[k]*.ordbeta_loglik(y[k], mu[k], theta, mode)
    out
  }
  Dd <- function(y, mu, theta, wt, level = 0) {
    d <- .ordbeta_Dd(y, mu, theta, wt, level, mode, bases, fbasis)
    if (stabilize_bam && level == 0) {
      # Ordinary bam takes a weighted least-squares step using EDeta2.
      # The Fisher denominator can give an oversize step when local observed
      # curvature is higher. Convert that curvature to mgcv's mu coordinates.
      observed <- d$Dmu2 + d$Dmu * (1 - 2 * mu) / (mu * (1 - mu))
      d$EDmu2 <- pmax(d$EDmu2, observed)
    }
    d
  }
  ls <- function(y, w, theta, scale) {
    if (length(scale) != 1L || !is.finite(scale) || scale != 1)
      stop("ordbeta likelihood scale is fixed at one")
    list(ls = 0, lsth1 = numeric(nt), lsth2 = matrix(0, nt, nt),
         LSTH1 = matrix(0, length(y), nt))
  }
  aic <- function(y, mu, theta = NULL, wt, dev) sum(dev.resids(y, mu, wt, theta))
  variance <- function(mu) {
    p <- .ordbeta_components(stats::qlogis(mu), getTheta(TRUE))
    mean <- p[, 3]+p[, 2]*mu
    p[,1]*mean^2+p[,2]*((mu-mean)^2+mu*(1-mu)/(1+getTheta(TRUE)[1]))+
      p[,3]*(1-mean)^2
  }
  saturated.ll <- function(y, wt, theta = NULL) {
    if (is.null(theta)) theta <- getTheta()
    term <- numeric(length(y)); interior <- which(y > 0 & y < 1 & wt > 0)
    # With fixed nuisance parameters the supremum for either endpoint is zero.
    for (j in interior) {
      objective <- function(eta) .ordbeta_loglik(y[j], stats::plogis(eta), theta, mode)
      opt <- stats::optimize(objective, c(-35, 35), maximum = TRUE, tol = 1e-10)
      if (opt$maximum < -34 || opt$maximum > 34)
        opt <- stats::optimize(objective, c(-700, stats::qlogis(1-.Machine$double.eps/2)),
                               maximum = TRUE, tol = 1e-10)
      term[j] <- opt$objective
    }
    list(f = sum(wt*term), term = term)
  }
  postproc <- function(family, y, prior.weights, fitted, linear.predictors, offset, intercept) {
    sat <- family$saturated.ll(y, prior.weights)
    objective <- function(b) sum(family$dev.resids(y, family$linkinv(b+offset), prior.weights))
    null <- if (intercept) {
      centre <- stats::weighted.mean(linear.predictors-offset, prior.weights)
      width <- 2
      repeat {
        opt <- stats::optimize(objective, centre+c(-width, width))
        if (abs(opt$minimum-centre) < .95*width || width >= 1024) break
        width <- width*2
      }
      opt$objective
    } else objective(0)
    label <- switch(mode, full = "full", lower = "no ones", upper = "no zeros", beta = "interior only")
    list(deviance = max(0, 2*sat$f+sum(family$dev.resids(y, fitted, prior.weights))),
         null.deviance = max(0, 2*sat$f+null),
         family = paste0("Ordered beta (", label, ")"))
  }
  residuals <- function(object, type = c("deviance", "working", "response", "pearson")) {
    type <- match.arg(type)
    if (type == "working") return(object$residuals)
    mu <- object$fitted.values
    p <- .ordbeta_components(object$linear.predictors, object$family$getTheta(TRUE))
    r <- object$y - (p[, 3]+p[, 2]*mu)
    if (type == "response") return(r)
    if (type == "pearson") return(r*sqrt(object$prior.weights/object$family$variance(mu)))
    sat <- object$family$saturated.ll(object$y, object$prior.weights)
    sign(r)*sqrt(pmax(0, 2*sat$term*object$prior.weights+
                       object$family$dev.resids(object$y, mu, object$prior.weights)))
  }
  predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                      beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta) || se) {
      location <- .inlaws_location_prediction(X, beta, off, if (se) Vb)
      if (is.null(eta)) eta <- location$eta
    }
    p <- .ordbeta_components(eta, family$getTheta(TRUE))
    mu <- stats::plogis(eta)
    fit <- cbind(mu = mu, p, mean = p[, 3]+p[, 2]*mu)
    colnames(fit) <- c("mu", "p0", "pI", "p1", "mean")
    if (!se) return(list(fit = fit))
    # Precision and cutpoints are fixed in this coefficient delta method.
    slopes <- .ordbeta_slopes(eta, p)
    slope <- cbind(slopes$conditional, slopes$probabilities, slopes$response)
    colnames(slope) <- colnames(fit)
    list(fit = fit, se.fit = abs(slope) * location$se)
  }
  cdf <- function(q, mu, wt = 1, scale = 1, logp = FALSE) {
    n <- max(length(q), length(mu)); q <- rep_len(q, n); mu <- rep_len(mu, n)
    lp <- .ordbeta_log_components(stats::qlogis(mu), getTheta(TRUE)); phi <- getTheta(TRUE)[1]
    a <- lp[, 1]; b <- lp[, 2]+stats::pbeta(q, mu*phi, (1-mu)*phi, log.p = TRUE)
    top <- pmax(a, b)
    out <- top+log1p(exp(pmin(a, b)-top))
    out[is.infinite(top) & top < 0] <- -Inf
    out[q < 0] <- -Inf; out[q >= 1] <- 0
    out[is.na(q) | is.na(mu)] <- NA_real_
    if (logp) out else exp(out)
  }
  qf <- function(p, mu, wt = 1, scale = 1) {
    n <- max(length(p), length(mu)); p <- rep_len(p, n); mu <- rep_len(mu, n)
    prob <- .ordbeta_components(stats::qlogis(mu), getTheta(TRUE)); phi <- getTheta(TRUE)[1]
    out <- rep(NA_real_, n); valid <- !is.na(p) & p >= 0 & p <= 1
    ok <- valid & !is.na(mu)
    out[ok & p <= prob[, 1]] <- 0
    out[ok & p >= 1-prob[, 3]] <- 1
    mid <- ok & p > prob[, 1] & p < 1-prob[, 3]
    out[mid] <- stats::qbeta((p[mid]-prob[mid,1])/prob[mid,2], mu[mid]*phi, (1-mu[mid])*phi)
    out[ok & p == 0] <- 0; out[ok & p == 1] <- 1
    if (any(!is.na(p) & !valid)) {
      out[!is.na(p) & !valid] <- NaN
      warning("NaNs produced")
    }
    out
  }
  rd <- function(mu, wt = 1, scale = 1) {
    n <- length(mu); p <- .ordbeta_components(stats::qlogis(mu), getTheta(TRUE))
    u <- stats::runif(n); out <- numeric(n)
    upper <- !is.na(mu) & u >= 1-p[, 3]; mid <- !is.na(mu) & u >= p[, 1] & !upper
    out[upper] <- 1; phi <- getTheta(TRUE)[1]
    out[mid] <- stats::rbeta(sum(mid), mu[mid]*phi, (1-mu[mid])*phi)
    out[is.na(mu)] <- NA_real_
    out
  }
  initialize <- expression({
    if (family$needs_setup()) family <- family$configure(family, y, weights, stabilize = TRUE)
    mustart <- (y + .1)/1.2
  })
  configure <- function(family, y, w, stabilize = FALSE) {
    if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) || any(y < 0 | y > 1))
      stop("ordbeta requires finite numeric responses in [0,1]")
    if (any(!is.finite(w)) || any(w < 0) || !any(w > 0) || !is.finite(sum(w)))
      stop("ordbeta requires finite nonnegative weights with positive total weight")
    active <- w > 0; z <- any(y[active] == 0); o <- any(y[active] == 1)
    mid <- y > 0 & y < 1 & active
    if (!any(mid)) stop("ordbeta requires interior observations: precision is unidentified")
    selected <- if (z && o) "full" else if (z) "lower" else if (o) "upper" else "beta"
    if (selected != "full") message("ordbeta: using reduced model (", switch(selected, lower="no ones", upper="no zeros", beta="interior only"), ")")
    m <- max(.Machine$double.eps, min(1-.Machine$double.eps, stats::weighted.mean(y[mid], w[mid])))
    v <- stats::weighted.mean((y[mid]-m)^2, w[mid])
    phi <- if (v > 0) max(.1, min(1000, m*(1-m)/v-1)) else 10
    eta <- stats::qlogis(m)
    total <- sum(w)
    cuts <- c(if (z) eta+stats::qlogis(sum(w[y == 0])/total) else -Inf,
              if (o) eta-stats::qlogis(sum(w[y == 1])/total) else Inf)
    if (!is.null(start$phi)) phi <- start$phi
    if (!is.null(start$cutpoints)) {
      live <- c(z, o)
      if (any(!is.finite(start$cutpoints[live]))) stop("active starting cutpoints must be finite")
      if (any(!live)) message("ordbeta: ignoring inactive starting cutpoints")
      cuts[live] <- start$cutpoints[live]
    }
    if (selected == "full" && !is.finite(cuts[2]-cuts[1]))
      stop("starting cutpoint spacing must be finite")
    init <- c(log(phi), switch(selected, full=c(cuts[1], log(cuts[2]-cuts[1])),
                              lower=cuts[1], upper=cuts[2], beta=numeric()))
    mode <<- selected
    theta <<- unname(init)
    needs_setup <<- FALSE
    stabilize_bam <<- stabilize
    if (nt != length(init)) {
      nt <<- length(init)
      bases <<- lapply(0:2, function(level) .ordbeta_basis(nt, level))
      fbasis <<- .ordbeta_basis(nt, first = TRUE)
    }
    family$mode <- mode
    family$n.theta <- nt
    family$ini.theta <- theta
    family
  }
  preinitialize <- function(y, family) {
    G <- attr(family, "G")
    if (is.null(G)) {
      # Ordinary bam provides the weights only when evaluating initialize.
      needs_setup <<- TRUE
      return(list())
    }
    attr(family, "G") <- NULL
    if (!isTRUE(G$discretize)) family <- mgcv::fix.family.link(family)
    list(family = family$configure(family, y, G$w))
  }
  attr(preinitialize, "needG") <- TRUE
  lk <- stats::make.link("logit")
  family <- structure(c(list(family = "Ordered beta", link = "logit", ordbeta = TRUE,
                   mode = mode, n.theta = nt, ini.theta = initial, scale = 1, no.r.sq = TRUE,
                   getTheta = getTheta, putTheta = putTheta, dev.resids = dev.resids,
                   Dd = Dd, ls = ls, aic = aic, variance = variance,
                   initialize = initialize, preinitialize = preinitialize,
                   needs_setup = function() needs_setup, configure = configure,
                   postproc = postproc, saturated.ll = saturated.ll, residuals = residuals,
                   predict = predict, rd = rd, qf = qf, cdf = cdf,
                   validmu = function(mu) all(is.finite(mu) & mu > 0 & mu < 1)),
              lk[c("linkfun", "linkinv", "mu.eta", "valideta")]),
            class = c("ordbeta.family", "extended.family", "family"))
  .inlaws_diagnostic_family(family, "ordbeta")
}

# Clone before bam captures a local family and passes the same state to its fitter.
#' @export
#' @importFrom mgcv fix.family.link
fix.family.link.ordbeta.family <- function(fam) {
  fam <- unserialize(serialize(fam, NULL))
  NextMethod()
}

.ordbeta_components <- function(eta, par) {
  exp(.ordbeta_log_components(eta, par))
}

.ordbeta_log_components <- function(eta, par) {
  z <- -.ordbeta_softplus(eta-par[2]); o <- -.ordbeta_softplus(par[3]-eta)
  lp <- -.ordbeta_softplus(par[2]-eta)-.ordbeta_softplus(eta-par[3])
  if (all(is.finite(par[2:3]))) lp <- lp+.ordbeta_log1mexp(par[3]-par[2])
  cbind(zero = z, interior = lp, one = o)
}

.ordbeta_slopes <- function(eta, p) {
  mu <- stats::plogis(eta)
  dz <- -p[,1]*(1-p[,1]); do <- p[,3]*(1-p[,3])
  di <- p[,2]*(p[,1]-p[,3])
  list(response = do+di*mu+p[,2]*mu*(1-mu), conditional = mu*(1-mu),
       probabilities = cbind(zero=dz, interior=di, one=do))
}

#' Predict components of an ordered beta GAM
#'
#' @param object A `gam` fitted with [ordbeta()].
#' @param newdata Optional prediction data; `NULL` uses the fitting data.
#' @param type `"response"` for the complete mean, `"conditional"` for the
#'   interior beta mean, or `"probabilities"` for zero/interior/one probabilities.
#' @param se.fit Return delta-method standard errors treating fitted precision
#'   and active cutpoints as fixed.
#' @param ... Passed to [mgcv::predict.gam()], including `unconditional`,
#'   `exclude`, `na.action`, and `block.size`.
#' @return A vector, or a three-column matrix named `p0`, `pI`, `p1` for
#'   probabilities of zero, an interior response, and one, respectively. With
#'   `se.fit = TRUE`, a list with `fit` and `se.fit` of matching dimensions.
#' @details Selects columns from `predict(object, type = "response")`, restoring
#'   column names if mgcv discards them. The standard response matrix is ordered
#'   `mu`, `p0`, `pI`, `p1`, `mean`; the helper's `"conditional"` and `"response"`
#'   types select `mu` and `mean`, respectively. See [ordbeta()] for interpretation.
#'   `unconditional = TRUE` can correct the coefficient covariance for smoothing
#'   parameter uncertainty, but does not add direct uncertainty in fitted
#'   precision or cutpoints.
#' @export
predict_ordbeta <- function(object, newdata = NULL,
                           type = c("response", "conditional", "probabilities"),
                           se.fit = FALSE, ...) {
  if (!inherits(object, "gam") || !isTRUE(object$family$ordbeta))
    stop("object must be a gam fitted with ordbeta()")
  type <- match.arg(type)
  args <- list(object = object, type = "response", se.fit = se.fit, ...)
  if (!is.null(newdata)) args$newdata <- newdata
  pred <- do.call(mgcv::predict.gam, args)
  columns <- switch(type, response = 5L, conditional = 1L, probabilities = 2:4)
  select <- function(x) {
    colnames(x) <- c("mu", "p0", "pI", "p1", "mean")
    x[, columns, drop = length(columns) == 1L]
  }
  if (se.fit) list(fit = select(pred$fit), se.fit = select(pred$se.fit)) else select(pred)
}
