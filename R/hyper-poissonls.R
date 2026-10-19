#' Hyper-Poisson mean and dispersion family
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Model the exact mean and dispersion of hyper-Poisson counts with two
#' additive predictors using [mgcv::gam()].
#' @param link Two links, for mean and dispersion; both must be `"log"`.
#' @param control Named numerical controls: `sum_tol = 1e-10`,
#'   `mean_tol = 1e-10`, `max_terms = 100000L`, `max_iter = 100L`.
#'   Tolerances must be between zero and one; limits are positive integers.
#' @return A `general.family` with distribution and native gratia diagnostic
#'   callbacks. Fitting does not require gratia.
#' @details Supply two formulae, with the response in the first. The predictors
#'   define `mu = exp(eta1)` and `sigma = exp(eta2)`. Both support offsets and
#'   shared coefficients using [mgcv::formula.gam()] syntax. Use `~ 1` for
#'   constant dispersion. Fitted values and response predictions have columns
#'   `mu`, `sigma` (where mgcv preserves names); link predictions contain their
#'   logarithms. Prediction standard errors use the delta method and the
#'   selected coefficient covariance, including shared coefficients.
#'
#'   The mass at integer y >= 0 is proportional to
#'   `lambda^y * Gamma(sigma) / Gamma(sigma + y)`. The internal rate lambda is
#'   solved numerically so that `E(Y) = mu`; it is not an additional parameter.
#'   sigma < 1 gives underdispersion, sigma = 1 gives Poisson, and sigma > 1
#'   gives overdispersion. This is the HYPERPO2 parameterization of the reference.
#'   Dispersion is not a free variance multiplier: at fixed mean its limiting
#'   distributions range from Bernoulli (mu <= 1) or a shifted Poisson (mu > 1)
#'   as sigma tends to zero, to geometric as sigma tends to infinity. Boundary
#'   dispersion estimates can therefore be weakly identified. Parameters are
#'   not artificially bounded.
#'
#'   REML supports outer Newton, outer BFGS, and extended Fellner--Schall
#'   (`optimizer = "efs"`). NCV and QNCV use BFGS, including default leave-one-out
#'   or explicit `nei` neighbourhoods. EFS selects REML and requires all smoothing
#'   parameters free and unlinked, or all fixed; use an outer optimizer for
#'   partly fixed or linked smoothing parameters. EFS uses gamma = 1 as in mgcv.
#'   ML, `bam()` and `gamm()` are unsupported. On mgcv 1.9-3, fixed smoothing
#'   parameters with a later intercept-only predictor may require
#'   `G <- gam(..., fit = FALSE)` followed by `gam(G = G, sp = ...)`.
#'
#'   Responses must be finite nonnegative integers. Weights must be finite and
#'   nonnegative, with positive total and at least one positive count among
#'   positive-weight observations. Weights multiply the likelihood and all its
#'   derivatives. Missing rows follow `gam()`'s `na.action`. Additional likelihood
#'   scale is fixed at one. On mgcv 1.9-4, fully unpenalized NCV/QNCV models
#'   require explicit `nei` with `jackknife = FALSE` to avoid an upstream
#'   empty-penalty covariance error.
#'
#'   Native `rd(mu, wt, scale)`,
#'   `cdf(q, mu, wt, scale, logp = FALSE, lower.tail = TRUE)`, and
#'   `qf(p, mu, wt, scale, lower.tail = TRUE, log.p = FALSE)` callbacks consume
#'   a two-column response-scale parameter matrix. Probability arguments are
#'   scalar or one per row. Missing inputs propagate; invalid inputs are
#'   rejected. Zero mean is allowed as a point mass at zero in these callbacks.
#'   For positive means, lower-tail quantiles at zero and one are zero and
#'   infinity. Weights and external scale do not alter the simulated distribution.
#'   The native `gratia_diagnostics()` factory supplies response moments,
#'   discrete CDF left limits, probabilities, quantiles, and simulation.
#'
#'   Deviance compares with the saturated mean y, holding dispersion fixed.
#'   Null deviance fits a constant mean predictor, retaining its offsets and
#'   fitted dispersion; without a mean intercept it uses the offset alone.
#'   These diagnostics are conditional on dispersion. Response residuals are
#'   `y - mu`; Pearson residuals use the response variance and square-root weights.
#'
#'   The pure-R implementation uses adaptive summation, safeguarded mean
#'   inversion and analytic Taylor derivatives through fourth order. Requested
#'   tails extend beyond the normalization support when necessary. Extreme
#'   parameters can exhaust the numerical controls and cause an explicit error.
#'   Mathematical details are in `inst/maths/hyper-poissonls.md`.
#' @references Hernandez-Barajas et al. (2026). Hyper-Poisson Regression
#'   Modeling via GAMLSS. The R Journal.
#'   <https://journal.r-project.org/articles/RJ-2026-004/>.
#' @seealso [cmpls()], [family_diagnostics]
#' @template ncv-neighbourhoods
#' @export
#' @examples
#' set.seed(104)
#' d <- data.frame(x = runif(120))
#' f <- hyper_poissonls()
#' d$y <- f$rd(cbind(exp(1 + d$x), exp(-0.3 + d$x)))
#' b <- mgcv::gam(list(y ~ x, ~ x), data = d, family = f, method = "REML")
#' head(predict(b, type = "response"))
#' head(b$family$cdf(d$y, fitted(b)))
hyper_poissonls <- function(link = list("log", "log"), control = list()) {
  control <- .hp_control(control)
  if (length(link) != 2L || !all(vapply(link, function(x)
      is.character(x) && length(x) == 1L && !is.na(x) && x == "log", logical(1))))
    stop('hyper_poissonls requires two "log" links')
  linfo <- lapply(link, function(x) {
    z <- stats::make.link(x)
    z$linkinv <- z$mu.eta <- function(eta) exp(eta)
    z
  })
  family <- structure(list(family = "hyper_poissonls", link = unlist(link),
    nlp = 2L, linfo = linfo, tri = mgcv::trind.generator(2), control = control,
    ll = .hp_ll, ncv = .hp_ncv, init = .hp_init,
    presetup = expression({
      if (missing(method)) method <- "REML"
      if (!method %in% c("REML", "NCV", "QNCV"))
        stop("hyper_poissonls supports REML, NCV and QNCV")
      if (method %in% c("NCV", "QNCV") || !is.null(nei))
        optimizer <- c("outer", "bfgs")
    }),
    preinitialize = .hp_preinitialize,
    initialize = expression({
      if (exists("penalized", inherits = FALSE) && !penalized &&
          scoreType == "NCV" && nei$jackknife > 2 &&
          utils::packageVersion("mgcv") == "1.9.4")
        stop("hyper_poissonls NCV/QNCV jackknife covariance for a fully unpenalized model is unavailable with mgcv 1.9-4; supply explicit nei with jackknife = FALSE, or use REML", call. = FALSE)
      n <- rep.int(1, nobs)
      start <- family$init(y, x, E, weights, offset, start)
    }),
    postproc = expression({
      colnames(object$fitted.values) <- c("mu", "sigma")
      object$deviance <- sum(object$family$dev.resids(object$y,
        object$fitted.values, object$prior.weights))
      object$null.deviance <- object$family$null.deviance(object, G$offset)
    }),
    predict = .hp_predict, residuals = .hp_residuals,
    null.deviance = .hp_null_deviance,
    dev.resids = function(y, mu, wt) .hp_deviance(y, mu, wt, control),
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    rd = function(mu, wt = 1, scale = 1)
      .hp_probability(stats::runif(nrow(mu)), mu, control, quantile = TRUE),
    cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
      .hp_probability(q, mu, control, lower.tail = lower.tail, log.p = logp),
    qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
      .hp_probability(p, mu, control, quantile = TRUE, lower.tail = lower.tail, log.p = log.p),
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
  .inlaws_ncv_family(.inlaws_diagnostic_family(family, "hyper_poissonls"))
}

.hp_preinitialize <- function(G) {
  # This hook runs in estimate.gam(), including for supplied G objects.
  # Read the actual requested optimizer there, rather than capturing one at setup.
  opt <- get0("optimizer", envir = parent.frame(), inherits = FALSE)
  method <- get0("method", envir = parent.frame(), inherits = FALSE)
  if (length(method) && !method %in% c("GCV.Cp", "REML", "NCV", "QNCV"))
    stop("hyper_poissonls supports REML, NCV and QNCV")
  if (length(opt) && opt[1] == "efs" && length(G$sp) &&
      ((!is.null(G$L) && !isTRUE(all.equal(unname(G$L),
        diag(length(G$sp)), check.attributes = FALSE))) || any(G$lsp0 != 0)))
    stop("hyper_poissonls EFS requires all smoothing parameters free and unlinked, or all fixed; use an outer optimizer")
  # mgcv 1.9-3's general-family EFS dispatch misspells this field.
  list(offxset = G$offset)
}

.hp_eta <- function(X, beta, offset) {
  if (!is.matrix(X)) stop("hyper_poissonls requires gam(), not bam()")
  jj <- attr(X, "lpi")
  if (length(jj) != 2L) stop("hyper_poissonls requires two linear predictors")
  eta <- matrix(0, nrow(X), 2L)
  for (j in 1:2) {
    eta[,j] <- drop(X[,jj[[j]],drop=FALSE] %*% beta[jj[[j]]])
    if (length(offset) >= j && !is.null(offset[[j]])) eta[,j] <- eta[,j]+offset[[j]]
  }
  eta
}

.hp_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                   d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                   D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (!is.matrix(X)) stop("hyper_poissonls requires gam(), not bam()")
  supplied <- !is.null(eta); jj <- attr(X,"lpi")
  if (!supplied) eta <- .hp_eta(X,coef,offset)
  wt <- rep_len(if (is.null(wt)) 1 else wt, length(y)); active <- wt > 0
  eta[!active,] <- 0; y[!active] <- 0
  pars <- exp(eta)
  bad <- function() list(l = -Inf, l0 = ifelse(active,-Inf,0), lb = rep(0,ncol(X)),
                        lbb = matrix(0,ncol(X),ncol(X)))
  if (any(!is.finite(eta)) || any(!is.finite(pars) | pars <= 0)) return(bad())
  order <- if (!deriv) 0L else if (deriv == 1L) 2L else if (deriv < 4L) 3L else 4L
  # A likelihood-only trial outside numerical limits must allow the fitter
  # to halve its step. Fail explicitly if an accepted derivative evaluation
  # or a public distribution callback exhausts the same controls.
  z <- if (deriv) .hp_moments(pars[,1],pars[,2],family$control) else
    tryCatch(.hp_moments(pars[,1],pars[,2],family$control), error = function(e) NULL)
  if (is.null(z)) return(bad())
  ell <- .hp_lljet(y,pars[,1],pars[,2],z,order)*wt
  ret <- list()
  if (deriv) {
    packed <- lapply(seq_len(order), function(k) {
      cols <- k*(k+1)/2+seq_len(k+1L)
      sweep(ell[,cols,drop=FALSE],2,factorial(k:0)*factorial(0:k),"*")
    })
    l3 <- if (order >= 3) packed[[3]] else 0
    l4 <- if (order >= 4) packed[[4]] else 0
    shared <- deriv == 4L && anyDuplicated(unlist(jj)) > 0L
    contracted <- deriv == 2L && any(lengths(jj) == 1L)
    ghderiv <- if (shared || contracted) 2L else deriv-1L
    if (!supplied && !is.null(jj)) {
      ret <- mgcv::gamlss.gH(X,jj,packed[[1]],packed[[2]],family$tri$i2,
        l3=l3,i3=family$tri$i3,l4=l4,i4=family$tri$i4,d1b=d1b,d2b=d2b,
        deriv=ghderiv,fh=fh,D=D,sandwich=sandwich)
      # Avoid mgcv's overlapping-column fourth-order trace overwrite.
      if (shared) ret$trHid2H <- .nbls_shared_trace(X,jj,
        list(l3=l3,l4=l4),family$tri,d1b,d2b,fh,D)
      if (contracted) ret$d1H <- vapply(ret$d1H,function(x) sum(x*fh),numeric(1))
    }
    if (ncv) {
      ret$l1 <- packed[[1]]; ret$l2 <- packed[[2]]; ret$l3 <- l3
      if (order >= 4L) ret$l4 <- l4
    }
  }
  ret$l0 <- ell[,1]; ret$l <- sum(ret$l0)
  ret
}

.hp_init <- function(y, x, E, wt, offset, start) {
  if (!is.matrix(x)) stop("hyper_poissonls requires gam(), not bam()")
  jj <- attr(x,"lpi")
  if (length(jj) != 2L) stop("hyper_poissonls requires two linear predictors")
  if (!is.numeric(y) || !is.null(dim(y)) || any(!is.finite(y) | y < 0 | y != floor(y)))
    stop("hyper_poissonls requires finite nonnegative integer responses")
  if (any(!is.finite(wt) | wt < 0) || !is.finite(sum(wt)) || !any(wt > 0))
    stop("hyper_poissonls requires finite nonnegative weights with positive total")
  keep <- wt > 0
  if (!any(y[keep] > 0)) stop("hyper_poissonls cannot fit an all-zero positive-weight response")
  if (!is.null(start)) {
    if (length(start) != ncol(x) || any(!is.finite(start))) stop("invalid hyper_poissonls starting coefficients")
    return(start)
  }
  target <- list(log(y+1/6),numeric(length(y)))
  design <- rhs <- vector("list",2)
  for (j in 1:2) {
    off <- if (length(offset) >= j && !is.null(offset[[j]]))
      rep_len(offset[[j]],length(y)) else numeric(length(y))
    xx <- matrix(0,sum(keep),ncol(x))
    xx[,jj[[j]]] <- x[keep,jj[[j]],drop=FALSE]
    design[[j]] <- xx*sqrt(wt[keep])
    rhs[[j]] <- (target[[j]][keep]-off[keep])*sqrt(wt[keep])
  }
  if (length(E) && sum(E^2) > 0) E <- E*.01*sqrt(sum(x[keep,,drop=FALSE]^2)/sum(E^2))
  ans <- qr.coef(qr(do.call(rbind,c(design,list(E)))),c(unlist(rhs),rep(0,nrow(E))))
  ans[!is.finite(ans)] <- 0
  ans
}

.hp_predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                        beta = NULL, off = NULL, Vb = NULL) {
  if (is.null(eta)) eta <- .hp_eta(X,beta,off) else se <- FALSE
  fit <- exp(eta); colnames(fit) <- c("mu","sigma")
  out <- list(fit=fit)
  if (se) {
    jj <- attr(X,"lpi"); s <- fit
    for (j in 1:2) {
      xx <- X[,jj[[j]],drop=FALSE]
      s[,j] <- fit[,j]*sqrt(pmax(0,rowSums((xx %*% Vb[jj[[j]],jj[[j]],drop=FALSE])*xx)))
    }
    out$se.fit <- s
  }
  out
}

.hp_deviance <- function(y, mu, wt, control) {
  out <- numeric(length(y)); wt <- rep_len(wt,length(y)); active <- which(wt > 0)
  sat <- .hp_moments(y[active],mu[active,2],control)
  fit <- .hp_moments(mu[active,1],mu[active,2],control)
  for (k in seq_along(active)) {
    i <- active[k]
    out[i] <- max(0,2*wt[i]*(.hp_logmass(y[i],sat[[k]])-.hp_logmass(y[i],fit[[k]])))
  }
  out
}

.hp_residuals <- function(object, type = c("deviance","pearson","response","scaled.pearson"), ...) {
  type <- match.arg(type); mu <- object$fitted.values
  r <- object$y-mu[,1]
  if (type == "deviance") return(sign(r)*sqrt(object$family$dev.resids(object$y,mu,object$prior.weights)))
  if (type != "response") {
    active <- object$prior.weights > 0; r[!active] <- 0
    z <- .hp_moments(mu[active,1],mu[active,2],object$family$control)
    v <- vapply(z,`[[`,numeric(1),"variance")
    r[active] <- r[active]*sqrt(object$prior.weights[active]/v)
  }
  r
}

.hp_null_deviance <- function(object, offset) {
  keep <- object$prior.weights > 0; w <- object$prior.weights[keep]
  y <- object$y[keep]; sigma <- object$fitted.values[keep,2]
  off <- if (length(offset) && !is.null(offset[[1]]))
    rep_len(offset[[1]],length(keep))[keep] else numeric(length(y))
  objective <- function(a) {
    mu <- exp(a+off)
    if (any(!is.finite(mu) | mu <= 0)) return(.Machine$double.xmax/100)
    z <- tryCatch(.hp_moments(mu,sigma,object$family$control),error=function(e) NULL)
    if (is.null(z)) return(.Machine$double.xmax/100)
    -sum(w*vapply(seq_along(y),function(i) .hp_logmass(y[i],z[[i]]),numeric(1)))
  }
  intercept <- attr(object$pterms[[1]],"intercept")
  if (isTRUE(intercept == 0)) value <- objective(0) else {
    initial <- stats::weighted.mean(log(object$fitted.values[keep,1])-off,w)
    fit <- stats::optim(initial,objective,method="BFGS",control=list(reltol=1e-10))
    if (fit$convergence) warning("hyper_poissonls null mean optimization did not converge")
    value <- fit$value
  }
  sat <- .hp_moments(y,sigma,object$family$control)
  max(0,2*(value+sum(w*vapply(seq_along(y),function(i) .hp_logmass(y[i],sat[[i]]),numeric(1)))))
}
