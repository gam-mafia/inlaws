#' Student-t location, scale and degrees-of-freedom family
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Three additive predictors model the location, scale and tail weight of a
#' symmetric Student-t response in [mgcv::gam()].
#' @param min.df Finite scalar strictly greater than 2, the lower bound on
#'   degrees of freedom. The default is 3, as in [mgcv::scat()].
#' @return An mgcv `general.family` with three predictors and distribution
#'   and native diagnostic callbacks.
#' @details Supply a list of three formulae, with the response in the first:
#' location `mu` uses the identity link, scale `sigma` uses the log link, and
#' degrees of freedom use `nu = min.df + exp(eta_nu)`. Smooths and offsets are
#' supported in all three predictors. Use `~ 1` to estimate a global parameter;
#' in particular the third formula `~ 1` estimates global degrees of freedom
#' when it has no offset. Separate coefficient blocks are required.
#'
#' The model is `Y = mu + sigma * T`, with `T` Student-t with `nu` degrees of
#' freedom. Its mean is `mu` and variance is `sigma^2 * nu / (nu - 2)`.
#' **Sigma is the Student-t scale, not the response standard deviation.**
#' The finite-variance constraint is structural; degrees of freedom are never
#' capped and the likelihood is never replaced by a Gaussian likelihood.
#' Analytic mixed derivatives through order four and stable gamma-ratio
#' expansions retain finite-df corrections near the Gaussian limit.
#' Large degrees of freedom can be weakly identified, and varying scale and
#' tail weight can be confounded. Extreme predictors can exceed numerical
#' range and weak identification can cause singular fits.
#'
#' REML supports outer Newton and BFGS optimization. NCV and QNCV support
#' BFGS, default leave-one-out and explicit `nei` neighbourhoods, including
#' different deletion and prediction subsets and offsets in all predictors.
#' NCV is a one-step deletion approximation; QNCV additionally uses a quadratic
#' held-out log-likelihood approximation. Neither is exact deletion/refitting.
#' With mgcv 1.9-4 the family uses signed-Hessian conjugate-gradient deletion
#' solves, local score/gradient assembly and local coefficient jackknife solves
#' to avoid upstream general-family NCV issues. The extended-family response
#' attribute-loss restriction does not apply to this ordinary numeric response.
#' These routes are tested on mgcv 1.9-4; fixes in future versions are not assumed.
#' Fully unpenalized NCV/QNCV models require explicit `nei` with
#' `jackknife = FALSE`: mgcv 1.9-4's jackknife reparameterization fails with
#' an empty penalty structure. Criterion evaluation remains supported;
#' requests for the affected covariance route fail explicitly.
#' ML, other criteria, shared coefficients, `bam()` and `gamm()` are
#' outside the supported scope.
#'
#' Responses must be finite numeric vectors. Missing rows and subsets follow
#' `gam()`'s model frame. Finite nonnegative weights multiply likelihood
#' contributions without changing the response distribution; at least one
#' weight must be positive. Zero-weight rows contribute zero derivatives.
#'
#' Fitted values and response predictions contain the distributional parameters
#' `mu`, `sigma`, `nu` in that order; column names are supplied where mgcv
#' preserves them. Link predictions contain `mu`, `log(sigma)`, and
#' `log(nu - min.df)`. Response standard errors use the delta method, with
#' multipliers `1`, `sigma`, and `nu - min.df`. They describe parameter estimation
#' uncertainty, not predictive standard deviations. They condition on smoothing
#' parameters by default; `unconditional = TRUE` requests mgcv's corrected
#' coefficient covariance when available. Predictions do not integrate over
#' coefficient uncertainty or random effects.
#'
#' Distribution callbacks consume a three-column response-scale parameter
#' matrix: `rd(mu, wt = 1, scale = 1)`,
#' `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)`, and
#' `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)`.
#' Probability arguments are scalar or one per parameter row. Missing inputs
#' propagate; invalid parameter matrices are rejected. Weights and external
#' scale are ignored by these response-distribution callbacks.
#'
#' Deviance saturates location at the observed response, keeping fitted scale
#' and df. Deviance residuals have sign `y - mu`; Pearson residuals use the
#' complete response variance and square-root case weights. Null deviance fits
#' a constant location predictor while retaining location offsets, scale and
#' df, or uses the location offset alone when there is no location intercept.
#' Native `gratia_diagnostics()` supplies Student-t probabilities, simulation,
#' and complete response moments without requiring gratia for fitting.
#'
#' Unlike `scatlss()`, [mgcv::scat()] estimates or fixes global scale and df
#' outside the location predictor. Equal parameter values give equal
#' likelihoods, but different parameter treatments need not give equal REML
#' estimates. [mgcv::shash()] is a different, four-parameter sinh-arcsinh
#' distribution with a skewness predictor; `scatlss()` is symmetric.
#' @seealso [mgcv::scat()], [mgcv::shash()], [family_diagnostics]
#' @template ncv-neighbourhoods
#' @export
#' @examples
#' set.seed(46)
#' dat <- data.frame(x = runif(600), z = runif(600))
#' pars <- cbind(sin(6 * dat$x), exp(-0.5 + 0.2 * dat$z),
#'               3 + exp(0.3 + dat$x))
#' dat$y <- scatlss()$rd(pars)
#' fit <- mgcv::gam(list(y ~ s(x, k = 6), ~ z, ~ x),
#'                  data = dat, family = scatlss(), method = "REML")
#' head(predict(fit, type = "response"))
#' # Global df, with location and scale still varying:
#' global <- mgcv::gam(list(y ~ s(x, k = 6), ~ z, ~ 1),
#'                     data = dat, family = scatlss(), method = "REML")
scatlss <- function(min.df = 3) {
  if (!is.numeric(min.df) || length(min.df) != 1L ||
      !is.finite(min.df) || min.df <= 2)
    stop("min.df must be a finite scalar strictly greater than 2")
  force(min.df)
  df_link <- stats::make.link("log")
  df_link$name <- "logshift"
  df_link$linkfun <- function(mu) log(mu - min.df)
  df_link$linkinv <- function(eta) min.df + exp(eta)
  df_link$mu.eta <- function(eta) exp(eta)
  scale_link <- stats::make.link("log")
  scale_link$linkinv <- scale_link$mu.eta <- function(eta) exp(eta)
  family <- structure(list(family = "scatlss", min.df = min.df,
    link = c("identity", "log", "logshift"), nlp = 3L,
    linfo = list(stats::make.link("identity"), scale_link, df_link),
    tri = mgcv::trind.generator(3), ll = .scatlss_ll, ncv = .scatlss_ncv,
    init = .scatlss_init, evaluate = .scatlss_eval,
    initialize = expression({
      y <- family$prepare.response(y)
      if (exists("penalized", inherits = FALSE) && !penalized &&
          scoreType == "NCV" && nei$jackknife > 2 &&
          utils::packageVersion("mgcv") == "1.9.4")
        stop("scatlss NCV/QNCV jackknife covariance for a fully unpenalized model is unavailable with mgcv 1.9-4 due to upstream empty-penalty reparameterization; supply explicit nei with jackknife = FALSE to compute the criterion, or use REML", call. = FALSE)
      n <- rep.int(1, nobs)
      start <- family$init(y, x, E, weights, offset, start, family$min.df)
    }),
    prepare.response = .scatlss_response,
    presetup = expression({
      if (method %in% c("NCV", "QNCV") || !is.null(nei))
        optimizer <- c("outer", "bfgs")
    }),
    preinitialize = function(G) list(y = .scatlss_response(G$y), offxset = G$offset),
    postproc = expression({
      colnames(object$fitted.values) <- c("mu", "sigma", "nu")
      object$deviance <- sum(object$family$dev.resids(object$y, object$fitted.values, object$prior.weights))
      object$null.deviance <- object$family$null.deviance(object, G$offset)
    }),
    predict = .scatlss_predict, residuals = .scatlss_residuals,
    dev.resids = .scatlss_deviance, null.deviance = .scatlss_null_deviance,
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    rd = function(mu, wt = 1, scale = 1) {
      .scatlss_parameters(mu)
      mu[, 1] + mu[, 2]*stats::rt(nrow(mu), mu[, 3])
    },
    cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
      .scatlss_probability(q, mu, lower.tail = lower.tail, log.p = logp),
    qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
      .scatlss_probability(p, mu, quantile = TRUE, lower.tail = lower.tail, log.p = log.p),
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
  .inlaws_ncv_family(.inlaws_diagnostic_family(family, "scatlss"))
}

.scatlss_response <- function(y) {
  if (!is.numeric(y) || !is.null(dim(y)) || any(!is.finite(y)))
    stop("scatlss requires a finite numeric response vector")
  y
}

.scatlss_predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                            beta = NULL, off = NULL, Vb = NULL) {
  if (is.null(eta)) eta <- .scatlss_eta(X, beta, off) else se <- FALSE
  fit <- cbind(eta[, 1], exp(eta[, 2]), family$min.df + exp(eta[, 3]))
  colnames(fit) <- c("mu", "sigma", "nu")
  out <- list(fit = fit)
  if (se) {
    jj <- attr(X, "lpi"); s <- fit
    multiplier <- cbind(1, fit[, 2], exp(eta[, 3]))
    for (k in 1:3) {
      xx <- X[, jj[[k]], drop = FALSE]
      s[, k] <- multiplier[, k]*sqrt(pmax(0,
        rowSums((xx %*% Vb[jj[[k]], jj[[k]], drop = FALSE])*xx)))
    }
    out$se.fit <- s
  }
  out
}

.scatlss_deviance <- function(y, mu, wt) {
  wt <- rep_len(wt, length(y)); keep <- wt > 0; out <- numeric(length(y))
  z <- (y[keep]-mu[keep, 1])/mu[keep, 2]
  nu <- mu[keep, 3]; z2 <- z^2; a <- z2/nu
  small <- a <= .25; dev <- numeric(length(z))
  # Retain tiny deviances even when z^2/nu underflows. In the tails,
  # evaluate the logarithm without requiring a representable squared z.
  ratio <- rep(1, sum(small)); nonzero <- a[small] > 0
  ratio[nonzero] <- log1p(a[small][nonzero])/a[small][nonzero]
  dev[small] <- (1 + 1/nu[small])*z2[small]*ratio
  loga <- 2*log(abs(z[!small])) - log(nu[!small])
  dev[!small] <- (nu[!small]+1)*(loga + log1p(exp(-loga)))
  out[keep] <- wt[keep]*pmax(0, dev)
  out
}

.scatlss_residuals <- function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
  type <- match.arg(type); mu <- object$fitted.values; w <- object$prior.weights
  if (type == "deviance") return(sign(object$y-mu[, 1]) *
    sqrt(.scatlss_deviance(object$y, mu, w)))
  ans <- object$y-mu[, 1]
  if (type != "response") {
    ans <- ans*sqrt(w/.scatlss_variance(mu)); ans[w == 0] <- 0
  }
  ans
}

.scatlss_null_deviance <- function(object, offset) {
  keep <- object$prior.weights > 0; y <- object$y[keep]; w <- object$prior.weights[keep]
  mu <- object$fitted.values[keep, , drop = FALSE]
  off <- if (length(offset) && !is.null(offset[[1]]))
    rep_len(offset[[1]], length(keep))[keep] else numeric(length(y))
  objective <- function(a) {
    mu[, 1] <- a + off
    sum(.scatlss_deviance(y, mu, w))
  }
  if (isTRUE(attr(object$pterms[[1]], "intercept") == 0)) return(objective(0))
  initial <- stats::weighted.mean(mu[, 1]-off, w)
  gradient <- function(a) {
    z <- (y-a-off)/mu[, 2]
    -2*sum(w*(1 + 1/mu[, 3])*z/(mu[, 2]*(1 + z^2/mu[, 3])))
  }
  ans <- stats::optim(initial, objective, gr = gradient, method = "BFGS")
  if (ans$convergence != 0L) warning("scatlss null-deviance optimization did not converge")
  max(0, ans$value)
}

.scatlss_eta <- function(X, beta, offset) {
  if (!is.matrix(X)) stop("scatlss requires ordinary model matrices")
  jj <- attr(X, "lpi")
  if (length(jj) != 3L) stop("scatlss requires three formulae")
  if (anyDuplicated(unlist(jj))) stop("scatlss does not support shared predictor coefficients")
  eta <- matrix(0, nrow(X), 3L)
  for (k in 1:3) {
    eta[, k] <- drop(X[, jj[[k]], drop = FALSE] %*% beta[jj[[k]]])
    if (length(offset) >= k && !is.null(offset[[k]])) eta[, k] <- eta[, k] + offset[[k]]
  }
  eta
}

.scatlss_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                          d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                          D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("scatlss does not support discrete model matrices")
  supplied <- !is.null(eta); jj <- attr(X, "lpi")
  if (!supplied) eta <- .scatlss_eta(X, coef, offset)
  wt <- rep_len(if (is.null(wt)) 1 else wt, nrow(eta)); active <- wt > 0
  eta[!active, ] <- 0; y[!active] <- 0
  bad <- function() list(l = -Inf, l0 = ifelse(active, -Inf, 0))
  parameters <- exp(eta[, 2:3, drop = FALSE])
  if (any(!is.finite(eta)) || any(!is.finite(parameters) | parameters <= 0) ||
      any(!is.finite(family$min.df + parameters[, 2]))) return(bad())
  order <- if (!deriv) 0L else if (deriv == 1L) 2L else if (deriv < 4L) 3L else 4L
  r <- family$evaluate(y, eta, order, family$min.df)
  ell <- r$ell * wt; ell[!active, ] <- 0
  if (any(!is.finite(ell))) return(bad())
  ret <- list()
  if (deriv) {
    # mgcv packs symmetric tensors lexicographically, not expand.grid order.
    packed <- lapply(seq_len(order), function(d) {
      p <- as.matrix(expand.grid(rep(list(1:3), d)))
      p <- p[apply(p, 1, function(x) all(diff(x) >= 0)), , drop = FALSE]
      p <- p[do.call(base::order, as.data.frame(p)), , drop = FALSE]
      ans <- matrix(0, nrow(eta), nrow(p))
      for (k in seq_len(nrow(p))) {
        powers <- tabulate(p[k, ], nbins = 3L)
        ix <- which(rowSums(abs(sweep(r$powers, 2, powers))) == 0L)
        ans[, k] <- ell[, ix]*prod(factorial(powers))
      }
      ans
    })
    l3 <- if (order >= 3L) packed[[3]] else 0
    l4 <- if (order >= 4L) packed[[4]] else 0
    contracted <- deriv == 2L && any(lengths(jj) == 1L)
    # mgcv 1.9-4 drops matrix dimensions in its trace-only route
    # when a predictor has a single coefficient. Contract the full Hessian
    # derivative instead, preserving intercept-only scale/shape predictors.
    gh_deriv <- if (contracted) 2L else deriv - 1L
    if (!supplied && !is.null(jj)) ret <- mgcv::gamlss.gH(
      X, jj, packed[[1]], packed[[2]], family$tri$i2,
      l3 = l3, i3 = family$tri$i3, l4 = l4, i4 = family$tri$i4,
      d1b = d1b, d2b = d2b, deriv = gh_deriv, fh = fh, D = D, sandwich = sandwich)
    if (contracted && !supplied) ret$d1H <- vapply(ret$d1H, function(z) sum(z*fh), numeric(1))
    if (ncv) { ret$l1 <- packed[[1]]; ret$l2 <- packed[[2]]; ret$l3 <- l3
      if (order >= 4L) ret$l4 <- l4 }
  }
  ret$l0 <- ell[, 1]; ret$l <- sum(ret$l0); ret
}

.scatlss_init <- function(y, x, E, wt, offset, start, min.df) {
  if (!is.matrix(x)) stop("scatlss requires gam(), not bam()")
  jj <- attr(x, "lpi")
  if (length(jj) != 3L) stop("scatlss requires three formulae")
  if (anyDuplicated(unlist(jj))) stop("scatlss does not support shared predictor coefficients")
  if (any(!is.finite(wt) | wt < 0) || !is.finite(sum(wt)) || !any(wt > 0))
    stop("scatlss requires finite nonnegative weights with positive total")
  if (!is.null(start)) {
    if (length(start) != ncol(x) || any(!is.finite(start))) stop("invalid scatlss starting coefficients")
    return(start)
  }
  keep <- wt > 0
  off <- lapply(1:3, function(k) if (length(offset) >= k && !is.null(offset[[k]]))
    rep_len(offset[[k]], length(y)) else numeric(length(y)))
  solve_start <- function(k, target) {
    xx <- x[keep, jj[[k]], drop = FALSE]*sqrt(wt[keep])
    ee <- E[, jj[[k]], drop = FALSE]
    if (length(ee) && sum(ee^2) > 0) ee <- ee*.01*sqrt(sum(xx^2)/sum(ee^2))
    ans <- qr.coef(qr(rbind(xx, ee)),
      c((target[keep]-off[[k]][keep])*sqrt(wt[keep]), rep(0, nrow(ee))))
    ans[!is.finite(ans)] <- 0; ans
  }
  start <- numeric(ncol(x)); start[jj[[1]]] <- solve_start(1, y)
  loc <- drop(x[, jj[[1]], drop = FALSE] %*% start[jj[[1]]]) + off[[1]]
  sigma <- sqrt(stats::weighted.mean((y[keep]-loc[keep])^2, wt[keep]))
  sigma <- max(sqrt(.Machine$double.eps), sigma * sqrt(1 - 2/(min.df + max(2, 5-min.df))))
  start[jj[[2]]] <- solve_start(2, rep(log(sigma), length(y)))
  start[jj[[3]]] <- solve_start(3, rep(log(max(2, 5 - min.df)), length(y)))
  start
}

# Three-predictor NCV/QNCV: compiled deletion updates, local score assembly
# for correct offsets, prediction subsets/reordering and gamma adjustment.
.scatlss_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL, Hi = NULL, R = NULL,
                      offset = NULL, dH = NULL, db = NULL, deriv = FALSE, nt = 1) {
  if (is.null(wt)) wt <- rep(1, length(y))
  n <- length(y)
  jj <- attr(X, "lpi")
  m <- length(nei$d)
  off <- lapply(1:3, function(j) {
    if (length(offset) >= j && !is.null(offset[[j]])) {
      rep_len(offset[[j]], n)
    } else {
      numeric(n)
    }
  })
  f <- family
  f$qapprox <- FALSE
  lf <- llf
  lf$gamma <- 1
  native <- utils::getFromNamespace("gamlss.ncv", "mgcv")
  # mgcv 1.9-4 Cholesky deletion updates can disagree with the actual
  # deleted Hessian for three predictors with negative cross-curvature.
  # Its conjugate-gradient route uses the full signed Hessian instead.
  if (!is.null(R)) { H <- crossprod(R); Hi <- chol2inv(R) }
  R <- NULL
  ret <- if (deriv < 0) .scatlss_jackknife(X, nei, beta, lf, H, off) else
    native(X, y, wt, nei, beta, f, lf, H = H, Hi = Hi, R = NULL,
           offset = off, dH = dH, db = db, deriv = deriv, nt = nt)
  cg_iterations <- ret$error
  if (any(!is.finite(cg_iterations) | cg_iterations >= 200L)) stop("scatlss NCV deletion solve failed to converge")
  # The CG helper returns iteration counts, whereas gam.fit5 expects a
  # scalar count of failed Cholesky downdates. Keep these distinct.
  ret$error <- 0L
  ecv <- attr(ret$NCV, "eta.cv")
  dcv <- attr(ret$NCV, "deta.cv")
  if (any(!is.finite(ecv)) || (deriv != 0 && any(!is.finite(dcv))))
    stop("scatlss NCV deletion solve returned non-finite predictions or derivatives")
  eta <- sapply(1:3, function(j) drop(X[, jj[[j]], drop = FALSE] %*% beta[jj[[j]]]) + off[[j]])
  base <- family$ll(y, X, beta, wt, family, eta = eta, deriv = if (deriv > 0) 3 else 1, ncv = TRUE)
  ix <- nei$d
  gamma <- llf$gamma
  if (is.null(gamma)) gamma <- 1
  score0 <- -sum(base$l0[ix])
  delta <- ecv - eta[ix, , drop = FALSE]
  tri <- family$tri
  if (isTRUE(family$qapprox)) {
    Hd <- matrix(0, m, 3)
    for (j in 1:3) {
      for (k in 1:3) {
        Hd[, j] <- Hd[, j] + base$l2[ix, tri$i2[j, k]] * delta[, k]
      }
    }
    score <- score0 - gamma * sum(base$l1[ix, , drop = FALSE] * delta + .5 * Hd * delta)
  } else {
    cv <- family$ll(y[ix], X[ix, , drop = FALSE], beta, wt[ix], family, eta = ecv, deriv = 1, ncv = TRUE)
    score <- -gamma * cv$l + (1 - gamma) * score0
  }
  grad <- Vg <- NULL
  if (deriv > 0) {
    ns <- ncol(db)
    g <- matrix(0, m, ns)
    d0 <- lapply(1:3, function(j) X[ix, jj[[j]], drop = FALSE] %*% db[jj[[j]], , drop = FALSE])
    for (j in 1:3) {
      dc <- dcv[seq_len(m) + (j - 1L) * m, , drop = FALSE]
      g <- g - (1 - gamma) * base$l1[ix, j] * d0[[j]]
      if (isTRUE(family$qapprox)) {
        g <- g - gamma * (base$l1[ix, j] + Hd[, j]) * dc
        td <- numeric(m)
        for (k in 1:3) {
          for (l in 1:3) {
            td <- td + base$l3[ix, tri$i3[j, k, l]] * delta[, k] * delta[, l]
          }
        }
        g <- g - .5 * gamma * td * d0[[j]]
      } else {
        g <- g - gamma * cv$l1[, j] * dc
      }
    }
    grad <- colSums(g)
    Vg <- crossprod(g)
  }
  if (!is.finite(score) || any(!is.finite(grad)))
    stop("scatlss NCV criterion or gradient is non-finite")
  attr(score, "eta.cv") <- ecv
  if (deriv != 0) attr(score, "deta.cv") <- dcv
  list(NCV = score, NCV1 = grad, error = ret$error, Vg = Vg, cg.iterations = cg_iterations)
}


# mgcv's CG helper in 1.9-4 does not initialize its negative-deriv output
# buffer. Assemble coefficient deletion updates directly for jackknife output.
.scatlss_jackknife <- function(X, nei, beta, llf, H, offset) {
  jj <- attr(X, "lpi"); p <- length(beta)
  ecv <- matrix(0, length(nei$d), 3L)
  delta <- matrix(0, length(nei$ma), p)
  previous_a <- previous_d <- 0L
  for (i in seq_along(nei$ma)) {
    a <- nei$a[previous_a + seq_len(nei$ma[i]-previous_a)]
    ix <- previous_d + seq_len(nei$md[i]-previous_d); pred <- nei$d[ix]
    hd <- H; score <- numeric(p); col <- 0L
    for (j in 1:3) {
      xj <- X[a, jj[[j]], drop = FALSE]
      score[jj[[j]]] <- drop(crossprod(xj, llf$l1[a, j]))
      for (k in j:3) {
        col <- col + 1L
        block <- crossprod(xj, X[a, jj[[k]], drop = FALSE]*llf$l2[a, col])
        hd[jj[[j]], jj[[k]]] <- hd[jj[[j]], jj[[k]], drop = FALSE] + block
        if (j != k) hd[jj[[k]], jj[[j]]] <- hd[jj[[k]], jj[[j]], drop = FALSE] + t(block)
      }
    }
    delta[i, ] <- solve(hd, score)
    for (j in 1:3) ecv[ix, j] <- drop(X[pred, jj[[j]], drop = FALSE] %*%
      (beta[jj[[j]]] - delta[i, jj[[j]]])) + offset[[j]][pred]
    previous_a <- nei$ma[i]; previous_d <- nei$md[i]
  }
  score <- 0
  attr(score, "eta.cv") <- ecv; attr(score, "deta.cv") <- delta
  list(NCV = score, error = 0L)
}
