#' Generalized gamma mean, scale and shape family
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Model the conditional arithmetic mean, scale and transformed shape of a
#' positive generalized gamma response with three additive predictors.
#' @param link Three links, all `"log"`, for mean, scale and transformed shape.
#' @return An mgcv `general.family` with distribution and diagnostic callbacks.
#' @details Supply three formulae to [mgcv::gam()], with the response in the
#' first and one-sided scale and shape formulae second and third. The predictors
#' are `log(mu)`, `log(sigma)` and `log(r)`, where `mu = E(Y | x)`, `sigma > 0`
#' is Prentice's scale, and `r` equals `1 + sigma * Q`, constrained to be
#' positive. Thus Prentice's shape is
#' `Q = (r - 1) / sigma`. All three predictors have unrestricted real domains,
#' and the distribution always has a finite mean. The third predictor models
#' transformed shape, not Q: constant r does not imply constant Q when sigma
#' varies. Smoothness penalties apply on the three predictor scales.
#'
#' In Prentice's parameterization, for `Q != 0`, draw `G ~ Gamma(Q^-2, 1)` and
#' set `Y = exp(a + sigma * log(Q^2 * G) / Q)`. The location `a` is adjusted
#' to give arithmetic mean mu. At `Q = 0` (\eqn{r=1}), this is log-normal with
#' `meanlog = log(mu) - sigma^2 / 2` and `sdlog = sigma`. The gamma occurs
#' when `Q = sigma` (\eqn{r=1+\sigma^2}); the Weibull when `Q = 1`
#' (\eqn{r=1+\sigma}). Negative Q permits heavy tails while retaining a finite
#' mean. Variance is finite only for \eqn{r>0.5}. Pearson residuals are `NA`
#' where variance is infinite; probability diagnostics remain available.
#'
#' Fitted values and response predictions have columns `mu`, `sigma`, `r`.
#' Link predictions give their logarithms. Response standard errors use the
#' delta method. Recover Q from response predictions as `(r - 1) / sigma`;
#' uncertainty in Q involves both predictors and their covariance. Predictions
#' condition on fitted coefficients, without integration over random effects
#' or coefficient uncertainty. Use `~ 1` for a constant parameter predictor.
#' To estimate constant scale and Prentice shape Q, use `~ 1` in both the
#' second and third formulae. This estimates constant sigma and r, so
#' `Q = (r - 1) / sigma` is constant too, while the mean can vary.
#'
#' Supports REML (outer Newton or BFGS) and NCV/QNCV (BFGS), including offsets
#' in all predictors and neighbourhoods supplied through `nei`. Separate
#' coefficient blocks are required. Shared coefficients, ML smoothing selection
#' and `gamm()` are outside the supported scope. mgcv currently rejects general
#' families in `bam()`, so neither ordinary nor discrete `bam()` is supported.
#'
#' Responses must be finite positive numeric vectors; censoring and zero
#' responses are not supported. Missing rows follow `gam()`'s `na.action`.
#' Finite nonnegative case weights multiply log likelihoods without changing
#' the response distribution; at least one weight must be positive. Extreme
#' predictors can exceed floating-point range and weakly identified scale or
#' shape can cause singular fits. Stable analytic fourth-order derivatives are
#' used through the log-normal limit; no fitted shape bound is imposed beyond
#' the finite-mean constraint.
#'
#' Distribution callbacks take a three-column response-scale matrix in the
#' same `mu`, `sigma`, `r` order: `rd(mu, wt = 1, scale = 1)`,
#' `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)`, and
#' `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)`.
#' Probability arguments are scalar or one per parameter row. Missing inputs
#' propagate; invalid parameter matrices are rejected. Weights and external
#' scale are ignored. Simulation returns one positive draw per parameter row.
#'
#' Deviance maximizes each observation's likelihood over mean, holding sigma
#' and r fixed. Its residual sign compares log response with Prentice location.
#' Null deviance retains fitted sigma and r and mean offsets, fitting a constant
#' mean predictor when that predictor includes an intercept. Log likelihood and
#' AIC include the response Jacobian and density constants.
#' @references Prentice, R. L. (1974). A log gamma model and its maximum
#' likelihood estimation. Biometrika 61, 539-544.
#' @export
#' @examples
#' set.seed(17)
#' dat <- data.frame(x = runif(250), z = runif(250))
#' pars <- cbind(exp(0.5 + sin(6 * dat$x)), exp(-0.7 + 0.2 * dat$z), 1.2)
#' dat$y <- gen_gammalss()$rd(pars)
#' fit <- mgcv::gam(list(y ~ s(x, k = 6), ~ z, ~ 1), data = dat,
#'                  family = gen_gammalss(), method = "REML")
#' head(predict(fit, type = "response"))
#'
#' # Estimate constant scale and shape, allowing the mean to vary:
#' fit_const <- mgcv::gam(list(y ~ s(x, k = 6), ~ 1, ~ 1), data = dat,
#'                        family = gen_gammalss(), method = "REML")
#' pars_const <- fitted(fit_const)
#' # Recover Prentice's Q from response-scale sigma and r:
#' Q <- (pars_const[, "r"] - 1) / pars_const[, "sigma"]
#' head(cbind(pars_const, Q = Q))
gen_gammalss <- function(link = list("log", "log", "log")) {
  if (!identical(unname(unlist(link)), rep("log", 3L)))
    stop('gen_gammalss requires three "log" links')
  linfo <- lapply(link, function(x) {
    lk <- stats::make.link(x)
    lk$linkinv <- lk$mu.eta <- function(eta) exp(eta)
    lk
  })
  family <- structure(list(family = "gen_gammalss", link = unlist(link), nlp = 3L,
    linfo = linfo, tri = mgcv::trind.generator(3), ll = .gen_gamma_ll,
    ncv = .gen_gamma_ncv, init = .gen_gamma_init, evaluate = .gen_gamma_eval,
    initialize = expression({
      y <- family$prepare.response(y)
      n <- rep.int(1, nobs)
      start <- family$init(y, x, E, weights, offset, start)
    }),
    prepare.response = .gen_gamma_response,
    preinitialize = function(G) list(y = .gen_gamma_response(G$y), offxset = G$offset),
    postproc = expression({
      colnames(object$fitted.values) <- c("mu", "sigma", "r")
      object$deviance <- sum(object$family$dev.resids(object$y, object$fitted.values, object$prior.weights))
      object$null.deviance <- object$family$null.deviance(object, G$offset)
    }),
    residuals = .gen_gamma_residuals, predict = .gen_gamma_predict,
    dev.resids = .gen_gamma_deviance, null.deviance = .gen_gamma_null_deviance,
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    rd = function(mu, wt = 1, scale = 1) {
      .gen_gamma_parameters(mu)
      .gen_gamma_probability(stats::runif(nrow(mu)), mu, quantile = TRUE)
    },
    cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
      .gen_gamma_probability(q, mu, lower.tail = lower.tail, log.p = logp),
    qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
      .gen_gamma_probability(p, mu, quantile = TRUE, lower.tail = lower.tail, log.p = log.p),
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
  .inlaws_diagnostic_family(family, "gen_gammalss")
}

.gen_gamma_response <- function(y) {
  if (!is.numeric(y) || !is.null(dim(y)) || any(!is.finite(y) | y <= 0))
    stop("gen_gammalss requires a finite positive numeric response vector")
  y
}

.gen_gamma_eta <- function(X, beta, offset) {
  if (!is.matrix(X)) stop("gen_gammalss requires ordinary model matrices")
  jj <- attr(X, "lpi")
  if (length(jj) != 3L) stop("gen_gammalss requires three formulae")
  if (anyDuplicated(unlist(jj))) stop("gen_gammalss does not support shared predictor coefficients")
  eta <- matrix(0, nrow(X), 3L)
  for (k in 1:3) {
    eta[, k] <- drop(X[, jj[[k]], drop = FALSE] %*% beta[jj[[k]]])
    if (length(offset) >= k && !is.null(offset[[k]])) eta[, k] <- eta[, k] + offset[[k]]
  }
  eta
}

.gen_gamma_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                          d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                          D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("gen_gammalss does not support discrete model matrices")
  supplied <- !is.null(eta); jj <- attr(X, "lpi")
  if (!supplied) eta <- .gen_gamma_eta(X, coef, offset)
  wt <- rep_len(if (is.null(wt)) 1 else wt, nrow(eta)); active <- wt > 0
  eta[!active, ] <- 0; y[!active] <- 1
  bad <- function() list(l = -Inf, l0 = ifelse(active, -Inf, 0))
  if (any(!is.finite(eta)) || any(!is.finite(exp(eta)) | exp(eta) <= 0)) return(bad())
  order <- if (!deriv) 0L else if (deriv == 1L) 2L else if (deriv < 4L) 3L else 4L
  r <- family$evaluate(y, eta, order)
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

.gen_gamma_init <- function(y, x, E, wt, offset, start) {
  if (!is.matrix(x)) stop("gen_gammalss requires gam(), not bam()")
  jj <- attr(x, "lpi")
  if (length(jj) != 3L) stop("gen_gammalss requires three formulae")
  if (anyDuplicated(unlist(jj))) stop("gen_gammalss does not support shared predictor coefficients")
  if (any(!is.finite(wt) | wt < 0) || !is.finite(sum(wt)) || !any(wt > 0))
    stop("gen_gammalss requires finite nonnegative weights with positive total")
  if (!is.null(start)) {
    if (length(start) != ncol(x) || any(!is.finite(start))) stop("invalid gen_gammalss starting coefficients")
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
  start <- numeric(ncol(x)); start[jj[[1]]] <- solve_start(1, log(y))
  loc <- drop(x[, jj[[1]], drop = FALSE] %*% start[jj[[1]]]) + off[[1]]
  sigma <- sqrt(stats::weighted.mean((log(y[keep])-loc[keep])^2, wt[keep]))
  sigma <- max(.05, min(sigma, 3)) # Initialization only.
  start[jj[[2]]] <- solve_start(2, rep(log(sigma), length(y)))
  start[jj[[3]]] <- solve_start(3, numeric(length(y)))
  start[jj[[1]]] <- solve_start(1, log(y) + sigma^2/2)
  start
}

.gen_gamma_predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                               beta = NULL, off = NULL, Vb = NULL) {
  if (is.null(eta)) eta <- .gen_gamma_eta(X, beta, off) else se <- FALSE
  fit <- exp(eta); colnames(fit) <- c("mu", "sigma", "r")
  out <- list(fit = fit)
  if (se) {
    jj <- attr(X, "lpi"); s <- fit
    for (k in 1:3) {
      xx <- X[, jj[[k]], drop = FALSE]
      s[, k] <- fit[, k]*sqrt(pmax(0, rowSums((xx %*% Vb[jj[[k]], jj[[k]], drop = FALSE])*xx)))
    }
    out$se.fit <- s
  }
  out
}

.gen_gamma_deviance <- function(y, mu, wt) {
  wt <- rep_len(wt, length(y)); keep <- wt > 0; out <- numeric(length(y))
  if (any(keep)) {
    r <- .gen_gamma_eval(y[keep], log(mu[keep, , drop = FALSE]))
    out[keep] <- pmax(0, 2*wt[keep]*(r$sat[, 1] - r$ell[, 1]))
  }
  out
}
.gen_gamma_residuals <- function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
  type <- match.arg(type); mu <- object$fitted.values; w <- object$prior.weights
  if (type == "deviance") {
    a <- .gen_gamma_distribution(mu)$location
    return(sign(log(object$y)-a)*sqrt(.gen_gamma_deviance(object$y, mu, w)))
  }
  ans <- object$y-mu[, 1]
  if (type != "response") {
    v <- .gen_gamma_variance(mu); ans <- ans*sqrt(w/v); ans[!is.finite(v)] <- NA_real_
    ans[w == 0] <- 0
  }
  ans
}
.gen_gamma_null_deviance <- function(object, offset) {
  keep <- object$prior.weights > 0; y <- object$y[keep]; w <- object$prior.weights[keep]
  eta <- log(object$fitted.values[keep, , drop = FALSE])
  off <- if (length(offset) && !is.null(offset[[1]]))
    rep_len(offset[[1]], length(keep))[keep] else numeric(length(y))
  objective <- function(a) {
    eta[, 1] <- a + off
    r <- .gen_gamma_eval(y, eta)
    v <- sum(2*w*(r$sat[, 1]-r$ell[, 1]))
    if (is.finite(v)) v else .Machine$double.xmax/100
  }
  if (isTRUE(attr(object$pterms[[1]], "intercept") == 0)) return(objective(0))
  initial <- stats::weighted.mean(eta[, 1]-off, w)
  max(0, stats::optim(initial, objective, method = "BFGS")$value)
}

# Three-predictor NCV/QNCV: compiled deletion updates, local score assembly
# for correct offsets, prediction subsets/reordering and gamma adjustment.
.gen_gamma_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL, Hi = NULL, R = NULL,
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
  ret <- if (deriv < 0) .gen_gamma_jackknife(X, nei, beta, lf, H, off) else
    native(X, y, wt, nei, beta, f, lf, H = H, Hi = Hi, R = NULL,
           offset = off, dH = dH, db = db, deriv = deriv, nt = nt)
  cg_iterations <- ret$error
  if (any(cg_iterations >= 200L)) stop("gen_gammalss NCV deletion solve failed to converge")
  # The CG helper returns iteration counts, whereas gam.fit5 expects a
  # scalar count of failed Cholesky downdates. Keep these distinct.
  ret$error <- 0L
  ecv <- attr(ret$NCV, "eta.cv")
  dcv <- attr(ret$NCV, "deta.cv")
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
  attr(score, "eta.cv") <- ecv
  if (deriv != 0) attr(score, "deta.cv") <- dcv
  list(NCV = score, NCV1 = grad, error = ret$error, Vg = Vg, cg.iterations = cg_iterations)
}


# mgcv's CG helper in 1.9-4 does not initialize its negative-deriv output
# buffer. Assemble coefficient deletion updates directly for jackknife output.
.gen_gamma_jackknife <- function(X, nei, beta, llf, H, offset) {
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
