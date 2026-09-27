#' Conway-Maxwell-Poisson mean and dispersion family
#'
#' Model the exact mean and dispersion of Conway-Maxwell-Poisson counts with
#' two additive predictors using [mgcv::gam()].
#' @param link Two links, for the mean and dispersion respectively. Both must
#'   be `"log"`.
#' @param control Numerical controls as in [cmp()]: `sum_tol`, `mean_tol`,
#'   `max_terms`, and `max_iter`.
#' @return A `general.family` object with simulation, CDF and quantile callbacks.
#' @details Supply two formulae: the response and mean formula first, followed
#'   by a one-sided dispersion formula. Both allow offsets. The predictors
#'   define `mu = exp(eta1)` and `nu = exp(eta2)`. The probability at integer
#'   y >= 0 is proportional to lambda^y / (y!)^nu; lambda is determined by
#'   numerical inversion so that the expectation is exactly mu. Thus both
#'   free distribution parameters are modelled. nu = 1 is Poisson, nu > 1
#'   is underdispersed, and 0 < nu < 1 is overdispersed.
#'
#'   REML supports outer Newton and BFGS optimization. NCV uses BFGS, with
#'   default leave-one-out or neighbourhoods supplied through `nei`.
#'   This family is intended for `gam()` only. ML, QNCV, EFS, `bam()` and
#'   `gamm()` are outside the supported scope. Separate coefficients are
#'   required for the predictors; covariates can appear in both formulae.
#'
#'   Fitted values and response predictions have columns `mu`, `nu`; link
#'   predictions give their logarithms. Response prediction standard errors
#'   use the delta method. The callbacks `rd(mu)`, `cdf(q, mu, logp = FALSE)`
#'   and `qf(p, mu)` take a two-column matrix of response-scale parameters.
#'   They condition on those parameters, independently of prior weights.
#'   A zero mean is a point mass at zero. For positive means, `qf(0, mu)`
#'   is zero and `qf(1, mu)` is infinity. Missing callback inputs return NA.
#'
#'   Responses must be finite nonnegative integers. Weights must be finite
#'   and nonnegative, with at least one positive count having positive weight.
#'   Weights multiply likelihood contributions. Missing rows follow `gam()`'s
#'   `na.action`. The additional likelihood scale is fixed at one.
#'
#'   Use `~ 1` to estimate constant dispersion. Its REML estimate need not
#'   match [cmp()], because dispersion is a regression coefficient here
#'   rather than an extra family parameter. Dispersion is not bounded;
#'   near-deterministic data can lack a finite estimate.
#'
#'   This pure-R implementation uses adaptive summation and analytic Taylor
#'   arithmetic through fourth order. Extreme parameters can exhaust the
#'   numerical limits and produce an error. Distribution callbacks use a
#'   truncated support, with tail accuracy limited by `sum_tol`.
#'   Deviance compares with the saturated mean at each fitted dispersion.
#'   Null deviance fits a constant mean predictor with fitted dispersions
#'   and mean offsets retained. These are conditional mean-fit diagnostics.
#' @seealso [cmp()]
#' @export
#' @examples
#' set.seed(17)
#' dat <- data.frame(x = runif(80))
#' pars <- cbind(exp(1 + dat$x), exp(0.5 + dat$x))
#' dat$y <- cmpls()$rd(pars)
#' fit <- mgcv::gam(list(y ~ x, ~ x), data = dat,
#'                  family = cmpls(), method = "REML")
#' head(predict(fit, type = "response"))
cmpls <- function(link = list("log", "log"), control = list()) {
  control <- .cmp_control(control)
  if (length(link) != 2L || !all(vapply(link, function(x)
      is.character(x) && length(x) == 1L && !is.na(x) && x == "log", logical(1))))
    stop('cmpls requires two "log" links')
  linfo <- lapply(link, function(x) {
    lk <- stats::make.link(x)
    lk$linkinv <- lk$mu.eta <- function(eta) exp(eta)
    lk
  })
  initialize <- expression({
    n <- rep.int(1, nobs)
    start <- family$init(y, x, E, weights, offset, start)
  })
  postproc <- expression({
    colnames(object$fitted.values) <- c("mu", "nu")
    object$deviance <- sum(object$family$dev.resids(object$y, object$fitted.values, object$prior.weights))
    object$null.deviance <- object$family$null.deviance(
      object, G$offset, intercept = attr(G$pterms[[1]], "intercept"))
  })
  structure(list(family = "cmpls", link = unlist(link), nlp = 2L,
    linfo = linfo, tri = mgcv::trind.generator(2), ll = .cmpls_ll, control = control,
    ncv = .cmpls_ncv, initialize = initialize, init = .cmpls_init,
    postproc = postproc, residuals = .cmpls_residuals, predict = .cmpls_predict,
    dev.resids = function(y, mu, wt) .cmpls_deviance(y, mu, wt, control),
    null.deviance = .cmpls_null_deviance,
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    rd = function(mu, wt = NULL, scale = NULL)
      .cmpls_probability(stats::runif(nrow(mu)), mu, control, quantile = TRUE),
    cdf = function(q, mu, wt = NULL, scale = NULL, logp = FALSE) {
      ans <- .cmpls_probability(q, mu, control)
      if (logp) log(ans) else ans
    },
    qf = function(p, mu, wt = NULL, scale = NULL)
      .cmpls_probability(p, mu, control, quantile = TRUE),
    d2link = 1, d3link = 1, d4link = 1, ls = 1,
    no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
}

.cmpls_eta <- function(X, beta, offset) {
  jj <- attr(X, "lpi")
  if (length(jj) != 2L) stop("cmpls requires two formulae")
  if (anyDuplicated(unlist(jj))) stop("cmpls does not support shared predictor coefficients")
  eta <- matrix(0, nrow(X), 2L)
  for (k in 1:2) {
    eta[, k] <- drop(X[, jj[[k]], drop = FALSE] %*% beta[jj[[k]]])
    if (length(offset) >= k && !is.null(offset[[k]])) eta[, k] <- eta[, k] + offset[[k]]
  }
  eta
}

.cmpls_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                      d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                      D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("cmpls requires gam(), not bam()")
  jj <- attr(X, "lpi"); supplied <- !is.null(eta)
  if (!supplied) eta <- .cmpls_eta(X, coef, offset)
  wt <- rep_len(if (is.null(wt)) 1 else wt, nrow(eta))
  active <- wt > 0
  eta[!active, ] <- 0
  y[!active] <- 0
  mu <- exp(eta[, 1]); t <- eta[, 2]
  bad <- function() list(l = -Inf, l0 = ifelse(active, -Inf, 0))
  if (any(!is.finite(eta)) || any(!is.finite(mu) | mu <= 0) ||
      any(!is.finite(exp(t)) | exp(t) <= 0)) return(bad())
  if (!deriv) {
    l0 <- wt * .cmpls_logpmf(y, mu, exp(t), family$control)
    return(list(l = sum(l0), l0 = l0))
  }
  z <- .cmp_moments(mu, exp(t), family$control)
  ell <- .cmp_lljet(y, mu, exp(t), z, predictor = TRUE)$ll
  ell <- ell * wt
  order <- if (deriv == 1) 2L else if (deriv < 4) 3L else 4L
  packed <- lapply(seq_len(order), function(k) {
    cols <- k * (k + 1) / 2 + seq_len(k + 1L)
    sweep(ell[, cols, drop = FALSE], 2, factorial(k:0) * factorial(0:k), `*`)
  })
  l3 <- if (order >= 3) packed[[3]] else 0
  l4 <- if (order >= 4) packed[[4]] else 0
  ret <- list()
  if (!supplied && !is.null(jj)) ret <- mgcv::gamlss.gH(
    X, jj, packed[[1]], packed[[2]], family$tri$i2,
    l3 = l3, i3 = family$tri$i3, l4 = l4, i4 = family$tri$i4,
    d1b = d1b, d2b = d2b, deriv = deriv - 1L, fh = fh, D = D, sandwich = sandwich)
  if (ncv) {
    ret$l1 <- packed[[1]]; ret$l2 <- packed[[2]]; ret$l3 <- l3
    if (order >= 4) ret$l4 <- l4
  }
  ret$l0 <- ell[, 1]; ret$l <- sum(ret$l0)
  ret
}

.cmpls_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                       Hi = NULL, R = NULL, offset = NULL, dH = NULL,
                       db = NULL, deriv = FALSE, nt = 1) {
  if (isTRUE(family$qapprox)) stop("cmpls does not support QNCV")
  if (is.null(wt)) wt <- rep(1, length(y))
  # Obtain the neighbourhood updates with gamma=1, then reconstruct the
  # gamma adjustment using the actual prediction rows and both offsets.
  # This also handles reordered/repeated nei$d and prediction subsets.
  lf <- llf; lf$gamma <- 1
  native <- utils::getFromNamespace("gamlss.ncv", "mgcv")
  ret <- native(X, y, wt, nei, beta, family, lf, H = H, Hi = Hi,
                R = R, offset = offset, dH = dH, db = db, deriv = deriv, nt = nt)
  ecv <- attr(ret$NCV, "eta.cv"); dcv <- attr(ret$NCV, "deta.cv")
  eta <- .cmpls_eta(X, beta, offset)
  base <- family$ll(y, X, beta, wt, family, eta = eta, deriv = 1, ncv = TRUE)
  ix <- nei$d; gamma <- if (is.null(llf$gamma)) 1 else llf$gamma
  cv <- family$ll(y[ix], X[ix, , drop = FALSE], beta, wt[ix], family,
                  eta = ecv, deriv = 1, ncv = TRUE)
  score <- -gamma * cv$l - (1 - gamma) * sum(base$l0[ix])
  grad <- Vg <- NULL
  if (deriv > 0) {
    jj <- attr(X, "lpi"); m <- length(ix)
    g <- matrix(0, m, ncol(db))
    for (j in 1:2) {
      d0 <- X[ix, jj[[j]], drop = FALSE] %*% db[jj[[j]], , drop = FALSE]
      dc <- dcv[seq_len(m) + (j - 1L) * m, , drop = FALSE]
      g <- g - gamma * cv$l1[, j] * dc - (1 - gamma) * base$l1[ix, j] * d0
    }
    grad <- colSums(g); Vg <- crossprod(g)
  }
  attr(score, "eta.cv") <- ecv
  if (deriv != 0) attr(score, "deta.cv") <- dcv
  list(NCV = score, NCV1 = grad, error = ret$error, Vg = Vg)
}

.cmpls_init <- function(y, x, E, wt, offset, start) {
  if (is.list(x)) stop("cmpls requires gam(), not bam()")
  jj <- attr(x, "lpi")
  if (length(jj) != 2L) stop("cmpls requires two formulae")
  if (anyDuplicated(unlist(jj))) stop("cmpls does not support shared predictor coefficients")
  if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y) | y < 0 | y != floor(y)))
    stop("cmpls requires finite nonnegative integer responses")
  if (any(!is.finite(wt) | wt < 0) || !is.finite(sum(wt)) || !any(wt > 0))
    stop("cmpls requires finite nonnegative weights with positive total")
  if (!any(y[wt > 0] > 0)) stop("cmpls cannot fit an all-zero positive-weight response")
  if (!is.null(start)) {
    if (length(start) != ncol(x) || any(!is.finite(start))) stop("invalid cmpls starting coefficients")
    return(start)
  }
  keep <- wt > 0
  off <- lapply(1:2, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]])) rep_len(offset[[k]], length(y)) else
      numeric(length(y))
  })
  solve_start <- function(k, target) {
    xx <- x[keep, jj[[k]], drop = FALSE] * sqrt(wt[keep])
    ee <- E[, jj[[k]], drop = FALSE]
    if (length(ee) && sum(ee^2) > 0) ee <- ee * (0.01 * sqrt(sum(xx^2)/sum(ee^2)))
    ans <- qr.coef(qr(rbind(xx, ee)),
      c((target[keep] - off[[k]][keep])*sqrt(wt[keep]), rep(0, nrow(ee))))
    ans[!is.finite(ans)] <- 0
    ans
  }
  start <- numeric(ncol(x))
  start[jj[[1]]] <- solve_start(1, log(y + 1/6))
  # Start at the Poisson submodel; both fitted predictors remain unbounded.
  start[jj[[2]]] <- solve_start(2, rep(0, length(y)))
  start
}

.cmpls_predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                           beta = NULL, off = NULL, Vb = NULL) {
  if (is.null(eta)) eta <- .cmpls_eta(X, beta, off) else se <- FALSE
  fit <- exp(eta); colnames(fit) <- c("mu", "nu")
  out <- list(fit = fit)
  if (se) {
    jj <- attr(X, "lpi"); s <- fit
    for (j in 1:2) {
      xj <- X[, jj[[j]], drop = FALSE]
      v <- rowSums((xj %*% Vb[jj[[j]], jj[[j]], drop = FALSE]) * xj)
      s[, j] <- fit[, j] * sqrt(pmax(0, v))
    }
    out$se.fit <- s
  }
  out
}

.cmpls_parameters <- function(mu) {
  if (!is.matrix(mu) || !is.numeric(mu) || ncol(mu) != 2L)
    stop("cmpls parameters must be a numeric matrix with mu and nu columns")
  if (any(!is.na(mu[, 1]) & (!is.finite(mu[, 1]) | mu[, 1] < 0)) ||
      any(!is.na(mu[, 2]) & (!is.finite(mu[, 2]) | mu[, 2] <= 0)))
    stop("cmpls requires finite nonnegative means and positive dispersions")
  invisible(NULL)
}

.cmpls_logpmf <- function(y, mu, nu, control) {
  z <- .cmp_moments(mu, nu, control)
  a <- vapply(z, `[[`, numeric(1), "a")
  ya <- y * a; ya[y == 0] <- 0
  ya - nu * lgamma(y + 1) - vapply(z, `[[`, numeric(1), "logZ")
}

.cmpls_probability <- function(q, mu, control, quantile = FALSE) {
  .cmpls_parameters(mu)
  v <- .cmp_recycle(q, mu[, 1], mu[, 2])
  q <- v[[1]]; m <- v[[2]]; nu <- v[[3]]
  if (quantile && any(!is.na(q) & (q < 0 | q > 1)))
    stop("CMP probabilities must lie in [0, 1]", call. = FALSE)
  out <- rep(NA_real_, length(q))
  keep <- which(!is.na(q) & !is.na(m) & !is.na(nu))
  z <- .cmp_moments(m[keep], nu[keep], control)
  for (j in seq_along(keep)) {
    i <- keep[j]; s <- z[[j]]
    if (quantile) {
      if (q[i] == 0 || m[i] == 0) out[i] <- 0 else if (q[i] == 1) out[i] <- Inf else {
        cs <- cumsum(s$p); cs[length(cs)] <- 1
        out[i] <- s$y[which(cs >= q[i])[1]]
      }
    } else out[i] <- min(1, sum(s$p[s$y <= floor(q[i])]))
  }
  out
}

.cmpls_deviance <- function(y, mu, wt, control) {
  wt <- rep_len(wt, length(y))
  out <- numeric(length(y)); active <- wt > 0
  if (any(active)) {
    nu <- mu[active, 2]; yy <- y[active]
    sat <- .cmpls_logpmf(yy, yy, nu, control)
    out[active] <- pmax(0, 2 * wt[active] *
      (sat - .cmpls_logpmf(yy, mu[active, 1], nu, control)))
  }
  out
}
.cmpls_residuals <- function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values[, 1]; nu <- object$fitted.values[, 2]
  r <- object$y - mu
  if (type == "deviance") return(sign(r) * sqrt(object$family$dev.resids(
    object$y, object$fitted.values, object$prior.weights)))
  if (type != "response") {
    active <- object$prior.weights > 0
    r[!active] <- 0
    v <- vapply(.cmp_moments(mu[active], nu[active], object$family$control),
                `[[`, numeric(1), "variance")
    r[active] <- r[active] * sqrt(object$prior.weights[active] / v)
  }
  r
}
.cmpls_null_deviance <- function(object, offset, intercept = attr(object$pterms[[1]], "intercept")) {
  active <- object$prior.weights > 0
  y <- object$y[active]; w <- object$prior.weights[active]
  nu <- object$fitted.values[active, 2]; control <- object$family$control
  off <- if (length(offset) && !is.null(offset[[1]])) rep_len(offset[[1]], length(active))[active] else
    numeric(length(y))
  objective <- function(a) {
    mu <- exp(a + off)
    if (any(!is.finite(mu) | mu <= 0)) return(.Machine$double.xmax/100)
    -sum(w * .cmpls_logpmf(y, mu, nu, control))
  }
  if (isTRUE(intercept == 0)) value <- objective(0) else {
    initial <- stats::weighted.mean(log(object$fitted.values[active, 1]) - off, w)
    width <- 1
    repeat {
      opt <- stats::optimize(objective, c(initial - width, initial + width))
      if (abs(opt$minimum - initial) < .95 * width) break
      width <- width * 2
      if (width > 128) stop("CMP null deviance optimisation failed", call. = FALSE)
    }
    value <- opt$objective
  }
  max(0, 2 * (sum(w * .cmpls_logpmf(y, y, nu, control)) + value))
}
