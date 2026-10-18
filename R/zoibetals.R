#' Beta regression with configurable endpoint inflation
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A normalized mixture with separate additive predictors for the interior
#' beta mean, precision, and endpoint membership, fitted by [mgcv::gam()].
#' @param inflation Which endpoints have mass: `"both"` (default), `"zero"`,
#'   or `"one"`.
#' @return A `general.family` object.
#' @details
#' Supply a list of formulas in the order `logit(mu)`, `log(phi)`, then
#' endpoint logits: zero followed by one for `inflation = "both"`, or the
#' selected endpoint for single inflation. Put the response in the first
#' formula; any predictor can be `~ 1`. Offsets and shared coefficients are
#' supported. The interior is the reference component: with endpoint logits
#' `z`, `o`, probabilities are proportional to `exp(z)`, `1`, `exp(o)`.
#' Single inflation uses a binary logit. Disabled endpoints have exactly zero
#' probability. Internally mgcv links are identity links because the likelihood
#' derivatives are evaluated directly on these predictor scales.
#'
#' Each endpoint has its own additive predictor: for example, the third
#' formula `~ z` models `log(p0/pI)` linearly in `z`, while the fourth
#' `~ w` models `log(p1/pI)` linearly in `w`. Smooth terms are also allowed.
#' In the both-endpoint model these are multinomial logits, so each endpoint
#' probability depends on both membership predictors through normalization.
#'
#' For `0 < y < 1` the density is `pI * dbeta(y, mu*phi, (1-mu)*phi)`;
#' endpoint masses are `p0` and `p1`. The complete mean is `p1 + pI*mu`,
#' and variance is `p1 + pI*(mu^2 + mu*(1-mu)/(1+phi)) - mean^2`.
#' `mu` and `phi` describe only the interior beta distribution.
#'
#' Responses must be finite numeric values in `[0,1]`; excluded endpoints
#' are rejected even in zero-weight rows, and observations are never clipped.
#' Missing rows follow `gam()`'s `na.action`. Finite nonnegative prior weights
#' multiply the likelihood, not precision or the individual draw distribution.
#' Each requested component must occur among positive-weight observations;
#' an absent active endpoint is rejected because its intercept has no finite
#' maximum. At least two distinct interior observations are required. Active
#' components are never silently disabled. Sparse components, separation,
#' interpolation of interior responses, or almost constant interior responses
#' can still give unidentified coefficients or unbounded precision.
#'
#' REML uses the extended Fellner--Schall optimizer (`optimizer = "efs"`),
#' selected automatically, with analytic first and second derivatives.
#' Smoothing parameters must be all free and unlinked, or all fixed; partial
#' fixing and shared smoothing IDs are rejected. In mgcv 1.9-3 use
#' `G <- gam(..., fit = FALSE)` followed by `gam(G = G, sp = ...)` for fixed
#' smoothing parameters. Newton/BFGS smoothing optimization, NCV, ML, `bam()`,
#' and `gamm()` are not supported interfaces. Check the returned fitting method:
#' mgcv can reset unsupported method requests to REML with a warning.
#'
#' `fitted()` and link predictions contain the three or four linear predictors.
#' Response predictions have seven columns, in this order: `mu`, `phi`, `p0`,
#' `pI`, `p1`, `mean`, `variance`. mgcv may drop column names; use
#' [predict_zoibetals()] to obtain named components or complete moments.
#' Standard errors use the full coefficient covariance and a multivariate delta
#' method, including covariance between predictors. By default they condition
#' on estimated smoothing parameters; `unconditional = TRUE` requests mgcv's
#' smoothing uncertainty correction when available. EFS does not currently
#' supply that correction. These are estimation standard errors, not predictive
#' standard deviations; no uncertainty from model selection is included.
#'
#' Native `rd`, `cdf`, `qf`, and `gratia_diagnostics` callbacks consume a matrix
#' of link-scale predictors (as returned by `fitted()`), not response columns.
#' Simulation conditions on fitted parameters. CDFs include endpoint jumps;
#' diagnostics supply left limits for randomized PIT and quantile residuals.
#' Response and Pearson residuals use complete moments. Deviance residuals,
#' deviance and null deviance are undefined because the unrestricted saturated
#' beta likelihood is unbounded. Use PIT, quantile, response or Pearson
#' diagnostics. Very small beta shapes can yield rounded endpoint draws in
#' floating-point arithmetic. Mathematical notes and a reproducible comparison
#' with [ordbeta()] are in `inst/maths/zoibetals.md` and
#' `inst/examples/zoibetals.R` in the source package.
#' @export
#' @examples
#' set.seed(18)
#' d <- data.frame(x = runif(250), z = runif(250))
#' f <- zoibetals("zero")
#' eta <- cbind(qlogis(.3 + .3*d$x), log(12), -1 + 2*d$z)
#' d$y <- f$rd(eta)
#' b <- mgcv::gam(list(y ~ s(x, k = 5), ~ 1, ~ z), data = d,
#'                family = f, method = "REML")
#' head(predict_zoibetals(b, type = "all"))
zoibetals <- function(inflation = c("both", "zero", "one")) {
  inflation <- match.arg(inflation)
  endpoints <- switch(inflation, both = c("p0", "p1"), zero = "p0", one = "p1")
  q <- 2L + length(endpoints)
  preinitialize <- function(G) {
    .bi_validate(G$y, G$w, endpoints)
    if (length(attr(G$X, "lpi")) != q)
      stop("zoibetals requires ", q, " formulas: mu, phi, then active endpoint logits")
    if (length(G$sp) > 0 && ((!is.null(G$L) &&
        !isTRUE(all.equal(unname(G$L), diag(length(G$sp)), check.attributes = FALSE))) ||
        any(G$lsp0 != 0)))
      stop("zoibetals EFS requires all smoothing parameters free and unlinked, or all fixed")
    # mgcv 1.9-3/1.9-4's EFS caller uses the misspelled G$offxset.
    list(offxset = G$offset)
  }
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    if (deriv > 1 || ncv) stop("zoibetals supports EFS REML fitting only")
    if (is.list(X)) stop("zoibetals does not support discrete matrices")
    if (is.null(eta)) eta <- .dirichlet_eta(X, coef, offset)
    z <- .bi_derivatives(y, eta, wt, endpoints, family$tri$i2, deriv > 0)
    if (!is.finite(z$l)) return(list(l = -Inf, lb = numeric(ncol(X)),
                                   lbb = matrix(0, ncol(X), ncol(X))))
    out <- if (deriv) mgcv::gamlss.gH(X, attr(X, "lpi"), z$l1, z$l2,
      family$tri$i2, deriv = 0, sandwich = sandwich) else list()
    out$l <- z$l; out$l0 <- z$l0
    out
  }
  initialize <- expression({
    n <- rep.int(1, nobs)
    family$validate(y, weights)
    if (is.null(start)) start <- family$initial(y, x, weights, offset, E)
  })
  predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                      beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta)) eta <- .dirichlet_eta(X, beta, off)
    z <- .bi_prediction(eta, endpoints)
    if (!se) return(list(fit = z$fit))
    jj <- attr(X, "lpi")
    ans <- z$fit * 0
    for (k in seq_len(ncol(ans))) {
      J <- matrix(0, nrow(X), ncol(X))
      for (j in seq_len(q)) J[, jj[[j]]] <- J[, jj[[j]], drop = FALSE] +
        X[, jj[[j]], drop = FALSE] * z$gradient[, k, j]
      ans[, k] <- sqrt(pmax(0, rowSums((J %*% Vb) * J)))
    }
    list(fit = z$fit, se.fit = ans)
  }
  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    if (type == "deviance") stop("Inflated-beta deviance residuals are undefined; use PIT, quantile, response or Pearson residuals")
    d <- .bi_distribution(object$linear.predictors, endpoints)
    r <- object$y - d$mean
    if (type == "pearson") r <- sqrt(object$prior.weights) * r / sqrt(d$variance)
    r
  }
  family <- structure(list(family = paste0("zoibetals (", inflation, ")"),
    zoibetals = TRUE, inflation = inflation, endpoints = endpoints,
    nlp = q, link = rep("identity", q),
    linfo = replicate(q, stats::make.link("identity"), simplify = FALSE),
    tri = mgcv::trind.generator(q), ll = ll, preinitialize = preinitialize,
    initialize = initialize,
    validate = function(y, wt) .bi_validate(y, wt, endpoints),
    initial = function(y, X, wt, off, E) .bi_initial(y, X, wt, off, E, endpoints),
    predict = predict, residuals = residuals,
    postproc = expression({
      object$deviance <- object$null.deviance <- NA_real_
      colnames(object$linear.predictors) <- c("logit.mu", "log.phi",
        paste0("log.", object$family$endpoints, ".over.pI"))
      object$fitted.values <- object$linear.predictors
    }),
    rd = function(mu, wt = 1, scale = 1) .bi_distribution(mu, endpoints)$rd(),
    cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE)
      .bi_distribution(mu, endpoints)$cdf(q, log_p = logp),
    qf = function(p, mu, wt = 1, scale = 1) .bi_distribution(mu, endpoints)$qf(p),
    gratia_diagnostics = function(mu, wt, scale, y, eta) {
      out <- .bi_distribution(if (is.null(eta)) mu else eta, endpoints)
      out$diagnostic_residuals <- function(y, type) {
        if (type == "response") return(y - out$mean)
        if (type == "pearson") return(sqrt(rep_len(wt, length(y))) * (y-out$mean)/sqrt(out$variance))
        stop("Inflated-beta deviance residuals are undefined; use PIT, quantile, response or Pearson residuals")
      }
      out
    },
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
    d2link = 1, d3link = 1, d4link = 1, ls = 1, available.derivs = 0,
    discrete.ok = FALSE), class = c("general.family", "extended.family", "family"))
  family
}

.bi_validate <- function(y, wt, endpoints) {
  if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) || any(y < 0 | y > 1))
    stop("zoibetals requires finite numeric responses in [0,1]")
  if ((!("p0" %in% endpoints) && any(y == 0)) ||
      (!("p1" %in% endpoints) && any(y == 1)))
    stop("response includes an endpoint excluded by inflation")
  if (is.null(wt)) wt <- rep(1, length(y))
  if (any(!is.finite(wt)) || any(wt < 0) || !any(wt > 0) || !is.finite(sum(wt)))
    stop("zoibetals requires finite nonnegative weights with positive total weight")
  mid <- y > 0 & y < 1 & wt > 0
  if (length(unique(y[mid])) < 2L)
    stop("zoibetals requires at least two distinct positive-weight interior observations: mu and phi are otherwise unidentified or unbounded")
  for (p in endpoints) if (!any(y[wt > 0] == if (p == "p0") 0 else 1))
    stop("active endpoint ", p, " has no positive-weight observations; no finite membership intercept exists")
}

.bi_parameters <- function(eta, endpoints) {
  if (!is.matrix(eta) || ncol(eta) != 2L + length(endpoints))
    stop("zoibetals callbacks require a matrix of link-scale predictors")
  z <- cbind(0, eta[, -(1:2), drop = FALSE])
  z <- z - apply(z, 1L, max)
  # Remove one maximal exponential (exactly one) before summing the rest.
  # log1p preserves tiny endpoint odds when the interior dominates.
  ez <- exp(z); ez[cbind(seq_len(nrow(z)), max.col(z, ties.method = "first"))] <- 0
  lp <- z - log1p(rowSums(ez))
  p <- matrix(0, nrow(eta), 3L, dimnames = list(rownames(eta), c("p0", "pI", "p1")))
  logp <- p - Inf
  p[, "pI"] <- exp(lp[, 1]); logp[, "pI"] <- lp[, 1]
  p[, endpoints] <- exp(lp[, -1, drop = FALSE]); logp[, endpoints] <- lp[, -1, drop = FALSE]
  list(mu = stats::plogis(eta[, 1]), phi = exp(eta[, 2]), p = p, logp = logp)
}

.bi_derivatives <- function(y, eta, wt = NULL, endpoints, i2, deriv = TRUE) {
  n <- length(y); q <- ncol(eta)
  if (is.null(wt)) wt <- rep(1, n)
  live <- wt > 0; eta[!live, ] <- 0
  p <- .bi_parameters(eta, endpoints)
  m <- p$mu; ph <- p$phi; a <- m*ph; b <- (1-m)*ph
  if (any(!is.finite(eta)) || any(!is.finite(a) | !is.finite(b) | a <= 0 | b <= 0))
    return(list(l = -Inf))
  component <- ifelse(y == 0, 1L, ifelse(y == 1, 3L, 2L))
  l0 <- p$logp[cbind(seq_len(n), component)]
  ii <- which(live & y > 0 & y < 1)
  l0[ii] <- l0[ii] + stats::dbeta(y[ii], a[ii], b[ii], log = TRUE)
  l0[!live] <- 0; l0 <- l0 * wt
  out <- list(l = sum(l0), l0 = l0)
  if (!deriv) return(out)
  g <- matrix(0, n, q); h <- matrix(0, n, q*(q+1)/2)
  if (length(ii)) {
    mm <- m[ii]; pp <- ph[ii]; aa <- a[ii]; bb <- b[ii]
    t <- log(y[ii])-log1p(-y[ii])-digamma(aa)+digamma(bb)
    lm <- pp*t
    lp <- digamma(pp)-mm*digamma(aa)-(1-mm)*digamma(bb)+
      mm*log(y[ii])+(1-mm)*log1p(-y[ii])
    lmm <- -pp^2*(trigamma(aa)+trigamma(bb))
    lmp <- t+pp*(-mm*trigamma(aa)+(1-mm)*trigamma(bb))
    lpp <- trigamma(pp)-mm^2*trigamma(aa)-(1-mm)^2*trigamma(bb)
    s <- mm*(1-mm)
    g[ii, 1] <- lm*s; g[ii, 2] <- lp*pp
    h[ii, i2[1, 1]] <- lmm*s^2+lm*s*(1-2*mm)
    h[ii, i2[1, 2]] <- lmp*s*pp
    h[ii, i2[2, 2]] <- lpp*pp^2+lp*pp
  }
  for (j in seq_along(endpoints)) {
    pj <- p$p[, endpoints[j]]
    g[, j+2L] <- (y == if (endpoints[j] == "p0") 0 else 1)-pj
    for (k in seq.int(j, length(endpoints)))
      h[, i2[j+2L, k+2L]] <- -pj*((j == k)-p$p[, endpoints[k]])
  }
  g[!live, ] <- 0; h[!live, ] <- 0
  out$l1 <- g*wt; out$l2 <- h*wt
  out
}

.bi_initial <- function(y, X, wt, off, E, endpoints) {
  if (is.null(wt)) wt <- rep(1, length(y))
  mid <- y > 0 & y < 1 & wt > 0
  m <- stats::weighted.mean(y[mid], wt[mid])
  v <- stats::weighted.mean((y[mid]-m)^2, wt[mid])
  phi <- max(.1, min(1000, m*(1-m)/v-1))
  target <- c(stats::qlogis(m), log(phi), vapply(endpoints, function(p)
    log(sum(wt[y == if (p == "p0") 0 else 1])/sum(wt[mid])), numeric(1)))
  jj <- attr(X, "lpi"); n <- length(y); q <- length(jj)
  A <- matrix(0, n*q, ncol(X)); rhs <- numeric(n*q)
  for (j in seq_len(q)) {
    rows <- (j-1L)*n+seq_len(n)
    A[rows, jj[[j]]] <- X[, jj[[j]], drop = FALSE]*sqrt(wt)
    oj <- if (length(off) >= j && !is.null(off[[j]])) off[[j]] else 0
    rhs[rows] <- (target[j]-oj)*sqrt(wt)
  }
  if (nrow(E)) { A <- rbind(A, E); rhs <- c(rhs, numeric(nrow(E))) }
  start <- qr.coef(qr(A), rhs); start[!is.finite(start)] <- 0
  start
}

.bi_prediction <- function(eta, endpoints) {
  z <- .bi_parameters(eta, endpoints); m <- z$mu; ph <- z$phi; p <- z$p
  mn <- p[, "p1"]+p[, "pI"]*m
  vb <- m*(1-m)/(1+ph); e2 <- m^2+vb
  # Nonnegative mixture form avoids cancellation near degenerate distributions.
  vr <- p[, "p0"]*mn^2+p[, "pI"]*(vb+(m-mn)^2)+p[, "p1"]*(1-mn)^2
  fit <- cbind(mu = m, phi = ph, p, mean = mn, variance = vr)
  g <- array(0, c(nrow(eta), 7L, ncol(eta)))
  g[, 1, 1] <- m*(1-m); g[, 2, 2] <- ph
  for (j in seq_len(ncol(eta))) {
    dp <- p*0
    if (j > 2L) for (k in colnames(p))
      dp[, k] <- p[, k]*((k == endpoints[j-2L])-p[, endpoints[j-2L]])
    dm <- g[, 1, j]; dph <- g[, 2, j]
    g[, 3:5, j] <- dp
    dmean <- dp[, "p1"]+dp[, "pI"]*m+p[, "pI"]*dm
    g[, 6, j] <- dmean
    g[, 7, j] <- dp[, "p1"]+dp[, "pI"]*e2+
      p[, "pI"]*((2*m+(1-2*m)/(1+ph))*dm-m*(1-m)/(1+ph)^2*dph)-2*mn*dmean
  }
  list(fit = fit, gradient = g)
}

.bi_logadd <- function(a, b) {
  top <- pmax(a, b)
  out <- top + log1p(exp(pmin(a, b)-top))
  out[is.infinite(top) & top < 0] <- -Inf
  pmin(0, out)
}

.bi_distribution <- function(eta, endpoints) {
  z <- .bi_parameters(eta, endpoints); m <- z$mu; ph <- z$phi; p <- z$p; lp <- z$logp
  n <- length(m); a <- m*ph; b <- (1-m)*ph
  moments <- .bi_prediction(eta, endpoints)$fit
  cdf <- function(q, log_p = FALSE, lower_tail = TRUE) {
    if (!length(q)) return(numeric())
    ix <- rep_len(seq_len(n), max(n, length(q)))
    q <- rep_len(q, length(ix))
    val <- if (lower_tail) .bi_logadd(lp[ix, "p0"], lp[ix, "pI"]+
      stats::pbeta(q, a[ix], b[ix], log.p = TRUE)) else
      .bi_logadd(lp[ix, "p1"], lp[ix, "pI"]+stats::pbeta(q, a[ix], b[ix], lower.tail = FALSE, log.p = TRUE))
    val[!is.na(q) & q < 0] <- if (lower_tail) -Inf else 0
    val[!is.na(q) & q >= 1] <- if (lower_tail) 0 else -Inf
    if (log_p) val else exp(val)
  }
  cdf_left <- function(q, log_p = FALSE, lower_tail = TRUE) {
    if (!length(q)) return(numeric())
    ix <- rep_len(seq_len(n), max(n, length(q)))
    q <- rep_len(q, length(ix))
    val <- cdf(q, log_p = TRUE, lower_tail = lower_tail)
    iz <- which(q == 0); io <- which(q == 1)
    val[iz] <- if (lower_tail) -Inf else 0
    val[io] <- if (lower_tail) .bi_logadd(lp[ix[io], "p0"], lp[ix[io], "pI"]) else lp[ix[io], "p1"]
    if (log_p) val else exp(val)
  }
  qf <- function(p, log_p = FALSE, lower_tail = TRUE) {
    if (!length(p)) return(numeric())
    ix <- rep_len(seq_len(n), max(n, length(p)))
    prob <- rep_len(p, length(ix)); p0 <- z$p[ix, "p0"]; pI <- z$p[ix, "pI"]
    if (log_p) prob <- if (lower_tail) exp(prob) else -expm1(prob)
    else if (!lower_tail) prob <- 1-prob
    out <- rep(NA_real_, length(ix)); ok <- !is.na(prob) & prob >= 0 & prob <= 1
    out[!is.na(prob) & !ok] <- NaN
    if (any(!is.na(prob) & !ok)) warning("NaNs produced")
    iz <- which(ok & prob <= p0)
    io <- which(ok & prob >= p0+pI)
    im <- which(ok & prob > p0 & prob < p0+pI)
    out[iz] <- 0; out[io] <- 1
    out[im] <- stats::qbeta((prob[im]-p0[im])/pI[im], a[ix[im]], b[ix[im]])
    out[which(prob == 0)] <- 0; out[which(prob == 1)] <- 1
    out
  }
  rd <- function() {
    u <- stats::runif(n); out <- numeric(n)
    io <- which(u >= p[, "p0"]+p[, "pI"])
    im <- which(u >= p[, "p0"] & u < p[, "p0"]+p[, "pI"])
    out[io] <- 1; out[im] <- stats::rbeta(length(im), a[im], b[im])
    out
  }
  list(mean = moments[, "mean"], variance = moments[, "variance"],
       cdf = cdf, cdf_left = cdf_left, qf = qf, rd = rd)
}

#' Predict inflated-beta parameters or complete response moments
#'
#' @param object A `gam` fitted with [zoibetals()].
#' @param newdata Optional prediction data.
#' @param type `"response"` for the complete response mean, `"variance"`,
#'   `"parameters"` for interior `mu`, `phi`, `"probabilities"` for `p0`, `pI`,
#'   `p1`, or `"all"` for all seven columns in that order.
#' @param se.fit Return delta-method standard errors.
#' @param ... Passed to [mgcv::predict.gam()], including `unconditional`,
#'   `exclude`, `na.action`, and `block.size`.
#' @return A vector or named matrix; with `se.fit = TRUE`, a list with matching
#'   `fit` and `se.fit`. See [zoibetals()] for uncertainty limitations.
#' @export
predict_zoibetals <- function(object, newdata = NULL,
                                 type = c("response", "variance", "parameters", "probabilities", "all"),
                                 se.fit = FALSE, ...) {
  if (!inherits(object, "gam") || !isTRUE(object$family$zoibetals))
    stop("object must be a gam fitted with zoibetals()")
  type <- match.arg(type)
  args <- list(object = object, type = "response", se.fit = se.fit, ...)
  if (!is.null(newdata)) args$newdata <- newdata
  pred <- do.call(mgcv::predict.gam, args)
  columns <- switch(type, response = 6L, variance = 7L, parameters = 1:2,
                    probabilities = 3:5, all = 1:7)
  select <- function(x) {
    colnames(x) <- c("mu", "phi", "p0", "pI", "p1", "mean", "variance")
    x[, columns, drop = length(columns) == 1L]
  }
  if (se.fit) list(fit = select(pred$fit), se.fit = select(pred$se.fit)) else select(pred)
}
