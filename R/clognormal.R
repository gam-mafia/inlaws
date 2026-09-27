# Independent lognormal likelihood. Taylor coefficients are derivatives divided
# by factorials, in (mu, log(sigma)); multiplication is polynomial convolution.
# No numerical differencing or censored-normal family implementation is used.
.clnormal_jets <- function(order, theta_order = 2L) {
  powers <- expand.grid(a = 0:order, b = 0:min(theta_order, order))
  powers <- as.matrix(powers[rowSums(powers) <= order, ])
  nc <- nrow(powers)
  pairs <- lapply(seq_len(nc), function(k) {
    z <- which(outer(powers[, 1], powers[, 1], "+") == powers[k, 1] &
               outer(powers[, 2], powers[, 2], "+") == powers[k, 2], arr.ind = TRUE)
    z
  })
  constant <- function(x) cbind(x, matrix(0, length(x), nc - 1L))
  multiply <- function(x, y) {
    out <- x * 0
    for (k in seq_len(nc)) {
      z <- pairs[[k]]
      out[, k] <- rowSums(x[, z[, 1], drop = FALSE] * y[, z[, 2], drop = FALSE])
    }
    out
  }
  compose <- function(x, coefficients) {
    dx <- x; dx[, 1] <- 0
    out <- constant(coefficients[[1]])
    power <- constant(rep(1, nrow(x)))
    if (order) for (k in seq_len(order)) {
      power <- multiply(power, dx)
      out <- out + power * coefficients[[k + 1L]]
    }
    out
  }
  exponential <- function(x) compose(x, lapply(0:order, function(k) exp(x[, 1])/factorial(k)))
  logarithm <- function(x) compose(x, c(list(log(x[, 1])),
    lapply(seq_len(order), function(k) (-1)^(k - 1)/(k * x[, 1]^k))))
  list(powers = powers, constant = constant, multiply = multiply,
       compose = compose, exp = exponential, log = logarithm)
}

.clnormal_response <- function(y) {
  if (!is.numeric(y)) stop("clognormal requires a numeric response")
  if (is.matrix(y)) {
    if (ncol(y) != 2L) stop("censored responses must have exactly two columns")
    bound <- y[, 2]; y <- y[, 1]
  } else {
    bound <- attr(y, "censor")
    if (is.null(bound)) bound <- y
  }
  if (length(bound) != length(y) || anyNA(y) || anyNA(bound))
    stop("missing responses or censoring bounds must be removed by na.action")
  if (any(!is.finite(y) | y <= 0))
    stop("responses and first-column censoring limits must be finite and positive")
  interval <- is.finite(bound) & bound != y
  if (any(interval & bound <= y))
    stop("finite interval upper limits must exceed positive lower limits")
  attr(y, "censor") <- bound
  y
}

# Log probabilities and their Taylor coefficients use scaled CDF derivatives.
# Narrow intervals use Gaussian quadrature on the log-response interval to
# avoid subtraction of almost equal CDFs (including their derivatives).
# In predictor mode, mu is log(mu), theta can vary by row, and all mixed
# derivatives through the requested total order are retained.
.clnormal_eval <- function(y, mu, theta, order = 0L, predictor = FALSE) {
  y <- .clnormal_response(y)
  bound <- attr(y, "censor")
  n <- length(y)
  mu <- rep_len(mu, n)
  j <- .clnormal_jets(order, if (predictor) order else 2L)
  C <- j$constant; M <- j$multiply
  t <- C(rep_len(theta, n))
  u <- C(mu)
  if (order) {
    t[, which(j$powers[, 1] == 0 & j$powers[, 2] == 1)] <- 1
    u[, which(j$powers[, 1] == 1 & j$powers[, 2] == 0)] <- 1
  }
  invsd <- j$exp(-t)
  location <- (if (predictor) u else j$log(u)) - j$exp(2*t)/2
  ly <- log(y)
  upper <- ly
  span <- rep(0, n)
  interval <- is.finite(bound) & bound != y
  span[interval] <- log1p((bound[interval] - y[interval])/y[interval])
  wide <- interval & !is.finite(span)
  span[wide] <- log(bound[wide]) - ly[wide]
  upper[interval] <- ly[interval] + span[interval]
  exact <- bound == y
  left <- bound == -Inf
  right <- bound == Inf
  # Integral over a short log interval, with the constant log density factored
  # out before exponentiation. Twelve nodes are ample under the width rule.
  quadrature <- function(lo, width, loc, ti, si) {
    q <- statmod::gauss.quad(12L, kind = "legendre")
    half <- width/2; mid <- lo + half
    zmid <- M(C(mid) - loc, si)
    ref <- -M(zmid, zmid)/2 - ti - C(rep(log(2*pi)/2, length(lo)))
    total <- C(rep(0, length(lo)))
    for (k in seq_along(q$nodes)) {
      z <- M(C(mid + half*q$nodes[k]) - loc, si)
      logd <- -M(z, z)/2 - ti - C(rep(log(2*pi)/2, length(lo)))
      total <- total + q$weights[k] * j$exp(logd - ref)
    }
    j$log(total) + ref + C(log(half))
  }
  probability <- function(lo, hi, loc, ti, si, width) {
    z0 <- M(C(lo) - loc, si); z1 <- M(C(hi) - loc, si)
    a <- z0[, 1]; b <- z1[, 1]
    # Use survival probabilities in the upper tail.
    useupper <- a > 0
    big <- stats::pnorm(b, log.p = TRUE)
    small <- stats::pnorm(a, log.p = TRUE)
    big[useupper] <- stats::pnorm(a[useupper], lower.tail = FALSE, log.p = TRUE)
    small[useupper] <- stats::pnorm(b[useupper], lower.tail = FALSE, log.p = TRUE)
    lp <- big + log(-expm1(small - big))
    out <- C(lp)
    narrow <- is.finite(a) & is.finite(b) &
      width * si[, 1] * (1 + abs((a + b)/2)) < 0.25
    regular <- !narrow
    if (order && any(regular)) {
      logcdf <- function(z) {
        v <- z[, 1]
        out <- C(stats::pnorm(v, log.p = TRUE))
        tail <- is.finite(v) & v < -12
        middle <- is.finite(v) & !tail
        if (any(middle)) {
          zz <- z[middle, , drop = FALSE]; vv <- v[middle]
          ratio <- exp(stats::dnorm(vv, log = TRUE) - out[middle, 1])
          hermite <- list(rep(1, length(vv)), -vv, vv*vv - 1, 3*vv - vv^3)
          coeff <- c(list(rep(1, length(vv))), lapply(seq_len(order),
            function(k) ratio * hermite[[k]]/factorial(k)))
          out[middle, ] <- j$log(j$compose(zz, coeff)) + C(out[middle, 1])
        }
        if (any(tail)) {
          # Mills' asymptotic series avoids cancellation of large terms in
          # high derivatives of log Phi. At z < -12, 16 terms give double
          # precision accuracy; compute the small correction as a jet.
          zz <- z[tail, , drop = FALSE]
          lz <- j$log(-zz)
          iz2 <- j$exp(-2*lz)
          term <- total <- C(rep(1, sum(tail)))
          for (k in 1:16) {
            term <- -(2*k - 1)*M(term, iz2)
            total <- total + term
          }
          out[tail, ] <- -M(zz, zz)/2 - C(rep(log(2*pi)/2, sum(tail))) - lz + j$log(total)
        }
        out
      }
      za <- z0[regular, , drop = FALSE]; zb <- z1[regular, , drop = FALSE]
      flip <- useupper[regular]
      za0 <- za[flip, , drop = FALSE]
      za[flip, ] <- -zb[flip, , drop = FALSE]
      zb[flip, ] <- -za0
      lb <- logcdf(zb); la <- logcdf(za)
      finite <- is.finite(la[, 1])
      result <- lb
      if (any(finite)) {
        delta <- la[finite, , drop = FALSE] - lb[finite, , drop = FALSE]
        difference <- -j$exp(delta)
        difference[, 1] <- -expm1(delta[, 1])
        result[finite, ] <- lb[finite, , drop = FALSE] + j$log(difference)
      }
      out[regular, ] <- result
    }
    if (any(narrow)) out[narrow, ] <- quadrature(lo[narrow], width[narrow],
      loc[narrow, , drop = FALSE], ti[narrow, , drop = FALSE], si[narrow, , drop = FALSE])
    out
  }
  ell <- sat <- C(rep(0, n))
  if (any(exact)) {
    z <- M(C(ly[exact]) - location[exact, , drop = FALSE], invsd[exact, , drop = FALSE])
    sat[exact, ] <- -t[exact, , drop = FALSE] - C(ly[exact] + log(2*pi)/2)
    ell[exact, ] <- sat[exact, , drop = FALSE] - M(z, z)/2
  }
  cens <- !exact
  if (any(cens)) {
    lo <- ly; hi <- upper
    lo[left] <- -Inf; hi[right] <- Inf
    width <- span; width[left | right] <- Inf
    ell[cens, ] <- probability(lo[cens], hi[cens], location[cens, , drop = FALSE],
      t[cens, , drop = FALSE], invsd[cens, , drop = FALSE], width[cens])
  }
  if (any(interval)) {
    sat[interval, ] <- probability(ly[interval], upper[interval],
      C((ly[interval] + upper[interval])/2), t[interval, , drop = FALSE],
      invsd[interval, , drop = FALSE], span[interval])
  }
  list(ell = ell, sat = sat, powers = j$powers)
}

#' Censored lognormal family with an arithmetic mean predictor
#'
#' An independent extended family for [mgcv::gam()] and [mgcv::bam()]. The
#' response and censoring limits are supplied on their original, positive scale.
#' With `mu = exp(eta)`, the model is
#' `log(Y) ~ N(log(mu) - sigma^2/2, sigma^2)`, so `mu = E[Y]`.
#'
#' @param theta Log-response standard deviation: `NULL` estimates it starting
#'   at one; a positive value fixes it; a negative value specifies the negative
#'   of its starting value for estimation. Zero is not allowed.
#' @param link Mean link; currently only `"log"` is supported.
#' @details A positive numeric vector specifies uncensored observations. For
#'   censoring, supply a two-column numeric matrix. Equal columns denote an
#'   exact observation; `(lower, upper)` denotes a finite interval;
#'   `(limit, -Inf)` denotes left censoring; `(limit, Inf)` denotes right
#'   censoring. All finite limits must be positive, with lower < upper.
#'   Limits can vary by observation. Missing rows are handled by the fitting
#'   function's `na.action`. Censoring is assumed non-informative.
#'
#'   Prior weights are non-negative case (likelihood) weights, not log-scale
#'   precision weights as in `mgcv::cnorm()`. The single global `sigma` is
#'   estimated alongside smoothing parameters; it has no separate predictor.
#'   `family$getTheta(TRUE)` retrieves `sigma` (`FALSE` gives `log(sigma)`).
#'
#'   Fitted values and response predictions are the arithmetic mean of the
#'   underlying, uncensored distribution, conditional on the model terms.
#'   Standard prediction errors condition on the estimated `sigma`, as usual
#'   for an extended family. There is no posterior bias correction or automatic
#'   integration over random effects. Link predictions are natural log means.
#'   This implementation does not call or wrap `cnorm()` or `mgcvUtils::clognorm()`.
#'
#'   Deviance compares each observation's likelihood with its maximum over the
#'   location at fixed `sigma`. AIC and log likelihood include the original-scale
#'   density Jacobian for exact observations. Response and Pearson residuals
#'   are `NA` for censored rows because the response is not observed; deviance
#'   residuals use censoring bounds and the fitted log-response location.
#'   The `rd`, `qf`, and `cdf` components describe the uncensored distribution.
#'
#'   The distribution components accept
#'   `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)` and
#'   `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)`.
#'   Set `lower.tail = FALSE` for upper-tail probabilities or quantiles.
#'   `log.p = TRUE` supplies log probabilities to `qf`; `logp = TRUE` returns
#'   log probabilities from `cdf`. These options retain accuracy for extreme
#'   tails without subtracting a CDF from one or logging an underflowed
#'   probability. Defaults use lower tails and ordinary probabilities.
#'   The `logp` spelling follows the mgcv CDF convention; `log.p` follows
#'   R's quantile functions.
#'   Case weights `wt` and the external `scale` are unused: the distribution
#'   uses the response mean `mu` and the family's current `sigma`.
#' @return An object of class `c("extended.family", "family")`.
#' @export
#' @examples
#' set.seed(12)
#' dat <- data.frame(x = runif(200))
#' mu <- exp(0.5 + sin(2*pi*dat$x))
#' y <- rlnorm(200, log(mu) - 0.6^2/2, 0.6)
#' dat$y <- cbind(pmax(y, 0.5), ifelse(y < 0.5, -Inf, y))
#' fit <- mgcv::gam(y ~ s(x), data = dat, family = clognormal(), method = "REML")
#' head(predict(fit, type = "response"))
#' fit$family$getTheta(TRUE)
clognormal <- function(theta = NULL, link = "log") {
  if (!identical(link, "log")) stop('clognormal supports only link = "log"')
  if (!is.null(theta) && (!is.numeric(theta) || length(theta) != 1L ||
      !is.finite(theta) || theta == 0)) stop("theta must be NULL or a finite non-zero number")
  estimated <- is.null(theta) || theta < 0
  state <- new.env(parent = emptyenv())
  state$theta <- if (is.null(theta)) 0 else log(abs(theta))
  getTheta <- function(trans = FALSE) if (trans) exp(state$theta) else state$theta
  putTheta <- function(theta) { state$theta <- theta }
  dev.resids <- function(y, mu, wt, theta = NULL) {
    if (is.null(theta)) theta <- getTheta()
    r <- .clnormal_eval(y, mu, theta)
    d <- -2*wt*(r$ell[, 1] - r$sat[, 1])
    d[wt == 0] <- 0
    d
  }
  Dd <- function(y, mu, theta, wt, level = 0) {
    r <- .clnormal_eval(y, mu, theta, if (level == 0) 2L else if (level == 1) 3L else 4L)
    d <- -2*wt*(r$ell - r$sat)
    d[wt == 0, ] <- 0
    take <- function(a, b) d[, which(r$powers[, 1] == a & r$powers[, 2] == b)] * factorial(a)*factorial(b)
    out <- list(Dmu = take(1, 0), Dmu2 = take(2, 0),
      EDmu2 = 2*wt/(mu^2*exp(2*theta)))
    if (level > 0) out <- c(out, list(Dth = take(0, 1), Dmuth = take(1, 1),
      Dmu3 = take(3, 0), Dmu2th = take(2, 1)))
    if (level > 1) out <- c(out, list(Dmu4 = take(4, 0), Dth2 = take(0, 2),
      Dmuth2 = take(1, 2), Dmu2th2 = take(2, 2), Dmu3th = take(3, 1)))
    out
  }
  ls <- function(y, w, theta, scale) {
    r <- .clnormal_eval(y, rep(1, length(if (is.matrix(y)) y[, 1] else y)), theta, 2L)
    take <- function(b) {
      v <- w*r$sat[, which(r$powers[, 1] == 0 & r$powers[, 2] == b)]*factorial(b)
      v[w == 0] <- 0
      v
    }
    list(ls = sum(take(0)), lsth1 = sum(take(1)),
         LSTH1 = matrix(take(1), ncol = 1), lsth2 = matrix(sum(take(2)), 1, 1))
  }
  aic <- function(y, mu, theta = NULL, wt, dev) {
    if (is.null(theta)) theta <- getTheta()
    v <- -2*wt*.clnormal_eval(y, mu, theta)$ell[, 1]
    v[wt == 0] <- 0
    sum(v)
  }
  initialize <- expression({
    y <- family$prepare.response(y)
    if (any(!is.finite(weights) | weights < 0) || !any(weights > 0))
      stop("clognormal requires finite non-negative weights with positive total")
    mustart <- as.numeric(y)
  })
  subsety <- function(y, ind) {
    if (is.matrix(y)) return(y[ind, , drop = FALSE])
    b <- attr(y, "censor")
    ans <- y[ind]
    if (!is.null(b)) attr(ans, "censor") <- b[ind]
    ans
  }
  postproc <- function(family, y, prior.weights, fitted, linear.predictors, offset, intercept) {
    null <- utils::getFromNamespace("find.null.dev", "mgcv")
    list(null.deviance = null(family, .clnormal_response(y), linear.predictors,
      offset, prior.weights), family = paste0("clognormal(", round(family$getTheta(TRUE), 3), ")"))
  }
  residuals <- function(object, type = "deviance", ...) {
    type <- match.arg(type, c("deviance", "pearson", "scaled.pearson", "working", "response"))
    if (type == "working") return(object$residuals)
    y <- .clnormal_response(object$y); b <- attr(y, "censor")
    mu <- object$fitted.values; w <- object$prior.weights
    if (type == "deviance") {
      loc <- log(mu) - object$family$getTheta(TRUE)^2/2
      target <- log(y)
      ii <- is.finite(b) & b != y
      target[ii] <- (log(y[ii]) + log(b[ii]))/2
      sgn <- sign(target - loc)
      sgn[b == -Inf] <- -1; sgn[b == Inf] <- 1
      return(sgn*sqrt(pmax(0, object$family$dev.resids(y, mu, w))))
    }
    ans <- as.numeric(y) - mu
    if (type != "response") ans <- ans*sqrt(w/object$family$variance(mu))
    ans[b != y] <- NA_real_
    ans
  }
  lk <- stats::make.link(link)
  structure(c(list(family = "clognormal", link = link), lk[c("linkfun", "linkinv", "mu.eta", "valideta")],
    list(validmu = function(mu) all(is.finite(mu) & mu > 0),
      dev.resids = dev.resids, Dd = Dd, ls = ls, aic = aic,
      initialize = initialize, prepare.response = .clnormal_response,
      subsety = subsety, postproc = postproc, residuals = residuals,
      n.theta = as.integer(estimated), ini.theta = state$theta,
      getTheta = getTheta, putTheta = putTheta, no.r.sq = TRUE,
      variance = function(mu) mu^2*expm1(getTheta(TRUE)^2),
      rd = function(mu, wt, scale) stats::rlnorm(length(mu), log(mu) - getTheta(TRUE)^2/2, getTheta(TRUE)),
      qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
        stats::qlnorm(p, log(mu) - getTheta(TRUE)^2/2, getTheta(TRUE),
                      lower.tail = lower.tail, log.p = log.p),
      cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
        stats::plnorm(q, log(mu) - getTheta(TRUE)^2/2, getTheta(TRUE),
                      lower.tail = lower.tail, log.p = logp))),
    class = c("extended.family", "family"))
}
