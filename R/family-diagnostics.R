#' Response distributions for fitted-model diagnostics
#'
#' Every inlaws family supplies a `gratia_diagnostics(mu, wt, scale, y, eta)`
#' function. It returns response moments and functions bound to the supplied
#' fitted parameters. No gratia installation is needed to construct or fit a
#' family. The native `rd`, `cdf` and `qf` functions keep their existing arguments.
#'
#' @details Count families supply CDF left limits, including at zero and between
#'   integers. Mixture and ordered-beta moments describe the complete response.
#'   Cumulative-link probability functions use the latent predictor `eta` and
#'   its moments describe numeric category values. GPD, generalized gamma and
#'   log-logistic families report whether moments exist; probability diagnostics
#'   remain available when a required moment is undefined.
#'
#'   Beta-binomial diagnostics describe success counts and require the fitted
#'   family's trial totals or a two-column success/failure response in `y`.
#'   Totals are recorded after fitting subsets and missing rows are removed.
#'   Its CDF and quantiles enumerate the relevant count support; very large trial
#'   totals can make these calculations expensive.
#'
#'   Dirichlet and Dirichlet-multinomial diagnostics provide a matrix of marginal
#'   response moments and CDFs, and simulate response columns jointly. They do
#'   not supply a joint inverse CDF. Componentwise Pearson residuals use marginal
#'   variances. Dirichlet deviance residuals are undefined; the componentwise
#'   Dirichlet-multinomial interface also requires PIT, quantile, response or
#'   Pearson residuals. Its existing native row-level deviance remains available
#'   through `residuals()`.
#'
#'   Censored families decode the observed intervals for probability diagnostics.
#'   Traditional simulated QQ and worm references additionally require an
#'   explicit `observe` rule, even when all realised observations are exact:
#'   exact values do not identify unobserved censoring thresholds. Use
#'   `censoring_process()` for common fixed limits, including `lower = 0` and
#'   `upper = Inf` for an uncensored positive response. Custom observation
#'   functions must apply a row-independent rule; stored vectors of row-specific
#'   thresholds are not aligned by this interface. Public response samplers and
#'   posterior draws continue to describe the underlying uncensored response.
#' @seealso [censoring_process()]
#' @name family_diagnostics
NULL

# These factories capture the constructor's distribution functions, whose
# parameter environments are shared with getTheta()/putTheta(). Rebuilding or
# serializing the family therefore keeps diagnostics and fitting in agreement.
.inlaws_diagnostic_family <- function(family, kind, observe = NULL) {
  force(family); force(kind); force(observe)
  family$gratia_diagnostics <- NULL
  if (!is.null(observe) && !is.function(observe))
    stop("observe must be a function or NULL", call. = FALSE)
  factory <- function(mu, wt, scale, y, eta) {
    out <- list()
    count <- kind %in% c("nbls", "cmp", "cmpls", "gpoisson", "pig", "pigls", "mixture")
    censored <- kind %in% c("cgamma", "cgammals", "clognormal", "clognormalls")
    m <- if (is.matrix(mu)) mu[, 1L] else mu
    if (count || censored || kind %in% c("log_logistic", "log_logisticls")) out$mean <- m
    if (count) {
      # ceil(q)-1 also gives the correct left limit between support values.
      out$cdf_left <- function(q, log_p = FALSE, lower_tail = TRUE) {
        lp <- family$cdf(ceiling(q) - 1, mu, wt, scale, logp = TRUE)
        .inlaws_log_probability(lp, log_p, lower_tail)
      }
      out$variance <- switch(kind,
        nbls = m + m^2 / mu[, 2L],
        cmp = family$variance(m),
        cmpls = vapply(.cmp_moments(m, mu[, 2L], family$control), `[[`, numeric(1), "variance"),
        gpoisson = family$variance(m), pig = family$variance(m),
        pigls = m + mu[, 2L] * m^2,
        mixture = {
          positive <- if (family$hurdle) -expm1(stats::dnbinom(0, mu = m, size = mu[, 3L], log = TRUE)) else 1
          out$mean <- mu[, 2L] * m / positive
          mu[, 2L] * (m + m^2 / mu[, 3L] + m^2) / positive - out$mean^2
        })
    }
    if (kind == "cumulative_link") {
      pars <- if (is.null(eta)) mu else eta
      out$cdf <- function(q, log_p = FALSE, lower_tail = TRUE)
        family$cdf(q, pars, wt, scale, log.p = log_p, lower.tail = lower_tail)
      out$cdf_left <- function(q, log_p = FALSE, lower_tail = TRUE)
        out$cdf(ceiling(q) - 1, log_p, lower_tail)
      out$qf <- function(p, log_p = FALSE, lower_tail = TRUE)
        family$qf(p, pars, wt, scale, log.p = log_p, lower.tail = lower_tail)
      out$rd <- function() as.numeric(family$rd(pars, wt, scale))
      R <- length(family$getTheta(TRUE)) + 1L
      survival <- matrix(vapply(seq_len(R - 1L), function(k)
        out$cdf(k, lower_tail = FALSE), numeric(length(pars))), nrow = length(pars))
      out$mean <- 1 + rowSums(survival)
      out$variance <- 1 + drop(survival %*% (2 * seq_len(R - 1L) + 1)) - out$mean^2
      out$diagnostic_residuals <- function(y, type) {
        if (type == "response") return(y - out$mean)
        if (type == "pearson") return(sqrt(rep_len(wt, length(y))) * (y - out$mean) / sqrt(out$variance))
        d <- family$dev.resids(y, pars, rep_len(wt, length(y)))
        attr(d, "sign") * sqrt(pmax(0, d))
      }
    }
    if (kind == "ordbeta") {
      th <- family$getTheta(TRUE)
      linear <- if (is.null(eta)) stats::qlogis(mu) else eta
      p <- .ordbeta_components(linear, th)
      out$mean <- p[, 3L] + p[, 2L] * mu
      out$variance <- p[, 3L] + p[, 2L] * (mu^2 + mu * (1 - mu) / (1 + th[1L])) - out$mean^2
      out$cdf_left <- function(q, log_p = FALSE, lower_tail = TRUE) {
        q <- rep_len(q, length(mu))
        lp <- family$cdf(q, mu, wt, scale, logp = TRUE)
        lp[q == 0] <- -Inf
        lp[q == 1] <- stats::plogis(th[3L] - linear[q == 1], log.p = TRUE)
        .inlaws_log_probability(lp, log_p, lower_tail)
      }
    }
    if (censored) {
      out$bounds <- .inlaws_censor_bounds
      out$reference_available <- is.function(observe)
      out$reference_reason <- "Supply the censoring observation process to generate reference residuals."
      if (is.function(observe)) out$observe <- observe
      nuisance <- if (is.matrix(mu)) mu[, 2L] else family$getTheta(TRUE)
      out$variance <- if (kind %in% c("cgamma", "cgammals")) m^2 * nuisance else m^2 * expm1(nuisance^2)
    }
    if (kind == "gen_gammalss") {
      out$continuous <- TRUE
      out$mean <- m
      out$variance <- .gen_gamma_variance(mu)
      out$variance_available <- all(mu[, 3L] > .5)
      out$moment_reason <- "Generalized gamma variance requires r > 0.5 (1 + 2*sigma*Q > 0)."
    }
    if (kind == "gpd") {
      shape <- mu[, 2L]
      out$continuous <- TRUE
      out$mean <- ifelse(shape < 1, m / (1 - shape), NA_real_)
      out$variance <- ifelse(shape < .5, m^2 / ((1 - shape)^2 * (1 - 2 * shape)), Inf)
      out$mean_available <- all(shape < 1)
      out$variance_available <- all(shape < .5)
      out$moment_reason <- "GPD means require shape < 1 and variances require shape < 0.5."
    }
    if (kind %in% c("log_logistic", "log_logisticls")) {
      s <- if (is.matrix(mu)) mu[, 2L] else family$getTheta(TRUE)
      out$continuous <- TRUE
      out$variance <- .llog_variance(m, s)
      out$variance_available <- all(s < .5)
      out$moment_reason <- "Log-logistic variance requires log-response scale < 0.5."
    }
    if (kind == "betabinomial") {
      trials <- if (is.matrix(y) && NCOL(y) == 2L) rowSums(y) else family$trials
      if (is.null(trials) || length(trials) != NROW(mu))
        stop("Beta-binomial diagnostics require trial totals in fitting order", call. = FALSE)
      out <- .inlaws_bb_distribution(mu[, 1L], mu[, 2L], trials)
    }
    if (kind %in% c("dirichlet", "dirmult")) {
      pars <- if (is.null(eta)) mu else eta
      pr <- .dirichlet_parameters(pars)
      a <- pr$mu * pr$phi
      b <- (1 - pr$mu) * pr$phi
      out$independent <- FALSE
      out$joint_available <- TRUE
      out$dependence_reason <- "Response columns share a fixed row sum."
      out$qf <- NULL # A multivariate distribution has no joint inverse CDF.
      out$mean <- pr$mu
      out$variance <- pr$mu * (1 - pr$mu) / (1 + pr$phi)
      if (kind == "dirichlet") {
        out$continuous <- TRUE
        out$cdf <- function(q, log_p = FALSE, lower_tail = TRUE)
          stats::pbeta(q, a, b, log.p = log_p, lower.tail = lower_tail)
        out$rd <- function() family$rd(pars, wt, scale)
      } else {
        if (!is.matrix(y) || NROW(y) != NROW(mu) || NCOL(y) != NCOL(pr$mu))
          stop("Dirichlet-multinomial diagnostics require the response count matrix", call. = FALSE)
        total <- rowSums(y)
        out$mean <- pr$mu * total
        out$variance <- pr$mu * (1 - pr$mu) * total * (total + pr$phi) / (1 + pr$phi)
        margins <- .inlaws_bb_distribution(as.vector(pr$mu), rep(pr$phi, NCOL(pr$mu)), rep(total, NCOL(pr$mu)))
        out$cdf <- function(q, log_p = FALSE, lower_tail = TRUE)
          matrix(margins$cdf(q, log_p, lower_tail), NROW(mu))
        out$cdf_left <- function(q, log_p = FALSE, lower_tail = TRUE)
          matrix(margins$cdf_left(q, log_p, lower_tail), NROW(mu))
        out$rd <- function() {
          probabilities <- .inlaws_dirichlet_draw(a)
          t(vapply(seq_along(total), function(i)
            as.numeric(stats::rmultinom(1L, total[i], probabilities[i, ])), numeric(NCOL(pr$mu))))
        }
      }
      out$diagnostic_residuals <- function(y, type) {
        if (type == "response") return(y - out$mean)
        if (type == "pearson") return((y - out$mean) / sqrt(out$variance))
        stop(if (kind == "dirichlet") "Dirichlet deviance residuals are undefined" else
          "Componentwise Dirichlet-multinomial deviance residuals are unavailable; use PIT, quantile, response or Pearson residuals", call. = FALSE)
      }
    }
    out
  }
  family$gratia_diagnostics <- factory
  family
}

.inlaws_log_probability <- function(lp, log_p, lower_tail) {
  if (!lower_tail) lp <- log(-expm1(lp))
  if (log_p) lp else exp(lp)
}

.inlaws_censor_bounds <- function(y) {
  if (is.matrix(y)) { value <- y[, 1L]; limit <- y[, 2L] } else {
    value <- as.numeric(y); limit <- attr(y, "censor")
    if (is.null(limit)) limit <- value
  }
  list(lower = pmin(value, limit), upper = pmax(value, limit))
}

# Sum the beta-binomial PMF directly, evaluating the requested tail rather than
# subtracting a lower-tail probability from one. Totals are observation-specific.
.inlaws_bb_distribution <- function(p, phi, trials) {
  force(p); force(phi); force(trials)
  n <- length(p)
  log_mass <- function(k, i) {
    if (p[i] == 0 || p[i] == 1) return(ifelse(k == trials[i] * p[i], 0, -Inf))
    lchoose(trials[i], k) +
      lbeta(k + p[i] * phi[i], trials[i] - k + (1 - p[i]) * phi[i]) -
      lbeta(p[i] * phi[i], (1 - p[i]) * phi[i])
  }
  phi <- rep_len(phi, n)
  if (length(trials) != n || anyNA(p) || anyNA(phi) || anyNA(trials) ||
      any(!is.finite(p) | p < 0 | p > 1) || any(!is.finite(phi) | phi <= 0) ||
      any(!is.finite(trials) | trials < 0 | trials != floor(trials)))
    stop("Beta-binomial diagnostics require valid probabilities, positive precision and integer trial totals", call. = FALSE)
  cdf <- function(q, log_p = FALSE, lower_tail = TRUE) {
    q <- rep_len(q, n)
    lp <- vapply(seq_len(n), function(i) {
      at <- floor(q[i]); size <- trials[i]
      if (is.na(at)) return(NA_real_)
      if (at < 0) return(if (lower_tail) -Inf else 0)
      if (at >= size) return(if (lower_tail) 0 else -Inf)
      k <- if (lower_tail) seq.int(0, at) else seq.int(at + 1, size)
      z <- log_mass(k, i); top <- max(z)
      if (top == -Inf) -Inf else min(0, top + log(sum(exp(z - top))))
    }, numeric(1))
    if (log_p) lp else exp(lp)
  }
  qf <- function(p, log_p = FALSE, lower_tail = TRUE) {
    p <- rep_len(p, n)
    if (anyNA(p) || any(if (log_p) p > 0 else p < 0 | p > 1))
      stop("Beta-binomial probabilities must lie in [0, 1]", call. = FALSE)
    vapply(seq_len(n), function(i) {
      prob <- if (log_p) exp(p[i]) else p[i]
      if (!lower_tail) prob <- if (log_p) -expm1(p[i]) else 1 - prob
      if (prob == 0) return(0)
      if (prob == 1) return(trials[i])
      k <- seq.int(0, trials[i]); cs <- cumsum(exp(log_mass(k, i)))
      cs[length(cs)] <- 1
      k[which(cs >= prob)[1L]]
    }, numeric(1))
  }
  list(mean = trials * p, variance = trials * p * (1 - p) * (trials + phi) / (1 + phi),
    cdf = cdf, cdf_left = function(q, log_p = FALSE, lower_tail = TRUE)
      cdf(ceiling(q) - 1, log_p, lower_tail), qf = qf,
    rd = function() stats::rbinom(n, trials, stats::rbeta(n, p * phi, (1 - p) * phi)))
}

.inlaws_dirichlet_draw <- function(a) {
  z <- matrix(log(stats::rgamma(length(a), a + 1)) + log(stats::runif(length(a))) / a, NROW(a))
  z <- exp(z - apply(z, 1L, max))
  z / rowSums(z)
}

#' A fixed censoring rule for diagnostic simulations
#'
#' @param lower,upper Scalar limits shared by all observations. Values below
#'   `lower` are left censored; values above `upper` are right censored.
#' @return A function encoding a simulated vector with a `censor` attribute.
#' @details Pass this function as `observe` to a censored family constructor.
#'   The same scalar limits are applied after subsetting and missing-row removal.
#'   Row-specific limits are not supported by this helper. Without a known rule,
#'   censored families retain interval PIT/quantile diagnostics but cannot generate
#'   traditional distribution-based QQ or worm reference residuals. The rule is
#'   used only for diagnostic references, not public or posterior response draws.
#' @export
#' @examples
#' f <- cgamma(theta = 0.6, observe = censoring_process(lower = 0.4))
#' f$gratia_diagnostics(rep(2, 3), 1, 1, c(1, 2, 3), NULL)$observe(c(.1, 1, 2))
censoring_process <- function(lower = 0, upper = Inf) {
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 1L ||
      length(upper) != 1L || is.na(lower) || is.na(upper) || lower < 0 || lower >= upper)
    stop("censoring limits must be ordered nonnegative scalars", call. = FALSE)
  force(lower); force(upper)
  function(y) {
    value <- pmin(upper, pmax(lower, y))
    structure(value, censor = ifelse(y < lower, -Inf, ifelse(y > upper, Inf, y)))
  }
}
