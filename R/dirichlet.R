#' Dirichlet mean and concentration additive models
#'
#' A general family for strictly positive compositional responses, fitted by
#' [mgcv::gam()] with REML smoothing parameter estimation.
#'
#' @param K Number of non-reference components: the response has `K + 1`
#'   columns and the model has `K + 1` linear predictors.
#' @return An object inheriting from `general.family`, `extended.family`, and
#'   `family`.
#' @details
#' Supply a list of `K + 1` formulas. The first formula specifies a matrix
#' response. The first `K` predictors model log mean ratios relative to the
#' first response column; the final predictor models log concentration, phi.
#' With `mu = softmax(c(0, eta[1:K]))`, the Dirichlet shapes are `phi * mu`.
#' The concentration formula can be intercept-only or contain smooths.
#'
#' Responses must be finite, strictly between zero and one, and sum to one
#' within `sqrt(.Machine$double.eps)`. Only rounding-level discrepancies are
#' normalized. Missing observations follow `gam()`'s `na.action`. Nonnegative
#' weights multiply each observation's log likelihood; they are not trial
#' counts or inverse variances. At least one weight must be positive.
#'
#' Response predictions are matrices of mean proportions, in response-column
#' order. mgcv currently drops prediction column names; component names are
#' retained in `colnames(fit$y)` and in simulated responses. Standard errors
#' use the full coefficient covariance. Link
#' predictions contain `K` log ratios followed by log concentration; exponentiate
#' the final column to obtain concentration. Internally, fitted values retain
#' these link-scale predictors for simulation. Use response predictions to
#' obtain fitted mean proportions.
#'
#' Response residuals are componentwise differences. Pearson residuals are
#' nonnegative distances, one per observation, whose square is
#' `weight * (phi + 1) * sum((y - mu)^2 / mu)`. They account for the simplex
#' covariance, but are not normally distributed. Use simulation-based checks
#' via `mgcv::qq.gam(fit, type = "pearson", rep = 100)`. `gam.check()` and
#' `k.check()` request deviance residuals internally and are not supported.
#' Deviance residuals are not
#' defined here: an unrestricted saturated Dirichlet likelihood is unbounded.
#' Deviance and null deviance are `NA`, so no deviance explained is reported.
#'
#' This family supports REML in `gam()`, not `bam()`, `gamm()`, or NCV.
#' Fourth-order derivative storage grows as `choose(K + 4, 4)` per observation
#' for the highest order alone, so modest component counts are recommended.
#' Changing the reference by reordering columns can change penalized fits.
#'
#' The mean/concentration parameterization agrees with `brms::dirichlet()`
#' using its default logit mean link, log concentration link, and first-column
#' reference. The fitting criteria and uncertainty calculations differ:
#' this family uses penalized likelihood and REML, not posterior sampling.
#' Use explicit namespaces when both packages are attached.
#' The family's `rd(mu, wt, scale)` callback takes link-scale predictors and
#' returns a response matrix. The compatibility arguments `wt` and `scale`
#' do not change the distribution. Very
#' small shape parameters can yield simulated zeros through floating-point
#' underflow; such draws cannot be used directly as fitting responses.
#'
#' @examples
#' set.seed(12)
#' n <- 150
#' dat <- data.frame(x = runif(n))
#' eta <- cbind(0, sin(2 * pi * dat$x), -dat$x)
#' mu <- exp(eta) / rowSums(exp(eta))
#' y <- matrix(rgamma(n * 3, shape = 12 * mu), n, 3)
#' y <- y / rowSums(y)
#' dat$y <- y
#' b <- mgcv::gam(list(y ~ s(x, k = 5), ~ s(x, k = 5), ~ 1),
#'                family = dirichlet(K = 2), data = dat, method = "REML")
#' head(predict(b, type = "response"))
#' head(exp(predict(b, type = "link")[, 3]))
#' @export
dirichlet <- function(K = 1) {
  if (length(K) != 1L || !is.numeric(K) || !is.finite(K) ||
      K < 1 || K != floor(K)) stop("K must be a positive integer")
  q <- K + 1L
  jets <- .dirichlet_jets(q)
  tri <- mgcv::trind.generator(q)

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, sandwich = FALSE) {
    if (is.null(eta)) eta <- .dirichlet_eta(X, coef, offset)
    order <- if (deriv == 0) 0L else if (deriv == 1) 2L else
      if (deriv < 4) 3L else 4L
    z <- .dirichlet_derivatives(y, eta, wt, jets, order)
    if (!deriv) return(list(l = z$l))
    if (!is.finite(z$l)) stop("non-finite Dirichlet derivatives at current coefficients")
    out <- mgcv::gamlss.gH(X, attr(X, "lpi"), z$l1, z$l2, tri$i2,
                          l3 = z$l3, i3 = tri$i3, l4 = z$l4, i4 = tri$i4,
                          d1b = d1b, d2b = d2b, deriv = deriv - 1,
                          fh = fh, D = D, sandwich = sandwich)
    out$l <- z$l
    out
  }

  preinitialize <- function(G) {
    y <- G$y
    if (!is.matrix(y) || !is.numeric(y) || ncol(y) != q)
      stop("Dirichlet response must be a numeric matrix with K + 1 columns")
    if (any(!is.finite(y)) || any(y <= 0 | y >= 1))
      stop("Dirichlet response entries must be finite and strictly between zero and one")
    if (any(abs(rowSums(y) - 1) > sqrt(.Machine$double.eps)))
      stop("Dirichlet response rows must sum to one")
    if (any(!is.finite(G$w)) || any(G$w < 0) || !any(G$w > 0))
      stop("Dirichlet weights must be finite, nonnegative, and include a positive weight")
    fam <- G$family
    fam$data <- list(component.names = colnames(y))
    list(y = y / rowSums(y), family = fam)
  }

  initialize <- expression({
    nobs <- nrow(y)
    n <- rep(1, nobs)
    if (is.null(start)) start <- family$initial(y, x, weights, offset, E)
  })

  predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                      beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta)) eta <- .dirichlet_eta(X, beta, off) else se <- FALSE
    mu <- .dirichlet_parameters(eta)$mu
    colnames(mu) <- family$data$component.names
    if (!se) return(list(fit = mu))
    jj <- attr(X, "lpi")
    ans <- mu * 0
    for (j in seq_len(ncol(mu))) {
      gradient <- matrix(0, nrow(X), ncol(X))
      for (k in seq_len(ncol(mu) - 1L)) {
        d <- mu[, j] * ((j == k + 1L) - mu[, k + 1L])
        gradient[, jj[[k]]] <- gradient[, jj[[k]], drop = FALSE] +
          d * X[, jj[[k]], drop = FALSE]
      }
      ans[, j] <- sqrt(pmax(0, rowSums((gradient %*% Vb) * gradient)))
    }
    list(fit = mu, se.fit = ans)
  }

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    if (type == "deviance")
      stop("Dirichlet deviance residuals are undefined; use type = 'pearson' or 'response'")
    p <- .dirichlet_parameters(object$linear.predictors)
    r <- object$y - p$mu
    if (type == "response") return(r)
    sqrt(object$prior.weights * (p$phi + 1) * rowSums(r^2 / p$mu))
  }

  rd <- function(mu, wt, scale) {
    p <- .dirichlet_parameters(mu)
    a <- p$mu * p$phi
    # Gamma(a + 1) * U^(1/a), evaluated on the log scale, also works
    # when Gamma(a) would underflow for small shape parameters.
    z <- matrix(log(stats::rgamma(length(a), shape = a + 1)) +
                  log(stats::runif(length(a))) / a, nrow(a), ncol(a))
    z <- exp(z - apply(z, 1L, max))
    z / rowSums(z)
  }

  structure(list(family = "dirichlet", link = NULL, nlp = q, ll = ll,
                 linfo = replicate(q, stats::make.link("identity"), simplify = FALSE),
                 tri = tri, initialize = initialize, initial = .dirichlet_initial,
                 preinitialize = preinitialize,
                 postproc = expression({
                   object$deviance <- object$null.deviance <- NA_real_
                 }),
                 predict = predict, residuals = residuals, rd = rd,
                 sandwich = function(y, X, coef, wt, family, offset = NULL)
                   ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 2, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}

.dirichlet_eta <- function(X, beta, offset) {
  jj <- attr(X, "lpi")
  eta <- matrix(0, nrow(X), length(jj))
  for (j in seq_along(jj)) {
    eta[, j] <- drop(X[, jj[[j]], drop = FALSE] %*% beta[jj[[j]]])
    if (length(offset) >= j && !is.null(offset[[j]]))
      eta[, j] <- eta[, j] + offset[[j]]
  }
  eta
}

.dirichlet_parameters <- function(eta) {
  q <- ncol(eta)
  z <- cbind(0, eta[, seq_len(q - 1L), drop = FALSE])
  z <- z - apply(z, 1L, max)
  logmu <- z - log(rowSums(exp(z)))
  list(mu = exp(logmu), logmu = logmu, phi = exp(eta[, q]))
}

# Taylor coefficients use multi-index factorial scaling. Ordered tuples match
# mgcv's symmetric derivative packing exactly. Cache convolution indices once.
.dirichlet_jets <- function(q) {
  tuples <- list(integer())
  extend <- function(prefix, lo, left) {
    if (!left) { tuples[[length(tuples) + 1L]] <<- prefix; return(invisible(NULL)) }
    for (j in seq.int(lo, q)) extend(c(prefix, j), j, left - 1L)
  }
  for (r in 1:4) extend(integer(), 1L, r)
  powers <- t(vapply(tuples, tabulate, integer(q), nbins = q))
  degree <- rowSums(powers)
  key <- function(x) paste(x, collapse = ",")
  keys <- apply(powers, 1L, key)
  pairs <- lapply(seq_along(tuples), function(i) {
    left <- which(apply(powers, 1L, function(a) all(a <= powers[i, ])))
    right <- vapply(left, function(j) match(key(powers[i, ] - powers[j, ]), keys), integer(1))
    list(left = left, right = right)
  })
  list(powers = powers, degree = degree, pairs = pairs,
       factorial = apply(powers, 1L, function(a) prod(factorial(a))))
}

.dirichlet_multiply <- function(a, b, jets) {
  ans <- a * 0
  for (j in seq_len(ncol(a))) {
    ix <- jets$pairs[[j]]
    ans[, j] <- rowSums(a[, ix$left, drop = FALSE] * b[, ix$right, drop = FALSE])
  }
  ans
}

.dirichlet_compose <- function(a, derivatives, jets, order) {
  delta <- a
  delta[, 1L] <- 0
  ans <- a * 0
  ans[, 1L] <- derivatives[[1L]]
  power <- delta
  for (r in seq_len(order)) {
    ans <- ans + power * (derivatives[[r + 1L]] / factorial(r))
    if (r < order) power <- .dirichlet_multiply(power, delta, jets)
  }
  ans
}

.dirichlet_derivatives <- function(y, eta, wt, jets, order = 4L) {
  n <- nrow(y); q <- ncol(y)
  if (is.null(wt)) wt <- rep(1, n)
  p <- .dirichlet_parameters(eta)
  alpha <- p$mu * p$phi
  active <- wt > 0
  if (any(!is.finite(alpha[active, , drop = FALSE])) ||
      any(alpha[active, , drop = FALSE] <= 0)) return(list(l = -Inf))
  # Zero-weight observations must not create 0 * Inf in the likelihood.
  alpha[!active, ] <- 1
  phi <- rowSums(alpha)
  logy <- log(y)
  logy[!active, ] <- 0
  l0 <- wt * (lgamma(phi) - rowSums(lgamma(alpha)) + rowSums((alpha - 1) * logy))
  if (!order) return(list(l = sum(l0)))
  cols <- which(jets$degree <= order)
  powers <- jets$powers[cols, , drop = FALSE]
  m <- length(cols)
  # Scale exponentials by their row maxima, held constant in differentiation.
  z <- cbind(0, eta[, seq_len(q - 1L), drop = FALSE])
  z[!active, ] <- 0
  z <- exp(z - apply(z, 1L, max))
  ej <- lapply(seq_len(q), function(j) {
    a <- matrix(0, n, m); a[, 1L] <- z[, j]
    if (j > 1L) {
      ix <- which(powers[, j - 1L] == jets$degree[cols])
      a[, ix] <- outer(z[, j], 1 / jets$factorial[ix])
    }
    a
  })
  den <- Reduce(`+`, ej)
  inv <- .dirichlet_compose(den,
    lapply(0:order, function(r) (-1)^r * factorial(r) / den[, 1L]^(r + 1L)), jets, order)
  phij <- matrix(0, n, m)
  ix <- which(powers[, q] == jets$degree[cols])
  phij[, ix] <- outer(phi, 1 / jets$factorial[ix])
  lg <- function(a) .dirichlet_compose(a,
    c(list(lgamma(a[, 1L])), lapply(seq_len(order), function(r)
      psigamma(a[, 1L], deriv = r - 1L))), jets, order)
  ans <- lg(phij)
  for (j in seq_len(q)) {
    a <- .dirichlet_multiply(.dirichlet_multiply(ej[[j]], inv, jets), phij, jets)
    ans <- ans - lg(a) + a * logy[, j]
  }
  ans <- ans * wt
  ans[!active, ] <- 0
  if (any(!is.finite(ans))) return(list(l = -Inf))
  out <- list(l = sum(l0), l3 = 0, l4 = 0)
  for (r in seq_len(order)) {
    ix <- which(jets$degree[cols] == r)
    out[[paste0("l", r)]] <- sweep(ans[, ix, drop = FALSE], 2L, jets$factorial[ix], `*`)
  }
  out
}

.dirichlet_initial <- function(y, X, wt, offset, E) {
  n <- nrow(y); q <- ncol(y)
  if (is.null(wt)) wt <- rep(1, n)
  mu <- colSums(y * wt) / sum(wt)
  v <- colSums(sweep(y, 2L, mu)^2 * wt) / sum(wt)
  phi <- stats::median(mu * (1 - mu) / pmax(v, .Machine$double.eps) - 1)
  phi <- min(1e4, max(0.1, phi))
  par <- log(mu * phi)
  logy <- colSums(log(y) * wt) / sum(wt)
  objective <- function(b) {
    a <- exp(b)
    val <- -(lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * logy))
    if (is.finite(val)) val else .Machine$double.xmax / 100
  }
  gradient <- function(b) {
    a <- exp(b)
    -a * (digamma(sum(a)) - digamma(a) + logy)
  }
  fit <- stats::optim(par, objective, gradient, method = "L-BFGS-B",
                      lower = -20, upper = 20)
  if (fit$convergence == 0 && all(is.finite(fit$par))) par <- fit$par
  target <- c(par[-1L] - par[1L], log(sum(exp(par))))
  jj <- attr(X, "lpi")
  design <- matrix(0, n * q, ncol(X))
  rhs <- numeric(n * q)
  for (j in seq_len(q)) {
    rows <- (j - 1L) * n + seq_len(n)
    design[rows, jj[[j]]] <- X[, jj[[j]], drop = FALSE] * sqrt(wt)
    off <- if (length(offset) >= j && !is.null(offset[[j]])) offset[[j]] else 0
    rhs[rows] <- (target[j] - off) * sqrt(wt)
  }
  if (nrow(E)) {
    design <- rbind(design, E)
    rhs <- c(rhs, rep(0, nrow(E)))
  }
  start <- qr.coef(qr(design), rhs)
  start[!is.finite(start)] <- 0
  start
}
