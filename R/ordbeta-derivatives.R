# Analytic, truncated multivariate Taylor arithmetic. Coefficients are
# derivatives divided by multi-index factorials; no finite differences are
# used in fitting. Only the derivative orders consumed by gam.fit4 are kept.
.ordbeta_basis <- function(nt, level = 2L, first = FALSE) {
  powers <- as.matrix(expand.grid(c(list(0:4), rep(list(0:2), nt))))
  keep <- rowSums(powers[, -1, drop = FALSE]) <= 2 & rowSums(powers) <= 4
  if (first) keep <- rowSums(powers) <= 1
  else if (level == 0) keep <- powers[, 1] <= 2 & rowSums(powers[, -1, drop = FALSE]) == 0
  else if (level == 1) keep <- rowSums(powers) <= 3 & rowSums(powers[, -1, drop = FALSE]) <= 1
  powers <- powers[keep, , drop = FALSE]
  key <- function(x) paste(x, collapse = ",")
  keys <- apply(powers, 1, key)
  pairs <- lapply(seq_len(nrow(powers)), function(k) {
    a <- which(apply(sweep(powers, 2, powers[k, ], "<="), 1, all))
    b <- match(apply(matrix(powers[k, ], length(a), ncol(powers), byrow = TRUE) -
                       powers[a, , drop = FALSE], 1, key), keys)
    cbind(a, b)
  })
  list(powers = powers, keys = keys, pairs = pairs,
       factorial = apply(factorial(powers), 1, prod), degree = max(rowSums(powers)))
}

.ordbeta_jet <- function(value, basis, variable = NULL) {
  ans <- matrix(0, length(value), nrow(basis$powers))
  ans[, 1] <- value
  if (!is.null(variable)) {
    e <- integer(ncol(basis$powers)); e[variable] <- 1L
    j <- match(paste(e, collapse = ","), basis$keys)
    if (!is.na(j)) ans[, j] <- 1
  }
  ans
}

.ordbeta_mul <- function(a, b, basis) {
  ans <- a * 0
  for (k in seq_along(basis$pairs)) {
    ij <- basis$pairs[[k]]
    ans[, k] <- rowSums(a[, ij[, 1], drop = FALSE] * b[, ij[, 2], drop = FALSE])
  }
  ans
}

.ordbeta_unary <- function(x, derivatives, basis) {
  out <- x * 0
  out[, 1] <- derivatives[[1]]
  dx <- x; dx[, 1] <- 0
  power <- dx
  for (k in seq_len(basis$degree)) {
    if (k > 1L) power <- .ordbeta_mul(power, dx, basis)
    out <- out + power * (derivatives[[k + 1L]] / factorial(k))
  }
  out
}

.ordbeta_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
.ordbeta_log1mexp <- function(x) {
  ans <- x
  small <- x <= log(2)
  ans[small] <- log(-expm1(-x[small]))
  ans[!small] <- log1p(-exp(-x[!small]))
  ans
}

.ordbeta_calculus <- function(basis, n) {
  constant <- function(x) .ordbeta_jet(rep_len(x, n), basis)
  multiply <- function(a, b) .ordbeta_mul(a, b, basis)
  unary <- function(x, kind) {
    v <- x[, 1]
    d <- switch(kind,
      exp = rep(list(exp(v)), 5),
      log = list(log(v), 1/v, -1/v^2, 2/v^3, -6/v^4),
      lgamma = list(lgamma(v), digamma(v), trigamma(v), psigamma(v, 2), psigamma(v, 3)),
      trigamma = lapply(1:5, function(k) psigamma(v, k)),
      softplus = {
        p <- stats::plogis(v); q <- stats::plogis(-v); pq <- p*q
        list(.ordbeta_softplus(v), p, pq, pq*(q-p), pq*(1-6*pq))
      },
      logistic = {
        p <- stats::plogis(v); q <- stats::plogis(-v); pq <- p*q
        list(p, pq, pq*(q-p), pq*(1-6*pq), pq*(q-p)*(1-12*pq))
      },
      log1mexp = {
        h <- 1/expm1(v)
        list(.ordbeta_log1mexp(v), h, -h*(1+h),
             h*(1+h)*(1+2*h), -h*(1+h)*(1+6*h+6*h^2))
      })
    .ordbeta_unary(x, d, basis)
  }
  list(c = constant, mul = multiply, f = unary)
}

.ordbeta_natural <- function(theta, mode) {
  cuts <- switch(mode, full = c(theta[2], theta[2] + exp(theta[3])),
                 lower = c(theta[2], Inf), upper = c(-Inf, theta[2]),
                 beta = c(-Inf, Inf))
  c(phi = unname(exp(theta[1])), cut_lower = unname(cuts[1]), cut_upper = unname(cuts[2]))
}

.ordbeta_loglik <- function(y, mu, theta, mode) {
  par <- .ordbeta_natural(theta, mode)
  eta <- stats::qlogis(mu)
  ans <- numeric(length(y))
  z <- y == 0; o <- y == 1; i <- !z & !o
  ans[z] <- -.ordbeta_softplus(eta[z] - par[2])
  ans[o] <- -.ordbeta_softplus(par[3] - eta[o])
  lp <- -.ordbeta_softplus(par[2] - eta[i]) - .ordbeta_softplus(eta[i] - par[3])
  if (mode == "full") lp <- lp + .ordbeta_log1mexp(exp(theta[3]))
  ans[i] <- lp + stats::dbeta(y[i], mu[i]*par[1], (1-mu[i])*par[1], log = TRUE)
  ans
}

.ordbeta_taylor <- function(y, mu, theta, mode, basis, fisher = FALSE) {
  n <- length(mu); calc <- .ordbeta_calculus(basis, n)
  C <- calc$c; M <- calc$mul; F <- calc$f
  m <- .ordbeta_jet(mu, basis, 1)
  t <- lapply(seq_along(theta), function(j) .ordbeta_jet(rep(theta[j], n), basis, j+1L))
  phi <- F(t[[1]], "exp"); om <- C(1)-m
  lm <- F(m, "log"); lom <- F(om, "log")
  a <- M(m, phi); b <- M(om, phi)
  lo <- hi <- C(0); nlo <- nhi <- C(1)
  lp0 <- lp1 <- NULL
  # Differentiate the gates in mu coordinates. Differentiating logit(mu)
  # first creates cancelling poles near the endpoints, losing accuracy in
  # the third/fourth derivatives of otherwise well-behaved endpoint terms.
  if (mode %in% c("full", "lower")) {
    dl <- F(M(F(t[[2]], "logistic"), om)+M(F(-t[[2]], "logistic"), m), "log")
    lp0 <- -F(-t[[2]], "softplus")+lom-dl
    lnlo <- -F(t[[2]], "softplus")+lm-dl
    lo <- F(lp0, "exp"); nlo <- F(lnlo, "exp")
  }
  if (mode %in% c("full", "upper")) {
    upper <- if (mode == "full") t[[2]] + F(t[[3]], "exp") else t[[2]]
    dh <- F(M(F(upper, "logistic"), om)+M(F(-upper, "logistic"), m), "log")
    lp1 <- -F(upper, "softplus")+lm-dh
    lnhi <- -F(-upper, "softplus")+lom-dh
    hi <- F(lp1, "exp"); nhi <- F(lnhi, "exp")
  }
  lpI <- switch(mode, beta=C(0), lower=lnlo, upper=lnhi,
                full=lm+lom-dl-dh-F(-upper,"softplus")-F(t[[2]],"softplus")+
                  F(F(t[[3]], "exp"), "log1mexp"))
  if (fisher) {
    pI <- F(lpI, "exp")
    cat <- M(lo, M(nlo, nlo)) + M(hi, M(nhi, nhi)) + M(pI, M(lo-hi, lo-hi))
    inv <- F(-2*(F(m, "log")+F(om, "log")), "exp")
    return(2*(M(cat, inv) + M(M(pI, M(phi, phi)), F(a, "trigamma")+F(b, "trigamma"))))
  }
  yy <- y; yy[y == 0 | y == 1] <- .5
  out <- lpI + F(phi, "lgamma")-F(a, "lgamma")-F(b, "lgamma") +
    (a-C(1))*log(yy) + (b-C(1))*log1p(-yy)
  if (any(y == 0)) out[y == 0, ] <- lp0[y == 0, , drop = FALSE]
  if (any(y == 1)) out[y == 1, ] <- lp1[y == 1, , drop = FALSE]
  -2*out
}

.ordbeta_Dd <- function(y, mu, theta, wt, level, mode, bases, fbasis) {
  wt <- rep_len(wt, length(y)); keep <- wt > 0
  basis <- bases[[level+1L]]
  all <- matrix(0, length(y), nrow(basis$powers))
  if (any(keep)) all[keep, ] <- .ordbeta_taylor(y[keep], mu[keep], theta, mode, basis)*wt[keep]
  extract <- function(dm, dt = integer(length(theta))) {
    j <- match(paste(c(dm, dt), collapse = ","), basis$keys)
    all[, j]*basis$factorial[j]
  }
  nt <- length(theta); unit <- diag(nt)
  pairs <- which(upper.tri(matrix(0, nt, nt), diag = TRUE), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
  first <- function(dm) matrix(vapply(seq_len(nt), function(j) extract(dm, unit[j, ]), numeric(length(y))), length(y), nt)
  second <- function(dm) matrix(vapply(seq_len(nrow(pairs)), function(j) extract(dm, unit[pairs[j,1], ]+unit[pairs[j,2], ]), numeric(length(y))), length(y), nrow(pairs))
  expected <- matrix(0, length(y), nrow(fbasis$powers))
  if (any(keep)) expected[keep, ] <- .ordbeta_taylor(y[keep], mu[keep], theta, mode, fbasis, TRUE)*wt[keep]
  r <- list(Dmu = extract(1), Dmu2 = extract(2), EDmu2 = expected[, 1])
  if (level > 0) {
    r$Dth <- first(0); r$Dmuth <- first(1); r$Dmu3 <- extract(3)
    r$Dmu2th <- first(2)
    r$EDmu3 <- expected[, match(paste(c(1, integer(nt)), collapse=","), fbasis$keys)]
    r$EDmu2th <- matrix(vapply(seq_len(nt), function(j) expected[, match(paste(c(0, unit[j, ]), collapse=","), fbasis$keys)], numeric(length(y))), length(y), nt)
  }
  if (level > 1) {
    r$Dmu4 <- extract(4); r$Dth2 <- second(0); r$Dmuth2 <- second(1)
    r$Dmu2th2 <- second(2); r$Dmu3th <- first(3)
  }
  r
}
