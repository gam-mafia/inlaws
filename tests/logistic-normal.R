library("inlaws")

near <- function(actual, expected, tol = 1e-6) {
  stopifnot(identical(dim(actual), dim(expected)),
    all(is.finite(actual)), all(is.finite(expected)),
    max(abs(actual - expected) / (1 + abs(expected))) < tol)
}

fails <- function(expr, pattern) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}

inverse_alr <- function(z, reference = 1L) {
  components <- seq_len(ncol(z) + 1L)
  logs <- matrix(0, nrow(z), length(components))
  logs[, components[-reference]] <- z
  values <- exp(logs - apply(logs, 1L, max))
  values / rowSums(values)
}

alr <- function(y, reference = 1L) log(y[, -reference, drop = FALSE] / y[, reference])

derivatives <- getFromNamespace(".logistic_normal_derivatives", "inlaws")

# The fully normalized weighted likelihood uses the density on the simplex.
y <- rbind(c(.26, .33, .41), c(.15, .64, .21), c(.08, .2, .72))
z <- alr(y, 2L)
jacobian <- -rowSums(log(y))
eta <- rbind(c(.2, -.3, log(.4), log(.75), .9),
             c(-.1, .7, log(.55), log(.9), -.4),
             c(.5, -.2, log(.7), log(.3), .25))
weights <- c(.7, 2, 0)
reference <- vapply(seq_len(nrow(y)), function(i) {
  sds <- exp(eta[i, 3:4])
  rho <- eta[i, 5L] / sqrt(1 + eta[i, 5L]^2)
  covariance <- outer(sds, sds) * matrix(c(1, rho, rho, 1), 2L)
  residual <- z[i, ] - eta[i, 1:2]
  weights[i] * (jacobian[i] - log(2 * pi) -
    as.numeric(determinant(covariance, logarithm = TRUE)$modulus) / 2 -
    drop(crossprod(residual, solve(covariance, residual))) / 2)
}, numeric(1L))
tri <- mgcv::trind.generator(5L)
value <- derivatives(z, jacobian, eta, weights, tri, 0L)
near(as.numeric(value$l0), reference)
near(as.numeric(value$l), sum(reference))

# The penalty derivatives must include mixed partials for entries in the same
# Cholesky row and vanish across distinct rows.
lkj_derivatives <- getFromNamespace(".logistic_normal_lkj", "inlaws")
theta <- c(.35, .8, -.55)
lkj <- 2.75
raw <- rbind(c(1, 0, 0), c(theta[1L], 1, 0), c(theta[2:3], 1))
correlation <- tcrossprod(raw / sqrt(rowSums(raw^2)))
log_det <- as.numeric(determinant(correlation, logarithm = TRUE)$modulus)
lkj_at <- function(theta, order) lkj_derivatives(theta, 3L, lkj, order)
penalty <- lkj_at(theta, 3L)
near(penalty$l, (lkj - 1) * log_det)
for (j in seq_along(theta)) {
  step <- numeric(length(theta))
  step[j] <- 1e-5
  a <- lkj_at(theta + step, 2L)
  b <- lkj_at(theta - step, 2L)
  near(penalty$gradient[j], (a$l - b$l) / (2 * step[j]))
  near(penalty$hessian[, j], (a$gradient - b$gradient) / (2 * step[j]))
  near(penalty$third[, , j], (a$hessian - b$hessian) / (2 * step[j]))
}

# Check all packed partials through third order, including interactions with
# both the correlation and scale predictors. Use unequal weights across rows.
set.seed(295)
for (m in 1:3) {
  q <- 2L * m + choose(m, 2L)
  tri <- mgcv::trind.generator(q)
  z <- matrix(rnorm(2L * m, sd = .8), 2L)
  eta <- matrix(rnorm(2L * q, sd = .35), 2L)
  jacobian <- -rowSums(log(inverse_alr(z)))
  evaluate <- function(eta, order) derivatives(z, jacobian, eta, c(.8, 1.5), tri, order)
  d <- evaluate(eta, 3L)
  h <- 1e-5
  for (j in seq_len(q)) {
    plus <- minus <- eta
    plus[, j] <- plus[, j] + h
    minus[, j] <- minus[, j] - h
    a <- evaluate(plus, 2L)
    b <- evaluate(minus, 2L)
    near(d$l1[, j], as.numeric(a$l0 - b$l0) / (2 * h))
    near(d$l2[, tri$i2[j, ], drop = FALSE], (a$l1 - b$l1) / (2 * h))
    for (k in seq_len(q)) {
      near(d$l3[, tri$i3[j, k, ], drop = FALSE],
        (a$l2[, tri$i2[k, ], drop = FALSE] - b$l2[, tri$i2[k, ], drop = FALSE]) / (2 * h))
    }
  }
}

# Fit smooth means and a varying log-SD with a named nonfirst reference.
set.seed(315)
n <- 180L
dat <- data.frame(x = runif(n), o1 = rnorm(n, sd = .08), o2 = rnorm(n, sd = .08),
  os = rnorm(n, sd = .04))
mu <- with(dat, cbind(sin(2 * pi * x) + o1, -.6 * cos(2 * pi * x) + o2))
sds <- with(dat, cbind(exp(-.9 + .3 * x + os), .5))
errors <- cbind(rnorm(n), rnorm(n)) %*% chol(matrix(c(1, .45, .45, 1), 2L))
dat$y <- inverse_alr(mu + errors * sds, 2L)
colnames(dat$y) <- c("a", "b", "c")
family <- logistic_normal(3L, refcat = "b")
simple_form <- logistic_normal_formula(y ~ x, family)
varying_form <- logistic_normal_formula(y ~ x, family, scale = ~ x)
stopifnot(identical(vapply(simple_form, deparse1, character(1L)),
  c("y ~ x", "~x", "~1", "~1", "~1")),
  identical(vapply(varying_form, deparse1, character(1L)),
    c("y ~ x", "~x", "~x", "~x", "~1")))
form <- list(y ~ s(x, k = 6) + offset(o1), ~ s(x, k = 6) + offset(o2),
  ~ x + offset(os), ~ 1, ~ 1)
fit <- gam(form, data = dat, family = family,
  method = "REML", optimizer = c("outer", "bfgs"))
lp <- predict(fit, type = "link")
center <- predict(fit, type = "response", se.fit = TRUE)
stopifnot(identical(fit$method, "REML"), identical(dim(lp), c(n, 5L)),
  sqrt(mean((lp[, 1:2] - mu)^2)) < .2, all(center$se.fit > 0),
  is.finite(as.numeric(logLik(fit))))
near(unname(center$fit), inverse_alr(lp[, 1:2], 2L))
near(unname(predict(fit, newdata = dat[1:5, ], type = "response")),
  unname(center$fit[1:5, , drop = FALSE]))
rho <- lp[1L, 5L] / sqrt(1 + lp[1L, 5L]^2)
near(fit$family$data$correlation, matrix(c(1, rho, rho, 1), 2L))
stopifnot(identical(colnames(fit$family$data$correlation), c("a", "c")),
  length(residuals(fit, type = "pearson")) == n)

diagnostics <- fit$family$gratia_diagnostics(lp, NULL, 1, fit$y, lp)
stopifnot(!diagnostics$mean_available, !diagnostics$independent,
  identical(dim(diagnostics$rd()), c(n, 3L)))
near(diagnostics$diagnostic_residuals(fit$y, "response"), fit$y - center$fit)

draw <- fit$family$rd(lp, NULL, 1)
stopifnot(identical(dim(draw), c(n, 3L)), all(draw > 0))
near(rowSums(draw), rep(1, n), 1e-13)

# A nonfirst reference is used consistently by both fitted and simulated values.
set.seed(842)
sim_eta <- matrix(rep(c(.4, -.7, log(.5), log(.8), .6), each = 5000L), 5000L)
log_ratios <- alr(fit$family$rd(sim_eta, NULL, 1), 2L)
near(colMeans(log_ratios), c(.4, -.7), .03)
near(cov(log_ratios), matrix(c(.25, .4 * .6 / sqrt(1 + .6^2),
  .4 * .6 / sqrt(1 + .6^2), .64), 2L), .03)

# The two-part reduction agrees with a Gaussian likelihood on the sole log
# ratio, and prior weights act as likelihood multipliers.
small <- dat
small$y <- dat$y[, c(1, 3)] / rowSums(dat$y[, c(1, 3)])
small$z <- drop(alr(small$y))
simple <- gam(list(y ~ x, ~ 1), data = small, family = logistic_normal(), method = "REML")
gaussian <- lm(z ~ x, data = small)
near(unname(coef(simple)[1:2]), unname(coef(gaussian)), 2e-5)
sigma_ml <- sqrt(mean(residuals(gaussian)^2))
near(unname(exp(coef(simple)[3L])), sigma_ml, 2e-5)
near(as.numeric(logLik(simple)), sum(dnorm(small$z, fitted(gaussian), sigma_ml, log = TRUE)) -
  sum(log(small$y)), 2e-5)

dat$w <- rep(c(0, 1, 2), length.out = n)
linear <- list(y ~ x + offset(o1), ~ x + offset(o2), ~ offset(os), ~ 1, ~ 1)
weighted <- gam(linear, data = dat, weights = w, family = logistic_normal(3L, 2L), method = "REML")
replicated <- gam(linear, data = dat[rep(seq_len(n), dat$w), ],
  family = logistic_normal(3L, "b"), method = "REML")
near(unname(coef(weighted)), unname(coef(replicated)), 2e-5)
near(as.numeric(logLik(weighted)), as.numeric(logLik(replicated)), 2e-6)

# A stronger LKJ penalty shrinks the residual correlation. The reported
# likelihood and AIC continue to use the normalized data likelihood.
regularized <- gam(form, data = dat, family = logistic_normal(3L, 2L, lkj = 30),
  method = "REML", optimizer = c("outer", "bfgs"))
stopifnot(identical(regularized$method, "REML"),
  abs(regularized$family$data$correlation[2L, 1L]) < abs(fit$family$data$correlation[2L, 1L]))
regularized_eta <- predict(regularized, type = "link")
regularized_data_loglik <- derivatives(alr(dat$y, 2L), -rowSums(log(dat$y)),
  regularized_eta, rep(1, n), NULL, 0L)$l
reported <- logLik(regularized)
near(as.numeric(reported), regularized_data_loglik)
near(as.numeric(AIC(regularized)), -2 * regularized_data_loglik + 2 * attr(reported, "df"))

# The neighbourhoods delete whole compositions, including reordered and
# repeated held-out rows. Validate criterion derivatives with offsets and gamma.
nei <- list(a = c(4, 5, 6, 20, 21, 22, 4, 5, 6, 39, 40, 41),
  ma = c(3, 6, 9, 12), d = c(6, 22, 21, 5, 40), md = c(1, 3, 4, 5), jackknife = FALSE)
X <- predict(fit, type = "lpmatrix")
beta <- coef(fit)
db <- matrix(seq(-.15, .15, length.out = length(beta)), ncol = 1L)
offset <- list(dat$o1, dat$o2, dat$os, numeric(n), numeric(n))
weights <- seq(.7, 1.3, length.out = n)
score <- function(step, want_gradient, family) {
  b <- beta + step * db[, 1L]
  ll <- family$ll(fit$y, X, b, weights, family, offset = offset,
    deriv = if (want_gradient) 3L else 1L, d1b = db, ncv = TRUE)
  ll$gamma <- 1.3
  H <- -ll$lbb + diag(10, length(b))
  family$ncv(X, fit$y, weights, nei, b, family, ll, R = chol(H),
    offset = offset, dH = ll$d1H, db = db, deriv = want_gradient)
}

# Assemble observation-level beta derivatives independently of the callback.
z <- alr(dat$y, 2L)
jacobian <- -rowSums(log(dat$y))
tri <- mgcv::trind.generator(5L)
d <- derivatives(z, jacobian, lp, weights, tri, 2L)
indices <- attr(X, "lpi")
design <- lapply(seq_len(n), function(i) {
  rows <- matrix(0, 5L, ncol(X))
  for (j in seq_len(nrow(rows))) rows[j, indices[[j]]] <- X[i, indices[[j]]]
  rows
})
hessians <- Map(function(A, i) {
  curvature <- matrix(d$l2[i, as.vector(tri$i2)], 5L)
  crossprod(A, curvature %*% A)
}, design, seq_len(n))
scores <- Map(function(A, i) crossprod(A, d$l1[i, ]), design, seq_len(n))
loglik_hessian <- Reduce(`+`, hessians)
full_score <- -sum(d$l0[nei$d])

# For D = 3, log(det(R)) = -log(1 + theta^2). The penalty is applied
# once to the full coefficient likelihood, with curvature retained when
# observations are deleted. Neither held-out score includes the penalty.
theta <- lp[1L, 5L]
theta_design <- design[[1L]][5L, ]
for (lkj in c(1, 30)) {
  family <- if (lkj == 1) fit$family else regularized$family
  strength <- lkj - 1
  penalty <- -strength * log1p(theta^2)
  penalty_gradient <- -2 * strength * theta / (1 + theta^2) * theta_design
  penalty_hessian <- -2 * strength * (1 - theta^2) / (1 + theta^2)^2 *
    tcrossprod(theta_design)
  ll <- family$ll(fit$y, X, beta, weights, family, offset = offset, deriv = 1L, ncv = TRUE)
  near(ll$l, d$l + penalty)
  near(ll$l0, d$l0)
  near(ll$l1, d$l1)
  near(ll$l2, d$l2)
  near(as.numeric(ll$lb), as.numeric(Reduce(`+`, scores)) + penalty_gradient)
  near(ll$lbb, loglik_hessian + penalty_hessian)
  supplied <- family$ll(fit$y, X, beta, weights, family, eta = lp, deriv = 1L, ncv = TRUE)
  near(supplied$l, d$l)
  near(supplied$l0, d$l0)
  near(supplied$l1, d$l1)
  near(supplied$l2, d$l2)

  H <- diag(10, ncol(X)) - loglik_hessian - penalty_hessian
  V <- chol2inv(chol(H))
  trace_derivative <- family$ll(fit$y, X, beta, weights, family, offset = offset,
    deriv = 2L, d1b = db, fh = V)$d1H
  curvature_at <- function(step) family$ll(fit$y, X, beta + step * db[, 1L],
    weights, family, offset = offset, deriv = 1L)$lbb
  numeric_trace <- sum(V * (curvature_at(1e-4) - curvature_at(-1e-4))) / 2e-4
  near(as.numeric(trace_derivative), numeric_trace)

  expected_eta <- matrix(0, length(nei$d), ncol(lp))
  for (i in seq_along(nei$ma)) {
    omitted <- nei$a[(c(0L, nei$ma)[i] + 1L):nei$ma[i]]
    pos <- (c(0L, nei$md)[i] + 1L):nei$md[i]
    step <- solve(H + Reduce(`+`, hessians[omitted]), Reduce(`+`, scores[omitted]))
    for (j in pos) {
      row <- nei$d[j]
      expected_eta[j, ] <- lp[row, ] - drop(design[[row]] %*% step)
    }
  }
  delta <- expected_eta - lp[nei$d, , drop = FALSE]
  quadratic_change <- vapply(seq_along(nei$d), function(i) {
    row <- nei$d[i]
    curvature <- matrix(d$l2[row, as.vector(tri$i2)], 5L)
    sum(delta[i, ] * d$l1[row, ]) + drop(crossprod(delta[i, ], curvature %*% delta[i, ])) / 2
  }, numeric(1L))
  held_score <- -derivatives(z[nei$d, , drop = FALSE], jacobian[nei$d], expected_eta,
    weights[nei$d], tri, 0L)$l

  for (method in c("NCV", "QNCV")) {
    family$qapprox <- method == "QNCV"
    value <- score(0, TRUE, family)
    difference <- (score(1e-4, FALSE, family)$NCV - score(-1e-4, FALSE, family)$NCV) / 2e-4
    near(as.numeric(value$NCV1), as.numeric(difference), 1e-5)
    near(attr(value$NCV, "eta.cv"), expected_eta, 1e-5)
    expected_score <- if (family$qapprox) full_score - 1.3 * sum(quadratic_change) else
      (1 - 1.3) * full_score + 1.3 * held_score
    near(as.numeric(value$NCV), expected_score, 1e-5)
    stopifnot(identical(dim(value$Vg), c(1L, 1L)))
    cv_fit <- gam(form, data = dat, family = logistic_normal(3L, 2L, lkj = lkj),
      method = method, nei = nei, gamma = 1.3)
    stopifnot(identical(cv_fit$method, method), is.finite(cv_fit$gcv.ubre),
      length(cv_fit$sp) == 2L, all(is.finite(coef(cv_fit))))
    if (lkj == 1) {
      loo_fit <- gam(form, data = dat, family = logistic_normal(3L, 2L), method = method)
      stopifnot(identical(loo_fit$method, method), is.finite(loo_fit$NCV),
        all(is.finite(vcov(loo_fit))))
    }
  }
}

for (bad in list(0, 1, 2.5, NA_real_, Inf, c(2, 3), "3"))
  fails(logistic_normal(bad), "integer|component|part")
for (bad in list(-1, .8, NA_real_, Inf, numeric(), c(1, 2), "2"))
  fails(logistic_normal(lkj = bad), "lkj")
fails(logistic_normal(3L, 4L), "refcat|reference")
fails(gam(linear, data = dat, family = logistic_normal(3L, "missing")), "refcat|reference")
fails(gam(form[-5L], data = dat, family = logistic_normal(3L)), "formula|predictor")
bad <- dat
bad$y[1L, ] <- c(0, .3, .7)
fails(gam(linear, data = bad, family = logistic_normal(3L)), "positive|between|zero")
bad$y <- .9 * dat$y
fails(gam(linear, data = bad, family = logistic_normal(3L)), "sum to one|sum to 1")
bad$y <- dat$y[, 1:2]
fails(gam(linear, data = bad, family = logistic_normal(3L)), "columns|components|parts")
bad <- linear
bad[[5L]] <- ~ x
fails(gam(bad, data = dat, family = logistic_normal(3L, 2L)), "intercept|correlation")
