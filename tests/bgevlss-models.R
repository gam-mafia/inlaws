library("inlaws")

near <- function(a, b, tol = 1e-7) {
  stopifnot(identical(dim(a), dim(b)), all(abs(a - b) <= tol*(1 + abs(b))))
}

set.seed(405)
n <- 240L
d <- data.frame(x = runif(n), z = runif(n), o = runif(n, -.12, .12))
f <- bgevlss()
d$y <- f$rd(cbind(1.5*sin(5*d$x), exp(d$o + .2*d$z), exp(-1.3 + d$o)))
form <- list(y ~ s(x, k = 6) + offset(o), ~ z + offset(o), ~ offset(o))

validate_prediction <- function(b, data) {
  mu <- fitted(b)
  stopifnot(identical(colnames(mu), c("location", "spread", "shape")),
            all(is.finite(mu)), all(mu[, 2:3] > 0))
  link <- predict(b, newdata = data, type = "link", se.fit = TRUE)
  response <- predict(b, newdata = data, type = "response", se.fit = TRUE)
  near(response$fit, mu)
  near(response$fit, cbind(link$fit[, 1], exp(link$fit[, 2:3])))
  near(response$se.fit, link$se.fit*cbind(1, mu[, 2:3]), 2e-6)
  lp <- f$ll(data$y, NULL, NULL, b$prior.weights, b$family, eta = link$fit)$l
  near(as.numeric(logLik(b)), lp, 5e-6)
}

# Reconstruct the Laplace REML criterion in the original coefficient basis.
setup <- gam(form, family = f, data = d, fit = FALSE)
reml_score <- function(b) {
  X <- predict(b, type = "lpmatrix")
  ll <- f$ll(d$y, X, coef(b), b$prior.weights, f, b$offset, deriv = 1)
  S <- matrix(0, ncol(X), ncol(X))
  ix <- setup$off[1L] + seq_len(ncol(setup$S[[1L]])) - 1L
  S[ix, ix] <- b$sp[1L]*setup$S[[1L]]
  eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  rank <- setup$rank[1L]
  log_det_S <- sum(log(eig[seq_len(rank)]))
  log_det_H <- 2*sum(log(diag(chol(S - ll$lbb))))
  penalty <- drop(crossprod(coef(b), S %*% coef(b)))
  -ll$l + (penalty + log_det_H - log_det_S - (ncol(X)-rank)*log(2*pi))/2
}

opts <- list(c("outer", "newton"), c("outer", "bfgs"), "efs")
reml <- lapply(opts, function(opt) {
  b <- gam(form, family = f, data = d, method = "REML", optimizer = opt)
  stopifnot(identical(b$method, "REML"), identical(b$optimizer, opt),
            identical(b$outer.info$conv, "full convergence"),
            is.finite(b$gcv.ubre), all(is.finite(b$sp) & b$sp > 0))
  validate_prediction(b, d)
  near(as.numeric(b$gcv.ubre), reml_score(b))
  b
})
near(as.numeric(reml[[1]]$gcv.ubre), as.numeric(reml[[2]]$gcv.ubre), 1e-6)
near(as.numeric(reml[[1]]$gcv.ubre), as.numeric(reml[[3]]$gcv.ubre), 1e-4)

for (opt in c("nlm", "optim")) {
  msg <- tryCatch(gam(form, family = f, data = d, method = "REML",
                     optimizer = c("outer", opt)), error = conditionMessage)
  stopifnot(is.character(msg), grepl("required objective scale", msg))
}

if (requireNamespace("gratia", quietly = TRUE) &&
    "quantile_residuals" %in% getNamespaceExports("gratia")) {
  fit <- reml[[1L]]
  pit <- f$cdf(d$y, fitted(fit))
  near(as.numeric(gratia::quantile_residuals(fit, type = "pit")), pit)
}

# Verify the actual held-out criterion using the returned deletion predictions.
# QNCV uses a quadratic approximation about each full-data predictor vector.
validate_cv <- function(b, method, idx, gamma) {
  stopifnot(identical(b$method, method),
            identical(b$optimizer, c("outer", "bfgs")),
            isTRUE(b$family$qapprox == (method == "QNCV")),
            is.finite(b$gcv.ubre))
  eta <- b$linear.predictors[idx, , drop = FALSE]
  cv <- attr(b$gcv.ubre, "eta.cv")
  w <- b$prior.weights[idx]
  stopifnot(identical(dim(cv), dim(eta)), all(is.finite(cv)))
  base <- f$ll(b$y[idx], NULL, NULL, w, f, eta = eta, deriv = 1, ncv = TRUE)
  if (method == "NCV") {
    exact <- f$ll(b$y[idx], NULL, NULL, w, f, eta = cv)$l
    expected <- -(gamma*exact + (1-gamma)*base$l)
  } else {
    delta <- cv - eta
    h_delta <- vapply(1:3, function(j)
      rowSums(base$l2[, f$tri$i2[j, ]] * delta), numeric(length(idx)))
    expected <- -base$l - gamma*sum(delta * (base$l1 + h_delta/2))
  }
  near(as.numeric(b$gcv.ubre), expected)
}

cv <- lapply(c("NCV", "QNCV"), function(method) {
  b <- gam(form, family = f, data = d, method = method)
  stopifnot(identical(b$outer.info$conv, "full convergence"))
  validate_cv(b, method, seq_len(n), 1)
  validate_prediction(b, d)
  b
})
stopifnot(abs(cv[[1]]$gcv.ubre - cv[[2]]$gcv.ubre) > .01)

# Explicit paired deletions, reordered predictions, fractional and zero weights.
# There is a fitted shape slope as well as offsets in all three predictors.
a <- seq_len(n)
idx <- c(seq(2, n, by = 2), seq(1, n, by = 2))
nei <- list(a = a, ma = seq(2, n, by = 2), d = idx,
            md = seq(2, n, by = 2), jackknife = FALSE)
d$w <- rep(c(0, .5, 2, 1, 1, 1), length.out = n)
form[[3L]] <- ~ x + offset(o)
for (method in c("NCV", "QNCV")) {
  b <- gam(form, family = f, data = d, method = method, nei = nei, gamma = 1.1,
           weights = w)
  validate_cv(b, method, idx, 1.1)
}
msg <- tryCatch(gam(form, family = f, data = d, method = "NCV", optimizer = "efs"),
                error = conditionMessage)
stopifnot(is.character(msg), grepl("requires BFGS", msg))

# Unpenalized NCV/QNCV evaluates the criterion with covariance disabled.
linear <- list(y ~ x + offset(o), ~ z + offset(o), ~ offset(o))
for (method in c("NCV", "QNCV")) {
  b <- gam(linear, family = f, data = d, method = method,
           nei = list(jackknife = FALSE))
  validate_cv(b, method, seq_len(n), 1)
  if (utils::packageVersion("mgcv") < "1.9.5") {
    msg <- tryCatch(gam(linear, family = f, data = d, method = method),
                    error = conditionMessage)
    stopifnot(is.character(msg), grepl("jackknife covariance", msg))
  }
}

# A plateau creates binding inequalities in the monotone location smooth.
if ("scasm" %in% getNamespaceExports("mgcv")) {
  set.seed(405)
  d <- data.frame(x = runif(n), z = runif(n), o = runif(n, -.12, .12))
  d$y <- f$rd(cbind(4*pmin(d$x, .45) + d$o, exp(d$o + .2*d$z), exp(-1.3 + d$o)))
  constrained <- list(y ~ s(x, bs = "sc", xt = "m+", k = 7) + offset(o),
                      ~ z + offset(o), ~ offset(o))
  b <- mgcv::scasm(constrained, family = f, data = d)
  slack <- drop(b$Ain %*% coef(b) - b$bin)
  stopifnot(identical(b$method, "REML"), identical(b$optimizer, "EFS"),
            isTRUE(b$outer.info$converged), b$outer.info$niter > 1,
            is.finite(b$gcv.ubre), identical(unname(b$gcv.ubre), -b$laml),
            all(slack >= -1e-8), any(slack < 1e-8), any(slack > .1),
            length(attr(coef(b), "active")) > 0L)
  grid <- data.frame(x = seq(0, 1, length.out = 200), z = .5, o = 0)
  fitted_location <- predict(b, newdata = grid, type = "response")[, 1L]
  stopifnot(min(diff(fitted_location)) > -1e-8, diff(range(fitted_location)) > 1)
  validate_prediction(b, d)

  for (method in c("NCV", "QNCV")) {
    msg <- tryCatch(mgcv::scasm(constrained, family = f, data = d, method = method),
                    error = conditionMessage)
    stopifnot(is.character(msg), grepl("REML/EFS", msg, fixed = TRUE))
  }
}
