library("inlaws")

near <- function(x, y, tol = 1e-7) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}
fails <- function(expr, pattern) {
  err <- tryCatch(expr, error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err)))
}
reference <- function(y, mu, k) dweibull(y, k, mu / gamma(1 + 1/k), log = TRUE)
check <- function(g, method, optimizer, general = TRUE, report_optimizer = TRUE) {
  stopifnot(identical(g$method, method), is.finite(g$gcv.ubre),
            all(is.finite(coef(g))), all(is.finite(g$Vp)), g$sig2 == 1,
            is.finite(g$null.deviance))
  if (report_optimizer) stopifnot(identical(g$optimizer, optimizer))
  if (identical(optimizer, "EFS")) stopifnot(isTRUE(g$outer.info$converged))
  else if (length(g$outer.info[["conv"]]))
    stopifnot(identical(g$outer.info$conv, "full convergence"))
  fit <- fitted(g)
  mu <- if (general) fit[, 1L] else fit
  k <- if (general) fit[, 2L] else g$family$getTheta(TRUE)
  near(logLik(g), sum(g$prior.weights * reference(g$y, mu, k)))
  near(g$deviance, sum(residuals(g, "deviance")^2))
  near(residuals(g, "response"), g$y - mu)
  variance <- mu^2 * (gamma(1 + 2/k) / gamma(1 + 1/k)^2 - 1)
  near(residuals(g, "pearson"), (g$y - mu) * sqrt(g$prior.weights / variance))
  lp <- predict(g, se.fit = TRUE)
  rp <- predict(g, type = "response", se.fit = TRUE)
  near(rp$fit, exp(lp$fit))
  near(rp$fit, fit)
  near(rp$se.fit, rp$fit * lp$se.fit)
}

set.seed(78)
n <- 220L
d <- data.frame(x = runif(n), z = runif(n), a = runif(n, -.1, .1),
                b = runif(n, -.1, .1), w = rep(c(0, .75, 1, 1.25), length.out = n))
f <- weibulls()
d$y <- f$rd(cbind(exp(sin(6*d$x) + d$a), exp(.4 + cos(6*d$z)*.5 + d$b)))
form <- y ~ s(x, k = 6) + offset(a)
forms <- list(form, ~ s(z, k = 6) + offset(b))
for (model in list(forms, list(form, ~ 1))) {
  for (optimizer in list(c("outer", "newton"), c("outer", "bfgs"), "efs")) {
    g <- gam(model, data = d, weights = w, family = f, method = "REML", optimizer = optimizer)
    check(g, "REML", optimizer)
    if (identical(optimizer, c("outer", "newton")) && identical(model, forms)) reml <- g
  }
}
stopifnot(identical(colnames(fitted(reml)), c("mu", "k")))

# Distinct deletion and prediction groups, repeated prediction rows, both
# offsets, non-unit gamma, and optimization of two smoothing parameters.
nei <- list(a = seq_len(n), ma = seq.int(4L, n, by = 4L),
  d = rep(seq.int(3L, n, by = 4L), each = 2L),
  md = seq.int(2L, n / 2L, by = 2L), jackknife = FALSE)
eval_weibull <- getFromNamespace(".weibull_eval", "inlaws")
for (method in c("NCV", "QNCV")) {
  loo <- gam(forms, data = d, weights = w, family = f, method = method)
  check(loo, method, c("outer", "bfgs"))
  g <- gam(forms, data = d, weights = w, family = f, method = method,
           nei = nei, gamma = 1.2)
  check(g, method, c("outer", "bfgs"))
  cv <- attr(g$gcv.ubre, "eta.cv")
  i <- nei$d
  eta <- g$linear.predictors[i, , drop = FALSE]
  ell <- eval_weibull(d$y[i], eta[, 1L], eta[, 2L], 2L)
  if (method == "NCV") {
    cv_ell <- reference(d$y[i], exp(cv[, 1L]), exp(cv[, 2L]))
  } else {
    delta <- cv - eta
    cv_ell <- ell[, 1L] + rowSums(ell[, 2:3] * delta) +
      .5 * rowSums(ell[, 4:6] * cbind(delta[, 1L]^2, 2*delta[, 1L]*delta[, 2L], delta[, 2L]^2))
  }
  near(g$gcv.ubre, -sum(d$w[i] * (ell[, 1L] + 1.2 * (cv_ell - ell[, 1L]))))
}

if (packageVersion("mgcv") == "1.9.4") {
  for (method in c("NCV", "QNCV")) {
    unpenalized <- list(y ~ x + offset(a), ~ z + offset(b))
    fails(gam(unpenalized, data = d, family = f, method = method), "jackknife = FALSE")
    g <- gam(unpenalized, data = d, family = f, method = method, nei = nei)
    stopifnot(identical(g$method, method), is.finite(g$gcv.ubre))
  }
}

# Integer replication, zero rows, and agreement with the single-predictor
# likelihood at matched mean and fixed shape. Compare REML only within fits
# that handle the shape parameter in the same way.
d$w <- rep(c(0, 1, 2, 3), length.out = n)
replicated <- d[rep(seq_len(n), d$w), ]
model <- list(y ~ x + offset(a), ~ 1)
g <- gam(model, data = d, weights = w, family = f)
g_rep <- gam(model, data = replicated, family = f)
near(coef(g), coef(g_rep), 1e-5)
p <- fitted(g)
fixed <- weibull(p[1L, 2L])
near(logLik(g), -fixed$aic(d$y, p[, 1L], wt = d$w, dev = 0) / 2)
near(g$deviance, sum(fixed$dev.resids(d$y, p[, 1L], d$w)))
near(fixed$cdf(d$y, p[, 1L]), f$cdf(d$y, p))

if ("scasm" %in% getNamespaceExports("mgcv")) {
  set.seed(34)
  m <- 260L
  dat <- data.frame(x = runif(m), z = runif(m), a = runif(m, -.1, .1),
                    b = runif(m, -.1, .1), w = rep(c(0, .75, 1, 1.25), length.out = m))
  # The data-generating curve turns down, so the increasing constraint binds.
  dat$y <- f$rd(cbind(exp(2*dat$x - 3*pmax(dat$x - .55, 0) + dat$a),
                     exp(.5 + .7*dat$z + dat$b)))
  constrained <- list(y ~ s(x, bs = "sc", xt = "m+", k = 8) + offset(a), ~ z + offset(b))
  g <- mgcv::scasm(constrained, data = dat, weights = w, family = f)
  check(g, "REML", "EFS")
  near(g$gcv.ubre, -g$laml)
  stopifnot(nrow(g$Ain) > 0L, length(attr(coef(g), "active")) > 0L,
            all(drop(g$Ain %*% coef(g)) >= g$bin - 1e-8))
  at <- data.frame(x = seq(0, 1, length.out = 401), z = .5, a = 0, b = 0)
  rp <- predict(g, newdata = at, type = "response", se.fit = TRUE)
  stopifnot(all(diff(rp$fit[, 1L]) >= -1e-8), all(is.finite(rp$se.fit)))
}

# Extended family routes and fitted shape state.
set.seed(81)
n <- 240L
d <- data.frame(x = runif(n), a = runif(n, -.15, .15),
                w = rep(c(0, .75, 1, 1.25), length.out = n))
d$y <- weibull(1.6)$rd(exp(.4 + sin(6*d$x) + d$a))
for (shape in list(NULL, 1.6)) {
  f <- weibull(shape)
  for (opt in list(c("outer", "newton"), c("outer", "bfgs"), "efs",
                   c("outer", "nlm"), c("outer", "optim"))) {
    g <- gam(form, data = d, weights = w, family = f, method = "REML", optimizer = opt)
    check(g, "REML", opt, general = FALSE)
    if (identical(opt, c("outer", "nlm"))) stopifnot(g$outer.info$code %in% 1:2)
    if (identical(opt, c("outer", "optim"))) stopifnot(g$outer.info$convergence == 0L)
    if (!is.null(shape)) near(g$family$getTheta(TRUE), shape)
  }
  g <- gam(form, data = d, weights = w, subset = w > 0, family = f, method = "ML")
  check(g, "ML", c("outer", "newton"), general = FALSE)
  grouped <- list(a = seq_len(n), ma = seq.int(3L, n, by = 3L),
    d = seq_len(n), md = seq.int(3L, n, by = 3L), jackknife = FALSE)
  for (method in c("NCV", "QNCV")) for (groups in list(NULL, grouped)) {
    g <- gam(form, data = d, weights = w, family = f, method = method,
             optimizer = c("outer", "bfgs"), nei = groups)
    check(g, method, c("outer", "bfgs"), general = FALSE)
    cv <- attr(g$gcv.ubre, "eta.cv")
    k <- g$family$getTheta(TRUE)
    if (method == "NCV") {
      ell <- reference(d$y, exp(cv), k)
    } else {
      eta <- g$linear.predictors
      at <- eval_weibull(d$y, eta, log(k), 2L)
      delta <- cv - eta
      ell <- at[, 1L] + at[, 2L] * delta + .5 * at[, 4L] * delta^2
    }
    near(g$gcv.ubre, -sum(d$w * ell))
  }
  for (method in c("REML", "ML", "fREML")) {
    opt <- if (method == "fREML") c("perf", "newton") else c("outer", "newton")
    g <- bam(form, data = d, weights = w, family = f, method = method, scale = 1)
    check(g, method, opt, general = FALSE)
    stopifnot(is.null(g$dinfo), !g$scale.estimated)
    if (method == "fREML") ordinary <- g
  }
  methods <- if (packageVersion("mgcv") >= "1.9.4") c("fREML", "NCV") else "fREML"
  for (method in methods) {
    g <- bam(form, data = d, weights = w, family = f, method = method,
             discrete = TRUE, scale = 1)
    check(g, method, c("perf", "chol"), general = FALSE)
    stopifnot(!is.null(g$dinfo), isTRUE(g$mgcv.conv), !g$scale.estimated)
    if (method == "fREML") near(fitted(g), fitted(ordinary), .005)
    default <- bam(form, data = d, weights = w, family = f, method = method, discrete = TRUE)
    near(coef(default), coef(g))
    near(default$family$getTheta(), g$family$getTheta())
    near(default$gcv.ubre, g$gcv.ubre)
    stopifnot(default$sig2 == 1, !is.null(default$dinfo), identical(default$method, method))
  }
  near(f$getTheta(TRUE), if (is.null(shape)) 1 else shape)
}
if (packageVersion("mgcv") <= "1.9.5")
  fails(gam(form, data = d, weights = w, family = weibull(), method = "ML"), "zero weights")
fails(bam(form, data = d, family = weibull(), scale = 2), "fixed at one")
fails(gam(form, data = d, family = weibull(), scale = 2), "fixed at one")
fails(gam(list(form, ~ 1), data = d, family = weibulls(), scale = 2), "fixed at one")
for (response in list(c(0, 1), c(-1, 1), c(1, Inf), cbind(1:3, 1:3))) {
  fails(gam(response ~ 1, family = weibull()), "positive, uncensored")
  fails(gam(list(response ~ 1, ~ 1), family = weibulls()), "positive, uncensored")
}

# At fixed shape one the single-predictor fit is the Gamma fit at scale one.
g <- gam(form, data = d, weights = w, family = weibull(1), method = "REML")
e <- gam(form, data = d, weights = w, family = Gamma(link = "log"), method = "REML", scale = 1)
near(coef(g), coef(e), 1e-5)
near(logLik(g), sum(d$w * dexp(d$y, rate = 1/fitted(g), log = TRUE)))

if (requireNamespace("gratia", quietly = TRUE)) {
  has_factory <- any(vapply(ls(asNamespace("gratia"), all.names = TRUE), function(nm) {
    obj <- get(nm, asNamespace("gratia"))
    is.function(obj) && any(grepl("gratia_diagnostics", deparse(body(obj)), fixed = TRUE))
  }, logical(1)))
  if (has_factory) for (fit in list(g, reml)) {
    pit <- gratia::quantile_residuals(fit, type = "pit", seed = 21)
    stopifnot(length(pit) == length(fit$y), all(pit >= 0 & pit <= 1))
    stopifnot(inherits(gratia::qq_plot(fit, method = "simulate", type = "deviance",
                                      n_simulate = 3, seed = 21), "ggplot"))
    stopifnot(inherits(gratia::worm_plot(fit, method = "simulate", type = "pearson",
                                        n_simulate = 3, seed = 21), "ggplot"))
  }
}
