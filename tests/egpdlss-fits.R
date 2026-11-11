library("inlaws")

near <- function(x, y, tolerance = 1e-7) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)))
}

# Real smoothing selection with offsets and nonuniform weights. Record the
# optimizer entry and actual criterion because gam() can echo the request.
f <- egpdlss()
set.seed(1203)
d <- data.frame(x = runif(400))
d$y <- f$rd(cbind(exp(2*d$x - d$x^2), .15, 1.4))
d$o <- .04 * sin(5*d$x)
d$w <- rep(c(.5, 1.5, 1, 1), 100)
form <- list(y ~ s(x, k = 6) + offset(o), ~ 1, ~ 1)
ns <- asNamespace("mgcv")
routes <- character()
trace("newton", tracer = quote(routes <<- c(routes, paste("Newton", scoreType))),
      where = ns, print = FALSE)
trace("bfgs", tracer = quote(routes <<- c(routes, paste("BFGS", scoreType,
      if (isTRUE(family$qapprox)) "quadratic" else "exact"))), where = ns, print = FALSE)
trace("efsud", tracer = quote(routes <<- c(routes, "EFS")), where = ns, print = FALSE)

check_fit <- function(b, data = d) {
  stopifnot(length(b$sp) > 0L, all(is.finite(b$sp)), all(b$sp > 0),
            all(is.finite(coef(b))), is.finite(b$gcv.ubre),
            identical(colnames(fitted(b)), c("scale", "shape", "kappa")))
  pr <- predict(b, type = "response", se.fit = TRUE)
  lp <- predict(b, type = "link", se.fit = TRUE)
  near(pr$fit, fitted(b))
  for (j in 1:3) {
    near(pr$fit[, j], f$linfo[[j]]$linkinv(lp$fit[, j]))
    near(pr$se.fit[, j], f$linfo[[j]]$mu.eta(lp$fit[, j]) * lp$se.fit[, j])
  }
  near(logLik(b), f$ll(data$y, NULL, NULL, data$w, f, eta = lp$fit)$l)
}
reml <- list()
for (opt in c("newton", "bfgs", "efs")) {
  routes <- character()
  b <- gam(form, family = f, data = d, weights = w, method = "REML",
           optimizer = if (opt == "efs") opt else c("outer", opt))
  stopifnot(identical(routes, switch(opt, newton = "Newton REML",
                                      bfgs = "BFGS REML exact", efs = "EFS")),
            identical(b$method, "REML"), identical(b$outer.info$conv, "full convergence"))
  check_fit(b)
  reml[[opt]] <- b
}
for (opt in c("bfgs", "efs"))
  near(reml[[opt]]$gcv.ubre, reml$newton$gcv.ubre, tolerance = 1e-5)

# Grouped deletions with distinct, reordered prediction rows and gamma > 1.
nei <- list(a = c(4, 5, 6, 10, 11, 12, 19, 20, 21), ma = c(3, 6, 9),
            d = c(5, 12, 11, 20), md = c(1, 3, 4), jackknife = FALSE)
scores <- numeric()
for (method in c("NCV", "QNCV")) {
  routes <- character()
  optimizer <- if (method == "NCV") c("outer", "bfgs") else "efs"
  b <- gam(form, family = f, data = d, weights = w, method = method,
           optimizer = optimizer, nei = nei, gamma = 1.2)
  stopifnot(identical(b$method, method), identical(routes,
    if (method == "NCV") "BFGS NCV exact" else "BFGS NCV quadratic"))
  check_fit(b)
  scores[method] <- b$gcv.ubre
  ix <- nei$d
  eta <- predict(b, type = "link")
  cv <- attr(b$gcv.ubre, "eta.cv")
  full <- f$ll(d$y[ix], NULL, NULL, d$w[ix], f, eta = eta[ix, ], deriv = 1, ncv = TRUE)
  if (method == "NCV") {
    rows <- f$ll(d$y[ix], NULL, NULL, d$w[ix], f, eta = cv)$l0
    score <- -sum(full$l0 + 1.2 * (rows - full$l0))
  } else {
    delta <- cv - eta[ix, ]
    correction <- vapply(seq_along(ix), function(i) {
      hessian <- matrix(full$l2[i, f$tri$i2], 3, 3)
      sum(full$l1[i, ] * delta[i, ]) + .5 * sum(delta[i, ] * (hessian %*% delta[i, ]))
    }, numeric(1))
    score <- -sum(full$l0 + 1.2 * correction)
  }
  near(b$gcv.ubre, score)
}
stopifnot(abs(diff(scores)) > 1e-7)
untrace("newton", where = ns)
untrace("bfgs", where = ns)
untrace("efsud", where = ns)

# Default leave-one-out also exercises the jackknife covariance path.
b <- gam(form, family = f, data = d, weights = w, method = "QNCV",
         optimizer = c("outer", "bfgs"))
stopifnot(identical(b$method, "QNCV"), is.finite(b$NCV), all(is.finite(b$Vj)))

# scasm() was added in mgcv 1.9-5. Check the actual constrained fitter,
# reported Laplace criterion and monotonicity on the response scale.
if ("scasm" %in% getNamespaceExports("mgcv")) {
  routes <- character()
  trace("scasm.fit", tracer = quote(routes <<- c(routes, "scasm.fit")),
        where = ns, print = FALSE)
  set.seed(22)
  for (constrained in c(FALSE, TRUE)) {
    ff <- if (constrained)
      list(y ~ s(x, bs = "sc", xt = "m+", k = 6) + offset(o), ~ 1, ~ 1) else form
    # Add a local decrease so some monotonicity inequalities bind.
    dd <- if (constrained) transform(d, y = y * exp(.5 * sin(2 * pi * x))) else d
    b <- mgcv::scasm(ff, family = f, data = dd, weights = w)
    stopifnot(identical(b$method, "REML"), identical(b$optimizer, "EFS"),
              isTRUE(b$outer.info$converged), b$outer.info$niter > 1)
    near(b$gcv.ubre, -b$laml)
    check_fit(b, dd)
    if (constrained) {
      slack <- drop(b$Ain %*% coef(b)) - b$bin
      stopifnot(all(slack >= -1e-8), any(abs(slack) < 1e-8), any(slack > 1e-8))
      nd <- data.frame(x = seq(0, 1, length.out = 101), o = 0)
      stopifnot(all(diff(predict(b, nd, type = "response")[, 1L]) >= -1e-8))
    }
  }
  stopifnot(identical(routes, rep("scasm.fit", 2L)))
  untrace("scasm.fit", where = ns)
}

# gratia can use the native distribution directly on a fitted model.
if (requireNamespace("gratia", quietly = TRUE) &&
    "quantile_residuals" %in% getNamespaceExports("gratia")) {
  b <- reml$newton
  near(gratia::quantile_residuals(b, type = "pit"), f$cdf(d$y, fitted(b)))
}
