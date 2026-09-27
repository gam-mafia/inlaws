library(inlaws)
near <- function(x, y, tol = 1e-6) stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)), max(abs(x - y) / pmax(1, abs(x), abs(y))) < tol)
error <- function(expr, pattern) {
  e <- tryCatch(expr, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
jet <- getFromNamespace(".zinb_logjet", "inlaws")
moments <- getFromNamespace(".zinb_moments", "inlaws")
set.seed(78)
for (h in c(FALSE, TRUE)) {
  fam <- if (h) zanb() else zinb()
  for (par in list(c(.01, .3, .1), c(3, .7, 2), c(20, .9, 10))) {
    yy <- 0:2000
    mu <- matrix(rep(par, each = length(yy)), ncol = 3)
    eta <- cbind(log(mu[, 1]), qlogis(mu[, 2]), log(mu[, 3]))
    pmf <- exp(jet(yy, eta, h)$l0)
    q0 <- dnbinom(0, mu = par[1], size = par[3])
    ref <- par[2] * dnbinom(yy, mu = par[1], size = par[3]) / (if (h) 1 - q0 else 1)
    ref[1] <- 1 - par[2] + if (h) 0 else par[2] * q0
    near(pmf, ref, 1e-10)
    near(sum(pmf), 1)
    mo <- moments(matrix(par, 1), h)
    near(sum(yy * pmf), mo$mean)
    near(sum(yy^2 * pmf) - mo$mean^2, mo$variance)
    near(fam$cdf(yy, mu), cumsum(pmf))
    near(exp(fam$cdf(yy, mu, logp = TRUE)), cumsum(pmf))
    pr <- seq(.001, .999, length.out = 99)
    pars <- matrix(rep(par, each = length(pr)), ncol = 3)
    quant <- fam$qf(pr, pars)
    stopifnot(all(fam$cdf(quant, pars) >= pr - 1e-12), all(fam$cdf(quant - 1, pars) < pr))
    sims <- fam$rd(matrix(rep(par, each = 20000), ncol = 3))
    stopifnot(abs(mean(sims == 0) - ref[1]) < .015)
  }
  tiny <- cbind(mu = rep(1e-16, 1000), p = rep(1, 1000), theta = rep(2, 1000))
  if (h) stopifnot(all(fam$rd(tiny) == 1))
  edge <- cbind(mu = c(1, 1, 1), p = c(0, 1, .5), theta = c(2, 2, 2))
  stopifnot(identical(fam$qf(c(1, 0, 1), edge), c(0, 0, Inf)))
  stopifnot(is.nan(fam$qf(NA_real_, matrix(c(1, .5, 2), 1))))
  n <- 650
  d <- data.frame(x = runif(n), z = runif(n), a = runif(n, -.15, .15), b = runif(n, -.1, .1), c = runif(n, -.1, .1))
  pars <- cbind(exp(1 + .6 * sin(6 * d$x) + d$a), plogis(.5 + .3 * d$x + d$b), exp(.6 + .5 * d$z + d$c))
  d$y <- fam$rd(pars)
  form <- list(y ~ s(x, k = 6) + offset(a), ~ x + offset(b), ~ z + offset(c))
  fits <- lapply(list(c("outer", "newton"), c("outer", "bfgs"), "efs"), function(opt) {
    gam(form, data = d, family = fam, method = "REML", optimizer = opt, control = gam.control(efs.tol = 1e-6))
  })
  f <- fits[[1]]
  for (g in fits) {
    stopifnot(g$outer.info$conv == "full convergence", all(is.finite(g$Vp)))
    near(g$gcv.ubre, f$gcv.ubre, 1e-4)
    near(g$fitted.values, f$fitted.values, .02)
  }
  for (method in c("NCV", "QNCV")) {
    g <- gam(form, data = d, family = fam, method = method)
    stopifnot(g$outer.info$conv == "full convergence", all(is.finite(g$gcv.ubre)))
  }
  near(as.numeric(logLik(f)), sum(jet(d$y, f$linear.predictors, h)$l0))
  near(AIC(f), -2 * as.numeric(logLik(f)) + 2 * attr(logLik(f), "df"))
  near(f$deviance, sum(residuals(f, "deviance")^2))
  mo <- moments(f$fitted.values, h)
  near(residuals(f, "response"), d$y - mo$mean)
  near(residuals(f, "pearson"), (d$y - mo$mean) / sqrt(mo$variance))
  stopifnot(identical(colnames(f$fitted.values), c("mu", "p", "theta")))
  pr <- predict(f, newdata = d[1:7, ], type = "response", se.fit = TRUE)
  lp <- predict(f, newdata = d[1:7, ], type = "link", se.fit = TRUE)
  near(pr$fit, cbind(exp(lp$fit[, 1]), plogis(lp$fit[, 2]), exp(lp$fit[, 3])))
  near(pr$se.fit, lp$se.fit * cbind(pr$fit[, 1], pr$fit[, 2] * (1 - pr$fit[, 2]), pr$fit[, 3]))
  # Draws and quantiles use the fitted family's predicted parameters.
  set.seed(1)
  a <- fam$rd(pr$fit)
  set.seed(1)
  near(a, fam$qf(runif(nrow(pr$fit)), pr$fit))
  stopifnot(length(a) == 7L, all(a >= 0), all(a == floor(a)),
            length(fam$rd(pr$fit[1, , drop = FALSE])) == 1L)
  for (p in c(.2, .8)) {
    q <- fam$qf(p, pr$fit)
    stopifnot(all(fam$cdf(q, pr$fit) >= p - 1e-12),
              all(fam$cdf(q - 1, pr$fit) < p + 1e-12))
  }
  # Independent constant-model MLE and hurdle component separation.
  constant <- gam(list(y ~ 1, ~1, ~1), data = d, family = fam, method = "REML")
  fn <- function(b) {
    mu <- exp(b[1])
    p <- plogis(b[2])
    t <- exp(b[3])
    q0 <- dnbinom(0, mu = mu, size = t)
    -sum(ifelse(d$y == 0, log(1 - p + if (h) 0 else p * q0), log(p) + dnbinom(d$y, mu = mu, size = t, log = TRUE) - if (h) log1p(-q0) else 0))
  }
  opt <- optim(coef(constant) + .05, fn, method = "BFGS", control = list(reltol = 1e-12))
  near(coef(constant), opt$par, 5e-4)
  if (h) {
    near(plogis(coef(constant)[2]), mean(d$y > 0))
    yp <- d$y[d$y > 0]
    pos <- optim(coef(constant)[c(1, 3)] + .05, function(b) {
      -sum(dnbinom(yp, mu = exp(b[1]), size = exp(b[2]), log = TRUE) - log1p(-dnbinom(0, mu = exp(b[1]), size = exp(b[2]))))
    },
    method = "BFGS", control = list(reltol = 1e-12)
    )
    near(coef(constant)[c(1, 3)], pos$par, 5e-4)
  }
  # Integer weights equal replication in a fixed parametric basis.
  d$w <- sample(0:3, n, replace = TRUE)
  ff <- list(y ~ x + offset(a), ~ x + offset(b), ~ z + offset(c))
  fw <- gam(ff, data = d, weights = w, family = fam)
  fr <- gam(ff, data = d[rep(seq_len(n), d$w), ], family = fam)
  near(coef(fw), coef(fr), 1e-4)
  d$xcopy <- d$x
  rankdef <- gam(list(y ~ x + xcopy + offset(a), ~ x + offset(b), ~ z + offset(c)), data = d, weights = w, family = fam)
  near(fw$fitted.values, rankdef$fitted.values, 1e-4)
  # mgcv 1.9-3 has an upstream fixed-sp setup bug when later predictors
  # have no smooths. Reusing the unfixed setup avoids that setup path.
  if (utils::packageVersion("mgcv") < "1.9.4") {
    G <- gam(form, data=d, family=fam, fit=FALSE)
    fs <- gam(G=G, start=coef(f), sp=f$sp)
  } else fs <- gam(form, data=d, family=fam, start=coef(f), sp=f$sp)

  near(fs$fitted.values, f$fitted.values, 1e-4)
  # Meaningful validation failures and boundary warnings.
  error(gam(list(y ~ 1, ~1, ~1), data = data.frame(y = c(0, .5)), family = fam), "integer counts")
  error(gam(list(y ~ 1, ~1, ~1), data = data.frame(y = c(0, 0)), family = fam), "all-zero")
  error(gam(y ~ x, data = d, family = fam), "three linear predictors")
  error(gam(list(y ~ 1, ~1, ~1), data = d, weights = rep(0, n), family = fam), "positive-weight")
  X <- matrix(1, 2, 3)
  attr(X, "lpi") <- list(1, 2, 3)
  warning_text <- character()
  withCallingHandlers(fam$init(c(1, 1), X, matrix(0, 0, 3), c(1, 1), NULL, c(0, 0, 0), fam),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(length(warning_text) == 2)
}
error(zinb(list("identity", "logit", "log")), "links")
error(zanb(list("log", "log")), "links")

# Poisson limits match the appropriate mixture and hurdle Poisson likelihoods.
zipll <- getFromNamespace("zipll", "mgcv")
y <- c(0, 1, 3, 8)
u <- c(-2, 0, 1, 2)
b <- c(-1, 0, 1, 2)
p <- plogis(b)
e <- cbind(u, b, rep(log(1e12), 4))
ref <- ifelse(y == 0, log(1 - p + p * exp(-exp(u))), log(p) + dpois(y, exp(u), log = TRUE))
near(jet(y, e, FALSE)$l0, ref, 1e-9)
near(jet(y, e, TRUE)$l0, zipll(y, u, log(-log1p(-p)))$l, 1e-9)
# Single-row parameter response predictions and offsets retain matrix shape.
stopifnot(identical(dim(predict(f, newdata = d[1, , drop = FALSE], type = "response")), c(1L, 3L)))
