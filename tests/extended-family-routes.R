library("inlaws")

near <- function(x, y, tol = 1e-6) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}
parameters <- c(cgamma = .4, clognormal = .5, cmp = 1.5, gpoisson = 2)
fixtures <- list()
for (name in names(parameters)) {
  set.seed(1843)
  n <- 96L
  d <- data.frame(x = runif(n), o = runif(n, -.2, .2),
    w = rep(c(.5, 1, 2, 1), length.out = n))
  d$y <- get(name)(parameters[[name]])$rd(exp(.5 + sin(6*d$x) + d$o), rep(1, n), 1)
  if (name %in% c("cgamma", "clognormal")) {
    y <- d$y; d$y <- cbind(y, y); kind <- seq_len(n) %% 4L
    left <- kind == 1L & y < 1; right <- kind == 2L & y > 2; interval <- kind == 3L
    d$y[left, ] <- cbind(rep(1, sum(left)), rep(-Inf, sum(left)))
    d$y[right, ] <- cbind(rep(2, sum(right)), rep(Inf, sum(right)))
    d$y[interval, ] <- cbind(.9*y[interval], 1.1*y[interval])
    stopifnot(any(left), any(right), any(interval), any(kind == 0L))
  }
  fixtures[[name]] <- d
}
check <- function(f, d, method, fixed, theta) {
  stopifnot(identical(f$method, method), isTRUE(f$converged),
    all(is.finite(coef(f))), all(is.finite(f$Vp)),
    all(is.finite(f$sp) & f$sp > 0), is.finite(f$gcv.ubre))
  if (!is.null(f$outer.info$conv)) stopifnot(f$outer.info$conv == "full convergence")
  pr <- predict(f, type = "response", se.fit = TRUE)
  near(pr$fit, fitted(f)); near(pr$fit, exp(predict(f, type = "link")))
  stopifnot(all(is.finite(pr$se.fit)), all(pr$se.fit >= 0))
  if (fixed) near(f$family$getTheta(TRUE), theta)
  near(f$deviance, sum(f$family$dev.resids(d$y, fitted(f), d$w)))
  ll <- -f$family$aic(d$y, fitted(f), wt = d$w, dev = 0)/2
  near(as.numeric(logLik(f)), ll)
  near(AIC(f), -2*ll + 2*attr(logLik(f), "df"))
}

# Validate the available routes; GAM NCV is explicitly blocked under mgcv 1.9-4.
# Convergence warnings are failures even though bam may set converged=TRUE.
options(warn = 2)
for (name in names(parameters)) {
  d <- fixtures[[name]]; constructor <- get(name); theta <- parameters[[name]]
  for (fixed in c(FALSE, TRUE)) {
    for (method in c("REML", "ML", "NCV")) {
      if (method == "NCV" && packageVersion("mgcv") == "1.9.4") next
      message(name, " gam ", method, " fixed=", fixed)
      fit <- gam(y ~ s(x, k = 5) + offset(o), data = d, weights = w,
        family = constructor(if (fixed) theta else NULL), method = method)
      check(fit, d, method, fixed, theta)
      if (method == "NCV") {
        singleton <- list(a = seq_len(nrow(d)), ma = seq_len(nrow(d)),
          d = seq_len(nrow(d)), md = seq_len(nrow(d)))
        explicit <- gam(y ~ s(x, k = 5) + offset(o), data = d, weights = w,
          family = constructor(if (fixed) theta else NULL), method = method, nei = singleton)
        near(coef(explicit), coef(fit)); near(explicit$sp, fit$sp)
        near(explicit$family$getTheta(TRUE), fit$family$getTheta(TRUE))
      }
    }
    for (method in c("REML", "ML", "fREML")) {
      message(name, " bam ", method, " fixed=", fixed)
      fit <- bam(y ~ s(x, k = 5) + offset(o), data = d, weights = w,
        family = constructor(if (fixed) theta else NULL), method = method,
        scale = 1, chunk.size = 31, nthreads = 1)
      check(fit, d, method, fixed, theta)
    }
    for (method in c("fREML", "NCV")) {
      if (method == "NCV" && packageVersion("mgcv") < "1.9.4") next
      message(name, " discrete bam ", method, " fixed=", fixed)
      fit <- bam(y ~ s(x, k = 5) + offset(o), data = d, weights = w,
        family = constructor(if (fixed) theta else NULL), method = method,
        discrete = TRUE, scale = 1, nthreads = 1)
      check(fit, d, method, fixed, theta)
    }
  }
  # A second fit cannot mutate the constructor or the first fitted family.
  for (discrete in c(FALSE, TRUE)) {
    template <- constructor()
    initial <- template$getTheta(TRUE)
    first <- bam(y ~ s(x, k = 5) + offset(o), data = d, weights = w,
      family = template, method = "fREML", discrete = discrete, scale = 1, nthreads = 1)
    saved <- first$family$getTheta(TRUE)
    second <- bam(y ~ s(x, k = 4), data = d, weights = w,
      family = template, method = "fREML", discrete = discrete, scale = 1, nthreads = 1)
    near(template$getTheta(TRUE), initial); near(first$family$getTheta(TRUE), saved)
    check(first, d, "fREML", FALSE, theta)
  }
}

# Independent Poisson identities, including actual binning and block processing.
# Use zero offsets for the mgcv NCV reference, whose native route omits them.
set.seed(371)
d <- data.frame(x = runif(90), w = rep(c(.5, 1, 2), 30))
d$y <- rpois(90, exp(.5 + sin(6*d$x)))
for (constructor in list(cmp, gpoisson)) {
  for (route in c("gam.ML", "gam.NCV", "bam.REML", "bam.ML", "bam.fREML", "discrete.fREML", "discrete.NCV")) {
    part <- strsplit(route, ".", fixed = TRUE)[[1]]
    if (route == "gam.NCV" && packageVersion("mgcv") == "1.9.4") next
    if (part[2] == "NCV" && part[1] == "discrete" && packageVersion("mgcv") < "1.9.4") next
    args <- list(formula = y ~ s(x, k = 5), data = d, weights = d$w, method = part[2])
    if (part[1] != "gam") args <- c(args, list(scale = 1, nthreads = 1,
      chunk.size = 31, discrete = if (part[1] == "discrete") 30 else FALSE))
    fun <- if (part[1] == "gam") gam else bam
    a <- do.call(fun, c(args, list(family = constructor(1))))
    b <- do.call(fun, c(args, list(family = poisson())))
    near(fitted(a), fitted(b), 5e-5)
    near(as.numeric(logLik(a)), as.numeric(logLik(b)), 1e-5)
  }
}

# Uncensored censored-family limits reproduce independent standard families.
d$o <- runif(nrow(d), -.2, .2)
mu <- exp(.5 + sin(6*d$x) + d$o)
for (name in c("cgamma", "clognormal")) {
  if (name == "cgamma") d$y <- rgamma(nrow(d), 1/.4, scale = mu*.4) else
    d$y <- rlnorm(nrow(d), log(mu) - .5^2/2, .5)
  d$logged <- log(d$y)
  for (method in c("REML", "ML")) {
    candidate <- bam(y ~ s(x, k = 5) + offset(o), data = d, weights = w,
      family = get(name)(parameters[[name]]), method = method, scale = 1, chunk.size = 31)
    if (name == "cgamma") {
      reference <- bam(y ~ s(x, k = 5) + offset(o), data = d, weights = w,
        family = Gamma("log"), method = method, scale = .4, chunk.size = 31)
      expected <- fitted(reference)
    } else {
      reference <- bam(logged ~ s(x, k = 5) + offset(o), data = d, weights = w,
        family = gaussian(), method = method, scale = .5^2, chunk.size = 31)
      expected <- exp(fitted(reference) + .5^2/2)
    }
    near(fitted(candidate), expected, 2e-5)
  }
}
