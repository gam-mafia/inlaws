library("inlaws")

near <- function(x, y, tol) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
    max(abs(x - y) / pmax(1, abs(x), abs(y))) < tol)
}

set.seed(1146)
n <- 240L
d <- data.frame(x = runif(n), z = runif(n), o = runif(n, -.2, .2),
                w = rep(c(.5, 2, 1, 1.5), length.out = n))
d$w[5L] <- 0
eta <- cbind(.9 * sin(2 * pi * d$x), .7 * cos(2 * pi * d$x),
             1.8 + .8 * sin(2 * pi * d$z))
p <- exp(cbind(0, eta[, 1:2]))
p <- p / rowSums(p)
total <- sample(18:35, n, replace = TRUE)
d$y <- t(vapply(seq_len(n), function(i) {
  alpha <- p[i, ] * exp(eta[i, 3L])
  drop(rmultinom(1, total[i], rgamma(3, shape = alpha)))
}, integer(3)))
form <- list(y ~ s(x, k = 5) + offset(o), ~s(x, k = 5) + offset(-o),
             ~s(z, k = 5) + offset(o / 2))

# Check the effective optimizer and criterion with free smoothing parameters.
ns <- asNamespace("mgcv")
routes <- character()
fourth_order_allocations <- 0L
trace("trind.generator", tracer = quote(fourth_order_allocations <<-
        fourth_order_allocations + 1L), where = ns, print = FALSE)
fam <- dirmult(2)
stopifnot(fourth_order_allocations == 0L)
trace("newton", tracer = quote(routes <<- c(routes, paste("newton", scoreType))),
      where = ns, print = FALSE)
trace("bfgs", tracer = quote(routes <<- c(routes, paste("bfgs", scoreType,
      isTRUE(family$qapprox)))), where = ns, print = FALSE)
trace("efsud", tracer = quote(routes <<- c(routes, "efs")), where = ns, print = FALSE)
newton <- gam(form, family = fam, data = d, weights = w, method = "REML")
bfgs <- gam(form, family = dirmult(2), data = d, weights = w, method = "REML",
            optimizer = c("outer", "bfgs"))
efs <- gam(form, family = dirmult(2), data = d, weights = w, method = "REML",
           optimizer = "efs")
stopifnot(identical(routes, c("newton REML", "bfgs REML FALSE", "efs")))
for (fit in list(newton, bfgs, efs)) {
  stopifnot(fit$method == "REML", length(fit$sp) == 3L,
            fit$outer.info$conv == "full convergence", all(is.finite(fit$sp)),
            all(is.finite(vcov(fit))))
}
near(newton$gcv.ubre, bfgs$gcv.ubre, 1e-7)
near(coef(newton), coef(bfgs), 5e-4)
near(newton$gcv.ubre, efs$gcv.ubre, 1e-5)

routes <- character()
for (method in c("NCV", "QNCV")) {
  fit <- gam(form, family = dirmult(2), data = d, weights = w, method = method)
  stopifnot(fit$method == method, length(fit$sp) == 3L,
            identical(fit$family$qapprox, method == "QNCV"),
            fit$outer.info$conv == "full convergence",
            is.finite(fit$gcv.ubre), is.finite(fit$NCV), all(is.finite(vcov(fit))))
}
stopifnot(identical(routes, c("bfgs NCV FALSE", "bfgs NCV TRUE")))
stopifnot(fourth_order_allocations == 1L)
untrace("trind.generator", where = ns)
untrace("newton", where = ns)
untrace("bfgs", where = ns)
untrace("efsud", where = ns)

# With no penalties, older mgcv versions can score NCV/QNCV but fail when
# computing default jackknife covariance through the empty smoothing map.
unpenalized <- list(y ~ 1 + offset(o), ~1, ~1)
for (method in c("NCV", "QNCV")) {
  fit <- gam(unpenalized, family = fam, data = d, weights = w, method = method,
             nei = list(jackknife = FALSE))
  stopifnot(fit$method == method, length(fit$sp) == 0L, is.finite(fit$gcv.ubre))
  default <- tryCatch(gam(unpenalized, family = fam, data = d, weights = w,
                         method = method), error = identity)
  if (packageVersion("mgcv") < "1.9.5") {
    stopifnot(inherits(default, "error"),
              grepl("supply nei with jackknife = FALSE", conditionMessage(default),
                    fixed = TRUE))
  } else {
    stopifnot(default$method == method, is.finite(default$gcv.ubre),
              is.finite(default$NCV), all(is.finite(vcov(default))))
  }
}

# Grouped deletions with differently sized, reordered and repeated prediction
# sets. Include the zero case weight and a non-unit gamma.
nei <- list(a = c(4, 5, 6, 11, 12, 13, 4, 5, 6, 19, 20, 21),
            ma = c(3, 6, 9, 12), d = c(5, 6, 13, 12, 6, 20),
            md = c(2, 4, 5, 6), jackknife = FALSE)

# Compare outer derivatives against finite differences with a singleton
# coefficient block and a coefficient block shared across two predictors.
shared <- list(y ~ offset(o), ~1 + offset(-o), ~1 + offset(o / 2),
               1 + 2 ~ s(x, k = 5) - 1)
G <- gam(shared, family = fam, data = d, weights = w, fit = FALSE)
nsp <- length(G$sp)
stopifnot(nsp == 1L)
Sl <- mgcv:::Sl.setup(G)
X <- mgcv:::Sl.initial.repara(Sl, G$X, both.sides = FALSE)
control <- gam.control(epsilon = 1e-10)
rho <- log(2)
start <- NULL
fit_at <- function(rho, score = "REML", deriv = 1L, qapprox = FALSE) {
  family <- fam
  family$qapprox <- qapprox
  mgcv:::gam.fit5(X, G$y, rho, Sl, weights = G$w, offset = G$offset,
    deriv = deriv, family = family, scoreType = score, control = control,
    start = start, gamma = 1.2, nei = nei)
}
base <- fit_at(rho, deriv = 2L)
start <- base$coefficients
h <- .01
stencil <- c(1, -8, 8, -1) / (12 * h)
reml <- lapply(rho + h * c(-2, -1, 1, 2), fit_at)
near(base$REML1, sum(stencil * vapply(reml, function(f) f$REML, numeric(1))), 5e-5)
near(base$REML2, sum(stencil * vapply(reml, function(f) f$REML1, numeric(1))), 1e-4)

for (quad in c(FALSE, TRUE)) {
  base <- fit_at(rho, "NCV", qapprox = quad)
  pp <- fit_at(rho + h, "NCV", qapprox = quad)
  mm <- fit_at(rho - h, "NCV", qapprox = quad)
  near(base$NCV1, (as.numeric(pp$NCV) - as.numeric(mm$NCV)) / (2 * h), 1e-4)
  if (!quad) {
    # Independently score the held-out count vectors, including normalization,
    # case weights and the weighted full-data correction for gamma.
    rows <- nei$d
    held_out <- attr(base$NCV, "eta.cv")
    logpmf <- function(eta) {
      y <- d$y[rows, , drop = FALSE]
      ex <- exp(cbind(0, eta[, 1:2, drop = FALSE]))
      phi <- exp(eta[, 3L])
      alpha <- phi * ex / rowSums(ex)
      N <- rowSums(y)
      lgamma(N + 1) - rowSums(lgamma(y + 1)) +
        lgamma(phi) - lgamma(phi + N) + rowSums(lgamma(alpha + y) - lgamma(alpha))
    }
    v <- d$w[rows] * (1.2 * logpmf(held_out) -
      .2 * logpmf(base$linear.predictors[rows, , drop = FALSE]))
    near(base$NCV, -sum(v), 1e-10)
  }
}
