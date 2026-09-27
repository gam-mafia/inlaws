library(inlaws)

near <- function(a, b, tol = 1e-6) {
  stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
}

# Keep all four observation types in a reproducible, small NCV fixture.
set.seed(43)
n <- 160L
dat <- data.frame(x = runif(n))
y <- rlnorm(n, sin(6 * dat$x), 0.5)
dat$y <- cbind(y, y)
kind <- seq_len(n) %% 4L
left <- kind == 1L & y < 1
right <- kind == 2L & y > 1
interval <- kind == 3L
dat$y[left, ] <- cbind(rep(1, sum(left)), rep(-Inf, sum(left)))
dat$y[right, ] <- cbind(rep(1, sum(right)), rep(Inf, sum(right)))
dat$y[interval, ] <- cbind(y[interval] * 0.9, y[interval] * 1.1)
stopifnot(any(left), any(right), any(interval), any(kind == 0L))

# Each neighbourhood drops and predicts four observations together.
grouped <- list(a = seq_len(n), ma = seq.int(4L, n, 4L),
                d = seq_len(n), md = seq.int(4L, n, 4L))
loo <- list(a = seq_len(n), ma = seq_len(n),
            d = seq_len(n), md = seq_len(n))

check_ncv <- function(fit, fixed = FALSE) {
  prediction <- predict(fit, type = "response")
  stopifnot(
    identical(fit$method, "NCV"), isTRUE(fit$converged),
    nobs(fit) == n,
    all(is.finite(coef(fit))),
    length(fit$sp) == 1L, all(is.finite(fit$sp)), all(fit$sp > 0),
    all(is.finite(fit$gcv.ubre)),
    length(prediction) == n, all(is.finite(prediction)), all(prediction > 0),
    is.finite(fit$family$getTheta(TRUE)), fit$family$getTheta(TRUE) > 0,
    fit$family$n.theta == if (fixed) 0L else 1L
  )
  near(prediction, fitted(fit))
  near(prediction, exp(predict(fit, type = "link")))
  if (fixed) near(fit$family$getTheta(TRUE), 0.5)
}

for (fixed in c(FALSE, TRUE)) {
  theta <- if (fixed) 0.5 else NULL
  # NULL nei requests leave-one-out NCV, including estimation of sigma when
  # theta is NULL. Grouped NCV exercises explicit censor-preserving subsets.
  for (neighbourhood in list(NULL, grouped)) {
    fit <- gam(y ~ s(x, k = 7), data = dat, family = clognormal(theta),
               method = "NCV", nei = neighbourhood)
    check_ncv(fit, fixed)
    if (is.null(neighbourhood) && !fixed) implicit_loo <- fit
  }

  # bam NCV optimizes the working-model criterion: do not require it to
  # reproduce gam's NCV optimum. Verify the discrete fitting route itself.
  fit <- bam(y ~ s(x, k = 7), data = dat, family = clognormal(theta),
             method = "NCV", discrete = TRUE, nei = grouped, nthreads = 1)
  check_ncv(fit, fixed)
}

# Explicit singleton neighbourhoods must reproduce default leave-one-out NCV.
explicit_loo <- gam(y ~ s(x, k = 7), data = dat, family = clognormal(),
                    method = "NCV", nei = loo)
check_ncv(explicit_loo)
near(coef(implicit_loo), coef(explicit_loo))
near(implicit_loo$sp, explicit_loo$sp)
near(implicit_loo$family$getTheta(), explicit_loo$family$getTheta())
near(implicit_loo$gcv.ubre, explicit_loo$gcv.ubre)
