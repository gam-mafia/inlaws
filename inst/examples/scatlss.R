library("inlaws")

# Symmetric Student-t data with covariate-dependent tail weight.
set.seed(46)
n <- 900L
dat <- data.frame(x = runif(n), z = runif(n))
truth <- cbind(mu = sin(6 * dat$x), sigma = exp(-.5 + .2 * dat$z),
               nu = 3 + exp(-.3 + 1.5 * dat$x))
dat$y <- scatlss()$rd(truth)
form <- list(y ~ s(x, k = 6), ~ z, ~ x)
fit <- gam(form, data = dat, family = scatlss(), method = "REML")
pars <- predict(fit, type = "response", se.fit = TRUE)
print(head(pars$fit))
# Third formula ~ 1 estimates global df, with varying location and scale.
global <- gam(list(y ~ s(x, k = 6), ~ z, ~ 1), data = dat,
              family = scatlss(), method = "REML")
print(range(fitted(global)[, "nu"]))
# Response standard deviations include the df adjustment.
print(head(pars$fit[, 2] * sqrt(pars$fit[, 3] / (pars$fit[, 3] - 2))))
