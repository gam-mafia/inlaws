# Separate predictors for endpoint membership and interior beta responses.
# Run after installing inlaws: source(system.file("examples", "zoibetals.R",
#                                            package = "inlaws"))
library("inlaws")
set.seed(217)
n <- 600
d <- data.frame(x = runif(n, -1, 1), z = runif(n, -1, 1))
f <- zoibetals("both")
# x changes the interior mean; z changes endpoint membership and precision.
eta <- cbind(.3 + 1.4*d$x, 2.5 + .5*d$z, -.8 + 1.8*d$z, -1.1 - d$z)
d$y <- f$rd(eta)
fit <- mgcv::gam(list(y ~ s(x, k = 5), ~ z, ~ z, ~ z),
                 data = d, family = f, method = "REML")
ordered <- mgcv::gam(y ~ s(x, k = 5) + s(z, k = 5),
                     data = d, family = ordbeta(), method = "REML")
# Ordered beta couples both endpoint probabilities and mu to one predictor.
# Inflated beta can change membership with z while retaining mu as a function of x.
grid <- expand.grid(x = seq(-1, 1, length.out = 25),
                    z = seq(-1, 1, length.out = 25))
truth_mu <- plogis(.3 + 1.4*grid$x)
truth_p <- exp(cbind(-.8 + 1.8*grid$z, 0, -1.1-grid$z))
truth_p <- truth_p/rowSums(truth_p)
p <- predict_zoibetals(fit, grid, type = "all")
op <- predict_ordbeta(ordered, grid, type = "probabilities")
om <- predict_ordbeta(ordered, grid, type = "conditional")
rmse <- function(a, b) sqrt(mean((a-b)^2))
print(rbind(inflated = c(interior_mean = rmse(p[, "mu"], truth_mu),
                         component_probabilities = rmse(p[, c("p0", "pI", "p1")], truth_p)),
            ordered = c(interior_mean = rmse(om, truth_mu),
                        component_probabilities = rmse(op, truth_p))))
predict_zoibetals(fit, head(grid), se.fit = TRUE) # complete mean with SE
head(predict_zoibetals(fit, grid, type = "variance")) # response variance
# Likelihood weights are not passed as beta precision or sampling multipliers.
simulated_response <- fit$family$rd(fitted(fit))
diagnostics <- fit$family$gratia_diagnostics(fitted(fit), fit$prior.weights,
                                            1, fit$y, fit$linear.predictors)
pit <- runif(n, diagnostics$cdf_left(d$y), diagnostics$cdf(d$y))
stopifnot(all(pit >= 0 & pit <= 1), length(simulated_response) == n)
# Other configurations use three formulas: mu, phi, then the active endpoint.
# zoibetals("zero") fixes p1=0; zoibetals("one") fixes p0=0.
# Active but unobserved endpoints cause an error; choose the model explicitly.
