# Optional developer validation; brms is suggested, not required at runtime.
# Run with an installed development version of inlaws and brms.
library(inlaws)
if (!requireNamespace("brms", quietly = TRUE)) stop("Install brms to run this comparison")

set.seed(71)
eta <- matrix(rnorm(60), 20, 3)
eta[, 3] <- eta[, 3] + 2
mu <- exp(cbind(0, eta[, 1:2]))
mu <- mu / rowSums(mu)
phi <- exp(eta[, 3])
y <- inlaws::dirichlet(2)$rd(eta, NULL, 1)
weights <- rep(c(0, .5, 1, 2), 5)
ll <- getFromNamespace(".dirichlet_derivatives", "inlaws")
jets <- getFromNamespace(".dirichlet_jets", "inlaws")(3)
ours <- ll(y, eta, weights, jets, 0)$l
reference <- sum(weights * brms::ddirichlet(y, alpha = mu * phi, log = TRUE))
stopifnot(isTRUE(all.equal(ours, reference, tolerance = 1e-12)))
fam <- inlaws::dirichlet(2)
stopifnot(isTRUE(all.equal(unname(fam$predict(fam, eta = eta)$fit),
                          unname(mu), tolerance = 1e-12)))

# Inspect the actual generated likelihood as well, without compiling or sampling.
# Explicitly model phi so that its log link is used, as in our final formula.
d <- data.frame(x = seq_len(nrow(y)))
colnames(y) <- c("a", "b", "c")
d$y <- y
code <- brms::make_stancode(brms::bf(y ~ x, phi ~ x), data = d,
                           family = brms::dirichlet(refcat = "a"))
stopifnot(grepl("softmax\\(mu\\) \\* phi", code),
          grepl("dirichlet_logit_lpdf", code))
cat("brms", as.character(packageVersion("brms")),
    ": likelihood, conditional means, and generated Stan parameterization agree\n")
cat("Weighted log likelihood:", ours, "\n")
