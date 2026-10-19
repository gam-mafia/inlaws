library("inlaws")

nanei <- getFromNamespace("nanei", "mgcv")
guard <- getFromNamespace(".inlaws_check_ncv_nei", "inlaws")
install_guard <- getFromNamespace(".inlaws_ncv_family", "inlaws")
version <- packageDescription("mgcv", fields = "Version")
near <- function(x, y, tol = 1e-7) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}
error_message <- function(expr) {
  err <- tryCatch(expr, error = identity)
  stopifnot(inherits(err, "error"))
  conditionMessage(err)
}
blocked <- function(expr) {
  msg <- error_message(expr)
  for (text in c(paste("installed mgcv", version), "neighbourhood preprocessing",
                 "response, predictor", "offset and weight", "resulting data's row indices"))
    stopifnot(grepl(text, msg, fixed = TRUE))
}

# Independent group-wise specification: filter each set, discard paired groups
# if either set is empty, then map surviving original indices to retained rows.
# No positional negative subscripts or cumulative-boundary subtraction.
expected <- function(nb, omit, n) {
  groups <- function(x, end) split(x, rep(seq_along(end), diff(c(0L, end))))
  a <- lapply(groups(nb$a, nb$ma), function(x) x[!x %in% omit])
  d <- lapply(groups(nb$d, nb$md), function(x) x[!x %in% omit])
  keep <- lengths(a) > 0L & lengths(d) > 0L
  a <- a[keep]; d <- d[keep]; rows <- setdiff(seq_len(n), omit)
  nb$a <- match(unlist(a, use.names = FALSE), rows)
  nb$d <- match(unlist(d, use.names = FALSE), rows)
  nb$ma <- cumsum(lengths(a)); nb$md <- cumsum(lengths(d))
  nb
}
same_nei <- function(x, y) {
  all(vapply(c("a", "ma", "d", "md"), function(nm)
    identical(as.numeric(x[[nm]]), as.numeric(y[[nm]])), logical(1)))
}
minimal <- list(a = 1:6, ma = c(3L, 6L), d = 1:6, md = c(3L, 6L))
want <- list(a = 1:5, ma = c(2L, 5L), d = 1:5, md = c(2L, 5L))
stopifnot(same_nei(expected(minimal, 2L, 6L), want))
actual <- nanei(minimal, 2L)
broken <- !same_nei(actual, want)
if (packageVersion("mgcv") == "1.9.4") {
  stopifnot(broken, identical(actual$a, 1:5), identical(actual$d, 1:5),
    identical(actual$ma, c(3L, 6L)), identical(actual$md, c(3L, 6L)))
}

# Includes zero matches in either set, partial groups, whole groups, cascading
# paired-group removal, and different ordered/repeated deletion/prediction sets.
cases <- list(
  both = minimal,
  deletion = list(a = 1:6, ma = c(3L, 6L), d = c(4L, 6L), md = 1:2),
  prediction = list(a = c(4L, 6L), ma = 1:2,
                    d = c(3L, 2L, 3L, 5L), md = c(3L, 4L)),
  neither = list(a = c(4L, 6L), ma = 1:2,
                 d = c(6L, 4L, 6L), md = c(1L, 3L)),
  whole = list(a = 1:6, ma = c(1L, 3L, 6L), d = 1:6, md = c(1L, 3L, 6L)),
  separate = list(a = c(1L, 2L, 3L, 1L, 2L, 3L, 8L, 10L), ma = c(3L, 6L, 8L),
                  d = c(3L, 3L, 2L, 9L, 8L, 9L, 10L), md = c(3L, 4L, 7L)),
  cascade = list(a = c(2L, 1L, 3L, 6L, 7L), ma = c(1L, 3L, 5L),
                 d = c(8L, 2L, 9L, 9L), md = c(1L, 2L, 4L)),
  mixed = list(a = 1:9, ma = c(1L, 4L, 9L),
               d = c(9L, 1L, 2L, 3L, 8L, 7L), md = c(1L, 4L, 6L)))
stopifnot(same_nei(expected(cases$deletion, 2L, 12L),
  list(a = 1:5, ma = c(2L, 5L), d = c(3L, 5L), md = 1:2)),
  same_nei(expected(cases$prediction, 2L, 12L),
  list(a = c(3L, 5L), ma = 1:2, d = c(2L, 2L, 4L), md = c(2L, 3L))),
  same_nei(expected(cases$neither, 2L, 12L),
  list(a = c(3L, 5L), ma = 1:2, d = c(5L, 3L, 5L), md = c(1L, 3L))))
for (nm in names(cases)) {
  omit <- if (nm == "whole") 1L else if (nm == "mixed") c(1L, 3L) else 2L
  nb <- cases[[nm]]; want <- expected(nb, omit, 12L); actual <- nanei(nb, omit)
  # Use gam itself, recording what the supported presetup hook receives.
  dat <- data.frame(y = seq(.1, .9, length.out = 12L), x = 1:12)
  dat$x[omit] <- NA
  fam <- betals(); seen <- NULL; omitted <- NULL
  fam$check.nei <- function(nei, missing) {
    seen <<- nei; omitted <<- missing
    guard(nei, missing)
  }
  result <- tryCatch(gam(list(y ~ x, ~ 1), data = dat, family = fam,
    method = "NCV", nei = nb, fit = FALSE), error = identity)
  stopifnot(same_nei(seen, actual), identical(as.integer(omitted), omit))
  if (!same_nei(actual, want)) {
    stopifnot(inherits(result, "error"))
    blocked(stop(result))
  } else {
    stopifnot(inherits(result, "gam.prefit"), same_nei(seen, want))
  }
}
if (packageVersion("mgcv") == "1.9.4") {
  # The no-match branches empty one or both otherwise valid index sets.
  # Keep these source-specific assertions separate from the partial-boundary
  # probe: a later mgcv may fix either defect without fixing the other.
  stopifnot(length(nanei(cases$deletion, 2L)$d) == 0L,
    length(nanei(cases$prediction, 2L)$a) == 0L,
    length(nanei(cases$neither, 2L)$a) == 0L,
    length(nanei(cases$neither, 2L)$d) == 0L)
}

# Every intended constructor installs the shared early hook. This checks public
# gam routes too, rather than equating a family field with integration coverage.
constructors <- list(betals, betabinomial, cgammals, clognormalls,
  log_logisticls, cmpls, pigls, gen_gammalss, cgen_gammalss, zinb, zanb)
for (constructor in constructors) {
  fam <- constructor()
  stopifnot(identical(fam$check.nei, guard),
    identical(install_guard(fam)$presetup, fam$presetup))
  dat <- data.frame(y = seq(.1, .9, length.out = 12L), x = 1:12)
  dat$x[2L] <- NA
  form <- c(list(y ~ x), rep(list(~ 1), fam$nlp - 1L))
  if (broken) for (method in c("NCV", "QNCV", "REML"))
    blocked(gam(form, data = dat, family = fam, method = method, nei = minimal))
}
# Preserve an existing supported presetup expression, once and in order.
fam <- betals(); fam$check.nei <- NULL
fam$presetup <- expression(presetup_ran <- TRUE)
fam <- install_guard(fam); presetup_ran <- FALSE
family <- fam; nei <- NULL; mf <- data.frame(y = 1:3)
eval(fam$presetup)
stopifnot(presetup_ran)

# Model-frame omission uses all formulae, offset arguments and case weights,
# after subset selection. Indices in nei refer to the selected pre-omission rows.
set.seed(441)
n <- 90L
dat <- data.frame(x = runif(n), z = runif(n), o = runif(n, -.1, .1),
                  o2 = runif(n, -.1, .1), w = 1)
p <- plogis(.3 + sin(6*dat$x))
dat$y <- rbeta(n, p*15, (1-p)*15)
form <- list(y ~ s(x, k = 4) + offset(o), ~ z + offset(o2))
grouped <- function(n) list(a = seq_len(n), ma = seq.int(3L, n, 3L),
  d = as.vector(apply(matrix(seq_len(n), 3L), 2L, rev)),
  md = seq.int(3L, n, 3L), jackknife = FALSE)
nb <- grouped(n)
if (broken) for (variable in c("y", "x", "z", "o", "o2", "w", "offset_argument")) {
  dd <- dat
  if (variable == "offset_argument") {
    dd$oa <- 0; dd$oa[4L] <- NA
    ff <- list(y ~ s(x, k = 4), ~ z)
  } else {
    dd[[variable]][4L] <- NA; dd$oa <- 0; ff <- form
  }
  for (method in c("NCV", "QNCV", "REML")) {
    blocked(gam(ff, data = dd, family = betals(), weights = w, offset = oa,
      method = method, nei = nb))
    blocked(gam(ff, data = dd, family = betals(), weights = w, offset = oa,
      subset = 3:62, method = method, nei = grouped(60L)))
  }
}
# Matrix responses and multiple simultaneous missing variables follow the same hook.
if (broken) {
  dd <- dat; dd$y <- cbind(rep(2L, n), rep(3L, n)); dd$y[2L, 2L] <- NA
  blocked(gam(form, data = dd, family = betabinomial(), method = "NCV", nei = nb))
  dd <- dat; dd$x[2L] <- NA; dd$o2[5L] <- NA; dd$w[8L] <- NA
  blocked(gam(form, data = dd, weights = w, family = betals(), method = "NCV", nei = nb))
}

# Successful end-to-end routes: method, convergence and actual NCV diagnostics.
check <- function(fit, method, ncv = TRUE) {
  stopifnot(identical(fit$method, method), all(is.finite(coef(fit))),
    all(is.finite(fit$Vp)), all(is.finite(fit$sp) & fit$sp > 0),
    length(fit$gcv.ubre) == 1L, is.finite(fit$gcv.ubre))
  if (!is.null(fit$converged)) stopifnot(isTRUE(fit$converged))
  if (!is.null(fit$outer.info$conv) && fit$outer.info$conv != "full convergence")
    stop(fit$family$family, " ", method, ": ", fit$outer.info$conv)
  if (ncv) stopifnot(length(fit$NCV) == 1L, is.finite(fit$NCV))
}
for (method in c("NCV", "REML")) {
  check(gam(form, data = dat, weights = w, family = betals(), method = method, nei = nb), method)
  dd <- dat; dd$x[1:3] <- NA
  whole <- gam(form, data = dd, weights = w, family = betals(), method = method, nei = nb)
  check(whole, method)
  filtered <- dd[complete.cases(dd), ]
  rebuilt <- gam(form, data = filtered, weights = w, family = betals(),
    method = method, nei = grouped(nrow(filtered)))
  check(rebuilt, method); near(coef(whole), coef(rebuilt)); near(whole$NCV, rebuilt$NCV)
  # Partial omission is usable after explicit filtering and rebuilding nei.
  dd <- dat; dd$x[2L] <- NA; dd$z[5L] <- NA; dd$w[8L] <- NA
  filtered <- dd[complete.cases(dd), ]
  check(gam(form, data = filtered, weights = w, family = betals(),
    method = method, nei = grouped(nrow(filtered))), method)
  # Default leave-one-out is constructed after omission; REML without nei is ordinary REML.
  dd <- dat; dd$x[2L] <- NA
  implicit <- gam(form, data = dd, weights = w, family = betals(), method = method)
  check(implicit, method, ncv = method == "NCV")
  filtered <- dd[complete.cases(dd), ]; m <- nrow(filtered)
  explicit <- gam(form, data = filtered, weights = w, family = betals(), method = method,
    nei = if (method == "NCV") list(a = 1:m, ma = 1:m, d = 1:m, md = 1:m) else NULL)
  check(explicit, method, ncv = method == "NCV")
  near(coef(implicit), coef(explicit)); near(implicit$gcv.ubre, explicit$gcv.ubre)
}
# Distinct deletion/prediction sets and repeated, reordered predictions retain
# their meaning on complete data (also covered numerically by existing NCV suites).
separate <- list(a = c(4:6, 10:12, 4:6, 19:21), ma = c(3L, 6L, 9L, 12L),
  d = c(5L, 11L, 12L, 5L, 20L), md = c(1L, 3L, 4L, 5L), jackknife = FALSE)
check(gam(form, data = dat, family = betals(), method = "NCV", nei = separate), "NCV")

# A precomputed G skips model-frame adjustment when it is fitted again. New
# neighbourhoods use G's retained rows, even when G$mf records past omission.
dd <- dat; dd$x[2L] <- NA
prefit <- gam(form, data = dd, family = betals(), fit = FALSE)
m <- nrow(prefit$mf)
check(gam(G = prefit, method = "NCV",
  nei = list(a = 1:m, ma = 1:m, d = 1:m, md = 1:m)), "NCV")

# Missing rows outside subset do not trigger adjustment. Explicit filtering
# after subset selection also leaves valid grouped neighbourhoods usable.
dd <- dat; dd$x[2L] <- NA
check(gam(form, data = dd, weights = w, family = betals(), subset = 4:90,
  method = "NCV", nei = grouped(87L)), "NCV")
dd <- dat[3:62, ]; dd$x[2L] <- NA; dd$o2[5L] <- NA; dd$w[8L] <- NA
dd <- dd[complete.cases(dd), ]
check(gam(form, data = dd, weights = w, family = betals(), method = "NCV",
  nei = grouped(nrow(dd))), "NCV")

# QNCV-capable general path: actual NCV/QNCV optimization with offsets, both
# whole-group omission and default leave-one-out after partial omission.
set.seed(23)
gg <- data.frame(x = runif(120L), z = runif(120L), a = runif(120L, -.1, .1),
                 b = runif(120L, -.1, .1), c = runif(120L, -.05, .05))
gg$y <- gen_gammalss()$rd(cbind(exp(.5 + .5*sin(6*gg$x) + gg$a),
  exp(-.7 + .3*gg$z + gg$b), exp(.2 + .1*gg$x + gg$c)))
gf <- list(y ~ s(x, k = 4) + offset(a), ~ z + offset(b), ~ x + offset(c))
for (method in c("NCV", "QNCV")) {
  check(gam(gf, data = gg, family = gen_gammalss(), method = method,
    nei = grouped(nrow(gg))), method)
  dd <- gg; dd$x[1:3] <- NA
  fit <- gam(gf, data = dd, family = gen_gammalss(), method = method,
    nei = grouped(nrow(dd)))
  check(fit, method)
  clean <- dd[complete.cases(dd), ]
  rebuilt <- gam(gf, data = clean, family = gen_gammalss(), method = method,
    nei = grouped(nrow(clean)))
  check(rebuilt, method); near(fit$NCV, rebuilt$NCV); near(coef(fit), coef(rebuilt))
  dd <- gg; dd$x[c(2L, 5L, 8L)] <- NA
  if (broken) blocked(gam(gf, data = dd, family = gen_gammalss(), method = method,
    nei = grouped(nrow(dd))))
  clean <- dd[complete.cases(dd), ]
  check(gam(gf, data = clean, family = gen_gammalss(), method = method,
    nei = grouped(nrow(clean))), method)
  dd <- gg; dd$x[2L] <- NA
  check(gam(gf, data = dd, family = gen_gammalss(), method = method), method)
}

# Supported extended-family routes remain available without the defective
# preprocessing: default LOO after omission and explicit rebuilt groups.
extended_nei <- function(n) {
  nb <- grouped(n); nb$d <- nb$a; nb
}
for (name in c("pig", "log_logistic", "cumulative_link")) {
  set.seed(912)
  fam <- switch(name, pig = pig(.6), log_logistic = log_logistic(.25),
                cumulative_link = cumulative_link(R = 3, link = "probit"))
  dd <- data.frame(x = runif(240L))
  eta <- .3 + sin(2*pi*dd$x)
  dd$y <- fam$rd(if (name == "cumulative_link") eta else exp(eta))
  for (method in c("NCV", "REML")) {
    fit <- gam(y ~ s(x, k = 6), data = dd, family = fam, method = method,
      nei = extended_nei(nrow(dd)))
    check(fit, method)
  }
  dd$x[c(2L, 5L, 8L)] <- NA
  check(gam(y ~ s(x, k = 6), data = dd, family = fam, method = "NCV"), "NCV")
  clean <- dd[complete.cases(dd), ]
  check(gam(y ~ s(x, k = 6), data = clean, family = fam, method = "NCV",
    nei = extended_nei(nrow(clean))), "NCV")
  if (packageVersion("mgcv") >= "1.9.4") {
    # BAM reports its working-model NCV criterion in gcv.ubre, not $NCV.
    check(bam(y ~ s(x, k = 6), data = dd, family = fam,
      method = "NCV", discrete = TRUE), "NCV", ncv = FALSE)
    check(bam(y ~ s(x, k = 6), data = clean, family = fam,
      method = "NCV", discrete = TRUE, nei = extended_nei(nrow(clean))), "NCV", ncv = FALSE)
  }
}

# Interception limits are tested separately from inlaws guard coverage.
# Extended families have no presetup; validation precedes preinitialize. BAM
# performs the same faulty adjustment, then validates before any family hook.
if (broken) {
  dd <- dat; dd$x[2L] <- NA
  extended <- list(pig(.4), log_logistic(.25), cumulative_link(R = 3),
                   cgamma(.4), clognormal(.5), cmp(1.5), gpoisson(2))
  for (fam in extended) {
    ran <- FALSE
    fam$preinitialize <- function(...) {ran <<- TRUE; stop("late hook reached")}
    for (method in c("NCV", "QNCV", "REML")) {
      msg <- error_message(gam(y ~ s(x, k = 4), data = dd, family = fam,
        method = method, nei = nb))
      stopifnot(!ran, identical(msg, "nei$ma faulty"))
    }
    msg <- error_message(bam(y ~ s(x, k = 4), data = dd, family = fam,
      method = "NCV", discrete = TRUE, nei = nb))
    stopifnot(!ran, identical(msg, "nei$ma does not match nei$a"))
  }
  # Incomplete lists fail inside nanei itself, even before general presetup.
  fam <- betals(); ran <- FALSE
  fam$check.nei <- function(...) {ran <<- TRUE; stop("guard reached")}
  msg <- error_message(gam(form, data = dd, family = fam, method = "NCV",
    nei = list(a = 1:n, ma = 1:n)))
  stopifnot(!ran, identical(msg, "full nei list needed if data incomplete"))
}

# The shared helper does not replace independent QNCV and extended-family guards.
for (fam in list(betals(), betabinomial())) {
  msg <- error_message(gam(form, data = dat, family = fam, method = "QNCV"))
  stopifnot(grepl("supports NCV but not QNCV", msg, fixed = TRUE))
}
if (packageVersion("mgcv") == "1.9.4") for (fam in list(cgamma(.4), clognormal(.5), cmp(1.5), gpoisson(2))) {
  for (method in c("NCV", "QNCV", "REML")) {
    msg <- error_message(gam(y ~ s(x, k = 4), data = dat, family = fam,
      method = method, nei = nb))
    stopifnot(grepl("upstream NCV bugs", msg, fixed = TRUE))
  }
}
