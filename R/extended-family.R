# Clone before bam retains separate local and G$family references. Replacing
# the family again in preinitialize would leave its final AIC using old theta.
.inlaws_extended_preinitialize <- function(y, family) {
  # mgcv sets qapprox for GAM NCV/QNCV, including NCV requested via nei.
  if (!is.null(family$qapprox) && utils::packageVersion("mgcv") == "1.9.4")
    stop("gam NCV/QNCV is unavailable with mgcv 1.9-4 due to upstream NCV bugs; ",
         "use REML or ML pending reassessment with the next mgcv release", call. = FALSE)
  if (!isTRUE(family$.inlaws_isolated))
    family <- unserialize(serialize(family, NULL))
  list(family = family)
}

#' @export
#' @importFrom mgcv fix.family.link
fix.family.link.inlaws_extended <- function(fam) {
  fam <- unserialize(serialize(fam, NULL))
  fam$.inlaws_isolated <- TRUE
  NextMethod()
}

.inlaws_extended_family <- function(family) {
  class(family) <- c("inlaws_extended", "extended.family", "family")
  family$preinitialize <- .inlaws_extended_preinitialize
  family
}
