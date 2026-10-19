# gam() calls nanei() after model.frame() and before general-family presetup.
# In mgcv 1.9-4 partial removals leave cumulative boundaries unchanged; an
# empty removal index also makes x[-integer(0)] discard the entire index set.
# Whole-group removal can be correct, so do not reject omission itself.
# Check the resulting structure rather than presuming a later version fixes it.
.inlaws_check_ncv_nei <- function(nei, omitted) {
  if (is.null(nei) || !length(omitted)) return(invisible(NULL))
  valid <- function(index, end) {
    is.numeric(index) && is.numeric(end) && length(end) > 0L &&
      all(is.finite(index)) && all(index >= 1 & index == round(index)) &&
      all(is.finite(end)) && all(end == round(end)) &&
      all(diff(c(0, end)) > 0) && utils::tail(end, 1L) == length(index)
  }
  if (length(nei$ma) != length(nei$md) ||
      !valid(nei$a, nei$ma) || !valid(nei$d, nei$md)) {
    version <- utils::packageDescription("mgcv", fields = "Version")
    stop("NCV/QNCV explicit nei neighbourhoods are unsafe after missing-row ",
      "omission with installed mgcv ", version, ": upstream neighbourhood ",
      "preprocessing left inconsistent indices or cumulative group boundaries. ",
      "Remove incomplete observations using all relevant response, predictor, ",
      "offset and weight variables (and apply subset selection), then construct ",
      "nei using the resulting data's row indices; or use REML without nei.",
      call. = FALSE)
  }
  invisible(NULL)
}

# Only general families have this supported early hook. preinitialize and
# initialize run after estimate.gam's neighbourhood validation. Extended
# families and bam therefore cannot use this integration (see the audit).
.inlaws_ncv_family <- function(family) {
  stopifnot(inherits(family, "general.family"), is.function(family$ncv))
  if (identical(family$check.nei, .inlaws_check_ncv_nei)) return(family)
  family$check.nei <- .inlaws_check_ncv_nei
  family$presetup <- as.expression(c(
    expression(family$check.nei(nei, attr(mf, "na.action"))),
    family$presetup))
  family
}
