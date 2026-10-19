#' @section NCV neighbourhood preprocessing:
#' In mgcv 1.9-4, model-frame omission can corrupt explicit `nei`
#' neighbourhoods before fitting. This is shared mgcv preprocessing, separate
#' from likelihood, deletion-update, QNCV and covariance restrictions. This
#' family checks the adjusted structure in `gam()`'s general-family `presetup`
#' hook and rejects inconsistent indices or cumulative group boundaries, also
#' when `nei` requests NCV calculations during REML. The check is not tied to
#' a presumed fixed mgcv version. Correctly adjusted whole-group omissions,
#' complete-data neighbourhoods and default leave-one-out (`nei = NULL`)
#' constructed after omission remain available within the family's supported
#' methods. Ordinary REML without `nei` is unaffected.
#'
#' Remove incomplete observations using all relevant response, predictor,
#' offset and weight variables, apply subset selection, then construct `nei`
#' using the resulting data's row indices. Deletion and prediction sets may
#' differ; their ordering, repetitions and grouping are not changed by inlaws.
#' If an incomplete `nei` list fails inside mgcv's adjustment itself, it fails
#' before any family hook and cannot receive the inlaws diagnostic. See the
#' family compatibility vignette for route and version qualifications.
