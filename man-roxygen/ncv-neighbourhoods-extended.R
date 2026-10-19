#' @section NCV neighbourhood preprocessing:
#' In mgcv 1.9-4, both `gam()` and `bam()` can corrupt explicit `nei`
#' neighbourhoods when the model frame omits incomplete observations. Partial
#' group removal can leave stale cumulative boundaries; a set containing no
#' omitted index can be emptied accidentally. This shared preprocessing issue
#' is separate from family-specific NCV restrictions. In `gam()`, extended
#' families have no `presetup` hook; `preinitialize` and `initialize` run after upstream
#' neighbourhood validation. In `bam()`, early link fixing runs before the
#' model frame exists, while hooks receiving the processed data run after
#' neighbourhood validation. Neither route offers a supported family hook
#' for this diagnostic. An upstream validation error is not inlaws guard
#' coverage.
#'
#' For otherwise supported NCV routes, remove incomplete observations using
#' all relevant response, predictor, offset and weight variables, apply subset
#' selection, then construct `nei` using the resulting data's row indices.
#' Complete-data neighbourhoods, correctly adjusted whole-group omissions,
#' default leave-one-out (`nei = NULL`) constructed after omission, and ordinary
#' REML without `nei` do not encounter this defect. Other documented method
#' restrictions still apply. In particular, complete-case filtering does not
#' lift an independent GAM NCV restriction. No fix in a later mgcv version is
#' presumed; see the family compatibility vignette for the tested versions.
