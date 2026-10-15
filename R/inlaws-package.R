#' Additional distribution families for mgcv
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Highly experimental and largely untested response families for models
#' fitted with [mgcv::gam()]. Interfaces and behaviour may change without
#' notice, and results require independent validation. See each family for its
#' supported fitting methods and parameterization.
#'
#' @importFrom mgcv gam
#' @importFrom stats coef fitted median model.matrix na.omit predict quantile
#'   residuals sd setNames var weighted.mean
#' @importFrom utils getFromNamespace
#' @keywords internal
"_PACKAGE"
