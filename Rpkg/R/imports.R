# Centralised import declarations.
# All @importFrom tags here are picked up by roxygen2 and written to NAMESPACE.

#' @importFrom dplyr `%>%` mutate filter case_when if_else inner_join
#'   bind_rows arrange left_join select relocate all_of add_count
#'   group_by group_split
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stats dist median var quantile p.adjust lowess pchisq
#'   coef vcov logLik sigma predict loess df.residual wilcox.test setNames
#' @importFrom utils write.table read.table
#' @importFrom graphics par segments points text lines polygon xspline xyinch
#' @importFrom stringr str_extract
#' @importFrom glmmTMB glmmTMB fixef
#' @importFrom VGAM vglm dirmultinomial
#' @importFrom igraph .from .to
#' @keywords internal
NULL

# Suppress R CMD check NOTEs for column names used in dplyr NSE expressions.
utils::globalVariables(c(
  # bipartition / n_choose_2 output column names
  "ref_ex_part", "ref_part_cnt", "setdiff1", "setdiff2",
  # shared event/gene columns
  "event", "gene", "n_samples",
  # multinomial / prec columns
  "count", "exon_part", "type",
  # rMATS mapping columns
  "dexseq_fragment", "position",
  # plot column
  "ex_or_in",
  # moderate_phi_trend columns
  "baseMean", "z", "var_z",
  # prec_estimate_vgam column
  "groups",
  # test_model_multinomial_vgam_wald_EB columns
  "coef_name", "comparison_option", "reference",
  "Estimate", "Std. Error", "z value", "Pr(>|z|)"
))
