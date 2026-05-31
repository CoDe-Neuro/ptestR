# Internal helpers — not exported

# Suppress R CMD check notes for dplyr column names used with NSE
utils::globalVariables(c(
  "term", "effect", "estimate", "statistic",
  "p.value", "p.perm", "perm_stat"
))

#' Compute permutation p-value
#'
#' @param perm_stat Numeric vector of test statistics from permuted models.
#' @param tref Observed test statistic.
#' @return Two-sided permutation p-value (proportion of |perm_stat| >= |tref|).
#' @noRd
get_newp <- function(perm_stat, tref) {
  length(perm_stat[abs(perm_stat) >= abs(tref)]) / length(perm_stat)
}
