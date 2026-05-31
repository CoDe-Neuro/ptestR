#' Permutation test for linear mixed-effects models
#'
#' @description `grouped_perm_glmm` is a nonparametric test for linear
#'   mixed-effects models. It assesses the significance of fixed-effect
#'   coefficients by permutation, computing the null distribution of the
#'   test statistic by randomly rearranging the outcome variable while
#'   preserving the random-effects structure.
#'
#' @param tbl A data frame or tibble containing all model variables.
#' @param formla A [formula] with both fixed- and random-effects parts;
#'   passed directly to [lme4::lmer()].
#' @param var_to_perm Character. Name of the column to permute (typically
#'   the outcome variable).
#' @param permNum Integer. Number of permutations to generate. Default `1000`.
#' @param seed Integer. Random seed for reproducibility; passed to
#'   [base::set.seed()]. Default `42`.
#'
#' @return A tibble with one row per fixed-effect term and columns:
#'   \describe{
#'     \item{`term`}{Name of the regression term.}
#'     \item{`effect`}{Always `"fixed"` (random-parameter rows are dropped).}
#'     \item{`estimate`}{Estimated coefficient.}
#'     \item{`statistic`}{Observed t-statistic.}
#'     \item{`p.perm`}{Permutation p-value: proportion of permuted |statistics|
#'       >= |observed statistic|. Replaces the conventional `p.value` because
#'       `lme4::lmer` does not compute degrees of freedom or p-values by
#'       default. A value of `0` means no permuted statistic was as extreme;
#'       report as `p < 1/permNum`.}
#'   }
#'
#' @export
#'
#' @importFrom modelr permute
#' @importFrom purrr map map_df map2
#' @importFrom broom.mixed tidy
#' @importFrom lme4 lmer
#' @importFrom dplyr filter group_by summarise ungroup inner_join mutate select
#' @importFrom tidyr all_of
#'
#' @examples
#' \dontrun{
#' library(sdamr)
#' data("anchoring")
#' grouped_perm_glmm(
#'   anchoring,
#'   everest_feet ~ anchor + sex + (1 | referrer),
#'   "everest_feet"
#' )
#' }
grouped_perm_glmm <- function(tbl, formla, var_to_perm,
                              permNum = 1000, seed = 42) {
  set.seed(seed)

  perms  <- modelr::permute(tbl, permNum, tidyr::all_of(var_to_perm))
  models <- purrr::map(perms$perm, ~ lme4::lmer(formla, data = .))

  tdy_idx <- purrr::map_df(models, broom.mixed::tidy, .id = "id")
  mod     <- broom.mixed::tidy(lme4::lmer(formla, data = tbl))

  tdy_stats <- tdy_idx |>
    dplyr::filter(effect == "fixed") |>
    dplyr::group_by(term, effect) |>
    dplyr::summarise(perm_stat = list(statistic), .groups = "drop")

  tdy_stats <- dplyr::inner_join(tdy_stats, mod, by = c("term", "effect"))

  tdy_stats |>
    dplyr::mutate(p.perm = unlist(purrr::map2(perm_stat, statistic, get_newp))) |>
    dplyr::select(term, effect, estimate, statistic, p.perm)
}
