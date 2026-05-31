#' Permutation test for generalised linear models
#'
#' @description `grouped_perm_glm` is a nonparametric test for generalised
#'   linear models. It assesses the significance of coefficients by
#'   permutation, computing the distribution of the test statistic by
#'   randomly rearranging the outcome variable.
#'
#' @param tbl A data frame or tibble containing all model variables.
#' @param formla A [formula] describing the regression model to fit;
#'   passed directly to [stats::glm()].
#' @param var_to_perm Character. Name of the column to permute (typically
#'   the outcome variable).
#' @param family A description of the error distribution and link function;
#'   passed to [stats::glm()]. Defaults to [stats::gaussian()].
#' @param permNum Integer. Number of permutations to generate. Default `1000`.
#' @param seed Integer. Random seed for reproducibility; passed to
#'   [base::set.seed()]. Default `42`.
#'
#' @return A tibble with one row per model term and columns:
#'   \describe{
#'     \item{`term`}{Name of the regression term.}
#'     \item{`estimate`}{Estimated coefficient.}
#'     \item{`statistic`}{Observed t-statistic.}
#'     \item{`p.value`}{Asymptotic two-sided p-value from the fitted model.}
#'     \item{`p.perm`}{Permutation p-value: proportion of permuted |statistics|
#'       >= |observed statistic|. A value of `0` means no permuted statistic
#'       was as extreme; report as `p < 1/permNum`.}
#'   }
#'
#' @export
#'
#' @importFrom stats glm gaussian
#' @importFrom modelr permute
#' @importFrom purrr map map_df map2
#' @importFrom broom tidy
#' @importFrom dplyr group_by summarise ungroup inner_join mutate select
#' @importFrom tidyr all_of
#'
#' @examples
#' counts    <- sample(1:100, 9, replace = TRUE)
#' outcomes  <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
#' treatment <- gl(3, 3)
#' TBL <- data.frame(counts, outcomes, treatment)
#' grouped_perm_glm(TBL, outcomes ~ counts + treatment, "outcomes")
grouped_perm_glm <- function(tbl, formla, var_to_perm, family = gaussian,
                             permNum = 1000, seed = 42) {
  set.seed(seed)

  perms  <- modelr::permute(tbl, permNum, tidyr::all_of(var_to_perm))
  models <- purrr::map(perms$perm, ~ glm(formla, family, data = .))

  tdy_idx <- purrr::map_df(models, broom::tidy, .id = "id")
  mod     <- broom::tidy(glm(formla, family, data = tbl))

  tdy_stats <- tdy_idx |>
    dplyr::group_by(term) |>
    dplyr::summarise(perm_stat = list(statistic)) |>
    dplyr::ungroup()

  tdy_stats <- dplyr::inner_join(tdy_stats, mod, by = "term")

  tdy_stats |>
    dplyr::mutate(p.perm = unlist(purrr::map2(perm_stat, statistic, get_newp))) |>
    dplyr::select(term, estimate, statistic, p.value, p.perm)
}
