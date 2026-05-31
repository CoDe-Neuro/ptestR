#' Permutation test for binomial logistic regression
#'
#' @description `grouped_perm_binoglm` is a nonparametric test for binomial
#'   logistic regression. It assesses the significance of coefficients by
#'   permutation, computing the null distribution of the z-statistic by
#'   randomly rearranging the binary outcome variable.
#'
#' @param tbl A data frame or tibble containing all model variables.
#' @param formla A [formula] describing the logistic regression model;
#'   passed to [stats::glm()] with `family = binomial`.
#' @param var_to_perm Character. Name of the binary outcome column to permute.
#' @param permNum Integer. Number of permutations to generate. Default `1000`.
#' @param seed Integer. Random seed for reproducibility; passed to
#'   [base::set.seed()]. Default `42`.
#'
#' @return A tibble with one row per model term and columns:
#'   \describe{
#'     \item{`term`}{Name of the regression term.}
#'     \item{`estimate`}{Estimated log-odds coefficient.}
#'     \item{`statistic`}{Observed z-statistic.}
#'     \item{`p.value`}{Asymptotic two-sided p-value from the fitted model.}
#'     \item{`p.perm`}{Permutation p-value: proportion of permuted |statistics|
#'       >= |observed statistic|. A value of `0` means no permuted statistic
#'       was as extreme; report as `p < 1/permNum`.}
#'   }
#'
#' @export
#'
#' @importFrom stats glm binomial
#' @importFrom modelr permute
#' @importFrom purrr map map_df map2
#' @importFrom broom tidy
#' @importFrom dplyr group_by summarise ungroup inner_join mutate select
#' @importFrom tidyr all_of
#'
#' @examples
#' counts <- c(10, 11, 8, 9, 6, 3, 5, 1)
#' gender <- c(0, 0, 0, 0, 0, 1, 1, 1)
#' TBL    <- data.frame(counts, gender)
#' grouped_perm_binoglm(TBL, gender ~ counts, "gender", permNum = 500, seed = 1)
grouped_perm_binoglm <- function(tbl, formla, var_to_perm,
                                 permNum = 1000, seed = 42) {
  set.seed(seed)

  var_to_perm <- as.factor(var_to_perm)
  perms  <- modelr::permute(tbl, permNum, tidyr::all_of(var_to_perm))
  models <- purrr::map(perms$perm, ~ glm(formla, family = binomial, data = .))

  tdy_idx <- purrr::map_df(models, broom::tidy, .id = "id")
  mod     <- broom::tidy(glm(formla, family = binomial, data = tbl))

  tdy_stats <- tdy_idx |>
    dplyr::group_by(term) |>
    dplyr::summarise(perm_stat = list(statistic)) |>
    dplyr::ungroup()

  tdy_stats <- dplyr::inner_join(tdy_stats, mod, by = "term")

  tdy_stats |>
    dplyr::mutate(p.perm = unlist(purrr::map2(perm_stat, statistic, get_newp))) |>
    dplyr::select(term, estimate, statistic, p.value, p.perm)
}
