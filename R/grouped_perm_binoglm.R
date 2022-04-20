#' Permutation test for binomial logistic regression
#'
#' @description \code{perm_test_binoglm} is a nonparametric test of binomial logistic
#' regression models. It assesses the significance of coefficients by
#' permutation tests, which calculate the distribution of the test statistic by
#' randomly rearranging the observed data.
#'
#' @param tbl A data frame or a data frame extension (e.g. a tibble) containing
#'   the variables in the model.
#' @param formla An object of class "formula" containing a symbolic description
#'   of the logistic regression model to be fitted. see
#'   \code{\link[stats]{formula}} for more details.
#' @param var_to_perm The columns of variables to permute.
#' @param permNum The number of permutations to generate.
#' @param seed A single value, interpreted as an integer for initializing a
#'   random permutation generator. see \code{\link[base]{set.seed}} for more
#'   details.
#'
#' @return It returns a tibble containing the following components:
#'   \tabular{ll}{\code{term}\tab  The name of the regression term.\cr
#'   \tab \cr
#'   \code{estimate} \tab The estimated value of the regression term.\cr
#'   \tab \cr
#'   \code{statistic} \tab The value of a Z-statistic to use in a hypothesis 
#'   that the regression term is non-zero.\cr
#'   \tab \cr
#'   \code{p.value} \tab The two-sided p-value associated with the observed 
#'   statistic.\cr
#'   \tab \cr
#'   \code{p.perm} \tab The likelihood of observing the test statistic of the 
#'   original data among that of permutations.}
#'   
#' @export
#'
#' @examples
#' counts <- c(10, 11, 8, 9, 6, 3, 5, 1)
#' gender <- c(0,0,0,0,0,1,1,1)
#' TBL <- data.frame(counts, gender)
#' grouped_perm_binoglm(TBL, gender ~ counts, "gender")
#' grouped_perm_binoglm(TBL, gender ~ counts, "gender", permNum = 500, seed = 1)
#'
grouped_perm_binoglm <- function(tbl, formla, var_to_perm, permNum = 1000,
                                 seed = 42){

  require(modelr)
  require(dplyr)
  require(purrr)
  require(tidyr)
  require(broom)
  require(tibble)
  set.seed(seed)
  var_to_perm = as.factor(var_to_perm)
  perms <- modelr::permute(tbl, permNum, all_of(var_to_perm))
  models <- purrr::map(perms$perm,
                       ~ glm(formla, family = binomial,
                             data = .))

  tdy_idx <- purrr::map_df(models, broom::tidy, .id = "id")
  mod <- broom::tidy(glm(formla, family = binomial,
                         data = tbl))
  tdy_stats <- tdy_idx %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(perm_stat = list(statistic)) %>%
    dplyr::ungroup()
  tdy_stats <- dplyr::inner_join(tdy_stats, mod)
  final <- tdy_stats %>%
    dplyr::mutate(p.perm = map2(perm_stat, statistic, get_newp)) %>%
    dplyr::mutate(p.perm = unlist(p.perm)) %>%
    dplyr::select(term, estimate, statistic, p.value, p.perm)

  return(final)
}

get_newp <- function(perm_stat, tref){

  p.perm <- length(perm_stat[abs(perm_stat) >= abs(tref)])/length(perm_stat)

  return(p.perm)
}



