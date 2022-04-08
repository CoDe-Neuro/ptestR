#' Permutation test for generalized linear models
#'
#' @description \code{perm_test_glm} is a nonparametric test of generalized
#' linear models. It assesses the significance of coefficients by
#' permutation tests, which calculate the distribution of the test statistic by
#' randomly rearranging the observed data.
#'
#' @param tbl A data frame or a data frame extension (e.g. a tibble) containing
#'   the variables in the model.
#' @param formla An object of class "formula" containing a symbolic description
#'   of the regression model to be fitted. See \code{\link[stats]{formula}} for 
#'   more details.
#' @param var_to_perm The columns of variables to permute.
#' @param family A description  of the error distribution and link function to 
#' be used in \code{glm}. See \code{\link[stats]{family}} for more details.
#' @param permNum The number of permutations to generate.
#' @param seed A single value interpreted as an integer for initializing a
#'   random permutation generator. See \code{\link[base]{set.seed}} for more
#'   details.
#'
#' @return It returns a tibble containing the following components:
#'   \tabular{ll}{\code{term}\tab  The name of the regression term.\cr
#'   \tab \cr
#'   \code{estimate} \tab The estimated value of the regression term.\cr
#'   \tab \cr
#'   \code{statistic} \tab The value of a test statistic to use in a hypothesis 
#'   that the regression term is non-zero.\cr
#'   \tab \cr
#'   \code{p.value} \tab The two-sided p-value associated with the observed 
#'   statistic.\cr
#'   \tab \cr
#'   \code{p.perm} \tab The likelihood of observing the test statistic of the 
#'   original data among that of permutations.}
#' @export
#'
#' @examples
#' counts <- sample(1:100, 9, replace=TRUE)
#' outcomes <- c(18,17,15,20,10,20,25,13,12)
#' treatment <- gl(3,3)
#' TBL <- data.frame(counts, outcomes, treatment)
#' grouped_perm_glm(TBL, outcomes ~ counts + treatment, "outcomes")
#'
#' counts <- sample(1:100, 9, replace=TRUE)
#' outcomes <- gl(2,1,3)
#' treatment <- gl(3,3)
#' TBL <- data.frame(counts, outcomes, treatment)
#' grouped_perm_glm(TBL, outcomes ~ counts + treatment, "outcomes",
#' family = binomial)
#'
grouped_perm_glm <- function(tbl, formla, var_to_perm, family = gaussian, permNum = 1000, seed = 42){

  require(modelr)
  require(dplyr)
  require(purrr)
  require(tidyr)
  require(broom)
  require(tibble)

  set.seed(seed)
  perms <- modelr::permute(tbl, permNum, all_of(var_to_perm))

  models <- purrr::map(perms$perm,
                       ~ glm(formla, family,
                             data = .))
  tdy_idx <- purrr::map_df(models, broom::tidy, .id = "id")

  mod <- broom::tidy(glm(formla, family,
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



