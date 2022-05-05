#' Permutation test for linear mixed-effects models
#' 
#' @description \code{grouped_perm_glmm} is a nonparametric test of linear
#' mixed-effects models. It assesses the significance of coefficients by
#' permutation tests, which calculate the distribution of the test statistic by
#' randomly rearranging the observed data.
#' 
#' @param tbl A data frame or a data frame extension (e.g. a tibble) containing
#'   the variables in the model.
#' @param formla An object of class "formula" containing a symbolic description
#'   of the linear mixed-effects model to be fitted. The object include both the 
#'   fixed-effects and random-effects part of the model.See 
#'   \code{\link[lme4]{lmer}} for more details.
#' @param var_to_perm The columns of variables to permute.
#' @param permNum The number of permutations to generate.
#' @param seed A single value interpreted as an integer for initializing a
#'   random permutation generator. See \code{\link[base]{set.seed}} for more
#'   details.
#'
#' @return It returns a tibble containing the following components:
#'   \tabular{ll}{\code{term}\tab  The name of the regression term.\cr
#'   \tab \cr
#'   \code{effect}\tab  The type of the tested effect of the 
#'   regression term.\cr
#'   \tab \cr
#'   \code{estimate} \tab The estimated value of the regression term.\cr
#'   \tab \cr
#'   \code{statistic} \tab The value of a test statistic to use in a hypothesis 
#'   that the regression term is non-zero.\cr
#'   \tab \cr
#'   \code{p.perm} \tab The likelihood of observing the test statistic of the 
#'   original data among that of permutations.}
#' 
#' @export
#'
#' @examples
#' library(sdamr)
#' data("anchoring")
#' grouped_perm_glmm(anchoring,everest_feet ~ anchor+sex + (1|referrer),  
#' "everest_feet")
#' grouped_perm_glmm(anchoring,everest_feet ~ anchor+sex + (1|referrer),  
#' "everest_feet", permNum = 500, seed = 2)
#' 

grouped_perm_glmm <- function(tbl,
                              formla,
                              var_to_perm,
                              permNum = 1000,
                              seed = 42){

  require(purrr)
  require(modelr)
  require(dplyr)
  require(broom.mixed)
  require(lme4)
  
  set.seed(seed)
  
  perms <- modelr::permute(tbl,
                           permNum,
                           all_of(var_to_perm))
  
  models <- purrr::map(perms$perm,
                ~ lme4::lmer(formla,
                       data = .))
  
  tdy_idx <- purrr::map_df(models,
                    broom.mixed::tidy,
                    .id = "id")
  
  mod <- broom.mixed::tidy(lme4::lmer(formla,
                                      data = tbl))
  
  tdy_stats <- tdy_idx %>%
    dplyr::filter(effect == 'fixed') %>%
    group_by(term, effect) %>%
    summarise(perm_stat = list(statistic)) %>%
    ungroup()
  tdy_stats <- dplyr::inner_join(tdy_stats,
                                 mod)
  
  
  final <- tdy_stats %>%
    dplyr::mutate(p.perm = map2(perm_stat,
                                statistic,
                                get_newp)) %>%
    dplyr::mutate(p.perm = unlist(p.perm)) %>%
    dplyr::select(term,
                  effect,
                  estimate,
                  statistic,
                  p.perm)
  
  return(final)
}

get_newp <- function(perm_stat, tref){
  p.perm <- length(perm_stat[abs(perm_stat) >= abs(tref)])/length(perm_stat)
  return(p.perm)
}
