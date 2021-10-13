get_newp <- function(perm_stat, tref){
  
  p.perm <- length(perm_stat[abs(perm_stat) >= abs(tref)])/length(perm_stat)
  return(p.perm)
}

grouped_perm_glm <- function(tbl, formla, var_to_perm, permNum = 1000, seed = 42){
  
  require(modelr)
  require(dplyr)
  require(purrr)
  require(tidyr)
  require(broom)
  require(tibble)
  
  set.seed(seed)
  
  perms <- modelr::permute(tbl, permNum, all_of(var_to_perm))
  models <- purrr::map(perms$perm, 
                       ~ glm(formla, 
                             data = .))
  tdy_idx <- purrr::map_df(models, broom::tidy, .id = "id")
  
  mod <- broom::tidy(glm(formla, 
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