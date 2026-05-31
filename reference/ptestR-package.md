# ptestR: Permutation-Based Significance Testing for Regression Models

ptestR wraps [`stats::glm()`](https://rdrr.io/r/stats/glm.html),
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html), and binomial
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) with a permutation
loop to produce nonparametric p-values. For each model, a null
distribution of the test statistic is built by randomly rearranging the
outcome variable (`permNum` times); `p.perm` is the proportion of
permuted statistics at least as extreme as the observed one.

This approach requires far fewer distributional assumptions than
standard Wald or likelihood-ratio tests, making it well-suited to
neuroimaging, EEG, and other biomedical datasets with repeated measures
and small samples.

### Main functions

|  |  |  |
|----|----|----|
| Function | Model class | Test statistic |
| [`grouped_perm_glm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glm.md) | Generalised linear model (`glm`) | t |
| [`grouped_perm_glmm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glmm.md) | Linear mixed-effects model (`lmer`) | t |
| [`grouped_perm_binoglm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_binoglm.md) | Binomial logistic regression (`glm`) | z |

All three share the same call signature and return a tidy tibble.

### Typical grouped usage

    df |>
      dplyr::group_by(State) |>
      tidyr::nest() |>
      dplyr::mutate(
        res = purrr::map(data, ~ grouped_perm_glmm(
          .x,
          formla      = feature ~ PMA + Sex + (1 | subject),
          var_to_perm = "feature",
          permNum     = 10000
        ))
      ) |>
      tidyr::unnest(res) |>
      dplyr::filter(term == "PMA")

## References

França LGS, et al. (2024). Neonatal brain dynamic functional
connectivity in term and preterm infants and its association with early
childhood neurodevelopment. *Nature Communications*, 15, 16.
<https://doi.org/10.1038/s41467-023-44050-z>

## See also

Useful links:

- <https://github.com/CoDe-Neuro/ptestR>

- Report bugs at <https://github.com/CoDe-Neuro/ptestR/issues>

## Author

**Maintainer**: Lucas G. S. França <lucas.franca@kcl.ac.uk>
([ORCID](https://orcid.org/0000-0003-0853-1319))

Authors:

- Yan Ge <yan.ge@kcl.ac.uk>

- Dafnis Batalle <dafnis.batalle@kcl.ac.uk>
  ([ORCID](https://orcid.org/0000-0003-2097-979X))
