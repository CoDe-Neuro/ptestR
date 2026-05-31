# Permutation test for linear mixed-effects models

`grouped_perm_glmm` is a nonparametric test for linear mixed-effects
models. It assesses the significance of fixed-effect coefficients by
permutation, computing the null distribution of the test statistic by
randomly rearranging the outcome variable while preserving the
random-effects structure.

## Usage

``` r
grouped_perm_glmm(tbl, formla, var_to_perm, permNum = 1000, seed = 42)
```

## Arguments

- tbl:

  A data frame or tibble containing all model variables.

- formla:

  A [formula](https://rdrr.io/r/stats/formula.html) with both fixed- and
  random-effects parts; passed directly to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).

- var_to_perm:

  Character. Name of the column to permute (typically the outcome
  variable).

- permNum:

  Integer. Number of permutations to generate. Default `1000`.

- seed:

  Integer. Random seed for reproducibility; passed to
  [`base::set.seed()`](https://rdrr.io/r/base/Random.html). Default
  `42`.

## Value

A tibble with one row per fixed-effect term and columns:

- `term`:

  Name of the regression term.

- `effect`:

  Always `"fixed"` (random-parameter rows are dropped).

- `estimate`:

  Estimated coefficient.

- `statistic`:

  Observed t-statistic.

- `p.perm`:

  Permutation p-value: proportion of permuted \|statistics\| \>=
  \|observed statistic\|. Replaces the conventional `p.value` because
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) does not
  compute degrees of freedom or p-values by default. A value of `0`
  means no permuted statistic was as extreme; report as `p < 1/permNum`.

## Examples

``` r
if (FALSE) { # \dontrun{
library(sdamr)
data("anchoring")
grouped_perm_glmm(
  anchoring,
  everest_feet ~ anchor + sex + (1 | referrer),
  "everest_feet"
)
} # }
```
