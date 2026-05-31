# Permutation test for binomial logistic regression

`grouped_perm_binoglm` is a nonparametric test for binomial logistic
regression. It assesses the significance of coefficients by permutation,
computing the null distribution of the z-statistic by randomly
rearranging the binary outcome variable.

## Usage

``` r
grouped_perm_binoglm(tbl, formla, var_to_perm, permNum = 1000, seed = 42)
```

## Arguments

- tbl:

  A data frame or tibble containing all model variables.

- formla:

  A [formula](https://rdrr.io/r/stats/formula.html) describing the
  logistic regression model; passed to
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
  `family = binomial`.

- var_to_perm:

  Character. Name of the binary outcome column to permute.

- permNum:

  Integer. Number of permutations to generate. Default `1000`.

- seed:

  Integer. Random seed for reproducibility; passed to
  [`base::set.seed()`](https://rdrr.io/r/base/Random.html). Default
  `42`.

## Value

A tibble with one row per model term and columns:

- `term`:

  Name of the regression term.

- `estimate`:

  Estimated log-odds coefficient.

- `statistic`:

  Observed z-statistic.

- `p.value`:

  Asymptotic two-sided p-value from the fitted model.

- `p.perm`:

  Permutation p-value: proportion of permuted \|statistics\| \>=
  \|observed statistic\|. A value of `0` means no permuted statistic was
  as extreme; report as `p < 1/permNum`.

## Examples

``` r
counts <- c(10, 11, 8, 9, 6, 3, 5, 1)
gender <- c(0, 0, 0, 0, 0, 1, 1, 1)
TBL    <- data.frame(counts, gender)
grouped_perm_binoglm(TBL, gender ~ counts, "gender", permNum = 500, seed = 1)
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> # A tibble: 2 × 5
#>   term        estimate statistic p.value p.perm
#>   <chr>          <dbl>     <dbl>   <dbl>  <dbl>
#> 1 (Intercept)    246.   0.000723   0.999      1
#> 2 counts         -44.7 -0.000728   0.999      1
```
