# Permutation test for generalised linear models

`grouped_perm_glm` is a nonparametric test for generalised linear
models. It assesses the significance of coefficients by permutation,
computing the distribution of the test statistic by randomly rearranging
the outcome variable.

## Usage

``` r
grouped_perm_glm(
  tbl,
  formla,
  var_to_perm,
  family = gaussian,
  permNum = 1000,
  seed = 42
)
```

## Arguments

- tbl:

  A data frame or tibble containing all model variables.

- formla:

  A [formula](https://rdrr.io/r/stats/formula.html) describing the
  regression model to fit; passed directly to
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html).

- var_to_perm:

  Character. Name of the column to permute (typically the outcome
  variable).

- family:

  A description of the error distribution and link function; passed to
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html). Defaults to
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html).

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

  Estimated coefficient.

- `statistic`:

  Observed t-statistic.

- `p.value`:

  Asymptotic two-sided p-value from the fitted model.

- `p.perm`:

  Permutation p-value: proportion of permuted \|statistics\| \>=
  \|observed statistic\|. A value of `0` means no permuted statistic was
  as extreme; report as `p < 1/permNum`.

## Examples

``` r
counts    <- sample(1:100, 9, replace = TRUE)
outcomes  <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
treatment <- gl(3, 3)
TBL <- data.frame(counts, outcomes, treatment)
grouped_perm_glm(TBL, outcomes ~ counts + treatment, "outcomes")
#> # A tibble: 4 × 5
#>   term        estimate statistic   p.value p.perm
#>   <chr>          <dbl>     <dbl>     <dbl>  <dbl>
#> 1 (Intercept)   22.4       11.5  0.0000887  0.012
#> 2 counts        -0.195     -4.58 0.00595    0.011
#> 3 treatment2     2.27       1.04 0.346      0.34 
#> 4 treatment3     6.23       2.47 0.0567     0.058
```
