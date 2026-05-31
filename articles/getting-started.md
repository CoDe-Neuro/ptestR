# Getting started with ptestR

## Why permutation tests?

Standard regression relies on asymptotic theory: p-values from Wald
tests or likelihood-ratio tests are valid when sample sizes are large
and residuals are approximately normal. In neuroscience and biomedical
research, these conditions are often not met — datasets are small,
outcomes are skewed, and repeated measures introduce complex dependency
structures.

Permutation tests sidestep these assumptions entirely. Instead of
comparing an observed statistic to a theoretical distribution, ptestR
builds the null distribution empirically: it randomly rearranges the
outcome variable many times, fits the model to each permuted dataset,
and asks how often the permuted statistic is at least as extreme as the
one we actually observed.

$`p_\text{perm} = \frac{\lvert\{|t_\text{perm}| \geq |t_\text{obs}|\}\rvert}{B}`$

where $`B`$ is the number of permutations. This value requires no
distributional assumptions and is exact in the sense that it converges
to the true p-value as $`B \to \infty`$.

------------------------------------------------------------------------

## The three functions

ptestR exports one function per model class, all sharing the same call
signature:

| Function | Model | Test statistic |
|----|----|----|
| [`grouped_perm_glm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glm.md) | Generalised linear model (`glm`) | t |
| [`grouped_perm_glmm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glmm.md) | Linear mixed-effects model (`lmer`) | t |
| [`grouped_perm_binoglm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_binoglm.md) | Binomial logistic regression | z |

This vignette covers
[`grouped_perm_glm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glm.md).
See
[`vignette("mixed-effects-models")`](https://code-neuro.github.io/ptestR/articles/mixed-effects-models.md)
for
[`grouped_perm_glmm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glmm.md)
and
[`vignette("grouped-analysis")`](https://code-neuro.github.io/ptestR/articles/grouped-analysis.md)
for the full nested workflow used in practice.

------------------------------------------------------------------------

## A simple example

We simulate a small dataset with a continuous outcome, a continuous
predictor, and a three-level grouping factor, then compare the
asymptotic and permutation p-values.

``` r

set.seed(42)
n  <- 30
df <- data.frame(
  outcome   = 2.5 * rnorm(n) + c(0, 1, 2)[rep(1:3, each = 10)],
  predictor = rnorm(n),
  group     = factor(rep(c("A", "B", "C"), each = 10))
)
head(df)
#>      outcome  predictor group
#> 1  3.4273961  0.4554501     A
#> 2 -1.4117454  0.7048373     A
#> 3  0.9078210  1.0351035     A
#> 4  1.5821565 -0.6089264     A
#> 5  1.0106708  0.5049551     A
#> 6 -0.2653113 -1.7170087     A
```

Fit a permutation GLM with 999 permutations:

``` r

res <- grouped_perm_glm(
  tbl         = df,
  formla      = outcome ~ predictor + group,
  var_to_perm = "outcome",
  permNum     = 999,
  seed        = 42
)
res
#> # A tibble: 4 × 5
#>   term        estimate statistic p.value p.perm
#>   <chr>          <dbl>     <dbl>   <dbl>  <dbl>
#> 1 (Intercept)    1.19      1.17    0.251  0.526
#> 2 groupB        -0.610    -0.430   0.671  0.670
#> 3 groupC         0.373     0.262   0.795  0.791
#> 4 predictor     -0.487    -0.864   0.395  0.395
```

The result is a tidy tibble with one row per model term. The key columns
are:

| Column      | Description                   |
|-------------|-------------------------------|
| `term`      | Regression term name          |
| `estimate`  | Fitted coefficient            |
| `statistic` | Observed t-statistic          |
| `p.value`   | Asymptotic p-value from `glm` |
| `p.perm`    | Permutation p-value           |

------------------------------------------------------------------------

## Interpreting `p.perm`

A `p.perm` of `0` does not mean the true p-value is zero — it means none
of the `permNum` permuted statistics were as extreme as the observed
one. Report it as:

    p < 1/permNum   (e.g. p < 0.001 for permNum = 1000)

The precision of `p.perm` depends on `permNum`. For exploratory work,
999 or 1 000 permutations are reasonable. For publication, 9 999 or 10
000 give better resolution, especially for small p-values.

``` r

# How many permutations do you need to resolve p = 0.05 with 10% relative error?
# Rule of thumb: permNum >= 10 / alpha
# For alpha = 0.05: permNum >= 200 (minimum); 1000+ recommended
# For alpha = 0.001: permNum >= 10000
```

------------------------------------------------------------------------

## The `family` argument

[`grouped_perm_glm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glm.md)
wraps [`stats::glm()`](https://rdrr.io/r/stats/glm.html), so you can
pass any family. The default is `gaussian` (ordinary least squares). For
count data you might use `poisson`:

``` r

set.seed(1)
counts_df <- data.frame(
  y = rpois(20, lambda = exp(0.5 * rnorm(20))),
  x = rnorm(20)
)

grouped_perm_glm(
  tbl         = counts_df,
  formla      = y ~ x,
  var_to_perm = "y",
  family      = poisson,
  permNum     = 499,
  seed        = 1
)
#> # A tibble: 2 × 5
#>   term        estimate statistic p.value p.perm
#>   <chr>          <dbl>     <dbl>   <dbl>  <dbl>
#> 1 (Intercept)  -0.0655    -0.277   0.782  0.361
#> 2 x             0.0951     0.278   0.781  0.737
```

------------------------------------------------------------------------

## Binomial logistic regression

For a binary outcome, use
[`grouped_perm_binoglm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_binoglm.md)
instead of passing `family = binomial` to
[`grouped_perm_glm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glm.md).
It handles the binary outcome correctly and returns z-statistics:

``` r

set.seed(7)
bin_df <- data.frame(
  outcome = rbinom(40, 1, prob = plogis(-0.5 + 0.8 * rnorm(40))),
  x1      = rnorm(40),
  x2      = rnorm(40)
)

grouped_perm_binoglm(
  tbl         = bin_df,
  formla      = outcome ~ x1 + x2,
  var_to_perm = "outcome",
  permNum     = 499,
  seed        = 7
)
#> # A tibble: 3 × 5
#>   term        estimate statistic p.value p.perm
#>   <chr>          <dbl>     <dbl>   <dbl>  <dbl>
#> 1 (Intercept)   -0.421    -1.19    0.236  0.222
#> 2 x1             0.327     0.884   0.377  0.359
#> 3 x2             0.323     0.820   0.412  0.437
```

------------------------------------------------------------------------

## Reproducibility

All three functions accept a `seed` argument passed to
[`base::set.seed()`](https://rdrr.io/r/base/Random.html). Always set
this in published analyses so results are exactly reproducible:

``` r

r1 <- grouped_perm_glm(df, outcome ~ predictor + group, "outcome",
                       permNum = 199, seed = 123)
r2 <- grouped_perm_glm(df, outcome ~ predictor + group, "outcome",
                       permNum = 199, seed = 123)

identical(r1$p.perm, r2$p.perm)
#> [1] TRUE
```

------------------------------------------------------------------------

## Multiple comparisons

`p.perm` values are raw — no multiple testing correction is applied
automatically. When running models across many features or brain
regions, apply FDR correction after collecting all results:

``` r

# After collecting results from many models into a single data frame `all_res`:
all_res$p.fdr <- p.adjust(all_res$p.perm, method = "fdr")
```

This matches the approach used in França et al. (2024), where 10 000
permutations and FDR correction at α = 5% were applied across brain
state features in 390 neonates.

------------------------------------------------------------------------

## Further reading

- [`vignette("mixed-effects-models")`](https://code-neuro.github.io/ptestR/articles/mixed-effects-models.md)
  — using
  [`grouped_perm_glmm()`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glmm.md)
  with `lmer`
- [`vignette("grouped-analysis")`](https://code-neuro.github.io/ptestR/articles/grouped-analysis.md)
  — nesting models across features or states
- [`?grouped_perm_glm`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glm.md),
  [`?grouped_perm_glmm`](https://code-neuro.github.io/ptestR/reference/grouped_perm_glmm.md),
  [`?grouped_perm_binoglm`](https://code-neuro.github.io/ptestR/reference/grouped_perm_binoglm.md)
