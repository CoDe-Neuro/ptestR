# 🧪 ptestR

An R package providing **permutation-based significance testing** for generalised linear models, linear mixed-effects models, and binomial logistic regression. Designed for small or non-normally distributed neuroscience datasets where standard asymptotic p-values are unreliable.

[![Nature Communications](https://img.shields.io/badge/Nature%20Comms-10.1038%2Fs41467--023--44050--z-4CAF50)](https://doi.org/10.1038/s41467-023-44050-z)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE.md)
[![R](https://img.shields.io/badge/R-≥4.0-276DC3.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.1.1-lightgrey)](https://github.com/CoDe-Neuro/ptestR/releases/tag/0.1.1)

---

## 📖 Overview

`ptestR` wraps `glm`, `lme4::lmer`, and binomial `glm` with a permutation loop: for each model, it generates a null distribution of the test statistic by randomly rearranging the outcome variable (`permNum` times), then computes `p.perm` as the proportion of permuted statistics that are at least as extreme as the observed one. This approach requires far fewer distributional assumptions than standard likelihood-ratio or Wald tests, making it well-suited to neuroimaging, EEG, and other biomedical datasets with repeated measures and small samples.

The package was developed for and is used in:

> França LGS, Ciarrusta J, Gale-Grant O, Fenn-Moltu S, Fitzgibbon S, Chew A, Falconer S, Dimitrova R, Cordero-Grande L, Price AN, Hughes E, O'Muircheartaigh J, Duff E, Tuulari JJ, Deco G, Counsell SJ, Hajnal JV, Nosarti C, Arichi T, Edwards AD, McAlonan G, Batalle D (2024).
> Neonatal brain dynamic functional connectivity in term and preterm infants and its association with early childhood neurodevelopment.
> *Nature Communications*, 15, 16. https://doi.org/10.1038/s41467-023-44050-z

In that study, `grouped_perm_glm` and `grouped_perm_glmm` were used with 10 000 permutations (two-sided, FDR-corrected at α = 5%) to test associations between neonatal fMRI brain state features and postmenstrual age, postnatal days, preterm birth, and neurodevelopmental outcomes (Bayley-III and Q-CHAT) at 18 months of age in 390 term- and preterm-born neonates from the developing Human Connectome Project (dHCP).

It is also used in:

> Ferré IBS, Corso G, dos Santos Lima GZ, Lopes SR, Leocadio-Miguel MA, França LGS, de Lima Prado T, Araújo JF (2024).
> Cycling reduces the entropy of neuronal activity in the human adult cortex.
> *PLOS ONE*, 19(10): e0298703. https://doi.org/10.1371/journal.pone.0298703

---

## 🚀 Installation

`ptestR` is not on CRAN. Install directly from GitHub:

```r
# install.packages("remotes")  # if needed
remotes::install_github("CoDe-Neuro/ptestR")
```

---

## 📦 Functions

The package exports three functions, one for each model class:

| Function | Model type | Underlying fit | Test statistic |
|----------|-----------|----------------|----------------|
| `grouped_perm_glm()` | Generalised linear model | `stats::glm` | t-statistic |
| `grouped_perm_glmm()` | Linear mixed-effects model | `lme4::lmer` | t-statistic |
| `grouped_perm_binoglm()` | Binomial logistic regression | `stats::glm(family = binomial)` | z-statistic |

All three share the same call signature and return a tidy tibble with one row per model term.

---

## 🔧 Function reference

### `grouped_perm_glm()`

Permutation test for generalised linear models.

```r
grouped_perm_glm(tbl, formla, var_to_perm, family = gaussian,
                 permNum = 1000, seed = 42)
```

| Argument | Type | Description |
|----------|------|-------------|
| `tbl` | data frame / tibble | Data containing all model variables |
| `formla` | formula | Model formula; passed to `glm` |
| `var_to_perm` | character | Name of the column to permute (typically the outcome) |
| `family` | family object | Error distribution and link function (default: `gaussian`) |
| `permNum` | integer | Number of permutations (default: `1000`) |
| `seed` | integer | Random seed for reproducibility (default: `42`) |

**Returns** a tibble with columns: `term`, `estimate`, `statistic`, `p.value`, `p.perm`.

**Example:**

```r
library(ptestR)

counts    <- sample(1:100, 9, replace = TRUE)
outcomes  <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
treatment <- gl(3, 3)
TBL <- data.frame(counts, outcomes, treatment)

grouped_perm_glm(TBL, outcomes ~ counts + treatment, "outcomes")
```

---

### `grouped_perm_glmm()`

Permutation test for linear mixed-effects models. The permutation is applied only to fixed effects; the random effects structure is preserved.

```r
grouped_perm_glmm(tbl, formla, var_to_perm, permNum = 1000, seed = 42)
```

| Argument | Type | Description |
|----------|------|-------------|
| `tbl` | data frame / tibble | Data containing all model variables |
| `formla` | formula | Mixed-effects formula with fixed and random parts; passed to `lme4::lmer` |
| `var_to_perm` | character | Name of the column to permute (typically the outcome) |
| `permNum` | integer | Number of permutations (default: `1000`) |
| `seed` | integer | Random seed for reproducibility (default: `42`) |

**Returns** a tibble with columns: `term`, `effect`, `estimate`, `statistic`, `p.perm`.

Note: `p.perm` replaces the conventional `p.value` because `lme4::lmer` does not compute degrees of freedom or p-values by default. The permutation p-value is a fully nonparametric alternative.

**Example:**

```r
library(ptestR)
library(sdamr)

data("anchoring")

grouped_perm_glmm(
  anchoring,
  everest_feet ~ anchor + sex + (1 | referrer),
  "everest_feet",
  permNum = 1000,
  seed = 42
)
```

**Typical usage pattern in a grouped analysis** (applied per brain state, as in the dHCP neonatal study):

```r
df %>%
  dplyr::group_by(State) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    mk = purrr::map(data, ~ grouped_perm_glmm(
      .x,
      formla      = feature ~ PMA + PND + Sex + motion_outliers + (1 | subject),
      var_to_perm = "feature",
      permNum     = 10000,
      seed        = 42
    ))
  ) %>%
  tidyr::unnest(mk) %>%
  dplyr::filter(term == "PMA")
```

---

### `grouped_perm_binoglm()`

Permutation test for binomial logistic regression.

```r
grouped_perm_binoglm(tbl, formla, var_to_perm, permNum = 1000, seed = 42)
```

| Argument | Type | Description |
|----------|------|-------------|
| `tbl` | data frame / tibble | Data containing all model variables |
| `formla` | formula | Logistic regression formula; passed to `glm(family = binomial)` |
| `var_to_perm` | character | Name of the binary outcome column to permute |
| `permNum` | integer | Number of permutations (default: `1000`) |
| `seed` | integer | Random seed for reproducibility (default: `42`) |

**Returns** a tibble with columns: `term`, `estimate`, `statistic`, `p.value`, `p.perm`.

**Example:**

```r
counts <- c(10, 11, 8, 9, 6, 3, 5, 1)
gender <- c(0, 0, 0, 0, 0, 1, 1, 1)
TBL    <- data.frame(counts, gender)

grouped_perm_binoglm(TBL, gender ~ counts, "gender", permNum = 500, seed = 1)
```

---

## 📐 Return value structure

All three functions return a tibble. The columns common to all functions are:

| Column | Description |
|--------|-------------|
| `term` | Name of the regression term (intercept or predictor) |
| `estimate` | Estimated coefficient from the fitted model |
| `statistic` | Observed test statistic (t for GLM/GLMM, z for binomial) |
| `p.value` | Asymptotic two-sided p-value (GLM and binoglm only; absent in GLMM) |
| `p.perm` | Permutation p-value: proportion of `permNum` permuted statistics ≥ \|observed\| |
| `effect` | `"fixed"` or `"ran_pars"` (GLMM only; results are filtered to `"fixed"`) |

`p.perm = 0` means no permuted statistic was as extreme as the observed one given `permNum` iterations; report as `p < 1/permNum` (e.g. `p < 0.0001` for `permNum = 10000`).

---

## 📦 Dependencies

| Package | Role |
|---------|------|
| `lme4` | Fits linear mixed-effects models (`lmer`) |
| `broom` | Tidies `glm` output |
| `broom.mixed` | Tidies `lmer` output |
| `modelr` | Generates permuted datasets (`permute`) |
| `dplyr` | Data manipulation |
| `purrr` | Functional iteration over permutations |
| `tidyr` | Reshaping outputs |
| `tibble` | Tidy data frames |

---

## 📄 Citation

If you use `ptestR` in your research, please cite the paper it was developed for:

```bibtex
@article{franca2024neonatal,
  title   = {Neonatal brain dynamic functional connectivity in term and preterm infants
             and its association with early childhood neurodevelopment},
  author  = {Fran{\c{c}}a, Lucas G. S. and Ciarrusta, Judit and Gale-Grant, Oliver and
             Fenn-Moltu, Sunniva and Fitzgibbon, Sean and Chew, Andrew and
             Falconer, Shona and Dimitrova, Ralica and Cordero-Grande, Lucilio and
             Price, Anthony N. and Hughes, Emer and O'Muircheartaigh, Jonathan and
             Duff, Eugene and Tuulari, Jetro J. and Deco, Gustavo and
             Counsell, Serena J. and Hajnal, Joseph V. and Nosarti, Chiara and
             Arichi, Tomoki and Edwards, A. David and McAlonan, Grainne and
             Batalle, Dafnis},
  journal = {Nature Communications},
  volume  = {15},
  pages   = {16},
  year    = {2024},
  doi     = {10.1038/s41467-023-44050-z}
}
```

---

## 👥 Authors

Lucas G. S. França, Yan Ge, and Dafnis Batallé — [CoDe-Neuro Lab](https://github.com/CoDe-Neuro), King's College London.

---

## 📄 Licence

[MIT](LICENSE.md) © 2022 Lucas G. S. França, Yan Ge, and Dafnis Batallé
