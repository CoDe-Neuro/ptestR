# ptestR (development version)

* Fixed the package-name wordmark rendering off-centre in some viewers;
  it's now rendered as vector artwork rather than live text.
* Fixed the hex sticker border being cropped at the image edges by adding
  canvas margin around the hex.
* Regenerated favicons from the corrected logo.

# ptestR 0.1.1

* Restructured package for CRAN-style hygiene: proper `DESCRIPTION`, `URL`, `BugReports`.
* Consolidated duplicated `get_newp()` helper into `R/utils.R`.
* Replaced `require()` calls with `@importFrom` declarations throughout.
* Added `@importFrom stats glm gaussian binomial` to silence R CMD check notes.
* Added `R/ptestR-package.R` with package-level documentation.
* Added `tests/testthat/` with 3 tests per exported function.
* Added `_pkgdown.yml` with Bootstrap 5 theme in the package colour palette.
* Added `.github/workflows/R-CMD-check.yaml` CI across ubuntu, macOS, and Windows.
* Added hex sticker logo (`man/figures/logo.svg`).

# ptestR 0.1.0

* Initial release.
* Exported `grouped_perm_glm()`, `grouped_perm_glmm()`, and `grouped_perm_binoglm()`.
* Used in França et al. (2024), *Nature Communications*, <https://doi.org/10.1038/s41467-023-44050-z>.
