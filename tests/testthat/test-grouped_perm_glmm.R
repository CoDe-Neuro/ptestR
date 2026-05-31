library(ptestR)

# ── grouped_perm_glmm ────────────────────────────────────────────────────────

test_that("grouped_perm_glmm returns a tibble with expected columns", {
  skip_if_not_installed("lme4")

  # Minimal repeated-measures dataset
  set.seed(42)
  tbl <- data.frame(
    y       = rnorm(30),
    x       = rnorm(30),
    subject = rep(1:10, each = 3)
  )

  res <- grouped_perm_glmm(tbl, y ~ x + (1 | subject), "y",
                           permNum = 49, seed = 42)

  expect_s3_class(res, "tbl_df")
  expect_true(all(c("term", "effect", "estimate", "statistic", "p.perm") %in% names(res)))
})

test_that("grouped_perm_glmm returns only fixed effects rows", {
  skip_if_not_installed("lme4")

  set.seed(42)
  tbl <- data.frame(
    y       = rnorm(30),
    x       = rnorm(30),
    subject = rep(1:10, each = 3)
  )

  res <- grouped_perm_glmm(tbl, y ~ x + (1 | subject), "y",
                           permNum = 49, seed = 42)

  expect_true(all(res$effect == "fixed"))
})

test_that("grouped_perm_glmm p.perm is between 0 and 1", {
  skip_if_not_installed("lme4")

  set.seed(42)
  tbl <- data.frame(
    y       = rnorm(30),
    x       = rnorm(30),
    subject = rep(1:10, each = 3)
  )

  res <- grouped_perm_glmm(tbl, y ~ x + (1 | subject), "y",
                           permNum = 49, seed = 42)

  expect_true(all(res$p.perm >= 0 & res$p.perm <= 1))
})
