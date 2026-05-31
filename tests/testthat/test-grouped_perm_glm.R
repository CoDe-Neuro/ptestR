library(ptestR)

# ── grouped_perm_glm ─────────────────────────────────────────────────────────

test_that("grouped_perm_glm returns a tibble with expected columns", {
  set.seed(1)
  counts    <- sample(1:100, 9, replace = TRUE)
  outcomes  <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  treatment <- gl(3, 3)
  tbl <- data.frame(counts, outcomes, treatment)

  res <- grouped_perm_glm(tbl, outcomes ~ counts + treatment, "outcomes",
                          permNum = 99, seed = 42)

  expect_s3_class(res, "tbl_df")
  expect_true(all(c("term", "estimate", "statistic", "p.value", "p.perm") %in% names(res)))
})

test_that("grouped_perm_glm p.perm is between 0 and 1", {
  tbl <- data.frame(x = rnorm(20), y = rnorm(20))
  res <- grouped_perm_glm(tbl, y ~ x, "y", permNum = 99, seed = 1)

  expect_true(all(res$p.perm >= 0 & res$p.perm <= 1))
})

test_that("grouped_perm_glm seed argument produces reproducible results", {
  tbl <- data.frame(x = rnorm(20), y = rnorm(20))
  r1  <- grouped_perm_glm(tbl, y ~ x, "y", permNum = 99, seed = 7)
  r2  <- grouped_perm_glm(tbl, y ~ x, "y", permNum = 99, seed = 7)

  expect_equal(r1$p.perm, r2$p.perm)
})
