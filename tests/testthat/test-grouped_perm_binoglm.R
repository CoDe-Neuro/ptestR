library(ptestR)

# ── grouped_perm_binoglm ─────────────────────────────────────────────────────

test_that("grouped_perm_binoglm returns a tibble with expected columns", {
  counts <- c(10, 11, 8, 9, 6, 3, 5, 1)
  gender <- c(0, 0, 0, 0, 0, 1, 1, 1)
  tbl    <- data.frame(counts, gender)

  res <- grouped_perm_binoglm(tbl, gender ~ counts, "gender",
                              permNum = 99, seed = 1)

  expect_s3_class(res, "tbl_df")
  expect_true(all(c("term", "estimate", "statistic", "p.value", "p.perm") %in% names(res)))
})

test_that("grouped_perm_binoglm p.perm is between 0 and 1", {
  counts <- c(10, 11, 8, 9, 6, 3, 5, 1)
  gender <- c(0, 0, 0, 0, 0, 1, 1, 1)
  tbl    <- data.frame(counts, gender)

  res <- grouped_perm_binoglm(tbl, gender ~ counts, "gender",
                              permNum = 99, seed = 1)

  expect_true(all(res$p.perm >= 0 & res$p.perm <= 1))
})

test_that("grouped_perm_binoglm seed argument produces reproducible results", {
  counts <- c(10, 11, 8, 9, 6, 3, 5, 1)
  gender <- c(0, 0, 0, 0, 0, 1, 1, 1)
  tbl    <- data.frame(counts, gender)

  r1 <- grouped_perm_binoglm(tbl, gender ~ counts, "gender",
                             permNum = 99, seed = 5)
  r2 <- grouped_perm_binoglm(tbl, gender ~ counts, "gender",
                             permNum = 99, seed = 5)

  expect_equal(r1$p.perm, r2$p.perm)
})
