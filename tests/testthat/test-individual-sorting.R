test_that("individual sorting uses the dominant ancestry component", {
  toy_dir <- system.file("extdata", package = "AIDmixture")
  fam <- data.table::fread(
    file.path(toy_dir, "toy.fam"),
    select = 1L,
    header = FALSE,
    col.names = "IDs"
  )
  q4 <- data.table::fread(file.path(toy_dir, "toy.4.Q"), header = FALSE)
  sort_ids <- unique(as.character(fam$IDs))

  expect_equal(
    AIDmixture:::.order_individuals(as.character(fam$IDs), sort_ids),
    seq_len(nrow(fam))
  )

  expect_equal(
    AIDmixture:::.order_individuals(as.character(fam$IDs), sort_ids, q4),
    c(5L, 1L, 3L, 2L, 4L, 8L, 6L, 10L, 7L, 9L, 15L, 13L, 11L, 12L, 14L)
  )
})

test_that("individual sorting respects population order", {
  toy_dir <- system.file("extdata", package = "AIDmixture")
  fam <- data.table::fread(
    file.path(toy_dir, "toy.fam"),
    select = 1L,
    header = FALSE,
    col.names = "IDs"
  )
  q4 <- data.table::fread(file.path(toy_dir, "toy.4.Q"), header = FALSE)

  expect_equal(
    AIDmixture:::.order_individuals(
      as.character(fam$IDs),
      c("POP_C", "POP_A", "POP_B"),
      q4
    ),
    c(15L, 13L, 11L, 12L, 14L, 5L, 1L, 3L, 2L, 4L, 8L, 6L, 10L, 7L, 9L)
  )
})

test_that("SortIndividuals is a binary flag", {
  expect_error(
    Admixture_ModPlot("prefix.", SortIndividuals = NA),
    "TRUE or FALSE"
  )
  expect_error(
    Admixture_ModPlot("prefix.", SortIndividuals = 1),
    "TRUE or FALSE"
  )
})
