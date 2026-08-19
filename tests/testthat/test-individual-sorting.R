test_that("individual sorting uses the population-major ancestry component", {
  toy_dir <- system.file("extdata", package = "AIDmixture")
  fam <- data.table::fread(
    file.path(toy_dir, "toy.fam"),
    header = FALSE,
    select = 1L,
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

test_that("individual sorting preserves requested population order", {
  toy_dir <- system.file("extdata", package = "AIDmixture")
  fam <- data.table::fread(
    file.path(toy_dir, "toy.fam"),
    header = FALSE,
    select = 1L,
    col.names = "IDs"
  )
  q4 <- data.table::fread(file.path(toy_dir, "toy.4.Q"), header = FALSE)

  order_index <- AIDmixture:::.order_individuals(
    as.character(fam$IDs),
    c("POP_C", "POP_A", "POP_B"),
    q4
  )

  expect_equal(
    as.character(fam$IDs)[order_index],
    c(rep("POP_C", 5L), rep("POP_A", 5L), rep("POP_B", 5L))
  )
})

test_that("SortIndividuals is a logical flag", {
  expect_error(
    Admixture_ModPlot("prefix.", SortIndividuals = NA),
    "`SortIndividuals` must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    Admixture_ModPlot("prefix.", SortIndividuals = 1),
    "`SortIndividuals` must be TRUE or FALSE",
    fixed = TRUE
  )
})
