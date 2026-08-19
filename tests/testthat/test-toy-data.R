test_that("bundled toy data is complete and valid", {
  toy_dir <- system.file("extdata", package = "AIDmixture")
  expect_true(nzchar(toy_dir))

  expected <- c("toy.fam", "toy.2.Q", "toy.3.Q", "toy.4.Q")
  source_files <- file.path(toy_dir, expected)
  expect_true(all(file.exists(source_files)))

  work <- tempfile("aidmixture-toy-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  expect_true(all(file.copy(source_files, work)))

  prefix <- file.path(work, "toy.")
  fam_file <- file.path(work, "toy.fam")
  prepared <- AIDmixture:::.prepare_input_files(prefix, fam_file, NULL, 2:4)
  expect_true(file.exists(prepared$Sort_file))

  fam <- data.table::fread(fam_file, header = FALSE)
  expect_equal(nrow(fam), 15L)
  expect_equal(length(unique(fam[[1L]])), 3L)

  for (K in 2:4) {
    q <- data.table::fread(paste0(prefix, K, ".Q"), header = FALSE)
    expect_equal(nrow(q), nrow(fam))
    expect_equal(ncol(q), K)
    expect_equal(rowSums(as.matrix(q)), rep(1, nrow(q)), tolerance = 1e-8)
  }
})
