test_that("scalar arguments are validated before file access", {
  expect_error(Admixture_ModPlot("prefix.", Kseq = c(2, 2)), "unique")
  expect_error(Admixture_ModPlot("prefix.", Kseq = 0), "positive integers")
  expect_error(Admixture_ModPlot("prefix.", Kseq = 2:3, KtoMod = 4), "one of the values")
  expect_error(Admixture_ModPlot("prefix.", Kseq = 2, lab.cex = 0), "positive number")
  expect_error(Admixture_ModPlot("prefix.", Kseq = 2, ToPDF = NA), "TRUE or FALSE")
  expect_error(Admixture_ModPlot("prefix.", AutoMatchColours = NA), "TRUE or FALSE")
})

test_that("ColourFile creates and extends deterministic mappings", {
  path <- tempfile("aidmixture-colours-")
  on.exit(unlink(path), add = TRUE)

  AIDmixture:::ColourFile(path, 2L)
  first <- data.table::fread(path, header = FALSE, fill = TRUE)
  expect_equal(as.integer(unlist(first[2, 1:2], use.names = FALSE)), 1:2)

  AIDmixture:::ColourFile(path, 2:4)
  extended <- data.table::fread(path, header = FALSE, fill = TRUE)
  expect_equal(nrow(extended), 4L)
  expect_equal(as.integer(unlist(extended[4, 1:4], use.names = FALSE)), 1:4)
})

test_that("Q dimensions must agree with the fam file and K", {
  dir <- tempfile("aidmixture-input-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  prefix <- file.path(dir, "run.")
  fam <- file.path(dir, "run.fam")
  writeLines(c("POP1 sample1 0 0 0 -9", "POP2 sample2 0 0 0 -9"), fam)

  data.table::fwrite(data.frame(V1 = c(0.6), V2 = c(0.4)), paste0(prefix, "2.Q"), col.names = FALSE)
  expect_error(
    Admixture_ModPlot(prefix, fam_file = fam, Kseq = 2L),
    "has 1 rows but the fam file has 2"
  )

  data.table::fwrite(
    data.frame(V1 = c(0.6, 0.2), V2 = c(0.4, 0.3), V3 = c(0, 0.5)),
    paste0(prefix, "2.Q"),
    col.names = FALSE
  )
  expect_error(
    Admixture_ModPlot(prefix, fam_file = fam, Kseq = 2L),
    "exactly 2 columns"
  )
})

test_that("sort files must include every population ID", {
  dir <- tempfile("aidmixture-sort-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  prefix <- file.path(dir, "run.")
  fam <- file.path(dir, "run.fam")
  sort <- file.path(dir, "order.txt")

  writeLines(c("POP1 sample1 0 0 0 -9", "POP2 sample2 0 0 0 -9"), fam)
  data.table::fwrite(
    data.frame(V1 = c(0.6, 0.2), V2 = c(0.4, 0.8)),
    paste0(prefix, "2.Q"),
    col.names = FALSE
  )
  writeLines("POP1", sort)

  expect_error(
    Admixture_ModPlot(prefix, fam_file = fam, Sort_file = sort, Kseq = 2L),
    "missing population IDs: POP2"
  )
})
