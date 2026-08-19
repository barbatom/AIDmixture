test_that("label options are validated before file access", {
  expect_error(
    Admixture_ModPlot("prefix.", Labels = "Samples"),
    "`Labels` must be one of",
    fixed = TRUE
  )
  expect_error(
    Admixture_ModPlot("prefix.", IndividualLabelCex = 0),
    "`IndividualLabelCex` must be a single positive number",
    fixed = TRUE
  )
})

test_that("fam metadata exposes PLINK population and individual IDs", {
  toy_dir <- system.file("extdata", package = "AIDmixture")
  fam <- AIDmixture:::.read_fam_metadata(
    file.path(toy_dir, "toy.fam"),
    Labels = "Both"
  )

  expect_equal(fam$Population[1:6], c(rep("POP_A", 5L), "POP_B"))
  expect_equal(fam$Individual[1:6], c("A1", "A2", "A3", "A4", "A5", "B1"))
})

test_that("individual label modes require fam column two", {
  fam <- tempfile("aidmixture-one-column-fam-")
  on.exit(unlink(fam), add = TRUE)
  writeLines(c("POP_A", "POP_B"), fam)

  expect_silent(AIDmixture:::.read_fam_metadata(fam, Labels = "Population"))
  expect_silent(AIDmixture:::.read_fam_metadata(fam, Labels = "None"))
  expect_error(
    AIDmixture:::.read_fam_metadata(fam, Labels = "Individual"),
    "second column with individual IDs",
    fixed = TRUE
  )
  expect_error(
    AIDmixture:::.read_fam_metadata(fam, Labels = "Both"),
    "second column with individual IDs",
    fixed = TRUE
  )
})

test_that("Both uses structurally separate label panels", {
  population <- AIDmixture:::.label_panel_spec("Population")
  individual <- AIDmixture:::.label_panel_spec("Individual")
  both <- AIDmixture:::.label_panel_spec("Both")
  none <- AIDmixture:::.label_panel_spec("None")

  expect_equal(population$types, "Population")
  expect_equal(individual$types, "Individual")
  expect_equal(both$types, c("Individual", "Population"))
  expect_length(both$heights, 2L)
  expect_length(none$types, 0L)
  expect_length(none$heights, 0L)
})

test_that("barplot axis names are explicitly disabled", {
  plot_body <- paste(deparse(body(AIDmixture:::plotAdmixture_MB2)), collapse = "\n")
  expect_match(plot_body, "axisnames = FALSE", fixed = TRUE)
})
