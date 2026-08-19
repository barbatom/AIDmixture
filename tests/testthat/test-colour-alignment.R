test_that("automatic colour matching follows label switching and a split", {
  dir <- tempfile("aidmixture-colour-match-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  prefix <- file.path(dir, "run.")
  colour_file <- paste0(prefix, "ADMXcolors")

  ancestry_a <- c(0.95, 0.85, 0.75, 0.65, 0.15, 0.10, 0.05, 0.20)
  ancestry_b <- 1 - ancestry_a

  q2 <- data.frame(ancestry_a, ancestry_b)
  q3 <- data.frame(
    old_b = ancestry_b,
    old_a_major_child = 0.7 * ancestry_a,
    old_a_new_child = 0.3 * ancestry_a
  )

  data.table::fwrite(q2, paste0(prefix, "2.Q"), col.names = FALSE)
  data.table::fwrite(q3, paste0(prefix, "3.Q"), col.names = FALSE)

  AIDmixture:::ColourFile(
    colour_file,
    2:3,
    Q_file = prefix,
    AutoMatchColours = TRUE,
    palette_length = 6L
  )

  mapping <- data.table::fread(colour_file, header = FALSE, fill = TRUE)
  expect_equal(
    as.integer(unlist(mapping[2, 1:2], use.names = FALSE)),
    c(1L, 2L)
  )
  expect_equal(
    as.integer(unlist(mapping[3, 1:3], use.names = FALSE)),
    c(2L, 1L, 3L)
  )
})

test_that("manual colour rows are preserved and seed later automatic matching", {
  dir <- tempfile("aidmixture-colour-manual-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  prefix <- file.path(dir, "run.")
  colour_file <- paste0(prefix, "ADMXcolors")

  ancestry_a <- c(0.95, 0.85, 0.75, 0.65, 0.15, 0.10, 0.05, 0.20)
  ancestry_b <- 1 - ancestry_a

  data.table::fwrite(
    data.frame(ancestry_a, ancestry_b),
    paste0(prefix, "2.Q"),
    col.names = FALSE
  )
  data.table::fwrite(
    data.frame(
      old_b = ancestry_b,
      old_a_major_child = 0.7 * ancestry_a,
      old_a_new_child = 0.3 * ancestry_a
    ),
    paste0(prefix, "3.Q"),
    col.names = FALSE
  )

  AIDmixture:::ColourFile(
    colour_file,
    2L,
    Q_file = prefix,
    AutoMatchColours = TRUE,
    palette_length = 6L
  )

  mapping <- data.table::fread(colour_file, header = FALSE, fill = TRUE)
  mapping[2, 1:2] <- list(2L, 1L)
  data.table::fwrite(mapping, colour_file, sep = "\t", col.names = FALSE)

  AIDmixture:::ColourFile(
    colour_file,
    2:3,
    Q_file = prefix,
    AutoMatchColours = TRUE,
    palette_length = 6L
  )

  rematched <- data.table::fread(colour_file, header = FALSE, fill = TRUE)
  expect_equal(
    as.integer(unlist(rematched[2, 1:2], use.names = FALSE)),
    c(2L, 1L)
  )
  expect_equal(
    as.integer(unlist(rematched[3, 1:3], use.names = FALSE)),
    c(1L, 2L, 3L)
  )
})

test_that("automatic colour matching can be disabled", {
  dir <- tempfile("aidmixture-colour-sequential-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  prefix <- file.path(dir, "run.")
  colour_file <- paste0(prefix, "ADMXcolors")

  ancestry_a <- c(0.9, 0.8, 0.2, 0.1)
  ancestry_b <- 1 - ancestry_a
  data.table::fwrite(
    data.frame(ancestry_a, ancestry_b),
    paste0(prefix, "2.Q"),
    col.names = FALSE
  )
  data.table::fwrite(
    data.frame(ancestry_b, 0.7 * ancestry_a, 0.3 * ancestry_a),
    paste0(prefix, "3.Q"),
    col.names = FALSE
  )

  AIDmixture:::ColourFile(
    colour_file,
    2:3,
    Q_file = prefix,
    AutoMatchColours = FALSE,
    palette_length = 6L
  )

  mapping <- data.table::fread(colour_file, header = FALSE, fill = TRUE)
  expect_equal(
    as.integer(unlist(mapping[3, 1:3], use.names = FALSE)),
    1:3
  )
})
