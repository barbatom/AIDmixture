test_that("default fresh colours maximize minimum Lab distance", {
  palette <- AIDmixture:::.default_colour_palette()
  used <- c(1L, 2L)

  selected <- AIDmixture:::.fresh_colour_indices(
    used,
    1L,
    length(palette)
  )

  canonical <- AIDmixture:::.canonical_colour_values(palette)
  lab <- AIDmixture:::.colours_to_lab(palette)
  candidates <- setdiff(seq_along(palette), used)
  candidates <- candidates[!canonical[candidates] %in% canonical[used]]
  candidates <- candidates[!duplicated(canonical[candidates])]

  min_distances <- vapply(
    candidates,
    function(candidate) {
      AIDmixture:::.lab_distance_to_set(lab, candidate, used)
    },
    numeric(1)
  )

  expect_equal(selected, candidates[which.max(min_distances)])
})

test_that("default fresh colours skip exact RGB duplicates", {
  palette <- AIDmixture:::.default_colour_palette()

  selected <- AIDmixture:::.fresh_colour_indices(
    8L,
    20L,
    length(palette)
  )

  selected_colours <- AIDmixture:::.canonical_colour_values(
    palette[c(8L, selected)]
  )

  expect_equal(anyDuplicated(selected_colours), 0L)
  expect_false(16L %in% selected)
})

test_that("custom palette lengths keep explicit index order", {
  expect_equal(
    AIDmixture:::.fresh_colour_indices(c(1L, 3L), 3L, 6L),
    c(2L, 4L, 5L)
  )
})
