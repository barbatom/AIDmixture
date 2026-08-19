# Perceptual allocation for new built-in palette colours.
#
# This file loads after R/AIDmixture.R and replaces the basic sequential
# `.fresh_colour_indices()` fallback. Existing ADMXcolors indices are never
# changed; only newly introduced ancestry components use this allocator.

#' @noRd
.canonical_colour_values <- function(colours) {
  rgb <- grDevices::col2rgb(colours)
  apply(
    rgb,
    2L,
    function(channel) {
      sprintf(
        "#%02X%02X%02X",
        channel[1L],
        channel[2L],
        channel[3L]
      )
    }
  )
}

#' @noRd
.colours_to_lab <- function(colours) {
  rgb <- t(grDevices::col2rgb(colours)) / 255
  grDevices::convertColor(
    rgb,
    from = "sRGB",
    to = "Lab",
    scale.in = 1
  )
}

#' @noRd
.lab_distance_to_set <- function(lab, candidate, reference) {
  reference_lab <- lab[reference, , drop = FALSE]
  candidate_lab <- matrix(
    lab[candidate, ],
    nrow = nrow(reference_lab),
    ncol = ncol(reference_lab),
    byrow = TRUE
  )
  min(sqrt(rowSums((reference_lab - candidate_lab)^2)))
}

#' @noRd
.fresh_colour_indices <- function(used, n, palette_length) {
  used <- unique(as.integer(used))
  n <- as.integer(n)

  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop("The number of requested fresh colours must be a non-negative integer.", call. = FALSE)
  }
  if (n == 0L) {
    return(integer())
  }
  if (anyNA(used) || any(used < 1L) || any(used > palette_length)) {
    stop("The existing colour mapping contains an invalid palette index.", call. = FALSE)
  }

  default_palette <- .default_colour_palette()

  # Custom palettes keep their explicit ordering. The perceptual allocator is
  # applied to the built-in palette, whose historical indices must remain
  # stable because ADMXcolors files persist those indices between sessions.
  if (palette_length != length(default_palette)) {
    available <- setdiff(seq_len(palette_length), used)
    if (length(available) < n) {
      stop("The colour palette does not contain enough unused colours for automatic matching.", call. = FALSE)
    }
    return(available[seq_len(n)])
  }

  canonical <- .canonical_colour_values(default_palette)
  lab <- .colours_to_lab(default_palette)

  available <- setdiff(seq_len(palette_length), used)
  if (length(used) > 0L) {
    available <- available[!canonical[available] %in% canonical[used]]
  }
  available <- available[!duplicated(canonical[available])]

  if (length(available) < n) {
    stop(
      "The built-in colour palette does not contain enough perceptually distinct unused colours for automatic matching.",
      call. = FALSE
    )
  }

  selected <- integer()

  if (length(used) == 0L) {
    selected <- available[1L]
    selected_key <- canonical[selected]
    available <- available[canonical[available] != selected_key]
  }

  while (length(selected) < n) {
    reference <- c(used, selected)

    if (length(reference) == 0L) {
      best <- available[1L]
    } else {
      min_distances <- vapply(
        available,
        function(candidate) {
          .lab_distance_to_set(lab, candidate, reference)
        },
        numeric(1)
      )
      best <- available[which.max(min_distances)]
    }

    selected <- c(selected, best)
    best_key <- canonical[best]
    available <- available[canonical[available] != best_key]

    if (length(selected) < n && length(available) == 0L) {
      stop(
        "The built-in colour palette does not contain enough perceptually distinct unused colours for automatic matching.",
        call. = FALSE
      )
    }
  }

  as.integer(selected[seq_len(n)])
}
