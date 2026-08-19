#' AIDmixture package imports
#'
#' @importFrom data.table fread fwrite
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline barplot layout locator par plot rect text
NULL

.default_colour_palette <- function() {
  c(
    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#E90C0C",
    "#7570B3", "#55D9F6", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
    "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3",
    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8",
    "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
    "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
  )
}

.validate_scalar_inputs <- function(
  Q_file,
  Kseq,
  ColourPalette,
  lab.cex,
  KtoMod,
  ToPDF,
  SortIndividuals,
  AutoMatchColours,
  Labels,
  IndividualLabelCex
) {
  if (!is.character(Q_file) || length(Q_file) != 1L || is.na(Q_file) || !nzchar(Q_file)) {
    stop("`Q_file` must be a non-empty character string.", call. = FALSE)
  }

  if (!is.numeric(Kseq) || length(Kseq) < 1L || anyNA(Kseq) ||
      any(!is.finite(Kseq)) || any(Kseq < 1) || any(Kseq != floor(Kseq))) {
    stop("`Kseq` must contain one or more positive integers.", call. = FALSE)
  }
  Kseq <- as.integer(Kseq)
  if (anyDuplicated(Kseq)) {
    stop("`Kseq` values must be unique.", call. = FALSE)
  }

  if (!is.numeric(lab.cex) || length(lab.cex) != 1L || is.na(lab.cex) ||
      !is.finite(lab.cex) || lab.cex <= 0) {
    stop("`lab.cex` must be a single positive number.", call. = FALSE)
  }

  if (!is.numeric(IndividualLabelCex) || length(IndividualLabelCex) != 1L ||
      is.na(IndividualLabelCex) || !is.finite(IndividualLabelCex) || IndividualLabelCex <= 0) {
    stop("`IndividualLabelCex` must be a single positive number.", call. = FALSE)
  }

  valid_labels <- c("Population", "Individual", "Both", "None")
  if (!is.character(Labels) || length(Labels) != 1L || is.na(Labels) || !Labels %in% valid_labels) {
    stop(
      "`Labels` must be one of: ", paste(valid_labels, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.numeric(KtoMod) || length(KtoMod) != 1L || is.na(KtoMod) ||
      !is.finite(KtoMod) || KtoMod < 0 || KtoMod != floor(KtoMod)) {
    stop("`KtoMod` must be a single non-negative integer.", call. = FALSE)
  }
  KtoMod <- as.integer(KtoMod)
  if (KtoMod > 0L && !KtoMod %in% Kseq) {
    stop("When greater than zero, `KtoMod` must be one of the values in `Kseq`.", call. = FALSE)
  }

  if (!is.logical(ToPDF) || length(ToPDF) != 1L || is.na(ToPDF)) {
    stop("`ToPDF` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(SortIndividuals) || length(SortIndividuals) != 1L || is.na(SortIndividuals)) {
    stop("`SortIndividuals` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(AutoMatchColours) || length(AutoMatchColours) != 1L || is.na(AutoMatchColours)) {
    stop("`AutoMatchColours` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(ColourPalette)) {
    if (!is.character(ColourPalette) || length(ColourPalette) < max(Kseq) || anyNA(ColourPalette)) {
      stop("`ColourPalette` must contain at least max(`Kseq`) valid colours.", call. = FALSE)
    }
    tryCatch(
      grDevices::col2rgb(ColourPalette),
      error = function(e) stop("`ColourPalette` contains an invalid R colour.", call. = FALSE)
    )
  }

  list(Kseq = Kseq, KtoMod = KtoMod)
}

#' @noRd
.read_fam_metadata <- function(fam_file, Labels = "Population") {
  fam_raw <- fread(fam_file, header = FALSE)
  if (nrow(fam_raw) < 1L || ncol(fam_raw) < 1L) {
    stop("The fam file is empty.", call. = FALSE)
  }

  population_ids <- as.character(fam_raw[[1L]])
  if (anyNA(population_ids) || any(!nzchar(population_ids))) {
    stop("The first column of the fam file contains missing or empty population IDs.", call. = FALSE)
  }

  need_individual_ids <- Labels %in% c("Individual", "Both")
  if (need_individual_ids && ncol(fam_raw) < 2L) {
    stop(
      "`fam_file` must contain a second column with individual IDs when `Labels` is ",
      sQuote(Labels), ".",
      call. = FALSE
    )
  }

  individual_ids <- rep(NA_character_, nrow(fam_raw))
  if (ncol(fam_raw) >= 2L) {
    individual_ids <- as.character(fam_raw[[2L]])
  }
  if (need_individual_ids && (anyNA(individual_ids) || any(!nzchar(individual_ids)))) {
    stop("The second column of the fam file contains missing or empty individual IDs.", call. = FALSE)
  }

  data.frame(
    Population = population_ids,
    Individual = individual_ids,
    stringsAsFactors = FALSE
  )
}

.prepare_input_files <- function(Q_file, fam_file, Sort_file, Kseq, Labels = "Population") {
  if (is.null(fam_file)) {
    fam_file <- paste0(Q_file, "fam")
  }
  if (!is.character(fam_file) || length(fam_file) != 1L || is.na(fam_file) || !file.exists(fam_file)) {
    stop("Could not open fam file: ", fam_file, call. = FALSE)
  }

  fam <- .read_fam_metadata(fam_file, Labels)

  for (K in Kseq) {
    q_path <- paste0(Q_file, K, ".Q")
    if (!file.exists(q_path)) {
      stop("Could not open Q file: ", q_path, call. = FALSE)
    }

    q <- fread(q_path, header = FALSE)
    if (nrow(q) != nrow(fam)) {
      stop(
        "Q file ", q_path, " has ", nrow(q), " rows but the fam file has ", nrow(fam), ".",
        call. = FALSE
      )
    }
    if (ncol(q) != K) {
      stop("Q file ", q_path, " must contain exactly ", K, " columns.", call. = FALSE)
    }
    if (!all(vapply(q, is.numeric, logical(1))) || any(!is.finite(as.matrix(q)))) {
      stop("Q file ", q_path, " must contain only finite numeric values.", call. = FALSE)
    }
  }

  if (is.null(Sort_file)) {
    Sort_file <- paste0(Q_file, "sort")
    if (!file.exists(Sort_file)) {
      fwrite(data.frame(IDs = unique(fam$Population)), Sort_file, col.names = FALSE)
    }
  }

  if (!is.character(Sort_file) || length(Sort_file) != 1L || is.na(Sort_file) || !file.exists(Sort_file)) {
    stop("Could not open sort file: ", Sort_file, call. = FALSE)
  }

  sort_table <- fread(Sort_file, header = FALSE, select = 1L)
  if (nrow(sort_table) < 1L) {
    stop("The sort file is empty.", call. = FALSE)
  }
  sort_ids <- as.character(sort_table[[1L]])
  if (anyNA(sort_ids) || any(!nzchar(sort_ids)) || anyDuplicated(sort_ids)) {
    stop("The sort file must contain unique, non-empty population IDs.", call. = FALSE)
  }

  missing_ids <- setdiff(unique(fam$Population), sort_ids)
  if (length(missing_ids) > 0L) {
    stop(
      "The sort file is missing population IDs: ", paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }

  list(fam_file = fam_file, Sort_file = Sort_file)
}

#' @noRd
.order_individuals <- function(fam_ids, sort_ids, q_ref = NULL) {
  if (is.null(q_ref)) {
    return(order(match(fam_ids, sort_ids), seq_along(fam_ids)))
  }

  q_ref <- as.matrix(q_ref)
  ordered <- lapply(sort_ids, function(population) {
    rows <- which(fam_ids == population)
    if (length(rows) == 0L) {
      return(integer())
    }

    population_q <- q_ref[rows, , drop = FALSE]
    major_component <- which.max(colMeans(population_q))
    rows[order(-population_q[, major_component], rows)]
  })

  as.integer(unlist(ordered, use.names = FALSE))
}

#' @noRd
.label_panel_spec <- function(Labels) {
  switch(
    Labels,
    Population = list(types = "Population", heights = 0.55),
    Individual = list(types = "Individual", heights = 0.75),
    Both = list(types = c("Individual", "Population"), heights = c(0.75, 0.55)),
    None = list(types = character(), heights = numeric())
  )
}

#' @noRd
.cosine_similarity <- function(x, y) {
  x_norm <- sqrt(sum(x * x))
  y_norm <- sqrt(sum(y * y))

  if (x_norm == 0 && y_norm == 0) {
    return(1)
  }
  if (x_norm == 0 || y_norm == 0) {
    return(0)
  }

  similarity <- sum(x * y) / (x_norm * y_norm)
  max(0, min(1, similarity))
}

#' @noRd
.component_similarity_matrix <- function(q_previous, q_current) {
  q_previous <- as.matrix(q_previous)
  q_current <- as.matrix(q_current)

  similarity <- matrix(
    0,
    nrow = ncol(q_previous),
    ncol = ncol(q_current)
  )

  for (i in seq_len(ncol(q_previous))) {
    for (j in seq_len(ncol(q_current))) {
      similarity[i, j] <- .cosine_similarity(q_previous[, i], q_current[, j])
    }
  }

  similarity
}

#' @noRd
.fresh_colour_indices <- function(used, n, palette_length) {
  available <- setdiff(seq_len(palette_length), unique(as.integer(used)))
  if (length(available) < n) {
    stop("The colour palette does not contain enough unused colours for automatic matching.", call. = FALSE)
  }
  available[seq_len(n)]
}

#' @noRd
.align_colour_row <- function(q_previous, q_current, previous_colours, palette_length) {
  q_previous <- as.matrix(q_previous)
  q_current <- as.matrix(q_current)
  previous_colours <- as.integer(previous_colours)

  previous_k <- ncol(q_previous)
  current_k <- ncol(q_current)

  if (nrow(q_previous) != nrow(q_current)) {
    stop("Q matrices used for colour matching must contain the same individuals.", call. = FALSE)
  }
  if (current_k <= previous_k) {
    stop("Automatic colour matching requires increasing K values.", call. = FALSE)
  }
  if (length(previous_colours) != previous_k || anyNA(previous_colours) || anyDuplicated(previous_colours)) {
    stop("The previous colour mapping is invalid.", call. = FALSE)
  }

  if (current_k == previous_k + 1L) {
    pairs <- utils::combn(seq_len(current_k), 2L)
    best_score <- -Inf
    best_assignment <- NULL
    best_groups <- NULL

    for (pair_index in seq_len(ncol(pairs))) {
      pair <- pairs[, pair_index]
      remaining <- setdiff(seq_len(current_k), pair)
      candidate <- cbind(
        rowSums(q_current[, pair, drop = FALSE]),
        q_current[, remaining, drop = FALSE]
      )
      groups <- c(list(pair), lapply(remaining, function(x) x))
      similarity <- .component_similarity_matrix(q_previous, candidate)
      assignment <- as.integer(clue::solve_LSAP(similarity, maximum = TRUE))
      score <- sum(similarity[cbind(seq_len(previous_k), assignment)])

      if (score > best_score + sqrt(.Machine$double.eps)) {
        best_score <- score
        best_assignment <- assignment
        best_groups <- groups
      }
    }

    mapping <- rep(NA_integer_, current_k)

    for (previous_component in seq_len(previous_k)) {
      group <- best_groups[[best_assignment[previous_component]]]
      colour <- previous_colours[previous_component]

      if (length(group) == 1L) {
        mapping[group] <- colour
      } else {
        child_similarity <- vapply(
          group,
          function(child) {
            .cosine_similarity(
              q_previous[, previous_component],
              q_current[, child]
            )
          },
          numeric(1)
        )
        keeper <- group[which.max(child_similarity)]
        mapping[keeper] <- colour
      }
    }

    new_components <- which(is.na(mapping))
    mapping[new_components] <- .fresh_colour_indices(
      previous_colours,
      length(new_components),
      palette_length
    )
    return(mapping)
  }

  similarity <- .component_similarity_matrix(q_previous, q_current)
  assignment <- as.integer(clue::solve_LSAP(similarity, maximum = TRUE))
  mapping <- rep(NA_integer_, current_k)

  for (previous_component in seq_len(previous_k)) {
    mapping[assignment[previous_component]] <- previous_colours[previous_component]
  }

  new_components <- which(is.na(mapping))
  mapping[new_components] <- .fresh_colour_indices(
    previous_colours,
    length(new_components),
    palette_length
  )
  mapping
}

#' @noRd
ColourFile <- function(
  colourADMX,
  Kseq,
  Q_file = NULL,
  AutoMatchColours = TRUE,
  palette_length = max(Kseq)
) {
  max_k <- max(Kseq)
  created <- !file.exists(colourADMX)

  if (created) {
    col_file <- as.data.frame(matrix(NA_integer_, max_k, max_k))
    col_file[1L, 1L] <- 1L
  } else {
    col_file <- as.data.frame(fread(colourADMX, header = FALSE, fill = TRUE))
  }

  changed <- created

  if (ncol(col_file) < max_k) {
    col_file <- cbind(
      col_file,
      as.data.frame(matrix(NA_integer_, nrow(col_file), max_k - ncol(col_file)))
    )
    names(col_file) <- paste0("V", seq_len(ncol(col_file)))
    changed <- TRUE
  }

  if (nrow(col_file) < max_k) {
    new_rows <- as.data.frame(
      matrix(NA_integer_, max_k - nrow(col_file), ncol(col_file))
    )
    names(new_rows) <- names(col_file)
    col_file <- rbind(col_file, new_rows)
    changed <- TRUE
  }

  ordered_k <- sort(as.integer(Kseq))

  for (position in seq_along(ordered_k)) {
    K <- ordered_k[position]
    existing <- as.integer(unlist(col_file[K, seq_len(K), drop = FALSE], use.names = FALSE))

    if (all(!is.na(existing))) {
      next
    }
    if (any(!is.na(existing))) {
      next
    }

    if (!AutoMatchColours || is.null(Q_file) || position == 1L) {
      mapping <- seq_len(K)
    } else {
      previous_candidates <- ordered_k[seq_len(position - 1L)]
      previous_candidates <- previous_candidates[vapply(
        previous_candidates,
        function(previous_k) {
          values <- as.integer(unlist(
            col_file[previous_k, seq_len(previous_k), drop = FALSE],
            use.names = FALSE
          ))
          all(!is.na(values))
        },
        logical(1)
      )]

      if (length(previous_candidates) == 0L) {
        mapping <- seq_len(K)
      } else {
        previous_k <- max(previous_candidates)
        previous_colours <- as.integer(unlist(
          col_file[previous_k, seq_len(previous_k), drop = FALSE],
          use.names = FALSE
        ))
        q_previous <- fread(paste0(Q_file, previous_k, ".Q"), header = FALSE)
        q_current <- fread(paste0(Q_file, K, ".Q"), header = FALSE)
        mapping <- .align_colour_row(
          q_previous,
          q_current,
          previous_colours,
          palette_length
        )
      }
    }

    col_file[K, seq_len(K)] <- as.list(as.integer(mapping))
    changed <- TRUE
  }

  if (changed) {
    fwrite(col_file, colourADMX, sep = "\t", col.names = FALSE)
  }

  invisible(col_file)
}

.validate_colour_file <- function(col_file, Kseq, palette_length) {
  for (K in Kseq) {
    if (nrow(col_file) < K || ncol(col_file) < K) {
      stop("The colour file does not contain enough rows or columns for K = ", K, ".", call. = FALSE)
    }

    indices <- as.integer(unlist(col_file[K, seq_len(K), with = FALSE], use.names = FALSE))
    if (anyNA(indices) || any(indices < 1L) || any(indices > palette_length) || anyDuplicated(indices)) {
      stop("The colour mapping for K = ", K, " is invalid. Delete the ADMXcolors file to regenerate it.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' @noRd
plotAdmixture_MB2 <- function(
  Q_is,
  Fam_is,
  Sort_is,
  Col_is,
  Kseq,
  toPDF = FALSE,
  colorPal,
  margin = 0.05,
  lab.cex = 1,
  modding = FALSE,
  KtoMod = 0L,
  SortIndividuals = FALSE,
  Labels = "Population",
  IndividualLabelCex = 0.7
) {
  if (toPDF) {
    pdf(paste0(Q_is, "ADMX.pdf"), paper = "a4r")
    on.exit(dev.off(), add = TRUE)
  }

  fam <- .read_fam_metadata(Fam_is, Labels)
  sort_table <- fread(Sort_is, header = FALSE, select = 1L)
  sort_ids <- as.character(sort_table[[1L]])

  q_ref <- NULL
  if (SortIndividuals) {
    q_ref <- fread(paste0(Q_is, max(Kseq), ".Q"), header = FALSE)
  }
  individual_order <- .order_individuals(fam$Population, sort_ids, q_ref)
  ordered_populations <- fam$Population[individual_order]
  ordered_individuals <- fam$Individual[individual_order]

  fam_change <- 0L
  if (length(ordered_populations) > 1L) {
    fam_change <- c(
      0L,
      which(ordered_populations[-1L] != ordered_populations[-length(ordered_populations)])
    )
  }

  n_plots <- length(Kseq)
  label_spec <- .label_panel_spec(Labels)
  extra_heights <- label_spec$heights
  if (modding) {
    extra_heights <- c(extra_heights, 0.5)
  }
  total_panels <- n_plots + length(extra_heights)
  layout(
    mat = matrix(seq_len(total_panels), ncol = 1L),
    heights = c(rep(1, n_plots), extra_heights)
  )
  par(mai = c(0, 0, margin, 0))

  for (K in Kseq) {
    q <- fread(paste0(Q_is, K, ".Q"), header = FALSE)
    q_plot <- as.data.frame(q)[individual_order, , drop = FALSE]
    row.names(q_plot) <- NULL

    colour_indices <- as.integer(unlist(Col_is[K, seq_len(K), with = FALSE], use.names = FALSE))
    barplot(
      t(as.matrix(q_plot)),
      width = 1,
      border = NA,
      col = colorPal[colour_indices],
      axes = FALSE,
      axisnames = FALSE,
      space = 0
    )
    text(0, 0.5, K, font = 2, pos = 2)
    abline(v = fam_change)
  }

  label_cex <- if (lab.cex == 1) n_plots^(-0.05) else lab.cex
  n_individuals <- length(ordered_populations)

  for (label_type in label_spec$types) {
    plot(
      NA,
      xlim = c(0, n_individuals),
      ylim = c(0, 1),
      axes = FALSE,
      xlab = "",
      ylab = ""
    )

    if (label_type == "Individual") {
      text(
        seq_len(n_individuals) - 0.5,
        0.95,
        labels = ordered_individuals,
        srt = 90,
        adj = c(1, 0.5),
        cex = IndividualLabelCex,
        xpd = FALSE
      )
    } else {
      starts <- fam_change
      ends <- c(fam_change[-1L], n_individuals)
      midpoints <- starts + (ends - starts) / 2
      population_labels <- ordered_populations[starts + 1L]
      text(
        midpoints,
        0.95,
        labels = population_labels,
        srt = 60,
        adj = c(1, 0.5),
        cex = label_cex,
        xpd = FALSE
      )
    }
  }

  if (modding) {
    plot(NA, xlim = c(0, KtoMod), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    text(0, 0.5, KtoMod, font = 2, pos = 2)
    rect(seq.int(0, KtoMod - 1L), 0, seq_len(KtoMod), 1, col = colorPal[seq_len(KtoMod)])
  }

  invisible(NULL)
}

#' Plot and modify ADMIXTURE results
#'
#' Plot ADMIXTURE Q-matrix results for one or more values of K. Population
#' ordering can be supplied with a one-column sort file. Colour assignments are
#' persisted in an `ADMXcolors` file and can be matched automatically between K
#' values or changed interactively by selecting two colour swatches.
#'
#' @param Q_file Character scalar giving the prefix shared by the ADMIXTURE
#'   files. For example, `"results/run."` expects `results/run.2.Q` for K = 2.
#' @param fam_file Path to a PLINK-style fam file. The first column is used as
#'   the population or group identifier. The second column is used as the
#'   individual identifier when `Labels = "Individual"` or `Labels = "Both"`.
#'   Defaults to `paste0(Q_file, "fam")`.
#' @param Sort_file Optional one-column file defining population plotting order.
#'   If omitted, `paste0(Q_file, "sort")` is created from the fam file when it
#'   does not already exist.
#' @param Kseq Positive integer vector containing the K values to plot.
#' @param ColourPalette Optional character vector of R colours. It must contain
#'   at least `max(Kseq)` entries. A built-in qualitative palette is used by
#'   default.
#' @param lab.cex Positive numeric scalar controlling population-label size.
#' @param KtoMod Non-negative integer. Set to one of the plotted K values to
#'   open the interactive colour-swap panel; use 0 to disable editing.
#' @param ToPDF Logical scalar. If `TRUE`, write the final plot to
#'   `paste0(Q_file, "ADMX.pdf")`.
#' @param SortIndividuals Logical scalar. If `TRUE`, identify the dominant
#'   ancestry component for each population at `max(Kseq)` from the population
#'   mean Q values, then sort individuals within that population from highest
#'   to lowest membership in that component. The same individual order is used
#'   for every plotted K. Defaults to `FALSE`.
#' @param AutoMatchColours Logical scalar. If `TRUE`, newly generated colour
#'   mappings align ancestry components from lower to higher K using
#'   individual-level Q profiles. For consecutive K values, the matcher
#'   explicitly evaluates one-component splits before applying a maximum-weight
#'   assignment. Existing `ADMXcolors` rows, including manual swaps, are never
#'   overwritten. Defaults to `TRUE`.
#' @param Labels Character scalar controlling labels below the plot. Use
#'   `"Population"` (default), `"Individual"`, `"Both"`, or `"None"`. When
#'   `"Both"` is selected, individual and population IDs are drawn in separate
#'   panels so they cannot overlap.
#' @param IndividualLabelCex Positive numeric scalar controlling individual-ID
#'   label size. Defaults to `0.7`.
#'
#' @return Invisibly, the path to the persistent `ADMXcolors` mapping file.
#' @export
#'
#' @examples
#' toy_dir <- system.file("extdata", package = "AIDmixture")
#' work <- tempfile("AIDmixture-example-")
#' dir.create(work)
#' invisible(file.copy(
#'   list.files(toy_dir, pattern = "^toy[.]", full.names = TRUE),
#'   work
#' ))
#'
#' Admixture_ModPlot(
#'   Q_file = file.path(work, "toy."),
#'   fam_file = file.path(work, "toy.fam"),
#'   Kseq = 2:4,
#'   SortIndividuals = TRUE
#' )
#'
#' unlink(work, recursive = TRUE)
#'
#' \dontrun{
#' # Plot your own ADMIXTURE results.
#' Admixture_ModPlot(
#'   Q_file = "results/run.",
#'   fam_file = "results/run.fam",
#'   Kseq = 2:6
#' )
#'
#' # Show both individual IDs and population IDs on separate label tiers.
#' Admixture_ModPlot(
#'   Q_file = "results/run.",
#'   fam_file = "results/run.fam",
#'   Kseq = 2:6,
#'   Labels = "Both"
#' )
#'
#' # Interactively swap two colours for K = 4.
#' Admixture_ModPlot(
#'   Q_file = "results/run.",
#'   fam_file = "results/run.fam",
#'   Kseq = 2:6,
#'   KtoMod = 4
#' )
#' }
Admixture_ModPlot <- function(
  Q_file,
  fam_file = NULL,
  Sort_file = NULL,
  Kseq = 2L,
  ColourPalette = NULL,
  lab.cex = 1,
  KtoMod = 0L,
  ToPDF = FALSE,
  SortIndividuals = FALSE,
  AutoMatchColours = TRUE,
  Labels = "Population",
  IndividualLabelCex = 0.7
) {
  checked <- .validate_scalar_inputs(
    Q_file,
    Kseq,
    ColourPalette,
    lab.cex,
    KtoMod,
    ToPDF,
    SortIndividuals,
    AutoMatchColours,
    Labels,
    IndividualLabelCex
  )
  Kseq <- checked$Kseq
  KtoMod <- checked$KtoMod

  files <- .prepare_input_files(Q_file, fam_file, Sort_file, Kseq, Labels)
  fam_file <- files$fam_file
  Sort_file <- files$Sort_file

  if (is.null(ColourPalette)) {
    ColourPalette <- .default_colour_palette()
  }
  if (length(ColourPalette) < max(Kseq)) {
    stop("The colour palette must contain at least max(`Kseq`) colours.", call. = FALSE)
  }

  colour_file_name <- paste0(Q_file, "ADMXcolors")
  ColourFile(
    colour_file_name,
    Kseq,
    Q_file = Q_file,
    AutoMatchColours = AutoMatchColours,
    palette_length = length(ColourPalette)
  )
  colour_file <- fread(colour_file_name, header = FALSE, fill = TRUE)
  .validate_colour_file(colour_file, Kseq, length(ColourPalette))

  mod_plot <- KtoMod > 0L
  plotAdmixture_MB2(
    Q_is = Q_file,
    Fam_is = fam_file,
    Sort_is = Sort_file,
    Col_is = colour_file,
    colorPal = ColourPalette,
    lab.cex = lab.cex,
    Kseq = Kseq,
    modding = mod_plot,
    KtoMod = KtoMod,
    toPDF = ToPDF && !mod_plot,
    SortIndividuals = SortIndividuals,
    Labels = Labels,
    IndividualLabelCex = IndividualLabelCex
  )

  if (mod_plot) {
    cols <- as.integer(unlist(colour_file[KtoMod, seq_len(KtoMod), with = FALSE], use.names = FALSE))

    message("Select two colours to swap.")
    coords <- locator(2)
    if (is.null(coords) || length(coords$x) < 2L || length(coords$y) < 2L) {
      warning("Colour selection was cancelled.", call. = FALSE)
      return(invisible(colour_file_name))
    }

    in_bounds <- coords$y[1:2] >= 0 & coords$y[1:2] <= 1 &
      coords$x[1:2] > 0 & coords$x[1:2] <= KtoMod
    if (!all(in_bounds)) {
      warning("Selected coordinates are outside the colour panel.", call. = FALSE)
      return(invisible(colour_file_name))
    }

    col1 <- as.integer(ceiling(coords$x[1]))
    col2 <- as.integer(ceiling(coords$x[2]))
    if (col1 == col2) {
      warning("The same colour was selected twice; nothing was changed.", call. = FALSE)
      return(invisible(colour_file_name))
    }

    message("Swapping ", ColourPalette[col1], " and ", ColourPalette[col2], ".")
    swapped <- cols
    swapped[cols == col1] <- col2
    swapped[cols == col2] <- col1

    for (v in seq_len(KtoMod)) {
      colour_file[KtoMod, v] <- swapped[v]
    }
    fwrite(colour_file, colour_file_name, sep = "\t", col.names = FALSE)

    plotAdmixture_MB2(
      Q_is = Q_file,
      Fam_is = fam_file,
      Sort_is = Sort_file,
      Col_is = colour_file,
      colorPal = ColourPalette,
      lab.cex = lab.cex,
      Kseq = Kseq,
      modding = mod_plot,
      KtoMod = KtoMod,
      toPDF = ToPDF,
      SortIndividuals = SortIndividuals,
      Labels = Labels,
      IndividualLabelCex = IndividualLabelCex
    )
  }

  invisible(colour_file_name)
}
