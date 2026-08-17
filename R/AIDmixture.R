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

.validate_scalar_inputs <- function(Q_file, Kseq, ColourPalette, lab.cex, KtoMod, ToPDF) {
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

.prepare_input_files <- function(Q_file, fam_file, Sort_file, Kseq) {
  if (is.null(fam_file)) {
    fam_file <- paste0(Q_file, "fam")
  }
  if (!is.character(fam_file) || length(fam_file) != 1L || is.na(fam_file) || !file.exists(fam_file)) {
    stop("Could not open fam file: ", fam_file, call. = FALSE)
  }

  fam <- fread(fam_file, select = 1L, header = FALSE, col.names = "IDs")
  if (nrow(fam) < 1L) {
    stop("The fam file is empty.", call. = FALSE)
  }
  fam$IDs <- as.character(fam$IDs)
  if (anyNA(fam$IDs) || any(!nzchar(fam$IDs))) {
    stop("The first column of the fam file contains missing or empty population IDs.", call. = FALSE)
  }

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
      fwrite(data.frame(IDs = unique(fam$IDs)), Sort_file, col.names = FALSE)
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

  missing_ids <- setdiff(unique(fam$IDs), sort_ids)
  if (length(missing_ids) > 0L) {
    stop(
      "The sort file is missing population IDs: ", paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }

  list(fam_file = fam_file, Sort_file = Sort_file)
}

#' @noRd
ColourFile <- function(colourADMX, Kseq) {
  max_k <- max(Kseq)

  if (!file.exists(colourADMX)) {
    col_matrix <- matrix(NA_integer_, max_k, max_k)
    for (r in seq_len(max_k)) {
      col_matrix[r, seq_len(r)] <- seq_len(r)
    }
    fwrite(as.data.frame(col_matrix), colourADMX, sep = "\t", col.names = FALSE)
  }

  col_file <- fread(colourADMX, header = FALSE, fill = TRUE)
  changed <- FALSE

  if (ncol(col_file) < max_k) {
    col_file <- as.data.frame(cbind(
      as.data.frame(col_file),
      matrix(NA_integer_, nrow(col_file), max_k - ncol(col_file))
    ))
    changed <- TRUE
  }

  if (nrow(col_file) < max_k) {
    old_n <- nrow(col_file)
    extra <- matrix(NA_integer_, max_k - old_n, max_k)
    for (r in seq_len(nrow(extra))) {
      k <- old_n + r
      extra[r, seq_len(k)] <- seq_len(k)
    }
    col_file <- as.data.frame(rbind(as.matrix(col_file), extra))
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
  KtoMod = 0L
) {
  if (toPDF) {
    pdf(paste0(Q_is, "ADMX.pdf"), paper = "a4r")
    on.exit(dev.off(), add = TRUE)
  }

  fam <- fread(Fam_is, select = 1L, header = FALSE, col.names = "IDs")
  fam$IDs <- as.character(fam$IDs)
  sort_table <- fread(Sort_is, header = FALSE, select = 1L)
  sort_ids <- as.character(sort_table[[1L]])

  n_plots <- length(Kseq)
  extra_panels <- if (!modding || toPDF) 1L else 2L
  layout(
    mat = matrix(seq_len(n_plots + extra_panels), ncol = 1L),
    heights = c(rep(1, n_plots), rep(0.5, extra_panels))
  )
  par(mai = c(0, 0, margin, 0))

  q_plot <- NULL
  fam_change <- NULL

  for (K in Kseq) {
    q <- fread(paste0(Q_is, K, ".Q"), header = FALSE)
    q_plot <- cbind(data.frame(IDs = fam$IDs), as.data.frame(q))
    q_plot <- q_plot[order(match(q_plot$IDs, sort_ids)), , drop = FALSE]

    fam_change <- 0
    if (nrow(q_plot) > 1L) {
      fam_change <- c(0, which(q_plot$IDs[-1L] != q_plot$IDs[-nrow(q_plot)]))
    }

    colour_indices <- as.integer(unlist(Col_is[K, seq_len(K), with = FALSE], use.names = FALSE))
    barplot(
      t(as.matrix(q_plot[, -1L, drop = FALSE])),
      width = 1,
      border = NA,
      col = colorPal[colour_indices],
      axes = FALSE,
      space = 0
    )
    text(0, 0.5, K, font = 2, pos = 2)
    abline(v = fam_change)
  }

  plot(NA, xlim = c(0, nrow(q_plot)), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  fam_unique <- data.frame(IDs = unique(q_plot$IDs), Start = fam_change)
  label_cex <- if (lab.cex == 1) n_plots^(-0.05) else lab.cex

  for (S in seq_len(nrow(fam_unique))) {
    next_start <- if (S != nrow(fam_unique)) fam_unique$Start[S + 1L] else nrow(q_plot)
    midpoint <- fam_unique$Start[S] + (next_start - fam_unique$Start[S]) / 2
    adjustment <- if (S %% 2L) 1 else 1.1
    text(
      midpoint, 1,
      adj = c(adjustment, adjustment),
      labels = fam_unique$IDs[S],
      srt = 60,
      cex = label_cex,
      xpd = TRUE
    )
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
#' persisted in an `ADMXcolors` file and can be changed interactively by
#' selecting two colour swatches.
#'
#' @param Q_file Character scalar giving the prefix shared by the ADMIXTURE
#'   files. For example, `"results/run."` expects `results/run.2.Q` for K = 2.
#' @param fam_file Path to a PLINK-style fam file. Only its first column is used
#'   as the population or group identifier. Defaults to `paste0(Q_file, "fam")`.
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
#'
#' @return Invisibly, the path to the persistent `ADMXcolors` mapping file.
#' @export
#'
#' @examples
#' \dontrun{
#' Admixture_ModPlot(
#'   Q_file = "results/run.",
#'   fam_file = "results/run.fam",
#'   Kseq = 2:6
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
  ToPDF = FALSE
) {
  checked <- .validate_scalar_inputs(Q_file, Kseq, ColourPalette, lab.cex, KtoMod, ToPDF)
  Kseq <- checked$Kseq
  KtoMod <- checked$KtoMod

  files <- .prepare_input_files(Q_file, fam_file, Sort_file, Kseq)
  fam_file <- files$fam_file
  Sort_file <- files$Sort_file

  if (is.null(ColourPalette)) {
    ColourPalette <- .default_colour_palette()
  }
  if (length(ColourPalette) < max(Kseq)) {
    stop("The colour palette must contain at least max(`Kseq`) colours.", call. = FALSE)
  }

  colour_file_name <- paste0(Q_file, "ADMXcolors")
  ColourFile(colour_file_name, Kseq)
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
    toPDF = ToPDF && !mod_plot
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
      toPDF = ToPDF
    )
  }

  invisible(colour_file_name)
}
