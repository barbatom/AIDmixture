# AIDmixture

AIDmixture is a small R package for plotting [ADMIXTURE](https://dalexander.github.io/admixture/) Q-matrix output across multiple values of K. It keeps population ordering and cluster colours consistent between plots and provides an interactive colour-swap mode.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("barbatom/AIDmixture")
```

## Expected input

For a prefix such as `results/run.`, AIDmixture expects one ADMIXTURE Q file per requested K:

```text
results/run.2.Q
results/run.3.Q
results/run.4.Q
...
```

The `fam_file` must contain the same number of rows as every Q file. Only its first column is used; that column is interpreted as the population or grouping identifier.

A sort file is optional. It is a one-column text file listing population IDs in the order they should appear. If `Sort_file` is omitted, AIDmixture uses `paste0(Q_file, "sort")`; if that file does not exist, it is created from the population order in the fam file.

## Basic plotting

```r
library(AIDmixture)

Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Kseq = 2:6
)
```

A persistent colour mapping is stored in `results/run.ADMXcolors`. Re-running the function reuses that mapping so cluster colours remain stable.

## Reorder populations

Create a one-column file such as:

```text
Population_C
Population_A
Population_B
```

Then pass it explicitly:

```r
Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Sort_file = "results/population-order.txt",
  Kseq = 2:6
)
```

Every population ID present in the fam file must be represented in the sort file.

## Change colours interactively

Set `KtoMod` to one of the K values being plotted:

```r
Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Kseq = 2:6,
  KtoMod = 4
)
```

A colour strip is added below the plot. Click two swatches to swap their assignments for that K. The change is written to the `ADMXcolors` file and is reused in later plots.

## Custom palette

Supply any R-compatible colour vector with at least `max(Kseq)` entries:

```r
my_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")

Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Kseq = 2:6,
  ColourPalette = my_palette
)
```

## PDF output

```r
Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Kseq = 2:6,
  ToPDF = TRUE
)
```

This writes `results/run.ADMX.pdf`.

## Validation

AIDmixture checks common input errors before plotting, including missing files, Q/fam row-count mismatches, incorrect Q column counts, incomplete population sort files, invalid K values, and palettes that are too short.

## Development

Tests use `testthat`. The GitHub Actions workflow runs `R CMD check` on pull requests and pushes to `master`.

```r
install.packages(c("devtools", "testthat"))
devtools::test()
devtools::check()
```

## License

GPL-3.
