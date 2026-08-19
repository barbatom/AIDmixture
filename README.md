# AIDmixture

AIDmixture is a small R package for plotting [ADMIXTURE](https://dalexander.github.io/admixture/) Q-matrix output across multiple values of K. It keeps population ordering and ancestry-component colours consistent between plots and provides an interactive colour-swap mode.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("barbatom/AIDmixture")
```

## Try it with bundled toy data

The package includes a small deterministic ADMIXTURE example with 15 individuals from three toy populations and Q matrices for K = 2, 3, and 4. Copy the files to a temporary writable directory and plot them directly:

```r
library(AIDmixture)

toy_dir <- system.file("extdata", package = "AIDmixture")
work <- file.path(tempdir(), "AIDmixture-toy")
unlink(work, recursive = TRUE)
dir.create(work)

invisible(file.copy(
  list.files(toy_dir, pattern = "^toy[.]", full.names = TRUE),
  work
))

Admixture_ModPlot(
  Q_file = file.path(work, "toy."),
  fam_file = file.path(work, "toy.fam"),
  Kseq = 2:4,
  SortIndividuals = TRUE
)
```

This produces a complete example plot immediately after installation. The copied files also let you inspect the expected `.fam` and `.Q` formats without supplying your own ADMIXTURE run first.

## Expected input

For a prefix such as `results/run.`, AIDmixture expects one ADMIXTURE Q file per requested K:

```text
results/run.2.Q
results/run.3.Q
results/run.4.Q
...
```

The `fam_file` must contain the same number of rows as every Q file. Its first column is interpreted as the population or grouping identifier. When individual labels are requested, the second column is interpreted as the individual identifier, matching the usual PLINK `.fam` convention.

A sort file is optional. It is a one-column text file listing population IDs in the order they should appear. If `Sort_file` is omitted, AIDmixture uses `paste0(Q_file, "sort")`; if that file does not exist, it is created from the population order in the fam file.

## Basic plotting with your own data

```r
library(AIDmixture)

Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Kseq = 2:6
)
```

A persistent colour mapping is stored in `results/run.ADMXcolors`. Re-running the function reuses that mapping so cluster colours remain stable.

## Keep ancestry colours consistent across K

`AutoMatchColours = TRUE` is the default. When an `ADMXcolors` row has not yet been generated, AIDmixture compares the individual-level Q profiles at lower and higher K values and assigns colours so the same ancestry component keeps its previous colour whenever possible.

For consecutive K values, AIDmixture explicitly evaluates every possible one-component split from K to K + 1. It merges each candidate pair temporarily, scores component similarity with cosine similarity, and uses a maximum-weight assignment to identify the best lineage. The child most similar to the parental component keeps the parental colour; the other child receives a new colour.

For the built-in palette, new colours are chosen by greedy maximin separation in CIE Lab space: among unused colours, AIDmixture selects the one whose nearest already-used colour is as far away as possible perceptually. Exact RGB duplicates in the built-in palette are skipped automatically. Existing saved palette indices are not changed. Custom palettes with their own length retain their explicit user-supplied ordering.

For gaps larger than one K, existing components are matched directly and unmatched components receive new colours.

```r
Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Kseq = 2:6,
  AutoMatchColours = TRUE
)
```

Existing `ADMXcolors` rows are never overwritten, including rows changed with `KtoMod`. This means manual colour choices remain authoritative and are inherited when a later K is added. To regenerate the automatic mapping from scratch, delete the corresponding `ADMXcolors` file before plotting again.

Set `AutoMatchColours = FALSE` to use the original column-based mapping (`1:K`) for newly created rows.

## Sort individuals within populations

Set `SortIndividuals = TRUE` to order individuals within each population by ancestry. AIDmixture uses the highest K in `Kseq` as the reference. For each population, it identifies the ancestry component with the highest mean Q value and orders individuals from highest to lowest membership in that component. The resulting individual order is then reused for every plotted K, so individuals remain aligned across panels.

```r
Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Kseq = 2:6,
  SortIndividuals = TRUE
)
```

The default is `SortIndividuals = FALSE`, which preserves the original row order within each population.

## Choose population and individual labels

`Labels = "Population"` is the default and draws only population IDs below the ancestry panels. The plotting code explicitly disables `barplot()` axis names, so data-frame row numbers or row names are never shown accidentally.

The available modes are:

- `Labels = "Population"`: population IDs only.
- `Labels = "Individual"`: individual IDs from column 2 of the `.fam` file only.
- `Labels = "Both"`: individual IDs and population IDs in two separate label panels, so the two label types cannot overlap.
- `Labels = "None"`: no labels below the ADMIXTURE panels.

Use `lab.cex` to control population-label size and `IndividualLabelCex` to control individual-label size.

```r
Admixture_ModPlot(
  Q_file = "results/run.",
  fam_file = "results/run.fam",
  Kseq = 2:6,
  Labels = "Both",
  IndividualLabelCex = 0.6
)
```

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

A colour strip is added below the plot. Click two swatches to swap their assignments for that K. The change is written to the `ADMXcolors` file and is reused in later plots. Automatic matching does not overwrite this manual row.

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

AIDmixture checks common input errors before plotting, including missing files, Q/fam row-count mismatches, incorrect Q column counts, incomplete population sort files, invalid K values, invalid label modes, missing individual IDs when requested, non-logical flags, and palettes that are too short.

## Development

Tests use `testthat`. GitHub Actions runs `R CMD check` on pushes to `master` and pull requests targeting `master`.

```r
install.packages(c("devtools", "testthat"))
devtools::test()
devtools::check()
```

## License

GPL-3.
