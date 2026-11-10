# longlineTypology

<!-- badges: start -->
[![R-CMD-check](https://github.com/PacificCommunity/longlineTypology/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PacificCommunity/longlineTypology/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

longlineTypology builds a typology of longline fishing practices: it groups
fisheries records into clusters of similar fishing strategy (species
composition, hooks between floats, fish length, area) with a PCA followed by
k-means, the number of clusters being chosen with the gap statistic. It is
designed for catch-effort and length-frequency data in the format produced by
the `seafisheries` package, though any data in that format can be used.

## Installation

Install the version from GitHub. Use `build_vignettes = TRUE` to
be able to open the vignette after installation:

``` r
# install.packages("remotes")
remotes::install_github("PacificCommunity/longlineTypology", build_vignettes = TRUE)
```

## Vignette

The vignette walks through the clustering step by step on a small simulated
dataset (`mockLonglineData()`):

``` r
browseVignettes("longlineTypology")                   # list, in the browser
vignette("clusters", package = "longlineTypology")    # open it directly
```

If R says the vignette is not found, the package was installed without
`build_vignettes = TRUE`. The source is also in the repository at
`vignettes/clusters.Rmd`; it can be rendered on your own data by changing its
parameters:

``` r
rmarkdown::render("vignettes/clusters.Rmd",
                  params = list(scenario  = "sp_hbf_latabs",
                                data_path = "path/to/LL_EC_LF.RDS",
                                outdir    = "cluster_outputs/"))
```

## Mock data

`mockLonglineData()` simulates operational longline data in the `seafisheries`
format, from four fishing strategies (tropical shallow, tropical deep,
subtropical albacore, temperate mixed). Column `mock_strategy` gives the
strategy behind each row, to compare with the clusters found.

``` r
library(longlineTypology)

mock    <- mockLonglineData()                    # small: 3 rows per strategy and month
mock_bp <- mockLonglineData(breakpoints = TRUE)  # 25 rows per strategy and month, with planted changes
attr(mock_bp, "events")                          # what was planted, where and when
```

With `breakpoints = TRUE`, six changes are planted, one for each label of the
breakpoint analysis (see `?mockLonglineData`):

| Change | Strategy | Start | Expected label |
|--------|----------|-------|----------------|
| CPUE drop in every strategy | all | 1995 | `common` |
| strong year class (small fish, more catch) | tropical deep | 2008 | `transient` in mean length (the CPUE increase is labelled `unexplained`) |
| mean length +15 cm over 6 years | subtropical albacore | 2001 | `gradual` |
| same flags fish 25 degrees further east, catch rate x 2.5 | temperate mixed | 2005 | `relocation` |
| CPUE drop to 35 % in one month | tropical shallow | 2010 | `unexplained` |
| three temperate flags switch to the albacore strategy | temperate mixed to subtropical albacore | 2014 | `switch` |

Cluster numbers follow the clustering (ordered by yellowfin catch), not the
order of the strategies. The analysis also finds a few changes that were not
planted: false alarms of the method on autocorrelated series.

## Basic example

Run the clustering on the small mock dataset, then make a PDF report of the
diagnostics. Small settings, so it runs in about a minute.

``` r
library(longlineTypology)

saveRDS(mockLonglineData(), "mock.RDS")

# 1. Clustering: writes mock_outputs/cluster_diagnostics.RDS
diag <- runClustering("mock.RDS", scenario = "sp_hbf_latabs", outdir = "mock_outputs",
                      sample_no = 500, max_k = 10)

# 2. Report: reloads the data from the path stored in the diagnostics
clusterReport("mock_outputs/cluster_diagnostics.RDS", k = 4, output_file = "mock_report.pdf")
```

## Breakpoint example

The same steps on the mock data with planted changes, then the breakpoint
report for K = 4: all cluster series with the changes found, and a detailed
section for each change labelled `switch`, `unexplained` or `ambiguous`.

``` r
saveRDS(mockLonglineData(breakpoints = TRUE), "mock_bp.RDS")

diag_bp <- runClustering("mock_bp.RDS", scenario = "sp_hbf_latabs", outdir = "mock_bp_outputs",
                         sample_no = 1000, kmeans_set = 5, max_k = 8, B = 10)

clusterReport(diag_bp, k = 4, output_file = "mock_bp_cluster_report.pdf")
breakpointReport(diag_bp, k = 4, output_file = "mock_bp_breakpoint_report.pdf")

# every change, with a detailed section for each
breakpointReport(diag_bp, k = 4, cards = "all", output_file = "mock_bp_all_changes.pdf")
```

## Notes

- `scenario` selects the clustering variables (see `?selectScenario`): e.g.
  `sp` = yellowfin, bigeye, albacore and other species fractions, `hbf` = hooks
  between floats, `latabs` = absolute latitude.
- The default `sample_no` (20000) is sized for real data; lower it for small
  datasets such as the mock data.
- Give `runClustering()` the **path** to the data file, so the reports can find
  the data again. If the data were given as a data frame, or the file has
  moved, pass the data to the report: `clusterReport(diag, data = "new/path.RDS")`.
- The reports need `rmarkdown` and a LaTeX installation (e.g.
  `tinytex::install_tinytex()`); the variable importance section also needs
  `randomForest`.
- Choose the report sections and values of K with
  `clusterReport(diag, sections = c("pca", "gap"), k = 4)`.
- `breakpointReport()` takes one K at a time; `manual =` adds an inspection of
  any cluster, variable and date, whether or not a change was found there.

The plots of the reports are also available as functions (`plotPCA()`,
`plotGapSet()`, `plotVarImportance()`, `plotClusterVar()`,
`plotClusterBreakpoints()`), e.g.:

``` r
dat <- clusterData(diag)                       # the clustered data
plotGapSet(diag$gap_stat)
plotClusterVar(dat, diag$pca_select, clusterLabels(diag, k = 4))
```
