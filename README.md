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

The vignette walks through the clustering step by step on a small bundled
mock dataset:

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

## Basic example

Run the clustering, then make a PDF report of the diagnostics. This example
uses the bundled mock data and small settings so it runs in a few seconds.

``` r
library(longlineTypology)

mock <- system.file("extdata", "mockdata", "LL_mock.RDS", package = "longlineTypology")

# 1. Clustering: writes mock_outputs/cluster_diagnostics.RDS
diag <- runClustering(mock, scenario = "sp_hbf_latabs", outdir = "mock_outputs",
                      sample_no = 100, max_k = 10)

# 2. Report: reloads the data from the path stored in the diagnostics
clusterReport("mock_outputs/cluster_diagnostics.RDS", output_file = "mock_report.pdf")
```

Notes:

- `scenario` selects the clustering variables (see `?selectScenario`): e.g.
  `sp` = yellowfin, bigeye, albacore and other species fractions, `hbf` = hooks
  between floats, `latabs` = absolute latitude.
- The default `sample_no` (20000) is sized for real data; lower it for small
  datasets such as the mock data.
- Give `runClustering()` the **path** to the data file, so the report can find
  the data again. If the data were given as a data frame, or the file has
  moved, pass the data to the report: `clusterReport(diag, data = "new/path.RDS")`.
- The report needs `rmarkdown` and a LaTeX installation (e.g.
  `tinytex::install_tinytex()`); the variable importance section also needs
  `randomForest`.
- Choose the report sections and values of K with
  `clusterReport(diag, sections = c("pca", "gap"), k = 4)`.

The plots of the report are also available as functions (`plotPCA()`,
`plotGapSet()`, `plotVarImportance()`, `plotClusterVar()`), e.g.:

``` r
dat <- clusterData(diag)                       # the clustered data
plotGapSet(diag$gap_stat)
plotClusterVar(dat, diag$pca_select, clusterLabels(diag, k = 4))
```
