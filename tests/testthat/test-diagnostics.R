## Small clustering of the mock data (few replicates, small samples).
mockDiagnostics <- function() {
	mock <- system.file("extdata", "mockdata", "LL_mock.RDS", package = "longlineTypology")
	diag <- runClustering(mock, "sp_hbf_latabs", sample_no = 100, kmeans_set = 3,
						  max_k = 8, B = 10, nstart = 5, verbose = FALSE)
	list(diagnostics = diag, data = clusterData(diag))
}

test_that("pcaScores reproduces the saved assignments and flags other data", {
	m  <- mockDiagnostics()
	sc <- expect_silent(pcaScores(m$diagnostics, m$data))
	expect_equal(dim(sc), c(nrow(m$data), as.integer(m$diagnostics$n_components)))
	shuffled <- m$data[rev(seq_len(nrow(m$data))), ]
	w <- character(0)
	withCallingHandlers(pcaScores(m$diagnostics, shuffled), warning = function(cnd) {
		w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
	})
	expect_true(length(w) > 0 && all(grepl("another cluster", w)))
	expect_error(pcaScores(m$diagnostics, m$data[-1, ]), "rows")
})

test_that("fitRows uses fit_years", {
	m <- mockDiagnostics()
	fr <- fitRows(m$diagnostics, m$data)
	fy <- m$diagnostics$fit_years
	yr <- as.integer(format(m$data$ymd, "%Y"))
	expect_equal(fr, which(yr >= fy[1] & yr <= fy[2]))
})

test_that("clusterLabels applies cluster_reorder", {
	m  <- mockDiagnostics()
	k  <- names(m$diagnostics$kmeans_full)[1]
	ro <- m$diagnostics$cluster_reorder[[k]]
	raw <- clusterLabels(m$diagnostics, k, reorder = FALSE)
	expect_equal(clusterLabels(m$diagnostics, k), ro$new[match(raw, ro$old)])
	expect_error(clusterLabels(m$diagnostics, 99), "available K")
})

test_that("plotPCA returns ggplots and subsamples large score sets", {
	m <- mockDiagnostics()
	p <- plotPCA(m$diagnostics$pca_res)
	expect_named(p, c("scree", "loadings"))
	sc <- pcaScores(m$diagnostics, m$data)
	p  <- plotPCA(m$diagnostics$pca_res, sc, max_points = 50)
	expect_s3_class(p$scatter, "ggplot")
	expect_equal(nrow(p$scatter$data), 50)
})

test_that("subsampleReplicates caps size and restores the RNG", {
	set.seed(5); a <- runif(1)
	set.seed(5); s <- subsampleReplicates(list(1:10, 1:200), size = 50, seed = 1); b <- runif(1)
	expect_equal(a, b)
	expect_equal(lengths(s), c(10L, 50L))
	expect_identical(s, subsampleReplicates(list(1:10, 1:200), size = 50, seed = 1))
})

test_that("varImportance returns one row per cluster, variable and metric", {
	skip_if_not_installed("randomForest")
	m  <- mockDiagnostics()
	fr <- fitRows(m$diagnostics, m$data)
	k  <- names(m$diagnostics$kmeans_full)[1]
	cl <- clusterLabels(m$diagnostics, k)[fr]
	vi <- varImportance(m$data[fr, ], m$diagnostics$pca_select, cl,
						samples = m$diagnostics$sample_indices, seed = 1, ntree = 50)
	expect_equal(nrow(vi), length(unique(cl)) * length(m$diagnostics$pca_select) * 3)
	expect_setequal(vi$metric, c("MeanDecreaseAccuracy", "MeanDecreaseAccuracyCluster",
								 "MeanDecreaseGini"))
	expect_identical(vi, varImportance(m$data[fr, ], m$diagnostics$pca_select, cl,
									   samples = m$diagnostics$sample_indices, seed = 1, ntree = 50))
	p <- plotVarImportance(vi)
	expect_length(p, 3)
	expect_error(varImportance(m$data[fr, ], "nope", cl), "missing")
})

test_that("plotClusterVar rescales to [0, 1] with optional limits", {
	m  <- mockDiagnostics()
	cl <- clusterLabels(m$diagnostics, names(m$diagnostics$kmeans_full)[1])
	p  <- plotClusterVar(m$data, m$diagnostics$pca_select, cl,
						 limits = list(mean_hbf = c(5, 15)))
	expect_s3_class(p, "ggplot")
	expect_true(all(p$data$value >= 0 & p$data$value <= 1, na.rm = TRUE))
})

test_that("clusterReport renders a PDF from the diagnostics alone", {
	skip_on_cran()
	skip_if_not_installed("rmarkdown")
	skip_if_not_installed("randomForest")
	skip_if(!rmarkdown::pandoc_available() || Sys.which("pdflatex") == "", "no pandoc / LaTeX")
	m   <- mockDiagnostics()
	out <- tempfile(fileext = ".pdf")
	clusterReport(m$diagnostics, k = names(m$diagnostics$kmeans_full)[1],
				  output_file = out, subsample_size = 50, null_silhouette = TRUE, n_null = 3)
	expect_true(file.exists(out))
})

test_that("runClustering stores the data source and restores the RNG", {
	mock <- system.file("extdata", "mockdata", "LL_mock.RDS", package = "longlineTypology")
	set.seed(3); a <- runif(1)
	set.seed(3)
	d <- runClustering(mock, "sp_hbf_latabs", sample_no = 50, kmeans_set = 2, max_k = 5,
					   B = 5, nstart = 2, verbose = FALSE)
	expect_equal(runif(1), a)
	expect_equal(d$data_path, normalizePath(mock))
	expect_equal(d$data_md5, unname(tools::md5sum(mock)))
	expect_equal(nrow(d$cluster_k), nrow(clusterData(d)))

	df <- runClustering(readRDS(mock), "sp_hbf_latabs", sample_no = 50, kmeans_set = 2,
						max_k = 5, B = 5, nstart = 2, verbose = FALSE)
	expect_null(df$data_path)
	expect_identical(df$cluster_k, d$cluster_k)   # same seed, same result
	expect_error(clusterData(df), "no data_path")
	expect_equal(nrow(clusterData(df, readRDS(mock))), nrow(df$cluster_k))
})

test_that("prepareClusterData is idempotent", {
	mock <- readRDS(system.file("extdata", "mockdata", "LL_mock.RDS", package = "longlineTypology"))
	once <- prepareClusterData(mock, "sp_hbf_latabs")
	expect_identical(prepareClusterData(once, "sp_hbf_latabs"), once)
	expect_true(all(once$hbf_outlier == "valid"))
})
