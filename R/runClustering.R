## R/runClustering.R -- the clustering pipeline of the clusters vignette as
## functions: prepareClusterData() (row filtering and species fractions),
## runClustering() (scaling, PCA, gap statistic replicates, kmeans per K,
## diagnostics object) and clusterData() (rebuild the clustered dataset
## from the diagnostics).


## Data given as a data frame or as the path to an RDS file.
.readClusterInput <- function(data) {
	if (is.character(data)) {
		if (length(data) != 1 || !file.exists(data))
			stop("data: no file at '", data, "'")
		return(list(df = readRDS(data),
					path = normalizePath(data),
					md5 = unname(tools::md5sum(data))))
	}
	if (!is.data.frame(data)) stop("data must be a data frame or the path to an RDS file")
	list(df = data, path = NULL, md5 = NULL)
}


#' Prepare a dataset for clustering
#'
#' The data processing done before the clustering, in one place so that the
#' clustering and the diagnostics use exactly the same rows:
#' \enumerate{
#'   \item recompute the species fractions of the scenario
#'     ([recalcSpeciesFractions()]);
#'   \item if the scenario uses HBF, keep rows with `hbf_outlier == "valid"`;
#'   \item if the scenario uses length, drop rows with
#'     `len_source == "imputed"`.
#' }
#' Running it on data it already prepared returns the same rows.
#'
#' @param data A data frame in the `seafisheries` output format.
#' @param scenario Scenario string, as passed to [selectScenario()].
#'
#' @return `data` with the fractions recomputed and the invalid rows
#'   removed. Rows with missing values in the clustering variables are not
#'   removed: [runClustering()] stops on them.
#'
#' @family clustering pipeline
#' @export
prepareClusterData <- function(data, scenario) {
	data <- recalcSpeciesFractions(as.data.frame(data), scenario)
	vars <- selectScenario(scenario)
	if (any(c("mean_hbf", "hbf_cat") %in% vars)) {
		if (!"hbf_outlier" %in% names(data))
			stop("Scenario '", scenario, "' uses HBF: data needs a column 'hbf_outlier'")
		data <- data[data$hbf_outlier %in% "valid", , drop = FALSE]
	}
	if ("mean_len" %in% vars) {
		if (!"len_source" %in% names(data))
			stop("Scenario '", scenario, "' uses length: data needs a column 'len_source'")
		data <- data[!data$len_source %in% "imputed", , drop = FALSE]
	}
	rownames(data) <- NULL
	data
}


#' Run the clustering
#'
#' Runs the clustering of the clusters vignette and returns the diagnostics
#' object the diagnostic functions and [clusterReport()] read:
#' \enumerate{
#'   \item prepares the data ([prepareClusterData()]);
#'   \item fits the z-score scaling and the PCA on the rows of the fit years,
#'     and projects every row on them;
#'   \item runs `kmeans_set` replicates of the gap statistic, each on
#'     `sample_no` rows drawn with replacement from the fit years
#'     ([customKmeans()]);
#'   \item for every K selected by at least one replicate, fits kmeans on
#'     the fit years and assigns every row to the nearest centroid;
#'   \item orders the clusters of each K ([orderClusters()]).
#' }
#'
#' @param data A data frame in the `seafisheries` output format, or the
#'   path to an RDS file of one. With a path, the path and the file's MD5
#'   checksum are stored in the diagnostics, so [clusterData()] and
#'   [clusterReport()] can reload the data without being given it.
#' @param scenario Scenario string, as passed to [selectScenario()].
#' @param outdir Folder where `cluster_diagnostics.RDS` is written (created
#'   if needed). `NULL` (default) writes nothing.
#' @param n_pc Number of principal components kept. `0` (default) uses
#'   `pca_method`.
#' @param pca_method `"n_components"` (number of species in the scenario
#'   minus one, [nSpeciesComponents()]) or `"threshold"` (enough components
#'   to explain `pca_variance_threshold` of the variance).
#' @param pca_variance_threshold Used when `pca_method = "threshold"`.
#' @param fit_years Years used to fit the scaling, the PCA and kmeans,
#'   `c(first, last)`. Rows outside are only assigned to the clusters.
#' @param sample_no Rows per gap statistic replicate.
#' @param kmeans_set Number of gap statistic replicates.
#' @param max_k Largest K tested by the gap statistic.
#' @param iter_max,nstart Passed to [stats::kmeans()].
#' @param B Reference datasets per gap statistic replicate.
#' @param ncores Cores of the gap statistic (Unix only).
#' @param order_by Criterion of [orderClusters()], or `NULL` to keep the
#'   kmeans labels.
#' @param seed Seed of the run (the caller's RNG state is restored).
#' @param verbose Print progress messages.
#'
#' @return The diagnostics list (invisibly when `outdir` is given). Besides
#'   the fitted objects, it holds the run settings (`settings`), the data
#'   source (`data_path`, `data_md5`, `NULL` when `data` was a data frame)
#'   and the package version.
#'
#' @examples
#' \dontrun{
#' mock <- tempfile(fileext = ".RDS")
#' saveRDS(mockLonglineData(), mock)
#' diag <- runClustering(mock, "sp_hbf_latabs", outdir = "mock_outputs",
#'                       sample_no = 100, max_k = 10)
#' clusterReport(diag, output_file = "mock_report.pdf")
#' }
#'
#' @family clustering pipeline
#' @export
runClustering <- function(data, scenario, outdir = NULL,
						  n_pc = 0, pca_method = c("n_components", "threshold"),
						  pca_variance_threshold = 0.7, fit_years = c(1990, 2022),
						  sample_no = 20000, kmeans_set = 15, max_k = 20,
						  iter_max = 1e6, nstart = 10, B = 25, ncores = 1,
						  order_by = "yft_w", seed = 123, verbose = TRUE) {
	pca_method <- match.arg(pca_method)
	say <- function(...) if (verbose) message(...)
	settings <- list(scenario = scenario, n_pc = n_pc, pca_method = pca_method,
					 pca_variance_threshold = pca_variance_threshold,
					 fit_years = fit_years, sample_no = sample_no,
					 kmeans_set = kmeans_set, max_k = max_k, iter_max = iter_max,
					 nstart = nstart, B = B, order_by = order_by, seed = seed)

	input <- .readClusterInput(data)
	df    <- prepareClusterData(input$df, scenario)
	vars  <- selectScenario(scenario)
	miss  <- setdiff(c(vars, "ymd"), names(df))
	if (length(miss)) stop("data is missing column(s): ", paste(miss, collapse = ", "))
	say("Variables: ", paste(vars, collapse = ", "), "; ", nrow(df), " valid rows")

	## Fit rows, scaling and PCA
	yrs     <- as.integer(format(df$ymd, "%Y"))
	idx_fit <- which(yrs >= fit_years[1] & yrs <= fit_years[2])
	if (length(idx_fit) < 2) stop("Fewer than 2 rows in the fit years")
	select_dat   <- df[, vars, drop = FALSE]
	scale_params <- attr(pcaScale(select_dat[idx_fit, , drop = FALSE], method = "zscore"),
						 "scale_params")
	pca_scaled   <- applyScaleParams(select_dat, scale_params)
	pca_res      <- customPCA(pca_scaled[idx_fit, , drop = FALSE],
							  variance_threshold = pca_variance_threshold, print_it = FALSE)
	max_pc <- ncol(pca_res$x)

	if (n_pc > 0) {
		if (n_pc != round(n_pc) || n_pc < 2) stop("n_pc must be an integer >= 2, got ", n_pc)
		if (n_pc > max_pc) {
			warning("n_pc = ", n_pc, " exceeds the ", max_pc, " PCs available; using ", max_pc)
			n_pc <- max_pc
		}
		pca_res$no_var <- as.integer(n_pc)
	} else if (pca_method == "n_components") {
		pca_res$no_var <- nSpeciesComponents(scenario, min_components = 2)
	}
	n_comp <- pca_res$no_var
	say("PCA: ", n_comp, " of ", max_pc, " components kept")

	pca_full <- stats::predict(pca_res, newdata = pca_scaled)[, seq_len(n_comp), drop = FALSE]
	bad <- !stats::complete.cases(pca_full)
	if (any(bad))
		stop(sum(bad), " rows have missing values in ", paste(vars, collapse = ", "),
			 "; filter them before the clustering")
	pca_fit <- pca_full[idx_fit, , drop = FALSE]

	run <- .withSeed(seed, {
		## Gap statistic replicates
		sample_indices <- lapply(seq_len(kmeans_set), function(i)
			sample(seq_len(nrow(pca_fit)), size = sample_no, replace = TRUE))
		kmeans_res <- vector("list", kmeans_set)
		gap_stat   <- vector("list", kmeans_set)
		for (i in seq_len(kmeans_set)) {
			say("Gap statistic replicate ", i, "/", kmeans_set)
			## customKmeans() always cat()s "Best K"; silenced unless verbose
			fit_rep <- function() customKmeans(pca_fit[sample_indices[[i]], , drop = FALSE],
											   max_k = max_k, random_set = B, iter_max = iter_max,
											   nstart = nstart, d.power = 2, ncores = ncores,
											   print_it = FALSE)
			if (verbose) res <- fit_rep() else utils::capture.output(res <- fit_rep())
			kmeans_res[[i]] <- res$kmeans
			gap_stat[[i]]   <- res$gap_stat
		}
		n_sel   <- vapply(kmeans_res, function(x) length(unique(x$cluster)), integer(1))
		ft      <- table(n_sel)
		freq_df <- data.frame(n_clusters = as.numeric(names(ft)), freq = as.numeric(ft))

		## kmeans per selected K, every row assigned to the nearest centroid
		kmeans_full <- list(); cluster_reorder <- list()
		cluster_k   <- data.frame(row.names = seq_len(nrow(df)))
		for (k in sort(unique(n_sel))) {
			say("kmeans on the fit years, K = ", k)
			km  <- stats::kmeans(pca_fit, centers = k, iter.max = iter_max, nstart = nstart)
			raw <- .nearestCenter(pca_full, km$centers)
			kmeans_full[[as.character(k)]] <- km
			cluster_k[[paste0("cluster_k", k)]] <- as.integer(raw)
			old_ids <- if (is.null(order_by)) seq_len(k) else
				attr(orderClusters(cbind(df[, setdiff(names(df), "cluster"), drop = FALSE],
										 cluster = raw), by = order_by), "old_ids")
			cluster_reorder[[as.character(k)]] <- data.frame(old = old_ids,
															 new = seq_along(old_ids))
		}
		list(sample_indices = sample_indices, kmeans_res = kmeans_res, gap_stat = gap_stat,
			 freq_df = freq_df, kmeans_full = kmeans_full, cluster_k = cluster_k,
			 cluster_reorder = cluster_reorder)
	})

	pca_res$x <- NULL   # rebuilt from the data with pcaScores()
	diagnostics <- c(list(
		scenario = scenario, pca_method = pca_method,
		pca_variance_threshold = if (pca_method == "threshold") pca_variance_threshold else NA,
		pca_select = vars, n_components = n_comp, n_pc_requested = settings$n_pc,
		n_pc_available = max_pc,
		train_nrow = nrow(df), train_nrow_1990 = length(idx_fit),
		train_date_range = range(df$ymd), fit_years = fit_years,
		scale_params = scale_params, pca_res = pca_res),
		run,
		list(data_path = input$path, data_md5 = input$md5, settings = settings,
			 package_version = as.character(utils::packageVersion("longlineTypology"))))

	if (!is.null(outdir)) {
		dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
		f <- file.path(outdir, "cluster_diagnostics.RDS")
		saveRDS(diagnostics, f)
		say("Saved ", f)
		return(invisible(diagnostics))
	}
	diagnostics
}


#' Rebuild the clustered dataset
#'
#' Returns the dataset the clustering was fitted on, prepared as in
#' [runClustering()] ([prepareClusterData()]), with rows in the same order
#' as the saved cluster assignments.
#'
#' @param diagnostics The diagnostics list from [runClustering()], or the
#'   path to `cluster_diagnostics.RDS`.
#' @param data The data given to [runClustering()] (data frame or path), or
#'   the already prepared data. `NULL` (default) reads the file at
#'   `diagnostics$data_path`, and warns if its checksum changed since the
#'   clustering.
#'
#' @return The prepared data frame. Stops if its number of rows differs from
#'   the clustering's.
#'
#' @family clustering pipeline
#' @export
clusterData <- function(diagnostics, data = NULL) {
	diagnostics <- .asDiagnostics(diagnostics)
	if (is.null(data)) {
		if (is.null(diagnostics$data_path))
			stop("diagnostics has no data_path (clustering run on a data frame, or by an ",
				 "older version): give the data.")
		if (!file.exists(diagnostics$data_path))
			stop("Data file of the clustering not found: ", diagnostics$data_path,
				 ". Give the data, or its new path.")
		data <- diagnostics$data_path
	}
	input <- .readClusterInput(data)
	if (!is.null(input$md5) && !is.null(diagnostics$data_md5) &&
		input$md5 != diagnostics$data_md5)
		warning("The data file changed since the clustering (different MD5 checksum): ",
				input$path)
	scenario <- diagnostics$scenario
	df <- prepareClusterData(input$df, scenario)
	if (nrow(df) != nrow(diagnostics$cluster_k))
		stop("The prepared data have ", nrow(df), " rows; the clustering had ",
			 nrow(diagnostics$cluster_k), ". This is not the clustered dataset.")
	df
}
