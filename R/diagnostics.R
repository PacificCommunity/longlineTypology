## R/diagnostics.R -- building blocks of the clustering diagnostics report.
##
## Every function works from the diagnostics object saved by the clustering
## vignette (cluster_diagnostics.RDS) and the processed dataset used to fit
## it. Computation and plotting are separate: compute functions return data
## frames, plot functions return ggplot objects (or draw with base graphics
## for plotGapSet()), and nothing is printed or written to disk.
## clusterReport() assembles them into a PDF.

utils::globalVariables(c("cluster", "variable", "value", "metric", "m"))


# ---- internal helpers ------------------------------------------------------

## Diagnostics given as a list or as the path to the saved RDS.
.asDiagnostics <- function(diagnostics) {
	if (is.character(diagnostics)) {
		if (length(diagnostics) != 1 || !file.exists(diagnostics))
			stop("diagnostics: no file at '", diagnostics, "'")
		diagnostics <- readRDS(diagnostics)
	}
	need <- c("pca_select", "n_components", "scale_params", "pca_res",
			  "kmeans_full", "cluster_k")
	miss <- setdiff(need, names(diagnostics))
	if (length(miss))
		stop("diagnostics is missing element(s): ", paste(miss, collapse = ", "),
			 ". Rerun the clustering vignette.")
	diagnostics
}

## Nearest centre by squared Euclidean distance, vectorised. Same result as
## assignClusters() but fast enough for the full dataset.
.nearestCenter <- function(X, centers) {
	X  <- as.matrix(X)
	C  <- as.matrix(centers)[, colnames(X), drop = FALSE]
	d2 <- outer(rowSums(X^2), rowSums(C^2), "+") - 2 * X %*% t(C)
	max.col(-d2, ties.method = "first")
}

## Evaluate `expr` with a fixed seed, then restore the caller's RNG state.
.withSeed <- function(seed, expr) {
	if (is.null(seed)) return(expr)
	old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
		get(".Random.seed", envir = globalenv()) else NULL
	on.exit({
		if (is.null(old)) {
			if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
				rm(".Random.seed", envir = globalenv())
		} else assign(".Random.seed", old, envir = globalenv())
	})
	set.seed(seed)
	expr
}


# ---- rows, labels and scores ---------------------------------------------

#' Rows of the fit window
#'
#' Positions of the rows of `data` whose year falls inside the fit window
#' used to fit the scaling, the PCA and kmeans (`diagnostics$fit_years`).
#' Diagnostics files saved before `fit_years` existed were fitted on 1990
#' onwards, which is used as the default.
#'
#' @param diagnostics The diagnostics list saved by the clustering vignette,
#'   or the path to `cluster_diagnostics.RDS`.
#' @param data The processed dataset used for the clustering (same rows, same
#'   order), with a `Date` column `ymd`.
#'
#' @return Integer vector of row positions in `data`.
#'
#' @family diagnostics
#' @export
fitRows <- function(diagnostics, data) {
	diagnostics <- .asDiagnostics(diagnostics)
	if (!"ymd" %in% names(data)) stop("fitRows(): data needs a column 'ymd'")
	fy  <- if (is.null(diagnostics$fit_years)) c(1990, Inf) else diagnostics$fit_years
	yrs <- as.integer(format(data$ymd, "%Y"))
	which(yrs >= fy[1] & yrs <= fy[2])
}


#' Cluster labels for one K
#'
#' Returns the cluster of each row for a given K, from the `cluster_k<K>`
#' columns of the diagnostics. By default the labels are the ordered ones
#' (`diagnostics$cluster_reorder`), i.e. the labels used in the reports.
#'
#' @inheritParams fitRows
#' @param k Number of clusters (one of `names(diagnostics$kmeans_full)`).
#' @param reorder Logical; use the ordered labels (`TRUE`, default) or the
#'   raw kmeans labels.
#'
#' @return Integer vector, one label per row of the clustered dataset.
#'
#' @family diagnostics
#' @export
clusterLabels <- function(diagnostics, k, reorder = TRUE) {
	diagnostics <- .asDiagnostics(diagnostics)
	col <- paste0("cluster_k", k)
	if (!col %in% names(diagnostics$cluster_k))
		stop("No ", col, " in diagnostics; available K: ",
			 paste(names(diagnostics$kmeans_full), collapse = ", "))
	raw <- as.integer(diagnostics$cluster_k[[col]])
	if (!reorder) return(raw)
	ro <- diagnostics$cluster_reorder[[as.character(k)]]
	if (is.null(ro)) stop("No cluster_reorder for K = ", k,
						  "; use reorder = FALSE or rerun the clustering vignette.")
	if (!setequal(ro$old, unique(raw)))
		stop("cluster_reorder for K = ", k, " does not match the cluster_k", k, " labels.")
	as.integer(ro$new[match(raw, ro$old)])
}


#' Project the data on the saved PCA
#'
#' Scales the clustering variables with the saved `scale_params` (fit-window
#' means and sds, not the data's own) and projects every row on the saved
#' PCA. The scores are then checked against the saved cluster assignments:
#' assigning each row to the nearest saved centroid must give back
#' `cluster_k<K>` for every K. A mismatch means `data` is not the dataset the
#' clustering was fitted on (other rows, other order or other processing).
#'
#' @inheritParams fitRows
#'
#' @return A matrix of PC scores (one row per row of `data`, the
#'   `diagnostics$n_components` retained PCs).
#'
#' @family diagnostics
#' @export
pcaScores <- function(diagnostics, data) {
	diagnostics <- .asDiagnostics(diagnostics)
	if (nrow(data) != nrow(diagnostics$cluster_k))
		stop("data has ", nrow(data), " rows; the clustering had ",
			 nrow(diagnostics$cluster_k), ". Use the same processed dataset.")
	miss <- setdiff(diagnostics$pca_select, names(data))
	if (length(miss)) stop("data is missing clustering variable(s): ",
						   paste(miss, collapse = ", "))

	scaled <- applyScaleParams(data[, diagnostics$pca_select, drop = FALSE],
							   diagnostics$scale_params)
	scores <- stats::predict(diagnostics$pca_res, newdata = scaled)[,
		seq_len(diagnostics$n_components), drop = FALSE]

	for (kk in names(diagnostics$kmeans_full)) {
		mis <- mean(.nearestCenter(scores, diagnostics$kmeans_full[[kk]]$centers) !=
					diagnostics$cluster_k[[paste0("cluster_k", kk)]])
		if (mis > 0)
			warning("K = ", kk, ": the projected scores assign ", round(100 * mis, 3),
					"% of rows to another cluster than saved: data probably differs ",
					"from the dataset the clustering was fitted on.")
	}
	scores
}


# ---- PCA ------------------------------------------------------------------

#' PCA diagnostic plots
#'
#' Scree plot, correlation circle and score scatter of a fitted PCA, as in
#' [customPCA()] but from a saved PCA (`diagnostics$pca_res`), without
#' refitting it.
#'
#' @param pca_res A `prcomp` object, e.g. `diagnostics$pca_res`.
#' @param scores Optional matrix of scores with columns `PC1` and `PC2`, e.g.
#'   from [pcaScores()]. `NULL` (default) skips the scatter plot.
#' @param max_points Maximum number of points drawn in the scatter plot; a
#'   random subset is drawn above that, to keep PDF files small. `Inf` draws
#'   them all.
#' @param seed Seed of the random subset (the caller's RNG state is
#'   restored).
#'
#' @return A named list of ggplot objects: `scree`, `loadings` and, if
#'   `scores` is given, `scatter`.
#'
#' @family diagnostics
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_text geom_segment geom_col geom_line geom_path scale_colour_identity coord_equal arrow unit labs theme_minimal
#' @export
plotPCA <- function(pca_res, scores = NULL, max_points = 1e4, seed = 1) {
	loadings      <- as.data.frame(pca_res$rotation)
	eigenvalues   <- pca_res$sdev^2
	variance_frac <- eigenvalues / sum(eigenvalues)

	scree <- data.frame(k = seq_along(variance_frac), pct = variance_frac * 100)
	p_scree <- ggplot(scree, aes(x = factor(k), y = pct, group = 1)) +
		geom_col(fill = "steelblue") +
		geom_point() +
		geom_line() +
		geom_hline(yintercept = 0) +
		labs(title = "Scree plot", x = "Principal component",
			 y = "% of explained variance") +
		theme_minimal()

	## Correlation circle; contrib = eigenvalue-weighted contribution of each
	## variable to the PC1-PC2 plane. Colours interpolated by hand and passed
	## through scale_colour_identity(), as in customPCA().
	e12  <- eigenvalues[1] + eigenvalues[2]
	vars <- data.frame(
		varname = rownames(loadings),
		PC1     = loadings$PC1 * sqrt(eigenvalues[1]),
		PC2     = loadings$PC2 * sqrt(eigenvalues[2]),
		contrib = (loadings$PC1^2 * eigenvalues[1] + loadings$PC2^2 * eigenvalues[2]) / e12 * 100,
		stringsAsFactors = FALSE)
	pal <- grDevices::colorRampPalette(c("#00AFBB", "#E7B800", "#FC4E07"))(100)
	rng <- range(vars$contrib)
	idx <- if (diff(rng) > 0)
		as.integer(round((vars$contrib - rng[1]) / diff(rng) * 99)) + 1L
	else rep(50L, nrow(vars))
	vars$col   <- pal[idx]
	vars$label <- sprintf("%s (%.0f%%)", vars$varname, vars$contrib)

	circ <- data.frame(x = cos(seq(0, 2 * pi, length.out = 200)),
					   y = sin(seq(0, 2 * pi, length.out = 200)))
	p_load <- ggplot(vars, aes(x = PC1, y = PC2)) +
		geom_path(data = circ, aes(x = x, y = y), colour = "grey70", inherit.aes = FALSE) +
		geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
		geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
		geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = col),
					 arrow = arrow(length = unit(0.2, "cm"))) +
		geom_text(aes(label = label), vjust = -0.6, size = 3) +
		scale_colour_identity() +
		coord_equal() +
		labs(title = "Variables - PCA",
			 x = sprintf("Dim1 (%.1f%%)", variance_frac[1] * 100),
			 y = sprintf("Dim2 (%.1f%%)", variance_frac[2] * 100)) +
		theme_minimal()

	out <- list(scree = p_scree, loadings = p_load)

	if (!is.null(scores)) {
		if (!all(c("PC1", "PC2") %in% colnames(scores)))
			stop("plotPCA(): scores needs columns PC1 and PC2")
		sc <- as.data.frame(scores[, c("PC1", "PC2"), drop = FALSE])
		n_drawn <- nrow(sc)
		if (nrow(sc) > max_points) {
			sc <- sc[.withSeed(seed, sample.int(nrow(sc), max_points)), , drop = FALSE]
			n_drawn <- max_points
		}
		out$scatter <- ggplot(sc, aes(x = PC1, y = PC2)) +
			geom_point(alpha = if (n_drawn > 2000) 0.3 else 1) +
			geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
			geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
			geom_text(data = loadings, aes(label = rownames(loadings)),
					  vjust = -0.5, hjust = -0.5) +
			geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
						 arrow = arrow(length = unit(0.2, "cm")), color = "red") +
			labs(title = "PCA scatter plot",
				 subtitle = if (n_drawn < nrow(scores))
				 	sprintf("%d of %d rows drawn at random", n_drawn, nrow(scores)),
				 x = "Principal Component 1", y = "Principal Component 2") +
			theme_minimal()
	}
	out
}


# ---- variable importance ----------------------------------------------------

#' Subsamples of the replicate samples
#'
#' Draws, without replacement, `size` positions from each element of
#' `sample_indices` (the replicate samples of the gap statistic,
#' `diagnostics$sample_indices`, which are positions in the fit-window rows).
#' Used to run costly diagnostics (random forests, distance matrices) on
#' smaller samples.
#'
#' @param sample_indices List of integer vectors, e.g.
#'   `diagnostics$sample_indices`.
#' @param size Number of positions drawn from each element; capped at the
#'   element's length.
#' @param seed Seed (the caller's RNG state is restored). `NULL` uses the
#'   current RNG state.
#'
#' @return A list of integer vectors, same length as `sample_indices`.
#'
#' @family diagnostics
#' @export
subsampleReplicates <- function(sample_indices, size = 3000, seed = NULL) {
	.withSeed(seed, lapply(sample_indices, function(s)
		s[sample.int(length(s), min(size, length(s)))]))
}


#' Variable importance of each cluster
#'
#' For each cluster, a one-vs-rest random forest (cluster against all other
#' clusters) is fitted on each sample of rows, and the importance of each
#' variable is averaged over the samples. Requires the `randomForest`
#' package.
#'
#' Three metrics are returned:
#' \describe{
#'   \item{`MeanDecreaseAccuracy`}{drop in out-of-bag accuracy over all rows
#'     when the variable is permuted.}
#'   \item{`MeanDecreaseAccuracyCluster`}{same, among the rows of the cluster
#'     only, so not capped by the cluster's share of the rows.}
#'   \item{`MeanDecreaseGini`}{total decrease in node impurity.}
#' }
#'
#' @param data A dataframe with the clustering variables.
#' @param variables Names of the variables (e.g. `diagnostics$pca_select`).
#' @param cluster Cluster label of each row of `data`.
#' @param samples List of integer vectors, positions in `data` of the rows of
#'   each sample, e.g. from [subsampleReplicates()]. `NULL` uses all rows
#'   once.
#' @param seed Seed of the random forests (the caller's RNG state is
#'   restored). `NULL` uses the current RNG state.
#' @param ... Passed to [randomForest::randomForest()], e.g. `ntree`.
#'
#' @return A dataframe with columns `cluster`, `variable`, `metric`, `m`
#'   (mean over samples) and `s` (sd over samples; `NA` with one sample).
#'
#' @family diagnostics
#' @export
varImportance <- function(data, variables, cluster, samples = NULL, seed = NULL, ...) {
	if (!requireNamespace("randomForest", quietly = TRUE))
		stop("varImportance() needs the 'randomForest' package.")
	if (length(cluster) != nrow(data))
		stop("cluster has ", length(cluster), " values; data has ", nrow(data), " rows")
	miss <- setdiff(variables, names(data))
	if (length(miss)) stop("data is missing variable(s): ", paste(miss, collapse = ", "))
	if (is.null(samples)) samples <- list(seq_len(nrow(data)))

	x  <- as.data.frame(data)[, variables, drop = FALSE]
	cl <- sort(unique(cluster))
	if (length(cl) < 2) stop("varImportance() needs at least 2 clusters")

	res <- .withSeed(seed, do.call(rbind, lapply(cl, function(c) {
		y <- factor(ifelse(cluster == c, "1", "0"), levels = c("0", "1"))
		do.call(rbind, lapply(seq_along(samples), function(i) {
			r   <- samples[[i]]
			imp <- randomForest::randomForest(x = x[r, , drop = FALSE], y = y[r],
											  importance = TRUE, ...)$importance
			metrics <- c(MeanDecreaseAccuracy        = "MeanDecreaseAccuracy",
						 MeanDecreaseAccuracyCluster = "1",
						 MeanDecreaseGini            = "MeanDecreaseGini")
			do.call(rbind, lapply(names(metrics), function(mt) data.frame(
				cluster = c, variable = rownames(imp), metric = mt,
				value = imp[, metrics[[mt]]], sample = i, row.names = NULL)))
		}))
	})))

	key <- interaction(res$cluster, res$variable, res$metric, drop = TRUE, lex.order = TRUE)
	agg <- do.call(rbind, lapply(split(res, key), function(d) data.frame(
		cluster = d$cluster[1], variable = d$variable[1], metric = d$metric[1],
		m = mean(d$value), s = if (nrow(d) > 1) stats::sd(d$value) else NA_real_)))
	rownames(agg) <- NULL
	agg
}


#' Plot variable importance per cluster
#'
#' One heat map per metric (clusters by variables), each with its own colour
#' scale since the metrics have different units.
#'
#' @param vi Output of [varImportance()].
#' @param metrics Metrics to plot, among those of `vi`.
#' @param variable_order Optional order of the variables (bottom to top);
#'   variables of `vi` not listed are added after it.
#'
#' @return A named list of ggplot objects, one per metric.
#'
#' @family diagnostics
#' @importFrom ggplot2 ggplot aes geom_tile labs theme_bw theme
#' @export
plotVarImportance <- function(vi, metrics = unique(vi$metric), variable_order = NULL) {
	vars   <- unique(vi$variable)
	levels <- c(intersect(variable_order, vars), setdiff(vars, variable_order))
	vi$variable <- factor(vi$variable, levels = levels)
	vi$cluster  <- factor(vi$cluster, levels = sort(unique(vi$cluster)))
	stats::setNames(lapply(metrics, function(mt) {
		ggplot(vi[vi$metric == mt, , drop = FALSE], aes(cluster, variable, fill = m)) +
			geom_tile() +
			labs(x = "Cluster", y = NULL, fill = "Importance", title = mt) +
			theme_bw() +
			theme(legend.position = "bottom", legend.direction = "horizontal")
	}), metrics)
}


# ---- variable distributions -------------------------------------------------

#' Distribution of the clustering variables by cluster
#'
#' Box plots (interquartile range and median) of each variable per cluster,
#' with every variable rescaled to \[0, 1\] so they share one axis. By
#' default a variable is rescaled over its range; `limits` sets another
#' range for some variables (values outside are clamped), e.g. to stop a few
#' extreme values from squashing the boxes.
#'
#' @param data A dataframe with the variables.
#' @param variables Names of the variables (e.g. `diagnostics$pca_select`).
#' @param cluster Cluster label of each row of `data`.
#' @param limits Optional named list of `c(min, max)` per variable.
#'
#' @return A ggplot object.
#'
#' @family diagnostics
#' @importFrom ggplot2 ggplot aes geom_boxplot labs scale_y_continuous guides guide_legend
#' @export
plotClusterVar <- function(data, variables, cluster, limits = NULL) {
	if (length(cluster) != nrow(data))
		stop("cluster has ", length(cluster), " values; data has ", nrow(data), " rows")
	miss <- setdiff(variables, names(data))
	if (length(miss)) stop("data is missing variable(s): ", paste(miss, collapse = ", "))

	long <- do.call(rbind, lapply(variables, function(v) {
		x   <- data[[v]]
		rng <- if (!is.null(limits[[v]])) limits[[v]] else range(x, na.rm = TRUE)
		x   <- pmin(pmax(x, rng[1]), rng[2])
		x   <- if (diff(rng) == 0) rep(0, length(x)) else (x - rng[1]) / diff(rng)
		data.frame(cluster = cluster, variable = v, value = x)
	}))
	long$variable <- factor(long$variable, levels = variables)

	ggplot(long, aes(x = factor(cluster), y = value, fill = variable)) +
		geom_boxplot(outlier.shape = NA, coef = 0) +
		labs(x = "Cluster", y = "Scaled value", fill = NULL) +
		scale_y_continuous(limits = c(0, 1)) +
		customTheme() +
		guides(fill = guide_legend(nrow = 3))
}


# ---- report -----------------------------------------------------------------

#' Clustering diagnostics report
#'
#' Renders a PDF report with the selected diagnostics, from the diagnostics
#' object of [runClustering()] (or of the clustering vignette) and the data
#' the clustering was fitted on. Needs `rmarkdown` and a LaTeX installation;
#' the `"varimp"` section also needs `randomForest`.
#'
#' Sections, in report order:
#' \describe{
#'   \item{`"pca"`}{loadings table, scree plot, correlation circle and score
#'     scatter of all rows ([plotPCA()]).}
#'   \item{`"gap"`}{gap profiles of the replicates ([plotGapSet()]) and the
#'     number of replicates selecting each K.}
#'   \item{`"centroids"`}{per K, variability of the replicate centroids
#'     around the full-fit centroids ([centroidVariability()]).}
#'   \item{`"jaccard"`}{per K, Jaccard stability of each cluster across the
#'     replicate fits, fit years ([jaccardStability()]).}
#'   \item{`"silhouette"`}{per K, silhouette widths on subsamples of the
#'     replicate samples, fit years ([silhouetteStats()]); with
#'     `null_silhouette = TRUE`, compared with a null reference.}
#'   \item{`"varimp"`}{per K, variable importance of each cluster
#'     ([varImportance()]) on the same subsamples as the silhouette.}
#'   \item{`"vardist"`}{per K, distribution of the clustering variables by
#'     cluster, all rows ([plotClusterVar()]).}
#'   \item{`"temporal"`}{per K, each cluster's share of each year's rows,
#'     fit years ([temporalStability()]).}
#' }
#'
#' The subsamples of the silhouette and the variable importance are drawn
#' once per K from the replicate samples of the gap statistic
#' (`diagnostics$sample_indices`, fit years) with [subsampleReplicates()],
#' seed `seed + K`, so both sections use the same rows.
#'
#' @inheritParams fitRows
#' @param data `NULL` (default) to reload the data from
#'   `diagnostics$data_path`, or the data given to [runClustering()] (data
#'   frame or path). Prepared with [clusterData()] in both cases.
#' @param sections Sections to include (see Details).
#' @param k Values of K reported in the per-K sections. `NULL` (default)
#'   reports every K of the diagnostics.
#' @param output_file Path of the PDF.
#' @param subsample_size Rows per subsample for the silhouette and the
#'   variable importance (distance matrices grow as the square of it).
#' @param seed Seed of the subsamples, random forests and null silhouette
#'   (`seed + K` for K).
#' @param null_silhouette Compare the silhouette with a null reference
#'   (slow: `n_null` kmeans fits per K). Default `FALSE`.
#' @param n_null Number of null datasets when `null_silhouette = TRUE`.
#' @param share_threshold Passed to [temporalStability()] as `threshold`.
#' @param var_limits Passed to [plotClusterVar()] as `limits`.
#' @param title Title of the report.
#' @param quiet Passed to [rmarkdown::render()].
#'
#' @return The path of the PDF, invisibly.
#'
#' @examples
#' \dontrun{
#' clusterReport("cluster_outputs/cluster_diagnostics.RDS",
#'               sections = c("pca", "gap", "varimp"),
#'               output_file = "report_sp_hbf_latabs.pdf")
#' }
#'
#' @family diagnostics
#' @export
clusterReport <- function(diagnostics, data = NULL,
						  sections = c("pca", "gap", "centroids", "jaccard", "silhouette",
									   "varimp", "vardist", "temporal"),
						  k = NULL, output_file = "cluster_report.pdf",
						  subsample_size = 3000, seed = 123,
						  null_silhouette = FALSE, n_null = 100, share_threshold = 0.01,
						  var_limits = NULL, title = "Clustering diagnostics report",
						  quiet = TRUE) {
	if (!requireNamespace("rmarkdown", quietly = TRUE))
		stop("clusterReport() needs the 'rmarkdown' package.")
	sections <- match.arg(sections, c("pca", "gap", "centroids", "jaccard", "silhouette",
									  "varimp", "vardist", "temporal"), several.ok = TRUE)
	if ("varimp" %in% sections && !requireNamespace("randomForest", quietly = TRUE))
		stop("The 'varimp' section needs the 'randomForest' package.")

	diagnostics <- .asDiagnostics(diagnostics)
	ks <- names(diagnostics$kmeans_full)
	if (!is.null(k)) {
		bad <- setdiff(as.character(k), ks)
		if (length(bad)) stop("K not in diagnostics: ", paste(bad, collapse = ", "),
							  " (available: ", paste(ks, collapse = ", "), ")")
		ks <- as.character(k)
	}

	## Data rebuilt and checked up front, so a wrong dataset fails before
	## rendering rather than halfway through the PDF.
	data     <- clusterData(diagnostics, data)
	scores   <- pcaScores(diagnostics, data)
	fit_rows <- fitRows(diagnostics, data)

	template <- system.file("reports", "cluster_report.Rmd", package = "longlineTypology")
	if (template == "") stop("Report template not found in the installed package.")
	## Rendered from a temporary copy: the installed package may not be writable.
	work <- tempfile("cluster_report_")
	dir.create(work)
	file.copy(list.files(dirname(template), full.names = TRUE), work)

	env <- new.env(parent = globalenv())
	env$diagnostics  <- diagnostics
	env$data         <- data
	env$scores       <- scores
	env$fit_rows     <- fit_rows
	env$sections     <- sections
	env$ks           <- ks
	env$subsample_size  <- subsample_size
	env$seed            <- seed
	env$null_silhouette <- null_silhouette
	env$n_null          <- n_null
	env$share_threshold <- share_threshold
	env$var_limits      <- var_limits

	output_file <- normalizePath(output_file, mustWork = FALSE)
	out <- rmarkdown::render(file.path(work, "cluster_report.Rmd"),
							 output_file = basename(output_file),
							 output_dir = dirname(output_file),
							 intermediates_dir = work,
							 params = list(title = title),
							 envir = env, quiet = quiet)
	invisible(out)
}
