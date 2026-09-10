## Base-R build: no dplyr, no rlang, no factoextra. See BASE-R-CHANGES.md.
utils::globalVariables(c("PC1", "PC2", "k", "gap", "sd", "pct", "x", "y",
                         "contrib", "varname", "col", "label"))

NULL


#' Select variables for PCA clustering from a scenario string
#'
#' Translates a scenario string into a character vector of column names to use
#' as PCA input. The scenario string is composed of tokens separated by `_`,
#' where each token maps to one or more column names. For example, `"yba_hbf"`
#' selects yellowfin, bigeye, and albacore fractions plus hooks between floats.
#'
#' @details Available tokens:
#' \describe{
#'   \item{Species fractions}{`yft`, `bet`, `alb`, `skj`, `oth`,
#'     `yba` (yft + bet + alb), `sp` (yft + bet + alb + oth)}
#'   \item{Length}{`len` (mean), `lens` (sd), `lenAll` (mean + sd)}
#'   \item{Gear}{`hbf` (mean hooks between floats), `hbfC` (categorised)}
#'   \item{Location}{`lat`, `latabs`, `lat7.5`, `lon`, `bat`, `batC`}
#'   \item{Time}{`month` (sine/cosine encoding), `qua` (quarter)}
#'   \item{Other}{`CPUE`, `flag`, `fleet`}
#' }
#'
#' @param scenario Character string of underscore-separated tokens
#'   (e.g. `"yba_hbf"`, `"yft_lat_month"`).
#'
#' @return A character vector of column names to extract from the data before
#'   scaling and PCA.
#'
#' @examples
#' selectScenario("yba_hbf")
#' selectScenario("yft_lat_month")
#' selectScenario("sp_hbf_lon_lat")
#'
#' @family clustering utilities
#' @export
selectScenario <- function(scenario) {
	token_map <- list(
		yft    = "yft_fraction",
		bet    = "bet_fraction",
		alb    = "alb_fraction",
		skj    = "skj_fraction",
		oth    = "oth_fraction",
		yba    = c("yft_fraction", "bet_fraction", "alb_fraction"),
		sp     = c("yft_fraction", "bet_fraction", "alb_fraction", "oth_fraction"),
		len    = "mean_len",
		lens   = "sd_len",
		lenAll = c("mean_len", "sd_len"),
		hbf    = "mean_hbf",
		hbfC   = "hbf_cat",
		lat    = "latitude",
		latabs = "latAbs",
		lat7.5 = "lat7.5",
		lon    = "longitude",
		bat    = "bat",
		batC   = "bat_cat",
		month  = c("cosmonth", "sinmonth"),
		qua    = "quarter",
		CPUE   = "CPUE",
		flag   = "flag",
		fleet  = "fleet_cat"
	)

	components <- unlist(strsplit(scenario, "_"))
	unknown <- setdiff(components, names(token_map))
	if (length(unknown) > 0)
		warning("Unknown scenario components: ", paste(unknown, collapse = ", "))

	unique(unlist(token_map[intersect(components, names(token_map))]))
}


# Z-score scaling with zero-variance handling.
# If all values are identical, centres to zero instead of dividing by zero.
custom_scale <- function(x) {
	if (stats::var(x, na.rm = TRUE) == 0) {
		return(x - mean(x, na.rm = TRUE))
	} else {
		return(as.numeric(scale(x)))
	}
}


#' Scale variables prior to PCA
#'
#' Scales a dataframe of numeric variables using one of two methods. Z-score
#' is the default and most broadly appropriate. Arcsine square root is suited
#' to proportion data bounded in \[0, 1\].
#'
#' @param select_dat A numeric dataframe of variables to scale (typically the
#'   subset of columns returned by [selectScenario()]).
#' @param method Character; one of `"zscore"` (default) or `"arcsine"`.
#'   \describe{
#'     \item{`zscore`}{Subtracts mean and divides by SD per column. Zero-variance
#'       columns are centred only.}
#'     \item{`arcsine`}{Applies arcsine square root transformation. Appropriate
#'       for proportions in \[0, 1\].}
#'   }
#'
#' @return A dataframe of the same dimensions as `select_dat` with scaled values.
#'
#' @examples
#' df <- data.frame(yft_fraction = c(0.1, 0.5, 0.9),
#'                  bet_fraction = c(0.3, 0.2, 0.1))
#' pcaScale(df)
#' pcaScale(df, method = "arcsine")
#'
#' @family clustering utilities
#' @export
pcaScale <- function(select_dat, method = "zscore") {
	f <- switch(method,
				"zscore"  = custom_scale,
				"arcsine" = function(x) asin(sqrt(x)),
				stop("Unknown scaling method: '", method,
					 "'. Use 'zscore' or 'arcsine'."))
	select_dat <- as.data.frame(select_dat)
	bad <- !vapply(select_dat, is.numeric, logical(1))
	if (any(bad))
		stop("pcaScale() needs numeric columns; these are not: ",
			 paste(names(select_dat)[bad], collapse = ", "))

	scale_params <- if (method == "zscore") {
		data.frame(
			variable = names(select_dat),
			mean     = vapply(select_dat, mean, numeric(1), na.rm = TRUE),
			sd       = vapply(select_dat, function(x) {
				s <- stats::sd(x, na.rm = TRUE)
				if (is.na(s) || s == 0) 1 else s   # matches custom_scale()'s zero-variance branch
			}, numeric(1)),
			row.names = NULL
		)
	} else NULL   # arcsine has no linear params to save

	select_dat[] <- lapply(select_dat, f)
	attr(select_dat, "scale_params") <- scale_params
	select_dat
}


#' Run PCA and select components by variance threshold
#'
#' Runs PCA on pre-scaled data, retains the minimum number of principal
#' components needed to explain at least `variance_threshold` of total
#' variance, and optionally prints diagnostic plots.
#'
#' @param data_df A numeric dataframe of scaled variables (output of
#'   [pcaScale()]).
#' @param variance_threshold Numeric in (0, 1); retain enough PCs to explain
#'   at least this fraction of total variance. Default is 0.7.
#' @param print_it Logical; if `TRUE`, prints the scree plot, variable
#'   contribution plot, and PC1/PC2 scatter. Default is `TRUE`.
#' @param return_print Logical; if `TRUE`, returns the diagnostic plots as a
#'   list instead of the PCA result. Default is `FALSE`.
#'
#' @return If `return_print = FALSE` (default): a `prcomp` object with an
#'   additional element `no_var` giving the number of retained components
#'   (minimum 2).
#'
#'   If `return_print = TRUE`: a named list with elements `loadings`,
#'   `variances_plot`, `loading_plot`, and `scatter_plot`.
#'
#' @examples
#' \dontrun{
#' pca_res <- customPCA(pca_scaled, variance_threshold = 0.7, print_it = FALSE)
#' pca_scores <- pca_res$x[, 1:pca_res$no_var]
#' }
#'
#' @family clustering utilities
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_text geom_segment geom_col geom_line geom_path scale_colour_identity coord_equal arrow unit labs theme_minimal
#' @export
customPCA <- function(data_df, variance_threshold = 0.7, print_it = TRUE,
				  return_print = FALSE) {
	pca_res    <- stats::prcomp(data_df, scale. = FALSE)
	pca_scores <- as.data.frame(pca_res$x)
	loadings   <- as.data.frame(pca_res$rotation)

	if (print_it) print(loadings)

	eigenvalues    <- (pca_res$sdev)^2
	variance_frac  <- eigenvalues / sum(eigenvalues)
	cumulative_var <- cumsum(variance_frac)
	no_pca         <- min(which(cumulative_var >= variance_threshold),
						  length(cumulative_var))
	pca_res$no_var <- max(no_pca, 2)

	## --- scree plot: replaces factoextra::fviz_eig() -------------------------
	scree <- data.frame(k = seq_along(variance_frac),
						pct = variance_frac * 100)
	p1 <- ggplot(scree, aes(x = factor(k), y = pct, group = 1)) +
		geom_col(fill = "steelblue") +
		geom_point() +
		geom_line() +
		geom_hline(yintercept = 0) +
		labs(title = "Scree plot", x = "Principal component",
			 y = "% of explained variance") +
		theme_minimal()

	## --- variable plot: replaces factoextra::fviz_pca_var() ------------------
	## Correlation circle. `contrib` is the contribution of each variable to the
	## PC1-PC2 plane, eigenvalue-weighted, as factoextra defines it.
	e12 <- eigenvalues[1] + eigenvalues[2]
	vars <- data.frame(
		varname = rownames(loadings),
		PC1     = loadings$PC1 * sqrt(eigenvalues[1]),
		PC2     = loadings$PC2 * sqrt(eigenvalues[2]),
		contrib = (loadings$PC1^2 * eigenvalues[1] +
				   loadings$PC2^2 * eigenvalues[2]) / e12 * 100,
		stringsAsFactors = FALSE)
	## Colours are resolved here with grDevices::colorRampPalette() and fed
	## through scale_colour_identity(), rather than scale_colour_gradientn().
	## A continuous colour scale forces ggplot2 to load `farver` at plot
	## CONSTRUCTION time, which fails on installs where farver predates R 4.0.
	## Doing the interpolation ourselves keeps this function usable there.
	pal <- grDevices::colorRampPalette(c("#00AFBB", "#E7B800", "#FC4E07"))(100)
	rng <- range(vars$contrib)
	idx <- if (diff(rng) > 0)
		as.integer(round((vars$contrib - rng[1]) / diff(rng) * 99)) + 1L
	else rep(50L, nrow(vars))
	vars$col   <- pal[idx]
	vars$label <- sprintf("%s (%.0f%%)", vars$varname, vars$contrib)

	circ <- data.frame(x = cos(seq(0, 2 * pi, length.out = 200)),
					   y = sin(seq(0, 2 * pi, length.out = 200)))
	p2 <- ggplot(vars, aes(x = PC1, y = PC2)) +
		geom_path(data = circ, aes(x = x, y = y),
				  colour = "grey70", inherit.aes = FALSE) +
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

	p3 <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
		geom_point() +
		geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
		geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
		geom_text(data = loadings, aes(label = rownames(loadings)),
				  vjust = -0.5, hjust = -0.5) +
		geom_segment(data = loadings,
					 aes(x = 0, y = 0, xend = PC1, yend = PC2),
					 arrow = arrow(length = unit(0.2, "cm")), color = "red") +
		labs(title = "PCA scatter plot",
			 x = "Principal Component 1",
			 y = "Principal Component 2") +
		theme_minimal()

	if (print_it) {
		print(p1)
		print(p2)
		print(p3)
	}

	if (return_print) {
		return(list(loadings       = loadings,
					variances_plot = p1,
					loading_plot   = p2,
					scatter_plot   = p3))
	} else {
		return(pca_res)
	}
}


#' Custom ggplot2 theme for fisheries plots
#'
#' A clean `theme_bw()`-based theme with sensible defaults: legend at the
#' bottom, wide legend keys, and blank facet strip backgrounds.
#'
#' Reused as-is from `seafisheries::customTheme()` (R/plotting.R) — copied
#' here rather than taken as a dependency, since this package is meant to
#' stand alone from `seafisheries`.
#'
#' @param text_size Numeric; base text size in points. Default is 11.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   customTheme()
#'
#' @family plotting utilities
#' @importFrom ggplot2 theme_bw theme unit element_text element_blank
#' @export
customTheme <- function(text_size = 11) {
	theme_bw() +
		theme(
			legend.position    = "bottom",
			legend.key.width   = unit(2.3, "cm"),
			text               = element_text(size = text_size),
			strip.background   = element_blank()
		)
}


#' Run k-means clustering with gap statistic to determine optimal K
#'
#' Uses [cluster::clusGap()] with the Tibshirani 2001 SEmax criterion to select
#' the optimal number of clusters, then runs k-means with that K. Returns the
#' k-means result, the gap statistic table, and a diagnostic plot.
#'
#' @param data_df A numeric dataframe or matrix of PCA scores (rows =
#'   observations, columns = retained principal components).
#' @param max_k Integer; maximum number of clusters to evaluate. Default is 15.
#' @param random_set Integer; number of bootstrap samples for the gap
#'   statistic (`B` in [cluster::clusGap()]). Default is 100.
#' @param iter_max Integer; maximum iterations passed to [kmeans()]. Default
#'   is 10.
#' @param nstart Integer; number of random starts passed to [kmeans()].
#'   Default is 1.
#' @param d.power Numeric; power of the Euclidean distance used in the gap
#'   statistic computation. Default is 2 (squared Euclidean).
#' @param print_it Logical; if `TRUE`, prints the gap statistic plot. Default
#'   is `TRUE`.
#' @param ncores Integer; number of cores passed through to
#'   [fastClusGap()]'s `ncores` argument for the reference-dataset bootstrap
#'   (embarrassingly parallel via `parallel::mclapply()`). Default 1
#'   (serial). Unix-alikes only, on Windows, values > 1 warn and fall back
#'   to serial execution.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{`kmeans`}{A `kmeans` object fit at the optimal K.}
#'   \item{`gap_stat`}{A dataframe with columns `logW`, `E.logW`, `gap`, and
#'     `SE.sim` for k = 1 to `max_k`.}
#'   \item{`plot`}{A ggplot2 object showing the gap statistic vs K with error
#'     bars and a vertical line at the selected K.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- customKmeans(pca_scores, max_k = 10, random_set = 50, print_it = FALSE)
#' res$kmeans$cluster  # cluster assignments
#' res$gap_stat        # gap statistic table
#' res$plot            # diagnostic plot
#' }
#'
#' @family clustering utilities
#' @importFrom cluster clusGap maxSE
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_vline
#' @export
customKmeans <- function(data_df, max_k = 15, random_set = 100, iter_max = 10,
						 nstart = 1, d.power = 2, print_it = TRUE, ncores = 1) {

	## --- amended gap statistic: exact O(n), no dist() ----------------------
	## The W_k identity holds only for squared Euclidean distance, so fall
	## back to cluster::clusGap() for any other d.power.
	if (d.power == 2) {
		Tab <- fastClusGap(data_df, K.max = max_k, B = random_set,
						   nstart = nstart, iter.max = iter_max,
						   ncores = ncores)
	} else {
		Tab <- cluster::clusGap(x = data_df, FUNcluster = stats::kmeans,
								K.max = max_k, B = random_set,
								d.power = d.power, nstart = nstart,
								iter.max = iter_max)$Tab
		attr(Tab, "k") <- seq_len(nrow(Tab))
	}

	## selectK(), NOT maxSE(): maxSE returns a POSITION in the vector, which
	## equals k only when the table starts at k = 1.
	nc <- selectK(Tab)


    plot_gap <- data.frame(k   = as.factor(attr(Tab, "k")),
						   gap = Tab[, "gap"],
						   sd  = Tab[, "SE.sim"])

	p1 <- ggplot(plot_gap, aes(x = .data$k, y = .data$gap)) +
		geom_point() +
		geom_errorbar(aes(ymin = .data$gap - .data$sd, ymax = .data$gap + .data$sd)) +
		geom_vline(xintercept = which(attr(Tab, "k") == nc)) +
		customTheme()

	if (print_it) print(p1)

	cat("Best K =", nc, "\n")
	if (nc == max_k) cat("Selected K is max K\n")

	kmeans_res <- stats::kmeans(data_df, centers = nc,
						 iter.max = iter_max, nstart = nstart)

	return(list(
		kmeans   = kmeans_res,
		gap_stat = as.data.frame(Tab),
		plot     = p1
	))
}


#' Assign observations to the nearest cluster centroid
#'
#' For each row in `df`, computes the Euclidean distance to every centroid and
#' returns the index of the closest one. Used to apply cluster labels from a
#' training subset to the full dataset.
#'
#' @param df A dataframe of PCA scores (rows = observations). Column names
#'   must match those of `centroids`.
#' @param centroids A dataframe of cluster centroids (rows = clusters, columns
#'   = PCA dimensions), typically `as.data.frame(kmeans_res$centers)`.
#'
#' @return An integer vector of length `nrow(df)` giving the cluster index
#'   (1-based row of the nearest centroid) for each observation.
#'
#' @examples
#' \dontrun{
#' EC_clean$cluster <- assignClusters(as.data.frame(pca_full),
#'                                    as.data.frame(kmeans_full$centers))
#' }
#'
#' @family clustering utilities
#' @export
assignClusters <- function(df, centroids) {
	# Rewritten from the original mutate()/pull() chain: that version relied
	# on the magrittr "." placeholder inside mutate(cluster = apply(., 1, ...))
	# to refer to the whole piped dataframe -- a feature specific to magrittr's
	# %>%, with no base |> equivalent. Base R subsetting + apply() sidesteps
	# the issue and drops the dplyr dependency for this function entirely.
	sub_df <- df[, names(centroids), drop = FALSE]
	apply(sub_df, 1, function(point) {
		which.min(apply(centroids, 1, function(cent) {
			sqrt(sum((point - cent)^2))
		}))
	})
}


#' Summarise key statistics per cluster
#'
#' Computes per-cluster summaries of catch, effort, CPUE, and optionally
#' length, hooks between floats, and latitude depending on which columns are
#' present. Always includes catch composition and date range.
#'
#' @param df A dataframe with at minimum columns `cluster`, `yft_n`, `bet_n`,
#'   `alb_n`, `E`, and `ymd`. Optional columns `mean_len`, `mean_hbf`, and
#'   `latitude` are included in the summary when present.
#'
#' @return A dataframe with one row per cluster containing:
#' \describe{
#'   \item{`CPUE_min`, `CPUE_max`}{Range of daily mean CPUE (yft_n / E).}
#'   \item{`len_min`, `len_max`}{Range of mean length (if `mean_len` present).}
#'   \item{`HBF_min`, `HBF_max`}{Range of mean HBF (if `mean_hbf` present).}
#'   \item{`lat_min`, `lat_max`}{Range of mean latitude (if `latitude` present).}
#'   \item{`ymd_min`, `ymd_max`}{Date range of observations in the cluster.}
#'   \item{`yft_n`}{Total yellowfin catch.}
#'   \item{`total_n`}{Total catch across yft, bet, and alb.}
#'   \item{`yft_catch`}{Yellowfin as percentage of total catch across all clusters.}
#'   \item{`yft_frac`}{Yellowfin as percentage of total catch within the cluster.}
#' }
#'
#' @examples
#' \dontrun{
#' cluster_summary <- summaryClusters(EC_clustered)
#' }
#'
#' @family clustering utilities
#' @export
summaryClusters <- function(df) {
	has_len <- "mean_len"  %in% names(df)
	has_hbf <- "mean_hbf"  %in% names(df)
	has_lat <- "latitude"  %in% names(df)

	## --- helpers (base equivalents of group_by()/summarise()) ---------------
	## Group keys are built once; rowsum() does the per-(cluster, date) sums in
	## one pass, which is what makes this usable on the full 2e6-row dataset.
	## reorder = FALSE keeps groups in order of first appearance, so the key
	## columns can be recovered with !duplicated() instead of a string split.
	key   <- paste(as.character(df$cluster), format(df$ymd), sep = "\r")
	first <- !duplicated(key)
	d_cluster <- df$cluster[first]      # keeps factor / integer class
	d_ymd     <- df$ymd[first]          # keeps Date class

	sums <- rowsum(as.matrix(df[, c("yft_n", "bet_n", "alb_n", "E")]),
				   group = key, reorder = FALSE)
	d_yft   <- sums[, "yft_n"]
	d_tot   <- sums[, "yft_n"] + sums[, "bet_n"] + sums[, "alb_n"]
	d_CPUE  <- sums[, "yft_n"] / sums[, "E"]

	## dplyr's summarise() returns groups in sorted order; match that.
	lev <- sort(unique(d_cluster))
	g   <- factor(d_cluster, levels = lev)
	num <- function(x, f) vapply(split(x, g), f, numeric(1), USE.NAMES = FALSE)
	## split() keeps the Date class, so min/max return Dates; c() reassembles.
	dat <- function(x, f) do.call(c, lapply(split(x, g), f))

	result <- data.frame(
		cluster  = lev,
		CPUE_min = num(d_CPUE, min),
		CPUE_max = num(d_CPUE, max),
		ymd_min  = dat(d_ymd, min),
		ymd_max  = dat(d_ymd, max),
		yft_n    = num(d_yft, sum),
		total_n  = num(d_tot, sum),
		stringsAsFactors = FALSE)

	## Optional blocks: per-(cluster, date) mean, then range across dates.
	## Computed on the same `lev` ordering, so they are cbind-ed directly
	## rather than joined; match() guards the ordering regardless.
	range_by_cluster <- function(values) {
		m  <- rowsum(cbind(v = values, n = 1), group = key, reorder = FALSE)
		mu <- m[, "v"] / m[, "n"]
		gg <- factor(d_cluster, levels = lev)
		lo <- vapply(split(mu, gg), min, numeric(1))
		hi <- vapply(split(mu, gg), max, numeric(1))
		list(lo = lo[match(lev, names(lo))], hi = hi[match(lev, names(hi))])
	}

	if (has_len) {
		r <- range_by_cluster(df$mean_len)
		result$len_min <- unname(r$lo); result$len_max <- unname(r$hi)
	}
	if (has_hbf) {
		r <- range_by_cluster(df$mean_hbf)
		result$HBF_min <- unname(r$lo); result$HBF_max <- unname(r$hi)
	}
	if (has_lat) {
		## Preserved verbatim from the dplyr version: the per-date statistic is
		## mean(quantile(latitude, c(0.01, 0.99))), i.e. a midrange of the 1st
		## and 99th percentiles, not a mean latitude.
		mid <- vapply(split(df$latitude, factor(key, levels = unique(key))),
					  function(v) mean(stats::quantile(v, c(0.01, 0.99))),
					  numeric(1))
		gg <- factor(d_cluster, levels = lev)
		lo <- vapply(split(mid, gg), min, numeric(1))
		hi <- vapply(split(mid, gg), max, numeric(1))
		result$lat_min <- unname(lo[match(lev, names(lo))])
		result$lat_max <- unname(hi[match(lev, names(hi))])
	}

	result$yft_catch <- result$yft_n / sum(result$yft_n) * 100
	result$yft_frac  <- result$yft_n / result$total_n * 100
	result
}


#' Reorder clusters by descending yellowfin fraction
#'
#' Renumbers cluster labels so that cluster 1 has the highest proportion of
#' yellowfin in its catch, cluster 2 the second highest, and so on. Ensures
#' consistent, interpretable ordering across runs.
#'
#' @param df A dataframe with columns `cluster`, `yft_n`, and `total_n`.
#'
#' @return The input dataframe with `cluster` replaced by a factor ordered
#'   by descending yellowfin fraction.
#'
#' @examples
#' \dontrun{
#' EC_clean$cluster <- assignClusters(as.data.frame(pca_full),
#'                                    as.data.frame(kmeans_full$centers))
#' EC_clustered <- orderClusters(EC_clean)
#' }
#'
#' @family clustering utilities
#' @export
orderClusters <- function(df) {
	lev  <- sort(unique(df$cluster))
	g    <- factor(df$cluster, levels = lev)
	yft  <- vapply(split(df$yft_n,   g), sum, numeric(1))
	tot  <- vapply(split(df$total_n, g), sum, numeric(1))
	frac <- yft / tot
	## order() is stable, matching dplyr::arrange(desc(.))
	old_ids <- lev[order(frac, decreasing = TRUE)]
	new_ids <- seq_along(old_ids)
	df$cluster <- factor(new_ids[match(df$cluster, old_ids)], levels = new_ids)
	attr(df, "old_ids") <- old_ids
	df
}

#' Number of PCA components to retain from species count in a scenario
#'
#' Resolves `scenario` via [selectScenario()] and counts how many of the
#' resulting columns are species-fraction variables (`yft_fraction`,
#' `bet_fraction`, `alb_fraction`, `skj_fraction`, `oth_fraction`). Returns
#' that count minus one, floored at `min_components`.
#'
#' @param scenario Character string, as passed to [selectScenario()].
#' @param min_components Integer; minimum value returned regardless of the
#'   species count. Default 2.
#'
#' @return Integer number of components to retain.
#'
#' @examples
#' nSpeciesComponents("yba_lat")   # 3 species -> 2
#' nSpeciesComponents("sp_hbf")    # 4 species -> 3
#'
#' @family clustering utilities
#' @export
nSpeciesComponents <- function(scenario, min_components = 2) {
	species_cols <- c("yft_fraction", "bet_fraction", "alb_fraction",
					  "skj_fraction", "oth_fraction")
	resolved  <- selectScenario(scenario)
	n_species <- length(intersect(resolved, species_cols))
	n_out     <- n_species - 1L
	if (n_out < min_components) {
		message("nSpeciesComponents(): n_species - 1 = ", n_out,
				" is below min_components = ", min_components,
				"; using ", min_components, " instead.")
		n_out <- min_components
	}
	n_out
}

#' Recompute species fractions scoped to the species used in a scenario
#'
#' Overwrites `total_n` and the relevant `*_fraction` columns so that,
#' within each row, the fractions of the species present in `scenario` sum
#' to 1. Species not referenced by `scenario` are left untouched.
#'
#' @param df A dataframe with `*_n` count columns for each species
#'   (`yft_n`, `bet_n`, `alb_n`, `skj_n`, `oth_n`).
#' @param scenario Character string, as passed to [selectScenario()].
#'
#' @return `df` with `total_n` and the scenario's `*_fraction` columns
#'   recomputed.
#'
#' @family clustering utilities
#' @export
recalcSpeciesFractions <- function(df, scenario) {
	species_cols <- c("yft_fraction", "bet_fraction", "alb_fraction",
					  "skj_fraction", "oth_fraction")
	resolved  <- selectScenario(scenario)
	frac_cols <- intersect(resolved, species_cols)
	if (length(frac_cols) == 0) {
		message("recalcSpeciesFractions(): no species tokens in scenario '",
				scenario, "'; leaving fractions and total_n unchanged.")
		return(df)
	}
	n_cols <- sub("_fraction$", "_n", frac_cols)
	missing_n <- setdiff(n_cols, names(df))
	if (length(missing_n) > 0)
		stop("recalcSpeciesFractions(): missing count column(s): ",
			 paste(missing_n, collapse = ", "))

	df$total_n <- rowSums(df[, n_cols, drop = FALSE])
	for (i in seq_along(frac_cols)) {
		df[[frac_cols[i]]] <- ifelse(df$total_n > 0,
									 df[[n_cols[i]]] / df$total_n, 0)
	}
	df
}
