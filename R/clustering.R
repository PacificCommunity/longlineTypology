#' @importFrom rlang .data
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
	if (var(x, na.rm = TRUE) == 0) {
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
	switch(method,
		   "zscore"  = select_dat |>
			   dplyr::mutate(dplyr::across(dplyr::everything(), custom_scale)),
		   "arcsine" = select_dat |>
			   dplyr::mutate(dplyr::across(dplyr::everything(), ~asin(sqrt(.x)))),
		   stop("Unknown scaling method: '", method, "'. Use 'zscore' or 'arcsine'.")
	)
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
#' @importFrom factoextra fviz_eig fviz_pca_var
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_text geom_segment arrow unit labs theme_minimal
#' @export
customPCA <- function(data_df, variance_threshold = 0.7, print_it = TRUE,
				  return_print = FALSE) {
	pca_res    <- prcomp(data_df, scale. = FALSE)
	pca_scores <- as.data.frame(pca_res$x)
	loadings   <- as.data.frame(pca_res$rotation)

	if (print_it) print(loadings)

	eigenvalues    <- (pca_res$sdev)^2
	variance_frac  <- eigenvalues / sum(eigenvalues)
	cumulative_var <- cumsum(variance_frac)
	no_pca         <- min(which(cumulative_var >= variance_threshold),
						  length(cumulative_var))
	pca_res$no_var <- max(no_pca, 2)

	p1 <- fviz_eig(pca_res)
	p2 <- fviz_pca_var(pca_res,
					   col.var      = "contrib",
					   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
					   repel        = TRUE)
	p3 <- ggplot(pca_scores, aes(x = .data$PC1, y = .data$PC2)) +
		geom_point() +
		geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
		geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
		geom_text(data = loadings, aes(label = rownames(loadings)),
				  vjust = -0.5, hjust = -0.5) +
		geom_segment(data = loadings,
					 aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2),
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
						 nstart = 1, d.power = 2, print_it = TRUE) {

	gap_stat <- cluster::clusGap(x = data_df, FUNcluster = kmeans,
								 K.max = max_k, B = random_set,
								 d.power = d.power, nstart = nstart,
								 iter.max = iter_max)

	nc <- cluster::maxSE(f    = gap_stat$Tab[, "gap"],
						 SE.f = gap_stat$Tab[, "SE.sim"],
						 method    = "Tibs2001SEmax",
						 SE.factor = 1)

	plot_gap <- data.frame(k   = as.factor(1:max_k),
						   gap = gap_stat$Tab[, 3],
						   sd  = gap_stat$Tab[, 4])

	p1 <- ggplot(plot_gap, aes(x = .data$k, y = .data$gap)) +
		geom_point() +
		geom_errorbar(aes(ymin = .data$gap - .data$sd, ymax = .data$gap + .data$sd)) +
		geom_vline(xintercept = nc) +
		customTheme()

	if (print_it) print(p1)

	cat("Best K =", nc, "\n")
	if (nc == max_k) cat("Selected K is max K\n")

	kmeans_res <- kmeans(data_df, centers = nc,
						 iter.max = iter_max, nstart = nstart)

	return(list(
		kmeans   = kmeans_res,
		gap_stat = as.data.frame(gap_stat$Tab),
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

	# Core per-date aggregation
	daily <- df |>
		dplyr::group_by(.data$cluster, .data$ymd) |>
		dplyr::summarise(
			yft_n   = sum(.data$yft_n),
			E       = sum(.data$E),
			total_n = sum(.data$yft_n) + sum(.data$bet_n) + sum(.data$alb_n),
			.groups = "drop"
		) |>
		dplyr::mutate(CPUE = .data$yft_n / .data$E)

	# Core cluster summary
	result <- daily |>
		dplyr::group_by(.data$cluster) |>
		dplyr::summarise(
			CPUE_min = min(.data$CPUE),
			CPUE_max = max(.data$CPUE),
			ymd_min  = min(.data$ymd),
			ymd_max  = max(.data$ymd),
			yft_n    = sum(.data$yft_n),
			total_n  = sum(.data$total_n),
			.groups  = "drop"
		)

	# Optional: length
	if (has_len) {
		len_summary <- df |>
			dplyr::group_by(.data$cluster, .data$ymd) |>
			dplyr::summarise(mean_len = mean(.data$mean_len), .groups = "drop") |>
			dplyr::group_by(.data$cluster) |>
			dplyr::summarise(len_min = min(.data$mean_len), len_max = max(.data$mean_len),
					  .groups = "drop")
		result <- dplyr::left_join(result, len_summary, by = "cluster")
	}

	# Optional: hooks between floats
	if (has_hbf) {
		hbf_summary <- df |>
			dplyr::group_by(.data$cluster, .data$ymd) |>
			dplyr::summarise(mean_hbf = mean(.data$mean_hbf), .groups = "drop") |>
			dplyr::group_by(.data$cluster) |>
			dplyr::summarise(HBF_min = min(.data$mean_hbf), HBF_max = max(.data$mean_hbf),
					  .groups = "drop")
		result <- dplyr::left_join(result, hbf_summary, by = "cluster")
	}

	# Optional: latitude
	if (has_lat) {
		lat_summary <- df |>
			dplyr::group_by(.data$cluster, .data$ymd) |>
			dplyr::summarise(mean_lat = mean(stats::quantile(.data$latitude, c(0.01, 0.99))),
					  .groups = "drop") |>
			dplyr::group_by(.data$cluster) |>
			dplyr::summarise(lat_min = min(.data$mean_lat), lat_max = max(.data$mean_lat),
					  .groups = "drop")
		result <- dplyr::left_join(result, lat_summary, by = "cluster")
	}

	result |>
		dplyr::mutate(
			yft_catch = .data$yft_n / sum(.data$yft_n) * 100,
			yft_frac  = .data$yft_n / .data$total_n * 100
		)
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
#' EC_clustered <- EC_clean |>
#'   dplyr::mutate(cluster = assignClusters(as.data.frame(pca_full),
#'                                          as.data.frame(kmeans_full$centers))) |>
#'   orderClusters()
#' }
#'
#' @family clustering utilities
#' @export
orderClusters <- function(df) {
	order <- df |>
		dplyr::group_by(.data$cluster) |>
		dplyr::summarise(yft_n   = sum(.data$yft_n),
				  total_n = sum(.data$total_n),
				  .groups = "drop") |>
		dplyr::mutate(yft_frac  = .data$yft_n / .data$total_n) |>
		dplyr::arrange(dplyr::desc(.data$yft_frac)) |>
		dplyr::mutate(new_clust = dplyr::row_number()) |>
		dplyr::select(.data$cluster, .data$new_clust)

	df |>
		dplyr::left_join(order, by = "cluster") |>
		dplyr::mutate(cluster = factor(.data$new_clust, levels = order$new_clust)) |>
		dplyr::select(-.data$new_clust)
}
