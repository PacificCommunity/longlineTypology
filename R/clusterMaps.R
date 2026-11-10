## R/clusterMaps.R -- cluster colours and maps: dominant cluster per cell,
## share of each cluster per cell, and maps of CPUE, mean length and mean
## HBF per cluster.

utils::globalVariables(c("lon_cell", "lat_cell", "share", "long", "lat", "group", "band", "value"))


#' Colours of the clusters
#'
#' Colourblind-safe colours (Okabe-Ito), ordered so that neighbouring
#' cluster numbers get the most distinct colours. Above 10 clusters, evenly
#' spaced hues are used instead.
#'
#' @param k Number of clusters.
#'
#' @return A character vector of `k` colours, named `1` to `k`.
#'
#' @examples
#' clusterColours(4)
#'
#' @family cluster maps
#' @export
clusterColours <- function(k) {
	pal <- c("#0072B2", "#882255", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7",
			 "#F0E442", "#999999", "#000000")
	cols <- if (k <= length(pal)) pal[seq_len(k)] else
		grDevices::hcl(h = seq(15, 375, length.out = k + 1)[seq_len(k)], c = 100, l = 65)
	stats::setNames(cols, seq_len(k))
}


# ---- internal helpers ------------------------------------------------------------

.squish <- function(x, range = c(0, 1), only.finite = TRUE) {
	ok <- is.finite(x)
	x[ok] <- pmin(pmax(x[ok], range[1]), range[2])
	x
}

.gridCell <- function(x, size) floor(x / size) * size + size / 2

## Coastlines (maps package, longitudes 0-360), or nothing if not installed
.landLayer <- function() {
	if (!requireNamespace("maps", quietly = TRUE)) return(NULL)
	ggplot2::geom_polygon(data = ggplot2::map_data("world2"), ggplot2::aes(long, lat, group = group),
						  fill = "grey80", inherit.aes = FALSE)
}

## Map extent: given limits, or the extent of the cells
.mapLims <- function(lim, x, size) if (!is.null(lim)) lim else range(x, na.rm = TRUE) + c(-1, 1) * size

.checkMapData <- function(data, cluster, need) {
	if (length(cluster) != nrow(data))
		stop("cluster has ", length(cluster), " values; data has ", nrow(data), " rows")
	miss <- setdiff(need, names(data))
	if (length(miss)) stop("data is missing column(s): ", paste(miss, collapse = ", "))
}

## Share of each cluster in the total of `weight` per cell
.cellShare <- function(data, cluster, weight, size) {
	lon <- .gridCell(data$longitude %% 360, size)
	lat <- .gridCell(data$latitude, size)
	w   <- data[[weight]]
	ok  <- !is.na(lon) & !is.na(lat) & !is.na(w)
	key <- paste(lon, lat, cluster, sep = "\r")[ok]
	tot <- rowsum(w[ok], key, reorder = FALSE)
	parts <- do.call(rbind, strsplit(rownames(tot), "\r", fixed = TRUE))
	out <- data.frame(lon_cell = as.numeric(parts[, 1]), lat_cell = as.numeric(parts[, 2]),
					  cluster = parts[, 3], w = as.numeric(tot), stringsAsFactors = FALSE)
	cell_tot <- stats::ave(out$w, out$lon_cell, out$lat_cell, FUN = sum)
	out$share <- out$w / cell_tot
	out[is.finite(out$share), , drop = FALSE]
}

## Rows used for each variable, and its value per group:
## CPUE = sum(yft_n) / sum(E); mean_len weighted by yft_n; mean_hbf row mean
.mapRows <- function(data, variable) {
	ok <- switch(variable,
		CPUE     = !is.na(data$yft_n) & !is.na(data$E),
		mean_len = !is.na(data$mean_len) & !is.na(data$yft_n) & data$yft_n > 0 &
			(if ("len_source" %in% names(data)) !data$len_source %in% "imputed" else TRUE),
		mean_hbf = !is.na(data$mean_hbf) &
			(if ("hbf_outlier" %in% names(data)) data$hbf_outlier %in% "valid" else TRUE),
		stop("variable must be CPUE, mean_len or mean_hbf"))
	ok & !is.na(data$longitude) & !is.na(data$latitude)
}

.groupValue <- function(data, variable, key) {
	num <- switch(variable,
		CPUE     = data$yft_n,
		mean_len = data$mean_len * data$yft_n,
		mean_hbf = data$mean_hbf)
	den <- switch(variable,
		CPUE     = data$E,
		mean_len = data$yft_n,
		mean_hbf = rep(1, nrow(data)))
	s <- rowsum(cbind(num, den), key, reorder = FALSE)
	stats::setNames(s[, 1] / s[, 2], rownames(s))
}

.splitKey <- function(keys, names) {
	parts <- do.call(rbind, strsplit(keys, "\r", fixed = TRUE))
	out <- as.data.frame(parts, stringsAsFactors = FALSE)
	names(out) <- names
	out
}


# ---- dominance -------------------------------------------------------------------

#' Map of the dominant cluster per cell
#'
#' For each grid cell, the cluster with the largest share of the cell's
#' total of `weight` (effort by default), all rows pooled.
#'
#' @param data A data frame with `longitude`, `latitude` and the `weight`
#'   column.
#' @param cluster Cluster label of each row of `data` (1 to K).
#' @param weight Column summed per cell and cluster. Default `"E"` (effort).
#' @param cell_size Cell size (degrees). Default 5.
#' @param xlim,ylim Map extent (longitude 0-360, latitude). `NULL`: extent of
#'   the cells.
#'
#' @return A ggplot object.
#'
#' @family cluster maps
#' @export
plotDominance <- function(data, cluster, weight = "E", cell_size = 5, xlim = NULL, ylim = NULL) {
	.checkMapData(data, cluster, c("longitude", "latitude", weight))
	k  <- max(as.integer(cluster), na.rm = TRUE)
	sh <- .cellShare(data, cluster, weight, cell_size)
	sh <- sh[order(sh$lon_cell, sh$lat_cell, -sh$share), , drop = FALSE]
	dom <- sh[!duplicated(sh[, c("lon_cell", "lat_cell")]), , drop = FALSE]
	dom$cluster <- factor(dom$cluster, levels = seq_len(k))
	ggplot2::ggplot(dom) + .landLayer() +
		ggplot2::geom_tile(ggplot2::aes(lon_cell, lat_cell, fill = cluster),
						   width = cell_size, height = cell_size) +
		ggplot2::scale_fill_manual(values = clusterColours(k), drop = FALSE) +
		ggplot2::coord_cartesian(xlim = .mapLims(xlim, dom$lon_cell, cell_size),
								 ylim = .mapLims(ylim, dom$lat_cell, cell_size), expand = FALSE) +
		ggplot2::labs(x = "Longitude", y = "Latitude", fill = "Dominant cluster") +
		customTheme() +
		ggplot2::guides(fill = ggplot2::guide_legend(title.position = "top", nrow = if (k > 5) 2 else 1)) +
		ggplot2::theme(legend.key.width = ggplot2::unit(0.5, "cm"))
}


#' Map of the share of each cluster per cell
#'
#' One panel per cluster: its share of the cell's total of `weight` (effort
#' by default), all rows pooled.
#'
#' @inheritParams plotDominance
#' @param cell_size Cell size (degrees). Default 1.
#' @param ncol Panels per row.
#'
#' @return A ggplot object.
#'
#' @family cluster maps
#' @export
plotClusterShare <- function(data, cluster, weight = "E", cell_size = 1, ncol = 2,
							 xlim = NULL, ylim = NULL) {
	.checkMapData(data, cluster, c("longitude", "latitude", weight))
	k  <- max(as.integer(cluster), na.rm = TRUE)
	sh <- .cellShare(data, cluster, weight, cell_size)
	sh$cluster <- factor(sh$cluster, levels = seq_len(k), labels = paste("cluster", seq_len(k)))
	ggplot2::ggplot(sh) + .landLayer() +
		ggplot2::geom_tile(ggplot2::aes(lon_cell, lat_cell, fill = share),
						   width = cell_size, height = cell_size) +
		ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
		ggplot2::facet_wrap(~ cluster, ncol = ncol, drop = FALSE) +
		ggplot2::coord_cartesian(xlim = .mapLims(xlim, sh$lon_cell, cell_size),
								 ylim = .mapLims(ylim, sh$lat_cell, cell_size), expand = FALSE) +
		ggplot2::labs(x = "Longitude", y = "Latitude",
					  fill = if (weight == "E") "Share of cell effort" else paste("Share of cell", weight)) +
		customTheme()
}


# ---- variable maps ---------------------------------------------------------------

#' Maps of CPUE, mean length or mean HBF per cluster
#'
#' One map per cluster of the variable per grid cell, all rows pooled, with
#' one colour scale shared by all clusters. With `profiles = TRUE` (needs the
#' `patchwork` package), each map has the same statistic per longitude band
#' above it and per latitude band on its right, computed from the rows of
#' each band.
#'
#' Values per cell: CPUE = `sum(yft_n) / sum(E)` (all rows); mean length
#' weighted by `yft_n` (rows with yellowfin catch, imputed lengths left out
#' when `len_source` is present); mean HBF = row mean (valid HBF only when
#' `hbf_outlier` is present).
#'
#' @inheritParams plotDominance
#' @param variable `"CPUE"`, `"mean_len"` or `"mean_hbf"`.
#' @param cell_size Cell size (degrees). Default 1.
#' @param limits Colour limits `c(min, max)`; values outside are drawn in
#'   the end colours. `NULL`: 2nd and 98th percentiles of the cells of all
#'   clusters.
#' @param trans `"identity"` or `"log10"` colour scale.
#' @param clusters Clusters drawn (all by default). The colour scale is
#'   still computed from all clusters, so maps drawn in several groups share
#'   one scale.
#' @param profiles Add longitude and latitude profiles (needs `patchwork`).
#' @param ncol Maps per row.
#'
#' @return A ggplot object, or a patchwork object with `profiles = TRUE`.
#'
#' @family cluster maps
#' @export
plotClusterMap <- function(data, cluster, variable = c("CPUE", "mean_len", "mean_hbf"),
						   cell_size = 1, limits = NULL, trans = "identity", clusters = NULL,
						   profiles = TRUE, ncol = 2, xlim = NULL, ylim = NULL) {
	variable <- match.arg(variable)
	need <- c("longitude", "latitude", switch(variable, CPUE = c("yft_n", "E"),
											  mean_len = c("mean_len", "yft_n"), mean_hbf = "mean_hbf"))
	.checkMapData(data, cluster, need)
	if (profiles && !requireNamespace("patchwork", quietly = TRUE)) {
		message("plotClusterMap(): profiles need the 'patchwork' package; drawing maps only.")
		profiles <- FALSE
	}
	ok  <- .mapRows(data, variable)
	d   <- data[ok, , drop = FALSE]
	cl  <- as.integer(cluster[ok])
	k   <- max(as.integer(cluster), na.rm = TRUE)
	if (is.null(clusters)) clusters <- seq_len(k)
	lon <- .gridCell(d$longitude %% 360, cell_size)
	lat <- .gridCell(d$latitude, cell_size)

	cv <- .groupValue(d, variable, paste(cl, lon, lat, sep = "\r"))
	cell <- cbind(.splitKey(names(cv), c("cluster", "lon_cell", "lat_cell")), value = unname(cv))
	cell[] <- lapply(cell, as.numeric)
	cell <- cell[is.finite(cell$value), , drop = FALSE]
	if (trans == "log10") cell$value[cell$value <= 0] <- NA
	if (is.null(limits)) limits <- unname(stats::quantile(cell$value, c(0.02, 0.98), na.rm = TRUE))
	xl <- .mapLims(xlim, cell$lon_cell, cell_size)
	yl <- .mapLims(ylim, cell$lat_cell, cell_size)
	lab <- c(CPUE = "CPUE", mean_len = "Mean length", mean_hbf = "Mean HBF")[[variable]]
	fill <- ggplot2::scale_fill_viridis_c(trans = trans, limits = limits, oob = .squish,
										  na.value = "grey30", name = lab)

	if (!profiles) {
		cd <- cell[cell$cluster %in% clusters, , drop = FALSE]
		cd$cluster <- factor(cd$cluster, levels = clusters, labels = paste("cluster", clusters))
		return(ggplot2::ggplot(cd) + .landLayer() +
			ggplot2::geom_tile(ggplot2::aes(lon_cell, lat_cell, fill = value),
							   width = cell_size, height = cell_size) +
			fill + ggplot2::facet_wrap(~ cluster, ncol = ncol, drop = FALSE) +
			ggplot2::coord_cartesian(xlim = xl, ylim = yl, expand = FALSE) +
			ggplot2::labs(x = "Longitude", y = "Latitude") + customTheme(text_size = 9))
	}

	prof <- function(band) {
		v <- .groupValue(d, variable, paste(cl, band, sep = "\r"))
		out <- cbind(.splitKey(names(v), c("cluster", "band")), value = unname(v))
		out[] <- lapply(out, as.numeric)
		out <- out[is.finite(out$value), , drop = FALSE]
		if (trans == "log10") out <- out[out$value > 0, , drop = FALSE]
		out[order(out$cluster, out$band), , drop = FALSE]
	}
	lonp <- prof(lon); latp <- prof(lat)
	y_scale <- if (trans == "log10") ggplot2::scale_y_log10() else ggplot2::scale_y_continuous()
	no_axis <- ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
	land <- .landLayer()

	panels <- lapply(clusters, function(cc) {
		p_map <- ggplot2::ggplot(cell[cell$cluster == cc, , drop = FALSE]) + land +
			ggplot2::geom_tile(ggplot2::aes(lon_cell, lat_cell, fill = value),
							   width = cell_size, height = cell_size) +
			fill + ggplot2::coord_cartesian(xlim = xl, ylim = yl, expand = FALSE) +
			ggplot2::labs(x = "Longitude", y = "Latitude") + customTheme(text_size = 9)
		p_lon <- ggplot2::ggplot(lonp[lonp$cluster == cc, , drop = FALSE], ggplot2::aes(band, value)) +
			ggplot2::geom_line() + y_scale + ggplot2::coord_cartesian(xlim = xl, expand = FALSE) +
			ggplot2::labs(x = NULL, y = NULL, title = paste("Cluster", cc)) +
			customTheme(text_size = 9) + no_axis
		p_lat <- ggplot2::ggplot(latp[latp$cluster == cc, , drop = FALSE], ggplot2::aes(band, value)) +
			ggplot2::geom_line() + y_scale + ggplot2::coord_flip(xlim = yl, expand = FALSE) +
			ggplot2::labs(x = NULL, y = NULL) + customTheme(text_size = 9) + no_axis
		p_lon + patchwork::plot_spacer() + p_map + p_lat +
			patchwork::plot_layout(widths = c(6, 1), heights = c(1, 4))
	})
	patchwork::wrap_plots(panels, ncol = ncol) + patchwork::plot_layout(guides = "collect") +
		patchwork::plot_annotation(theme = ggplot2::theme(legend.position = "bottom"))
}
