## R/breakpoints.R -- breakpoint detection and ecological validation,
## ported from myLibrary.r's breakpointAnalysis()/addBreakpointsToPlot()/
## summaryBP() and associated helpers.
##
## Design: analysis and plotting are fully separate. analyseBreakpoints()
## returns a plain data frame (no copies of the input series) meant to be
## saved; plotBreakpoints() reconstructs the figure from that data frame plus
## the original time series. Three previously-fused concerns -- statistical
## detection, age validation, spatial validation -- are now three functions
## that can be tested and called independently.
##
## Needs `changepoint` added to DESCRIPTION Imports (cpt.mean/cpts).
##
## Known open items, not resolved here:
##  - Linf/K/t0 defaults below need a literature citation (Baldrech to add).
##  - The age formula's t0 term is left un-scaled to days, exactly as in
##    myLibrary.r -- see check_age_formula_units.R for what that does to
##    results. Not changed pending a decision.
##  - calculateDistance()'s formal args are (lat1, lon1, lat2, lon2) but
##    every caller below passes (lon, lat, lon, lat), unchanged from the
##    original. Because both points get the same swapped treatment, this is
##    self-consistent (baseline and observed values are biased the same way),
##    but the distances are not true Haversine km and "Centroid_Distance" in
##    any output table is mislabelled. Not fixed here -- flagging it.


# ---- age validation ---------------------------------------------------

#' Convert mean fork length to implied age via a von Bertalanffy growth curve
#'
#' Inverts the von Bertalanffy growth equation to estimate age (in days) from
#' mean fork length. Used to check whether the time between two detected
#' breakpoints is consistent with a single cohort ageing through the
#' population, rather than an artefact of the changepoint detector.
#'
#' @param length_cm Numeric vector of fork lengths, in cm.
#' @param Linf Numeric; asymptotic length (cm). Default 150.3 -- TODO:
#'   citation needed.
#' @param K Numeric; growth coefficient (per year). Default 0.442.
#' @param t0 Numeric; theoretical age (years) at length zero. Default -0.244.
#'
#' @return Numeric vector of implied ages in days, same length as
#'   `length_cm`. `NA` where `length_cm >= Linf` (undefined under the model).
#'
#' @family breakpoint analysis utilities
#' @export
ageFromLength <- function(length_cm, Linf = 150.3, K = 0.442, t0 = -0.244) {
	over <- length_cm >= Linf
	if (any(over, na.rm = TRUE))
		warning("ageFromLength(): ", sum(over, na.rm = TRUE),
				" length(s) >= Linf (", Linf, "); returning NA for those.")
	age <- t0 - log(1 - length_cm / Linf) / K * 365   # t0 left un-scaled -- see file header
	age[over] <- NA
	age
}

#' Test whether breakpoints are consistent with cohort ageing
#'
#' For each candidate breakpoint, compares the time until the next breakpoint
#' (or the end of the series) against the age implied by mean length around
#' that breakpoint ([ageFromLength()]). A breakpoint passes if the fleet
#' could plausibly still be fishing the same cohort at that point, i.e. the
#' phase length does not exceed the implied age.
#'
#' @param bp_dates Date vector of candidate breakpoint dates, sorted ascending.
#' @param dates Date vector for the full analysis series (aligned with
#'   `mean_len`).
#' @param mean_len Numeric vector of mean length values, aligned to `dates`.
#' @param series_end Date; end of the analysis series, used as the phase end
#'   for the last breakpoint.
#' @param window_months Numeric; months either side of a breakpoint used to
#'   compute its implied age. Default 6.
#' @param skip Logical; if `TRUE`, bypasses the test -- every breakpoint
#'   passes with `age_limit_days = Inf` (for variables with no length-age
#'   relationship, e.g. `mean_hbf`). Default `FALSE`.
#' @param Linf,K,t0 Growth parameters passed to [ageFromLength()].
#'
#' @return A data frame, one row per element of `bp_dates`, with columns
#'   `bp_date`, `age_limit_days`, `phase_days`, `passed`.
#'
#' @family breakpoint analysis utilities
#' @export
validateBreakpointAge <- function(bp_dates, dates, mean_len, series_end,
								  window_months = 6, skip = FALSE,
								  Linf = 150.3, K = 0.442, t0 = -0.244) {
	n <- length(bp_dates)
	if (n == 0)
		return(data.frame(bp_date = as.Date(character()), age_limit_days = numeric(),
						  phase_days = numeric(), passed = logical()))

	phase_days <- c(as.numeric(diff(bp_dates)), as.numeric(series_end - bp_dates[n]))

	if (skip)
		return(data.frame(bp_date = bp_dates, age_limit_days = Inf,
						  phase_days = phase_days, passed = TRUE))

	## seq.Date(by = "n months") mirrors adding whole calendar months, EXCEPT
	## at month-end overflow (e.g. Jan 31 + 1 month rolls to Mar 3 rather than
	## clamping to Feb 28/29, unlike lubridate::months()). Not handled here.
	age_limit_days <- vapply(bp_dates, function(bp) {
		window_start <- seq(bp, by = paste(-window_months, "months"), length.out = 2)[2]
		window_end   <- seq(bp, by = paste( window_months, "months"), length.out = 2)[2]
		in_window <- dates >= window_start & dates <= window_end
		stats::median(ageFromLength(mean_len[in_window], Linf, K, t0), na.rm = TRUE)
	}, numeric(1))

	data.frame(bp_date = bp_dates, age_limit_days = age_limit_days,
			  phase_days = phase_days, passed = phase_days <= age_limit_days)
}


# ---- spatial validation -------------------------------------------------

#' Great-circle (Haversine) distance between two points
#'
#' @param lat1,lon1,lat2,lon2 Numeric; coordinates in decimal degrees.
#' @return Numeric distance in km (see file header re: caller argument order).
#' @family breakpoint analysis utilities
#' @export
calculateDistance <- function(lat1, lon1, lat2, lon2) {
	lat1r <- lat1 * pi / 180; lon1r <- lon1 * pi / 180
	lat2r <- lat2 * pi / 180; lon2r <- lon2 * pi / 180
	dlat <- lat2r - lat1r; dlon <- lon2r - lon1r
	a <- sin(dlat / 2)^2 + cos(lat1r) * cos(lat2r) * sin(dlon / 2)^2
	6371 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

#' Spread of fishing positions around a centroid
#'
#' @param longitude,latitude Numeric vectors of positions.
#' @param centroid Numeric length-2 vector, `c(longitude, latitude)`.
#' @return Numeric; SD of distances (km) from each point to `centroid`. `0`
#'   if fewer than 2 points.
#' @family breakpoint analysis utilities
#' @export
calculateFleetSpread <- function(longitude, latitude, centroid) {
	if (length(longitude) <= 1) return(0)
	d <- mapply(function(lon, lat) calculateDistance(lon, lat, centroid[1], centroid[2]),
			   longitude, latitude)
	stats::sd(d, na.rm = TRUE)
}

#' Spatial (Jaccard) overlap between two sets of points
#'
#' Grids both point sets onto a common lon/lat grid and returns the Jaccard
#' index (grid-cell overlap) between them.
#'
#' @param points1,points2 Data frames with `longitude`/`latitude` columns.
#' @param cell_size Numeric; grid cell size in degrees. Default 1.
#' @return Numeric in \[0, 1\]; 0 if the two sets share no grid cell.
#' @family breakpoint analysis utilities
#' @export
calculateSpatialOverlap <- function(points1, points2, cell_size = 1) {
	min_lon <- min(c(points1$longitude, points2$longitude), na.rm = TRUE)
	max_lon <- max(c(points1$longitude, points2$longitude), na.rm = TRUE)
	min_lat <- min(c(points1$latitude, points2$latitude), na.rm = TRUE)
	max_lat <- max(c(points1$latitude, points2$latitude), na.rm = TRUE)

	lon_range <- max_lon - min_lon; lat_range <- max_lat - min_lat
	min_lon <- min_lon - 0.1 * lon_range; max_lon <- max_lon + 0.1 * lon_range
	min_lat <- min_lat - 0.1 * lat_range; max_lat <- max_lat + 0.1 * lat_range

	lon_breaks <- seq(min_lon, max_lon, by = cell_size)
	lat_breaks <- seq(min_lat, max_lat, by = cell_size)

	assign_to_grid <- function(points) {
		n_cells <- length(lon_breaks) * length(lat_breaks)
		if (nrow(points) == 0) return(integer(n_cells))
		lon_idx <- findInterval(points$longitude, lon_breaks)
		lat_idx <- findInterval(points$latitude, lat_breaks)
		valid <- lon_idx > 0 & lon_idx < length(lon_breaks) &
				 lat_idx > 0 & lat_idx < length(lat_breaks)
		if (!any(valid)) return(integer(n_cells))
		cell_idx <- (lat_idx[valid] - 1) * (length(lon_breaks) - 1) + lon_idx[valid]
		tab <- table(cell_idx)
		result <- integer(n_cells)
		result[as.integer(names(tab))] <- as.integer(tab)
		result
	}

	g1 <- assign_to_grid(points1) > 0
	g2 <- assign_to_grid(points2) > 0
	union <- sum(g1 | g2)
	if (union == 0) return(0)
	sum(g1 & g2) / union
}

#' Empirical baseline of period-to-period fleet movement for one cluster
#'
#' Computes centroid distance, spatial overlap, and spread change between
#' every pair of consecutive observation dates within a cluster, and
#' summarises each as a distribution (mean/SD/quantiles) used as the
#' reference against which a candidate breakpoint's movement is judged in
#' [validateBreakpointSpatial()]. Recomputed on every call (no memoisation).
#'
#' @param df Data frame with `cluster`, `date`, `longitude`, `latitude`
#'   columns.
#' @param cluster_id Cluster identifier to subset to.
#'
#' @return A list with `raw_metrics` (the per-period-pair values) and
#'   `thresholds` (mean/SD/quantiles per metric).
#'
#' @family breakpoint analysis utilities
#' @export
calculateClusterMovementBaseline <- function(df, cluster_id) {
	cluster_data <- df[df$cluster == cluster_id, , drop = FALSE]
	dates <- sort(unique(cluster_data$date))

	centroid_distances <- spatial_overlaps <- spread_changes <- numeric(0)

	for (i in seq_len(max(length(dates) - 1, 0))) {
		data_start <- cluster_data[cluster_data$date == dates[i], , drop = FALSE]
		data_end   <- cluster_data[cluster_data$date == dates[i + 1], , drop = FALSE]
		if (nrow(data_start) == 0 || nrow(data_end) == 0) next

		start_centroid <- c(mean(data_start$longitude, na.rm = TRUE), mean(data_start$latitude, na.rm = TRUE))
		end_centroid   <- c(mean(data_end$longitude, na.rm = TRUE),   mean(data_end$latitude, na.rm = TRUE))

		centroid_distance <- calculateDistance(start_centroid[1], start_centroid[2],
											   end_centroid[1], end_centroid[2])
		start_spread <- calculateFleetSpread(data_start$longitude, data_start$latitude, start_centroid)
		end_spread   <- calculateFleetSpread(data_end$longitude, data_end$latitude, end_centroid)
		spread_change <- if (start_spread > 0) abs(end_spread - start_spread) / start_spread else 0
		overlap <- calculateSpatialOverlap(data_start[, c("longitude", "latitude")],
										   data_end[,   c("longitude", "latitude")])

		centroid_distances <- c(centroid_distances, centroid_distance)
		spatial_overlaps   <- c(spatial_overlaps, overlap)
		spread_changes      <- c(spread_changes, spread_change)
	}

	summarise_metric <- function(x) list(
		mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE),
		q10 = stats::quantile(x, 0.10, na.rm = TRUE), q25 = stats::quantile(x, 0.25, na.rm = TRUE),
		q75 = stats::quantile(x, 0.75, na.rm = TRUE), q90 = stats::quantile(x, 0.90, na.rm = TRUE)
	)

	list(
		raw_metrics = list(centroid_distances = centroid_distances,
						  spatial_overlaps = spatial_overlaps, spread_changes = spread_changes),
		thresholds = list(centroid_distance = summarise_metric(centroid_distances),
						 spatial_overlap = summarise_metric(spatial_overlaps),
						 spread_change = summarise_metric(spread_changes))
	)
}

#' Decide whether a breakpoint's spatial movement is significant
#'
#' Compares observed centroid distance / spatial overlap / spread change
#' against a cluster's empirical baseline ([calculateClusterMovementBaseline()])
#' at a chosen strictness level.
#'
#' @param centroid_distance,spatial_overlap,spread_change Numeric; observed
#'   values at the candidate breakpoint.
#' @param baseline_thresholds The `thresholds` element of
#'   [calculateClusterMovementBaseline()]'s return value.
#' @param strictness One of `"lenient"`, `"moderate"`, `"strict"`. Default
#'   `"moderate"`.
#'
#' @return A list: `is_significant`, `criteria_met`, `z_scores`,
#'   `thresholds_used`, `raw_values`.
#'
#' @family breakpoint analysis utilities
#' @export
evaluateBreakpointMovement <- function(centroid_distance, spatial_overlap, spread_change,
									   baseline_thresholds, strictness = "moderate") {
	thresholds <- baseline_thresholds
	if (strictness == "lenient") {
		distance_threshold <- thresholds$centroid_distance$q75
		overlap_threshold  <- thresholds$spatial_overlap$q25
		spread_threshold   <- thresholds$spread_change$q75
		min_criteria <- 1
	} else if (strictness == "moderate") {
		distance_threshold <- thresholds$centroid_distance$q90
		overlap_threshold  <- thresholds$spatial_overlap$q10
		spread_threshold   <- thresholds$spread_change$q90
		min_criteria <- 1
	} else if (strictness == "strict") {
		distance_threshold <- thresholds$centroid_distance$mean + 2 * thresholds$centroid_distance$sd
		overlap_threshold  <- max(0, thresholds$spatial_overlap$mean - 2 * thresholds$spatial_overlap$sd)
		spread_threshold   <- thresholds$spread_change$mean + 2 * thresholds$spread_change$sd
		min_criteria <- 2
	} else {
		## Original had no matching else branch: an unrecognised value left
		## distance_threshold etc. undefined, failing later with a cryptic
		## "object not found". Made explicit -- a deliberate change, not a
		## silent one.
		stop("strictness must be one of 'lenient', 'moderate', 'strict'")
	}

	criteria_met <- c(
		distance_significant = centroid_distance > distance_threshold,
		overlap_significant  = spatial_overlap < overlap_threshold,
		spread_significant   = spread_change > spread_threshold
	)
	z_scores <- c(
		distance_z = (centroid_distance - thresholds$centroid_distance$mean) / thresholds$centroid_distance$sd,
		overlap_z  = (thresholds$spatial_overlap$mean - spatial_overlap) / thresholds$spatial_overlap$sd,
		spread_z   = (spread_change - thresholds$spread_change$mean) / thresholds$spread_change$sd
	)

	list(is_significant = sum(criteria_met) >= min_criteria, criteria_met = criteria_met, z_scores = z_scores,
		thresholds_used = list(distance = distance_threshold, overlap = overlap_threshold, spread = spread_threshold),
		raw_values = list(centroid_distance = centroid_distance, spatial_overlap = spatial_overlap, spread_change = spread_change))
}

#' Test whether breakpoints coincide with the fleet relocating
#'
#' For each candidate breakpoint, compares the fleet's position just before
#' and just after against the cluster's empirical movement baseline. A
#' breakpoint passes (is spatially valid) when the fleet did *not* move
#' significantly -- a compositional shift explained by relocation isn't an
#' independent regime change.
#'
#' @param bp_dates Date vector of candidate breakpoint dates.
#' @param df Data frame with `cluster`, `date`, `longitude`, `latitude`
#'   columns (the raw, un-summarised fisheries data).
#' @param cluster_id Cluster identifier.
#' @param strictness Passed to [evaluateBreakpointMovement()]. Default
#'   `"moderate"`.
#' @param baseline Optional, pre-computed
#'   [calculateClusterMovementBaseline()] result; if `NULL` (default) it is
#'   computed here, once per call (no cross-call caching).
#'
#' @return A data frame, one row per element of `bp_dates` that had both a
#'   before- and after-period in `df` (breakpoints with no following
#'   observation are dropped -- matches the original, which never added such
#'   a breakpoint to `valid_bp_dates`), with columns `bp_date`,
#'   `centroid_distance`, `spatial_overlap`, `spread_change_ratio`, `passed`.
#'
#' @family breakpoint analysis utilities
#' @export
validateBreakpointSpatial <- function(bp_dates, df, cluster_id, strictness = "moderate", baseline = NULL) {
	if (is.null(baseline)) baseline <- calculateClusterMovementBaseline(df, cluster_id)

	rows <- lapply(bp_dates, function(bp_date) {
		data_start <- df[df$cluster == cluster_id & df$date == bp_date, , drop = FALSE]
		future_dates <- unique(df$date[df$cluster == cluster_id & df$date > bp_date])
		if (nrow(data_start) == 0 || length(future_dates) == 0) return(NULL)

		next_date <- min(future_dates)
		data_end <- df[df$cluster == cluster_id & df$date == next_date, , drop = FALSE]
		if (nrow(data_end) == 0) return(NULL)

		start_centroid <- c(mean(data_start$longitude, na.rm = TRUE), mean(data_start$latitude, na.rm = TRUE))
		end_centroid   <- c(mean(data_end$longitude, na.rm = TRUE),   mean(data_end$latitude, na.rm = TRUE))

		centroid_distance <- calculateDistance(start_centroid[1], start_centroid[2], end_centroid[1], end_centroid[2])
		start_spread <- calculateFleetSpread(data_start$longitude, data_start$latitude, start_centroid)
		end_spread   <- calculateFleetSpread(data_end$longitude, data_end$latitude, end_centroid)
		spread_change_ratio <- if (start_spread > 0) abs(end_spread - start_spread) / start_spread else 0
		overlap <- calculateSpatialOverlap(data_start[, c("longitude", "latitude")], data_end[, c("longitude", "latitude")])

		eval <- evaluateBreakpointMovement(centroid_distance, overlap, spread_change_ratio, baseline$thresholds, strictness)

		data.frame(bp_date = bp_date, centroid_distance = centroid_distance, spatial_overlap = overlap,
				  spread_change_ratio = spread_change_ratio, passed = !eval$is_significant)
	})

	rows <- rows[!vapply(rows, is.null, logical(1))]
	if (length(rows) == 0)
		return(data.frame(bp_date = as.Date(character()), centroid_distance = numeric(),
						  spatial_overlap = numeric(), spread_change_ratio = numeric(), passed = logical()))
	do.call(rbind, rows)
}


# ---- input prep ------------------------------------------------------------

#' Aggregate fisheries data into per-(cluster, date) time series
#'
#' Minimal, purpose-built aggregation for breakpoint analysis: catch per unit
#' effort (CPUE), mean length, and mean hooks-between-floats, one row per
#' (cluster, date). Not a general-purpose cluster summary -- see
#' [summaryClusters()] for per-cluster (not per-date) statistics, and note
#' this does none of `summarizeClusterData()`'s (myLibrary.r) extra work:
#' no lat/lon centroid, no fleet-movement distance, no imputation-flag
#' handling.
#'
#' @param df Data frame with `cluster`, `ymd`, `yft_n`, `E`, `mean_len`,
#'   `mean_hbf` columns. `mean_len`/`mean_hbf` are assumed NA-free (filter
#'   upstream, as `report.qmd` already does) -- this does not `na.rm`.
#'
#' @return A data frame with one row per (cluster, date) present in `df`,
#'   columns `cluster`, `date`, `CPUE`, `mean_len`, `mean_hbf`, sorted by
#'   cluster then date.
#'
#' @family breakpoint analysis utilities
#' @export
aggregateClusterTimeSeries <- function(df) {
	key   <- paste(as.character(df$cluster), format(df$ymd), sep = "\r")
	first <- !duplicated(key)
	d_cluster <- df$cluster[first]
	d_date    <- df$ymd[first]

	sums <- rowsum(as.matrix(df[, c("yft_n", "E")]), group = key, reorder = FALSE)
	CPUE <- sums[, "yft_n"] / sums[, "E"]

	mean_by_key <- function(x) {
		m <- rowsum(cbind(v = x, n = 1), group = key, reorder = FALSE)
		m[, "v"] / m[, "n"]
	}

	out <- data.frame(cluster = d_cluster, date = d_date, CPUE = CPUE,
					  mean_len = mean_by_key(df$mean_len), mean_hbf = mean_by_key(df$mean_hbf))
	out[order(out$cluster, out$date), ]
}


# ---- statistical detection ------------------------------------------------

#' Detect and classify candidate regime-change breakpoints in a time series
#'
#' Smooths the series with LOESS, runs `changepoint::cpt.mean()` (PELT) to
#' find raw changepoints, classifies the segment between each pair of
#' consecutive changepoints as increasing/decreasing/stable (via segment
#' slope significance and magnitude), and keeps only transitions between
#' stable segments whose means differ by more than `change_threshold`
#' (proportionally). Breakpoints within one year of the series end are
#' dropped. No domain (age/spatial) validation happens here.
#'
#' @param dates Date vector.
#' @param x Numeric vector, same length as `dates`; the variable to detect
#'   changepoints in.
#' @param min_segment_length Integer; minimum segment length for `cpt.mean()`
#'   and for a series to be analysed at all (needs > `2 * min_segment_length`
#'   points). Default 12.
#' @param penalty_factor Passed to `cpt.mean(penalty = )`. Default `"BIC"`.
#' @param smooth_span LOESS span. Default 0.2.
#' @param change_threshold Numeric; minimum proportional change between two
#'   stable segments to count as a candidate breakpoint. Default 0.03.
#' @param slope_p_threshold,slope_magnitude_threshold Thresholds for
#'   classifying a segment as increasing/decreasing vs stable. Defaults 0.1,
#'   0.05.
#'
#' @return `NULL` if there isn't enough data; otherwise a list with
#'   `bp_dates` (Date vector of candidate breakpoints, sorted) and
#'   `segment_info` (data frame describing every segment between raw
#'   changepoints -- kept for [plotBreakpoints()]'s segment-classification
#'   overlay, see note there).
#'
#' @family breakpoint analysis utilities
#' @importFrom changepoint cpt.mean cpts
#' @export
detectChangepoints <- function(dates, x, min_segment_length = 12, penalty_factor = "BIC",
							   smooth_span = 0.2, change_threshold = 0.03,
							   slope_p_threshold = 0.1, slope_magnitude_threshold = 0.05) {
	ord <- order(dates)
	dates <- dates[ord]; x <- x[ord]
	if (length(x) <= min_segment_length * 2) return(NULL)

	time_num <- as.numeric(dates - min(dates))
	loess_fit <- stats::loess(x ~ time_num, span = smooth_span, degree = 1, na.action = stats::na.exclude)
	smoothed <- stats::predict(loess_fit, newdata = data.frame(time_num = time_num))

	keep <- !is.na(smoothed)
	if (sum(keep) <= min_segment_length * 2) return(NULL)
	dates_k <- dates[keep]; ts_data <- smoothed[keep]

	cpt_result <- changepoint::cpt.mean(ts_data, method = "PELT", penalty = penalty_factor,
										minseglen = min_segment_length, test.stat = "Normal")
	changepoints <- changepoint::cpts(cpt_result)
	if (length(changepoints) == 0) return(list(bp_dates = as.Date(character()), segment_info = NULL))

	bp_indices <- changepoints[changepoints < length(ts_data)]
	segment_starts <- c(1, bp_indices + 1)
	segment_ends   <- c(bp_indices, length(ts_data))

	segment_info <- do.call(rbind, lapply(seq_along(segment_starts), function(i) {
		seg <- ts_data[segment_starts[i]:segment_ends[i]]
		seg_type <- "stable"; seg_slope <- 0; seg_p <- NA
		if (length(seg) > 3) {
			seg_time <- segment_starts[i]:segment_ends[i]
			fit <- stats::lm(seg ~ seg_time)
			seg_slope <- stats::coef(fit)[2]
			seg_p <- summary(fit)$coefficients[2, 4]
			if (seg_p < slope_p_threshold && abs(seg_slope) > stats::sd(seg) * slope_magnitude_threshold)
				seg_type <- if (seg_slope > 0) "increasing" else "decreasing"
		}
		data.frame(segment = i, start_idx = segment_starts[i], end_idx = segment_ends[i],
				  mean_value = mean(seg, na.rm = TRUE), start_value = seg[1], end_value = seg[length(seg)],
				  type = seg_type, slope = seg_slope, p_value = seg_p)
	}))

	candidate_idx <- c()
	for (i in seq_len(nrow(segment_info) - 1)) {
		cur <- segment_info[i, ]; nxt <- segment_info[i + 1, ]
		if (cur$type == "stable" && nxt$type == "stable") {
			pct <- abs(nxt$mean_value - cur$mean_value) / abs(cur$mean_value)
			if (pct > change_threshold) candidate_idx <- c(candidate_idx, bp_indices[i])
		} else if (cur$type == "stable" && i + 1 < nrow(segment_info)) {
			future <- segment_info[(i + 2):nrow(segment_info), , drop = FALSE]
			future_stable <- future[future$type == "stable", , drop = FALSE][1, ]
			if (!is.na(future_stable$segment)) {
				pct <- abs(future_stable$mean_value - cur$mean_value) / abs(cur$mean_value)
				if (pct > change_threshold) candidate_idx <- c(candidate_idx, bp_indices[i])
			}
		}
	}

	bp_dates <- sort(dates_k[unique(candidate_idx)])
	bp_dates <- bp_dates[bp_dates <= max(dates_k) - 365]
	list(bp_dates = bp_dates, segment_info = segment_info)
}


# ---- orchestrator (this is what should be saved) ---------------------------

#' Full breakpoint analysis for one (cluster, variable) time series
#'
#' Orchestrates [detectChangepoints()], [validateBreakpointAge()], and
#' [validateBreakpointSpatial()] into a single tidy result. This -- not a
#' `bp_results`-style list with embedded copies of the series -- is what
#' should be saved: it holds only breakpoint dates and validation outcomes.
#' Reconstructing a plot needs this data frame plus the original time
#' series, via [plotBreakpoints()].
#'
#' @param dates Date vector for the analysis series.
#' @param x Numeric vector, the variable to detect breakpoints in.
#' @param mean_len Numeric vector of mean length, aligned to `dates`; used
#'   for age validation. Ignored if `skip_age = TRUE`.
#' @param df Data frame with `cluster`, `date`, `longitude`, `latitude`
#'   (raw, un-summarised fisheries data); used for spatial validation.
#' @param cluster_id Cluster identifier, for spatial validation.
#' @param skip_age Logical; skip age validation for this variable (e.g.
#'   `mean_hbf`). Default `FALSE`.
#' @param window_months,Linf,K,t0 Passed to [validateBreakpointAge()].
#' @param movement_strictness Passed to [validateBreakpointSpatial()].
#' @param ... Passed to [detectChangepoints()].
#'
#' @return `NULL` if [detectChangepoints()] finds too little data; a 0-row
#'   data frame if it runs but finds no candidate breakpoints; otherwise one
#'   row per candidate breakpoint with columns `bp_date`, `age_limit_days`,
#'   `phase_days`, `passed_age`, `centroid_distance`, `spatial_overlap`,
#'   `spread_change_ratio`, `passed_spatial`, `valid`. A breakpoint with no
#'   spatial test result (dropped by [validateBreakpointSpatial()]) gets
#'   `valid = FALSE`, matching the original's behaviour of only counting a
#'   breakpoint valid when the spatial test explicitly passed.
#'
#' @family breakpoint analysis utilities
#' @export
analyseBreakpoints <- function(dates, x, mean_len, df, cluster_id, skip_age = FALSE,
							   window_months = 6, Linf = 150.3, K = 0.442, t0 = -0.244,
							   movement_strictness = "moderate", ...) {
	detected <- detectChangepoints(dates, x, ...)
	if (is.null(detected) || length(detected$bp_dates) == 0) return(detected)

	age <- validateBreakpointAge(detected$bp_dates, dates, mean_len, max(dates),
								 window_months = window_months, skip = skip_age, Linf = Linf, K = K, t0 = t0)
	spatial <- validateBreakpointSpatial(detected$bp_dates, df, cluster_id, strictness = movement_strictness)

	out <- merge(age, spatial, by = "bp_date", all = TRUE, suffixes = c("_age", "_spatial"))
	out$valid <- out$passed_age & !is.na(out$passed_spatial) & out$passed_spatial
	out
}


# ---- plotting (no analysis logic here) -------------------------------------

#' Overlay breakpoint analysis results on a time series plot
#'
#' Draws the raw series, its LOESS-smoothed curve, and a vertical line at
#' each candidate breakpoint (solid = valid, dashed = rejected), from the
#' output of [analyseBreakpoints()]. No analysis is performed here.
#'
#' @param dates Date vector for the series (same data passed to
#'   [analyseBreakpoints()]).
#' @param x Numeric vector, the plotted variable.
#' @param bp_table Data frame from [analyseBreakpoints()] (or `NULL`/0-row,
#'   plotted as a plain time series with no breakpoint lines).
#' @param smooth_span LOESS span for the overlaid smoothed curve; should
#'   match what [detectChangepoints()] was called with. Default 0.2.
#' @param y_lab Y-axis label. Default `NULL`.
#'
#' @return A ggplot2 object.
#'
#' @family breakpoint plotting utilities
#' @importFrom ggplot2 ggplot aes geom_line geom_vline scale_linetype_manual labs
#' @export
plotBreakpoints <- function(dates, x, bp_table, smooth_span = 0.2, y_lab = NULL) {
	plot_df <- data.frame(date = dates, value = x)
	p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = date, y = value)) +
		ggplot2::geom_line() +
		ggplot2::labs(x = NULL, y = y_lab) +
		customTheme()

	if (sum(!is.na(x)) > 5) {
		time_num <- as.numeric(dates - min(dates))
		fit <- stats::loess(x ~ time_num, span = smooth_span, degree = 1, na.action = stats::na.exclude)
		plot_df$smoothed <- stats::predict(fit, newdata = data.frame(time_num = time_num))
		p <- p + ggplot2::geom_line(data = plot_df, ggplot2::aes(y = smoothed), colour = "blue")
	}

	if (!is.null(bp_table) && nrow(bp_table) > 0) {
		bp_table$status <- ifelse(bp_table$valid, "valid", "rejected")
		p <- p +
			ggplot2::geom_vline(data = bp_table,
								ggplot2::aes(xintercept = as.numeric(bp_date), linetype = status),
								colour = "firebrick") +
			ggplot2::scale_linetype_manual(name = NULL, values = c(valid = "solid", rejected = "dashed"))
	}
	p
}
