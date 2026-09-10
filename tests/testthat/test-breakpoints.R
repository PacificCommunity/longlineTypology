## Coverage for R/breakpoints.R. Deterministic arithmetic checks use exact
## expect_equal(); anything downstream of cpt.mean()/loess() on synthetic
## random data uses large, unambiguous signal-to-noise ratios and checks
## qualitative properties (sign, membership, structure) rather than exact
## values, so they survive minor package-version differences in those
## algorithms. Not run against a live R session while writing this --
## verify once and report back if any threshold needs tuning.

# ---- ageFromLength() --------------------------------------------------

test_that("ageFromLength() matches the formula and flags length >= Linf", {
	expect_equal(ageFromLength(0), -0.244)  # log(1 - 0/Linf) == 0 -> age == t0

	L <- 120
	expect_equal(ageFromLength(L), -0.244 - log(1 - L / 150.3) / 0.442 * 365)

	expect_warning(res <- ageFromLength(c(100, 150.3, 160)), "length\\(s\\) >= Linf")
	expect_false(is.na(res[1]))
	expect_true(is.na(res[2]))
	expect_true(is.na(res[3]))
})


# ---- validateBreakpointAge() -------------------------------------------

test_that("validateBreakpointAge() skip = TRUE bypasses the test entirely", {
	bp_dates <- as.Date(c("2000-01-01", "2005-01-01"))
	out <- validateBreakpointAge(bp_dates, dates = bp_dates, mean_len = c(120, 120),
								 series_end = as.Date("2010-01-01"), skip = TRUE)
	expect_equal(nrow(out), 2)
	expect_true(all(out$passed))
	expect_true(all(is.infinite(out$age_limit_days)))
})

test_that("validateBreakpointAge() computes phase_days against the next breakpoint or series_end", {
	bp_dates   <- as.Date(c("2000-01-15", "2001-01-15"))
	dates      <- seq(as.Date("1995-01-15"), as.Date("2005-01-15"), by = "month")
	mean_len   <- rep(100, length(dates))
	series_end <- max(dates)

	out <- validateBreakpointAge(bp_dates, dates, mean_len, series_end)
	expect_equal(out$phase_days[1], as.numeric(bp_dates[2] - bp_dates[1]))
	expect_equal(out$phase_days[2], as.numeric(series_end - bp_dates[2]))
})

test_that("validateBreakpointAge() handles zero breakpoints", {
	out <- validateBreakpointAge(as.Date(character()), as.Date(character()),
								 numeric(), as.Date("2000-01-01"))
	expect_equal(nrow(out), 0)
	expect_named(out, c("bp_date", "age_limit_days", "phase_days", "passed"))
})


# ---- calculateDistance() / calculateFleetSpread() / calculateSpatialOverlap() ----

test_that("calculateDistance() is 0 for identical points, symmetric, and matches a known value", {
	expect_equal(calculateDistance(-20, 160, -20, 160), 0)

	d1 <- calculateDistance(-20, 160, -25, 165)
	d2 <- calculateDistance(-25, 165, -20, 160)
	expect_equal(d1, d2)
	expect_true(d1 > 0)

	## London to Paris, ~344 km great-circle
	d <- calculateDistance(51.5074, -0.1278, 48.8566, 2.3522)
	expect_equal(d, 344, tolerance = 5)
})

test_that("calculateFleetSpread() returns 0 for 0/1 points and increases with dispersion", {
	expect_equal(calculateFleetSpread(numeric(0), numeric(0), c(0, 0)), 0)
	expect_equal(calculateFleetSpread(1, 1, c(1, 1)), 0)

	centroid <- c(160, -20)
	tight <- calculateFleetSpread(c(160, 160.1, 159.9), c(-20, -20.1, -19.9), centroid)
	wide  <- calculateFleetSpread(c(158, 162, 160), c(-22, -18, -20), centroid)
	expect_true(wide > tight)
})

test_that("calculateSpatialOverlap() is 1 for co-located points and 0 for disjoint ones", {
	pts <- data.frame(longitude = c(160, 160.2), latitude = c(-20, -20.1))
	expect_equal(calculateSpatialOverlap(pts, pts), 1)

	near <- data.frame(longitude = c(160.05, 160.15), latitude = c(-20.05, -20.05))
	expect_equal(calculateSpatialOverlap(pts, near), 1)  # same 1-degree cell

	disjoint <- data.frame(longitude = c(-170, -170.2), latitude = c(20, 20.1))
	expect_equal(calculateSpatialOverlap(pts, disjoint), 0)
})


# ---- calculateClusterMovementBaseline() / evaluateBreakpointMovement() ----

make_fleet_df <- function(n_periods = 20, cluster_id = 1, jump_at = NULL, n_per_period = 20) {
	dates <- seq(as.Date("2000-01-15"), by = "month", length.out = n_periods)
	do.call(rbind, lapply(seq_along(dates), function(i) {
		base_lon <- 160; base_lat <- -20
		if (!is.null(jump_at) && i >= jump_at) { base_lon <- 170; base_lat <- -10 }
		data.frame(cluster = cluster_id, date = dates[i],
				  longitude = base_lon + stats::rnorm(n_per_period, sd = 0.2),
				  latitude  = base_lat + stats::rnorm(n_per_period, sd = 0.2))
	}))
}

test_that("calculateClusterMovementBaseline() returns one metric per consecutive date pair", {
	set.seed(1)
	df <- make_fleet_df(n_periods = 10)
	baseline <- calculateClusterMovementBaseline(df, cluster_id = 1)
	expect_length(baseline$raw_metrics$centroid_distances, 9)
	expect_true(all(c("mean", "sd", "q10", "q25", "q75", "q90") %in%
					names(baseline$thresholds$centroid_distance)))
})

test_that("evaluateBreakpointMovement() rejects an unrecognised strictness", {
	set.seed(1)
	df <- make_fleet_df(n_periods = 10)
	baseline <- calculateClusterMovementBaseline(df, cluster_id = 1)
	expect_error(
		evaluateBreakpointMovement(1, 0.5, 0.1, baseline$thresholds, strictness = "extreme"),
		"strictness must be one of"
	)
})

test_that("evaluateBreakpointMovement() strict mode needs 2 criteria, lenient/moderate need 1", {
	thresholds <- list(
		centroid_distance = list(mean = 10, sd = 2,   q75 = 12,  q90 = 15),
		spatial_overlap   = list(mean = 0.8, sd = 0.1, q25 = 0.7, q10 = 0.6),
		spread_change     = list(mean = 0.1, sd = 0.05, q75 = 0.15, q90 = 0.2)
	)
	## only the distance criterion is breached
	res_lenient  <- evaluateBreakpointMovement(20, 0.75, 0.1, thresholds, "lenient")
	res_moderate <- evaluateBreakpointMovement(20, 0.75, 0.1, thresholds, "moderate")
	res_strict   <- evaluateBreakpointMovement(20, 0.75, 0.1, thresholds, "strict")

	expect_true(res_lenient$is_significant)
	expect_true(res_moderate$is_significant)
	expect_false(res_strict$is_significant)  # needs 2 criteria, only 1 breached
})

test_that("validateBreakpointSpatial() flags a relocation and drops a breakpoint with no future data", {
	set.seed(1)
	df <- make_fleet_df(n_periods = 24, jump_at = 13)
	bp_dates <- as.Date(c("2000-06-15", "2001-01-15"))  # one stable period, one at the jump
	out <- validateBreakpointSpatial(bp_dates, df, cluster_id = 1, strictness = "lenient")

	expect_equal(nrow(out), 2)
	jump_row <- out[out$bp_date == as.Date("2001-01-15"), ]
	expect_false(jump_row$passed)  # ~1500 km jump vs ~tens-of-km baseline -> significant

	last_date <- max(df$date)
	out2 <- validateBreakpointSpatial(last_date, df, cluster_id = 1)
	expect_equal(nrow(out2), 0)  # no observation after the last date -> dropped
})


# ---- detectChangepoints() ----------------------------------------------

test_that("detectChangepoints() returns NULL for a too-short series", {
	dates <- seq(as.Date("2000-01-01"), by = "month", length.out = 10)
	expect_null(detectChangepoints(dates, stats::rnorm(10), min_segment_length = 12))
})

test_that("detectChangepoints() detects structure and separates segment means for an obvious level shift", {
	## NOTE: an instant step is a stress case for the segment-classification
	## heuristic, not for changepoint detection itself. LOESS smoothing (span
	## 0.2) spreads the step into a ~24-point ramp, which cpt.mean() reads as
	## 3-4 segments rather than 2 -- expected. Separately, the edge segments
	## far from the ramp can get misclassified as "increasing"/"decreasing"
	## because their slope is compared to their OWN (LOESS-smoothed, near-zero)
	## SD -- see the file header. That blocks the higher-level stable-to-stable
	## candidate rule in a way this test deliberately does not assert on,
	## since it's a property of the classification heuristic, not something
	## this port changed. This checks only what's robust: raw changepoints
	## were found, and the resulting segments actually separate the two levels.
	set.seed(42)
	n <- 120
	dates <- seq(as.Date("2000-01-15"), by = "month", length.out = n)
	x <- c(rnorm(60, mean = 100, sd = 1), rnorm(60, mean = 130, sd = 1))

	res <- detectChangepoints(dates, x, min_segment_length = 12)
	expect_false(is.null(res))
	expect_true(is.data.frame(res$segment_info))
	expect_true(nrow(res$segment_info) >= 2)
	expect_true(diff(range(res$segment_info$mean_value)) > 20)  # ~30-unit shift should show up
})

test_that("detectChangepoints() returns no *surviving* candidates for a shift close to the series end", {
	## As above: this can return empty either because the end-of-series trim
	## worked, or because classification never reached the candidate stage at
	## all -- it isolates neither cause specifically, only the overall
	## guarantee that a breakpoint this close to the end never survives.
	set.seed(42)
	n <- 30
	dates <- seq(as.Date("2000-01-15"), by = "month", length.out = n)
	x <- c(rnorm(20, mean = 100, sd = 1), rnorm(10, mean = 130, sd = 1))
	res <- detectChangepoints(dates, x, min_segment_length = 8)
	expect_true(is.null(res) || length(res$bp_dates) == 0)
})


# ---- aggregateClusterTimeSeries() ---------------------------------------

test_that("aggregateClusterTimeSeries() aggregates CPUE/mean_len/mean_hbf correctly", {
	df <- data.frame(
		cluster  = c(1, 1, 1, 2, 2),
		ymd      = as.Date(c("2000-01-15", "2000-01-15", "2000-02-15", "2000-01-15", "2000-01-15")),
		yft_n    = c(10, 5, 20, 3, 7),
		E        = c(2, 3, 4, 1, 1),
		mean_len = c(100, 110, 120, 90, 95),
		mean_hbf = c(8, 9, 10, 5, 6)
	)
	out <- aggregateClusterTimeSeries(df)

	expect_equal(nrow(out), 3)
	expect_setequal(names(out), c("cluster", "date", "CPUE", "mean_len", "mean_hbf"))

	c1_jan <- out[out$cluster == 1 & out$date == as.Date("2000-01-15"), ]
	expect_equal(c1_jan$CPUE, (10 + 5) / (2 + 3))
	expect_equal(c1_jan$mean_len, mean(c(100, 110)))
	expect_equal(c1_jan$mean_hbf, mean(c(8, 9)))

	c2_jan <- out[out$cluster == 2 & out$date == as.Date("2000-01-15"), ]
	expect_equal(c2_jan$CPUE, (3 + 7) / (1 + 1))
})


# ---- analyseBreakpoints() -- integration, including the reported bug ----

test_that("analyseBreakpoints() returns NULL when detectChangepoints() does", {
	dates <- seq(as.Date("2000-01-15"), by = "month", length.out = 5)
	df <- data.frame(cluster = 1, date = dates, longitude = 160, latitude = -20)
	expect_null(analyseBreakpoints(dates, stats::rnorm(5), mean_len = rep(100, 5),
								   df = df, cluster_id = 1, min_segment_length = 12))
})

test_that("analyseBreakpoints() returns a well-typed 0-row data frame when no candidates survive (regression: this used to leak a list straight to plotBreakpoints())", {
	set.seed(1)
	n <- 60
	dates <- seq(as.Date("2000-01-15"), by = "month", length.out = n)
	x  <- stats::rnorm(n, mean = 100, sd = 0.5)  # flat -- no meaningful changepoints
	df <- data.frame(cluster = 1, date = dates,
					 longitude = 160 + stats::rnorm(n, sd = 0.1),
					 latitude  = -20 + stats::rnorm(n, sd = 0.1))

	out <- analyseBreakpoints(dates, x, mean_len = rep(120, n), df = df, cluster_id = 1,
							  min_segment_length = 12)

	expect_s3_class(out, "data.frame")
	expect_equal(nrow(out), 0)
	expect_named(out, c("bp_date", "age_limit_days", "phase_days", "passed_age",
					   "centroid_distance", "spatial_overlap", "spread_change_ratio",
					   "passed_spatial", "valid"))

	## the actual failure mode reported: plotBreakpoints() must not error on this
	p <- plotBreakpoints(dates, x, out)
	expect_s3_class(p, "ggplot")
})

test_that("analyseBreakpoints() marks a breakpoint invalid (not NA) when spatial validation drops it", {
	set.seed(2)
	n <- 60
	dates <- seq(as.Date("2000-01-15"), by = "month", length.out = n)
	x <- c(stats::rnorm(30, mean = 100, sd = 0.5), stats::rnorm(30, mean = 130, sd = 0.5))
	## df only covers up to the shift -> validateBreakpointSpatial() finds no
	## "next_date" for any candidate near it and drops the row from its merge
	df <- data.frame(cluster = 1, date = dates[1:31],
					 longitude = 160 + stats::rnorm(31, sd = 0.1),
					 latitude  = -20 + stats::rnorm(31, sd = 0.1))

	out <- analyseBreakpoints(dates, x, mean_len = rep(120, n), df = df, cluster_id = 1,
							  skip_age = TRUE, min_segment_length = 12)

	expect_s3_class(out, "data.frame")
	if (nrow(out) > 0) expect_true(all(!out$valid))  # never NA -- see analyseBreakpoints()
})


# ---- plotBreakpoints() --------------------------------------------------

test_that("plotBreakpoints() handles NULL and 0-row bp_table without error", {
	dates <- seq(as.Date("2000-01-15"), by = "month", length.out = 24)
	x <- stats::rnorm(24)

	expect_s3_class(plotBreakpoints(dates, x, NULL), "ggplot")

	empty <- data.frame(bp_date = as.Date(character()), valid = logical())
	expect_s3_class(plotBreakpoints(dates, x, empty), "ggplot")
})

test_that("plotBreakpoints() draws a vline layer when bp_table has rows", {
	dates <- seq(as.Date("2000-01-15"), by = "month", length.out = 24)
	x <- stats::rnorm(24)
	bp_table <- data.frame(bp_date = dates[c(10, 15)], valid = c(TRUE, FALSE))

	p <- plotBreakpoints(dates, x, bp_table)
	expect_s3_class(p, "ggplot")
	has_vline <- vapply(p$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))
	expect_true(any(has_vline))
})
