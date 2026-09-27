## Tests for R/breakpointsV2.R. Synthetic series with large signal-to-noise
## so outcomes do not depend on package versions of EnvCpt / changepoint.

mkSeries <- function(n = 420, seed = 1, sd = 1, shift = 0, at = 200, duration = 1) {
	set.seed(seed)
	dates <- seq(as.Date("1990-01-15"), by = "month", length.out = n)
	r <- pmin(pmax((seq_len(n) - at + 1) / duration, 0), 1)
	list(dates = dates, x = 100 + shift * r + stats::rnorm(n, sd = sd))
}

test_that("vbAgeMonths() uses t0 in years and returns NA at or above Linf", {
	expect_equal(vbAgeMonths(0), -0.244 * 12)
	L <- 117
	expect_equal(vbAgeMonths(L), (-0.244 - log(1 - L / 150.3) / 0.442) * 12)
	expect_true(is.na(vbAgeMonths(150.3)))
	expect_true(is.na(vbAgeMonths(NA)))
})

test_that("analyseBreakpointsV2() returns a typed 0-row frame on a flat series", {
	skip_if_not_installed("EnvCpt")
	s <- mkSeries(seed = 2)
	out <- analyseBreakpointsV2(s$dates, s$x, "mean_len", rep(117, 420))
	expect_s3_class(out, "data.frame")
	expect_true(all(c("start_date", "d_hat", "label") %in% names(out)))
	expect_false(any(out$label %in% c("abrupt", "unexplained")))
})

test_that("a large abrupt step is found at the right time and classified abrupt", {
	skip_if_not_installed("EnvCpt")
	s <- mkSeries(shift = 8, seed = 3)
	out <- analyseBreakpointsV2(s$dates, s$x, "mean_len", rep(117, 420))
	expect_gte(nrow(out), 1)
	best <- out[which.max(out$lr), ]
	expect_lte(abs(as.numeric(best$start_date - s$dates[200])), 62)   # within ~2 months
	expect_equal(best$class, "abrupt")
	expect_equal(best$label, "abrupt")   # no spatial table -> not split further
})

test_that("a slow ramp is not classified abrupt", {
	skip_if_not_installed("EnvCpt")
	s <- mkSeries(shift = 8, seed = 4, duration = 96, at = 150)
	out <- analyseBreakpointsV2(s$dates, s$x, "mean_len", rep(117, 420))
	expect_false(any(out$class == "abrupt"))
})

test_that("CPUE is analysed on the log scale (a x2 step is found whatever the level)", {
	skip_if_not_installed("EnvCpt")
	set.seed(5)
	dates <- seq(as.Date("1990-01-15"), by = "month", length.out = 420)
	x <- 0.05 * exp(c(rep(0, 200), rep(log(2), 220)) + stats::rnorm(420, sd = 0.1))
	out <- analyseBreakpointsV2(dates, x, "CPUE", rep(117, 420))
	expect_gte(sum(out$class == "abrupt"), 1)
	expect_true(any(abs(out$delta_rel - 1) < 0.3))   # ~ +100 %
})

test_that("a composition change labels an abrupt change 'composition'", {
	skip_if_not_installed("EnvCpt")
	set.seed(6)
	dates <- seq(as.Date("1990-01-15"), by = "month", length.out = 420)
	rows <- do.call(rbind, lapply(seq_along(dates), function(i) {
		fl <- if (i < 200) c("A", "B") else c("A", "C")      # flag B replaced by C
		data.frame(cluster = 1, date = dates[i], flag = rep(fl, each = 10),
				   longitude = 160 + stats::rnorm(20), latitude = -10 + stats::rnorm(20))
	}))
	tab <- spatialMonthTable(rows, 1)
	x <- 100 + c(rep(0, 199), rep(8, 221)) + stats::rnorm(420)
	out <- analyseBreakpointsV2(dates, x, "mean_len", rep(117, 420), spatial_table = tab)
	best <- out[which.max(out$lr), ]
	expect_equal(best$label, "composition")
})

test_that("flagCommonBreaks() needs the share of clusters and at least 2", {
	bp <- data.frame(cluster = c(1, 2, 3, 1), variable = c("CPUE", "CPUE", "CPUE", "mean_len"),
					 start_date = as.Date(c("2005-01-15", "2005-04-15", "2010-01-15", "2005-01-15")))
	out <- flagCommonBreaks(bp, n_clusters = 3, tol_months = 6, min_share = 0.6)
	expect_equal(out$common, c(TRUE, TRUE, FALSE, FALSE))
	expect_false(any(flagCommonBreaks(bp, n_clusters = 1)$common))
})

test_that("analyseClusterBreakpoints() and plotBreakpointsV2() handle series with no change", {
	skip_if_not_installed("EnvCpt")
	set.seed(7)
	dates <- seq(as.Date("1990-01-15"), by = "month", length.out = 120)
	ts <- data.frame(cluster = rep(1:2, each = 120), date = rep(dates, 2),
					 CPUE = exp(stats::rnorm(240, sd = 0.05)), mean_len = 117 + stats::rnorm(240),
					 mean_hbf = 20 + stats::rnorm(240))
	out <- analyseClusterBreakpoints(ts)
	expect_true(all(c("cluster", "variable", "label", "common") %in% names(out)))
	expect_s3_class(plotBreakpointsV2(dates, ts$mean_len[1:120], out[0, ]), "ggplot")
})

test_that("flagCommonBreaks() never marks or counts composition / switch changes", {
	bp <- data.frame(cluster = c(1, 2), variable = "CPUE",
					 start_date = as.Date(c("2008-01-15", "2008-02-15")),
					 label = c("switch", "switch"))
	expect_false(any(flagCommonBreaks(bp, n_clusters = 2)$common))
	bp$label <- c("unexplained", "unexplained")
	expect_true(all(flagCommonBreaks(bp, n_clusters = 2)$common))
})

test_that("an ambiguous change can still be labelled composition, never unexplained", {
	skip_if_not_installed("EnvCpt")
	set.seed(8)
	dates <- seq(as.Date("1990-01-15"), by = "month", length.out = 420)
	rows <- do.call(rbind, lapply(seq_along(dates), function(i) {
		fl <- if (i < 200) c("A", "B") else c("A", "C")
		data.frame(cluster = 1, date = dates[i], flag = rep(fl, each = 10),
				   longitude = 160 + stats::rnorm(20), latitude = -10 + stats::rnorm(20))
	}))
	tab <- spatialMonthTable(rows, 1)
	x <- 100 + c(rep(0, 199), rep(8, 221)) + stats::rnorm(420)
	## timescale set so that any estimated duration > 1 month straddles it
	out <- analyseBreakpointsV2(dates, x, "mean_len", rep(117, 420), spatial_table = tab, abrupt_months = 1)
	expect_false(any(out$class == "ambiguous" & out$label == "unexplained"))
	expect_true(all(out$label[out$class == "ambiguous"] %in% c("ambiguous", "composition", "relocation")))
})

test_that("rows moving between two clusters are labelled 'switch', a new fleet stays 'composition'", {
	skip_if_not_installed("EnvCpt")
	set.seed(9)
	dates <- seq(as.Date("1990-01-15"), by = "month", length.out = 420)
	mk <- function(move) do.call(rbind, lapply(seq_along(dates), function(i) {
		after <- i >= 200
		## cluster 1: flags A, B; cluster 2: flags C, D. From month 200, flag B
		## either moves to cluster 2 (move = TRUE) or leaves cluster 1 and a new
		## flag E appears in cluster 2 (move = FALSE).
		fl1 <- if (after) "A" else c("A", "B")
		fl2 <- if (!after) c("C", "D") else if (move) c("C", "D", "B") else c("C", "D", "E")
		rbind(data.frame(cluster = 1, date = dates[i], flag = rep(fl1, each = 10)),
			  data.frame(cluster = 2, date = dates[i], flag = rep(fl2, each = 10)))
	}))
	res <- lapply(c(TRUE, FALSE), function(move) {
		rows <- mk(move)
		rows$longitude <- 160 + stats::rnorm(nrow(rows)); rows$latitude <- -10 + stats::rnorm(nrow(rows))
		ts <- data.frame(cluster = rep(1:2, each = 420), date = rep(dates, 2),
						 CPUE = exp(stats::rnorm(840, sd = 0.05)),
						 mean_len = c(117 + stats::rnorm(420), 110 + c(rep(0, 199), rep(8, 221)) + stats::rnorm(420)),
						 mean_hbf = 20 + stats::rnorm(840))
		out <- analyseClusterBreakpoints(ts, rows, variables = "mean_len")
		out[out$cluster == 2, ][which.max(out$lr[out$cluster == 2]), ]
	})
	expect_equal(res[[1]]$label, "switch")
	expect_equal(res[[1]]$mirror_cluster, "1")
	expect_equal(res[[2]]$label, "composition")
	expect_false(res[[1]]$common)
})
