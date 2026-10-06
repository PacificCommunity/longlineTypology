## Tests for R/breakpointInspect.R. Two clusters, two flags each; a large
## CPUE step in cluster 1 and flag "B" moving from cluster 2 to cluster 1.

mkRows <- function(seed = 1) {
	set.seed(seed)
	dates <- seq(as.Date("1990-01-15"), by = "month", length.out = 360)
	one <- function(flag, cluster, lon, cpue, keep = rep(TRUE, 360)) {
		n <- 20
		d <- data.frame(flag = flag, cluster = cluster, date = rep(dates, each = n),
						longitude = rep(lon, 360 * n) + stats::rnorm(360 * n),
						latitude = stats::rnorm(360 * n), E = 100,
						yft_n = stats::rpois(360 * n, rep(cpue, each = n) * 100),
						mean_len = stats::rnorm(360 * n, 117, 2), mean_hbf = stats::rnorm(360 * n, 20, 1))
		d[rep(keep, each = n), ]
	}
	step <- ifelse(seq_len(360) >= 180, 0.2, 1)
	moved <- seq_len(360) >= 240
	rows <- rbind(one("A", 1, 170, step), one("B", 1, 175, 0.5, moved),
				  one("B", 2, 175, 0.5, !moved), one("C", 2, 200, 0.5))
	rows$ymd <- rows$date
	rows
}

test_that("rankBreakpoints() puts unexplained first and common last", {
	bp <- data.frame(cluster = 1:4, variable = "CPUE", start_date = as.Date("2000-01-15") + 0:3,
					 label = c("gradual", "unexplained", "ambiguous", "unexplained"),
					 common = c(FALSE, FALSE, FALSE, TRUE))
	r <- rankBreakpoints(bp)
	expect_equal(r$cluster, c(2, 3, 1, 4))
	expect_equal(r$priority, c(1, 3, 7, 8))
})

test_that(".bpEmpThreshold() is the value an observation must exceed for p < alpha", {
	v <- c(1:99, NA)
	thr <- longlineTypology:::.bpEmpThreshold(v, 0.05)
	empP <- function(o) (1 + sum(v >= o, na.rm = TRUE)) / (1 + sum(!is.na(v)))
	expect_gte(empP(thr), 0.05)
	expect_lt(empP(thr + 1e-9), 0.05)
	expect_true(is.na(longlineTypology:::.bpEmpThreshold(1:10, 0.05)))   # too short: never significant
})

test_that("inspectBreakpoint() reproduces every detected change", {
	skip_if_not_installed("EnvCpt")
	rows <- mkRows()
	ts <- aggregateClusterTimeSeries(rows)
	ctx <- breakpointContext(ts, rows)
	expect_gte(nrow(ctx$bp), 1)
	for (i in seq_len(nrow(ctx$bp))) {
		ins <- expect_silent(inspectBreakpoint(ctx, row = ctx$bp[i, ], envcpt = FALSE))
		expect_equal(ins$label, ctx$bp$label[i])
		expect_equal(ins$char$d_lo, ctx$bp$d_lo[i])
		expect_equal(ins$char$d_hi, ctx$bp$d_hi[i])
		expect_equal(ins$spatial$p_comp, ctx$bp$p_comp[i])
		expect_equal(ins$common, ctx$bp$common[i])
		expect_equal(sum(ins$flags$comp_part), ins$spatial$comp_change)
	}
})

test_that("manual inspection works where nothing was detected and on a range", {
	skip_if_not_installed("EnvCpt")
	rows <- mkRows()
	ts <- aggregateClusterTimeSeries(rows)
	ctx <- breakpointContext(ts, rows)
	ins <- inspectBreakpoint(ctx, cluster = 2, variable = "mean_len", at = "2000-06-01")
	expect_s3_class(ins, "bpInspection")
	expect_equal(ins$mode, "manual")
	ins2 <- inspectBreakpoint(ctx, cluster = 1, variable = "CPUE", at = c("2003-01-01", "2006-12-31"))
	expect_equal(nrow(ins2$scan), 48)
	expect_lte(abs(as.numeric(ins2$char$start_date - as.Date("2004-12-15"))), 100)
	pl <- plotInspection(ins2, ctx)
	expect_true(all(vapply(pl[c("series", "local", "profile", "shares", "context")], inherits, logical(1), "ggplot")))
	expect_error(inspectBreakpoint(ctx, cluster = 1, variable = "CPUE"), "Give `row`")
})

test_that("the move of flag B is mirrored by cluster 2", {
	skip_if_not_installed("EnvCpt")
	rows <- mkRows()
	ts <- aggregateClusterTimeSeries(rows)
	ctx <- breakpointContext(ts, rows)
	ins <- inspectBreakpoint(ctx, cluster = 1, variable = "CPUE", at = "2009-12-15", envcpt = FALSE)
	expect_equal(as.character(ins$mirror$cluster), "2")
	expect_gte(ins$mirror$share, 0.5)
	expect_gt(ins$flags$d_share[ins$flags$flag == "B"], 0.3)   # B appears in cluster 1
})
