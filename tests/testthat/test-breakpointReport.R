## Tests for R/breakpointReport.R.

test_that("decomposeChange() splits a change exactly into composition and within", {
	before <- data.frame(flag = c("A", "A", "B"), mean_len = c(100, 110, 130))
	after  <- data.frame(flag = c("A", "B", "B", "C"), mean_len = c(105, 130, 140, 90))
	dc <- decomposeChange(before, after, "mean_len")
	expect_equal(dc$total, mean(after$mean_len) - mean(before$mean_len))
	expect_equal(dc$composition + dc$within, dc$total)
	tab <- dc$table
	expect_equal(tab$within[tab$flag == "C"], 0)        # only after: all composition
	expect_equal(tab$share_before[tab$flag == "A"], 2 / 3)
	expect_equal(tab$value_after[tab$flag == "B"], 135)

	## CPUE: effort shares and sum(yft_n) / sum(E)
	b <- data.frame(flag = c("A", "B"), yft_n = c(10, 30), E = c(10, 10))
	a <- data.frame(flag = c("A", "B"), yft_n = c(30, 30), E = c(30, 10))   # A: more effort, same CPUE
	dc <- decomposeChange(b, a, "CPUE")
	expect_equal(dc$total, 60 / 40 - 40 / 20)
	expect_equal(dc$within, 0)                           # flag CPUEs unchanged
	expect_equal(dc$composition, dc$total)
})

test_that("breakpointWindows() and changeAcrossClusters() use the windows of the change", {
	skip_if_not_installed("EnvCpt")
	rows <- mkRows()
	rows$year <- as.integer(format(rows$date, "%Y"))
	ctx <- breakpointContext(aggregateClusterTimeSeries(rows), rows)
	r   <- ctx$bp[ctx$bp$cluster == 1 & ctx$bp$variable == "CPUE", , drop = FALSE][1, ]
	ins <- inspectBreakpoint(ctx, row = r, envcpt = FALSE)
	w   <- breakpointWindows(rows, ins)
	mi  <- longlineTypology:::.bpMonthIdx
	s   <- ins$char$s_hat
	expect_true(all(mi(w$before$date) >= s - 12 & mi(w$before$date) < s))
	expect_true(all(w$before$cluster == 1))
	expect_equal(w$per_month$months[w$per_month$window == "before"], 12)
	expect_equal(w$per_month$empty_months[w$per_month$window == "before"], 0)
	expect_equal(nrow(breakpointWindows(rows, ins, cluster = NULL)$before),
				 sum(mi(rows$date) >= s - 12 & mi(rows$date) < s))

	ac <- changeAcrossClusters(rows, ins)
	expect_equal(ac$cluster, c("1", "2", "all"))
	expect_lt(ac$rel_change[1], -0.3)                    # the CPUE step of cluster 1
	expect_equal(ac$before[1], sum(w$before$yft_n) / sum(w$before$E))
})

test_that("plotClusterBreakpoints() draws with and without changes", {
	ts <- aggregateClusterTimeSeries(mkRows())
	expect_s3_class(plotClusterBreakpoints(ts, NULL), "ggplot")
	bp <- data.frame(cluster = 1, variable = "CPUE", start_date = as.Date("2005-01-15"),
					 end_date = as.Date("2005-03-15"), label = "unexplained", common = FALSE)
	p <- plotClusterBreakpoints(ts, bp)
	expect_s3_class(p, "ggplot")
	expect_silent(ggplot2::ggplot_build(p))
})

## Clustering of the mock data with planted changes, for the report-level tests
mock_bp_run <- function() {
	runClustering(mockPath(breakpoints = TRUE), "sp_hbf_latabs", sample_no = 1000, kmeans_set = 5, max_k = 8,
				  B = 10, nstart = 5, verbose = FALSE)
}

test_that("breakpointData() keeps the fit years, orders the clusters and drops imputed values", {
	d  <- mock_bp_run()
	k  <- "4"
	bd <- breakpointData(d, k)
	fy <- d$fit_years
	expect_true(all(bd$year >= fy[1] & bd$year <= fy[2]))
	expect_identical(bd$date, bd$ymd)
	expect_setequal(unique(bd$cluster), seq_len(as.integer(k)))
	expect_true(all(is.na(bd$mean_len[bd$len_source %in% "imputed"])))
	kept <- breakpointData(d, k, use_imputed = TRUE)
	expect_false(all(is.na(kept$mean_len[kept$len_source %in% "imputed"])))
	expect_error(breakpointReport(d, k = 99), "k must be one of")
})

test_that("breakpointReport() renders, with detailed and manual sections", {
	skip_on_cran()
	skip_if_not_installed("EnvCpt")
	skip_if_not_installed("rmarkdown")
	skip_if(!rmarkdown::pandoc_available() || Sys.which("pdflatex") == "", "no pandoc / LaTeX")
	d   <- mock_bp_run()
	out <- tempfile(fileext = ".pdf")
	breakpointReport(d, k = 4, output_file = out, detail = TRUE,
					 manual = list(list(cluster = 1, variable = "CPUE", at = "2003-06-15")))
	expect_true(file.exists(out))
})
