## Tests for R/mockData.R.

test_that("mockLonglineData() has the seafisheries columns and is reproducible", {
	set.seed(9); a <- runif(1)
	set.seed(9)
	m <- mockLonglineData()
	expect_equal(runif(1), a)                            # caller's RNG restored
	need <- c("flag", "latitude", "longitude", "ymd", "E", "mean_hbf", "mean_len",
			  paste0(c("yft", "bet", "alb", "skj", "oth"), "_n"), "yft_w",
			  "latAbs", "hbf_outlier", "hbf_cat", "len_source", "mock_strategy")
	expect_true(all(need %in% names(m)))
	expect_s3_class(m$ymd, "Date")
	expect_equal(nrow(m), 4 * 12 * length(1985:2022) * 3)
	expect_identical(m, mockLonglineData())
	expect_false(identical(m, mockLonglineData(seed = 5)))
	expect_null(attr(m, "events"))
	expect_true(all(m$longitude >= 0 & m$longitude < 360))
	expect_equal(m$latAbs, abs(m$latitude))
	expect_setequal(unique(m$hbf_outlier), c("valid", "outlier"))
	expect_true(all(grepl("^[A-Z]{2}_V[1-5]$", m$vesselname)))
})

test_that("planted changes act on the rows they target", {
	m  <- mockLonglineData(breakpoints = TRUE, n_per_month = 10)
	ev <- attr(m, "events")
	expect_identical(ev, mockEvents())
	expect_setequal(ev$label, c("unexplained", "transient", "gradual", "relocation", "switch", "common"))

	## switch: the switching flags fish the albacore strategy afterwards
	sw <- ev[ev$type == "switch", ]
	fl <- strsplit(sw$flags, ",")[[1]]
	after <- m$flag %in% fl & m$ymd >= sw$start
	expect_true(all(m$mock_strategy[after] == "subtropical_albacore"))
	expect_true(all(m$mock_strategy[m$flag %in% fl & m$ymd < sw$start - 31] == "temperate_mixed"))

	## relocation: the relocated flags fish further east
	rl <- ev[ev$type == "relocation", ]
	fl <- strsplit(rl$flags, ",")[[1]]
	keep <- m$flag %in% fl & m$mock_strategy == "temperate_mixed"
	before <- mean(m$longitude[keep & m$ymd < rl$start])
	after  <- mean(m$longitude[keep & m$ymd >= rl$start])
	expect_gt(after - before, 20)

	## step: CPUE of the targeted strategy drops by about the planted factor
	st <- ev[ev$type == "step", ]
	x  <- m[m$mock_strategy == "tropical_shallow" & abs(as.numeric(m$ymd - st$start)) < 365, ]
	cpue <- tapply(x$yft_n / x$E, x$ymd >= st$start, mean)
	expect_equal(unname(cpue[["TRUE"]] / cpue[["FALSE"]]), exp(st$size), tolerance = 0.25)
})

test_that("the planted changes are found by the clustering and breakpoint analysis", {
	skip_on_cran()
	skip_if_not_installed("EnvCpt")
	d  <- runClustering(mockPath(breakpoints = TRUE), "sp_hbf_latabs", sample_no = 1000,
						kmeans_set = 5, max_k = 8, B = 10, nstart = 5, verbose = FALSE)
	expect_true("4" %in% names(d$kmeans_full))
	bd <- breakpointData(d, 4)
	## each strategy is recovered as one cluster
	tab <- table(bd$cluster, bd$mock_strategy)
	expect_true(all(apply(tab, 1, max) / rowSums(tab) > 0.95))
	bp <- analyseClusterBreakpoints(aggregateClusterTimeSeries(bd), bd)
	bp$lab <- ifelse(bp$common, "common", bp$label)
	strat_cluster <- setNames(as.integer(rownames(tab))[apply(tab, 2, which.max)], colnames(tab))
	ev <- mockEvents()
	for (i in seq_len(nrow(ev))) {
		e <- ev[i, ]
		cl <- if (e$strategy == "all") unname(strat_cluster) else
			if (e$type == "switch") strat_cluster[["subtropical_albacore"]] else strat_cluster[[e$strategy]]
		## a change overlapping the planted period, give or take a year
		from <- e$start - 365; to <- e$start + 30.5 * e$duration + 365
		near <- bp[bp$cluster %in% cl & bp$end_date >= from & bp$start_date <= to, , drop = FALSE]
		expect_true(e$label %in% near$lab, info = paste(e$type, "->", e$label))
	}
})
