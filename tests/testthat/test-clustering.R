## Tests for R/clustering.R.

two_blobs <- function(n_per_blob = 30) {
	rbind(matrix(rnorm(n_per_blob * 2, sd = 0.3), ncol = 2),
		  matrix(rnorm(n_per_blob * 2, mean = 4, sd = 0.3), ncol = 2))
}

# ---- selectScenario() / nSpeciesComponents() ---------------------------------

test_that("selectScenario() maps tokens to columns, without duplicates", {
	expect_equal(selectScenario("yba_hbf"),
				 c("yft_fraction", "bet_fraction", "alb_fraction", "mean_hbf"))
	expect_equal(selectScenario("sp_hbf_latabs"),
				 c("yft_fraction", "bet_fraction", "alb_fraction", "oth_fraction",
				   "mean_hbf", "latAbs"))
	expect_equal(selectScenario("yft_yba"), c("yft_fraction", "bet_fraction", "alb_fraction"))
	expect_equal(selectScenario("lat_month"), c("latitude", "cosmonth", "sinmonth"))
	expect_warning(out <- selectScenario("yft_nope"), "nope")
	expect_equal(out, "yft_fraction")
})

test_that("nSpeciesComponents() is the species count minus one, floored", {
	expect_equal(nSpeciesComponents("sp_hbf"), 3)
	expect_equal(nSpeciesComponents("yba_lat"), 2)
	expect_message(n <- nSpeciesComponents("yft_bet_hbf"), "min_components")
	expect_equal(n, 2)
	expect_equal(suppressMessages(nSpeciesComponents("hbf_lat", min_components = 3)), 3)
})

# ---- recalcSpeciesFractions() ---------------------------------------------------

test_that("recalcSpeciesFractions() rescales to the scenario's species", {
	df <- data.frame(yft_n = c(2, 0), bet_n = c(1, 0), alb_n = c(1, 0),
					 skj_n = c(4, 0), oth_n = c(2, 0),
					 skj_fraction = c(0.9, 0.9))
	yb <- recalcSpeciesFractions(df, "yft_bet_hbf")
	expect_equal(yb$total_n, c(3, 0))
	expect_equal(yb$yft_fraction, c(2 / 3, 0))   # 0 when no catch, not NaN
	expect_equal(yb$bet_fraction, c(1 / 3, 0))
	expect_equal(yb$skj_fraction, c(0.9, 0.9))   # untouched

	## sp: skipjack is counted with "other"
	sp <- recalcSpeciesFractions(df, "sp_hbf")
	expect_equal(sp$total_n[1], 10)
	expect_equal(sp$oth_fraction[1], (2 + 4) / 10)
	expect_equal(sp$yft_fraction[1] + sp$bet_fraction[1] + sp$alb_fraction[1] +
				 	sp$oth_fraction[1], 1)
	expect_equal(sp$oth_n, df$oth_n)             # counts untouched

	expect_error(recalcSpeciesFractions(df[, -1], "yft_bet"), "yft_n")
	expect_message(out <- recalcSpeciesFractions(df, "hbf_lat"), "no species")
	expect_identical(out, df)
})

# ---- pcaScale() / applyScaleParams() -------------------------------------------

test_that("pcaScale() z-scores, stores its parameters and centres constant columns", {
	df <- data.frame(a = c(1, 2, 3, 6), b = c(5, 5, 5, 5))
	out <- pcaScale(df)
	expect_equal(out$a, as.numeric(scale(df$a)))
	expect_equal(out$b, rep(0, 4))
	sp <- attr(out, "scale_params")
	expect_equal(sp$variable, c("a", "b"))
	expect_equal(sp$mean, c(3, 5))
	expect_equal(sp$sd, c(stats::sd(df$a), 1))

	arc <- pcaScale(data.frame(p = c(0, 0.25, 1)), method = "arcsine")
	expect_equal(arc$p, asin(sqrt(c(0, 0.25, 1))))
	expect_null(attr(arc, "scale_params"))

	expect_error(pcaScale(df, method = "minmax"), "Unknown scaling method")
	expect_error(pcaScale(data.frame(a = 1:3, f = letters[1:3])), "f")
})

test_that("applyScaleParams() reproduces pcaScale() on the fit rows and reuses it on new rows", {
	fit <- data.frame(a = c(1, 2, 3, 6), b = c(10, 20, 30, 40))
	sc  <- pcaScale(fit)
	sp  <- attr(sc, "scale_params")
	expect_equal(applyScaleParams(fit, sp)$a, sc$a)
	new <- applyScaleParams(data.frame(a = 3, b = 0, other = "x"), sp)
	expect_equal(new$a, 0)
	expect_equal(new$b, (0 - 25) / stats::sd(fit$b))
	expect_equal(new$other, "x")
	expect_error(applyScaleParams(data.frame(a = 1), sp), "missing column")
	expect_error(applyScaleParams(fit, NULL), "scale_params")
})

# ---- customPCA() -----------------------------------------------------------------

test_that("customPCA() keeps enough components for the threshold, at least 2", {
	set.seed(1)
	z <- rnorm(200)
	df <- pcaScale(data.frame(a = z + rnorm(200, sd = 0.01), b = z + rnorm(200, sd = 0.01),
							  c = rnorm(200), d = rnorm(200)))
	res <- customPCA(df, variance_threshold = 0.99, print_it = FALSE)
	expect_s3_class(res, "prcomp")
	cum <- cumsum(res$sdev^2) / sum(res$sdev^2)
	expect_equal(res$no_var, max(2, which(cum >= 0.99)[1]))
	expect_equal(customPCA(df, variance_threshold = 0.1, print_it = FALSE)$no_var, 2)

	plots <- customPCA(df, print_it = FALSE, return_print = TRUE)
	expect_named(plots, c("loadings", "variances_plot", "loading_plot", "scatter_plot"))
	expect_s3_class(plots$loading_plot, "ggplot")
})

# ---- customKmeans() ----------------------------------------------------------------

test_that("customKmeans() returns the documented structure on toy data", {
	set.seed(1)
	x <- two_blobs()
	expect_output(res <- customKmeans(x, max_k = 5, random_set = 10, iter_max = 50,
									  nstart = 2, print_it = FALSE), "Best K")
	expect_named(res, c("kmeans", "gap_stat", "plot"))
	expect_s3_class(res$kmeans, "kmeans")
	expect_s3_class(res$plot, "ggplot")
	expect_equal(nrow(res$gap_stat), 5)
	expect_equal(nrow(res$kmeans$centers), selectK(as.matrix(res$gap_stat)))
})

test_that("customKmeans() falls back to cluster::clusGap() when d.power != 2", {
	set.seed(1)
	expect_output(res <- customKmeans(two_blobs(15), max_k = 3, random_set = 5, d.power = 1,
									  print_it = FALSE), "Best K")
	expect_true(all(c("logW", "E.logW", "gap", "SE.sim") %in% names(res$gap_stat)))
})

# ---- assignClusters() --------------------------------------------------------------

test_that("assignClusters() returns the row of the nearest centroid", {
	centroids <- data.frame(PC1 = c(0, 10), PC2 = c(0, 10))
	pts <- data.frame(PC1 = c(1, 9, 4, 6), PC2 = c(-1, 11, 4, 6), other = 1:4)
	expect_equal(unname(assignClusters(pts, centroids)), c(1, 2, 1, 2))

	set.seed(2)
	x  <- two_blobs()
	km <- stats::kmeans(x, 2, nstart = 5)
	d  <- setNames(as.data.frame(x), c("PC1", "PC2"))
	expect_equal(unname(assignClusters(d, setNames(as.data.frame(km$centers), c("PC1", "PC2")))),
				 unname(km$cluster))
})

# ---- summaryClusters() / orderClusters() -----------------------------------------

test_that("summaryClusters() summarises catch, CPUE and dates per cluster", {
	df <- data.frame(
		cluster  = c(2, 2, 2, 1),
		ymd      = as.Date(c("2000-01-15", "2000-01-15", "2000-02-15", "2001-05-15")),
		yft_n    = c(10, 10, 30, 5),
		bet_n    = c(0, 10, 0, 5),
		alb_n    = c(0, 0, 10, 0),
		E        = c(1, 1, 2, 1),
		mean_len = c(100, 110, 120, 90),
		mean_hbf = c(10, 20, 30, 40),
		latitude = c(-10, -10, -20, 5))
	s <- summaryClusters(df)
	expect_equal(s$cluster, c(1, 2))
	c2 <- s[s$cluster == 2, ]
	expect_equal(c2$CPUE_min, 10)          # January: 20 / 2
	expect_equal(c2$CPUE_max, 15)          # February: 30 / 2
	expect_equal(c2$yft_n, 50)
	expect_equal(c2$total_n, 70)
	expect_equal(c2$len_min, 105)          # mean of the January rows
	expect_equal(c2$HBF_max, 30)
	expect_equal(c2$lat_min, -20)
	expect_equal(c2$ymd_min, as.Date("2000-01-15"))
	expect_equal(c2$ymd_max, as.Date("2000-02-15"))
	expect_equal(sum(s$yft_catch), 100)
	expect_equal(c2$yft_frac, 50 / 70 * 100)
	expect_false(any(c("len_min", "HBF_min", "lat_min") %in%
					 	names(summaryClusters(df[, 1:6]))))
})

test_that("orderClusters() ranks clusters by yellowfin weight or fraction", {
	df <- data.frame(cluster = c(1, 1, 2, 3),
					 yft_n = c(10, 5, 40, 2), total_n = c(20, 20, 50, 10),
					 yft_w = c(1, 1, 3, 9))
	by_w <- orderClusters(df)
	expect_equal(attr(by_w, "old_ids"), c(3, 2, 1))
	expect_equal(as.integer(as.character(by_w$cluster)), c(3, 3, 2, 1))
	expect_s3_class(by_w$cluster, "factor")
	expect_equal(attr(orderClusters(df, by = "yft_frac"), "old_ids"), c(2, 1, 3))
	expect_error(orderClusters(df[, c("cluster", "yft_n")]), "yft_w")
})

# ---- customTheme() ----------------------------------------------------------------

test_that("customTheme() is a ggplot2 theme with the legend at the bottom", {
	th <- customTheme(text_size = 9)
	expect_s3_class(th, "theme")
	expect_equal(th$legend.position, "bottom")
	expect_equal(th$text$size, 9)
})
