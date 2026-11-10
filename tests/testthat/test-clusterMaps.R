mapData <- function() {
	set.seed(1)
	n <- 400
	d <- data.frame(longitude = runif(n, 140, 200), latitude = runif(n, -30, 10),
					E = runif(n, 1, 10), yft_n = rpois(n, 5),
					mean_len = runif(n, 90, 140), mean_hbf = runif(n, 5, 30))
	list(data = d, cluster = sample(1:3, n, replace = TRUE))
}

test_that("clusterColours gives k distinct colours", {
	expect_length(clusterColours(4), 4)
	expect_named(clusterColours(3), as.character(1:3))
	expect_length(unique(clusterColours(12)), 12)
})

test_that("plotDominance shows one cluster per cell", {
	m <- mapData()
	p <- plotDominance(m$data, m$cluster, cell_size = 10)
	expect_s3_class(p, "ggplot")
	expect_silent(ggplot2::ggplot_build(p))
	cells <- p$data[, c("lon_cell", "lat_cell")]
	expect_false(anyDuplicated(cells) > 0)
})

test_that("plotClusterShare gives shares summing to 1 per cell", {
	m <- mapData()
	p <- plotClusterShare(m$data, m$cluster, cell_size = 10)
	expect_s3_class(p, "ggplot")
	tot <- tapply(p$data$share, paste(p$data$lon_cell, p$data$lat_cell), sum)
	expect_equal(unname(as.numeric(tot)), rep(1, length(tot)))
})

test_that("plotClusterMap draws maps with and without profiles", {
	m <- mapData()
	p <- plotClusterMap(m$data, m$cluster, "CPUE", cell_size = 10, profiles = FALSE)
	expect_s3_class(p, "ggplot")
	expect_silent(ggplot2::ggplot_build(p))
	skip_if_not_installed("patchwork")
	for (v in c("CPUE", "mean_len", "mean_hbf")) {
		pp <- plotClusterMap(m$data, m$cluster, v, cell_size = 10, clusters = 1:2)
		expect_s3_class(pp, "patchwork")
	}
	expect_error(plotClusterMap(m$data[, -3], m$cluster, "CPUE"), "E")
})
