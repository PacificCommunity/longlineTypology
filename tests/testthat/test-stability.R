## Small clustering of the mock data, shared by the tests of this file.
mock_run <- local({
	mock <- system.file("extdata", "mockdata", "LL_mock.RDS", package = "longlineTypology")
	diag <- runClustering(mock, "sp_hbf_latabs", sample_no = 100, kmeans_set = 4,
						  max_k = 8, B = 10, nstart = 5, verbose = FALSE)
	dat  <- clusterData(diag)
	fr   <- fitRows(diag, dat)
	sc   <- pcaScores(diag, dat)[fr, , drop = FALSE]
	k    <- names(diag$kmeans_full)[1]
	list(diag = diag, data = dat, fit_rows = fr, scores = sc, k = k,
		 cluster = clusterLabels(diag, k)[fr])
})

test_that(".hungarian finds the optimal assignment", {
	set.seed(1)
	perms <- function(n) if (n == 1) list(1L) else
		do.call(c, lapply(seq_len(n), function(i) lapply(perms(n - 1), function(p)
			c(i, setdiff(seq_len(n), i)[p]))))
	for (n in 1:6) {
		cost <- matrix(runif(n * n), n)
		best <- min(vapply(perms(n), function(p) sum(cost[cbind(seq_len(n), p)]), numeric(1)))
		perm <- longlineTypology:::.hungarian(cost)
		expect_setequal(perm, seq_len(n))
		expect_equal(sum(cost[cbind(seq_len(n), perm)]), best)
	}
})

test_that("centroidVariability matches replicate centroids to the full fit", {
	m  <- mock_run
	cv <- centroidVariability(m$diag, m$k, m$scores)
	n_pc <- m$diag$n_components
	if (cv$n_replicates > 0) {
		expect_equal(nrow(cv$stats), as.integer(m$k) * n_pc)
		expect_true(all(cv$stats$nrmse >= 0))
		expect_true(all(diff(cv$stats$nrmse) <= 0))
	}
	expect_equal(sum(cv$points$source == "full fit"), as.integer(m$k))
	expect_s3_class(plotCentroidVariability(cv), "ggplot")

	## A replicate identical to the full fit, with shuffled rows, has rmse 0
	d <- m$diag
	ref <- d$kmeans_full[[m$k]]$centers
	d$kmeans_res <- list(list(centers = ref[rev(seq_len(nrow(ref))), , drop = FALSE]))
	cv0 <- centroidVariability(d, m$k, m$scores)
	expect_equal(max(cv0$stats$rmse), 0)
})

test_that("jaccardStability is 1 when the replicates are the full fit", {
	m <- mock_run
	js <- jaccardStability(m$diag, m$scores, m$cluster)
	expect_equal(nrow(js$values), length(unique(m$cluster)) * length(m$diag$kmeans_res))
	expect_true(all(js$summary$jaccard >= js$summary$jaccard_min))
	expect_true(all(js$values$jaccard > 0 & js$values$jaccard <= 1))

	d <- m$diag
	d$kmeans_res <- list(list(centers = d$kmeans_full[[m$k]]$centers))
	expect_equal(jaccardStability(d, m$scores, m$cluster)$summary$jaccard,
				 rep(1, length(unique(m$cluster))))
	expect_s3_class(plotJaccard(js), "ggplot")
})

test_that("temporalStability counts coverage and internal gaps", {
	year    <- c(2000, 2000, 2001, 2003, 2003, 2004, 2004)
	cluster <- c(1,    2,    1,    1,    2,    1,    2)
	ts <- temporalStability(year, cluster)
	s2 <- ts$summary[ts$summary$cluster == 2, ]
	## cluster 2: present 2000, 2003, 2004; absent 2001, 2002 (2002 has no rows)
	expect_equal(s2$coverage, 3 / 5)
	expect_equal(s2$max_gap, 2)
	s1 <- ts$summary[ts$summary$cluster == 1, ]
	expect_equal(s1$max_gap, 1)   # 2002 only
	expect_equal(sum(ts$share$share[ts$share$year == 2001]), 1)
	expect_equal(sum(ts$share$share[ts$share$year == 2002]), 0)
	expect_s3_class(plotTemporalShare(ts), "ggplot")
})

test_that("silhouetteStats matches cluster::silhouette and adds the null on request", {
	m  <- mock_run
	ss <- silhouetteStats(m$scores, m$cluster)
	ref <- cluster::silhouette(as.integer(factor(m$cluster)), stats::dist(m$scores))
	expect_equal(ss$overall, mean(ref[, "sil_width"]))
	expect_null(ss$null)

	set.seed(4); a <- runif(1)
	set.seed(4)
	smp <- subsampleReplicates(m$diag$sample_indices, 50, seed = 1)
	ssn <- silhouetteStats(m$scores, m$cluster, smp, n_null = 3, seed = 1)
	expect_equal(runif(1), a)
	expect_named(ssn$null, c("mean", "sd", "z", "n", "size"))
	expect_equal(ssn$null$size, 50)
	expect_s3_class(plotSilhouette(ssn), "ggplot")
})
