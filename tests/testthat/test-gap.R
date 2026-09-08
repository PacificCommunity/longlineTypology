## Toy two-blob data: two well-separated clusters, so the true K = 2 is
## unambiguous and these tests stay fast and deterministic.
two_blobs <- function(n_per_blob = 30) {
	rbind(
		matrix(rnorm(n_per_blob * 2, sd = 0.3), ncol = 2),
		matrix(rnorm(n_per_blob * 2, mean = 4, sd = 0.3), ncol = 2)
	)
}

test_that("fastClusGap() returns a table selectK() can read", {
	set.seed(1)
	x   <- two_blobs()
	Tab <- fastClusGap(x, K.max = 5, B = 10, verbose = FALSE)

	expect_true(all(c("logW", "E.logW", "gap", "SE.sim") %in% colnames(Tab)))
	expect_identical(attr(Tab, "k"), 1:5)

	nc <- selectK(Tab)
	expect_true(nc >= 1 && nc <= 5)
})

test_that("selectK() reads k from the table's attribute, not row position", {
	## Regression test for the bug documented in fastClusGap.R: with
	## k.min > 1, a row's POSITION in Tab no longer equals its k, so
	## selectK() must use attr(Tab, "k") rather than cluster::maxSE()'s
	## return value directly.
	set.seed(1)
	x   <- two_blobs()
	Tab <- fastClusGap(x, K.max = 6, k.min = 3, B = 10, verbose = FALSE)

	expect_identical(attr(Tab, "k"), 3:6)
	nc <- selectK(Tab)
	expect_true(nc %in% 3:6)
})

test_that("customKmeans() returns the documented structure on toy data", {
	set.seed(1)
	x   <- two_blobs()
	res <- customKmeans(x, max_k = 5, random_set = 10, iter_max = 50,
						nstart = 2, print_it = FALSE)

	expect_named(res, c("kmeans", "gap_stat", "plot"))
	expect_s3_class(res$kmeans, "kmeans")
	expect_s3_class(res$plot, "ggplot")
	expect_true(nrow(res$kmeans$centers) >= 1)
})

test_that("plotGapSet() runs without error and reports one selected k per replicate", {
	set.seed(1)
	gap_stat <- lapply(1:3, function(i) {
		customKmeans(two_blobs(), max_k = 5, random_set = 10, iter_max = 50,
					nstart = 2, print_it = FALSE)$gap_stat
	})

	out <- plotGapSet(gap_stat, file = NULL)
	expect_length(out$selected, 3)
})
