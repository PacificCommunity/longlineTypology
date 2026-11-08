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
	## With k.min > 1, a row's POSITION in Tab no longer equals its k, so
	## selectK() must use attr(Tab, "k") rather than cluster::maxSE()'s
	## return value directly.
	set.seed(1)
	x   <- two_blobs()
	Tab <- fastClusGap(x, K.max = 6, k.min = 3, B = 10, verbose = FALSE)

	expect_identical(attr(Tab, "k"), 3:6)
	nc <- selectK(Tab)
	expect_true(nc %in% 3:6)
})

test_that("fastClusGap() matches the W_k of cluster::clusGap() at k = 1", {
	## At k = 1 there is no clustering, so logW is deterministic and must
	## equal clusGap()'s value exactly (same 0.5 factor).
	set.seed(2)
	x   <- two_blobs()
	Tab <- fastClusGap(x, K.max = 2, B = 5, verbose = FALSE)
	ref <- cluster::clusGap(x, FUNcluster = stats::kmeans, K.max = 2, B = 5,
							d.power = 2, verbose = FALSE)$Tab
	expect_equal(unname(Tab[1, "logW"]), unname(ref[1, "logW"]))
})

test_that("fastClusGap() gives the same result whatever ncores and restores the RNG", {
	skip_on_os("windows")
	x <- two_blobs()
	set.seed(3); a <- fastClusGap(x, K.max = 4, B = 6, verbose = FALSE); r1 <- runif(1)
	set.seed(3); b <- fastClusGap(x, K.max = 4, B = 6, ncores = 2, verbose = FALSE); r2 <- runif(1)
	expect_equal(a[, "E.logW"], b[, "E.logW"])
	expect_equal(r1, r2)
})

test_that("selectK() falls back to rownames when the k attribute is missing", {
	set.seed(1)
	Tab <- fastClusGap(two_blobs(), K.max = 5, k.min = 2, B = 10, verbose = FALSE)
	k_attr <- selectK(Tab)
	attr(Tab, "k") <- NULL
	expect_identical(selectK(Tab), k_attr)
})

test_that("checkWkIdentity() confirms the identity and Wk_fast() matches it", {
	set.seed(4)
	x <- matrix(rnorm(300), ncol = 3)
	expect_output(expect_true(checkWkIdentity(x, k = 3)), "identity holds")
	expect_equal(longlineTypology:::Wk_fast(x, 1L), 0.5 * sum(sweep(x, 2, colMeans(x))^2))
})
