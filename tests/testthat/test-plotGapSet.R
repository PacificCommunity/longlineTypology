## Gap tables of three replicates on toy data with two well-separated groups.
gapTables <- function() {
	set.seed(1)
	x <- rbind(matrix(rnorm(60, sd = 0.3), ncol = 2),
			   matrix(rnorm(60, mean = 4, sd = 0.3), ncol = 2))
	lapply(1:3, function(i) as.data.frame(fastClusGap(x, K.max = 5, B = 10, verbose = FALSE)))
}

test_that("plotGapSet() reports one selected k per replicate, as selectK() does", {
	gs <- gapTables()
	grDevices::pdf(NULL)
	on.exit(grDevices::dev.off())
	out <- plotGapSet(gs)
	expect_length(out$selected, 3)
	expect_equal(sum(out$table), 3)
	expect_equal(out$selected, vapply(gs, function(t) selectK(as.matrix(t)), integer(1)))
})

test_that("plotGapSet() writes a png when given a file", {
	f <- tempfile(fileext = ".png")
	plotGapSet(gapTables(), file = f)
	expect_true(file.exists(f))
	expect_gt(file.size(f), 0)
})
