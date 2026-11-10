## Path to an RDS file of mockLonglineData(...), written once per test run.
mockPath <- local({
	cache <- list()
	function(breakpoints = FALSE) {
		key <- as.character(breakpoints)
		if (is.null(cache[[key]])) {
			f <- tempfile(fileext = ".RDS")
			saveRDS(mockLonglineData(breakpoints = breakpoints), f)
			cache[[key]] <<- f
		}
		cache[[key]]
	}
})
