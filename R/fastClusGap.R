## fastClusGap.R -- exact O(n) replacement for cluster::clusGap when d.power = 2
##
## WHY. seafisheries::customKmeans() calls cluster::clusGap(..., d.power = 2).
## clusGap computes W_k inside its internal W.k() as
##     0.5 * sum_r ( sum(dist(X_r)^d.power) / n_r )
## i.e. it builds the FULL PAIRWISE DISTANCE MATRIX of every cluster, for every
## k in 1:K.max, for the real data and for each of the B reference datasets.
## That is O(n^2) in both time and memory. The k = 1 call is the worst: one
## dist() over all n rows, which is n(n-1)/2 doubles -- 1.5 GB at n = 20,000 and
## 149 GB at n = 200,000. This, not the clustering, is what makes the gap step
## take days, and it is why a 1% subsample was the largest that would run.
##
## THE IDENTITY. For squared Euclidean distance (d.power = 2),
##     sum_{i<j in r} ||x_i - x_j||^2  ==  n_r * sum_{i in r} ||x_i - mu_r||^2
## so
##     W_k  ==  0.5 * sum_r SS_r  ==  0.5 * kmeans(...)$tot.withinss
## which kmeans() has already computed and returned for free. Verified to 12
## significant digits against clusGap's own W.k(). The factor 0.5 is common to
## logW and E.logW and cancels in gap = E.logW - logW, so results are identical
## up to K-means random starts.
##
## NOTE. The identity holds ONLY for d.power = 2. clusGap's default is
## d.power = 1 (plain Euclidean), which has no closed form. d.power = 2 is the
## version that matches Tibshirani et al. (2001) for K-means, and is what
## customKmeans() already passes, so nothing changes statistically.
##
## Cost, measured (4 dims, K.max = 20, one core):
##   cluster::clusGap  n = 20,000  ~20 min per call at B = 100 ; O(n^2) beyond
##   fastClusGap       n = 20,000  ~2.3 min ; n = 100,000  ~12 min ; linear in n
##   speed-up at n = 20,000: 9x, and it grows proportionally to n thereafter.

Wk_fast <- function(X, k, nstart = 1L, iter.max = 30L) {
	if (k == 1L) 0.5 * sum(sweep(X, 2L, colMeans(X))^2)
	else 0.5 * stats::kmeans(X, k, nstart = nstart, iter.max = iter.max)$tot.withinss
}

#' Drop-in replacement for cluster::clusGap() when d.power = 2
#'
#' @param x        numeric matrix or data frame (e.g. the retained PCA scores)
#' @param K.max    largest k to test
#' @param B        number of uniform reference datasets (clusGap's B / random_set)
#' @param nstart   passed to kmeans()
#' @param iter.max passed to kmeans()
#' @param k.min    smallest k to test. clusGap always starts at 1; k = 1 is the
#'                 single most expensive call there and is usually uninformative,
#'                 so k.min = 2 is a cheap saving. SEE THE WARNING ON selectK().
#' @param ncores   number of cores for the reference-dataset loop, which is
#'                 embarrassingly parallel. Default 1 (serial). Uses
#'                 parallel::mclapply(), so Unix-alikes only; on Windows it
#'                 warns and runs serially. Results do NOT depend on ncores:
#'                 each reference sample draws from its own pre-generated seed,
#'                 so serial and parallel runs are bit-identical.
#' @param verbose  progress tracing in the style of cluster::clusGap(), plus
#'                 an ETA printed after the first bootstrap sample. Defaults
#'                 to interactive(), matching clusGap's own default.
#' @return matrix with columns logW, E.logW, gap, SE.sim (as clusGap()$Tab) plus
#'         an attribute "k" giving the k value of each row.
fastClusGap <- function(x, K.max, B = 100L, nstart = 1L, iter.max = 30L,
						k.min = 1L, ncores = 1L, verbose = interactive()) {
	x <- as.matrix(x); n <- nrow(x)
	stopifnot(k.min >= 1L, K.max >= k.min, n > K.max)

	## Capture the caller's RNG state FIRST -- the main k loop below already
	## consumes randomness through kmeans(nstart = ), so capturing later would
	## restore to the wrong point. One seed per reference sample is drawn here,
	## which is what makes the result independent of `ncores`: worker b builds
	## the same reference dataset serially or in a forked child.
	if (!exists(".Random.seed", envir = .GlobalEnv)) stats::runif(1)
	old_seed <- get(".Random.seed", envir = .GlobalEnv)
	on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
	seeds <- sample.int(.Machine$integer.max, B)
	ks <- k.min:K.max
	el <- function() proc.time()[["elapsed"]]
	t0 <- el()
	say <- function(...) if (verbose) { cat(..., sep = ""); utils::flush.console() }
	## seconds below two minutes, minutes above; hours above two hours
	dur <- function(sec) if (sec < 120) sprintf("%.0f s", sec)
						 else if (sec < 7200) sprintf("%.1f min", sec / 60)
						 else sprintf("%.1f h", sec / 3600)

	## Progress format follows cluster::clusGap(), with two additions: one dot
	## per k in the main loop rather than a single ".." for all of them, and an
	## ETA after the first bootstrap sample.
	say("Clustering k = ", k.min, ",", min(k.min + 1L, K.max),
		",..., K.max (= ", K.max, "): ")
	logW <- vapply(ks, function(k) {
		v <- log(Wk_fast(x, k, nstart, iter.max)); say("."); v
	}, 0)
	t_main <- el() - t0
	say(" done (", dur(t_main), ")\n")

	## reference distribution: clusGap's spaceH0 = "scaledPCA", reproduced exactly
	xs  <- scale(x, center = TRUE, scale = FALSE)
	m.x <- rep(attr(xs, "scaled:center"), each = n)
	V   <- svd(xs, nu = 0)$v
	rng <- apply(xs %*% V, 2L, range)

	one_ref <- function(b) {
		set.seed(seeds[b])
		z1 <- apply(rng, 2L, function(M) stats::runif(n, M[1L], M[2L]))
		z  <- tcrossprod(z1, V) + m.x
		vapply(ks, function(k) log(Wk_fast(z, k, nstart, iter.max)), 0)
	}

	ncores <- max(1L, as.integer(ncores))
	if (ncores > 1L && .Platform$OS.type == "windows") {
		warning("ncores > 1 needs fork(); running serially on Windows.")
		ncores <- 1L
	}

	say("Bootstrapping, b = 1,2,..., B (= ", B, ")",
		if (ncores > 1L) paste0(" on ", ncores, " cores") else
			"  [one \".\" per sample]", ":\n")
	t1 <- el()

	if (ncores == 1L) {
		logWks <- matrix(0, B, length(ks))
		for (b in seq_len(B)) {
			logWks[b, ] <- one_ref(b)
			say(".")
			if (b %% 50 == 0) say(" ", b, "\n")
			if (b == 1L) {
				eta <- (el() - t1) * B
				if (eta > 20) say("\n  [~", dur(eta), " to go, ~", dur(t_main + eta),
								  " total for this call]\n  ")
			}
		}
		if (B %% 50 != 0) say(" ", B, "\n")
	} else {
		## Sample 1 runs here so its cost can drive the ETA, then the rest fork.
		first <- one_ref(1L)
		t_one <- el() - t1
		eta   <- t_one * ceiling((B - 1L) / ncores)
		say("  sample 1 took ", dur(t_one), "; ~", dur(eta), " for the other ",
			B - 1L, " across ", ncores, " cores\n")
		rest <- parallel::mclapply(seq_len(B)[-1L], one_ref, mc.cores = ncores)
		bad  <- vapply(rest, function(r) inherits(r, "try-error") ||
									   length(r) != length(ks), logical(1))
		if (any(bad)) stop("mclapply: ", sum(bad), " reference sample(s) failed; ",
						   "re-run with ncores = 1 to see the error.")
		logWks <- rbind(first, do.call(rbind, rest))
		say("  done\n")
	}
	say("Total ", dur(el() - t0), "  (n = ", n, ", d = ", ncol(x),
		", K.max = ", K.max, ", B = ", B, ", nstart = ", nstart,
		 if (ncores > 1L) paste0(", ncores = ", ncores) else "", ")\n")
	E.logW <- colMeans(logWks)
	SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2L, stats::var))
	out <- cbind(logW = logW, E.logW = E.logW, gap = E.logW - logW, SE.sim = SE.sim)
	rownames(out) <- ks
	attr(out, "k") <- ks
	out
}

#' Apply the Tibshirani 2001 SEmax rule and return the NUMBER OF CLUSTERS
#'
#' WARNING, and the reason this wrapper exists: cluster::maxSE() returns a
#' POSITION in the vector it is given, not a k. customKmeans() gets away with
#'   nc <- cluster::maxSE(gap_stat$Tab[, "gap"], ...)
#' only because clusGap always starts at k = 1, so position == k. If you set
#' k.min = 2 to skip the expensive k = 1 call, position and k differ by one and
#' that line silently returns the wrong number of clusters. Always go through
#' selectK() rather than calling maxSE() on the table directly.
selectK <- function(Tab, SE.factor = 1) {
	ks <- attr(Tab, "k"); if (is.null(ks)) ks <- as.integer(rownames(Tab))
	i  <- cluster::maxSE(f = Tab[, "gap"], SE.f = Tab[, "SE.sim"],
						 method = "Tibs2001SEmax", SE.factor = SE.factor)
	ks[i]
}

#' Self-test: confirm the identity on YOUR data before trusting any of this.
#' Returns TRUE invisibly, or stops. Cheap only for small n (it calls dist()).
checkWkIdentity <- function(X, k = 5L, tol = 1e-8) {
	X <- as.matrix(X)
	if (nrow(X) > 5000L) X <- X[sample.int(nrow(X), 5000L), , drop = FALSE]
	km <- stats::kmeans(X, k, nstart = 1L, iter.max = 50L)
	ii <- seq_len(nrow(X))
	ref <- 0.5 * sum(vapply(split(ii, km$cluster), function(I) {
		xs <- X[I, , drop = FALSE]; sum(stats::dist(xs)^2 / nrow(xs))
	}, 0))
	stopifnot(isTRUE(all.equal(ref, 0.5 * km$tot.withinss, tolerance = tol)))
	cat("W_k identity holds: clusGap W.k =", ref,
		"  0.5*tot.withinss =", 0.5 * km$tot.withinss, "\n")
	invisible(TRUE)
}

## --- drop-in use in customKmeans() -------------------------------------------
## replace
##   gap_stat <- cluster::clusGap(x = data_df, FUNcluster = kmeans, K.max = max_k,
##                               B = random_set, d.power = d.power,
##                               nstart = nstart, iter.max = iter_max)
##   nc <- cluster::maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"],
##                        method = "Tibs2001SEmax", SE.factor = 1)
## with
##   Tab <- fastClusGap(data_df, K.max = max_k, B = random_set,
##                      nstart = nstart, iter.max = iter_max)
##   nc  <- selectK(Tab)
## and use `Tab` wherever `gap_stat$Tab` was used (same four columns). Note the
## gap-curve plot in customKmeans() builds its k axis as `1:max_k`; with k.min
## kept at 1 that is still correct, otherwise use attr(Tab, "k").
