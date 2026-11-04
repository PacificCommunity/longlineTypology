## R/stability.R -- cluster stability and separation diagnostics:
## centroid variability, Jaccard stability, temporal stability and
## silhouette. Compute functions return data frames (no printing), plot
## functions return ggplot objects. Cluster labels are the ordered ones
## (clusterLabels()) throughout.

utils::globalVariables(c("source", "panel", "jaccard", "share", "year",
						 "sil_avg", "sil_sd", "replicate"))


# ---- internal helpers ------------------------------------------------------

## Optimal one-to-one assignment (Hungarian algorithm) for a square cost
## matrix: returns perm with perm[i] = column assigned to row i, minimising
## sum(cost[cbind(i, perm[i])]). Same result as clue::solve_LSAP(), in base
## R. Vectors are 1-offset: element [j + 1] stands for index j (0 = dummy).
.hungarian <- function(cost) {
	cost <- as.matrix(cost)
	n <- nrow(cost)
	if (n != ncol(cost)) stop(".hungarian(): cost must be square")
	u <- v <- numeric(n + 1)
	p <- way <- integer(n + 1)
	for (i in seq_len(n)) {
		p[1] <- i
		j0   <- 0L
		minv <- rep(Inf, n + 1)
		used <- rep(FALSE, n + 1)
		repeat {
			used[j0 + 1] <- TRUE
			i0    <- p[j0 + 1]
			delta <- Inf
			j1    <- 0L
			for (j in seq_len(n)) {
				if (!used[j + 1]) {
					cur <- cost[i0, j] - u[i0 + 1] - v[j + 1]
					if (cur < minv[j + 1]) { minv[j + 1] <- cur; way[j + 1] <- j0 }
					if (minv[j + 1] < delta) { delta <- minv[j + 1]; j1 <- j }
				}
			}
			for (j in 0:n) {
				if (used[j + 1]) {
					u[p[j + 1] + 1] <- u[p[j + 1] + 1] + delta
					v[j + 1]        <- v[j + 1] - delta
				} else minv[j + 1] <- minv[j + 1] - delta
			}
			j0 <- j1
			if (p[j0 + 1] == 0) break
		}
		repeat {
			j1 <- way[j0 + 1]
			p[j0 + 1] <- p[j1 + 1]
			j0 <- j1
			if (j0 == 0) break
		}
	}
	perm <- integer(n)
	for (j in seq_len(n)) perm[p[j + 1]] <- j
	perm
}

## Raw kmeans label -> ordered label, for one K.
.labelMap <- function(diagnostics, k) {
	ro <- diagnostics$cluster_reorder[[as.character(k)]]
	if (is.null(ro)) return(function(x) x)
	function(x) as.integer(ro$new[match(x, ro$old)])
}

## Hennig-style stability: for each cluster of `ref`, the best Jaccard
## overlap with any cluster of `other` (which can have another K).
.maxJaccard <- function(ref, other) {
	tab   <- table(ref, other)
	n_ref <- rowSums(tab)
	n_oth <- colSums(tab)
	jac   <- tab / (outer(n_ref, n_oth, "+") - tab)
	stats::setNames(apply(jac, 1, max), rownames(tab))
}


# ---- centroid variability ---------------------------------------------------

#' Variability of the centroids across replicates
#'
#' Compares the centroids of the gap statistic replicates that selected K
#' (`diagnostics$kmeans_res`) with the centroids of the full fit
#' (`diagnostics$kmeans_full`). Each replicate's centroids are matched one
#' to one with the full-fit centroids by minimising the total distance
#' (Hungarian algorithm).
#'
#' @inheritParams fitRows
#' @param k Number of clusters.
#' @param scores PC scores of the fit-year rows (e.g.
#'   `pcaScores(diagnostics, data)[fitRows(diagnostics, data), ]`), used to
#'   express the deviations in PC standard deviations.
#'
#' @return A list with:
#' \describe{
#'   \item{`stats`}{one row per cluster and PC: `mean` and `sd` of the
#'     matched replicate centroids, `reference` (full-fit centroid), `rmse`
#'     (root-mean-square deviation of the replicate centroids from the
#'     reference, i.e. spread and bias) and `nrmse` (`rmse` divided by the
#'     sd of the PC scores). Sorted by decreasing `nrmse`. `NULL` when no
#'     replicate selected K.}
#'   \item{`points`}{the matched replicate centroids and the reference
#'     centroids (`source` = `"replicate"` or `"full fit"`), for
#'     [plotCentroidVariability()].}
#'   \item{`n_replicates`, `n_total`}{replicates that selected K, and all
#'     replicates.}
#' }
#'
#' @family stability diagnostics
#' @export
centroidVariability <- function(diagnostics, k, scores) {
	diagnostics <- .asDiagnostics(diagnostics)
	reference <- diagnostics$kmeans_full[[as.character(k)]]$centers
	if (is.null(reference)) stop("No kmeans fit for K = ", k, " in diagnostics")
	relabel <- .labelMap(diagnostics, k)
	n_total <- length(diagnostics$kmeans_res)
	reps    <- which(vapply(diagnostics$kmeans_res, function(x) nrow(x$centers),
							integer(1)) == nrow(reference))
	ref_pts <- data.frame(cluster = relabel(seq_len(nrow(reference))),
						  source = "full fit", replicate = NA_integer_, reference,
						  row.names = NULL)
	if (!length(reps))
		return(list(stats = NULL, points = ref_pts, n_replicates = 0L, n_total = n_total))

	pcs <- colnames(reference)
	matched <- do.call(rbind, lapply(reps, function(i) {
		cen  <- diagnostics$kmeans_res[[i]]$centers[, pcs, drop = FALSE]
		cost <- outer(seq_len(nrow(reference)), seq_len(nrow(cen)), Vectorize(function(a, b)
			sqrt(sum((reference[a, ] - cen[b, ])^2))))
		perm <- .hungarian(cost)
		data.frame(cluster = relabel(seq_len(nrow(reference))), source = "replicate",
				   replicate = i, cen[perm, , drop = FALSE], row.names = NULL)
	}))

	pc_sd <- apply(as.matrix(scores)[, pcs, drop = FALSE], 2, stats::sd)
	stats <- do.call(rbind, lapply(seq_len(nrow(reference)), function(r) {
		cl <- relabel(r)
		do.call(rbind, lapply(pcs, function(pc) {
			x    <- matched[[pc]][matched$cluster == cl]
			rmse <- sqrt(mean((x - reference[r, pc])^2))
			data.frame(cluster = cl, PC = pc, mean = mean(x), sd = stats::sd(x),
					   reference = unname(reference[r, pc]), rmse = rmse,
					   nrmse = rmse / pc_sd[[pc]])
		}))
	}))
	stats <- stats[order(-stats$nrmse), ]
	rownames(stats) <- NULL
	list(stats = stats, points = rbind(matched, ref_pts),
		 n_replicates = length(reps), n_total = n_total)
}


#' Plot the centroid variability
#'
#' Replicate centroids (crosses) and full-fit centroids (dots) of each
#' cluster, in the planes PC1 against each other PC.
#'
#' @param cv Output of [centroidVariability()].
#'
#' @return A ggplot object.
#'
#' @family stability diagnostics
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_shape_manual scale_alpha_manual labs theme_bw theme
#' @export
plotCentroidVariability <- function(cv) {
	pts <- cv$points
	pcs <- grep("^PC", names(pts), value = TRUE)
	if (length(pcs) < 2) stop("plotCentroidVariability(): needs at least 2 PCs")
	long <- do.call(rbind, lapply(pcs[-1], function(pc) data.frame(
		cluster = factor(pts$cluster, levels = sort(unique(pts$cluster))),
		source = factor(pts$source, levels = c("replicate", "full fit")),
		x = pts$PC1, y = pts[[pc]], panel = paste("PC1 vs", pc))))
	long$panel <- factor(long$panel, levels = paste("PC1 vs", pcs[-1]))
	long <- long[order(long$source), ]   # full-fit points drawn on top
	ggplot(long, aes(x, y, colour = cluster, shape = source, alpha = source)) +
		geom_point(size = 2.5) +
		facet_wrap(~ panel, scales = "free") +
		scale_shape_manual(values = c("replicate" = 4, "full fit" = 16), drop = FALSE) +
		scale_alpha_manual(values = c("replicate" = 0.4, "full fit" = 1), guide = "none") +
		labs(x = "PC1", y = NULL, colour = "Cluster", shape = NULL) +
		theme_bw() +
		theme(legend.position = "bottom")
}


# ---- Jaccard stability ------------------------------------------------------

#' Jaccard stability of the clusters
#'
#' Assigns the fit-year rows to the nearest centroid of each replicate fit
#' (`diagnostics$kmeans_res`, every replicate whatever K it selected). For
#' each cluster of the full fit, the stability in one replicate is the best
#' Jaccard overlap |A and B| / |A or B| with any replicate cluster (Hennig
#' 2007: below 0.6 the cluster dissolves, above 0.75 it is stable).
#'
#' @inheritParams centroidVariability
#' @param scores PC scores of the fit-year rows.
#' @param cluster Ordered cluster labels of the same rows (e.g.
#'   `clusterLabels(diagnostics, k)[fit_rows]`).
#'
#' @return A list with `summary` (one row per cluster: `jaccard`, the mean
#'   over replicates, and `jaccard_min`, the worst replicate) and `values`
#'   (one row per cluster and replicate).
#'
#' @family stability diagnostics
#' @export
jaccardStability <- function(diagnostics, scores, cluster) {
	diagnostics <- .asDiagnostics(diagnostics)
	if (length(cluster) != nrow(scores))
		stop("cluster has ", length(cluster), " values; scores has ", nrow(scores), " rows")
	if (length(unique(cluster)) < 2) stop("jaccardStability() needs at least 2 clusters")
	cl  <- sort(unique(cluster))
	jac <- vapply(diagnostics$kmeans_res, function(fit)
		.maxJaccard(cluster, .nearestCenter(scores, fit$centers))[as.character(cl)],
		numeric(length(cl)))
	jac <- matrix(jac, nrow = length(cl))
	values <- data.frame(cluster = rep(cl, ncol(jac)),
						 replicate = rep(seq_len(ncol(jac)), each = length(cl)),
						 jaccard = as.numeric(jac))
	summary <- data.frame(cluster = cl, jaccard = rowMeans(jac),
						  jaccard_min = apply(jac, 1, min))
	list(summary = summary, values = values)
}


#' Plot the Jaccard stability
#'
#' Distribution over replicates of each cluster's Jaccard stability, with
#' Hennig's thresholds (0.6 and 0.75).
#'
#' @param js Output of [jaccardStability()].
#'
#' @return A ggplot object.
#'
#' @family stability diagnostics
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline scale_y_continuous labs
#' @export
plotJaccard <- function(js) {
	ggplot(js$values, aes(factor(cluster), jaccard)) +
		geom_hline(yintercept = c(0.6, 0.75), linetype = "dashed", colour = "grey50") +
		geom_boxplot() +
		scale_y_continuous(limits = c(0, 1)) +
		labs(x = "Cluster", y = "Jaccard stability") +
		customTheme()
}


# ---- temporal stability -----------------------------------------------------

#' Temporal stability of the clusters
#'
#' Share of each year's rows in each cluster. Gradual trends are expected;
#' what this flags is a cluster that disappears and reappears, or jumps
#' from year to year.
#'
#' @param year Year of each row.
#' @param cluster Cluster label of each row.
#' @param threshold A cluster is present in a year when its share is at
#'   least `threshold`.
#' @param eps Added to the shares before taking logs, for years with share 0.
#'
#' @return A list with:
#' \describe{
#'   \item{`summary`}{one row per cluster: `mean_share`; `coverage`, the
#'     fraction of years where the cluster is present; `max_gap`, the
#'     longest run of years where it is absent between its first and last
#'     year present (a late start or an early end is not counted; a cluster
#'     never present gets all years); `jump`, the median of
#'     |log(share in year y+1 / share in year y)| (0.1 is about 10% change
#'     per year).}
#'   \item{`share`}{one row per year and cluster, for [plotTemporalShare()].
#'     Years without rows inside the range get share 0.}
#' }
#'
#' @family stability diagnostics
#' @export
temporalStability <- function(year, cluster, threshold = 0.01, eps = 0.001) {
	if (length(year) != length(cluster))
		stop("year and cluster must have the same length")
	year <- as.integer(year)
	cl   <- sort(unique(cluster))
	tab  <- table(factor(year, levels = seq(min(year), max(year))), factor(cluster, levels = cl))
	share <- sweep(tab, 1, pmax(rowSums(tab), 1), "/")

	longestRun <- function(x) {
		r <- rle(x)
		if (any(r$values)) max(r$lengths[r$values]) else 0L
	}
	internalGap <- function(present) {
		if (!any(present)) return(length(present))
		idx <- which(present)
		longestRun(!present[min(idx):max(idx)])
	}
	summary <- do.call(rbind, lapply(seq_along(cl), function(j) {
		s <- as.numeric(share[, j])
		present <- s >= threshold
		data.frame(cluster = cl[j], mean_share = mean(s), coverage = mean(present),
				   max_gap = internalGap(present),
				   jump = if (length(s) > 1) stats::median(abs(diff(log(s + eps)))) else NA_real_)
	}))
	long <- data.frame(year = as.integer(rownames(share)[row(share)]),
					   cluster = cl[col(share)], share = as.numeric(share))
	list(summary = summary, share = long, threshold = threshold)
}


#' Plot the cluster shares through time
#'
#' @param ts Output of [temporalStability()].
#'
#' @return A ggplot object (stacked areas of each cluster's share of the
#'   year's rows).
#'
#' @family stability diagnostics
#' @importFrom ggplot2 ggplot aes geom_area scale_y_continuous labs
#' @export
plotTemporalShare <- function(ts) {
	d <- ts$share
	d$cluster <- factor(d$cluster, levels = sort(unique(d$cluster)))
	ggplot(d, aes(year, share, fill = cluster)) +
		geom_area(colour = "white", linewidth = 0.2) +
		scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
						   expand = c(0, 0)) +
		labs(x = "Year", y = "Share of the year's rows", fill = "Cluster") +
		customTheme() +
		ggplot2::theme(legend.key.width = ggplot2::unit(0.6, "cm"))
}


# ---- silhouette ---------------------------------------------------------------

#' Silhouette of the clusters
#'
#' Silhouette widths ([cluster::silhouette()]) on samples of rows, since the
#' full distance matrix is too large for real datasets. Optionally compared
#' with a null reference: the same number of points drawn uniformly in the
#' bounding box of the scores (no cluster structure, as in the gap
#' statistic) and clustered by kmeans at the same K.
#'
#' @param scores PC scores (e.g. of the fit-year rows).
#' @param cluster Cluster label of each row of `scores`.
#' @param samples List of integer vectors, positions in `scores` of the rows
#'   of each sample (e.g. [subsampleReplicates()]). `NULL` uses all rows
#'   once (only for small data: the distance matrix has n^2 cells).
#' @param n_null Number of null datasets; `0` (default) skips the null
#'   reference.
#' @param null_size Points per null dataset. Default: size of the first
#'   sample.
#' @param seed Seed of the null datasets (the caller's RNG state is
#'   restored). `NULL` uses the current RNG state.
#'
#' @return A list with:
#' \describe{
#'   \item{`per_cluster`}{`cluster`, `sil_avg`, `sil_sd` (over all points
#'     of all samples) and `n_points`.}
#'   \item{`overall`}{mean silhouette width over all points.}
#'   \item{`null`}{`NULL`, or a list with `mean` and `sd` of the null mean
#'     silhouettes, `z` = (`overall` - `mean`) / `sd`, `n` and `size`.}
#'   \item{`n_samples`}{samples with at least 2 clusters.}
#' }
#'
#' @family stability diagnostics
#' @export
silhouetteStats <- function(scores, cluster, samples = NULL, n_null = 0,
							null_size = NULL, seed = NULL) {
	scores <- as.matrix(scores)
	if (length(cluster) != nrow(scores))
		stop("cluster has ", length(cluster), " values; scores has ", nrow(scores), " rows")
	if (is.null(samples)) samples <- list(seq_len(nrow(scores)))

	sil <- do.call(rbind, lapply(seq_along(samples), function(i) {
		r  <- samples[[i]]
		cl <- cluster[r]
		if (length(unique(cl)) < 2) return(NULL)
		s  <- cluster::silhouette(as.integer(factor(cl)), stats::dist(scores[r, , drop = FALSE]))
		data.frame(cluster = cl, sil_width = s[, "sil_width"], sample = i)
	}))
	if (is.null(sil)) stop("silhouetteStats(): every sample has fewer than 2 clusters")

	g <- factor(sil$cluster, levels = sort(unique(sil$cluster)))
	per_cluster <- data.frame(
		cluster  = as.integer(levels(g)),
		sil_avg  = vapply(split(sil$sil_width, g), mean, numeric(1)),
		sil_sd   = vapply(split(sil$sil_width, g), stats::sd, numeric(1)),
		n_points = as.integer(table(g)), row.names = NULL)
	overall <- mean(sil$sil_width)

	null <- NULL
	if (n_null > 0) {
		k <- length(unique(cluster))
		if (is.null(null_size)) null_size <- length(samples[[1]])
		rng <- apply(scores, 2, range)
		sil_null <- .withSeed(seed, vapply(seq_len(n_null), function(b) {
			X  <- apply(rng, 2, function(r) stats::runif(null_size, r[1], r[2]))
			km <- stats::kmeans(X, centers = k, nstart = 10, iter.max = 100)
			mean(cluster::silhouette(km$cluster, stats::dist(X))[, "sil_width"])
		}, numeric(1)))
		null <- list(mean = mean(sil_null), sd = stats::sd(sil_null),
					 z = (overall - mean(sil_null)) / stats::sd(sil_null),
					 n = n_null, size = null_size)
	}
	list(per_cluster = per_cluster, overall = overall, null = null,
		 n_samples = length(unique(sil$sample)))
}


#' Plot the silhouette per cluster
#'
#' Mean silhouette width of each cluster with +/- 1 sd; solid line: mean
#' over all points; dashed line: mean of the null reference, if computed.
#'
#' @param ss Output of [silhouetteStats()].
#'
#' @return A ggplot object.
#'
#' @family stability diagnostics
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline labs
#' @export
plotSilhouette <- function(ss) {
	p <- ggplot(ss$per_cluster, aes(factor(cluster), sil_avg)) +
		geom_hline(yintercept = ss$overall, colour = "grey40") +
		geom_errorbar(aes(ymin = sil_avg - sil_sd, ymax = sil_avg + sil_sd), width = 0.2) +
		geom_point(size = 2.5) +
		labs(x = "Cluster", y = "Silhouette width") +
		customTheme()
	if (!is.null(ss$null))
		p <- p + geom_hline(yintercept = ss$null$mean, linetype = "dashed", colour = "firebrick")
	p
}
