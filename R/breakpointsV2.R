## R/breakpointsV2.R -- breakpoint analysis, version 2.
##
## Replaces the detection / validation logic of R/breakpoints.R (kept for
## backward compatibility and comparison). Design, calibration and limits are
## documented in breakpoint_report.qmd.
##
## Per (cluster, variable) monthly series:
##   1. candidates: EnvCpt (default; mean/trend x iid/AR1 x with/without
##      changepoints, best by AIC) or PELT on the noise-scaled series.
##   2. characterisation of each candidate on +/- `window` months:
##        ramp  y = a + b t + delta * ramp(t; start, duration)
##        pulse y = a + b t + delta * exp(-(t - start) / tau), t >= start
##      against no change y = a + b t, over grids of start / duration / tau.
##      Likelihood ratio with an AR(1)-corrected effective sample size;
##      profile interval [d_lo, d_hi] on the transition duration.
##   3. classification against a timescale T_c: cohort timescale (von
##      Bertalanffy age at the local mean length) or a fixed number of months
##      (HBF: 12 by default in analyseClusterBreakpoints()):
##        pulse fits better and tau <= transient_tau_frac * T_c -> "transient"
##        d_lo >  T_c   -> "gradual"
##        d_hi <= T_c   -> "abrupt"; d_lo <= T_c < d_hi -> "ambiguous"
##      then, for abrupt AND ambiguous changes (if spatial data given):
##        flag mix changed more than usual     -> "composition"
##        else same flags moved more than usual -> "relocation"
##        else abrupt -> "unexplained" (ambiguous stays "ambiguous")
##   4. analyseClusterBreakpoints() (needs all clusters):
##      - "composition" changes whose flag gains / losses are mirrored by
##        another cluster at the same time -> "switch" (rows moving between
##        clusters: the clustering-relevant signal). Not mirrored = the
##        cluster's membership changed for another reason (fleet entering /
##        leaving the fishery): stays "composition".
##      - changes found at the same time in most clusters -> common = TRUE
##        (data / reporting events). Composition / switch changes are never
##        common and do not count towards it.

utils::globalVariables(c("x", "start_date", "end_plot", "lab"))

.BP_ORIGIN <- as.Date("1900-01-15")

.bpMonthIdx <- function(d, origin = .BP_ORIGIN) {
	d <- as.POSIXlt(d); o <- as.POSIXlt(origin)
	(d$year - o$year) * 12 + (d$mon - o$mon) + 1
}

.bpIdxToDate <- function(i, origin = .BP_ORIGIN) {
	out <- rep(as.Date(NA), length(i)); ok <- !is.na(i)
	if (any(ok)) out[ok] <- seq(origin, by = "month", length.out = max(i[ok]))[i[ok]]
	out
}

.bpRamp <- function(t, s, d) pmin(pmax((t - s + 1) / d, 0), 1)

.bpHaversine <- function(lon1, lat1, lon2, lat2) {
	rad <- pi / 180
	dlat <- (lat2 - lat1) * rad; dlon <- (lon2 - lon1) * rad
	a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
	6371 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

#' Age (months) at length from the von Bertalanffy growth curve
#'
#' Unlike [ageFromLength()], t0 is treated in years throughout.
#'
#' @param len Numeric, fork length (cm).
#' @param Linf,K,t0 Growth parameters (cm, per year, years).
#' @return Numeric, age in months; `NA` where `len >= Linf`.
#' @family breakpoint analysis v2
#' @export
vbAgeMonths <- function(len, Linf = 150.3, K = 0.442, t0 = -0.244) {
	out <- rep(NA_real_, length(len)); ok <- !is.na(len) & len < Linf
	out[ok] <- (t0 - log(1 - len[ok] / Linf) / K) * 12
	out
}

.bpCandidatesPELT <- function(t, y, smooth_span = 0.05, minseglen = 12, penalty = "BIC") {
	sigma <- stats::mad(diff(y)) / sqrt(2)
	z <- if (is.finite(sigma) && sigma > 0) y / sigma else y
	if (!is.null(smooth_span) && !is.na(smooth_span))
		z <- as.numeric(stats::predict(stats::loess(z ~ t, span = smooth_span, degree = 1)))
	if (length(z) <= 2 * minseglen) return(integer(0))
	cp <- changepoint::cpts(changepoint::cpt.mean(z, method = "PELT", penalty = penalty,
												  minseglen = minseglen, test.stat = "Normal"))
	cp[cp < length(z)]
}

.bpCandidatesEnvCpt <- function(y, models = c("mean", "meancpt", "meanar1", "meanar1cpt",
											  "trend", "trendcpt", "trendar1", "trendar1cpt"),
								minseglen = 12) {
	fit <- EnvCpt::envcpt(y, models = models, minseglen = minseglen, verbose = FALSE)
	aic <- stats::AIC(fit)
	if (is.null(names(aic))) names(aic) <- models
	best <- names(aic)[which.min(aic)]
	if (!grepl("cpt", best)) return(integer(0))
	cp <- changepoint::cpts(fit[[best]])
	cp[cp < length(y)]
}

.bpCharacterise <- function(t, y, cand_t, T_c, window = 72,
							durations = c(1, 3, 6, 12, 18, 24, 36, 48, 72, 96),
							pulse_taus = c(3, 6, 9, 12, 18, 24, 36, 48, 72, 120),
							start_slack = 12, margin = 6) {
	in_w <- t >= cand_t - window & t <= cand_t + window
	tw <- t[in_w]; yw <- y[in_w]; n <- length(yw)
	if (n < 24) return(NULL)
	X0   <- cbind(1, tw - cand_t)
	rss0 <- sum(stats::.lm.fit(X0, yw)$residuals^2)

	prof <- stats::setNames(rep(NA_real_, length(durations)), durations)
	best <- list(rss = Inf)
	for (d in durations) {
		starts <- seq(cand_t - d - start_slack, cand_t + start_slack)
		starts <- starts[starts > min(tw) + margin & starts + d < max(tw) - margin]
		if (!length(starts)) next
		rss_d <- vapply(starts, function(s)
			sum(stats::.lm.fit(cbind(X0, .bpRamp(tw, s, d)), yw)$residuals^2), numeric(1))
		k <- which.min(rss_d)
		prof[as.character(d)] <- rss_d[k]
		if (rss_d[k] < best$rss) best <- list(rss = rss_d[k], s = starts[k], d = d)
	}
	if (!is.finite(best$rss)) return(NULL)

	fit   <- stats::.lm.fit(cbind(X0, .bpRamp(tw, best$s, best$d)), yw)
	delta <- fit$coefficients[3]
	res   <- fit$residuals
	phi   <- suppressWarnings(stats::cor(res[-1], res[-n]))
	phi   <- if (is.finite(phi)) min(max(phi, 0), 0.95) else 0
	n_eff <- n * (1 - phi) / (1 + phi)

	ok    <- is.finite(prof)
	in_ci <- n_eff * log(prof[ok] / best$rss) <= stats::qchisq(0.95, 1)
	d_ok  <- durations[ok]

	p_starts <- seq(cand_t - 24, cand_t + start_slack)
	p_starts <- p_starts[p_starts > min(tw) + margin & p_starts < max(tw) - margin]
	best_p <- list(rss = Inf, s = NA_real_, tau = NA_real_)
	for (tau in pulse_taus) {
		if (!length(p_starts)) break
		rss_p <- vapply(p_starts, function(s)
			sum(stats::.lm.fit(cbind(X0, ifelse(tw >= s, exp(-(tw - s) / tau), 0)), yw)$residuals^2),
			numeric(1))
		k <- which.min(rss_p)
		if (rss_p[k] < best_p$rss) best_p <- list(rss = rss_p[k], s = p_starts[k], tau = tau)
	}
	rss_best <- min(best$rss, best_p$rss)
	lr <- n_eff * log(rss0 / rss_best)

	pre <- y[t >= best$s - 12 & t < best$s]
	data.frame(cand_t = cand_t, s_hat = best$s, d_hat = best$d,
			   d_lo = min(d_ok[in_ci]), d_hi = max(d_ok[in_ci]),
			   T_c = T_c, delta = unname(delta), pre_level = if (length(pre)) mean(pre) else NA_real_,
			   lr = lr, p_nominal = stats::pchisq(lr, df = 3, lower.tail = FALSE),
			   phi = phi, n_eff = n_eff,
			   pulse_tau = best_p$tau,
			   pulse_gain = if (is.finite(best_p$rss)) n_eff * log(best$rss / best_p$rss) else NA_real_)
}

.bpDedupe <- function(ch, gap = 6) {
	ch <- ch[order(-ch$lr), , drop = FALSE]
	keep <- logical(nrow(ch))
	for (i in seq_len(nrow(ch))) {
		a1 <- ch$s_hat[i] - gap; a2 <- ch$s_hat[i] + ch$d_hat[i] + gap
		clash <- keep & (ch$s_hat - gap <= a2) & (ch$s_hat + ch$d_hat + gap >= a1)
		keep[i] <- !any(clash)
	}
	ch[keep, , drop = FALSE]
}

#' Per-month, per-group position sums for one cluster
#'
#' Pre-computation for the composition / relocation test of
#' [analyseBreakpointsV2()]; build once per cluster.
#'
#' @param df Row-level data with `cluster`, `longitude`, `latitude`, the date
#'   column and the group column.
#' @param cluster_id Cluster to tabulate.
#' @param group Column defining the fleet units (default `"flag"`).
#' @param date_col Date column (default `"date"`).
#' @return A list (`mi`, `n`, `slon`, `slat`).
#' @family breakpoint analysis v2
#' @export
spatialMonthTable <- function(df, cluster_id, group = "flag", date_col = "date") {
	d <- df[df$cluster == cluster_id & !is.na(df$longitude) & !is.na(df$latitude), , drop = FALSE]
	mi <- .bpMonthIdx(d[[date_col]])
	g  <- as.character(d[[group]])
	n    <- tapply(rep(1, nrow(d)), list(mi, g), sum)
	slon <- tapply(d$longitude, list(mi, g), sum)
	slat <- tapply(d$latitude,  list(mi, g), sum)
	n[is.na(n)] <- 0; slon[is.na(slon)] <- 0; slat[is.na(slat)] <- 0
	list(mi = as.integer(rownames(n)), n = n, slon = slon, slat = slat)
}

.bpWindowStats <- function(tab, a1, a2, b1, b2, min_rows = 5) {
	ia <- tab$mi >= a1 & tab$mi < a2
	ib <- tab$mi >= b1 & tab$mi < b2
	na <- colSums(tab$n[ia, , drop = FALSE]); nb <- colSums(tab$n[ib, , drop = FALSE])
	if (sum(na) == 0 || sum(nb) == 0) return(c(comp = NA_real_, reloc = NA_real_))
	comp <- 0.5 * sum(abs(na / sum(na) - nb / sum(nb)))
	both <- na >= min_rows & nb >= min_rows
	reloc <- NA_real_
	if (any(both)) {
		lon_a <- colSums(tab$slon[ia, both, drop = FALSE]) / na[both]
		lat_a <- colSums(tab$slat[ia, both, drop = FALSE]) / na[both]
		lon_b <- colSums(tab$slon[ib, both, drop = FALSE]) / nb[both]
		lat_b <- colSums(tab$slat[ib, both, drop = FALSE]) / nb[both]
		w <- pmin(na[both], nb[both])
		reloc <- sum(w * .bpHaversine(lon_a, lat_a, lon_b, lat_b)) / sum(w)
	}
	c(comp = comp, reloc = reloc)
}

.bpSpatialTest <- function(s, d, tab, window = 12, null_step = 3, min_rows = 5) {
	obs <- .bpWindowStats(tab, s - window, s, s + d, s + d + window, min_rows)
	lo <- min(tab$mi) + window; hi <- max(tab$mi) - d - window
	p_comp <- p_reloc <- NA_real_
	if (lo <= hi) {
		u <- seq(lo, hi, by = null_step)
		nul <- vapply(u, function(ui) .bpWindowStats(tab, ui - window, ui, ui + d, ui + d + window, min_rows),
					  numeric(2))
		empP <- function(o, v) { v <- v[!is.na(v)]; if (is.na(o) || !length(v)) NA_real_ else (1 + sum(v >= o)) / (1 + length(v)) }
		p_comp  <- empP(obs[["comp"]],  nul[1, ])
		p_reloc <- empP(obs[["reloc"]], nul[2, ])
	}
	list(comp_change = obs[["comp"]], p_comp = p_comp, reloc_km = obs[["reloc"]], p_reloc = p_reloc)
}

.bpEmpty <- function() {
	data.frame(start_date = as.Date(character()), end_date = as.Date(character()),
			   d_hat = numeric(), d_lo = numeric(), d_hi = numeric(), T_c = numeric(),
			   delta = numeric(), delta_rel = numeric(), lr = numeric(), p_nominal = numeric(),
			   phi = numeric(), pulse_tau = numeric(), pulse_gain = numeric(), class = character(),
			   comp_change = numeric(), p_comp = numeric(), reloc_km = numeric(), p_reloc = numeric(),
			   label = character(), stringsAsFactors = FALSE)
}

#' Detect and classify changes in one cluster time series
#'
#' See the file header of R/breakpointsV2.R and breakpoint_report.qmd for the
#' method, its calibration and detection limits.
#'
#' @param dates Date vector (monthly).
#' @param x Numeric series (`CPUE` is analysed on the log scale).
#' @param variable `"CPUE"`, `"mean_len"` or `"mean_hbf"`.
#' @param mean_len Mean length aligned to `dates`; sets the cohort timescale.
#' @param detector `"envcpt"` (default) or `"pelt"`.
#' @param detector_args List of extra arguments for the candidate generator.
#' @param deseason Remove the monthly climatology first. Default `TRUE`.
#' @param window Half-width (months) of the local fit. Default 72.
#' @param alpha Significance level of the change test. Default 0.1 (EnvCpt's
#'   own model choice is the main false-positive filter; see report).
#' @param abrupt_months `NULL` = cohort timescale from `mean_len`; a number =
#'   fixed threshold (months).
#' @param transient_tau_frac A pulse decaying with time constant <= this x
#'   T_c is labelled "transient". Default 0.5.
#' @param spatial_table [spatialMonthTable()] for this cluster, or `NULL` to
#'   skip the composition / relocation test.
#' @param spatial_window,spatial_alpha Window (months) and significance level
#'   of that test.
#' @param Linf,K,t0 Growth parameters for [vbAgeMonths()].
#' @return Data frame, one row per significant change (0 rows if none).
#' @family breakpoint analysis v2
#' @export
analyseBreakpointsV2 <- function(dates, x, variable, mean_len,
								 detector = c("envcpt", "pelt"), detector_args = list(),
								 deseason = TRUE, window = 72, alpha = 0.1,
								 abrupt_months = NULL, transient_tau_frac = 0.5,
								 spatial_table = NULL, spatial_window = 12, spatial_alpha = 0.05,
								 Linf = 150.3, K = 0.442, t0 = -0.244) {
	detector <- match.arg(detector)
	empty <- .bpEmpty()
	ok <- !is.na(x) & is.finite(x)
	dates <- dates[ok]; x <- x[ok]; ml <- mean_len[ok]
	ord <- order(dates); dates <- dates[ord]; x <- x[ord]; ml <- ml[ord]
	if (length(x) < 48) return(empty)

	t <- .bpMonthIdx(dates)
	y <- x
	if (variable == "CPUE") {
		pos <- x[x > 0]
		y <- log(pmax(x, if (length(pos)) min(pos) / 2 else 1e-6))
	}
	if (deseason) {
		moy <- (t - 1) %% 12
		clim <- tapply(y, moy, mean)
		y <- as.numeric(y - clim[as.character(moy)] + mean(y))
	}

	cand <- tryCatch(switch(detector,
		envcpt = do.call(.bpCandidatesEnvCpt, c(list(y = y), detector_args)),
		pelt   = do.call(.bpCandidatesPELT,   c(list(t = t, y = y), detector_args))),
		error = function(e) structure(integer(0), error = conditionMessage(e)))
	if (!length(cand)) { attr(empty, "error") <- attr(cand, "error"); return(empty) }

	rows <- lapply(unique(t[cand]), function(ct) {
		T_c <- if (is.null(abrupt_months))
			suppressWarnings(stats::median(vbAgeMonths(ml[abs(t - ct) <= 12], Linf, K, t0), na.rm = TRUE))
		else abrupt_months
		.bpCharacterise(t, y, ct, T_c = T_c, window = window)
	})
	rows <- rows[!vapply(rows, is.null, logical(1))]
	if (!length(rows)) return(empty)
	ch <- do.call(rbind, rows)
	ch <- ch[ch$p_nominal < alpha, , drop = FALSE]
	if (nrow(ch) == 0) return(empty)
	ch <- .bpDedupe(ch)

	ch$start_date <- .bpIdxToDate(ch$s_hat)
	ch$end_date   <- .bpIdxToDate(ch$s_hat + ch$d_hat)
	ch$delta_rel  <- if (variable == "CPUE") expm1(ch$delta) else ch$delta / abs(ch$pre_level)
	transient <- !is.na(ch$pulse_gain) & ch$pulse_gain > 0 & !is.na(ch$T_c) &
				 ch$pulse_tau <= transient_tau_frac * ch$T_c
	ch$class <- ifelse(is.na(ch$T_c), "no_timescale",
				ifelse(transient,         "transient",
				ifelse(ch$d_hi <= ch$T_c, "abrupt",
				ifelse(ch$d_lo >  ch$T_c, "gradual", "ambiguous"))))

	ch$comp_change <- ch$p_comp <- ch$reloc_km <- ch$p_reloc <- NA_real_
	ch$label <- ch$class
	if (!is.null(spatial_table)) {
		for (i in seq_len(nrow(ch))) {
			sp <- .bpSpatialTest(ch$s_hat[i], ch$d_hat[i], spatial_table, window = spatial_window)
			ch$comp_change[i] <- sp$comp_change; ch$p_comp[i]  <- sp$p_comp
			ch$reloc_km[i]    <- sp$reloc_km;    ch$p_reloc[i] <- sp$p_reloc
		}
		ab   <- ch$class == "abrupt"
		fast <- ab | ch$class == "ambiguous"   # abrupt, or possibly abrupt
		comp_sig  <- !is.na(ch$p_comp)  & ch$p_comp  < spatial_alpha
		reloc_sig <- !is.na(ch$p_reloc) & ch$p_reloc < spatial_alpha
		ch$label[fast & comp_sig]               <- "composition"
		ch$label[fast & !comp_sig & reloc_sig]  <- "relocation"
		ch$label[ab   & !comp_sig & !reloc_sig] <- "unexplained"
	}
	out <- ch[order(ch$s_hat), names(empty)]
	rownames(out) <- NULL
	out
}

#' Mark changes found at the same time in most clusters
#'
#' A change is "common" if, for the same variable, at least
#' `ceiling(min_share * n_clusters)` clusters (and at least 2) have a change
#' starting within `tol_months` of it. Such changes point to data or
#' reporting events rather than to the clustering. Changes labelled
#' "composition" or "switch" (the cluster's membership changed) are never
#' common and are not counted: a reporting event does not change which flags
#' are in a cluster, and rows moving between clusters change several
#' clusters at once, which would otherwise pass as common (always, at K = 2).
#'
#' @param bp Output of [analyseClusterBreakpoints()] (columns `cluster`,
#'   `variable`, `start_date`).
#' @param n_clusters Number of clusters analysed.
#' @param tol_months Tolerance on start dates. Default 6.
#' @param min_share Share of clusters required. Default 0.75.
#' @return `bp` with a logical column `common`.
#' @family breakpoint analysis v2
#' @export
flagCommonBreaks <- function(bp, n_clusters, tol_months = 6, min_share = 0.75) {
	bp$common <- logical(nrow(bp))
	need <- max(2, ceiling(min_share * n_clusters))
	if (nrow(bp) == 0 || n_clusters < 2) return(bp)
	mi <- .bpMonthIdx(bp$start_date)
	eligible <- if (is.null(bp$label)) rep(TRUE, nrow(bp)) else !(bp$label %in% c("composition", "switch"))
	for (i in which(eligible)) {
		same <- eligible & bp$variable == bp$variable[i] & abs(mi - mi[i]) <= tol_months
		bp$common[i] <- length(unique(bp$cluster[same])) >= need
	}
	bp
}

## Rows per flag after minus before the change, in one cluster.
.bpFlagDelta <- function(tab, s, d, window) {
	ia <- tab$mi >= s - window & tab$mi < s
	ib <- tab$mi >= s + d & tab$mi < s + d + window
	colSums(tab$n[ib, , drop = FALSE]) - colSums(tab$n[ia, , drop = FALSE])
}

## Is a composition change in cluster `cl` (start month index s, duration d)
## mirrored by another cluster? Per flag, rows gained by `cl` and lost by the
## other cluster (and the reverse) count as moved (the smaller of the two).
## Mirror share = moved / total per-flag change in `cl` (0-1). Overall fleet
## growth or decline (same sign in both clusters) does not count.
## Returns list(share, cluster) for the best-matching other cluster.
.bpMirror <- function(tabs, cl, s, d, window = 12) {
	dc <- .bpFlagDelta(tabs[[cl]], s, d, window)
	tot <- sum(abs(dc))
	others <- setdiff(names(tabs), cl)
	if (!length(others) || tot == 0) return(list(share = NA_real_, cluster = NA))
	shares <- vapply(others, function(o) {
		d2 <- .bpFlagDelta(tabs[[o]], s, d, window)
		fl <- union(names(dc), names(d2))
		a <- stats::setNames(numeric(length(fl)), fl); b <- a
		a[names(dc)] <- dc; b[names(d2)] <- d2
		sum(pmin(pmax(a, 0), pmax(-b, 0)) + pmin(pmax(-a, 0), pmax(b, 0))) / tot
	}, numeric(1))
	list(share = max(shares), cluster = others[which.max(shares)])
}

#' Breakpoint analysis for all clusters and variables
#'
#' @param ts Output of [aggregateClusterTimeSeries()] (`cluster`, `date`,
#'   `CPUE`, `mean_len`, `mean_hbf`).
#' @param df Row-level data (for the composition / relocation test), or
#'   `NULL` to skip it.
#' @param variables Variables to analyse.
#' @param group Fleet-unit column of `df`. Default `"flag"`.
#' @param hbf_abrupt_months Fixed timescale for `mean_hbf` (gear changes are
#'   not tied to cohort turnover). Default 12. `NULL` = cohort timescale.
#' @param common_tol_months,common_min_share Passed to [flagCommonBreaks()].
#' @param mirror_min Share of a composition change mirrored by another
#'   cluster (rows gained here lost there, or the reverse) above which it is
#'   relabelled "switch". Default 0.5.
#' @param spatial_window Months before / after for the composition,
#'   relocation and mirror tests. Default 12.
#' @param ... Passed to [analyseBreakpointsV2()].
#' @return Data frame with `cluster`, `variable`, the columns of
#'   [analyseBreakpointsV2()], `mirror_share`, `mirror_cluster` and `common`.
#' @family breakpoint analysis v2
#' @export
analyseClusterBreakpoints <- function(ts, df = NULL, variables = c("CPUE", "mean_len", "mean_hbf"),
									  group = "flag", hbf_abrupt_months = 12,
									  common_tol_months = 6, common_min_share = 0.75,
									  mirror_min = 0.5, spatial_window = 12, ...) {
	clusters <- sort(unique(ts$cluster))
	tabs <- if (is.null(df)) NULL else
		stats::setNames(lapply(clusters, function(cl) spatialMonthTable(df, cl, group = group)),
						as.character(clusters))
	rows <- list()
	for (cl in clusters) {
		d <- ts[ts$cluster == cl, , drop = FALSE]
		tab <- if (is.null(tabs)) NULL else tabs[[as.character(cl)]]
		for (v in variables) {
			res <- analyseBreakpointsV2(d$date, d[[v]], v, d$mean_len, spatial_table = tab,
										spatial_window = spatial_window,
										abrupt_months = if (v == "mean_hbf") hbf_abrupt_months else NULL, ...)
			if (nrow(res)) rows[[length(rows) + 1]] <- cbind(cluster = cl, variable = v, res,
															 stringsAsFactors = FALSE)
		}
	}
	out <- if (length(rows)) do.call(rbind, rows) else
		cbind(data.frame(cluster = ts$cluster[0], variable = character()), .bpEmpty())

	out$mirror_share   <- rep(NA_real_, nrow(out))
	out$mirror_cluster <- rep(NA_character_, nrow(out))
	if (!is.null(tabs)) for (i in which(out$label == "composition")) {
		m <- .bpMirror(tabs, as.character(out$cluster[i]), .bpMonthIdx(out$start_date[i]), out$d_hat[i],
					   window = spatial_window)
		out$mirror_share[i]   <- m$share
		out$mirror_cluster[i] <- as.character(m$cluster)
		if (!is.na(m$share) && m$share >= mirror_min) out$label[i] <- "switch"
	}
	flagCommonBreaks(out, length(clusters), common_tol_months, common_min_share)
}

#' Plot one series with the changes found by analyseBreakpointsV2()
#'
#' @param dates,x The series.
#' @param bp Rows of [analyseBreakpointsV2()] / [analyseClusterBreakpoints()]
#'   for this series (may be 0 rows).
#' @param y_lab Y-axis label.
#' @return A ggplot object.
#' @family breakpoint analysis v2
#' @export
plotBreakpointsV2 <- function(dates, x, bp, y_lab = NULL) {
	cols <- c(unexplained = "#d7301f", abrupt = "#d7301f", switch = "#ae017e", composition = "#7a0177",
			  relocation = "#2b8cbe", transient = "#41ab5d", gradual = "#969696",
			  ambiguous = "#fdae61", common = "#000000")
	p <- ggplot2::ggplot(data.frame(date = dates, x = x), ggplot2::aes(date, x)) +
		ggplot2::geom_line(linewidth = 0.3)
	if (!is.null(bp) && nrow(bp) > 0) {
		bp$lab <- if ("common" %in% names(bp)) ifelse(bp$common, "common", bp$label) else bp$label
		bp$end_plot <- pmax(bp$end_date, bp$start_date + 20)
		p <- p +
			ggplot2::geom_rect(data = bp, ggplot2::aes(xmin = start_date, xmax = end_plot, ymin = -Inf, ymax = Inf,
													   fill = lab), alpha = 0.35, inherit.aes = FALSE) +
			ggplot2::geom_vline(data = bp, ggplot2::aes(xintercept = start_date, colour = lab), linewidth = 0.4) +
			ggplot2::scale_fill_manual(values = cols, name = NULL) +
			ggplot2::scale_colour_manual(values = cols, guide = "none")
	}
	p + ggplot2::labs(x = NULL, y = y_lab) + ggplot2::theme_bw()
}
