## R/breakpointInspect.R -- zoom on one change of the breakpoint analysis v2.
##
## For a change found by analyseClusterBreakpoints(), or for any date / date
## range given by hand, recompute every quantity the classification uses and
## the thresholds it is compared with:
##   - local fit: ramp / pulse / no change, LR test, duration profile and
##     [d_lo, d_hi], T_c and the T_c values that would change the class;
##   - spatial test: comp_change and reloc_km with their null distributions,
##     the value each must exceed to be significant, sensitivity to the
##     transition duration and to the spatial window;
##   - per-flag shares and centroids before / after, mirror shares;
##   - context: other variables of the cluster, same variable in the other
##     clusters (common rule), EnvCpt model comparison.
## Uses the same internals as analyseBreakpointsV2() (.bpPrepare,
## .bpTimescale, .bpCharacterise, .bpClassify, .bpSpatialTest,
## .bpSpatialLabel, .bpMirror), so a detected change is reproduced exactly.

utils::globalVariables(c("date", "value", "fit", "d", "stat", "in_ci", "share", "when", "flag",
						 "lon_before", "lat_before", "lon_after", "lat_after", "stat_name", "obs",
						 "thr", "xmin", "xmax", "series", "y"))

.BP_LABEL_PRIORITY <- c(unexplained = 1, switch = 2, ambiguous = 3, composition = 4, relocation = 5,
						transient = 6, gradual = 7, no_timescale = 9)

#' Set up the breakpoint inspection for one clustering
#'
#' Runs [analyseClusterBreakpoints()] once (unless `bp` is given) and keeps
#' everything [inspectBreakpoint()] needs: the cluster series, the per-cluster
#' spatial tables, the changes found and the settings used. Settings have
#' the same names and defaults as in [analyseClusterBreakpoints()] and
#' [analyseBreakpointsV2()].
#'
#' @param ts Output of [aggregateClusterTimeSeries()].
#' @param df Row-level data (`cluster`, `date`, `longitude`, `latitude` and
#'   the `group` column).
#' @param bp Optional output of [analyseClusterBreakpoints()] obtained with
#'   the same settings; recomputed if `NULL`. It must contain `cand_date`.
#' @param variables,group,hbf_abrupt_months,common_tol_months,common_min_share,mirror_min,spatial_window
#'   As in [analyseClusterBreakpoints()].
#' @param detector,detector_args,deseason,window,alpha,transient_tau_frac,spatial_alpha,Linf,K,t0
#'   As in [analyseBreakpointsV2()].
#' @return A list of class `bpContext`.
#' @family breakpoint analysis v2
#' @export
breakpointContext <- function(ts, df, bp = NULL, variables = c("CPUE", "mean_len", "mean_hbf"),
							  group = "flag", hbf_abrupt_months = 12,
							  common_tol_months = 6, common_min_share = 0.75, mirror_min = 0.5,
							  spatial_window = 12, detector = "envcpt", detector_args = list(),
							  deseason = TRUE, window = 72, alpha = 0.1, transient_tau_frac = 0.5,
							  spatial_alpha = 0.05, Linf = 150.3, K = 0.442, t0 = -0.244) {
	settings <- list(variables = variables, group = group, hbf_abrupt_months = hbf_abrupt_months,
					 common_tol_months = common_tol_months, common_min_share = common_min_share,
					 mirror_min = mirror_min, spatial_window = spatial_window, detector = detector,
					 detector_args = detector_args, deseason = deseason, window = window, alpha = alpha,
					 transient_tau_frac = transient_tau_frac, spatial_alpha = spatial_alpha,
					 Linf = Linf, K = K, t0 = t0)
	if (is.null(bp))
		bp <- analyseClusterBreakpoints(ts, df, variables = variables, group = group,
										hbf_abrupt_months = hbf_abrupt_months,
										common_tol_months = common_tol_months,
										common_min_share = common_min_share, mirror_min = mirror_min,
										spatial_window = spatial_window, detector = detector,
										detector_args = detector_args, deseason = deseason, window = window,
										alpha = alpha, transient_tau_frac = transient_tau_frac,
										spatial_alpha = spatial_alpha, Linf = Linf, K = K, t0 = t0)
	if (!"cand_date" %in% names(bp))
		stop("`bp` has no cand_date column: it was made with an older breakpointsV2.R. Leave bp = NULL.")
	clusters <- sort(unique(ts$cluster))
	tabs <- stats::setNames(lapply(clusters, function(cl) spatialMonthTable(df, cl, group = group)),
							as.character(clusters))
	structure(list(ts = ts, tabs = tabs, bp = bp, clusters = clusters, settings = settings),
			  class = "bpContext")
}

#' Order changes from most to least concerning
#'
#' Order: `unexplained`, `switch`, `ambiguous`, `composition`, `relocation`,
#' `transient`, `gradual`; changes flagged `common` (data or fleet-wide
#' events) come after all of these. Ties by cluster, variable and date.
#'
#' @param bp Output of [analyseClusterBreakpoints()].
#' @return `bp` reordered, with a `priority` column (1 = most concerning).
#' @family breakpoint analysis v2
#' @export
rankBreakpoints <- function(bp) {
	if (!nrow(bp)) { bp$priority <- numeric(0); return(bp) }
	pr <- unname(.BP_LABEL_PRIORITY[bp$label])
	pr[is.na(pr)] <- 9
	if ("common" %in% names(bp)) pr[bp$common %in% TRUE] <- 8
	bp$priority <- pr
	out <- bp[order(pr, bp$cluster, match(bp$variable, c("CPUE", "mean_len", "mean_hbf")), bp$start_date), ,
			  drop = FALSE]
	rownames(out) <- NULL
	out
}

## Smallest value an observed statistic must EXCEED for the empirical p of
## .bpSpatialTest() to fall below alpha. -Inf: any value passes; NA: no value
## can pass (too few null positions, see breakpoint_report.qmd).
.bpEmpThreshold <- function(v, alpha) {
	v <- sort(v[!is.na(v)], decreasing = TRUE); N <- length(v)
	if (!N) return(NA_real_)
	c_max <- ceiling(alpha * (1 + N) - 1) - 1
	if (c_max < 0) return(NA_real_)
	if (c_max >= N) return(-Inf)
	v[c_max + 1]
}

## Per-flag rows, shares and centroids in the windows before / after a change.
.bpFlagTable <- function(tab, s, d, window = 12, min_rows = 5) {
	ia <- tab$mi >= s - window & tab$mi < s
	ib <- tab$mi >= s + d & tab$mi < s + d + window
	na <- colSums(tab$n[ia, , drop = FALSE]); nb <- colSums(tab$n[ib, , drop = FALSE])
	lon_a <- colSums(tab$slon[ia, , drop = FALSE]) / na; lat_a <- colSums(tab$slat[ia, , drop = FALSE]) / na
	lon_b <- colSums(tab$slon[ib, , drop = FALSE]) / nb; lat_b <- colSums(tab$slat[ib, , drop = FALSE]) / nb
	used <- na >= min_rows & nb >= min_rows
	out <- data.frame(flag = colnames(tab$n), n_before = na, n_after = nb,
					  share_before = if (sum(na)) na / sum(na) else NA_real_,
					  share_after  = if (sum(nb)) nb / sum(nb) else NA_real_,
					  lon_before = ifelse(na > 0, lon_a, NA), lat_before = ifelse(na > 0, lat_a, NA),
					  lon_after  = ifelse(nb > 0, lon_b, NA), lat_after  = ifelse(nb > 0, lat_b, NA),
					  stringsAsFactors = FALSE)
	out$d_share   <- out$share_after - out$share_before
	out$comp_part <- 0.5 * abs(out$d_share)   # sums to comp_change
	out$dist_km   <- ifelse(na > 0 & nb > 0, .bpHaversine(out$lon_before, out$lat_before,
														   out$lon_after, out$lat_after), NA)
	out$used_reloc  <- used
	out$reloc_weight <- ifelse(used, pmin(na, nb) / sum(pmin(na, nb)[used]), 0)   # weights of reloc_km
	out <- out[out$n_before > 0 | out$n_after > 0, , drop = FALSE]
	out <- out[order(-abs(out$d_share)), , drop = FALSE]
	rownames(out) <- NULL
	out
}

## Spatial test + label for one (s, d, window); used for the sensitivity tables.
.bpSpatialRow <- function(tab, s, d, window, alpha, class) {
	sp <- .bpSpatialTest(s, d, tab, window = window, detail = TRUE)
	data.frame(d = d, window = window,
			   comp_change = sp$comp_change, comp_threshold = .bpEmpThreshold(sp$null$comp, alpha),
			   p_comp = sp$p_comp,
			   reloc_km = sp$reloc_km, reloc_threshold = .bpEmpThreshold(sp$null$reloc, alpha),
			   p_reloc = sp$p_reloc, n_null = nrow(sp$null),
			   label = .bpSpatialLabel(class, sp$p_comp, sp$p_reloc, alpha))
}

#' Inspect one change of the breakpoint analysis
#'
#' Recomputes, for one cluster and variable, everything that decides the
#' class and label of a change, with the thresholds each quantity is
#' compared with. Give either `row` (a row of `ctx$bp`, i.e. a detected
#' change, reproduced exactly) or `at` (manual inspection, whether or not a
#' change was detected there).
#'
#' With `at` a single date, the local fit is centred on that month (the
#' start of the change is still searched within about 12 months of it, and
#' up to `d` months before for long ramps). With `at` two dates, every month
#' of the range is tried as centre and the one with the largest likelihood
#' ratio is kept (`scan` in the result).
#'
#' A manual change gets the label it would have had if it had been
#' detected; `local_test_pass` says whether it would have passed the local
#' test, and `envcpt` shows why EnvCpt did or did not propose it.
#'
#' @param ctx Output of [breakpointContext()].
#' @param row One row of `ctx$bp`.
#' @param cluster,variable Series to inspect (ignored when `row` is given).
#' @param at A date or a range of two dates (anything `as.Date()` accepts).
#' @param spatial_windows Spatial windows (months) for the sensitivity table.
#' @param envcpt Run EnvCpt on the series to report its model comparison.
#' @return A list of class `bpInspection`; see [plotInspection()].
#' @family breakpoint analysis v2
#' @export
inspectBreakpoint <- function(ctx, row = NULL, cluster = NULL, variable = NULL, at = NULL,
							  spatial_windows = c(6, 12, 24), envcpt = TRUE) {
	stopifnot(inherits(ctx, "bpContext"))
	s <- ctx$settings
	if (!is.null(row)) {
		if (nrow(row) != 1) stop("`row` must be one row of ctx$bp.")
		cluster <- row$cluster; variable <- row$variable
		mode <- "detected"
	} else {
		if (is.null(cluster) || is.null(variable) || is.null(at))
			stop("Give `row`, or `cluster`, `variable` and `at`.")
		mode <- "manual"
	}
	if (!cluster %in% ctx$clusters) stop("Unknown cluster: ", cluster)
	if (!variable %in% c("CPUE", "mean_len", "mean_hbf")) stop("Unknown variable: ", variable)
	cl_key <- as.character(cluster)
	d0 <- ctx$ts[ctx$ts$cluster == cluster, , drop = FALSE]
	pr <- .bpPrepare(d0$date, d0[[variable]], variable, d0$mean_len, s$deseason)
	if (length(pr$y) < 48) stop("Series shorter than 48 months: not analysed by the pipeline either.")
	abrupt_months <- if (variable == "mean_hbf") s$hbf_abrupt_months else NULL
	tcOf <- function(ct) .bpTimescale(pr$t, pr$ml, ct, abrupt_months, s$Linf, s$K, s$t0)

	## ---- centre of the local fit ----
	scan <- NULL
	if (mode == "detected") {
		ct <- .bpMonthIdx(row$cand_date)
	} else {
		at <- as.Date(at)
		if (length(at) == 1) {
			ct <- .bpMonthIdx(at)
		} else if (length(at) == 2) {
			cts <- seq(.bpMonthIdx(min(at)), .bpMonthIdx(max(at)))
			scan <- do.call(rbind, lapply(cts, function(c1) {
				r <- .bpCharacterise(pr$t, pr$y, c1, T_c = tcOf(c1), window = s$window)
				if (is.null(r)) NULL else data.frame(cand_date = .bpIdxToDate(c1), lr = r$lr,
													 p_nominal = r$p_nominal, s_hat = r$s_hat, d_hat = r$d_hat)
			}))
			if (is.null(scan)) stop("No month of the range can be characterised (too close to the series ends?).")
			ct <- .bpMonthIdx(scan$cand_date[which.max(scan$lr)])
		} else stop("`at` must be one date or two dates.")
	}

	## ---- local fit, timescale, class ----
	T_c <- tcOf(ct)
	ch <- .bpCharacterise(pr$t, pr$y, ct, T_c = T_c, window = s$window, detail = TRUE)
	if (is.null(ch)) stop("Cannot characterise around ", .bpIdxToDate(ct),
						  ": fewer than 24 months in the window (too close to the series ends?).")
	det <- attr(ch, "detail")
	ch$cand_date  <- .bpIdxToDate(ch$cand_t)
	ch$start_date <- .bpIdxToDate(ch$s_hat)
	ch$end_date   <- .bpIdxToDate(ch$s_hat + ch$d_hat)
	ch$delta_rel  <- if (variable == "CPUE") expm1(ch$delta) else ch$delta / abs(ch$pre_level)
	ch$class <- .bpClassify(ch, s$transient_tau_frac)
	ch$local_test_pass <- ch$p_nominal < s$alpha
	len_local <- if (is.null(abrupt_months)) stats::median(pr$ml[abs(pr$t - ct) <= 12], na.rm = TRUE) else NA_real_

	## ---- spatial test (at the estimated d and spatial window) ----
	tab <- ctx$tabs[[cl_key]]
	sp <- .bpSpatialTest(ch$s_hat, ch$d_hat, tab, window = s$spatial_window, detail = TRUE)
	comp_thr  <- .bpEmpThreshold(sp$null$comp,  s$spatial_alpha)
	reloc_thr <- .bpEmpThreshold(sp$null$reloc, s$spatial_alpha)
	label <- .bpSpatialLabel(ch$class, sp$p_comp, sp$p_reloc, s$spatial_alpha)

	## ---- mirror (computed for every change, used for the label only if composition) ----
	mir <- .bpMirror(ctx$tabs, cl_key, ch$s_hat, ch$d_hat, window = s$spatial_window)
	if (label == "composition" && !is.na(mir$share) && mir$share >= s$mirror_min) label <- "switch"
	all_flags <- sort(unique(unlist(lapply(ctx$tabs, function(x) colnames(x$n)))))
	flag_deltas <- do.call(cbind, lapply(names(ctx$tabs), function(k)
		unname(.bpFlagDelta(ctx$tabs[[k]], ch$s_hat, ch$d_hat, s$spatial_window)[all_flags])))
	flag_deltas[is.na(flag_deltas)] <- 0
	dimnames(flag_deltas) <- list(all_flags, names(ctx$tabs))
	flag_deltas <- flag_deltas[rowSums(abs(flag_deltas)) > 0, , drop = FALSE]
	flag_deltas <- flag_deltas[order(-abs(flag_deltas[, cl_key])), , drop = FALSE]

	## ---- common rule ----
	bp_all <- ctx$bp
	mi_start <- .bpMonthIdx(ch$start_date)
	near_other <- bp_all[bp_all$variable == variable & bp_all$cluster != cluster &
						 bp_all$label != "switch" &
						 abs(.bpMonthIdx(bp_all$start_date) - mi_start) <= s$common_tol_months, , drop = FALSE]
	need <- max(2, ceiling(s$common_min_share * length(ctx$clusters)))
	n_common <- length(unique(near_other$cluster)) + (label != "switch")
	common <- label != "switch" && length(ctx$clusters) >= 2 && n_common >= need

	## a detected change must come out identical
	if (mode == "detected" && (label != row$label || !isTRUE(all.equal(ch$d_hat, row$d_hat))))
		warning("Inspection of ", cluster, "/", variable, " at ", row$start_date, " gives label '", label,
				"', the pipeline gave '", row$label, "': settings of ctx differ from those used for bp?")

	## ---- sensitivity ----
	d_grid <- det$profile$d[det$profile$in_ci]
	sens_d <- do.call(rbind, lapply(d_grid, function(dd)
		.bpSpatialRow(tab, ch$s_hat, dd, s$spatial_window, s$spatial_alpha, ch$class)))
	sens_window <- do.call(rbind, lapply(spatial_windows, function(w)
		.bpSpatialRow(tab, ch$s_hat, ch$d_hat, w, s$spatial_alpha, ch$class)))
	vbLen <- function(age_m) s$Linf * (1 - exp(-s$K * (age_m / 12 - s$t0)))
	tc_bounds <- data.frame(
		class = c("abrupt", "ambiguous", "gradual"),
		T_c_range = c(sprintf("T_c >= %g", ch$d_hi), sprintf("%g <= T_c < %g", ch$d_lo, ch$d_hi),
					  sprintf("T_c < %g", ch$d_lo)),
		mean_len_range = if (is.null(abrupt_months))
			c(sprintf(">= %.1f cm", vbLen(ch$d_hi)), sprintf("%.1f to %.1f cm", vbLen(ch$d_lo), vbLen(ch$d_hi)),
			  sprintf("< %.1f cm", vbLen(ch$d_lo))) else rep("fixed T_c", 3),
		stringsAsFactors = FALSE)
	transient_tc <- if (!is.na(ch$pulse_gain) && ch$pulse_gain > 0) ch$pulse_tau / s$transient_tau_frac else NA_real_

	## ---- context ----
	same_cluster <- bp_all[bp_all$cluster == cluster &
						   abs(.bpMonthIdx(bp_all$start_date) - mi_start) <= s$window, , drop = FALSE]
	same_variable <- bp_all[bp_all$variable == variable & bp_all$cluster != cluster &
							abs(.bpMonthIdx(bp_all$start_date) - mi_start) <= 24, , drop = FALSE]
	nearest_detected <- {
		b <- bp_all[bp_all$cluster == cluster & bp_all$variable == variable, , drop = FALSE]
		if (nrow(b)) b[which.min(abs(.bpMonthIdx(b$start_date) - mi_start)), , drop = FALSE] else b
	}
	env <- NULL
	if (envcpt && s$detector == "envcpt") {
		cp <- tryCatch(do.call(.bpCandidatesEnvCpt, c(list(y = pr$y), s$detector_args)),
					   error = function(e) structure(integer(0), error = conditionMessage(e)))
		aic <- attr(cp, "aic")
		if (!is.null(aic)) {
			cps <- attr(cp, "cpts")
			env <- data.frame(model = names(aic), AIC = unname(aic), dAIC = unname(aic) - min(aic, na.rm = TRUE),
							  n_cpts = vapply(names(aic), function(m) if (m %in% names(cps)) length(cps[[m]]) else 0L,
											  integer(1)),
							  nearest_cpt = .bpIdxToDate(vapply(names(aic), function(m) {
								  v <- if (m %in% names(cps)) pr$t[cps[[m]]] else integer(0)
								  if (length(v)) v[which.min(abs(v - ct))] else NA_real_ }, numeric(1))),
							  stringsAsFactors = FALSE)
			env <- env[!is.na(env$AIC), , drop = FALSE]
			env <- env[order(env$AIC), , drop = FALSE]; rownames(env) <- NULL
		} else env <- data.frame(error = attr(cp, "error"))
	}

	## ---- the rules and their thresholds, in the order they are applied ----
	f <- function(x, d = 2) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f"))
	pulse_better <- !is.na(ch$pulse_gain) && ch$pulse_gain > 0
	rules <- data.frame(
		step = c("local test", "transient", "abrupt", "gradual", "composition", "relocation", "switch", "common"),
		statistic = c("p (LR vs no change)", "pulse decay tau (months)", "d_hi (months)", "d_lo (months)",
					  "comp_change", "reloc_km", "mirror_share", "clusters changing within tol"),
		value = c(f(ch$p_nominal, 3), f(ch$pulse_tau, 0), f(ch$d_hi, 0), f(ch$d_lo, 0),
				  f(sp$comp_change, 3), f(sp$reloc_km, 0), f(mir$share, 2), as.character(n_common)),
		threshold = c(if (is.null(scan)) sprintf("< %g", s$alpha) else
						  sprintf("< %g (best of %d centres: optimistic)", s$alpha, nrow(scan)),
					  sprintf("<= %s (%g x T_c) and pulse fits better", f(s$transient_tau_frac * ch$T_c, 1),
							  s$transient_tau_frac),
					  sprintf("<= T_c = %s", f(ch$T_c, 1)), sprintf("> T_c = %s", f(ch$T_c, 1)),
					  sprintf("> %s (p = %s < %g)", f(comp_thr, 3), f(sp$p_comp, 3), s$spatial_alpha),
					  sprintf("> %s (p = %s < %g)", f(reloc_thr, 0), f(sp$p_reloc, 3), s$spatial_alpha),
					  if (is.na(mir$share)) sprintf(">= %g (no flag change)", s$mirror_min) else
						  sprintf(">= %g (best: cluster %s)", s$mirror_min, as.character(mir$cluster)),
					  sprintf(">= %d of %d", need, length(ctx$clusters))),
		met = c(ch$local_test_pass,
				pulse_better && ch$pulse_tau <= s$transient_tau_frac * ch$T_c,
				ch$d_hi <= ch$T_c, ch$d_lo > ch$T_c,
				!is.na(sp$p_comp) && sp$p_comp < s$spatial_alpha,
				!is.na(sp$p_reloc) && sp$p_reloc < s$spatial_alpha,
				!is.na(mir$share) && mir$share >= s$mirror_min,
				common),
		stringsAsFactors = FALSE)

	structure(list(
		mode = mode, cluster = cluster, variable = variable, at = at, scan = scan,
		char = ch, T_c = T_c, len_local = len_local, hbf_fixed = !is.null(abrupt_months),
		label = label, common = common, rules = rules,
		profile = det$profile, local = det,
		series = data.frame(date = .bpIdxToDate(pr$t), raw = pr$x, y = pr$y),
		spatial = list(comp_change = sp$comp_change, p_comp = sp$p_comp, comp_threshold = comp_thr,
					   reloc_km = sp$reloc_km, p_reloc = sp$p_reloc, reloc_threshold = reloc_thr,
					   null = sp$null, window = s$spatial_window),
		sens_d = sens_d, sens_window = sens_window, tc_bounds = tc_bounds, transient_tc = transient_tc,
		flags = .bpFlagTable(tab, ch$s_hat, ch$d_hat, s$spatial_window),
		mirror = list(share = mir$share, cluster = mir$cluster, all = mir$all, flag_deltas = flag_deltas),
		context = list(same_cluster = same_cluster, same_variable = same_variable,
					   n_common = n_common, need = need, nearest_detected = nearest_detected),
		envcpt = env, settings = s),
		class = "bpInspection")
}

#' Plots of one breakpoint inspection
#'
#' @param ins Output of [inspectBreakpoint()].
#' @param ctx Output of [breakpointContext()], for the context panel (other
#'   variables of the cluster). `NULL` skips it.
#' @return A named list of ggplot objects: `series` (whole series),
#'   `local` (window of the local fit with the three fits), `profile`
#'   (duration profile against T_c), `null` (spatial null distributions),
#'   `shares` (flag shares before / after), `centroids` (flag centroid
#'   shifts) and `context` (all variables of the cluster around the change).
#' @family breakpoint analysis v2
#' @export
plotInspection <- function(ins, ctx = NULL) {
	stopifnot(inherits(ins, "bpInspection"))
	gg <- ggplot2::ggplot; a <- ggplot2::aes; th <- ggplot2::theme_bw(base_size = 9)
	ch <- ins$char
	y_lab <- if (ins$variable == "CPUE") "log CPUE (deseasoned)" else paste(ins$variable, "(deseasoned)")
	sw <- ins$spatial$window
	rects <- data.frame(xmin = c(ch$start_date, .bpIdxToDate(ch$s_hat - sw), .bpIdxToDate(ch$s_hat + ch$d_hat)),
						xmax = c(ch$end_date, ch$start_date, .bpIdxToDate(ch$s_hat + ch$d_hat + sw)),
						what = c("transition (d_hat)", "spatial: before", "spatial: after"))
	fills <- c(`transition (d_hat)` = "#fdae6b", `spatial: before` = "#c6dbef", `spatial: after` = "#6baed6")
	## plausible end of the transition, start + [d_lo, d_hi]
	d_rng <- data.frame(xmin = .bpIdxToDate(ch$s_hat + ch$d_lo), xmax = .bpIdxToDate(ch$s_hat + ch$d_hi),
						y = max(ins$local$y, na.rm = TRUE))
	win <- range(.bpIdxToDate(ins$local$t))

	p_series <- gg(ins$series, a(date, y)) +
		ggplot2::annotate("rect", xmin = win[1], xmax = win[2], ymin = -Inf, ymax = Inf, alpha = 0.12) +
		ggplot2::geom_line(linewidth = 0.3) +
		ggplot2::geom_vline(xintercept = ch$start_date, colour = "#e6550d") +
		ggplot2::labs(x = NULL, y = y_lab, title = "Whole series (grey: local window)") + th

	loc <- data.frame(date = .bpIdxToDate(ins$local$t), y = ins$local$y)
	fits <- rbind(data.frame(date = loc$date, value = ins$local$fitted_null, fit = "no change"),
				  data.frame(date = loc$date, value = ins$local$fitted_ramp, fit = "ramp"),
				  if (!is.null(ins$local$fitted_pulse))
					  data.frame(date = loc$date, value = ins$local$fitted_pulse, fit = "pulse"))
	p_local <- gg(loc, a(date, y)) +
		ggplot2::geom_rect(data = rects, a(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = what),
						   inherit.aes = FALSE, alpha = 0.3) +
		ggplot2::geom_line(linewidth = 0.3, colour = "grey30") +
		ggplot2::geom_line(data = fits, a(date, value, colour = fit), linewidth = 0.6) +
		ggplot2::geom_vline(xintercept = ch$cand_date, linetype = "dotted") +
		ggplot2::geom_segment(data = d_rng, a(x = xmin, xend = xmax, y = y, yend = y), inherit.aes = FALSE,
							  linewidth = 1.2, colour = "#e6550d") +
		ggplot2::scale_fill_manual(values = fills, breaks = names(fills), name = NULL) +
		ggplot2::scale_colour_manual(values = c(`no change` = "grey50", ramp = "#d7301f", pulse = "#41ab5d"),
									 name = NULL) +
		ggplot2::labs(x = NULL, y = y_lab, title = "Local fit") +
		th + ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

	prof <- ins$profile[is.finite(ins$profile$stat), , drop = FALSE]
	p_profile <- gg(prof, a(d, stat)) +
		ggplot2::geom_hline(yintercept = ins$local$crit, linetype = "dashed") +
		ggplot2::geom_vline(xintercept = ch$T_c, colour = "#3182bd") +
		ggplot2::geom_line() + ggplot2::geom_point(a(colour = in_ci), size = 2) +
		ggplot2::scale_colour_manual(values = c(`TRUE` = "#d7301f", `FALSE` = "grey60"),
									 name = "in [d_lo, d_hi]") +
		ggplot2::scale_x_log10(breaks = prof$d) +
		ggplot2::labs(x = "transition duration d (months, log)", y = "LR against best ramp",
					  title = sprintf("Duration profile (T_c = %.1f)", ch$T_c)) +
		th + ggplot2::theme(legend.position = "bottom")

	nl <- ins$spatial$null
	p_null <- NULL
	if (nrow(nl)) {
		nd <- rbind(data.frame(stat_name = "comp_change", value = nl$comp),
					data.frame(stat_name = "reloc_km", value = nl$reloc))
		nd <- nd[!is.na(nd$value), , drop = FALSE]
		ref <- data.frame(stat_name = c("comp_change", "reloc_km"),
						  obs = c(ins$spatial$comp_change, ins$spatial$reloc_km),
						  thr = c(ins$spatial$comp_threshold, ins$spatial$reloc_threshold))
		ref$thr[!is.finite(ref$thr)] <- NA
		p_null <- gg(nd, a(value)) +
			ggplot2::geom_histogram(bins = 30, fill = "grey70", colour = "white") +
			ggplot2::geom_vline(data = ref[!is.na(ref$thr), ], a(xintercept = thr), linetype = "dashed") +
			ggplot2::geom_vline(data = ref[!is.na(ref$obs), ], a(xintercept = obs), colour = "#d7301f") +
			ggplot2::facet_wrap(~ stat_name, scales = "free") +
			ggplot2::labs(x = NULL, y = "null positions",
						  title = "Spatial null distributions") + th
	}

	fl <- ins$flags
	sh <- rbind(data.frame(flag = fl$flag, when = "before", share = fl$share_before),
				data.frame(flag = fl$flag, when = "after", share = fl$share_after))
	sh$flag <- factor(sh$flag, levels = rev(fl$flag)); sh$when <- factor(sh$when, c("before", "after"))
	p_shares <- gg(sh, a(share, flag, fill = when)) +
		ggplot2::geom_col(position = "dodge") +
		ggplot2::scale_fill_manual(values = c(before = "#9ecae1", after = "#3182bd"), name = NULL) +
		ggplot2::labs(x = "share of the cluster's rows", y = NULL,
					  title = sprintf("Flag mix (comp_change = %.3f)", ins$spatial$comp_change)) + th

	fc <- fl[fl$used_reloc, , drop = FALSE]
	p_centroids <- NULL
	if (nrow(fc)) p_centroids <- gg(fc) +
		ggplot2::geom_segment(a(x = lon_before, y = lat_before, xend = lon_after, yend = lat_after),
							  arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"))) +
		ggplot2::geom_point(a(lon_before, lat_before), colour = "#9ecae1", size = 2) +
		ggplot2::geom_text(a(lon_after, lat_after, label = flag), size = 2.5, vjust = -0.7) +
		ggplot2::labs(x = "longitude", y = "latitude",
					  title = sprintf("Centroid shifts (reloc_km = %.0f)", ins$spatial$reloc_km)) + th

	p_context <- NULL
	if (!is.null(ctx)) {
		d0 <- ctx$ts[ctx$ts$cluster == ins$cluster, , drop = FALSE]
		d0 <- d0[d0$date >= win[1] & d0$date <= win[2], , drop = FALSE]
		cx <- do.call(rbind, lapply(c("CPUE", "mean_len", "mean_hbf"), function(v)
			data.frame(date = d0$date, value = d0[[v]], series = v)))
		cx$series <- factor(cx$series, c("CPUE", "mean_len", "mean_hbf"))
		b <- ins$context$same_cluster
		p_context <- gg(cx, a(date, value)) +
			ggplot2::geom_line(linewidth = 0.3) +
			ggplot2::geom_vline(xintercept = ch$start_date, colour = "#e6550d") +
			ggplot2::facet_wrap(~ series, ncol = 1, scales = "free_y") +
			ggplot2::labs(x = NULL, y = NULL,
						  title = "All variables of the cluster (raw scale)") + th
		if (nrow(b)) {
			b$series <- factor(b$variable, c("CPUE", "mean_len", "mean_hbf"))
			b <- b[b$start_date != ch$start_date | b$variable != ins$variable, , drop = FALSE]
			if (nrow(b)) p_context <- p_context +
				ggplot2::geom_vline(data = b, a(xintercept = start_date), linetype = "dashed")
		}
	}
	list(series = p_series, local = p_local, profile = p_profile, null = p_null,
		 shares = p_shares, centroids = p_centroids, context = p_context)
}
