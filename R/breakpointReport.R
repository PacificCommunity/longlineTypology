## R/breakpointReport.R -- breakpoint analysis of a clustering: data set-up,
## the tables behind one change (rows per window, flag decomposition, same
## windows in the other clusters), the overview plot, and the PDF report.

utils::globalVariables(c("label", "end_plot", "start_date"))


# ---- data ----------------------------------------------------------------------

#' Rows of a clustering prepared for the breakpoint analysis
#'
#' Rebuilds the clustered data ([clusterData()]), keeps the rows of the fit
#' years ([fitRows()]), and adds the columns the breakpoint functions use:
#' `cluster` (ordered labels of [clusterLabels()] for K = `k`), `date`
#' (copy of `ymd`) and `year`.
#'
#' Imputed values are set to `NA` (unless `use_imputed = TRUE`): lengths with
#' `len_source == "imputed"` and HBF with `hbf_outlier != "valid"`. These are
#' the values [aggregateClusterTimeSeries()] leaves out of the series, so
#' every table computed from these rows describes the values behind the
#' series. The rows themselves are kept: CPUE, flags and positions use them.
#'
#' @inheritParams fitRows
#' @param k Number of clusters (one of `names(diagnostics$kmeans_full)`).
#' @param data `NULL` (default) to reload the data from
#'   `diagnostics$data_path`, or the data given to [runClustering()].
#' @param use_imputed Keep imputed lengths and HBF. Default `FALSE`.
#'
#' @return A data frame: the fit-year rows with `cluster`, `date` and `year`.
#'
#' @family breakpoint report
#' @export
breakpointData <- function(diagnostics, k, data = NULL, use_imputed = FALSE) {
	diagnostics <- .asDiagnostics(diagnostics)
	data <- clusterData(diagnostics, data)
	data$cluster <- clusterLabels(diagnostics, k)
	data <- data[fitRows(diagnostics, data), , drop = FALSE]
	data$date <- data$ymd
	data$year <- as.integer(format(data$ymd, "%Y"))
	if (!use_imputed) {
		if ("len_source" %in% names(data))  data$mean_len[data$len_source %in% "imputed"] <- NA
		if ("hbf_outlier" %in% names(data)) data$mean_hbf[!data$hbf_outlier %in% "valid"] <- NA
	}
	rownames(data) <- NULL
	data
}


## Value of a variable over a set of rows, computed as in the series:
## CPUE = sum(yft_n) / sum(E); mean_len and mean_hbf = mean over rows.
.bpValue <- function(d, v) {
	if (!nrow(d)) return(NA_real_)
	if (v == "CPUE") sum(d$yft_n, na.rm = TRUE) / sum(d$E, na.rm = TRUE) else mean(d[[v]], na.rm = TRUE)
}


# ---- one change ------------------------------------------------------------------

#' Rows before, during and after a change
#'
#' Splits the rows of one cluster (or of all clusters) into the windows used
#' by the breakpoint analysis: `before` (the `window` months before the
#' change starts), `during` (the transition) and `after` (the `window`
#' months after it ends). The window is the spatial window of the analysis
#' (`spatial_window`, 12 months by default).
#'
#' @param data Output of [breakpointData()].
#' @param ins Output of [inspectBreakpoint()].
#' @param cluster Cluster whose rows are used; `NULL` uses all clusters.
#'   Default: the cluster of the change.
#'
#' @return A list with the rows of each window (`before`, `during`, `after`)
#'   and `per_month`: one row per window with the number of rows, the
#'   number of months, and the minimum, median and maximum rows per month
#'   and the number of empty months, counting only rows with a value of the
#'   variable (all rows for CPUE).
#'
#' @family breakpoint report
#' @export
breakpointWindows <- function(data, ins, cluster = ins$cluster) {
	stopifnot(inherits(ins, "bpInspection"))
	s <- ins$char$s_hat; d <- ins$char$d_hat; w <- ins$settings$spatial_window
	dd <- if (is.null(cluster)) data else data[data$cluster == cluster, , drop = FALSE]
	mi <- .bpMonthIdx(dd$date)
	win <- list(before = c(s - w, s), during = c(s, s + d), after = c(s + d, s + d + w))
	rows <- lapply(win, function(r) dd[mi >= r[1] & mi < r[2], , drop = FALSE])
	v <- ins$variable
	per_month <- do.call(rbind, lapply(names(win), function(nm) {
		r <- win[[nm]]
		x <- rows[[nm]]
		has <- if (v == "CPUE") rep(TRUE, nrow(x)) else !is.na(x[[v]])
		cnt <- tabulate(.bpMonthIdx(x$date[has]) - r[1] + 1, nbins = r[2] - r[1])
		data.frame(window = nm, rows = nrow(x), rows_with_value = sum(has), months = r[2] - r[1],
				   min_per_month = if (length(cnt)) min(cnt) else NA,
				   median_per_month = if (length(cnt)) stats::median(cnt) else NA,
				   max_per_month = if (length(cnt)) max(cnt) else NA,
				   empty_months = sum(cnt == 0))
	}))
	c(rows, list(per_month = per_month))
}


#' Split a change between groups: composition and within-group parts
#'
#' Decomposes the difference of a variable between two sets of rows (e.g.
#' the windows before and after a change, [breakpointWindows()]) into a
#' composition part (the shares of the groups changing) and a within-group
#' part (the values of the groups changing):
#' `total = composition + within`, exactly.
#'
#' For each group, `w` is its share (of effort for CPUE, of rows with a value
#' otherwise) and `m` its value (CPUE, or mean). Then
#' `composition = (w_after - w_before) * (m_before + m_after) / 2` and
#' `within = (w_before + w_after) / 2 * (m_after - m_before)`. A group
#' present in one window only counts entirely as composition.
#'
#' @param before,after Data frames of rows, with the variable and the group
#'   column (and `yft_n` and `E` for CPUE).
#' @param variable `"CPUE"`, `"mean_len"` or `"mean_hbf"`.
#' @param group Grouping column. Default `"flag"`.
#'
#' @return A list: `table` (one row per group: rows, shares and values
#'   before and after, composition and within parts, sorted by the size of
#'   their contribution), `total`, `composition` and `within`.
#'
#' @family breakpoint report
#' @export
decomposeChange <- function(before, after, variable, group = "flag") {
	groups <- sort(unique(c(as.character(before[[group]]), as.character(after[[group]]))))
	groups <- groups[!is.na(groups)]
	st <- function(d) {
		g <- factor(as.character(d[[group]]), levels = groups)
		if (variable == "CPUE") {
			wt <- tapply(d$E, g, sum, na.rm = TRUE)
			m  <- tapply(d$yft_n, g, sum, na.rm = TRUE) / wt
		} else {
			ok <- !is.na(d[[variable]])
			m  <- tapply(d[[variable]][ok], g[ok], mean)
			wt <- tapply(ok, g, sum)
		}
		wt[is.na(wt)] <- 0
		m[!is.finite(m)] <- NA
		data.frame(rows = as.numeric(table(g)), w = as.numeric(wt) / max(sum(wt), 1e-12),
				   m = as.numeric(m))
	}
	b <- st(before); a <- st(after)
	in_b <- b$w > 0 & !is.na(b$m); in_a <- a$w > 0 & !is.na(a$m)
	m_bar  <- ifelse(in_b & in_a, (b$m + a$m) / 2, ifelse(in_a, a$m, b$m))
	comp   <- (a$w - b$w) * m_bar
	comp[is.na(comp)] <- 0
	within <- ifelse(in_b & in_a, (b$w + a$w) / 2 * (a$m - b$m), 0)
	tab <- data.frame(group = groups, rows_before = b$rows, rows_after = a$rows,
					  share_before = b$w, share_after = a$w, value_before = b$m, value_after = a$m,
					  composition = comp, within = within, stringsAsFactors = FALSE)
	tab <- tab[order(-(abs(tab$composition) + abs(tab$within))), , drop = FALSE]
	rownames(tab) <- NULL
	names(tab)[1] <- group
	list(table = tab, total = .bpValue(after, variable) - .bpValue(before, variable),
		 composition = sum(comp), within = sum(within))
}


#' The windows of a change in every cluster
#'
#' Value of the variable in the windows before and after a change
#' ([breakpointWindows()]), in every cluster and in all clusters pooled.
#'
#' @inheritParams breakpointWindows
#'
#' @return A data frame: `cluster` (`"all"` for the pooled rows), rows and
#'   value before and after, `change` and `rel_change` (change relative to
#'   the value before).
#'
#' @family breakpoint report
#' @export
changeAcrossClusters <- function(data, ins) {
	v <- ins$variable
	cls <- sort(unique(data$cluster))
	do.call(rbind, lapply(c(as.list(cls), list(NULL)), function(cl) {
		w <- breakpointWindows(data, ins, cluster = cl)
		vb <- .bpValue(w$before, v); va <- .bpValue(w$after, v)
		data.frame(cluster = if (is.null(cl)) "all" else as.character(cl),
				   rows_before = nrow(w$before), rows_after = nrow(w$after),
				   before = vb, after = va, change = va - vb, rel_change = (va - vb) / abs(vb))
	}))
}


# ---- overview plot -----------------------------------------------------------------

#' Plot every cluster series with the changes found
#'
#' One panel per cluster and variable, the changes of
#' [analyseClusterBreakpoints()] shaded from their start to their end (at least 9
#' months wide, so that short changes stay visible) and
#' coloured by label (black for changes marked `common`).
#'
#' @param ts Output of [aggregateClusterTimeSeries()].
#' @param bp Output of [analyseClusterBreakpoints()] (may have 0 rows).
#' @param variables Variables shown.
#'
#' @return A ggplot object.
#'
#' @family breakpoint report
#' @export
plotClusterBreakpoints <- function(ts, bp, variables = c("CPUE", "mean_len", "mean_hbf")) {
	long <- do.call(rbind, lapply(variables, function(v)
		data.frame(cluster = ts$cluster, variable = v, date = ts$date, value = ts[[v]])))
	long$variable <- factor(long$variable, variables)
	long$cluster  <- factor(long$cluster, sort(unique(ts$cluster)))
	p <- ggplot2::ggplot(long, ggplot2::aes(date, value))
	if (!is.null(bp) && nrow(bp) > 0) {
		bp <- bp[bp$variable %in% variables, , drop = FALSE]
		bp$lab <- if ("common" %in% names(bp)) ifelse(bp$common, "common", bp$label) else bp$label
		bp$end_plot <- pmax(bp$end_date, bp$start_date + 270)   # short changes drawn 9 months wide
		bp$variable <- factor(bp$variable, variables)
		bp$cluster  <- factor(bp$cluster, levels(long$cluster))
		p <- p + ggplot2::geom_rect(data = bp, ggplot2::aes(xmin = start_date, xmax = end_plot,
															ymin = -Inf, ymax = Inf, fill = lab),
									alpha = 0.5, inherit.aes = FALSE) +
			ggplot2::scale_fill_manual(values = .BP_COLS, name = "Label")
	}
	p + ggplot2::geom_line(linewidth = 0.25, na.rm = TRUE) +
		## facet_wrap, not facet_grid: each panel gets its own y scale
		ggplot2::facet_wrap(~ cluster + variable, ncol = length(variables), scales = "free_y",
							labeller = function(lab) list(paste0("cluster ", lab$cluster, ": ", lab$variable))) +
		ggplot2::labs(x = NULL, y = NULL) +
		ggplot2::theme_bw(base_size = 8) +
		ggplot2::theme(legend.position = "bottom")
}


# ---- report ------------------------------------------------------------------------

#' Breakpoint report of a clustering
#'
#' Renders a PDF report of the breakpoint analysis for one value of K: all
#' cluster series with the changes found, a table of the changes, and for
#' each selected change the plots and tables behind its label (data in each
#' window, timing, flag mix and positions, rows moving between clusters,
#' the same windows in the other clusters). A section explains how to read
#' each plot and table; the report itself draws no conclusion. Needs
#' `rmarkdown` and a LaTeX installation.
#'
#' @inheritParams breakpointData
#' @param cards Labels of the changes given a detailed section, or `"all"`.
#'   Default `c("switch", "unexplained", "ambiguous")`.
#' @param include_common Also give a detailed section to changes marked
#'   `common`. Default `FALSE`.
#' @param manual Optional list of changes to inspect by hand, each a list
#'   with `cluster`, `variable` and `at` (one date, or two dates for a range;
#'   see [inspectBreakpoint()]). Works whether or not a change was found.
#' @param detail Add the EnvCpt model comparison to each section (slower).
#'   Default `FALSE`.
#' @param top_n Flags listed in the tables (the others are summed).
#' @param output_file Path of the PDF.
#' @param title Title of the report.
#' @param quiet Passed to [rmarkdown::render()].
#' @param ... Settings of the analysis, passed to [breakpointContext()]
#'   (e.g. `window`, `alpha`, `spatial_window`).
#'
#' @return The path of the PDF, invisibly.
#'
#' @examples
#' \dontrun{
#' breakpointReport("cluster_outputs/cluster_diagnostics.RDS", k = 7,
#'                  output_file = "breakpoint_report_k7.pdf")
#' }
#'
#' @family breakpoint report
#' @export
breakpointReport <- function(diagnostics, k, data = NULL,
							 cards = c("switch", "unexplained", "ambiguous"), include_common = FALSE,
							 manual = NULL, detail = FALSE, top_n = 8,
							 output_file = "breakpoint_report.pdf",
							 title = "Breakpoint report", quiet = TRUE, ...) {
	diagnostics <- .asDiagnostics(diagnostics)
	k <- as.character(k)
	if (length(k) != 1 || !k %in% names(diagnostics$kmeans_full))
		stop("k must be one of: ", paste(names(diagnostics$kmeans_full), collapse = ", "))
	if (as.numeric(k) < 2) stop("The breakpoint report needs K >= 2.")
	for (m in manual)
		if (!all(c("cluster", "variable", "at") %in% names(m)))
			stop("Each element of `manual` needs cluster, variable and at.")

	bp_data <- breakpointData(diagnostics, k, data)
	ts  <- aggregateClusterTimeSeries(bp_data)
	ctx <- breakpointContext(ts, bp_data, ...)
	bp  <- ctx$bp
	bp  <- bp[order(bp$cluster, match(bp$variable, c("CPUE", "mean_len", "mean_hbf")), bp$start_date), ,
			  drop = FALSE]
	rownames(bp) <- NULL
	sel <- if (identical(cards, "all")) rep(TRUE, nrow(bp)) else bp$label %in% cards
	if (!include_common) sel <- sel & !bp$common

	.renderReport("breakpoint_report.Rmd", output_file, title, quiet, list(
		diagnostics = diagnostics, k = k, bp_data = bp_data, ts = ts, ctx = ctx, bp = bp,
		card_rows = which(sel), cards = cards, include_common = include_common,
		manual = manual, detail = detail, top_n = top_n))
}


## Render a report template of inst/reports with the given objects in its
## environment. Rendered from a temporary copy of the folder, since the
## installed package may not be writable.
.renderReport <- function(template, output_file, title, quiet, objects) {
	if (!requireNamespace("rmarkdown", quietly = TRUE))
		stop("Reports need the 'rmarkdown' package.")
	src <- system.file("reports", template, package = "longlineTypology")
	if (src == "") stop("Report template not found in the installed package: ", template)
	work <- tempfile("report_")
	dir.create(work)
	file.copy(list.files(dirname(src), full.names = TRUE), work)
	env <- list2env(objects, envir = new.env(parent = globalenv()))
	output_file <- normalizePath(output_file, mustWork = FALSE)
	out <- rmarkdown::render(file.path(work, template),
							 output_file = basename(output_file),
							 output_dir = dirname(output_file),
							 intermediates_dir = work,
							 params = list(title = title),
							 envir = env, quiet = quiet)
	invisible(out)
}
