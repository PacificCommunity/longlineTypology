## R/mockData.R -- simulated longline data in the seafisheries output format,
## for examples, tests and the README.


## Fishing strategies of the mock data: species composition of the catch,
## catch rate (fish per unit of effort), hooks between floats, latitude and
## longitude ranges, mean fork length, and the flags that fish them.
.MOCK_STRATEGIES <- data.frame(
	strategy = c("tropical_shallow", "tropical_deep", "subtropical_albacore", "temperate_mixed"),
	yft = c(0.55, 0.20, 0.08, 0.10),
	bet = c(0.15, 0.55, 0.07, 0.25),
	alb = c(0.05, 0.05, 0.70, 0.30),
	skj = c(0.15, 0.05, 0.02, 0.05),
	oth = c(0.10, 0.15, 0.13, 0.30),
	rate = c(2.0, 1.5, 4.0, 6.0),
	hbf = c(6, 14, 24, 30),
	lat_min = c(-2, -12, -32, -46), lat_max = c(12, 2, -20, -36),
	lon_min = c(140, 190, 160, 150), lon_max = c(185, 235, 200, 190),
	len = c(118, 122, 105, 112),
	stringsAsFactors = FALSE)

.MOCK_FLAGS <- data.frame(
	flag = c("JP", "KR", "TW", "CN", "US", "KI", "FJ", "VU", "NC", "WS", "NZ", "AU", "ES", "SB"),
	strategy = rep(1:4, c(3, 3, 4, 4)),
	weight = c(0.5, 0.3, 0.2, 0.4, 0.4, 0.2, 0.3, 0.3, 0.2, 0.2, 0.3, 0.3, 0.2, 0.2),
	stringsAsFactors = FALSE)


#' Planted changes of the mock data
#'
#' The changes [mockLonglineData()] plants with `breakpoints = TRUE`, one for
#' each label of the breakpoint analysis, spread across strategies and years.
#'
#' @return A data frame, one row per change: `type`, `strategy`, `flags`
#'   (comma-separated), `variable` (series it shows in), `start`, `duration`
#'   (months), `size` (see Details of [mockLonglineData()]) and
#'   `label` (the label the breakpoint analysis is expected to give).
#'
#' @family mock data
#' @export
mockEvents <- function() {
	data.frame(
		type     = c("step", "year_class", "trend", "relocation", "switch", "common_step"),
		strategy = c("tropical_shallow", "tropical_deep", "subtropical_albacore", "temperate_mixed",
					 "temperate_mixed", "all"),
		flags    = c("JP,KR,TW", "CN,US,KI", "FJ,VU,NC,WS", "NZ,AU,ES,SB", "AU,ES,SB", "all"),
		variable = c("CPUE", "mean_len", "mean_len", "CPUE", "CPUE", "CPUE"),
		start    = as.Date(c("2010-01-15", "2008-01-15", "2001-01-15", "2005-01-15",
							 "2014-01-15", "1995-01-15")),
		duration = c(1, 1, 72, 1, 1, 1),
		size     = c(log(0.35), 0.6, 15, log(2.5), NA, log(0.4)),
		label    = c("unexplained", "transient", "gradual", "relocation", "switch", "common"),
		stringsAsFactors = FALSE)
}


#' Simulated longline data
#'
#' Simulates operational longline data in the `seafisheries` output format
#' (one row per set), from four fishing strategies with distinct species
#' composition, hooks between floats (HBF), latitude band and fish length,
#' so that the clustering has structure to find.
#'
#' Each strategy has monthly series of catch rate, length and HBF with
#' autocorrelated noise (AR(1), phi = 0.85) and a seasonal cycle; rows add
#' their own noise (negative binomial catch, lognormal effort, normal length,
#' HBF and position).
#'
#' With `breakpoints = TRUE`, the changes of [mockEvents()] are added, one
#' for each label of the breakpoint analysis:
#' \describe{
#'   \item{step}{catch rate of every species multiplied by `exp(size)` in
#'     one month (species fractions unchanged, so rows stay in their
#'     cluster).}
#'   \item{year_class}{a strong year class: a share `size` of the catch is
#'     small fish (age 1.5 years) that grow and fade out (total mortality
#'     1.2 per year), with a matching increase in catch.}
#'   \item{trend}{mean length increased by `size` cm over `duration` months.}
#'   \item{relocation}{the flags fish 25 degrees further east, with catch rate
#'     multiplied by `exp(size)`.}
#'   \item{switch}{the flags adopt the species composition, HBF and area of
#'     the subtropical albacore strategy, keeping their own catch rate and
#'     fish length: the clustering moves them to that cluster.}
#'   \item{common_step}{catch rate multiplied by `exp(size)` in every
#'     strategy at once.}
#' }
#' The breakpoint analysis then finds them in the clusters that the
#' clustering recovers. Cluster numbers follow the clustering (ordered by
#' yellowfin catch), not the strategy order; column `mock_strategy` gives the
#' strategy behind each row.
#'
#' @param breakpoints Add the planted changes of [mockEvents()]. Default
#'   `FALSE`.
#' @param n_per_month Rows per strategy and month. Default 3 without
#'   breakpoints (a small dataset), 40 with them (the breakpoint analysis
#'   needs many rows per month).
#' @param years Years simulated. Rows before the clustering's fit years are
#'   only assigned to the clusters.
#' @param seed Seed (the caller's RNG state is restored).
#'
#' @return A data frame, one row per set, with the columns of the
#'   `seafisheries` output used by the package (`flag`, `vesselname`,
#'   `latitude`, `longitude`, `ymd`, `E`, `mean_hbf`, `<sp>_w`, `<sp>_n`
#'   for yft, bet, alb, skj and oth, `latCent`, `lonCent`, `latAbs`,
#'   `hbf_outlier`, `hbf_cat`, `mean_len`, `len_source`) and `mock_strategy`.
#'   With `breakpoints = TRUE`, attribute `"events"` holds [mockEvents()].
#'
#' @examples
#' mock <- mockLonglineData()
#' table(mock$mock_strategy)
#'
#' \dontrun{
#' ## With planted changes: see the README for the clustering and the reports
#' mock_bp <- mockLonglineData(breakpoints = TRUE)
#' attr(mock_bp, "events")
#' }
#'
#' @family mock data
#' @export
mockLonglineData <- function(breakpoints = FALSE, n_per_month = if (breakpoints) 25 else 3,
							 years = 1985:2022, seed = 4) {
	.withSeed(seed, .mockRows(breakpoints, n_per_month, years))
}


.mockRows <- function(breakpoints, n_per_month, years) {
	S  <- .MOCK_STRATEGIES; Fl <- .MOCK_FLAGS
	sp <- c("yft", "bet", "alb", "skj", "oth")
	dates <- seq(as.Date(paste0(min(years), "-01-15")), as.Date(paste0(max(years), "-12-15")),
				 by = "month")
	nt <- length(dates); ns <- nrow(S); t <- seq_len(nt)
	moy <- as.POSIXlt(dates)$mon

	## Monthly latent series per strategy: AR(1) + seasonal cycle
	ar1 <- function(sd, phi = 0.85) {
		x <- numeric(nt); x[1] <- stats::rnorm(1, 0, sd)
		e <- stats::rnorm(nt, 0, sd * sqrt(1 - phi^2))
		for (i in 2:nt) x[i] <- phi * x[i - 1] + e[i]
		x
	}
	seas <- function(amp) amp * sin(2 * pi * moy / 12 + stats::runif(1, 0, 2 * pi))
	lrate <- sapply(seq_len(ns), function(s) log(S$rate[s]) + ar1(0.10) + seas(0.08))
	len   <- sapply(seq_len(ns), function(s) S$len[s] + ar1(1.0) + seas(0.4))
	hbf   <- sapply(seq_len(ns), function(s) S$hbf[s] + ar1(0.6) + seas(0.25))

	## Rows: n_per_month per strategy and month, flags drawn by weight
	g  <- expand.grid(t = t, s = seq_len(ns))
	g  <- g[rep(seq_len(nrow(g)), each = n_per_month), ]
	fl <- vapply(g$s, function(s) {
		f <- Fl[Fl$strategy == s, ]
		f$flag[sample.int(nrow(f), 1, prob = f$weight)]
	}, character(1))
	n  <- nrow(g)
	s_home <- g$s; s_eff <- g$s           # strategy of the flag / strategy it fishes
	lrate_row <- lrate[cbind(g$t, s_home)]
	len_row   <- len[cbind(g$t, s_home)]
	lon_shift <- numeric(n)
	yc_share  <- numeric(n); yc_len <- rep(NA_real_, n)

	ev <- if (breakpoints) mockEvents() else NULL
	for (i in seq_len(if (is.null(ev)) 0 else nrow(ev))) {
		e  <- ev[i, ]
		s0 <- match(e$start, dates)
		on <- if (e$flags == "all") rep(TRUE, n) else fl %in% strsplit(e$flags, ",")[[1]]
		ramp <- pmin(pmax((g$t - s0 + 1) / e$duration, 0), 1)
		after <- g$t >= s0
		switch(e$type,
			step = , common_step = {
				lrate_row[on] <- lrate_row[on] + e$size * ramp[on]
			},
			trend = {
				len_row[on] <- len_row[on] + e$size * ramp[on]
			},
			relocation = {
				lrate_row[on & after] <- lrate_row[on & after] + e$size
				lon_shift[on & after] <- 25
			},
			switch = {
				s_eff[on & after] <- match("subtropical_albacore", S$strategy)
			},
			year_class = {
				dt <- (g$t - s0) / 12
				p  <- ifelse(after, e$size * exp(-1.2 * pmax(dt, 0)), 0)
				yc_share[on] <- p[on]
				yc_len[on]   <- 150.3 * (1 - exp(-0.442 * (1.5 + pmax(dt[on], 0) + 0.244)))
				lrate_row[on] <- lrate_row[on] + log1p(p[on])
			})
	}

	## Effort, catch by species, length, HBF, positions
	E <- stats::rlnorm(n, 4, 0.3)
	frac <- as.matrix(S[s_eff, sp])
	frac <- frac * matrix(stats::rlnorm(n * length(sp), 0, 0.15), n)
	frac <- frac / rowSums(frac)
	mu <- exp(lrate_row) * E
	catch <- matrix(stats::rnbinom(n * length(sp), size = 20, mu = mu * frac), n,
					dimnames = list(NULL, sp))

	mean_len <- (1 - yc_share) * len_row + ifelse(yc_share > 0, yc_share * yc_len, 0) +
		stats::rnorm(n, 0, 1.5)
	mean_hbf <- pmax(hbf[cbind(g$t, s_eff)] + stats::rnorm(n, 0, 1), 3)
	latitude  <- stats::runif(n, S$lat_min[s_eff], S$lat_max[s_eff])
	longitude <- (stats::runif(n, S$lon_min[s_eff], S$lon_max[s_eff]) + lon_shift) %% 360

	## Weight (tonnes) from length: yellowfin length-weight relationship
	## for every species, enough for the mock data.
	w1 <- 1.96e-5 * mean_len^2.98 / 1000
	out <- data.frame(flag = fl,
					  vesselname = paste0(fl, "_V", sample.int(5, n, replace = TRUE)),
					  latitude = latitude, longitude = longitude,
					  ymd = dates[g$t] + sample(-14:13, n, replace = TRUE),
					  E = E, mean_hbf = mean_hbf, stringsAsFactors = FALSE)
	for (x in sp) {
		out[[paste0(x, "_w")]] <- catch[, x] * w1
		out[[paste0(x, "_n")]] <- catch[, x]
		out[[paste0(x, "_n_est")]] <- NA_integer_
	}
	out$latCent <- floor(out$latitude) + 0.5
	out$lonCent <- floor(out$longitude) + 0.5
	out$latAbs  <- abs(out$latitude)

	## Data quality columns: a few HBF outliers, length sources
	outlier <- stats::runif(n) < 0.02
	out$mean_hbf[outlier] <- stats::runif(sum(outlier), 40, 60)
	out$hbf_outlier <- ifelse(outlier, "outlier", "valid")
	out$hbf_cat <- as.integer(cut(out$mean_hbf, c(-Inf, 10, 20, Inf)))
	out$mean_len <- mean_len
	src <- sample(c("operational", "aggregated", "imputed"), n, replace = TRUE, prob = c(0.8, 0.1, 0.1))
	imp <- src == "imputed"
	out$mean_len[imp] <- S$len[s_home[imp]] + stats::rnorm(sum(imp), 0, 4)
	out$len_source <- src
	out$mock_strategy <- S$strategy[s_eff]

	out <- out[order(out$ymd), , drop = FALSE]
	rownames(out) <- NULL
	if (breakpoints) attr(out, "events") <- ev
	out
}
