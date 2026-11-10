## Row-level data for the breakpoint tests: two clusters, two flags each; a
## large CPUE step in cluster 1 (month 180) and flag "B" moving from cluster 2
## to cluster 1 (month 240).

mkRows <- function(seed = 1) {
	set.seed(seed)
	dates <- seq(as.Date("1990-01-15"), by = "month", length.out = 360)
	one <- function(flag, cluster, lon, cpue, keep = rep(TRUE, 360)) {
		n <- 20
		d <- data.frame(flag = flag, cluster = cluster, date = rep(dates, each = n),
						longitude = rep(lon, 360 * n) + stats::rnorm(360 * n),
						latitude = stats::rnorm(360 * n), E = 100,
						yft_n = stats::rpois(360 * n, rep(cpue, each = n) * 100),
						mean_len = stats::rnorm(360 * n, 117, 2), mean_hbf = stats::rnorm(360 * n, 20, 1))
		d[rep(keep, each = n), ]
	}
	step <- ifelse(seq_len(360) >= 180, 0.2, 1)
	moved <- seq_len(360) >= 240
	rows <- rbind(one("A", 1, 170, step), one("B", 1, 175, 0.5, moved),
				  one("B", 2, 175, 0.5, !moved), one("C", 2, 200, 0.5))
	rows$ymd <- rows$date
	rows
}
