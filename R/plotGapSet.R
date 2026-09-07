#' Overlay the gap profiles from a set of replicates and show the selection rule
#'
#' @param gap_stat list of gap tables (the `gap_stat` list built by the
#'   clustering vignette); each element needs columns `gap` and `SE.sim`
#'   with rownames giving k.
#' @param SE.factor as passed to cluster::maxSE(); 1 matches customKmeans().
#' @param file optional png path; if NULL, draws to the current device.
#' @return invisibly, the selected k per replicate and their frequency table.
#' @export
## ---- plotGapSet(): overlay the 15 gap profiles and show the selection rule ----
## gap_stat is the list built by the vignette loop; each element has columns
## logW, E.logW, gap, SE.sim with rownames giving k. Base graphics only --
## no continuous colour scale, so it is unaffected by the broken farver.
plotGapSet <- function(gap_stat, SE.factor = 1, file = NULL) {
  if (!is.null(file)) png(file, 1100, 520, res = 110)
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 1))
  ks <- as.integer(rownames(gap_stat[[1]]))
  G  <- sapply(gap_stat, function(d) d[, "gap"])
  S  <- sapply(gap_stat, function(d) d[, "SE.sim"])

  ## panel 1 -- the gap curves themselves
  matplot(ks, G, type = "l", lty = 1, col = "#00000030", lwd = 1.5,
          xlab = "k", ylab = "Gap(k)", main = "Gap profiles, 15 replicates")
  lines(ks, rowMeans(G), col = "firebrick", lwd = 2.5)
  arrows(ks, rowMeans(G) - rowMeans(S), ks, rowMeans(G) + rowMeans(S),
         angle = 90, code = 3, length = 0.03, col = "firebrick")

  ## panel 2 -- the Tibshirani rule made explicit.
  ## Selected k = FIRST k with  Gap(k) - (Gap(k+1) - SE.factor*s_{k+1})  >= 0
  D <- G[-nrow(G), , drop = FALSE] -
       (G[-1, , drop = FALSE] - SE.factor * S[-1, , drop = FALSE])
  matplot(ks[-length(ks)], D, type = "l", lty = 1, col = "#00000030", lwd = 1.5,
          xlab = "k", ylab = "Gap(k) - [Gap(k+1) - s(k+1)]",
          main = "Selection criterion (first k above 0 wins)")
  abline(h = 0, col = "firebrick", lwd = 2, lty = 2)
  sel <- apply(D, 2, function(z) { i <- which(z >= 0); if (length(i)) ks[i[1]] else NA })
  points(sel, rep(0, length(sel)), pch = 19, col = "firebrick", cex = 1.1)
  par(op); if (!is.null(file)) dev.off()
  invisible(list(selected = sel, table = table(sel, useNA = "ifany")))
}

