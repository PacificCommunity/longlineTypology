library(ggplot2)
for (f in list.files(pattern=".R", full.names = TRUE)) source(f)

EC_clean <- readRDS("~/WORK/Team/Romain/Fisheries-data/2026/LL_1x1_imp_raised_customLat.RDS")

stat_list <- list()
stat_list$nrow_total <- nrow(EC_clean)



scenario <- "yba_lat"
scenario <- "yba_hbf_lat"
EC_clean$total_n      <- EC_clean$yft_n + EC_clean$alb_n + EC_clean$bet_n
EC_clean$yft_fraction <- ifelse(EC_clean$total_n > 0, EC_clean$yft_n / EC_clean$total_n, 0)
EC_clean$alb_fraction <- ifelse(EC_clean$total_n > 0, EC_clean$alb_n / EC_clean$total_n, 0)
EC_clean$bet_fraction <- ifelse(EC_clean$total_n > 0, EC_clean$bet_n / EC_clean$total_n, 0)
names(EC_clean)[names(EC_clean) == "hbf"] <- "mean_hbf"
EC_clean <- EC_clean[!is.na(EC_clean$mean_hbf), , drop = FALSE]
summary(rowSums(EC_clean[, c("yft_fraction","bet_fraction","alb_fraction")]))

#scenario <- "sp_lat"
#scenario <- "sp"
#scenario <- "sp_hbf_lat"
#all_n <- with(EC_clean, yft_n + bet_n + alb_n + skj_n + oth_n)
#EC_clean$yft_fraction <- ifelse(all_n > 0, EC_clean$yft_n / all_n, 0)
#EC_clean$bet_fraction <- ifelse(all_n > 0, EC_clean$bet_n / all_n, 0)
#EC_clean$alb_fraction <- ifelse(all_n > 0, EC_clean$alb_n / all_n, 0)
#EC_clean$oth_fraction <- ifelse(all_n > 0, (EC_clean$oth_n + EC_clean$skj_n) / all_n, 0)
# sanity: should be exactly 1 everywhere
#summary(rowSums(EC_clean[, c("yft_fraction","bet_fraction","alb_fraction","oth_fraction")]))


# PCA and kmeans parameters
pca_variance_threshold <- 0.7   # Retain PCs explaining 70% of variance
sample_no              <- 2000 # Sample size per replicate for gap statistic
kmeans_set             <- 15    # Number of kmeans replicates
max_k                  <- 10    # Maximum number of clusters to test
iter_max               <- 1e6   # Maximum iterations for kmeans
nstart                 <- 3     # Number of random starts for kmeans
B_var                  <- 25    # Bootstrap replicates for gap statistic
nb_cores			   <- 5

# Cluster merging parameters
catch_threshold <- 3            # Merge clusters with < 3% of total catch

# Select and prepare variables
pca_var <- pca_variance_threshold * 100
pca_select <- selectScenario(scenario)

cat("Variables selected for clustering:\n")
print(pca_select)

select_dat <- EC_clean[, pca_select]
pca_scaled <- pcaScale(select_dat, method = "zscore")

pca_res <- customPCA(pca_scaled,
				 variance_threshold = pca_variance_threshold,
				 print_it = FALSE)

cat("PCA retained", pca_res$no_var, "components explaining >=",
	pca_variance_threshold * 100, "% of variance\n\n")

pca_full <- pca_res$x[, 1:pca_res$no_var]

indices_90  <- which(as.integer(format(EC_clean$date, "%Y")) >= 1990)
pca_full_90 <- pca_full[indices_90, ]

cat("Using", nrow(pca_full_90), "observations from 1990 onwards for clustering\n")
stat_list$nrow_clustering <- nrow(pca_full_90)


#Gap-statistic
set.seed(123)
sample_indices <- lapply(1:kmeans_set, function(x) {
	sample(1:nrow(pca_full_90), size = sample_no, replace = TRUE)
})

cat("\nRunning", kmeans_set, "kmeans replicates on samples of",
	sample_no, "observations...\n")

kmeans_res <- list()
gap_stat   <- list()

for (i in 1:kmeans_set) {
	message("i=",i)
	pca_sub <- pca_full_90[sample_indices[[i]], 1:pca_res$no_var]

	res <- customKmeans(pca_sub,
						max_k      = max_k,
						random_set = B_var,
						iter_max   = iter_max,
						nstart     = nstart,
						d.power    = 2,
						print_it   = FALSE)

	kmeans_res[[i]] <- res$kmeans
	gap_stat[[i]]   <- res$gap_stat
	plotGapSet(gap_stat, file = paste0("figs/gap_",scenario,"_pca",pca_var,"_",sample_no,".png"))
}

# Determine most common K across replicates
cluster_counts <- data.frame(
	cluster_n = sapply(kmeans_res, function(x) length(unique(x$cluster))),
	replicate = 1:kmeans_set
)

freq_table   <- table(cluster_counts$cluster_n)
main_clust_n <- as.numeric(names(freq_table)[which.max(freq_table)])

# Store as data frame for use in diagnostic report
freq_df <- data.frame(
	n_clusters = as.numeric(names(freq_table)),
	freq       = as.numeric(freq_table)
)

cat("\nCluster number distribution across replicates:\n")
print(freq_table)
cat("\nSelected K =", main_clust_n, "clusters\n")
stat_list$selected_k <- main_clust_n


