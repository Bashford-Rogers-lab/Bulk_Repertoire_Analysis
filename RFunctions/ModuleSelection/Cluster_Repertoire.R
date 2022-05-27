## Function to Cluster filtered matrix 
## Lauren Overend and Rachael Bashford-Rogers
## lauren.overend@oriel.ox.ac.uk
## Jan 2022 
# used tutorial from https://uc-r.github.io/kmeans_clustering#silo and https://uc-r.github.io/hc_clustering
library(dendextend)
library(cluster)
library(Rfast)
options(bitmapType='cairo-png')

cluster_features <- function(outputdir, correlation_matrix, type, iso_type, data_type){
	mat_feature_corr <- correlation_matrix
	
	w=10
	pdf(paste0(outputdir, "Plots/CLUSTERING_isotyper_metrics_FINAL_METRICS__CorrelatedvsAnticorrelated", type, "_", iso_type, ".pdf"), height=w*1*2, width=w*4.5*1.5)
	par(mfrow= c(2,1), mar = c(5,5,3,3))
	#### build distance tree
	## option 1: anticorrelated metrics will be more distant than correlated
	dist1 = as.dist(2-(mat_feature_corr))
	plot(hclust(dist1), main = "option 1", hang = -1)
	## option 2: anticorrelated metrics will be same as correlated (i.e. distance based on abs(correlation))
	dist2 = as.dist(2-abs(mat_feature_corr))
	plot(hclust(dist2), main = "option 2", hang = -1)
	dev.off()
	

	## Cluster features
	## Using anti-correlated tree 
	## this is the one we will take forward 
	dist1 = as.dist(2-(mat_feature_corr))
	
	## What aglomeration method to use for clustering 
	## default of hclus is complete 
	m <- c( "average", "single", "complete", "ward", "weighted")
	names(m) <- c( "average", "single", "complete", "ward", "weighted")

	print(paste0("Determining which cluster aglomeration method to use"))
	# function to compute coefficient
	ac <- function(x) {
	  agnes(dist1, method = x)$ac
	}
	# Best ac score is best agglomerative coefficient 
	# Amoung of clustering found!
	ac_scores <- data.frame(map_dbl(m, ac))
	colnames(ac_scores) <- "agglomerative_coefficient"
	ac_scores$Method <- rownames(ac_scores)
	
	## Plotting cluster scores per method
	pdf(paste0(outputdir, "Plots/CLUSTERING_isotyper_metrics_FINAL_METRICS__AC_score", type, "_", iso_type, ".pdf"), height=5, width=5)
	plot(ggplot(ac_scores, aes(x=Method, y=agglomerative_coefficient, fill=Method)) +geom_col() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Cluster Agglomeration Method")+ ylab("Agglomerative Coefficient"))
	dev.off()
	
	pdf(paste0(outputdir, "Plots/CLUSTERING_isotyper_metrics_FINAL_METRICS__AC_score", type, "_", iso_type, ".pdf"), height=5, width=5)
	plot(ggplot(ac_scores, aes(x=Method, y=agglomerative_coefficient, fill=Method)) +geom_col() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Cluster Agglomeration Method")+ ylab("Agglomerative Coefficient"))
	dev.off()
	
	## Select Method
	method_ac_use <- ac_scores$Method[ac_scores$agglomerative_coefficient ==max(ac_scores$agglomerative_coefficient)]
	
	##rename for using hclust
	if(method_ac_use=="ward"){
		method_ac_use="ward.D2"
	}
	print(paste0("Cluster Agglomeration Method: ",  method_ac_use))
	
	
	print(paste0("Determining Optimal Number of clusters:"))
	## manual silhouett
	avg_sil <- function(k) {
		  km.res <- kmeans(mat_feature_corr, centers = k, nstart = 20)
		  ss <- silhouette(km.res$cluster, dist(mat_feature_corr))
		  mean(ss[, 3])
	}
	
	## Specify number of K 
	## Must be up to 1-number of branches of data 
	## Will use this in the paramater later
	k.max_use <- (dim(mat_feature_corr)[1]-1)
	if(k.max_use >40){
		k.max_use <- 40
	}
	k.values <- 2:k.max_use

	
	# extract avg silhouette for 2-15 clusters
	avg_sil_values <- map_dbl(k.values, avg_sil)
	
	pdf(paste0(outputdir, "Plots/CLUSTERING_isotyper_metrics_FINAL_METRICS__Number_Clusters", type, "_", iso_type, ".pdf"), height=5, width=10)
	## Calculating silhoute plotts
	x <- fviz_nbclust(mat_feature_corr, FUN = hcut, method = "silhouette", k.max=k.max_use)
	avg_sil_values <- x$data
	avg_sil_values$clusters <- as.numeric(avg_sil_values$clusters)
	
		
	if(colnames(mat_feature_corr) %like% "IGHV|TRAV|TRBV|TRGV|TRDV"){
		min_use <- 10
	} else { 
		min_use <- 10
	}
	max1 <- Rfast::nth(avg_sil_values$y[avg_sil_values$clusters >=min_use], 1, descending = T)
	max1 <- avg_sil_values$clusters[avg_sil_values$y==max1]
	x <- x + geom_vline(xintercept=max1, col="red")
	plot(x)
	# Other plots 	
	plot(fviz_nbclust(mat_feature_corr, FUN = hcut, method = "wss", k.max=k.max_use))
	plot(fviz_nbclust(mat_feature_corr, FUN = hcut, method = "gap_stat", k.max=k.max_use))
	dev.off()
	print(paste0("Optimal Number of Clusters ", max1))
	
	## Do the actual clustering 
	w=10
	pdf(paste0(outputdir, "Plots/CLUSTERING_isotyper_metrics_FINAL_METRICS__FINAL", type, "_", iso_type, ".pdf"), width=40, height=10)
	hc = hclust(dist1, method=method_ac_use)
	par(mfrow= c(1,1), mar = c(5,5,3,3))
	plot(hc, cex=0.8, hang = -1)
	##settimg limits for trimming tree
	if(type %like% "BCR"){
		limit <- 0.5
	} 
	if(type %like% "TCR"){
		limit <- 0.1
	} 
	if(type %like% "V"){
		limit <- 0.2
	} 
	#threshold = quantile(dist1,limit) ## define threshold for clustering features
	rect.hclust(hc, border = "red", k=max1)
	dev.off()
	
	################################
	# Identify features/clusters 
	# Get clusters of dendogram (e.g. group by type)
	ct = cutree(hc, k=max1)
	clusters = sort(unique(ct))
	list_cluster = NULL
	for(c in c(1:length(clusters))){
	  w = which(ct==clusters[c])
	  list_cluster = c(list_cluster, list(names(ct[w])))
	}
	# How many variables are in each module...
	## save list to txt file 
	df_cluster_assignement <- data.frame()
	for(i in 1:length(colnames(mat_feature_corr))){
		feature_find <- colnames(mat_feature_corr)[i]
		op <- lapply(list_cluster, function(ch) grep(paste0("\\b", feature_find,"\\b"),ch, perl=TRUE))
		op1 <- op > 0
		cluster_member <- which(op1=="TRUE")
		row_df <- c(feature_find, cluster_member)
		df_cluster_assignement <- rbind(row_df, df_cluster_assignement)
		}
	colnames(df_cluster_assignement) <- c("Feature", "Cluster")
	df_cluster_assignement$Datatype <- data_type
	write.table(df_cluster_assignement, paste0(outputdir, "Summary/Clustered_Features_", type, "_", iso_type, "_", data_type, ".txt"), row.names=FALSE)
	return(list_cluster)

}





