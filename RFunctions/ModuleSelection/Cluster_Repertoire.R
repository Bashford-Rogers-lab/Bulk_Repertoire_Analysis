## Function to Cluster filtered matrix 
## Lauren Overend and Rachael Bashford-Rogers
## lauren.overend@oriel.ox.ac.uk
## Jan 2022 
# used tutorial from https://uc-r.github.io/kmeans_clustering#silo and https://uc-r.github.io/hc_clustering
library(dendextend)
library(cluster)
library(Rfast)
options(bitmapType='cairo-png')
library(NbClust)
library(dynamicTreeCut)
library(moduleColor)
library(WGCNA)
library(fastcluster)
library(gridExtra)
library(grid)


cluster_features <- function(outputdir, correlation_matrix, type, iso_type, data_type){
	mat_feature_corr <- correlation_matrix
	
	w=10
	pdf(paste0(outputdir, "Plots/CLUSTERING_CorrelatedvsAnticorrelated_dendrogram_", type, "_", iso_type, ".pdf"), height=w*1*2, width=w*4.5*1.5)
	par(mfrow= c(2,1), mar = c(5,5,3,3))
	#### build distance tree
	## option 1: anticorrelated metrics will be more distant than correlated
	#https://bioinformatics.mdanderson.org/Software/OOMPA/ClassDiscovery/html/distanceMatrix.html#:~:text=The%20most%20common%20metric%20used,%2Dcor(dataset))%2F2%20.&text=The%20spearman%20metric%20used%20the,correlation%20for%20the%20Pearson%20correlation.
	dist1 =as.dist(1-mat_feature_corr)/2
	plot(hclust(dist1), main = "option 1: correlation", hang = -1)
	## option 2: anticorrelated metrics will be same as correlated (i.e. distance based on abs(correlation))
	dist2 = as.dist(1-abs(mat_feature_corr))
	plot(hclust(dist2), main = "option 2: absolute correlation", hang = -1)
	dev.off()
	
	## Cluster features
	## Using anti-correlated tree 
	## this is the one we will take forward 
	print("We will use correlation not absolute correlation")
	dist1 = as.dist((1-mat_feature_corr)/2)
	## Plot a map of distance between measures
	heights=dim(mat_feature_corr)[1]/4
	pdf(paste0(outputdir, "Plots/Heamap_of_distance_", type, "_", iso_type, ".pdf"), height=10, width=10)
	plot(fviz_dist(dist1, show_labels = FALSE))
	dev.off()
 
 
	################################
	## What aglomeration method to use for clustering 
	## default of hclus is complete 
	m <- c( "average", "single", "complete", "ward", "weighted")
	names(m) <- c( "average", "single", "complete", "ward", "weighted")
	print(paste0("Determining which cluster aglomeration method to use"))
	# function to compute coefficient
	## agnes will do the heirachical clustering if we provide a df, but as we provide distance matrix it just shows best method
	ac <- function(x) {
	  agnes(dist1, method = x)$ac
	}
	# Best ac score is best agglomerative coefficient 
	# Indicates amount of clustering structure found!
	ac_scores <- data.frame(map_dbl(m, ac))
	colnames(ac_scores) <- "agglomerative_coefficient"
	ac_scores$Method <- rownames(ac_scores)
	## Plotting cluster scores per method
	pdf(paste0(outputdir, "Plots/CLUSTER_AGLOMERATION_METHOD_", type, "_", iso_type, ".pdf"), height=5, width=5)
	plot(ggplot(ac_scores, aes(x=Method, y=agglomerative_coefficient, fill=Method)) +geom_col() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Cluster Agglomeration Method")+ ylab("Agglomerative Coefficient"))
	dev.off()
	## Select Method
	method_ac_use <- ac_scores$Method[ac_scores$agglomerative_coefficient ==max(ac_scores$agglomerative_coefficient)]
	print("Maximum AC score")
	print(max(ac_scores$agglomerative_coefficient))
	
	##rename for using hclust
	if(method_ac_use=="ward"){
		method_ac_use="ward.D2"
	}
	print(paste0("Hericical Cluster Method: ",  method_ac_use))
	################################
	
	################################
	print(paste0("Determining Optimal Number of clusters with NB clust:"))
	if(type %like% "TCRAB") {
		print("Type Receptor: TCRAB")
		top_clust <- 18
		bot_clust <- 10
		min_size <- 5
		dimension <- 40
	} 
	if(type %like% "TCRGD") {
		print("Type Receptor: TCRGD")
		top_clust <- 18
		bot_clust <- 8
		min_size <- 3
		dimension <- 15
	} 
	if(type %like% "BCR"){
		print("Type Receptor: BCR")
		top_clust <- 25
		bot_clust <- 10
		min_size <- 10
		dimension <- 70
	} 
	
	#################################################
	## Lets run a load of iteractions
	cluster_iterate <- c()
	print("Running Fixed Height Cutting Iteration min cluster 2")
	for(x in c((2+2):40)){
		if(x > dim(mat_feature_corr)[1]){
			break 
		}
		resd <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = 2, max.nc = x, method = "ward.D2", index = 'dunn'))
		ress <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = 2, max.nc = x, method = "ward.D2", index = 'silhouette'))
		resm <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = 2, max.nc = x, method = "ward.D2", index = 'mcclain'))
		resc <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = 2, max.nc = x, method = "ward.D2", index = 'cindex'))
		resf <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = 2, max.nc = x, method = "ward.D2", index = 'frey'))
		result <- rbind(resd$Best.nc, ress$Best.nc, resm$Best.nc, resc$Best.nc, resf$Best.nc)
		rownames(result) <- c("Dunn", "Silhouette", "Mcclain", "Cindex", "Frey")
		result <- data.frame(result)
		maj <- data.frame(table(result$Number_clusters))
		colnames(maj) <- c("NoClusters", "CountMethods")
		number_clusters <- maj[maj$CountMethods==max(maj$CountMethods),]
		###Get the maximum number of clusters 
		if(dim(number_clusters)[1]==1){
			number_k <- as.numeric(as.character(number_clusters$NoClusters))
		} else {
			number_clusters$NoClusters <- as.numeric(as.character(number_clusters$NoClusters))
			number_k <- number_clusters$NoClusters[number_clusters$NoClusters==max(number_clusters$NoClusters)]
		}
		row_clust <- c(x, number_k)
		cluster_iterate <- rbind(cluster_iterate,row_clust)
	}		
	cluster_iterate <- data.frame(cluster_iterate)
	cluster_iterate[,1] <- as.numeric(cluster_iterate[,1])
	cluster_iterate[,2] <- as.numeric(cluster_iterate[,2])
	pdf(paste0(outputdir, "Plots/SELECT_OPTIMAL_CLUSTERS_Iteraction_Min2_", type, "_", iso_type, ".pdf"), height=10, width=10)
	p_it <- ggplot(cluster_iterate, aes(x=X1, y=X2)) + geom_point() +geom_line() + theme_bw()+xlab("User Specified Max Number of Clusters") +ylab("Optimal Cluster Number")+labs(colour="none")+ggtitle("Fixed Height Cutting: Min Cluster 2")
	plot(p_it)
	dev.off()
	
	###################################################################
	cluster_iterate <- c()
	print(paste0("Running Fixed Height Cutting Iteration min cluster ", min_size))
	for(x in c((min_size+2):40)){
		if(x > dim(mat_feature_corr)[1]){
			break 
		}
		resd <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = min_size, max.nc = x, method = "ward.D2", index = 'dunn'))
		ress <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = min_size, max.nc = x, method = "ward.D2", index = 'silhouette'))
		resm <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = min_size, max.nc = x, method = "ward.D2", index = 'mcclain'))
		resc <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = min_size, max.nc = x, method = "ward.D2", index = 'cindex'))
		resf <- suppressMessages(NbClust(diss=dist1, distance = NULL, min.nc = min_size, max.nc = x, method = "ward.D2", index = 'frey'))
		result <- rbind(resd$Best.nc, ress$Best.nc, resm$Best.nc, resc$Best.nc, resf$Best.nc)
		rownames(result) <- c("Dunn", "Silhouette", "Mcclain", "Cindex", "Frey")
		result <- data.frame(result)
		maj <- data.frame(table(result$Number_clusters))
		colnames(maj) <- c("NoClusters", "CountMethods")
		number_clusters <- maj[maj$CountMethods==max(maj$CountMethods),]
		###Get the maximum number of clusters 
		if(dim(number_clusters)[1]==1){
			number_k <- as.numeric(as.character(number_clusters$NoClusters))
		} else {
			number_clusters$NoClusters <- as.numeric(as.character(number_clusters$NoClusters))
			number_k <- number_clusters$NoClusters[number_clusters$NoClusters==max(number_clusters$NoClusters)]
		}
		row_clust <- c(x, number_k)
		cluster_iterate <- rbind(cluster_iterate,row_clust)
	}		
	cluster_iterate <- data.frame(cluster_iterate)
	cluster_iterate[,1] <- as.numeric(cluster_iterate[,1])
	cluster_iterate[,2] <- as.numeric(cluster_iterate[,2])
	pdf(paste0(outputdir, "Plots/SELECT_OPTIMAL_CLUSTERS_Iteraction_", type, "_", iso_type, ".pdf"), height=10, width=10)
	p_it2 <- ggplot(cluster_iterate, aes(x=X1, y=X2)) + geom_point() +geom_line() + theme_bw()+xlab("User Specified Max Number of Clusters") +ylab("Optimal Cluster Number")+labs(colour="none")+ggtitle(paste0("Fixed Height Cutting: Min Cluster ", min_size))
	plot(p_it)
	dev.off()
	
	##########################################################
	### Run properly 	
	resd <- NbClust(diss=dist1, distance = NULL, min.nc = bot_clust, max.nc = top_clust, method = "ward.D2", index = 'dunn')
	ress <- NbClust(diss=dist1, distance = NULL, min.nc = bot_clust, max.nc = top_clust, method = "ward.D2", index = 'silhouette')
	resm <- NbClust(diss=dist1, distance = NULL, min.nc = bot_clust, max.nc = top_clust, method = "ward.D2", index = 'mcclain')
	resc <- NbClust(diss=dist1, distance = NULL, min.nc = bot_clust, max.nc = top_clust, method = "ward.D2", index = 'cindex')
	resf <- NbClust(diss=dist1, distance = NULL, min.nc = bot_clust, max.nc = top_clust, method = "ward.D2", index = 'frey')
	## Number of Clusters
	result <- rbind(resd$Best.nc, ress$Best.nc, resm$Best.nc, resc$Best.nc, resf$Best.nc)
	rownames(result) <- c("Dunn", "Silhouette", "Mcclain", "Cindex", "Frey")
	result <- data.frame(result)
	## Score change
	result2 <- rbind(resd$All.index, ress$All.index, resm$All.index, resc$All.index, resf$All.index)
	rownames(result2) <- c("Dunn", "Silhouette", "Mcclain", "Cindex", "Frey")
	result2 <- data.frame(result2)
	result2$Method <- rownames(result2)
	result2 <- reshape2::melt(result2, id.vars = c("Method"))
	colnames(result2) <- c("Method", "Number_of_Clusters", "Score")
	result2$Number_of_Clusters <- as.numeric(gsub("X", "", result2$Number_of_Clusters))
	
	## Plot these
	pdf(paste0(outputdir, "Plots/SELECT_OPTIMAL_CLUSTERS_", type, "_", iso_type, ".pdf"), height=8, width=10)
	pall <- ggplot(result2, aes(x=Number_of_Clusters, y=Score, colour=Method)) + geom_point() +geom_line() + theme_bw() +facet_wrap(~Method, scales="free")+xlab("Number of Clusters") +ylab("Method Specific Score")+guides(colour="none")+ggtitle("Fixed Height Cutting")
	plot(pall)
	plot(ggarrange(ggarrange(p_it, p_it2, nrow=2), pall, ncol=2, labels="AUTO"))
	dev.off()
	
	## Majority rule
	maj <- data.frame(table(result$Number_clusters))
	colnames(maj) <- c("NoClusters", "CountMethods")
	number_clusters <- maj[maj$CountMethods==max(maj$CountMethods),]
	
	if(dim(number_clusters)[1]==1){
		number_k <- as.numeric(as.character(number_clusters$NoClusters))
		print("Majoriy Rules")
		print(paste0("Optimum Number of clusters is ", number_k))
	} else {
		number_clusters$NoClusters <- as.numeric(as.character(number_clusters$NoClusters))
		number_k <- number_clusters$NoClusters[number_clusters$NoClusters==max(number_clusters$NoClusters)]
		print("Tied Optimal Number of Clusters - taking largest number")
		print(paste0("Optimum Number of clusters is ", number_k))
	}
	
	#########################################
	## Partition the tree and plot 
	w=10
	pdf(paste0(outputdir, "Plots/CLUSTERING_PARTITION_", type, "_", iso_type, ".pdf"), width=dimension, height=10)
	hc = hclust(dist1, method=method_ac_use)
	par(mfrow= c(1,1), mar = c(5,5,3,3))
	plot(hc, cex=0.8, hang = -1)
	rect.hclust(hc, border = "red", k=number_k)
	dev.off()
	
	#########################################
	### Lets try using dynamic tree cut allows for cutting at different heights of the dendrogram
	print("Attempting Dynamic Tree Cutting: deep Split")
	pdf(paste0(outputdir, "Plots/CLUSTERING_PARTITION_DynamicTreeCut_", type, "_", iso_type, ".pdf"), width=dimension, height=10)
	dynamic_cut <-  cutreeDynamic(hc, minClusterSize = min_size, method="hybrid", deepSplit = TRUE, distM=as.matrix(dist1), maxDistToLabel = 0, verbose = 0)
	dynamic_cut2 <- labels2colors(dynamic_cut)
	plotDendroAndColors(dendro=hc,colors=dynamic_cut2, main="Feature Dendrogram: Dynamic Tree Cut", cex.rowText=0.3)
	dev.off()
	
	pdf(paste0(outputdir, "Plots/CLUSTERING_PARTITION_DynamicTreeCut_NoLabels_", type, "_", iso_type, ".pdf"), width=15, height=7.5)
	plotDendroAndColors(dendro=hc,colors=dynamic_cut2, main=paste0("Feature Dendrogram: Dynamic Tree Cut ", type, " Repertoire"), cex.rowText=0.3, dendroLabels=FALSE, addGuide=TRUE)
	dev.off()
	
	## How many clusters have we assigned
	no_detected <- length(unique(dynamic_cut[dynamic_cut !=0]))
	print(paste0("Dynamic Tree Cut Proposes ", 	no_detected, " Clusters"))
	no_detected <- length(unique(dynamic_cut[dynamic_cut !=0]))
	if(no_detected < number_k){
		print("This is less than the number proposed by fixed height branch cut")
		print("Proceeding with Dynamic Tree Cut")
	} else if (no_detected== number_k) {
		print("This is the same as the number proposed by fixed height branch cut")
		print("Proceeding with Dynamic Tree Cut")
	} else {
		print("This is more than the number proposed by fixed height branch cut")
		print("Proceeding with Dynamic Tree Cut")
	}
	
	## Lets get the clusters 
	features_x <- data.frame(rownames(as.matrix(dist1)))
	features_x$cluster <- dynamic_cut
	features_xx <- features_x[order(features_x$cluster),]
	colnames(features_xx) <- c("feature", "cluster")
    ## For saving later 
	df_cluster_assignement <- features_xx
	df_cluster_assignement$cluster[df_cluster_assignement$cluster==0] <- "Unassigned"

	## Unassigned objects are labelled 0 by dynamic tree cut 
	## We will want to remove these so they are not put into module reduction as one level!
	if(any(features_xx$cluster==0)){
		print("WARNING!")
		print("Unassigned features present!!!!!")
		print("Will need to store seperately")
		unassigned <- features_xx[features_xx$cluster==0,]
		## only want assigned clusters for PCA otherwise this will cause confusion
		features_xx <- features_xx[features_xx$cluster !=0,]
		print("Filtered unassigned features")
		## Save output so easy to check if we have any unassigned and will want to explore seperately!
		colnames(unassigned) <- c("Feature", "Cluster")
		unassigned$Datatype <- data_type
		write.table(unassigned, paste0(outputdir, "Summary/UNASSIGNED_FEATURES_", type, "_", iso_type, "_", data_type, ".txt"), row.names=FALSE)
	}
	## Lets turn it into a sorted list for PCA 
	clusters = sort(unique(features_xx$cluster))
	list_cluster = NULL
	for(c in c(1:length(clusters))){
	  w = features_xx$feature[features_xx$cluster==clusters[c]]
	  list_cluster = c(list_cluster, list(w))
	}
	
	
	## Old code for running with HC clust static clustering
	################################
	# Identify features/clusters 
	# Get clusters of dendogram (e.g. group by type)
	#ct = cutree(hc, k=number_k)
	#clusters = sort(unique(ct))
	#list_cluster = NULL
	#for(c in c(1:length(clusters))){
	 # w = which(ct==clusters[c])
	 # list_cluster = c(list_cluster, list(names(ct[w])))
	#}
	# How many variables are in each module...
	## save list to txt file 
	#df_cluster_assignement <- data.frame()
	#for(i in 1:length(colnames(mat_feature_corr))){
		#feature_find <- colnames(mat_feature_corr)[i]
		#op <- lapply(list_cluster, function(ch) grep(paste0("\\b", feature_find,"\\b"),ch, perl=TRUE))
		#op1 <- op > 0
		#cluster_member <- which(op1=="TRUE")
		#row_df <- c(feature_find, cluster_member)
		#df_cluster_assignement <- rbind(row_df, df_cluster_assignement)
		#}
	
	colnames(df_cluster_assignement) <- c("Feature", "Cluster")
	#### Save and write to file!!!
	df_cluster_assignement$Datatype <- data_type
	write.table(df_cluster_assignement, paste0(outputdir, "Summary/Clustered_Features_assignment_", type, "_", iso_type, "_", data_type, ".txt"), row.names=FALSE)
	saveRDS(list_cluster, paste0(outputdir, "Summary/ClusterLIST.rds"))

	print("Done Feature Clustering")
	return(list_cluster)
}





