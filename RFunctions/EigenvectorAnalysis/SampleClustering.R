library(reshape2)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(stringr)
library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(data.table)
library(foreach)
library(doParallel) 
library(ggforce)
library(plot3D)
library(missForest)
library(mice)
library(VIM)
library(foreach)
library(doParallel)
library(umap)
library(psych)
library(corrplot)
library(lavaan)
library(factoextra)
library(cowplot)
library(ggpubr)
library(matrixStats)
library(RColorBrewer)
library(factoextra)
library(dqrng)
library(dynamicTreeCut)
library(moduleColor)
library(WGCNA)
library(fastcluster)
library(cluster) 

#new_eigenvectors <- "/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_BCR_PRODUCTIVE.txt"
# Lets use the imputed data  - see what is driving the seperation
#new_eigenvectors <- "/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_BCR_PRODUCTIVE.txt"
#outputdir <- "/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/"
#type_use <- "BCR"
#iso_type <- "PRODUCTIVE"
#meta <- '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/sepsis_meta_health.txt'
#minClusterSizex <- 20

cluster_samples <- function(new_eigenvectors, outputdir, type_use, iso_type, meta, minClusterSizex, thresh=1.4){
	#-----------------------------------------------------------------------------
   	# Can we cluster samples based on their module eigenvectors
	## First generate distance tree
	## It is common to scale before hand 
	
	new_eigenvectors <- read.delim(new_eigenvectors)
	
	print("Removing Samples where RNAseq suggests a sample mixup")
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	new_eigenvectors <- new_eigenvectors[!rownames(new_eigenvectors) %in% bad_ids,]
	
	### 1 Scale!!!
	new_eigenvectors_scaled <- scale(new_eigenvectors)
	res.dist <- dist(new_eigenvectors_scaled, method = "euclidean")
	
	## Now which cluster aglomeration method 
	m <- c( "average", "single", "complete", "ward", "weighted")
	names(m) <- c( "average", "single", "complete", "ward", "weighted")
	print(paste0("Determining which cluster aglomeration method to use"))
	# function to compute coefficient
	ac <- function(x) {
	  agnes(res.dist, method = x)$ac
	}
	# Best ac score is best agglomerative coefficient 
	# Amount of clustering found!
	ac_scores <- data.frame(map_dbl(m, ac))
	colnames(ac_scores) <- "agglomerative_coefficient"
	ac_scores$Method <- rownames(ac_scores)
	method_ac_use <- ac_scores$Method[ac_scores$agglomerative_coefficient ==max(ac_scores$agglomerative_coefficient)]
	##rename for using hclust
	if(method_ac_use=="ward"){
		method_ac_use="ward.D2"
	}
	print(paste0("Cluster Agglomeration Method: ",  method_ac_use))
	
	## Plotting cluster scores per method
	pdf(paste0(outputdir, "Plots/CLUSTER_AGLOMERATION_METHOD_SAMPLECLUSTERING_", type_use, "_", iso_type, ".pdf"), height=5, width=5)
	plot(ggplot(ac_scores, aes(x=Method, y=agglomerative_coefficient, fill=Method)) +geom_col() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Cluster Agglomeration Method")+ ylab("Agglomerative Coefficient"))
	dev.off()
	
	###############################################################
	## lets get optimal number of clusters
	## If you go higher than this it starts causing issues with NB clust - may need tweaking based on scenario
	print(paste0("Determining Optimal Number of Sample Clusters :"))
	## Now lets plot the dedrogram using euclidean distance and the best aglomeration method
	w=10
	## Lets make the dendrogram and plot partition
	hc = hclust(res.dist, method=method_ac_use)
	
	######################
	### Now lets try dynamic tree cut 
	print("Attempting Dynamic Tree Cutting: deep Split")
	pdf(paste0(outputdir, "Plots/CLUSTERING_PARTITION_Samples_DynamicTreeCut_", type_use, "_", iso_type, ".pdf"), width=40, height=10)
	
	## This is the smallest group (health n=12)
	dynamic_cut <-  cutreeDynamic(hc, minClusterSize = minClusterSizex, method="hybrid", deepSplit = FALSE, distM=as.matrix(res.dist), maxDistToLabel = 0, verbose = 0)
	dynamic_cut2 <- labels2colors(dynamic_cut)
	plotDendroAndColors(dendro=hc,colors=dynamic_cut2, main="Feature Dendrogram: Dynamic Tree Cut")
	dev.off()
	
	
	#########	
	no_detected <- length(unique(dynamic_cut[dynamic_cut !=0]))
	print(paste0("Dynamic Tree Cut Proposes ", 	no_detected, " Clusters"))
	no_detected <- length(unique(dynamic_cut[dynamic_cut !=0]))
	
	
	### Lets get assignments 
	features_x <- data.frame(rownames(as.matrix(res.dist)))
	features_x$cluster <- dynamic_cut
	features_xx <- features_x[order(features_x$cluster),]
	rownames(features_xx) <- features_xx[,1]
	features_xx[,1] <- NULL
	
	## If any is assigned a 0 it means it couldnt be assigned and we want to annotate that!	
	if(any(features_xx$cluster==0)){
		print("Unassigned features present")
		print("Will need to store seperately")
		unassigned <- data.frame(features_xx[features_xx$cluster==0,])
		rownames(unassigned) <- rownames(features_xx)[features_xx$cluster==0]
		print("Filtered unassigned features")
		## Save output so easy to check if we have any unassigned and will want to explore seperately!
		colnames(unassigned) <- c("Cluster")
		unassigned$Datatype <- "Imputed"
		unassigned$Sample <- rownames(unassigned)
		write.table(unassigned, paste0(outputdir, "Summary/UNASSIGNED_Samples_", type_use, "_", iso_type, ".txt"), row.names=FALSE)
		
	}
	colnames(features_xx) <- "Cluster"
    features_xx$Cluster[features_xx$Cluster ==0] <- "X"
	features_xx$Cluster <- as.factor(features_xx$Cluster)
	
	## Sample assignments 
	##assignment <- data.frame(e$Best.partition)
	##colnames(assignment) <- c("Cluster")
	##assignment$Cluster <- as.factor(assignment$Cluster)
	
	## Save the assignments 
	#assignment_write <- assignment
	assignment <- features_xx
	assignment_write <- features_xx
	assignment_write$Sample <- rownames(assignment_write)
	write.table(assignment_write, paste0(outputdir, "Summary/Clustered_Samples_assignment_", type_use, "_", iso_type, ".txt"), row.names=FALSE)
	
    ########################################
	## now lets merge with pcs 
	
	print("Performing Sample Clustering on Module Scores")
	mtcars.pca <- prcomp(new_eigenvectors_scaled, center = TRUE,scale. = TRUE)
	pcs <- data.frame(mtcars.pca$x)
	pcs$DAY <- NA
	pcs$DAY[rownames(pcs) %like% "_1"] <- "Day1"
	pcs$DAY[rownames(pcs) %like% "_3"] <- "Day3"
	pcs$DAY[rownames(pcs) %like% "_5"] <- "Day5"
	rownames(pcs) <- gsub("_productive", "", rownames(pcs))
	rownames(assignment) <- gsub("_productive", "", rownames(assignment))
	loadingsx <- data.frame(mtcars.pca$rotation)
	loadingsx$feature <- rownames(loadingsx)
	meta_datause <- read.delim(meta, sep='\t', header=TRUE)
	rownames(meta_datause) <- meta_datause$SampleID_alternative
	pcs <- merge(pcs, meta_datause, by = 0, all.x=TRUE)
	rownames(pcs) <- pcs$Row.names
	pcs$Row.names <- NULL
	
	pcs <- merge(pcs, assignment, by=0)
	pcs$Mortality2[is.na(pcs$Mortality2)] <- "TECHNICAL"
	pcs$DAY[pcs$Mortality2=="TECHNICAL"] <- "TECHNICAL"
	pcs$DISEASE[pcs$Mortality2=="TECHNICAL"] <- "TECHNICAL"
	pc1_sd <- sd(loadingsx$PC1)
	pc1_mean<- mean(loadingsx$PC1)
	up <- pc1_mean+(pc1_sd*thresh)
	down <- pc1_mean-(pc1_sd*thresh)
	pc2_sd <- sd(loadingsx$PC2)
	pc2_mean<- mean(loadingsx$PC2)
	up1 <- pc2_mean+(pc2_sd*thresh)
	down1 <- pc2_mean-(pc2_sd*thresh)
	loadingsx$feature <- gsub("BCR_READS_", " ",loadingsx$feature)
	loadingsx$feature <- gsub("__", " ",loadingsx$feature)
	loadingsx$feature <- gsub("_", " ",loadingsx$feature)
	loadingsx$feature <- str_wrap(loadingsx$feature, width = 20)
	
	
	variance.explained <-t(data.frame(summary(mtcars.pca)$importance))
	variance.explained <- data.frame(variance.explained)
	variance.explained$PC <- rownames(variance.explained)
	variance.explained$PC <- gsub("PC", "", variance.explained$PC)
	variance.explained$PC <- as.numeric(variance.explained$PC)
	
	## Going to remove the technical for this plot 
	pc2 <- pcs
	pc2 <- pc2[pc2$DISEASE != "TECHNICAL" & pc2$Cluster !="X",]
	props <- pc2 %>%
	group_by(Cluster, DAY, Mortality2) %>%
	  summarise(n = n()) %>%
	  mutate(freq = n / sum(n))
	props <- data.frame(props)
  
	 pie_chart <- ggplot(props, aes(x = "", y = freq, fill = Mortality2)) +
	  geom_bar(width = 1, stat = "identity") +
	  coord_polar(theta = "y") + theme_classic()+theme(axis.text = element_blank(),
	        axis.ticks = element_blank())+
	  facet_grid(cols=vars(Cluster), rows=vars(DAY))+ylab("Cluster")+xlab("Frequency") +ggtitle(paste0("Min Cluster Size ", minClusterSizex))+labs(fill="Mortality")
  
  
	## Now lets plot 
	pdf(paste0(outputdir, "Plots/PCAEigenvectors_",  type_use, "_", iso_type,".pdf"), height=11, width=10)
	#p <- ggplot(pcs, aes(x = PC1, y = PC2, shape=Mortality2, colour=Cluster)) +   geom_point()+theme_classic()+labs(shape="Mortality", colour="Cluster")
	p1 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=Cluster)) +   geom_point()+theme_classic() + facet_wrap(~DAY, ncol=4)+labs(colour="Cluster")
	p2 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=DISEASE, shape=Mortality2)) +   geom_point()+theme_classic() + facet_wrap(~DAY, ncol=4)+labs(colour="DISEASE",shape="Mortality")
	p3 <- ggplot(loadingsx, aes(x = PC1, y = PC2)) +   geom_point(colour="red")+theme_classic()+ gghighlight((PC1> up |  PC1<down)& (PC2> up1 |  PC2<down1), label_key = feature)+xlab("PC1 Loadings") + ylab("PC2 Loadings") +ggtitle(paste0("Mean +/- 1SD*", thresh))
	p4 <- ggplot(variance.explained, aes(x = PC, y = Cumulative.Proportion)) +   geom_point(colour="red")+theme_classic()+ gghighlight(Proportion.of.Variance>0.05)+ylab("Cumulative Proportion of Variance") +xlab("Principal Component")
	#p5 <- ggplot(props, aes(x = Cluster, y = freq)) +   geom_col()+theme_classic()+facet_grid(rows=vars(Mortality2), cols=vars(DAY))

	plot(plot_grid(p2, p1, pie_chart, ncol=1, rel_heights=c(1, 1,2), align = "hv", axis="tblr", labels="AUTO"))
	plot(plot_grid(p3, p4, ncol=1, rel_heights=c(1, 1,1), align = "hv", axis="tblr", labels="AUTO"))
	dev.off()	
	
	write.table(pcs, paste0(outputdir, "Summary/SampleClustering_PCS_assignment_", type_use, "_", iso_type, ".txt"), row.names=FALSE)
	write.table(loadingsx, paste0(outputdir, "Summary/SampleClustering_PCS_LOADINGS_", type_use, "_", iso_type, ".txt"), row.names=FALSE)
	write.table(variance.explained, paste0(outputdir, "Summary/SampleClustering_PCS_VARIANCE_", type_use, "_", iso_type, ".txt"), row.names=FALSE)


	print("Finished Module Reduction Pipeline")
}