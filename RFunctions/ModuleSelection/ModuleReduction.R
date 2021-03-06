## Code for preforming module reduction and clustering of isotyper metrics
## Based on WGCNA type analysis 
## Lauren Overend lauren.overend@oriel.ox.ac.uk
## Rachael Bashford Rogers
## Welcome Trust Centre for Human Genetics 
## 07/01/2022 

## Packages 
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

cluster_immunerep <- function(type_use, subsampled_deptha,subsampled_depthb=NA, productivity, outputdir1, outputdir2=NA, plot_outputdir, meta_data){
## Functions for Script 
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/ModuleImputationMethods.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/CorrelationMatrix.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/Get_Subsample_Depth.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/Cluster_Repertoire.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/Eigengenes.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/PlotModules.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/Other.R')

	## Variable Names!!!
	iso_type <- productivity
	type_use <- type_use
	outputdir <- plot_outputdir
	meta_data <- meta_data
	#################################
	## Reading in the Final Matrix 
	if(outputdir1 %like% "BCR"){
		mat_filtered <- read.delim(paste0(outputdir1, "Summary/isotyper_metrics_filtered_FINAL_METRICS_", subsampled_deptha, "_", iso_type, ".txt"), sep="\t")
	}
	## Reading in the Final Matrix 
	if(outputdir1 %like% "TCR"){
		mat_filtered1 <- read.delim(paste0(outputdir1, "Summary/isotyper_metrics_filtered_FINAL_METRICS_", subsampled_deptha , "_", iso_type, ".txt"), sep="\t")
		mat_filtered1$sample <- rownames(mat_filtered1)
		mat_filtered2 <- read.delim(paste0(outputdir2, "Summary/isotyper_metrics_filtered_FINAL_METRICS_", subsampled_depthb, "_", iso_type, ".txt"), sep="\t")
		mat_filtered2$sample <- rownames(mat_filtered2)
		mat_filtered <- merge(mat_filtered1, mat_filtered2, by="sample", all=TRUE)
		rownames(mat_filtered) <- mat_filtered$sample
		mat_filtered$sample <- NULL
	}
	
	colnames(mat_filtered) <- gsub("_PRODUCTIVE", "", colnames(mat_filtered))
	colnames(mat_filtered) <- gsub("_UNPRODUCTIVE", "", colnames(mat_filtered))
	
	## We will use this later
	mat_filtered_vgenes <- mat_filtered[,colnames(mat_filtered) %in% grep("IGHV|TRAV|TRBV|TRGV|TRDV", colnames(mat_filtered), value=TRUE)]
	## Split into Features and VDJ usage 
	mat_filtered <- mat_filtered[,!colnames(mat_filtered) %in% grep("IGHV|TRAV|TRBV|TRGV|TRDV", colnames(mat_filtered), value=TRUE)]
	print("Split into Features and V gene usage dataframe")
	

	############################
	## Compare Imputation Methods for selected features 
	## Missing-ness threshold was set in the Summary Isotyper as >60% 
	method_try="NO"
	p <- mat_filtered
	features <- colnames(mat_filtered)
			
	if(method_try=="YES"){
			
			imputation_methods <- c()
			if(type_use %like% "BCR"){
				for(i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)){ #, 0.3, 0.4, 0.5, 0.6
					u <- optimise_imputing(mat_filtered, iso_type, i, outputdir, type_use)
					imputation_methods <- rbind(imputation_methods,u)
				} 
			} else {
			  for(i in c(0.1, 0.2, 0.3, 0.4)){ #, 0.3, 0.4, 0.5, 0.6
					u <- optimise_imputing(mat_filtered, iso_type, i, outputdir, type_use)
					imputation_methods <- rbind(imputation_methods,u)
				}		
			}
			
			method_means <- imputation_methods %>% dplyr::group_by(method, NAseedfrequency) %>% dplyr::summarise(mean_error = mean(rsme))
			method_means <- data.frame(method_means)
			method_use <- method_means$method[method_means$mean_error==min(method_means$mean_error)]
			print(paste0("Best imputation method: ", method_use))
			nrmse_values <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
			nrmse_files <- grep("nrmse", nrmse_values, value=TRUE)
			nrmse_files <- grep(type_use, nrmse_files, value=TRUE)
			nrmse_list <- lapply(nrmse_files, fread, header = TRUE, sep="\t", fill=TRUE)
			nrmse_df <- rbindlist(nrmse_list, idcol = "na_freq")	
			nrmse_df2<- reshape2::melt(nrmse_df, id.vars=c("na_freq"))
			nrmse_df2$na_freq <- nrmse_df2$na_freq/10
			nrmse_df2$variable <- gsub("nrmse_", "", nrmse_df2$variable)
			colnames(nrmse_df2) <- c("NAseedfrequency", "Method", "NRMSE")

			## Plot the result of different levels of missingness per method: 
			pdf(paste0(outputdir, "Plots/Imputation1_Summary_", type_use, "_", iso_type, ".pdf"), width=13, height=13)
			plot(ggplot(imputation_methods, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_violin(alpha=0.5) + stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red") +ggtitle("Wilcox test comparing RSME against base-mean") +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE"))
			plot(ggplot(imputation_methods, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_boxplot(alpha=0.5)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red")+ggtitle("Wilcox test comparing RSME against base-mean") +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE"))
			plot(ggplot(imputation_methods, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_boxplot(alpha=0.5)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red")+ggtitle("Wilcox test comparing RSME against base-mean") +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE"))
			plot(ggplot(nrmse_df2, aes(x=Method, y=NRMSE, fill=Method)) +geom_col() +theme_bw() +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("NRMSE"))
			dev.off()
	}

	################################
	# GENERATE CORRELATION MATRIX!!!!:
	## We do this only on MEASURED VALUES!
	# First we must identify depth to subsample 
	features <- colnames(mat_filtered)
	min_sample_input <- get_subsample_depth(mat_filtered, features)
	## Calculate Correlation Between different Measures 
	## Output will be saved to the Summary Directory
	mat_feature_corr <- Get_subsample_corr_matrix(features, p, min_sample_input, outputdir, iso_type, type_use)
	## Cluster the correlation Matrix and Plot 
	list_cluster <- cluster_features(outputdir, mat_feature_corr, type_use, iso_type, "IMPUTED")
	# How many variables are in each module...
	size_cluster <- lengths(list_cluster)

	## Original data(used_later) (we have to reasign it here becuase we add categorical columns which arent for imputation)
	imputed_df <- mat_filtered
	imputed_df$Type <- "Original"
	imputed_df$sample <- rownames(imputed_df)
	imputed_df$Imputed <- "NO"
	imputed_df <- reshape2::melt(imputed_df, id.vars=c("sample", "Imputed", "Type"))


	###########################################################
	## MISSING DATA IMPUTATION for final Matrix 
	## Using the missForest non-parametric missing value imputation for mixed-type data
	###################
	imputed_Data <- missForest(mat_filtered, maxiter=20)
	Error <- (imputed_Data$OOBerror)
	print(paste0("NRMSE for missForest Imputed Data: ", Error))
	completeData = data.frame(imputed_Data$ximp)
	completeData$Type <- "Imputed"
	completeData$sample <- rownames(completeData)
	completeData$Imputed <- "NO"
	completeData <- reshape2::melt(completeData, id.vars=c("sample", "Imputed", "Type"))
	#Bind orignal and imputed data for visualising 
	final_data <- rbind(imputed_df, completeData)
	final_data$variable <- as.character(final_data$variable)

	###########################################################
	## Plotting the Imputed values for missForest
	imputed_samples <- final_data[is.na(final_data$value),]
	final_data$Type <- factor(final_data$Type, levels=c("Original", "Imputed"))
	for(i in 1:length(imputed_samples$sample)){
		#print(i)
		final_data$Imputed[final_data$sample==imputed_samples$sample[i] & final_data$Type=="Imputed" & final_data$variable==imputed_samples$variable[i]] <- "YES"
	} 
	plot_list <- list()
	for(i in unique(final_data$variable)){ 
		#print(i)
		datax <- final_data[final_data$variable==i,]
		datax$value <- as.numeric(datax$value)
		if(length(unique(datax$Imputed))>1){
			x <- ggplot(data = datax, aes_string(x="Type", y="value")) +geom_boxplot(aes(group = Type)) +geom_point(aes(color = as.factor(Imputed)), alpha = 0.5, position = "jitter", set.seed(1)) +scale_color_manual(values = c("red", "blue"),  drop=FALSE) +ylab(i) +theme_bw() +xlab("Imputation Type") + guides(color=guide_legend(title="Imputed Value")) + ggtitle("missForest Imputation")  + stat_compare_means(aes(label = ..p.signif..), col="darkgreen", label.x = 1.5)
			plot_list[[i]] <- x 
		}
	}
	pdf(paste0(outputdir, "Plots/ImputedValues_",  type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=20, width=20)
	for(i in 1:(round(length(plot_list)/16)+1)){
		stop_val <- i*16
		start_val <- stop_val-15
		#print(start_val)
		#print(stop_val)
		if(start_val > length(plot_list)){
			break
		}
		new_list <- plot_list[c(start_val:stop_val)]
		multiplot(plotlist = new_list, cols = 4)
	}
	dev.off()

	###########################################################
	## Calculate Eigen Vectors!!!!!
	full_data <- data.frame(imputed_Data$ximp)
	## scale data
	full_data <- scale(full_data) 
	## Get and plot eigen genes
	mat_eigenvectors <- get_eigengenes(full_data, list_cluster, outputdir,  type_use, subsampled_deptha, iso_type) 
	### Save our mat eigen vectors for plotting!!!
	#write.table(mat_eigenvectors, paste0(outputdir, "Eigenvectors_", type_use, "_", iso_type, ".txt"), sep='\t')
	
	## Want to repeat on the V genes MATRIX for which we dont need to do imputation!!
	print("Performing Module Analysis on V genes")
	
	##########################################################
	################################
	## Rerruning for V genes!!
	p <- mat_filtered_vgenes
	features <- colnames(mat_filtered_vgenes)
	
	# First we must identify depth to subsample 
	min_sample_input <- get_subsample_depth(mat_filtered_vgenes, features)
	## Calculate Correlation Between different Measures 
	## Output will be saved to the Summary Directory
	mat_feature_corr <- Get_subsample_corr_matrix(features, p, min_sample_input, outputdir, iso_type, paste0(type_use, "_Vgene"))
	## Cluster the correlation Matrix and Plot 
	list_cluster <- cluster_features(outputdir, mat_feature_corr, paste0(type_use, "_Vgene"), iso_type, "IMPUTED")
	# How many variables are in each module...
	size_cluster <- lengths(list_cluster)
	
	full_datav <- data.frame(mat_filtered_vgenes)
	## scale data
	full_datav <- scale(full_datav) 
	## Get and plot eigen genes
	mat_eigenvectors_vgenes <- get_eigengenes(full_datav, list_cluster, outputdir,  paste0(type_use, "_Vgene"), subsampled_deptha, iso_type) 
	colnames(mat_eigenvectors_vgenes) <- paste0("Vgene_", colnames(mat_eigenvectors_vgenes))
	
	## Merge Feature modules with V gene Modules. 
	print("Done Module Analysis on V genes")
	
	mat_eigenvectors_features <- mat_eigenvectors
	new_eigenvectors <- merge(mat_eigenvectors, mat_eigenvectors_vgenes, by=0)
	rownames(new_eigenvectors) <- new_eigenvectors$Row.names
	new_eigenvectors$Row.names <- NULL
	print("Calculated Eigenvectors and saved")
	write.table(new_eigenvectors, paste0(outputdir, "Eigenvectors_", type_use, "_", iso_type, ".txt"), sep='\t')
	
	## Plot 
	print("Plotting Modules")
	plot_eigenvectors(new_eigenvectors, meta_data, iso_type, outputdir, type_use, subsampled_deptha)
	print("Done and Saved Module Reduction")
	
	########################################
	# Can we cluster samples based on their module eigenvectors
	res.dist <- dist(new_eigenvectors, method = "euclidean")
	
	m <- c( "average", "single", "complete", "ward", "weighted")
	names(m) <- c( "average", "single", "complete", "ward", "weighted")
	print(paste0("Determining which cluster aglomeration method to use"))
	# function to compute coefficient
	ac <- function(x) {
	  agnes(res.dist, method = x)$ac
	}
	# Best ac score is best agglomerative coefficient 
	# Amoung of clustering found!
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
	pdf(paste0(outputdir, "Plots/CLUSTERING_isotyper_metrics_SAMPLES__AC_score", type_use, "_", iso_type, ".pdf"), height=5, width=5)
	plot(ggplot(ac_scores, aes(x=Method, y=agglomerative_coefficient, fill=Method)) +geom_col() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Cluster Agglomeration Method")+ ylab("Agglomerative Coefficient"))
	dev.off()
	
	print(paste0("Determining Optimal Number of clusters:"))
	## manual silhouett
	avg_sil <- function(k) {
		  km.res <- kmeans(new_eigenvectors, centers = k, nstart = 20)
		  ss <- silhouette(km.res$cluster, dist(new_eigenvectors))
		  mean(ss[, 3])
	}
	
	## Specify number of K 
	## Must be up to 1-number of branches of data 
	## Will use this in the paramater later
	k.max_use <- (dim(new_eigenvectors)[2]-1)
	if(k.max_use >40){
		k.max_use <- 40
	}
	k.values <- 2:k.max_use
		
	pdf(paste0(outputdir, "Plots/CLUSTERING_isotyper_metrics_SAMPLES__Number_Clusters", type_use, "_", iso_type, ".pdf"), height=5, width=10)
	## Calculating silhoute plotts
	x <- fviz_nbclust(new_eigenvectors, FUN = hcut, method = "silhouette", k.max=k.max_use)
	avg_sil_values <- x$data
	avg_sil_values$clusters <- as.numeric(avg_sil_values$clusters)
	max1 <- Rfast::nth(avg_sil_values$y[avg_sil_values$clusters >=4], 1, descending = T)
	max1 <- avg_sil_values$clusters[avg_sil_values$y==max1]
	x <- x + geom_vline(xintercept=max1, col="red")
	plot(x)
	# Other plots 	
	plot(fviz_nbclust(new_eigenvectors, FUN = hcut, method = "wss", k.max=k.max_use))
	plot(fviz_nbclust(new_eigenvectors, FUN = hcut, method = "gap_stat", k.max=k.max_use))
	dev.off()
	print(paste0("Optimal Number of Clusters ", max1))
	
	
	w=10
	pdf(paste0(outputdir, "Plots/CLUSTERING_isotyper_metrics_SAMPLES__FINAL", type_use, "_", iso_type, ".pdf"), width=20, height=10)
	hc = hclust(res.dist, method=method_ac_use)
	par(mfrow= c(1,1), mar = c(5,5,3,3))
	plot(hc, cex=0.5, hang = -1)
	##settimg limits for trimming tree
	#threshold = quantile(dist1,limit) ## define threshold for clustering features
	rect.hclust(hc, border = "red", k=max1)
	dev.off()
}

