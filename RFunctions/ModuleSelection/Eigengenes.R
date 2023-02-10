## Function to get Eigen vectors from scaled and filtered matrix 
## Lauren Overend and Rachael Bashford-Rogers
## lauren.overend@oriel.ox.ac.uk
## Jan 2022 

## For development
#scaled_matrix <- full_data
#cluster_assignment <- list_cluster
#outputdir
#type <- type_use
#subsampled_deptha
#iso_type

library(factoextra)

get_eigengenes <- function(scaled_matrix, cluster_assignment, outputdir, type, subsampled_deptha, iso_type){
		
		full_data <- scaled_matrix
		list_cluster <- cluster_assignment
		full_data[is.na(full_data)] <- 0 
		## set up plot shape
		nclus = length(list_cluster)
		ncols= ceiling(sqrt(nclus))
		w=2.5
		pdf(paste0(outputdir, "Plots/TRY_EIGEN_scaleTRUE_",  type, "_", subsampled_deptha, "_", iso_type,".pdf"), height=w*1*5, width=w*1*5)
		par(mfrow= c(ncols,ncols), mar = c(5,5,3,3))
		list_eigenvectors = NULL
		columns_single <- c()
		
		#### Extra
		variance_explained <- c()
		### Loadings 
		feature_loadings <- NULL
		plot_list <- list()
		plot_list2 <- list()
		
		for(k in 1:length(list_cluster)){ 
			if(length(list_cluster[[k]])>=2){
				cols_use <- list_cluster[[k]]
				data_cluster <- full_data[, c(cols_use)]
				eigenvector = apply(data_cluster, 1, mean)
				blah <- prcomp(data_cluster)
				### THIS IS THE ACTUAL EIGENVECTOR 
				eigenvector2 <- blah$x[,1]
				
				## Get the loading score of each feature 
				loadings_pca1 <- blah$rotation[,1]
				
				### Do we need to invert them to correlate with mean?
				if(cor(eigenvector, eigenvector2)<0){
					print(paste0("Module ", k, " mean and eigenvector are anticorrelated... converting"))
					loadings_pca1 <- loadings_pca1*-1
				}
				
				loadings_pca1 <- sort(loadings_pca1)
				loadings_pca11 <- list(loadings_pca1)
				names(loadings_pca11) <- paste0("Module_", k)
				feature_loadings <- c(feature_loadings, loadings_pca11)
				
				## One feature worth of loadings. 
				ft_cutoff <- sqrt(1/length(cols_use))
				
				colourx <- c("orange", rep("blue", 9))
				p1 <- fviz_eig(blah, main= paste0("Module ", k), ylab="Proportion of Variance", xlab="Principle Component",  ylim=c(0,100)) 
				plot_list[[k]] <- p1
				
				loadings_pca1<- data.frame(loadings_pca1)
				colnames(loadings_pca1) <- "Loading_Score"
				loadings_pca1$feature <- rownames(loadings_pca1)
				orderx <- rownames(loadings_pca1)
				loadings_pca1$feature  <- factor(loadings_pca1$feature, levels=c(orderx))
				loadings_pca1$Loading_Score <- as.numeric(loadings_pca1$Loading_Score)
				loadings_pca1$important <- NA
				loadings_pca1$important[abs(loadings_pca1$Loading_Score)>=ft_cutoff] <- "True"
				loadings_pca1$important[abs(loadings_pca1$Loading_Score)< ft_cutoff] <- "False"

				p2x <- ggplot(loadings_pca1, aes(x=feature, y=Loading_Score, colour=important))+ggtitle(paste0("Module ", k))+xlab("Feature")+ylab("PCA1 Loading Score") +geom_hline(yintercept=0, col="lightblue", type="dashed") + geom_point()+geom_hline(yintercept=(-1*ft_cutoff), col="red") +geom_hline(yintercept=ft_cutoff, col="red")+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ labs(colour="Important\nContributor") +theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
				plot_list2[[k]] <- p2x
				
				###########################################################
				variance <- data.frame(t(data.frame(as.matrix(summary(blah)$importance[,1]))))
				variance$Module <- paste0("Module_", k)
				rownames(variance) <- variance$Module
				variance_explained <- rbind(variance, variance_explained)
				## Make sure the eigenvector and mean are in the same direction!
				if(cor(eigenvector, eigenvector2)<0){
					print(paste0("Module ", k, " mean and eigenvector are anticorrelated... converting"))
					eigenvector2 <- eigenvector2*-1
				}
				
				plist <-  list(eigenvector2)
				names(plist) <- paste0("Module_", k)
				list_eigenvectors = c(list_eigenvectors, plist)
				# order samples and plot
				order_samples = order(rowSums(data_cluster))
				mat_scale_ordered = data_cluster[order_samples,]
				eigenvector_ordered = eigenvector[order_samples]
				eigenvector_ordered2 = eigenvector2[order_samples]
				range_y = range(data_cluster)
				range_y[1] <- range_y[1]-3
				range_y[2] <- range_y[2]+3
				range_x = c(1,length(mat_scale_ordered[,1]))
				plot(range_x, range_y, col = "white", pch = 21, main = paste0(c("Module ", k)), xlab = "ordered sample", ylab = "value")
				col = "grey"
				x = c(1:length(mat_scale_ordered[,1]))
				for(i in c(1:length(mat_scale_ordered[1,]))){
					points(x,mat_scale_ordered[,i], type = "l", col = col, lwd = 0.5, lty = 1)
				 }
				 col1 = add.alpha("red", alpha = 0.5)
				 points(x,eigenvector_ordered, type = "l", col = col1, lwd = 2, lty = 1)
				 col1 = add.alpha("blue", alpha = 0.5)
				 points(x,eigenvector_ordered2, type = "l", col = col1, lwd = 2, lty = 1)
				} else {
				columns_single <- c(columns_single, list_cluster[[k]])
				}
		}
	dev.off()
	
	variance_explained$Module <- factor(variance_explained$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25", "Module_26", "Module_27", "Module_28", "Module_29", "Module_30",  "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40" ))
	##Lets plot how much variation is explained by each module 
	variance_explained$Module <- factor(variance_explained$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25", "Module_26", "Module_27", "Module_28", "Module_29", "Module_30",  "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40" ))
	pdf(paste0(outputdir, "Plots/VarianceExplained_",  type, "_", subsampled_deptha, "_", iso_type,".pdf"), height=5, width=6)
	plot(ggplot(variance_explained, aes(x=Module, y=Proportion.of.Variance, fill=Module))+geom_col(colour="black") + theme_classic() +guides(fill="none") + xlab("Module") +ylab("Proportion of Variance")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	ggtitle("Variance explained by PC1 for each repertoire cluster"))
	dev.off()
	##save these stats 
	write.table(variance_explained, paste0(outputdir, "Summary/VarianceExplainedPC1_",  type, "_", subsampled_deptha, "_", iso_type,".txt"), sep="\t")

	## Lets plot the variance explained 
	col_number <- 5
	pdf(paste0(outputdir, "Plots/Module_PCAS_VARIANCE_",  type, "_", subsampled_deptha, "_", iso_type,".pdf"), height=11, width=12)
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
	multiplot(plotlist = plot_list, cols = ncols)
	dev.off()
	
	## Lets plot the loadings
	pdf(paste0(outputdir, "Plots/Module_PCAS_LOADINGS_",  type, "_", subsampled_deptha, "_", iso_type,".pdf"), height=15, width=16)
	multiplot(plotlist = plot_list2, cols = ncols)
	dev.off()
	
	## Lets save the list of loading scores
	saveRDS(feature_loadings, paste0(outputdir, "Summary/Module_FeaturePCALoadings_",  type, "_", subsampled_deptha, "_", iso_type,".rds"))
	
	#####
	## Tidy up to return
	mat_eigenvectors = NULL
	for(c in c(1:length(list_eigenvectors))){
	  if(length(mat_eigenvectors)==0){mat_eigenvectors = list_eigenvectors[[c]]
	  }else{mat_eigenvectors = cbind(mat_eigenvectors, list_eigenvectors[[c]])}
	}
	colnames(mat_eigenvectors) = names(list_eigenvectors)
	rownames(mat_eigenvectors) <- gsub("_productive", "", rownames(mat_eigenvectors))
	rownames(mat_eigenvectors) <- gsub("_unproductive", "", rownames(mat_eigenvectors))
	## Columns that werent included in the PCA!!!
	if(length(columns_single)>0){
		other_columns <- full_data[, c(columns_single), drop = FALSE]
		colnames(other_columns) <- paste0("Scale_", colnames(other_columns))
		#colnames(other_columns) <- gsub("_PRODUCTIVE", "", colnames(other_columns))
		rownames(other_columns) <- gsub("_productive", "", rownames(other_columns))
		mat_eigenvectors<- merge(mat_eigenvectors, other_columns, by=0)
		rownames(mat_eigenvectors) <- mat_eigenvectors$Row.names 
		mat_eigenvectors$Row.names <- NULL
	}
	return(mat_eigenvectors)

	}

		