## Function to get Eigen vectors from scaled and filtered matrix 
## Lauren Overend and Rachael Bashford-Rogers
## lauren.overend@oriel.ox.ac.uk
## Jan 2022 

get_eigengenes <- function(scaled_matrix, cluster_assignment, outputdir, type, subsampled_deptha, iso_type){
		
		full_data <- scaled_matrix
		list_cluster <- cluster_assignment
		
		w=2.5
		pdf(paste0(outputdir, "Plots/TRY_EIGEN_scaleTRUE_",  type, "_", subsampled_deptha, "_", iso_type,".pdf"), height=w*1*5, width=w*1*6)
		par(mfrow= c(5,6), mar = c(5,5,3,3))
		list_eigenvectors = NULL
		columns_single <- c()
		for(k in 1:length(list_cluster)){ 
			if(length(list_cluster[[k]])>=2){
				cols_use <- list_cluster[[k]]
				data_cluster <- full_data[, c(cols_use)]
				eigenvector = apply(data_cluster, 1, mean)
				blah <- prcomp(data_cluster)
				eigenvector2 <- blah$x[,1]
				
				## Make sure the eigenvector and mean are in the same direction!
				if(cor(eigenvector, eigenvector2)<0){
					eigenvector2 <- eigenvector2*-1
				}
				
				plist <-  list(eigenvector)
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
				plot(range_x, range_y, col = "white", pch = 21, main = paste0(c("module ", k)), xlab = "ordered sample", ylab = "value")
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

		