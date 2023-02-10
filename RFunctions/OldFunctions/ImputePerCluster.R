list_scaled_matrices = list()
clusters_use = NULL
imputed_df <- data.frame()
features_used <- c()
Error_percluster <- data.frame()
Error1 <- c()
full_data <- data.frame( x=rep(NA, dim(mat_filtered)[1]))
full_data$x <- NULL
# Only considering clusters >= 1 features for imputatation
for(c in c(1:length(clusters))){
  features_sub = list_cluster[[c]]
  if(length(features_sub)>=1){
    print(paste0("One or more features present in cluster ", c))
    mat = p[,features_sub]
    mat_full = mat
	# Impute missing values
    if(any(is.na(mat))==TRUE){
      print(paste0("NAS present in cluster", c))
	  print("Comensing Imputation")
	  imputed_Data <- missForest(mat, variablewise = TRUE)
	  ## Get errors  per individual variable (MSE)
	  Error_clus <- imputed_Data$OOBerror
	  Error_clus <- data.frame(imputed_Data$OOBerror)
	  rownames(Error_clus) <- colnames(mat)
	  Error_clus$variable <- rownames(Error_clus)
	  Error_percluster <- rbind(Error_clus, Error_percluster) 
      ## Get errors for the whole analysis (Normalised Value allowing comparison across methods)
	  imputed_Data <- missForest(mat)
	  Error_1 <- imputed_Data$OOBerror
	  Error1 <- c(Error1, Error_1)
	  #print(paste0("Error ", Error_1))
	  completeData1= imputed_Data$ximp
	  features_used <- c(features_used, colnames(mat))
      ##########
	  comparitor1<- data.frame(mat)
	  comparitor1$Type <- "Original"
	  comparitor1$sample <- rownames(comparitor1)
	  comparitor2 <- data.frame(completeData1)
	  comparitor2$Type <- "PerCluster"
	  comparitor2$sample <- rownames(comparitor2)
	  comparitor <- rbind(comparitor1, comparitor2)
	  ###########
	  for(i in names(comparitor)[1:(length(comparitor)-2)]){
		comparitor$Imputed <- "NO" 
		imputed_samples <- rownames(comparitor[is.na(comparitor[, i]),])
		comparitor$Imputed[comparitor$sample %in% imputed_samples & comparitor$Type=="PerCluster"] <- "YES" 
		comparitor$Imputed <- factor(comparitor$Imputed, levels=c("YES", "NO"))
		comparitor$Type <- factor(comparitor$Type,levels=c("Original", "PerCluster"))
		## We are going to combine the methods so we can compare imputation 
		new <- comparitor[, c(i, "Type", "Imputed", "sample")]
		new <- reshape2::melt(new, id.vars=c("sample", "Imputed", "Type"))
		imputed_df <- rbind(imputed_df, new)
		}
		#Get data frame of all the Imputed Data!
		full_data <- cbind(full_data, completeData1)
		clusters_use = c(clusters_use, c)
		} else {
		# if there are no NaS there will be no imputation
		# we still have to bind the values
		full_data <- cbind(full_data, mat)
		}
		
}}
## Now we want to scale the data for PCA!!!
full_data <- scale(full_data)
