make_matrices1 <- function(file=file, ids_all=ids_all){
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.delim(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.Id"]) %in% ids_all),]
		if(class(p)[1]=="character"){
			p <- as.matrix(t(p))
		} 
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
		if(class(p)[1]=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		#calculate normalised scores
		p$vrenyi <- 1-(as.numeric(p[,"Vertex.Reyni"])/log(as.numeric(p[,"N.vertices"])))
		p$crenyi <- 1-(as.numeric(p[,"Cluster_Renyi"])/log(as.numeric(p[,"N_clusters"])))
		a <- -0.001
		p$vrenyi[which(p$vrenyi<0 & p$vrenyi >=a)] <- 0
		p$crenyi[which(p$crenyi<0 & p$crenyi >=a)] <- 0
		p$vrenyi[which(p$vrenyi< a)] <- NA
		p$crenyi[which(p$crenyi<a)] <- NA
		## Rename columns 
		colnames(p) <- c("X.Id", "Isotype", "N_repeats", "N_reads", "N_vertices", "N_cluster", "Vertex.Gini.Index", "Cluster_Gini_Index", "mean_total_cluster_size", "mean_vertex_size", "max_cluster_size", "max_vertex_size", "Vertex_Renyi", "Cluster_Renyi", "ThielC1", "ThielC2", "D5", "D10", "D50", "vrenyi", "crenyi")
		p <- p[,c("X.Id", "Isotype", "Vertex.Gini.Index", "Cluster_Gini_Index", "mean_total_cluster_size", "mean_vertex_size", "max_cluster_size", "max_vertex_size", "Vertex_Renyi", "Cluster_Renyi", "vrenyi", "crenyi")]
		colnames(p) <- c("X.Id", "Isotype", "Vertex_Gini_Index", "Cluster_Gini_Index", "Mean_Total_Cluster_Size", "Mean_Vertex_Size", "Maximum_Cluster_Size", "Maximum_Vertex_Size", "Vertex_Renyi", "Cluster_Renyi", "Normalised_Vertex_Renyi", "Normalised_Cluster_Renyi")
		## Make wideformat
		colnames(p) <- gsub("\\.", "_", colnames(p))
		l <- reshape(p, idvar = "X_Id", timevar = "Isotype", direction = "wide")
		colnames(l) <- gsub("\\.", "__", colnames(l))
		rownames(l) <- l$X_Id
		l$X_Id <- NULL
		colnames(l) <- gsub("all", "ALL", colnames(l))
		# Convert to matrix
		#######
		q <- list(l)
		names(q) <- "Clustering"
		analysis_matrices = q	
		## Check for -1 values which need to be replaced with NA 
		for(i in 1:length(analysis_matrices)){
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		}
		## if there is no file
		} else { 
				print(paste0("File: ", file, " IS EMPTY"))
				analysis_matrices <- vector(mode = "list", length = 0)
		} 
	print("DONE 1")
	return(analysis_matrices)
}