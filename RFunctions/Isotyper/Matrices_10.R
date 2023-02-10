make_matrices10 <- function(file=file, ids_all=ids_all){
	info = file.info(file)	
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		if(class(p)[1]=="character"){
				p <- as.matrix(t(p))
		} 
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1)] <- apply(p[ , c(-1)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)[1]=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[,  c(-1)] <- apply(p[ , c(-1)], 2, as.numeric)
		############################
		
		p <- p[, c("X.sample", "d5_norm", "d5_secondary","mean_clone_size_norm", "X_mean_clone_size_secondary")]
		colnames(p) <- c("X.sample", "d5_norm__ALL", "d5_secondary__ALL", "mean_clone_size_norm__ALL", "mean_clone_size_secondary__ALL")
		rownames(p) <- p$X.sample
		p$X.sample <- NULL

		analysis_matrices = list(p)
		analysis_names = c("V_gene_replacement_clonal_expansion")
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "_", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

		} 
		analysis_matrices10 = analysis_matrices
		} else { 
				print(paste0("File: ", file, " IS EMPTY"))
				analysis_matrices10 <- vector(mode = "list", length = 0)	
		}
		print("DONE 10")
		return(analysis_matrices10)
}

