make_matrices12 <- function(file=file, ids_all=ids_all){
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		if(class(p)[1]=="character"){
				p <- as.matrix(t(p))
		} 
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
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
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[,  c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		############################
		p <- p[, c("X.sample", "isotype", "mean_CDR2.3_R_K_residues")]
		colnames(p) <- c("Sample", "isotype", "mean_CDR2&3_R_K_residues_CHARGE")
		s <- reshape(p, idvar = "Sample", timevar = "isotype", direction = "wide")
		colnames(s) <- gsub("\\.", "__", colnames(s))
		rownames(s) <- s$Sample
		s$Sample <- NULL	
		analysis_matrices = list(s)
		analysis_names = c("CDR_Charge")
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices12 = analysis_matrices
	} else { 
		print(paste0("File: ", file, " IS EMPTY"))
		analysis_matrices12 <- vector(mode = "list", length = 0)	
	} 
	print("DONE 12")
	return(analysis_matrices12)
}