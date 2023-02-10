make_matrices4 <- function(file=file, ids_all=ids_all){
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.Sample"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"isotype"]))),]
		p=p[which(is.na(as.numeric(p[,"mean.mutations"]))==F),]
		p=p[which(as.numeric(p[,"number.of.unique.sequences"])>10),]
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
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
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		############################
		p <- p[, c("X.Sample", "isotype", "mean.mutations", "perc_unumtated")]
		colnames(p) <- c("X_sample", "isotype", "mean_mutations_per_vdj", "percentage_vdj_unmutated")
		l <- reshape(p, idvar = "X_sample", timevar = "isotype", direction = "wide")
		colnames(l) <- gsub("\\.", "__", colnames(l))
		rownames(l) <- l$X_sample
		l$X_sample <- NULL
		# convert to matrices
		analysis_matrices = list(l)
		analysis_names = c("Mutation")
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices4 = analysis_matrices
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices4 <- vector(mode = "list", length = 0)
		}
	print("DONE 4")
	return(analysis_matrices4)
	
}
