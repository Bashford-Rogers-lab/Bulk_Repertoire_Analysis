make_matrices7 <- function(file=file, ids_all=ids_all){
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		p=p[grep("IGH",as.character(p[,"chain"])),]
		w = setdiff(p[,"chain"], p[,"chain"][grep("IGHV",as.character(p[,"chain"]))])
		p=p[which(as.character(p[,"chain"]) %in% w),]
		w = setdiff(p[,"chain"], p[,"chain"][grep("mut",as.character(p[,"chain"]))])
		p=p[which(as.character(p[,"chain"]) %in% w),]
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		p <- p[p$total>=50,]
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
		p[,  c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		############################
		p <- p[p$total >50,]
		p <- p[, c("X.sample", "chain", "percentage")]
		colnames(p) <- c("X.sample", "chain", "V_gene_replacement_frequency")
		# reshape
		colnames(p) <- gsub("\\.", "_", colnames(p))
		l <- reshape(p, idvar = "X_sample", timevar = "chain", direction = "wide")
		colnames(l) <- gsub("\\.", "__", colnames(l))
		rownames(l) <- l$X_sample
		l$X_sample <- NULL
		
		# Reformat 
		analysis_matrices = list(l)
		analysis_names = c("V_gene_replacement_frequency")
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
				analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices7 = analysis_matrices
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices7 <- vector(mode = "list", length = 0)
		}
	print("DONE 7")
	return(analysis_matrices7)
}