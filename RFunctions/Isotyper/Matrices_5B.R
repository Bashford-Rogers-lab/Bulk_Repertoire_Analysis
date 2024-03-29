make_matrices5b <- function(file=file){
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("IGHG4",as.character(p[,"iso1"]))),]
		p=p[setdiff(c(1:length(p[,1])), grep("IGHG4",as.character(p[,"iso2"]))),]

		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1, -4, -5)] <- apply(p[ , c(-1, -4, -5)], 2, as.numeric)
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
		p[, c(-1, -4, -5)] <- apply(p[ , c(-1, -4, -5)], 2, as.numeric)
		############################
		p$Class <- paste0(p$iso1, "_", p$iso2)
		p <- p[c("X.sample", "Class", "mean_overlap")]
		colnames(p) <- c("X.sample", "Class", "Isotype_normalised_overlap_frequencies")
		l <- reshape(p, idvar = "X.sample", timevar = "Class", direction = "wide")
		colnames(l) <- gsub("\\.", "__", colnames(l))
		rownames(l) <- l$X__sample
		l$X__sample <- NULL
		l$Isotype_normalised_overlap_frequencies__all_0__0 <- NULL
		
		analysis_names = "Isotype_normalised_overlap_frequencies"
		analysis_matrices = list(l)
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices5 = analysis_matrices
		print("Done5A")
	} else { 
		print(paste0("File: ", file, " IS EMPTY"))
		analysis_matrices5 <- vector(mode = "list", length = 0)	
	}  
} 