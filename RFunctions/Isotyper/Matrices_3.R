make_matrices3 <- function(file=file, ids_all=ids_all){
	info = file.info(file)
	if(info$size != 0) {	
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		if(class(p)[1]=="character"){
			p <- as.matrix(t(p))
		} 
		#p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"isotype"]))),]
		if(class(p)[1]=="character"){
			p <- as.matrix(t(p))
		} 
		p=p[which(as.numeric(p[,"Number_of_BCRs"])>10),]
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
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		## Reshape into wide format 
		value2 = reshape(p, idvar = "X.sample", timevar = "isotype", direction = "wide")
		colnames(value2) <- gsub("\\.", "__", colnames(value2))
		cols_keep <- colnames(value2)[ colnames(value2) %like% "CDR3"]
		value2 <- value2[, c("X__sample", cols_keep)]
		rownames(value2) <- value2$X__sample
		value2$X__sample <- NULL
		############################
		analysis_matrices = list(value2)
		analysis_names = c("Mean_CDR3_Lengths")
		names(analysis_matrices) <- analysis_names
		analysis_matrices3 = analysis_matrices
		print("DONE 3")
		return(analysis_matrices3)
		} else { 
				print(paste0("File: ", file, " IS EMPTY"))
				analysis_matrices3 <- vector(mode = "list", length = 0)	
		}
}