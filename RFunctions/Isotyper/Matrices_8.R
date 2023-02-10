make_matrices8 <- function(file=file, ids_all=ids_all){
info = file.info(file)
if(info$size != 0) {
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	if(class(p)[1]=="character"){
		p <- as.matrix(t(p))
	} 
	p=p[which(as.numeric(p[,"d5"]) !=-1),]
	if(class(p)[1]=="character"){
		p <- as.matrix(t(p))
	} 
	#######################
	p <- data.frame(p)
	## Want to make the relevant columns numeric
	p[, c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
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
	colnames(p) <- c("X.sample", "isotype", "D5", "D10", "D50")
	colnames(p) <- gsub("\\.", "_", colnames(p))
	l <- reshape(p, idvar = "X_sample", timevar = "isotype", direction = "wide")
	colnames(l) <- gsub("\\.", "__", colnames(l))
	rownames(l) <- l$X_sample
	l$X_sample <- NULL	
	colnames(l) <- gsub("all", "ALL", colnames(l))
	analysis_matrices = list(l)
	analysis_names = "Clonal"
	## Replace -1 with NA (not enough sample depth for these samples!)
	for(i in 1:length(analysis_matrices)){
		analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
	}
	analysis_matrices8 = analysis_matrices
	names(analysis_matrices8) <- analysis_names	
	} else { 
		print(paste0("File: ", file, " IS EMPTY"))
		analysis_matrices8 <- vector(mode = "list", length = 0)	
	}
	print("DONE 8")
	return(analysis_matrices8)
}