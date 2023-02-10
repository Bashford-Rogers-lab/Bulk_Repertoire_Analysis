get_productivity <- function(path=path, ids_all=ids_all, chain_vdj=chain_vdj, counts_used=counts_used){	
	files <- list.files(path, full.name=TRUE)
	files <- grep("_1_Summary", files, value=TRUE)
	files <- grep("productive", files, value=TRUE, invert=TRUE)	
	df_list <- lapply(files, fread, header = FALSE, sep="\t", fill=TRUE, colClasses='character')
	# identify the sample they came from 
	d <- str_split(files, 'IMGT_SPLIT/') 
	d <- sapply(d, "[[", 2)  
	d <- gsub("_1_Summary.txt", "", d)
	d <- gsub("IMGT_", "", d)
	d <- gsub("BCR_", "", d)
	 
	# Name each dataframe with the run and filename
	names(df_list) <- d
	# Create combined dataframe  
	df <- df_list %>% bind_rows(.id = 'Sample')
	df <- data.frame(df)
			
	FullData <- df[!df[,4]=="rearranged sequence (but no junction found) (see comment)",]
	FullData <- FullData[!FullData[,4]=="No rearrangement found",]
	FullData <- FullData[!FullData[,4]=="No results",]		
	FullData[,4][FullData[,4]=="unproductive (see comment)"] <- "unproductive"
	FullData[,4][FullData[,4]=="productive (see comment)"] <- "productive"
	functionality <- unique(FullData[,4])
	functionality <- prop.table(table(FullData$Sample, FullData[,4]),1)

	## We must add on the chain info!
	if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
		all_class <- c("ALL")
	} else {
		counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
		receptor_type <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
		all_class <- receptor_type
	}
	colnames(functionality) <- paste0("prop_", colnames(functionality), "__", all_class)
	functionality <- as.data.frame.matrix(functionality)
	
	if(any(ids_all %like% "_productive")){
		rownames(functionality) <- paste0(rownames(functionality), "_productive")
	}
	
	if(any(ids_all %like% "_unproductive")){
		rownames(functionality) <- paste0(rownames(functionality), "_unproductive")
	}
	
	functionality <- functionality[rownames(functionality) %in% ids_all,]	
	functionality <- as.matrix(functionality)
	analysis_matrices16 = list(functionality)
	names(analysis_matrices16) <- "VDJ_Functionality"
	return(analysis_matrices16)
	print("DONE 16: VDJ Functionality")	
}
