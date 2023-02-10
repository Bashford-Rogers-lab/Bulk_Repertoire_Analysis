eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	
	if(!is.na(modules_keep)){
		eigenvectors <- eigenvectors[, c( modules_keep,  "sample", "DAY", "DISEASE")]
	}
		
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	## Include for just samples that are in the dataset (after filtering etc)
	metadata <- metadata[metadata$SampleID_alternative %in% rownames(eigenvectors),]
	
	## We also want to add read depths for the analysis!!!
	if(type_receptor %like% "BCR"){
		read_depths <- read.delim(paste0(outputdir, "/Summary/Read_Depths_per_isotype_allPRODUCTIVE.txt"))
		read_depths <- read_depths[,-2]
		cols_keep <- c("sample")
		features <- read.delim(paste0(outputdir, "/Summary/Clustered_Features_assignment_BCR_PRODUCTIVE_NON_IMPUTED.txt"))
		for(i in 2:length(colnames(read_depths))){
			if(any(features[,1] %like% colnames(read_depths)[i])){
				cols_keep <- c(cols_keep,  colnames(read_depths)[i])
			}
		}
		read_depths <- read_depths[,c(cols_keep)]
	}
	
	if(type_receptor %like% "TCRAB"){
		outdir1 <- gsub("TCRAB", "TCRA", outputdir)
		read_depths1 <- read.delim(paste0(outdir1, "/Summary/Read_Depths_per_isotype_allPRODUCTIVE.txt"))
		read_depths1 <- read_depths1[,c("sample", "TRAC")]
		outdir2 <- gsub("TCRAB", "TCRB", outputdir)
		read_depths2 <- read.delim(paste0(outdir2, "/Summary/Read_Depths_per_isotype_allPRODUCTIVE.txt"))
		read_depths2 <- read_depths2[,c("sample", "TRBC")]
		read_depths <- merge(read_depths1, read_depths2, by="sample")
	}
	
	if(type_receptor %like% "TCRGD"){
		outdir1 <- gsub("TCRGD", "TCRG", outputdir)
		read_depths1 <- read.delim(paste0(outdir1, "/Summary/Read_Depths_per_isotype_allPRODUCTIVE.txt"))
		read_depths1 <- read_depths1[,c("sample", "TRGC")]
		outdir2 <- gsub("TCRGD", "TCRD", outputdir)
		read_depths2 <- read.delim(paste0(outdir2, "/Summary/Read_Depths_per_isotype_allPRODUCTIVE.txt"))
		read_depths2 <- read_depths2[,c("sample", "TRDC")]
		read_depths <- merge(read_depths1, read_depths2, by="sample")
	}
	
	colnames(read_depths) <- paste0("Total_Read_Depth.", colnames(read_depths))
	colnames(read_depths)[1] <- "Total_Read_Depth.sample"
	metadata <- merge(metadata, read_depths, by.x="SampleID_alternative", by.y="Total_Read_Depth.sample")