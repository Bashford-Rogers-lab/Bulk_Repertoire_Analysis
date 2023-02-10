file <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/ORIENTATED_SEQUENCES/ISOTYPER/All_Isotype_overlapping_frequencies_ALL.txt'


make_matrices5b <- function(file=file, ids_all=ids_all){
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.ID"]) %in% ids_all),]
		#######################
		p <- data.frame(p)
		#### We only want downsample to 50 
		p <- p[p$sample_depth==50,]
		
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		############################
		p$Class <- paste0(p$isotype1, "_", p$isotype2)
		p <- p[c("X.ID", "Class", "mean_overlap")]
		colnames(p) <- c("X.sample", "Class", "Isotype_subsampled_50_overlap_frequencies")
		l <- reshape(p, idvar = "X.sample", timevar = "Class", direction = "wide")
		colnames(l) <- gsub("\\.", "__", colnames(l))
		rownames(l) <- l$X__sample
		l$X__sample <- NULL
		#############################################################
		analysis_names = "Isotype_subsampled_overlap_frequencies"
		analysis_matrices = list(l)
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices5b = analysis_matrices
	} else { 
		print(paste0("File: ", file, " IS EMPTY"))
		analysis_matrices5b <- vector(mode = "list", length = 0)	
	}  
	print("Done 5")
	return(analysis_matrices5b)
} 


