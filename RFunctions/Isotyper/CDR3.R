get_CDR3 <- function(all_files=all_files, ids_all=ids_all){
	cdr3_summary <- data.frame()
	for(c in 1:length(all_files)){
		file_use <- read.delim(all_files[c])
		cdr3_stats <- data.frame()
		for(i in 1:length(unique(file_use$isotype))){
			subset_data <- file_use[file_use$isotype==unique(file_use$isotype)[i],]
			length_dist <- rep(subset_data$CDR3_length, times=subset_data$Number_of_BCRs)
			# calculate kurotis per isotype
			## Kurtosis 3 - roughly normal, >3 it has a sharp peak, less than 3 it is very wide 
			kurtosis_new <- kurtosis(length_dist)
			if(kurtosis_new=="NaN"){
				kurtosis_new <- NA
			} 
			## Skewness 0 - symetric normally distributed, > 0 positively skewed, <0 negatively skewed (values >/< mean)
			skewness_new <- skewness(length_dist)
			if(skewness_new=="NaN"){
				skewness_new <- NA
			} 
			bimodality_new <- bimodality_coefficient(length_dist)
			row_stats <- c(unique(file_use$isotype)[i], kurtosis_new, skewness_new, bimodality_new)
			cdr3_stats <- rbind(cdr3_stats, row_stats)
		}
		colnames(cdr3_stats) <- c("class", "CDR3_Length_kurtosis", "CDR3_Length_skewness", "CDR3_Length_bimodality")
		## Need to do for all 
		all_use <- rep(file_use$CDR3_length, times=file_use$Number_of_BCRs)
		kurtosis_all <- kurtosis(all_use)
		skewness_all <- skewness(all_use)
		bimodality_all <- bimodality_coefficient(all_use)
		row_all <- c("ALL", kurtosis_all, skewness_all,bimodality_all)
		cdr3_stats <- rbind(cdr3_stats, row_all)
		cdr3_stats$sample <- unique(file_use$X.sample)
		kurtosis_df <- reshape(cdr3_stats, idvar = "sample", timevar = "class", direction = "wide")
		kurtosis_df[is.na(kurtosis_df)] <- NA
		cdr3_summary <- suppressMessages(plyr::join(kurtosis_df, cdr3_summary, type="full"))
	}

	rownames(cdr3_summary) <- cdr3_summary$sample
	cdr3_summary$sample <- NULL
	#cdr3_summary <- as.numeric(cdr3_summary)
	cdr3_summary <- as.matrix(cdr3_summary)
	cdr3_summary <- cdr3_summary[rownames(cdr3_summary) %in% ids_all,]	
	colnames(cdr3_summary) <- gsub("\\.", "__", colnames(cdr3_summary))
	storage.mode(cdr3_summary) <- "numeric"
	analysis_matrices15.b = list(cdr3_summary)
	names(analysis_matrices15.b) <- "CDR3_Summary"
	print("DONE 15b: CDR3 Kurtosis and Skewness")
	return(analysis_matrices15.b)
}

	