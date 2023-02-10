get_hydrophobicity <- function(all_files=all_files, ids_all=ids_all){	
	mean_charges <- data.frame()
	for(c in 1:length(all_files)){
		file_use <- read.delim(all_files[c])
		file_use$Hydrophobicity <- hydrophobicity(file_use$CDR3, scale = "KyteDoolittle")
		if(any(file_use$all_classes %like% ',')){
			#print(paste0("Double Assignment Present in ", all_files[c], " - will duplicate bad row"))
			bad_rows <- file_use[file_use$all_classes %like% ',',]
			bad_rows1 <- bad_rows
			bad_rows2 <- bad_rows
			two_classes <- str_split_fixed(bad_rows$all_classes, ",", 2)
			bad_rows1$all_classes <- two_classes[,1]
			bad_rows2$all_classes <- two_classes[,2]
			rownames(bad_rows2) <- paste0(rownames(bad_rows1), "_duplicate")
			all_bad_rows <- rbind(bad_rows1, bad_rows2)
			## Remove the bad row 
			file_use <- file_use[-c(as.numeric(rownames(bad_rows))),]
			## add in the dubplication row 
			file_use <- rbind(file_use, all_bad_rows)
		}
		mean_hydrophobicity <- file_use %>%group_by(all_classes) %>%summarise_at(vars(Hydrophobicity), list(CDR3_Hydrophobicity = mean, CDR3_Hydrophobicity_kurtosis=kurtosis, CDR3_Hydrophobicity_skewness=skewness, CDR3_Hydrophobicity_bimodality=bimodality_coefficient))
		all_mean <- mean(file_use$Hydrophobicity)
		all_skew <- skewness(file_use$Hydrophobicity)
		all_kurtosis <- kurtosis(file_use$Hydrophobicity)
		all_bimodality <- bimodality_coefficient(file_use$Hydrophobicity)
		new_row <- c("ALL", all_mean, all_kurtosis, all_skew, all_bimodality)
		## ALSO need to do for all!
		mean_hydrophobicity <- data.frame(mean_hydrophobicity)
		mean_hydrophobicity <- rbind(mean_hydrophobicity, new_row)
		mean_hydrophobicity$sample <- unique(file_use$X.ID)
		mean_hydrophobicity <- reshape(mean_hydrophobicity, idvar = "sample", timevar = "all_classes", direction = "wide")
		mean_hydrophobicity[mean_hydrophobicity=="NaN"] <- NA
		mean_charges <- suppressMessages(plyr::join(mean_hydrophobicity, mean_charges, type="full"))
	} 
	rownames(mean_charges) <- mean_charges$sample
	mean_charges$sample <- NULL
	mean_charges <- as.matrix(mean_charges)
	mean_charges <- mean_charges[rownames(mean_charges) %in% ids_all,]	
	colnames(mean_charges) <- gsub("\\.", "__", colnames(mean_charges))
	storage.mode(mean_charges) <- "numeric"
	analysis_matrices15 = list(mean_charges)
	names(analysis_matrices15) <- "CDR3_Hydrophobicity"
	print("DONE 15: CDR3 Hydrophobicity")
	return(analysis_matrices15)
}

