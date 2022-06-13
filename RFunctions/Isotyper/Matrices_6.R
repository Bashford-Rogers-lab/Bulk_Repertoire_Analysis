make_matrices6 <- function(file=file, ids_all=ids_all){
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"chain"]))),]
		p=p[which(as.numeric(p[,"total_BCRs_counted"])>10),]
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1,-2, -3)] <- apply(p[ , c(-1,-2, -3)], 2, as.numeric)
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
		p[,  c(-1,-2, -3)] <- apply(p[ , c(-1,-2, -3)], 2, as.numeric)
		############################
		p <- p[, c(1:4)]
		
		# Calculate total mutations_per_VDJ
		blank <- NULL
		blank <- data.frame("X.sample"=as.character(), "count_type"=as.character(), "isotype"=as.character(), "mean.value"=as.numeric()) 
		for(i in 1:length(unique(p$X.sample))){
			samples <- unique(p$X.sample)[i]
			for (c in 1:length(unique(p$chain))){
				isotype <- unique(p$chain)[c]
				total <- p$mean.value[p$X.sample==samples & p$chain==isotype & p$count_type =="mean CDR_mm per BCR"] + p$mean.value[p$X.sample==samples & p$chain==isotype & p$count_type =="mean FWR_mm per BCR"]
				if(length(total)==0){
					total <- NA
				}
				row <- c(samples, "mutations_per_VDJ", isotype, total)
				blank <- rbind(blank, row)
			}
		}
		colnames(blank) <- c("X.sample", "count_type", "chain", "mean.value")
		p <- rbind(p, blank)
		
		p$count_type <- gsub(" ", "_", p$count_type)
		p$count_type[p$count_type=="mean_CDR_mm_per_BCR"] <- "CDR_mutations_per_VDJ"
		p$count_type[p$count_type=="mean_FWR_mm_per_BCR"] <- "FRAMEWORK_mutations_per_VDJ"
		p$count_type[p$count_type=="mean_CDR_FWR_ratio"] <- "CDR/FRAMEWORK_mutation_ratio_per_VDJ"
		p$count_type[p$count_type=="mean_silent_per_BCR"] <- "silent_mutations_per_VDJ"
		p$count_type[p$count_type=="mean_nonsilent_per_BCR"] <- "nonsilent_mutations_per_VDJ"
		p$count_type[p$count_type=="mean_silent_nonsilent_ratio"] <- "nonsilent/silent_mutation_ratio_per_VDJ"
		p$count_type[p$count_type=="FWR3_mm"] <- "FRAMEWORK3_mutations_per_VDJ"
		p$type <- paste0(p$count_type, "__", p$chain)
		
		## Reshape
		p <- p[, c("X.sample", "type", "mean.value")]
		colnames(p) <- c("X.sample", "type", "mean")
		l <- reshape(p, idvar = "X.sample", timevar = c("type"), direction = "wide")
		
		colnames(l) <- gsub("\\.", "_", colnames(l))
		rownames(l) <- l$X_sample
		l$X_sample <- NULL

		analysis_matrices = list(l)
		names(analysis_matrices) <- c("SHM")
		for(i in 1:length(analysis_matrices)){
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices6 = analysis_matrices
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices6 <- vector(mode = "list", length = 0)	
		}
		print("DONE 6")
		return(analysis_matrices6)
} 

