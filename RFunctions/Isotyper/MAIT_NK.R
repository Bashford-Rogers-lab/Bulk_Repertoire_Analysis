#file <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRA/ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_isotype_frequency_ALL.txt'
get_mait <- function(file, ids_all){
	info = file.info(file)
	if(info$size != 0) {
		data_all <- read.delim(file)
		data_all <- data_all[data_all$class=="TRAC",]
		data_all <- data_all[data_all$X.sample %in% ids_all,]
		data_all$frequency <- as.numeric(data_all$frequency)
		data_all$total <- NA

		for(i in 1:length(unique(data_all$X.sample))){
			sampler <- unique(data_all$X.sample)[i]
			data_all$total[data_all$X.sample==sampler] <- sum(data_all$frequency[data_all$X.sample==sampler])
			}
		data_all$prop <- data_all$frequency/data_all$total

		## Invariant

		data_all$Type <- "Standard T cell"
		data_all$Type[data_all$V.gene=="TRAV1-2" & data_all$J.gene=="TRAJ33"] <- "PossibleMAIT"
		data_all$Type[data_all$V.gene=="TRAV1-2" & data_all$J.gene=="TRAJ12"] <- "PossibleMAIT"
		data_all$Type[data_all$V.gene=="TRAV1-2" &  data_all$J.gene=="TRAJ20"] <- "PossibleMAIT"

		data_all$Type[data_all$V.gene=="TRAV1-2" & data_all$J.gene=="TRAJ9"] <- "PossibleGEM"
		data_all$Type[data_all$V.gene=="TRAV10" & data_all$J.gene=="TRAJ18"] <- "PossibleNKT"
		data_all$Type[data_all$V.gene=="TRAV24" & data_all$J.gene=="TRAJ18"] <- "PossibleNKT"

		## make a matrix of proporiton per individual
		
		l <- data.frame(xtabs(data_all$frequency ~data_all$X.sample+data_all$Type, data_all)) 
		colnames(l) <- c("Sample", "Type", "Count")

		l$prop <- NA
		for(i in 1:length(unique(l$Sample))){
			sample <-unique(l$Sample)[i] 
			total <- sum(l$Count[l$Sample==sample])
			l$prop[l$Sample==sample] <- l$Count[l$Sample==sample] / total *100
		}
		l <- l[, c("Sample", "Type", "prop")]
		z <- reshape(l, idvar = "Sample", timevar = c("Type"), direction="wide")
		rownames(z) <- z$Sample
		z$Sample <- NULL
		colnames(z) <- gsub("\\.", "_", colnames(z))
		colnames(z) <- paste0(colnames(z), "__TRAC")
		
		analysis_names = "TcellType"
		analysis_matrices = list(z)
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices15 = analysis_matrices
		} else {
		print(paste0("File: ", file, " IS EMPTY"))
		analysis_matrices15 <- vector(mode = "list", length = 0)
		}
		print("DONE NK/MAIT detection")
		return(analysis_matrices15)
}
	
		

