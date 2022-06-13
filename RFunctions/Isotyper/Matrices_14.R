make_matrices14 <- function(file=file, chain_vdj=chain_vdj, ids_all=ids_all, counts_used=counts_used, iso_type=iso_type){
	info = file.info(file)
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	p <- data.frame(p)
	p$uniq_read_freq <- as.numeric(p$uniq_read_freq)		
	#####
	isotypes <- c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHEP2", "IGHGP", "IGHE")
	other_class <- c("Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "IGHD,IGHM_unmutated_singleton")
	expansion_class <- c("unexpanded", "expanded")
	####	
	if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
		all_class <- c("ALL")
	} else {
		counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
		receptor_type <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
		all_class <- receptor_type
	}
	## Total V gene usage 
	#####
	p_all <- p[p$class %in% all_class,]
	p_iso <-p[p$class %in% isotypes,] 
			
	## Calculate a percentage of repertoire which is each read 
	p_all$percent_repertoire <- NA
	for(i in unique(p_all[, "X.sample"])){
		sample_id <- i 
		sum_frequency <- sum(p_all$uniq_read_freq[p_all$X.sample==sample_id])
		p_all$percent_repertoire[p_all$X.sample==sample_id] <- ((p_all$uniq_read_freq[p_all$X.sample==sample_id])/sum_frequency)*100
	} 
	p_allx <- p_all[, c("X.sample", "percent_repertoire", "V.gene")]

	## We fill in any missing combinations with 0 as 0 percent of repertoire 
	a <- spread(p_allx, key = V.gene, value = percent_repertoire, fill=0)
	rownames(a) <- a$X.sample
	a$X.sample <- NULL
	colnames(a) <- gsub("-", "_", colnames(a))
	colnames(a) <- paste0(colnames(a), "__", all_class)
	a <- as.matrix(a)
	
	### Look at V gene Family per isotype 
	### Only relevant for BCRs!!!!
	
	if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
		p$family <- str_split_fixed(p$V.gene, "-", 2)[,1]
		p1 <- p[p$class %in% isotypes,] 
		q <- aggregate(p1$uniq_read_freq, list(p1$X.sample, p1$family, p1$class), FUN=sum) 
		colnames(q) <- c("Sample", "V.Family", "isotype", "Unique_Frequency")
		
		q$percent_repertoire <- NA
		for(i in unique(q[, "Sample"])){
			sample_id <- i 
			sum_frequency <- sum(q$Unique_Frequency[q$Sample==sample_id])
			q$percent_repertoire[q$Sample==sample_id] <- ((q$Unique_Frequency[q$Sample==sample_id])/sum_frequency)*100
		} 
		q <- q[, c("Sample", "V.Family", "isotype", "percent_repertoire")]
		q$type <- paste0(q$V.Family, "__", q$isotype)
		q <- q[, c("Sample", "type", "percent_repertoire")]
		b <- spread(q, key = type, value = percent_repertoire, fill=0)
		rownames(b) <- b$Sample
		b$Sample <- NULL
		b <- as.matrix(b)
		a <- merge(a, b, by=0)
		rownames(a) <- a$Row.names
		a$Row.names <- NULL
	}	
	analysis_matrices14 = list(a)
	names(analysis_matrices14) <- "V_GENE_USAGE"
	## Save the dataframe!
	write.table(a, paste0(outputdir, "Summary/V_Gene_usage_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	print("DONE 14: V GENE Usages")
	return(analysis_matrices14)	
} 

