get_read_depths <- function(outputdir=outputdir, productivity=productivity, type=type){

	files <- list.files(paste0(outputdir,'ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT'), full.name=TRUE)
	files <- grep("_Summary.txt", files, value=TRUE)
	if(productivity=="PRODUCTIVE"){
		files <- grep("_productive_1_Summary", files, value=TRUE)
	}
	if(productivity=="UNPRODUCTIVE"){
		files <- grep("_unproductive_1_Summary", files, value=TRUE)
	}
	if(productivity=="ALL"){
		files <- grep("_Summary.txt", files, value=TRUE)
		files <- grep("unproductive|productive", files, value=TRUE, invert=TRUE)
	}
	
	cl <- parallel::makeCluster(10)
	doParallel::registerDoParallel(cl)
	depths_all <- foreach(i = 1:length(files), .combine = rbind, .packages=c("stringr", "data.table")) %dopar% {

		## Get the sample id
		sample <- files[i]
		sample <- gsub(paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT/IMGT_"), "", sample)
		sample <- gsub("_unproductive_1_Summary.txt", "", sample)
		sample <- gsub("_productive_1_Summary.txt", "", sample)
		sample <- gsub("_1_Summary.txt", "", sample)
	
		## locate counts file 
		s <- read.delim(paste0(outputdir, 'ORIENTATED_SEQUENCES/NETWORKS/Att_', sample, '.txt'), header=FALSE)
		s <- s[,1]
		s <- gsub("\\|", "__", s)
		s <- str_split_fixed(s, "__", 3)
		## get TCRB counts
		cols <- s[,3]
		col_number <- (sum(str_count(cols, "_"))/length(cols))+1
		col_names <- as.vector(str_split_fixed(cols[1], "_", col_number))

		s <- s[,1:2]
		t <- str_split_fixed(s[,2], "_", col_number)
		t <- data.frame(t)
		colnames(t) <- col_names
		
		t <- data.frame(sapply(t, as.numeric ))
		t$ALL <- rowSums(t)
		rownames(t) <- s[,1]
		
		## Now we need to only use ids in imgt file 
		data <- read.delim(files[i], header=FALSE)
		data <- data[(data$V3=="productive (see comment)" | data$V3=="productive" | data$V3=="unproductive (see comment)" | data$V3=="unproductive"), ]

		
		seq_id <- data[, 2]
		seq_id2 <- str_split_fixed(seq_id, "__", 2)[,2]
		ids <- str_split_fixed(seq_id, "__", 2)[,1]
		
		#### For Total reads
		t <- t[rownames(t) %in% ids,]
		depths <- colSums(t)
		depths <- as.data.frame(t(data.frame(depths)))
		depths$type <- "TOTAL"
		
		## For Unique Reads 
		
		dunique <- colSums(t !=0)
		dunique <- as.data.frame(t(data.frame(dunique)))
		dunique$type <- "UNIQUE"
		depths <- rbind(dunique, depths)
		 
		 ## Want to save ids and counts ## can use this as look up table for other analysis
		t$Sample <- paste0(sample, "_", productivity)
		t$seqid <- rownames(t)
		
		if (!dir.exists(paste0(outputdir, "Summary/SummarySequences/"))){
			dir.create(paste0(outputdir, "Summary/SummarySequences/"))
		} 
		write.table(t, paste0(outputdir, "Summary/SummarySequences/", sample, "_", productivity, "_allsequences.txt"), sep="\t", row.names=FALSE)
		####################
		
		
		if(type =="IGH"){
			all_types <- c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHEP2",  "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM", "ALL")
			found <- colnames(depths) 
			not_found <- all_types[!all_types %in% found]
			if(length(not_found)>0){
				#print(i)
				depths[, c(not_found)] <- 0
			}
			depths <- depths[,c(all_types, "type")]
		}
		
		if(type %like% "TR"){
			all_types <- c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "ALL")
			found <- colnames(depths) 
			not_found <- all_types[!all_types %in% found]
			if(length(not_found)>0){
				#print(i)
				depths[, c(not_found)] <- 0
			}
			depths <- depths[,c(all_types, "type")]
		}
		depths$sample <- sample
		return(depths)
		}
	
		
	depths_allx <- depths_all[, c((dim(depths_all)[2]), (dim(depths_all)[2]-1), 1:(dim(depths_all)[2]-2))]
	
	if(type %like% "TR"){
		depths_allx$TRBC <- depths_allx$TRBC1 +  depths_allx$TRBC2
		depths_allx$TRGC <- depths_allx$TRGC1 +  depths_allx$TRGC2
	}

	depths_all_unique <- depths_allx[ depths_allx$type=="UNIQUE",]
	depths_all_total <- depths_allx[depths_allx$type=="TOTAL",]
	write.table(depths_all_unique, paste0(outputdir, "Summary/Read_Depths_per_isotype_unique", productivity, ".txt"), sep="\t", row.names=FALSE)
	write.table(depths_all_total, paste0(outputdir, "Summary/Read_Depths_per_isotype_all", productivity, ".txt"), sep="\t", row.names=FALSE)
	stopCluster(cl)
	
	allx <- rbind(depths_all_unique, depths_all_total)
	return(allx)
}
	