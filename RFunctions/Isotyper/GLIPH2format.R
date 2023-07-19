## Function to take IMGT output and convert it into format for GLIPH2 for TCRB
## Use with RBR Preprocessing pipeline
## Lauren Overend Jul 2022
## lauren.overend@oriel.ox.ac.uk

suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
#suppressMessages(library(ShortRead))


#productivity <- "ALL"
#outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRB/'

makegliph2format <- function(outputdir=outputdir, productivity=productivity){
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
	print("Extracting CDR3 information from file")
	data_all <- foreach(i = 1:length(files), .combine = rbind,  .packages=c("stringr", "data.table")) %dopar% {
		print(i)
		
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
		
		## Subset for just TCRB counts (as we have counts for all the chains)
		t <- t[,c("TRBC1", "TRBC2")]
		
		## count incidence
		t <- data.frame(sapply(t, as.numeric ))
		t$freq <- rowSums(t, na.rm=TRUE)
		
		s <- data.frame(cbind(s[,1], t$freq))
		colnames(s) <- c("ids", "count")
		rownames(s) <- s$ids
		s$ids <- NULL
		s$count <- as.numeric(s$count)

		## Get frequencies
		data <- read.delim(files[i], header=FALSE)		
		seq_id <- data[, 2]
		seq_id2 <- str_split_fixed(seq_id, "__", 2)[,2]
		ids <- str_split_fixed(seq_id, "__", 2)[,1]
		
		## Get cdr3
		cdr3 <- data.frame(data[,21])
		rownames(cdr3) <- ids
		colnames(cdr3) <- "CDR3b"
		
		# Get V gene 
		vgene <- data[, 4]
		vgene <- gsub("Homsap ", "", vgene)
		vgene <- str_split_fixed(vgene, " ", 2)[,1]
		vgene <- str_split_fixed(vgene, "\\*", 2)[,1]
		vgene <- data.frame(vgene)
		rownames(vgene) <- ids
		
		# Get J gene 
		jgene <- data[, 10]
		jgene <- gsub("Homsap ", "", jgene)
		jgene <- str_split_fixed(jgene, " ", 2)[,1]
		jgene <- str_split_fixed(jgene, "\\*", 2)[,1]
		jgene <- data.frame(jgene)
		rownames(jgene) <- ids
		
		## Put together 
		data_new <- cbind(cdr3, vgene, jgene)
		data_new$CDR3a <- "NA"
		
		## Sample 
		sample <- gsub("HV_", "HV", sample)
		subject <- str_split_fixed(sample, "_", 2)[,1]
		condition <- str_split_fixed(sample, "_", 2)[,2]
		if(outputdir %like% "SEPSIS"){
			condition <- paste0("DAY", condition)
		}
		data_new[,'subject:condition'] <- paste0(subject, ":", condition)
		## add counts
		data_new <- merge(data_new, s, by=0)
		
		## Remove those with missing cdr3
		data_new <- data_new[data_new$CDR3b!="",]
		## check for missing
		any(data_new$vgene=="")
		any(data_new$jgene=="")
		
		data_new$seqid <- data_new$Row.names
		data_new$Row.names <- NULL
		## reformatcolnames
		colnames(data_new) <- c("#CDR3b", "TRBV", "TRBJ", "CDR3a", "subject:condition", "count", "seqid")
		## Need to look for and replace this v gene functionality warnining
		data_new[, '#CDR3b'] <- gsub("\\s*\\([^\\)]+\\)","", data_new[, '#CDR3b'])

		return(data_new)
		}
		###########################
		stopCluster(cl)
		print("Concatenated Results")
		data_all <- data_all[data_all$TRBV %like% "TRBV",]
		## Some final checks
		any(data_all$vgene=="")
		any(data_all$jgene=="")
		## Checking that we have no 0 counts 
		any(data_all$count==0)
		
		## want to remove the rogue 0 counts - this will be an annotated constant region that doesnt match the V genes annotated
		data_all <- data_all[data_all$count!=0,]
		any(data_all$count==0)
		## All done and ready to go!
		data_all$seqid <- NULL
		write.table(data_all, paste0(outputdir, 'Summary/TCRBforGliph2_', productivity, '.txt'), sep="\t", row.names=FALSE, quote=FALSE)
		print("DataFrame Made and Saved")
		print("Done")
}