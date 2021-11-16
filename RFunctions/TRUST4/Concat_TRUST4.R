#module purge
#module use -a /apps/eb/dev/ivybridge/modules/all
#module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0

# Function to convert output of TRUST4 into isotype compatible format
# Also trims reads to location of best matching VDJ prior to uploading to IMGT (constant region encoded in gene name)
# Splits into BCR/TCR 
# Lauren Overend lauren.overend@oriel.ox.ac.uk
# November 2021
# Welcome Trust Centre for Human Genetics


##############################################################################################
convert_trust4 <- function(outputdir=outputdir, sampleid=sample_id) {

	# load Recquired packages
	# all availible on rescomp Bioconductor 
	suppressMessages(library("Biostrings"))
	suppressMessages(library("stringr"))
	suppressMessages(library("tidyr"))
	suppressMessages(library("seqinr"))
	suppressMessages(library("gridExtra"))
	suppressMessages(library("grid"))
	suppressMessages(library("ggplot2"))
	suppressMessages(library("data.table"))

	# read fastq file 
	path <- paste0(outputdir, sampleid, "_annot.fa")
	fastaFile <- readDNAStringSet(path)
	# extract names of samples 
	ids <- names(fastaFile)
	all_seq <- as.character(fastaFile)

	## extracting coordinates of CDR3 region
	ids2 <- data.frame(do.call('rbind', strsplit(ids,' ',fixed=TRUE)))
	a <- unlist(regmatches(ids2$X10, gregexpr("(?<=\\().*(?=\\))", ids2$X10, perl=TRUE)))
	ids2$LocationCDR3 <- a
	b <- str_split_fixed(ids2$LocationCDR3, "-", 2)
	ids2$CDR3Start <- b[,1]
	ids2$CDR3Ends <- b[,2]

	## Give each read a unique id based on the sample name
	ids2$sampleid <- paste0(sampleid, ids2$X1)
	ids2$sampleid <- gsub("assemble", "_", ids2$sampleid)

	#Finding out isotype/constant usage 
	c <- str_split_fixed(ids2$X7, "\\*", 2)
	c <- c[,1]
	ids2$isotype <- c

	## assign just the assembly name to the sequence
	## THIS WILL BE USED LATER!!!!!!!!! 
	names(all_seq) <- ids2$X1
	sequences <- data.frame("ID"=names(all_seq), "sequence"=all_seq, row.names=NULL)

	## Finding out counts of each assembly 
	counts <- read.delim(paste0(outputdir, sampleid, "_cdr3.out"), header=FALSE)
	countsx <- data.frame(table(counts$V1))

	## Merge to get count info 
	## this will also have the benefit of removing 'incomplete' transcripts 
	ids3 <- merge(ids2, countsx, by.x="X1", by.y="Var1")
	ids3$isotype[ids3$isotype==""] <- "IX"


	## Assigning Isotype to name as in RBR isotyper pipeline
	## BCR
	matrixiso <- matrix(,nrow=dim(ids3)[1], ncol=12)
	colnames(matrixiso) <-  c("sample", c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHEP2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM"))
	for(i in 1:length(ids3$sampleid)){
		matrixiso[, "sample"][i] <- ids3$sampleid[i]
		for(c in 2:length(colnames(matrixiso))){
		#print(colnames(matrixiso)[c])
			if (ids3$isotype[i]==colnames(matrixiso)[c]){
				#print(ids3$isotype[i])
				matrixiso[i,c] <- ids3$Freq[i]
			} else {
				#print("NO")
				matrixiso[i,c] <- 0
			}
		}
	}
	matrixiso <- data.frame(matrixiso)
	# Example ID >M01488331000000000JFPJG12115213403353__0_2_0_0_0_0_0_0_0_0_0|IGHA1_IGHA2_IGHD_IGHE_IGHEP2_IGHG1_IGHG2_IGHG3_IGHG4_IGHGP_IGHM	
	FULLIDisotype <- paste0(matrixiso$sample, "__", matrixiso$IGHA1, "_",  matrixiso$IGHA2, "_", matrixiso$IGHD, "_", matrixiso$IGHE, "_", matrixiso$IGHEP2, "_", matrixiso$IGHG1, "_", matrixiso$IGHG2, "_", matrixiso$IGHG3, "_", matrixiso$IGHG4, "_", matrixiso$IGHGP, "_", matrixiso$IGHM, "|IGHA1_IGHA2_IGHD_IGHE_IGHEP2_IGHG1_IGHG2_IGHG3_IGHG4_IGHGP_IGHM")
	ids3$BCRID <- FULLIDisotype

	## TCR 
	matrixtcr <- matrix(,nrow=dim(ids3)[1], ncol=7)
	colnames(matrixtcr) <-  c("sample", c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2"))
	for(i in 1:length(ids3$sampleid)){
		matrixtcr[, "sample"][i] <- ids3$sampleid[i]
		for(c in 2:length(colnames(matrixtcr))){
		#print(colnames(matrixtcr)[c])
			if (ids3$isotype[i]==colnames(matrixtcr)[c]){
				#print(ids3$isotype[i])
				matrixtcr[i,c] <- ids3$Freq[i]
			} else {
				#print("NO")
				matrixtcr[i,c] <- 0
			}
		}
	}
	matrixtcr <- data.frame(matrixtcr)
	# Example ID >M00113396000000000JVVJY1111177789892__0_0_0_0_0_1|TRAC_TRBC1_TRBC2_TRDC_TRGC1_TRGC2
	FULLIDtcr <- paste0(matrixtcr$sample, "__", matrixtcr$TRAC, "_",  matrixtcr$TRBC1, "_", matrixtcr$TRBC2, "_", matrixtcr$TRDC, "_", matrixtcr$TRGC1, "_",  matrixtcr$TRGC2, "|TRAC_TRBC1_TRBC2_TRDC_TRGC1_TRGC2")
	ids3$TCRID <- FULLIDtcr

	## Merge so we have new name and sequence all in one dataframe 
	new <- merge(ids3, sequences, by.x="X1", by.y="ID")
	new$readlength <- (nchar(new$sequence)*-1)
	
	# Trim the sequence so that it ends at the CDR3 region
	new$trimmed <- substr(new$sequence, 0, new$CDR3Ends)
	new$trimmedlength <- (nchar(new$trimmed)*-1)
	
	# Create a directory to save summary plots
	if (!dir.exists(paste0(outputdir, "Plots"))){
		dir.create(paste0(outputdir, "Plots"))
	} else {
		print("Dir already exists!")
	}
  
	## getting start points
	## There are some sequences with no V gene annotated - lets remove those 
	new1 <- new[!new$X4=="*",]

	## String split by "'" to get the highest matching V gene and then take that column
	## Note this might mean exluding those that have different earlier start points 
	s <- data.frame(str_split_fixed(new1$X4, ",", 2)[,1])
	row.names(s) <- new1$X1
	colnames(s) <- c("VGENEloc")
	s$col2 <- unlist(regmatches(s$VGENEloc, gregexpr("(?<=\\):\\().*(?=\\):\\()", s$VGENEloc, perl=TRUE)))
	
	# Append coordinates to dataframe 
	new1$LocVstart <- s$col2
	b <- str_split_fixed(new1$LocVstart, "-", 2)
	new1$VStart <- b[,1]
	new1$VEnds <- b[,2]
	new1$CDR3Ends <- as.numeric(new1$CDR3Ends)
	new1$VStart <- as.numeric(new1$VStart)
	
	## Now we have everything we need to trim down to the V gene
	new1$Finaltrimmed <- substr(new1$sequence, new1$VStart, new1$CDR3Ends)
	new1$Finaltrimmedlength <- (nchar(new1$Finaltrimmed)*-1)

	# Plot reads with trim1, trim2, etc 
	pdf(paste0(outputdir, "Plots/", "TRUST4lengths_Trim_", sampleid, ".pdf"), width=10, height=15)
	c <- ggplot(data=new, aes(readlength)) + geom_histogram()+ xlab("Read Length") +theme_bw() +ylab("Count") +ggtitle(paste0("Pre-Trim: ", sampleid))
	c1 <- ggplot(new,aes(readlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50) + xlab("Read Length") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0("Pre-Trim: ", sampleid))
	a <- ggplot(data=new, aes(trimmedlength)) + geom_histogram()+ xlab("TRIMMED V-CDR3 LENGTH") +theme_bw() +ylab("Count") +ggtitle(paste0("Trim 1: ", sampleid))
	a1 <- ggplot(new,aes(trimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50) + xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0("Trim 1: ", sampleid))
	b <- ggplot(data=new1, aes(Finaltrimmedlength)) + geom_histogram() + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw()+ylab("Count") +ggtitle(paste0("Trim 2: ", sampleid))
	b1 <- ggplot(new1,aes(Finaltrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0("Trim 2: ", sampleid))
	grid.arrange(c, c1, a,a1, b, b1, ncol=2)
	dev.off()


	## Set up directory to output fasta files to 
	if (!dir.exists(paste0(outputdir, "TCR"))){
	dir.create(paste0(outputdir, "TCR"))
	} else {
		print("Dir already exists!")
	}

	if (!dir.exists(paste0(outputdir, "BCR"))){
	dir.create(paste0(outputdir, "BCR"))
	} else {
		print("Dir already exists!")
	}


	## Final dataset 
	## NOTE WE WILL ALSO ADD THOSE WITH NO CONSTANT REGION ASSIGNED TO THIS DATASET
	bcr_only <- new1[new1$isotype %like% "IGH" | new1$isotype == "IX",] 
	final_data_bcr <- bcr_only[, c("BCRID", "TCRID", "Finaltrimmed")]
	final_data_bcr <- as.list(c(bcr_only$Finaltrimmed))
	names(final_data_bcr) <- bcr_only$BCRID

	outfile_name <- paste0(outputdir, "BCR/", sampleid, "_BCR.fa")
	write.fasta(final_data_bcr, names(final_data_bcr), outfile_name, open = "w", nbchar = 100, as.string = FALSE)


	## Final dataset 
	## NOTE WE WILL ALSO ADD THOSE WITH NO CONSTANT REGION ASSIGNED TO THIS DATASET
	tcr_only <- new1[new1$isotype %like% "TR", ] 
	final_data_tcr <- tcr_only[, c("BCRID", "TCRID", "Finaltrimmed")]
	final_data_tcr <- as.list(c(tcr_only$Finaltrimmed))
	names(final_data_tcr) <- tcr_only$TCRID

	outfile_name <- paste0(outputdir, "TCR/", sampleid, "_TCR.fa")
	write.fasta(final_data_tcr, names(final_data_tcr), outfile_name, open = "w", nbchar = 100, as.string = FALSE)

	# Final Plot comparing TCR and BCR 
	pdf(paste0(outputdir, "Plots/TRUST4lengths_FinalTrim_bychain_", sampleid, ".pdf"))
	a <- ggplot(data=bcr_only, aes(Finaltrimmedlength)) + geom_histogram() + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("BCR: ", sampleid)) +ylab("Count")
	a1 <- ggplot(bcr_only,aes(Finaltrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("BCR: ", sampleid))
	b <- ggplot(data=tcr_only, aes(Finaltrimmedlength)) + geom_histogram() + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("TCR: ", sampleid))+ylab("Count")
	b1 <- ggplot(tcr_only,aes(Finaltrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("TCR: ", sampleid))
	grid.arrange(a,a1, b, b1, ncol=2)
	dev.off()

	# done 
}
