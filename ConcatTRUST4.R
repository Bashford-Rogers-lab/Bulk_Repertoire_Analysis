# Function to reduce the number of files from the TRUST4 convertor 
# Then upload these to IMGT 
# Lauren Overend 
# Lauren.overend@oriel.ox.ac.uk
# 16/11/2021


suppressMessages(library("Biostrings"))
suppressMessages(library("stringr"))
suppressMessages(library("tidyr"))
suppressMessages(library("seqinr"))
suppressMessages(library("gridExtra"))
suppressMessages(library("grid"))
suppressMessages(library("ggplot2"))
suppressMessages(library("data.table"))
########################################################
## Path to output 
outputdir <- '/well/immune-rep/shared/MISEQ/TRUST4_GAINS/'
BCRs <- paste0(outputdir, "BCR")
TCRAs <- paste0(outputdir, "TCRA")
TCRBs <- paste0(outputdir, "TCRB")
TCRGs <- paste0(outputdir, "TCRG")
TCRDs <- paste0(outputdir, "TCRD")

# list files
if (!dir.exists(paste0(outputdir, "/BCR/Reduced/"))){
dir.create(paste0(outputdir, "/BCR/Reduced/"))
} else {
    print("Dir already exists!")
}

BCR_files <- list.files(BCRs, full.names=TRUE, include.dirs=FALSE)
BCR_files <- grep(".fa", BCR_files, value=TRUE)
blank_fasta <-  DNAStringSet()
full_data <-  DNAStringSet()
count <- 0
other_count <- 0 
for(i in 1:length(BCR_files)){
	print(i)
	fastaFile <- readDNAStringSet(BCR_files[i])
	length_new <- length(fastaFile)
	if((count+length_new) < (1000000-1)){
		blank_fasta <- c(blank_fasta, fastaFile)
		count <- count+length_new
	} else {
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/BCR/Reduced/Fully_reduced_BCR_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		full_data <- c(full_data, blank_fasta)
		blank_fasta <- DNAStringSet()
		count <- length_new
		blank_fasta <- c(blank_fasta, fastaFile)
		print(paste0("New Fully Reduced File generated: ", other_count))
		}
	if(i==length(BCR_files)){
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/BCR/Reduced/Fully_reduced_BCR_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		print(paste0("New Fully Reduced File generated: ", other_count))
		print("Finished Concatenating files")
		full_data <- c(full_data, blank_fasta)
		}
}
sequence_widths <- data.frame(width(full_data)*-1)
colnames(sequence_widths) <- "Length"	
pdf(paste0(outputdir, "/BCR/Reduced/Read_Length_Full_Data.pdf"), width=10, height=5)
a <- ggplot(data=sequence_widths, aes(Length)) + geom_histogram(binwidth=10) + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("BCR")) +ylab("Count")
a1 <- ggplot(sequence_widths,aes(Length))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("BCR"))
grid.arrange(a,a1, ncol=2)
dev.off()

## Read depth per sample 
depths <- names(full_data)
depths <- sub("_.*", "", depths)
depths <- data.frame(table(depths)) 
pdf(paste0(outputdir, "/BCR/Reduced/Read_Depth_Scatter.pdf"), width=90, height=5)
a <- ggplot(depths, aes(x=depths, y=Freq)) + geom_point() + xlab("Sample")+theme_bw() +ggtitle(paste0("BCR")) +ylab("Count") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(a)
dev.off()
pdf(paste0(outputdir, "/BCR/Reduced/Read_Depth_Histo.pdf"), width=10, height=5)
a <- ggplot(depths, aes(Freq)) + geom_histogram(binwidth=50, color="black", fill="red") +theme_bw() +ggtitle(paste0("BCR")) +ylab("Count") +xlab("Read Depth") + geom_vline(aes(xintercept=mean(Freq)),color="blue", size=1)
plot(a)
dev.off()

##########################################################			
## Run for TCRA_files
if (!dir.exists(paste0(outputdir, "/TCRA/Reduced/"))){
dir.create(paste0(outputdir, "/TCRA/Reduced/"))
} else {
    print("Dir already exists!")
}
TCR_files <- list.files(TCRAs, full.names=TRUE, include.dirs=FALSE)
TCR_files <- grep(".fa", TCR_files, value=TRUE)
blank_fasta <-  DNAStringSet()
full_data <-  DNAStringSet()
count <- 0
other_count <- 0 
for(i in 1:length(TCR_files)){
	print(i)
	fastaFile <- readDNAStringSet(TCR_files[i])
	length_new <- length(fastaFile)
	if((count+length_new) < (1000000-1)){
		blank_fasta <- c(blank_fasta, fastaFile)
		count <- count+length_new
	} else {
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCRA/Reduced/Fully_reduced_TCRA_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		full_data <- c(full_data, blank_fasta)
		blank_fasta <- DNAStringSet()
		count <- length_new
		blank_fasta <- c(blank_fasta, fastaFile)
		print(paste0("New Fully Reduced File generated: ", other_count))
		
		}
	if(i==length(TCR_files)){
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCRA/Reduced/Fully_reduced_TCRA_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		print(paste0("New Fully Reduced File generated: ", other_count))
		print("Finished Concatenating files")
		full_data <- c(full_data, blank_fasta)
		}
}

sequence_widths <- data.frame(width(full_data)*-1)
colnames(sequence_widths) <- "Length"	
pdf(paste0(outputdir, "/TCRA/Reduced/Read_Length_Full_Data.pdf"), width=10, height=5)
a <- ggplot(data=sequence_widths, aes(Length)) + geom_histogram(binwidth=10) + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("TCRA")) +ylab("Count")
a1 <- ggplot(sequence_widths,aes(Length))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("TCRA"))
grid.arrange(a,a1, ncol=2)
dev.off()

## Read depth per sample 
depths <- names(full_data)
depths <- sub("_.*", "", depths)
depths <- data.frame(table(depths)) 
pdf(paste0(outputdir, "/TCRA/Reduced/Read_Depth_Scatter.pdf"), width=90, height=5)
a <- ggplot(depths, aes(x=depths, y=Freq)) + geom_point() + xlab("Sample")+theme_bw() +ggtitle(paste0("TCRA")) +ylab("Count") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(a)
dev.off()
pdf(paste0(outputdir, "/TCRA/Reduced/Read_Depth_Histo.pdf"), width=10, height=5)
a <- ggplot(depths, aes(Freq)) + geom_histogram(binwidth=10, color="black", fill="red") +theme_bw() +ggtitle(paste0("TCRA")) +ylab("Count") +xlab("Read Depth") + geom_vline(aes(xintercept=mean(Freq)),color="blue", size=1)
plot(a)
dev.off()


###################################################################
## Run for TCRB_files
if (!dir.exists(paste0(outputdir, "/TCRB/Reduced/"))){
dir.create(paste0(outputdir, "/TCRB/Reduced/"))
} else {
    print("Dir already exists!")
}
TCR_files <- list.files(TCRBs, full.names=TRUE, include.dirs=FALSE)
TCR_files <- grep(".fa", TCR_files, value=TRUE)
blank_fasta <-  DNAStringSet()
full_data <-  DNAStringSet()
count <- 0
other_count <- 0 
for(i in 1:length(TCR_files)){
	print(i)
	fastaFile <- readDNAStringSet(TCR_files[i])
	length_new <- length(fastaFile)
	if((count+length_new) < (1000000-1)){
		blank_fasta <- c(blank_fasta, fastaFile)
		count <- count+length_new
	} else {
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCRB/Reduced/Fully_reduced_TCRB_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		full_data <- c(full_data, blank_fasta)
		blank_fasta <- DNAStringSet()
		count <- length_new
		blank_fasta <- c(blank_fasta, fastaFile)
		print(paste0("New Fully Reduced File generated: ", other_count))
		
		}
	if(i==length(TCR_files)){
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCRB/Reduced/Fully_reduced_TCRB_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		print(paste0("New Fully Reduced File generated: ", other_count))
		print("Finished Concatenating files")
		full_data <- c(full_data, blank_fasta)
		}
}

sequence_widths <- data.frame(width(full_data)*-1)
colnames(sequence_widths) <- "Length"	
pdf(paste0(outputdir, "/TCRB/Reduced/Read_Length_Full_Data.pdf"), width=10, height=5)
a <- ggplot(data=sequence_widths, aes(Length)) + geom_histogram(binwidth=10) + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("TCRB")) +ylab("Count")
a1 <- ggplot(sequence_widths,aes(Length))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("TCRB"))
grid.arrange(a,a1, ncol=2)
dev.off()

## Read depth per sample 
depths <- names(full_data)
depths <- sub("_.*", "", depths)
depths <- data.frame(table(depths)) 
pdf(paste0(outputdir, "/TCRB/Reduced/Read_Depth_Scatter.pdf"), width=90, height=5)
a <- ggplot(depths, aes(x=depths, y=Freq)) + geom_point() + xlab("Sample")+theme_bw() +ggtitle(paste0("TCRB")) +ylab("Count") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(a)
dev.off()
pdf(paste0(outputdir, "/TCRB/Reduced/Read_Depth_Histo.pdf"), width=10, height=5)
a <- ggplot(depths, aes(Freq)) + geom_histogram(binwidth=10, color="black", fill="red") +theme_bw() +ggtitle(paste0("TCRB")) +ylab("Count") +xlab("Read Depth") + geom_vline(aes(xintercept=mean(Freq)),color="blue", size=1)
plot(a)
dev.off()


###################################################################
## Run for TCRG_files
if (!dir.exists(paste0(outputdir, "/TCRG/Reduced/"))){
dir.create(paste0(outputdir, "/TCRG/Reduced/"))
} else {
    print("Dir already exists!")
}
TCR_files <- list.files(TCRGs, full.names=TRUE, include.dirs=FALSE)
TCR_files <- grep(".fa", TCR_files, value=TRUE)
blank_fasta <-  DNAStringSet()
full_data <-  DNAStringSet()
count <- 0
other_count <- 0 
for(i in 1:length(TCR_files)){
	print(i)
	fastaFile <- readDNAStringSet(TCR_files[i])
	length_new <- length(fastaFile)
	if((count+length_new) < (1000000-1)){
		blank_fasta <- c(blank_fasta, fastaFile)
		count <- count+length_new
	} else {
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCRG/Reduced/Fully_reduced_TCRG_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		full_data <- c(full_data, blank_fasta)
		blank_fasta <- DNAStringSet()
		count <- length_new
		blank_fasta <- c(blank_fasta, fastaFile)
		print(paste0("New Fully Reduced File generated: ", other_count))
		
		}
	if(i==length(TCR_files)){
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCRG/Reduced/Fully_reduced_TCRG_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		print(paste0("New Fully Reduced File generated: ", other_count))
		print("Finished Concatenating files")
		full_data <- c(full_data, blank_fasta)
		}
}

sequence_widths <- data.frame(width(full_data)*-1)
colnames(sequence_widths) <- "Length"	
pdf(paste0(outputdir, "/TCRG/Reduced/Read_Length_Full_Data.pdf"), width=10, height=5)
a <- ggplot(data=sequence_widths, aes(Length)) + geom_histogram(binwidth=10) + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("TCRG")) +ylab("Count")
a1 <- ggplot(sequence_widths,aes(Length))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("TCRG"))
grid.arrange(a,a1, ncol=2)
dev.off()

## Read depth per sample 
depths <- names(full_data)
depths <- sub("_.*", "", depths)
depths <- data.frame(table(depths)) 
pdf(paste0(outputdir, "/TCRG/Reduced/Read_Depth_Scatter.pdf"), width=90, height=5)
a <- ggplot(depths, aes(x=depths, y=Freq)) + geom_point() + xlab("Sample")+theme_bw() +ggtitle(paste0("TCRG")) +ylab("Count") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(a)
dev.off()
pdf(paste0(outputdir, "/TCRG/Reduced/Read_Depth_Histo.pdf"), width=10, height=5)
a <- ggplot(depths, aes(Freq)) + geom_histogram(binwidth=10, color="black", fill="red") +theme_bw() +ggtitle(paste0("TCRG")) +ylab("Count") +xlab("Read Depth") + geom_vline(aes(xintercept=mean(Freq)),color="blue", size=1)
plot(a)
dev.off()


###################################################################
## Run for TCRD_files
if (!dir.exists(paste0(outputdir, "/TCRD/Reduced/"))){
dir.create(paste0(outputdir, "/TCRD/Reduced/"))
} else {
    print("Dir already exists!")
}
TCR_files <- list.files(TCRDs, full.names=TRUE, include.dirs=FALSE)
TCR_files <- grep(".fa", TCR_files, value=TRUE)
blank_fasta <-  DNAStringSet()
full_data <-  DNAStringSet()
count <- 0
other_count <- 0 
for(i in 1:length(TCR_files)){
	print(i)
	fastaFile <- readDNAStringSet(TCR_files[i])
	length_new <- length(fastaFile)
	if((count+length_new) < (1000000-1)){
		blank_fasta <- c(blank_fasta, fastaFile)
		count <- count+length_new
	} else {
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCRD/Reduced/Fully_reduced_TCRD_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		full_data <- c(full_data, blank_fasta)
		blank_fasta <- DNAStringSet()
		count <- length_new
		blank_fasta <- c(blank_fasta, fastaFile)
		print(paste0("New Fully Reduced File generated: ", other_count))
		
		}
	if(i==length(TCR_files)){
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCRD/Reduced/Fully_reduced_TCRD_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		print(paste0("New Fully Reduced File generated: ", other_count))
		print("Finished Concatenating files")
		full_data <- c(full_data, blank_fasta)
		}
}

sequence_widths <- data.frame(width(full_data)*-1)
colnames(sequence_widths) <- "Length"	
pdf(paste0(outputdir, "/TCRD/Reduced/Read_Length_Full_Data.pdf"), width=10, height=5)
a <- ggplot(data=sequence_widths, aes(Length)) + geom_histogram(binwidth=10) + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("TCRD")) +ylab("Count")
a1 <- ggplot(sequence_widths,aes(Length))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("TCRD"))
grid.arrange(a,a1, ncol=2)
dev.off()

## Read depth per sample 
depths <- names(full_data)
depths <- sub("_.*", "", depths)
depths <- data.frame(table(depths)) 
pdf(paste0(outputdir, "/TCRD/Reduced/Read_Depth_Scatter.pdf"), width=90, height=5)
a <- ggplot(depths, aes(x=depths, y=Freq)) + geom_point() + xlab("Sample")+theme_bw() +ggtitle(paste0("TCRD")) +ylab("Count") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(a)
dev.off()
pdf(paste0(outputdir, "/TCRD/Reduced/Read_Depth_Histo.pdf"), width=10, height=5)
a <- ggplot(depths, aes(Freq)) + geom_histogram(binwidth=10, color="black", fill="red") +theme_bw() +ggtitle(paste0("TCRD")) +ylab("Count") +xlab("Read Depth") + geom_vline(aes(xintercept=mean(Freq)),color="blue", size=1)
plot(a)
dev.off()
