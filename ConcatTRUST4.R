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
TCRs <- paste0(outputdir, "TCR")

# list files
if (!dir.exists(paste0(outputdir, "/BCR/Reduced/"))){
dir.create(paste0(outputdir, "/BCR/Reduced/"))
} else {
    print("Dir already exists!")
}

BCR_files <- list.files(BCRs, full.names=TRUE, include.dirs=FALSE)
BCR_files <- grep(".fa", BCR_files, value=TRUE)

blank_fasta <-  DNAStringSet()
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
		print("New Fully Reduced File generated")
		}
	if(i==length(BCR_files)){
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/BCR/Reduced/Fully_reduced_BCR_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		print("New Fully Reduced File generated")
		print("Finished Concatenating files")
		full_data <- c(full_data, blank_fasta)
		}
}
sequence_widths <- data.frame(width(full_data)*-1)
colnames(sequence_widths) <- "Length"	
pdf(paste0(outputdir, "/BCR/Reduced/Read_Length_Full_Data.pdf"), width=10, height=5)
a <- ggplot(data=sequence_widths, aes(Length)) + geom_histogram(bins=50) + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("BCR")) +ylab("Count")
a1 <- ggplot(sequence_widths,aes(Length))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("BCR"))
grid.arrange(a,a1, ncol=2)
dev.off()

##########################################################			
## Run for TCR_files
if (!dir.exists(paste0(outputdir, "/TCR/Reduced/"))){
dir.create(paste0(outputdir, "/TCR/Reduced/"))
} else {
    print("Dir already exists!")
}
TCR_files <- list.files(TCRs, full.names=TRUE, include.dirs=FALSE)
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
		outfile_name <- paste0(outputdir, "/TCR/Reduced/Fully_reduced_TCR_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		full_data <- c(full_data, blank_fasta)
		blank_fasta <- DNAStringSet()
		count <- length_new
		blank_fasta <- c(blank_fasta, fastaFile)
		print("New Fully Reduced File generated")
		
		}
	if(i==length(TCR_files)){
		other_count <- other_count +1
		outfile_name <- paste0(outputdir, "/TCR/Reduced/Fully_reduced_TCR_", other_count, ".fa")
		writeXStringSet(blank_fasta, outfile_name, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
		print("New Fully Reduced File generated")
		print("Finished Concatenating files")
		full_data <- c(full_data, blank_fasta)
		}
}

sequence_widths <- data.frame(width(full_data)*-1)
colnames(sequence_widths) <- "Length"	
pdf(paste0(outputdir, "/TCR/Reduced/Read_Length_Full_Data.pdf"), width=10, height=5)
a <- ggplot(data=sequence_widths, aes(Length)) + geom_histogram(bins=50) + xlab("TRIMMED V-CDR3 LENGTH")+theme_bw() +ggtitle(paste0("TCR")) +ylab("Count")
a1 <- ggplot(sequence_widths,aes(Length))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("TCR"))
grid.arrange(a,a1, ncol=2)
dev.off()
###################################################################

