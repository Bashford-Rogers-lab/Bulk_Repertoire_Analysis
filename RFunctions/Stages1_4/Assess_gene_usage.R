# Function to visualise the V gene usage across the dataset to assess if any are missing
# Author: Lauren Overend 
# lauren.overend@oriel.ox.ac.uk 

#technical_replicates_file<-'/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/LEO_SEPSIS_BCR_ALL_technicals.txt'
#path_to_outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG'


visualise_vj_QC <- function(path_to_outputdir = path_to_outputdir, cluster_nodes = 5){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	library(seqinr)
	path <- path_to_outputdir
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Gene_frequencies', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	VJ_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=FALSE)
				sample <- output[1,1]
				genes <- strsplit(output$V2, "|", fixed=TRUE)
				df  <- matrix(unlist(genes), ncol=2, byrow=TRUE)
				df <- data.frame(df)
				output$V5 <- df$X1
				output$V6 <- df$X2
			    colnames(output) <- c("Sample", "VJ_Segement", "Count", "V_family", "V_Gene", "J_gene")
				return(output)
	}
	## Check VJ usage 
	frequency_variable_genes_id <- aggregate(VJ_Results$Count, by=list(V_Gene=VJ_Results$V_Gene), FUN=sum)
	colnames(frequency_variable_genes_id)[2] <- "count"
	frequency_variable_genes_family <- aggregate(VJ_Results$Count, by=list(V_family=VJ_Results$V_family), FUN=sum)
	colnames(frequency_variable_genes_family)[2] <- "count"
	frequency_variable_genes_J <- aggregate(VJ_Results$Count, by=list(J_gene=VJ_Results$J_gene), FUN=sum)
	colnames(frequency_variable_genes_J)[2] <- "count"
	
	## BCR!!!
	## Get Known V genes from library file for BCR
	if(any(frequency_variable_genes_id$V_Gene %like% "IGH")){
		Known_V_genes <- read.fasta('RFunctions/Stages1_4/Reference_nn_HOMO_SAPIENS_IGHV.fasta')
		Known_V_genes <- names(Known_V_genes)
		Known_V_genes <- str_split_fixed(Known_V_genes, "\\*", 2)[,1]
		Known_V_genes <- data.frame(unique(Known_V_genes))
		colnames(Known_V_genes) <- "Vgene"
		## See if they are all detected generate table with all known vs counts for V genes!
		#Known_V_genes$Vgene[!Known_V_genes$Vgene %in% frequency_variable_genes_id$V_Gene]
		Known_V_genes <- merge(Known_V_genes, frequency_variable_genes_id, by.x="Vgene", by.y="V_Gene", all=TRUE)
		Known_V_genes$count[is.na(Known_V_genes$count)] <- 0
		Known_V_genes$colour <- NA
		Known_V_genes$colour[Known_V_genes$count >0] <- "Black"
		Known_V_genes$colour[Known_V_genes$count ==0] <- "Red"
		pdf(paste0(path_to_outputdir, '/Plots/V_Gene_Usage_Detection.pdf'), width=20, height=10)
		p1 <- ggplot(Known_V_genes, aes(x=Vgene, y=count)) +geom_bar(stat="identity") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=Known_V_genes$colour)) +xlab("All Known V genes") + ylab("Count") + ggtitle("V gene detection in full dataset")
		plot(p1)
		dev.off()
	}
	
	## TCR!!!
	# If it is a TCR we need to use different file 
	if(any(frequency_variable_genes_id$V_Gene %like% "TR")){
		TCRusded <- frequency_variable_genes_id$V_Gene[frequency_variable_genes_id$count==max(frequency_variable_genes_id$count)]
		if(TCRusded %like% "TRG"){
			Known_V_genes <- read.fasta('RFunctions/Stages1_4/Reference_nn_HOMO_SAPIENS_TRGV.fasta')
			TCRusded <- "TRG"
		}
		if(TCRusded %like% "TRD"){
			Known_V_genes <- read.fasta('RFunctions/Stages1_4/Reference_nn_HOMO_SAPIENS_TRDV.fasta')
			TCRusded <- "TRD"
		}
		if(TCRusded %like% "TRA"){
			Known_V_genes <- read.fasta('RFunctions/Stages1_4/Reference_nn_HOMO_SAPIENS_TRAV.fasta')
			TCRusded <- "TRA"
		}
		if(TCRusded %like% "TRB"){
			Known_V_genes <- read.fasta('RFunctions/Stages1_4/Reference_nn_HOMO_SAPIENS_TRBV.fasta')
			TCRusded <- "TRB"
		}
		Known_V_genes <- names(Known_V_genes)
		Known_V_genes <- str_split_fixed(Known_V_genes, "\\*", 2)[,1]
		Known_V_genes <- data.frame(unique(Known_V_genes))
		colnames(Known_V_genes) <- "Vgene"
		## See if they are all detected generate table with all known vs counts for V genes!
		#Known_V_genes$Vgene[!Known_V_genes$Vgene %in% frequency_variable_genes_id$V_Gene]
		Known_V_genes <- merge(Known_V_genes, frequency_variable_genes_id, by.x="Vgene", by.y="V_Gene", all=TRUE)
		Known_V_genes$count[is.na(Known_V_genes$count)] <- 0
		Known_V_genes$colour <- NA
		Known_V_genes$colour[Known_V_genes$count >0] <- "Black"
		Known_V_genes$colour[Known_V_genes$count ==0] <- "Red"
		# This is for leaky V genes!
		if(!TCRusded == "TRD"){
			Known_V_genes$colour[!Known_V_genes$Vgene %like% TCRusded] <- "LightGrey"
		} else {
			Known_V_genes$colour[!Known_V_genes$Vgene %like% "TCRusded" & !Known_V_genes$Vgene %like% "DV"] <- "LightGrey"
		} 
		
		pdf(paste0(path_to_outputdir, '/Plots/V_Gene_Usage_Detection.pdf'), width=20, height=10)
		p1 <- ggplot(Known_V_genes, aes(x=Vgene, y=count)) +geom_bar(stat="identity") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=Known_V_genes$colour)) +xlab("All known V genes + V genes identified in data") + ylab("Count") + ggtitle("V gene detection in full dataset")
		plot(p1)
		dev.off()
	}
}
		