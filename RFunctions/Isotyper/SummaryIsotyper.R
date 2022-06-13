#------------------------------------------------------------------------------------------------
## FUNCTION FOR PLOTTING THE RESULTS OF RBR ISOTYPER SCRIPT
## ASSESSING CORRELATION WITH READ DEPTH AND MISSINGNESS OF DATA 
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk

# Recquired Packages
library(reshape2)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(stringr)
library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(data.table)
library(foreach)
library(doParallel) 
library(ggforce)
library(plot3D)
library(Peptides)
library(plyr)
library(moments)
library(mousetrap)

## Function 
summary_isotyper <- function(outputdir=outputdir, samplesfilepost=samplesfilepost, iso_type=iso_type){

iso_type=iso_type
#Getting SampleIDs
ids1 <- read.delim(samplesfilepost, sep='\t', header=FALSE)
ids <- as.character(ids1$V1)
chain_vdj <- unique(as.character(ids1$V4))
ids_all <- as.character(ids)
if(iso_type == "UNPRODUCTIVE"){
	ids_all <- paste0(ids_all, "_unproductive")
}
if(iso_type=="PRODUCTIVE"){
	ids_all <- paste0(ids_all, "_productive")
}


## Getting Read Depths for Each Sample based on RBR pipeline: 
path <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
source('RFunctions/Isotyper/ReadDepths.R')
read_depths_all <- getdepths(path, iso_type)
write.table(read_depths_all, paste0(outputdir, "Summary/Read_Depths_", iso_type, ".txt"), sep="\t", row.names=TRUE)

## Getting subsample depths which were used for isotyper script 
counts_used <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS")
all_files <- list.files(counts_used, full.name=TRUE)
all_files <- grep("depth_per_isotype", all_files, value=TRUE)
counts_used <- read.delim(all_files[1], sep="\t", header=TRUE)
counts_used <- counts_used[counts_used$type=="UNIQ",]
subsampled_depth_all <- counts_used$min[counts_used$X.isotype=="all"]
## Extract receptor information (e.g. chain usage if TCR
if(chain_vdj=="TCR" | chain_vdj=="TR"){
	counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
	receptor <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
}

## Plot read depths
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	IGHA2 = counts_used$min[counts_used$X.isotype=="IGHA2"]
	IGHA1 = counts_used$min[counts_used$X.isotype=="IGHA1"]
	IGHGP = counts_used$min[counts_used$X.isotype=="IGHGP"]
	IGHG1 = counts_used$min[counts_used$X.isotype=="IGHG1"]
	IGHG2 = counts_used$min[counts_used$X.isotype=="IGHG2"]
	IGHG3 = counts_used$min[counts_used$X.isotype=="IGHG3"]
	IGHG4 = counts_used$min[counts_used$X.isotype=="IGHG4"]
	IGHE = counts_used$min[counts_used$X.isotype=="IGHE"]
	IGHD = counts_used$min[counts_used$X.isotype=="IGHD"]
	IGHEP2 = counts_used$min[counts_used$X.isotype=="IGHEP2"]
	IGHM = counts_used$min[counts_used$X.isotype=="IGHM"]
	
	pdf(paste0(outputdir, "Plots/Read_Depths_", iso_type, ".pdf"), width=5, height=5)
	s <- ggplot(data=read_depths_all, aes(ReadDepth)) + geom_histogram() + theme_bw() + geom_vline(aes(xintercept=subsampled_depth_all, col="ALL"), show.legend=TRUE) 
	s <- s +  geom_vline(aes(xintercept=IGHA2, col="IGHA2"), show.legend=TRUE) 
	s <- s + geom_vline(aes(xintercept=IGHA1, col="IGHA1"), show.legend=TRUE, linetype = "longdash") 
	s <- s + geom_vline(aes(xintercept=IGHGP, col="IGHGP"), show.legend=TRUE)
	s <- s + geom_vline(aes(xintercept=IGHG1, col="IGHG1"), show.legend=TRUE) 
	s <- s + geom_vline(aes(xintercept=IGHG2, col="IGHG2"), show.legend=TRUE, linetype = "longdash") 
	s <- s + geom_vline(aes(xintercept=IGHG3, col="IGHG3"), show.legend=TRUE, linetype = "longdash") 
	s <- s + geom_vline(aes(xintercept=IGHG4, col="IGHG4"), show.legend=TRUE, linetype = "dotted") 
	s <- s + geom_vline(aes(xintercept=IGHE, col="IGHE"), show.legend=TRUE) 	
	s <- s + geom_vline(aes(xintercept=IGHM, col="IGHM"), show.legend=TRUE, linetype = "longdash") 
	s <- s + geom_vline(aes(xintercept=IGHD, col="IGHD"), show.legend=TRUE, linetype = "dotted") 
	s <- s + geom_vline(aes(xintercept=IGHEP2, col="IGHEP2"), show.legend=TRUE, linetype = "twodash") 
	s <- s + scale_color_manual(name = "Thresholds", values = c(ALL="red", IGHA2="BLUE", IGHA1="ORANGE", IGHGP="DARKGREEN", IGHG1="BROWN4", IGHG2="YELLOW", IGHG3="GREY", IGHG4="PURPLE", IGHE="BLACK", IGHM="LIGHTBLUE", IGHEP2="LIGHTGREEN",IGHD="DEEPPINK"))
	s <- s + facet_zoom(xlim=c(0, (subsampled_depth_all+50)))
	plot(s)
	dev.off()
} else {
	pdf(paste0(outputdir, "Plots/Read_Depths_", iso_type, ".pdf"), width=10, height=5)
	s <- ggplot(data=read_depths_all, aes(ReadDepth)) + geom_histogram() + theme_bw() + geom_vline(aes(xintercept=subsampled_depth_all, col="ALL"), show.legend=TRUE) 
	s <- s + scale_color_manual(name = "Thresholds", values = c(ALL="red"))
	plot(s)
	dev.off()
}

## Correct subsampled depth for productive/nonproductive 
if(iso_type == "UNPRODUCTIVE"){ 
	subsampled_depth_allx <- subsampled_depth_all/10
} else {
	subsampled_depth_allx <- subsampled_depth_all
} 

# Low read depth samples:
# Probably want to exlude them
samples_to_low_all <- read_depths_all$SampleIDforDepths[read_depths_all$ReadDepth < subsampled_depth_allx]

## Read Depths of excluded samples :D
read_exclude <- read_depths_all[read_depths_all$SampleIDforDepths %in% samples_to_low_all, ]

## Plot how many samples are being excluded!
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	pdf(paste0(outputdir, "Plots/Read_Depths_excluded_samples_", iso_type, ".pdf"), width=5, height=5)
	s <- ggplot(data=read_exclude, aes(ReadDepth)) + geom_histogram() + theme_bw() + geom_vline(aes(xintercept=subsampled_depth_allx, col="ALL"), show.legend=TRUE) 
	s <- s + scale_color_manual(name = "Thresholds", values = c(ALL="red"))
	plot(s)
	dev.off()
} else {
	pdf(paste0(outputdir, "Plots/Read_Depths_excluded_samples_", iso_type, ".pdf"), width=10, height=5)
	s <- ggplot(data=read_exclude, aes(ReadDepth)) + geom_histogram() + theme_bw() + geom_vline(aes(xintercept=subsampled_depth_allx, col="ALL"), show.legend=TRUE) 
	s <- s + scale_color_manual(name = "Thresholds", values = c(ALL="red"))
	plot(s)
	dev.off()
}

print("Begining Matrix Compiling")	

##Begining Data compiling:
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
## File Number 1: GINI: checked BCR!!!
source('RFunctions/Isotyper/Matrices_1.R')
file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_", iso_type, ".txt")
analysis_matrices1 <- make_matrices1(file, ids_all)

##---------------------------------------------------------------------------------------------------------------------
## File Number 2: Class Switching Summary : checked BCR!!!
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	source('RFunctions/Isotyper/Matrices_2.R')
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_per_sequence_network_parameters_", iso_type, ".txt")
	analysis_matrices2 <- make_matrices2(file, ids_all)
}

##---------------------------------------------------------------------------------------------------------------------
## File Number 3: CDR3 lengths  (averaged): checked BCR!!!
source('RFunctions/Isotyper/Matrices_3.R')
file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_lengths_overall_", iso_type, ".txt")
analysis_matrices3 <- make_matrices3(file, ids_all)

##---------------------------------------------------------------------------------------------------------------------
## File Number 4: ALL SHM: checked BCR!!!
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	source('RFunctions/Isotyper/Matrices_4.R')
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Unmutated_sequences_", iso_type, ".txt")
	analysis_matrices4 <- make_matrices4(file, ids_all)
}
	
##----------------------------------------------------------------------
## File Number 5: Isotype overlapping frequencies:  checked BCR!!!
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	source('RFunctions/Isotyper/Matrices_5A.R')
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Isotype_normalised_overlap_frequencies_uniq_", iso_type, ".txt")
	analysis_matrices5a <- make_matrices5a(file, ids_all)
}

##---------------------------------------------------------------------------------------------------------------------
## File Number 6: SHM  split across different regions :  checked BCR!!!
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	source('RFunctions/Isotyper/Matrices_6.R')
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Mutation_summmary_selection_", iso_type, ".txt")
	analysis_matrices6 <- make_matrices6(file, ids_all)
}
	
##---------------------------------------------------------------------------------------------------------------------
## File Number 7: Secondary Rearangments by isotype:  checked BCR!!!
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	source('RFunctions/Isotyper/Matrices_7.R')
	file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_", iso_type, ".txt")
	analysis_matrices7 <- make_matrices7(file, ids_all)
}

##---------------------------------------------------------------------------------------------------------------------
## File Number 8: Cluster Expansion (d5,d10, d50):  checked BCR!!!
## Corrected by Lauren!!
source('RFunctions/Isotyper/Matrices_8.R')
file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_expansion_isotype_", iso_type, ".txt")
analysis_matrices8 <- make_matrices8(file, ids_all)

##---------------------------------------------------------------------------------------------------------------------
## File Number 9: ALL V gene IGHV4_34_quantification:  checked BCR!!!
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	source('RFunctions/Isotyper/Matrices_9.R')
	file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_IGHV4_34_quantification_", iso_type, ".txt")
	analysis_matrices9 <- make_matrices9(file, ids_all)
}

##---------------------------------------------------------------------------------------------------------------------
## File Number 10: Secondary Rearangments on clone size:  checked BCR!!!
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	source('RFunctions/Isotyper/Matrices_10.R')
	file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_clone_sizes_", iso_type, ".txt")
	analysis_matrices10 <- make_matrices10(file, ids_all)
}

##---------------------------------------------------------------------------------------------------------------------
## File Number 12 CDR Charge :  checked BCR!!!
source('RFunctions/Isotyper/Matrices_12.R')
file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_CDR_charge_", iso_type, ".txt")
analysis_matrices12 <- make_matrices12(file, ids_all)

##---------------------------------------------------------------------------------------------------------------------
## File Number 13 CDR3 Charge :  checked BCR!!!
source('RFunctions/Isotyper/Matrices_13.R')
file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_charge_", iso_type, ".txt")
analysis_matrices13 <- make_matrices13(file, ids_all)

##---------------------------------------------------------------------------------------------------------------------
## File Number 14 V gene Usages:  checked BCR!!!  
source('RFunctions/Isotyper/Matrices_14.R')
file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_grouped_isotype_frequency_", iso_type, ".txt")
analysis_matrices14 <- make_matrices14(file, chain_vdj, ids_all, counts_used, iso_type)
	
##---------------------------------------------------------------------------------------------------------------------
## Get Hydrophobicity:  checked BCR!!! 
## ignore the warning I think there is an error in the package?? I've checked to see that the CDR3 contain the listed amino acids
source('RFunctions/Isotyper/Hydrophobicity.R')
all_files <- list.files(paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/", iso_type, '/Classification_per_sequence'), full.name=TRUE)
analysis_matrices15 <- get_hydrophobicity(all_files, ids_all)

##---------------------------------------------------------------------------------------------------------------------
##Skewness and kurtosis (statistics describing the CDR3 lengths):  checked BCR!!! 
source('RFunctions/Isotyper/CDR3.R')
all_files <- list.files(paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/", iso_type, '/CDR3_length_distribution'), full.name=TRUE)
analysis_matrices15b <- get_CDR3(all_files, ids_all)

##---------------------------------------------------------------------------------------------------------------------
## Get proportion productive:  checked BCR!!! 
source('RFunctions/Isotyper/Productivity.R')
path <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
analysis_matrices16 <- get_productivity(path, ids_all, chain_vdj, counts_used)

##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
## PART 2!!!!!!
## COMPOSING THE OVERALL MATRIX 

if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	print_info = c(analysis_matrices1, analysis_matrices2,  analysis_matrices3, analysis_matrices4, analysis_matrices5a,  analysis_matrices6, analysis_matrices7, analysis_matrices8, analysis_matrices9, analysis_matrices10,  analysis_matrices12, analysis_matrices13, analysis_matrices14, analysis_matrices15,analysis_matrices15b, analysis_matrices16 )
} else {
	print_info = c(analysis_matrices1, analysis_matrices3, analysis_matrices8, analysis_matrices12, analysis_matrices13, analysis_matrices14, analysis_matrices15, analysis_matrices15b, analysis_matrices16)
}

## Combining all the matrices into one big dataframe
for(i in c(1:length(print_info))){
	print_info[[i]] <- data.frame(print_info[[i]])
	print_info[[i]]$sample <- row.names(print_info[[i]])
} 
overall_matrix <- print_info %>% purrr::reduce(full_join, by="sample")
colnames(overall_matrix) <- gsub("\\.", "_", colnames(overall_matrix))
overall_matrix <- data.frame(overall_matrix)
# Merge with readdepths 
read_depths_all$SampleIDforDepths <- gsub("BCR_", "", read_depths_all$SampleIDforDepths)
read_depths_all$SampleIDforDepths <- gsub("TCR_", "", read_depths_all$SampleIDforDepths)
overall_matrix <- merge(overall_matrix, read_depths_all, by.x="sample", by.y="SampleIDforDepths")
rownames(overall_matrix) <- overall_matrix$sample

##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
## EDIT TO FILTER FOR DIFFERENT READ DEPTHS 
###############################################################################################
## GOING TO REMOVE SAMPLES FROM WHICH THE READ DEPTH IS LOWER THAN SPECIFIED
overall_matrix <- overall_matrix[overall_matrix$ReadDepth >=  subsampled_depth_allx, ]
print("###############################################################################################")
print("###############################################################################################")
print("ATTENTION!")
print(paste0("REMOVED Samples with read depth less than: ", subsampled_depth_allx))
print(paste0("REMOVED ", length(samples_to_low_all), " Samples"))
print("###############################################################################################")
print("###############################################################################################")
###############################################################################################

##---------------------------------------------------------------------------------------------------------------------
## For TCRs we must filter for the chain of interest!
if(chain_vdj %like% "T"){
	counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
	receptor <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
	if(receptor=="TRBC"){
		receptor <- c("TRBC1", "TRBC2", "TRBC")
	} else if (receptor=="TRCG" | receptor=="TRGC"  ){
		receptor <- c("TRGC1", "TRGC2", "TRGC", "TCRCG1", "TCRCG2", "TCRCG" )
	} 
	keeping_columns <- c()
	for(i in 1:length(receptor)){
		keep <- grep(receptor[i], colnames(overall_matrix), value=TRUE)
		keeping_columns <- c(keeping_columns, keep)
	}
	keeping_columns<- c("sample", keeping_columns, "ReadDepth")
	keep <- colnames(overall_matrix)[colnames(overall_matrix) %in% keeping_columns]
	overall_matrix <- overall_matrix[, c(keep)]
}
## this will be used later for naming the columns of the final matrix!
if(chain_vdj %like% "T"){
	counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
	receptor <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
}
print("Filtered for valid Metrics from file")
##---------------------------------------------------------------------------------------------------------------------
## Add in the proportions of different isotypes across groups!
## Only relevant for BCRs
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
		file <- paste0(outputdir, "Summary/IMGT/IMGT_Prop_SHM_", iso_type, ".txt")
		proportions_file <- read.delim(file, sep="\t", header=TRUE)
		colnames(proportions_file) <- gsub("\\.", "_", colnames(proportions_file))
		oo <- data.frame(str_split_fixed(colnames(proportions_file), "__", 2))
		oo <- oo$X1
		oo <- paste0(oo, "__", "ALL")
		colnames(proportions_file) <- oo
		if(iso_type=="UNPRODUCTIVE"){
			proportions_file$Sample__ALL <- paste0(proportions_file$Sample__ALL, "_unproductive")
		}
		if(iso_type=="PRODUCTIVE"){
			proportions_file$Sample__ALL <- paste0(proportions_file$Sample__ALL, "_productive")
		}	
		overall_matrix <- merge(overall_matrix, proportions_file, by.x="sample", by.y=paste0("Sample__ALL"))
		overall_matrix$Sample__ALL <- NULL
		rownames(overall_matrix) <- overall_matrix$sample
} 
overall_matrix <- overall_matrix %>% select(ReadDepth, everything())

## REPLACE ANY INSTANCES OF '-1' (meaning to small sample size so not calculated with NA!" 
################We want to add in percentage each group  (SHM - relevant for BCRS) 	
## IF we keep this we will NA the CDR3 hydrophobicity but I think its okay as I did this in each step 
empty <- apply(overall_matrix, 2, function(x){length(which(x==0 | x=="NA"))})
empty_cols <- empty[empty==dim(overall_matrix)[1]]
overall_matrix<- overall_matrix[, c(!colnames(overall_matrix) %in% names(empty_cols))]
#return(overall_matrix)
##---------------------------------------------------------------
# Save the raw UNFILTERED output!!! 
################################
out_file_table=paste0(outputdir, "Summary/All_raw_values_unfiltered_", subsampled_depth_all, "_", iso_type, ".txt")
write.table(overall_matrix, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
#################################################
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
## PART 3
## Calculating the relationship between read depth and the metric 

new <- overall_matrix

## ***** Default of corelation test will be to omit na values *****
## Calculate correlation and missingness 
# note that missingness is really percentage of samples which are not na! (so high is good!)
values <- c("Metric", "percentage_present", "pval", "correlation")
for(i in 3:(length(colnames(new)))){
	variable <- new[,i]
	id <- colnames(new)[i]
	depths <- new$ReadDepth
	data <- data.frame(cbind(variable, depths))
	#Hydrophobicity can be less than 0
	if(!id %like% "Hydrophobicity"){
		data[data=="-1"] <- "NA"
	}
	suppressWarnings(data$variable <- as.numeric(data$variable))
	suppressWarnings(data$depths <- as.numeric(data$depths))
	missingness <- data$variable[!is.na(data$variable)]
	if(length(missingness) >= (0.25*length(ids_all))){
			missingness <- length(missingness) / length(ids_all) * 100
			rval <- cor.test(data$variable,data$depths, method="spearman", exact = FALSE)
			pval <- rval$p.value
			corval <- rval$estimate
			names(corval) <- NULL
	} else {
			missingness <- length(missingness) / length(ids_all) * 100
			rval <- NA
			pval <- NA
			corval <- NA
	}
	result <- c(id, missingness, pval, corval)
	values <- rbind(values, result)
}

# Reformat data frame to annotate with significance p<0.05
# Reformat to warn with high missingness  
colnames(values) <- values[1,]
values <- values[2:length(values[,1]),]
values <- data.frame(values)
values$sig <- ""
values$pval <- as.numeric(as.character(values$pval))
values$sig[values$pval > 0.05] <- "ns"
values$sig[values$pval <= 0.05 & values$pval >= 0.001 ] <- "*"
values$sig[values$pval <= 0.001 & values$pval >= 0.0005 ] <- "**"
values$sig[values$pval <= 0.0005] <- "***"
values$sig[values$missingness < 0.6] <- "warning >60% missingness"
values$sig[is.na(values$pval)] <- "WARNING MISSINGNESS >75% STAT not Calculated"

## Fill in with subampled depth for those metrics that were subsampled. 
subsampled_data <- grep("Subsampled", values$Metric, value=TRUE)
data_metrics <- values$Metric 
data_metrics <- str_split(data_metrics, "__")
iso <-  data.frame(matrix(unlist(data_metrics), ncol=2, byrow=TRUE))
iso <- iso[,2]
values$isotype <- iso

## Identify any samples with very low detection rate!
values$percentage_present <- as.numeric(as.character(values$percentage_present))
## Assess columns with high missingness
## Use 60% as missingness threshold 
bad_columns <- values$Metric[values$percentage_present<50]
bad_columns2 <- values$Metric[values$percentage_present<25]
bad_columns3 <- values$Metric[values$percentage_present<60]
## Save Significance to a txt file 
write.table(values, paste0(outputdir, "Summary/isotyper_metrics_summary_stats", subsampled_depth_all,  "_", iso_type, ".txt"), sep="\t")
#Extract just those that are significant 
#Save file
## If there are many that are significant we may want to consider rerunning with higher
## redepth 
significant_samples <- values[values$sig=="significant",]
write.table(significant_samples, paste0(outputdir, "Summary/isotyper_metrics_summary_stats_sig_correlated_readdepth_", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t")


##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
## PART 3 
## Assessing which metrics to keeep 
## Baed on Rachaels script 
p <- as.matrix(overall_matrix)
#return(p)
## Remove read depth column 
p1 = p[,which(colnames(p)!="ReadDepth" )]
samples = p1[, "sample"]
## Converting to a numerical matrix 
p1 = p1[,c(2:length(p1[1,]))]
headers = colnames(p1)
mat = matrix(data = -1, nrow = length(samples), ncol = length(headers), dimnames = c(list(samples), list(headers)))
for(i in c(1:length(headers))){
	mat[,headers[i]] = as.numeric(p1[, headers[i]])
}
## remove all incomplete rows then columns
nz_cols = apply(mat, 2, function(x){length(which(!is.na(x)))})
nz_rows = apply(mat, 1, function(x){length(which(!is.na(x)))})
## Create plots directory if not present
if (!file.exists(paste0(outputdir, "Plots/ISOTYPER"))){
    dir.create(paste0(outputdir, "Plots/ISOTYPER"))
}
###############
# Doing plot 
pdf(paste0(outputdir, "Plots/ISOTYPER/Metrics_Missingness_Summary_", subsampled_depth_all, "_", iso_type,".pdf"), width=6, height=6)
par(mfrow= c(2,2), mar = c(5,5,3,3))
##1
nz_cols = apply(mat, 2, function(x){length(which(!is.na(x)))})
threshold = quantile(nz_cols, 0.25)
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main=paste0("run1 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat1 = mat[,which(nz_cols>= threshold)]
##2
nz_rows = apply(mat1, 1, function(x){length(which(!is.na(x)))})
threshold = quantile(nz_rows, 0.1)
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main=paste0("run2 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat2 = mat1[which(nz_rows>= threshold),]
##3
nz_cols = apply(mat2, 2, function(x){length(which(!is.na(x)))})
threshold = quantile(nz_cols, 0.25)
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main=paste0("run3 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat3 = mat2[,which(nz_cols>= threshold)]
##4
nz_rows = apply(mat3, 1, function(x){length(which(!is.na(x)))})
threshold = quantile(nz_rows, 0.1)
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main=paste0("run4 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat4 = mat3[which(nz_rows>= threshold),]
dev.off()
print("plot1")


##----------------------------------------------------------------------------------------
#################################
## FILTERING 
##----------------------------------------------------------------------------------------
## Removing repetitive annotations 
## Keeping those that are biologically informative 
## Removing those with high misingness! 
if(iso_type=="UNPRODUCTIVE"){
	threshold_new <- 90 
} else {
	threshold_new <- 60
}

o_matrix_old <- overall_matrix
## Filter those with < 60% presence 
keep <- values$Metric[values$percentage_present >= threshold_new]
overall_matrix <- overall_matrix[, c(keep)]
print(paste0("Removed ", (dim(o_matrix_old)[2]-length(keep)-2), " Metrics"))

mat_filtered <- overall_matrix 

###########
## Making plot of final filtered metrics! 

if((dim(mat_filtered))[2] >0){
	pdf(paste0(outputdir, "Plots/ISOTYPER/Metrics_Analysis_Filtered_Matrix", subsampled_depth_all, "_", iso_type,".pdf"), height=6, width=6)
	par(mfrow= c(2,2), mar = c(5,5,3,3))
	nz_rows = apply(mat_filtered, 1, function(x){length(which(!is.na(x)))})
	plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main ="run2.1")
	threshold = quantile(nz_rows, 0.1)
	segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
	mat_filtered1 = mat_filtered[which(nz_rows>= threshold),]
	nz_cols = apply(mat_filtered1, 2, function(x){length(which(!is.na(x)))})
	plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main ="run2.2")
	dev.off()
	print("plot2!")
} else {
	print("No metrics passed filtering")
}
##---------------------------------------------------------------------------------------------
## Appending chain type to columns to ensure we can merge them later

if(chain_vdj %like% "T" & (dim(mat_filtered))[2] >0){
	colnames(mat_filtered) <- paste0(receptor, "_",iso_type, "_", colnames(mat_filtered))
	colnames(mat_filtered) <- gsub(paste0("__", receptor), "", colnames(mat_filtered))
	colnames(mat_filtered) <- gsub("[.]", "_", colnames(mat_filtered))
	print("colnames edited")
} 

if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
	if((dim(mat_filtered))[2] >0){
		colnames(mat_filtered) <- paste0("BCR_", iso_type, "_READS_", colnames(mat_filtered))
		print("colnames edited")
	}
} 

## Code to replace technical replicate with mean of all technical replicates (may need to edit if has diferent name)
print("Collapsing Technical replicate by name")
control_row <-data.frame(mat_filtered[rownames(mat_filtered) %like% "JR1795_1003",])
if(dim(control_row)[1]>=1){
	control_rowid <- rownames(mat_filtered) %like% "JR1795_1003"
	control_row<- mutate_all(control_row, function(x) as.numeric(as.character(x)))
	new_row <- t(data.frame(colSums(control_row, na.rm=TRUE)/dim(control_row)[1]))
	rownames(new_row) <- "JR1795_1003_POSITIVE_MEAN"
	mat_filtered <- mat_filtered[!control_rowid,]
	mat_filtered <- rbind(mat_filtered, new_row)
}


#---------------------------------------------------------------------------------------------		
#---------------------------------------------------------------------------------------------	

## Final filtering 
## Remove features relating to pseudo genes 
print("Final Filtering...")
p <- mat_filtered
features = colnames(p)
## remove all features with very few unique values
## This likely means there are a lot of 0 values!
diff_values = apply(p[,features], 2, function(x){length(unique(x[which(is.na(x)==F)]))})
features = intersect(features, names(diff_values)[which(diff_values>=8)])
print("Filtering Performed")
xx <- length(colnames(p))-length(features)
print(paste0(xx, " Features Removed from Matrix"))
feature_columns <-   mat_filtered[,!colnames(mat_filtered) %in% grep("IGHV|TRAV|TRBV|TRGV|TRDV", colnames(mat_filtered), value=TRUE)]
## We also want to check samples for greater than 40% missingness
mat_filtered <- mat_filtered[, c(features)]
## Missingness per sample (must be less than 40%!): 
## AFter this point imputation gets messy.....
missing_threshold <- 0.4*length(feature_columns)
if(any(rowSums(is.na(mat_filtered))>missing_threshold)){
	samples_exclude <- rownames(mat_filtered)[rowSums(is.na(mat_filtered))>missing_threshold]
	mat_filtered <- mat_filtered[!rownames(mat_filtered) %in% samples_exclude,]
} 
### 
print("Done filtering on unique values filtering and missingness filtering!")
	
#---------------------------------------------------------------------------------------------	
#---------------------------------------------------------------------------------------------	
## Write the final filtered matric to an outs file - used for weighted network correlation analysis!
write.table(mat_filtered, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)

## We also want to do one with no rownames or headers for tensor decomposition analysis!!!!  
write.table(mat_filtered, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_TENSOR_FORMAT_", subsampled_depth_all, "_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
print("DONE PART 1")

## get sample and feature order which are needed for tensor analysis
sample_order <- rownames(mat_filtered)
feature_order <- colnames(mat_filtered)

write.table(sample_order, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_TENSOR_FORMAT_SAMPLE_ORDER", subsampled_depth_all, "_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(feature_order, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_TENSOR_FORMAT_FEATURE_ORDER", subsampled_depth_all, "_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)

print("Finished")
}