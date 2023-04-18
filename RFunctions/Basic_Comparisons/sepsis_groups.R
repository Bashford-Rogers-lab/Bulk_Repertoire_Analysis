library(ggplot2) 
library(data.table) 
library(stringr)
library(ggpubr)
library(rmcorr) 
library(lme4)
library(lmerTest)
library(ggplot2) 
library(M3C) 
library(lmerTest)

cohort1_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/logcpm_181_20600_cohort2.txt'
cohort2_file <- "/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/Logcpm_864_20416.txt"
cohort_assignment <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/RepertoireRNAseqCohorts.txt'
eigenvectors_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/'
genes_use <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/GeneExpPlots/ALL_GENES/GeneTables/Module_19All_genes_LME_significant.txt'
module_info <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_SCALED_BCR_PRODUCTIVE.txt'
	
function_umap <- function( cohort1_file, cohort2_file, eigenvectors_file, outputdir, cohort_assignment, module_use))
		
		## this will give us ensembl id and gene name 
		gtf <- read.delim('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Homo_sapiens.GRCh38.94_gene_annotation_table.txt', sep="\t")
	    gtf_lookup <- gtf[, c("gene_id", "GeneSymbol")]

		if(!dir.exists(outputdir)) {dir.create(outputdir)}
		plot_dir <- paste0(outputdir, "ALL_GENES/")
		if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
		###############################################################################
		### Read in files of interest (gene expression, sample cohort assignment etc 
		cohort2 <- t(read.delim(cohort1_file, header=TRUE))
		cohort1 <- t(read.delim(cohort2_file, header=TRUE))
		repertoire <- read.delim(cohort_assignment)
		## Clean UP 
		repertoire$RNACohort[repertoire$RNACohort =="Cohort1 "] <- "Cohort1"
		repertoire$RNACohort[repertoire$RNACohort =="Cohort2 "] <- "Cohort2"
		repertoire$RNACohort[repertoire$RNACohort =="NA "] <- NA
		repertoire$RNACohort[repertoire$RNACohort =="Exclude "] <- "Exclude"
		repertoire <- repertoire[!is.na(repertoire$RNACohort),]

		### Subset to get gene expression
		#### Cohort1 gene expression!!!!!!!!!
		cohort1_samples <- repertoire[repertoire$RNACohort=="Cohort1",]
		cohort1_geneexp <- cohort1[rownames(cohort1) %in% cohort1_samples$SampleID,]
		## Cohort2 gene expression
		cohort2_samples <- repertoire[repertoire$RNACohort=="Cohort2",]
		cohort2_geneexp <- cohort2[rownames(cohort2) %in% cohort2_samples$SampleID,]

		## Now we want to merge it together 
		common_col_names <- intersect(colnames(cohort2_geneexp), colnames(cohort1_geneexp))
		gene_expression <- rbind(subset(cohort2_geneexp, select = common_col_names), subset(cohort1_geneexp, select = common_col_names))
		
		###########################################################################################
		## Now lets get the module scores and subset to remove the bad ids 
		eigenvectors <- read.delim(eigenvectors_file, sep="\t", header=TRUE)
		bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
		eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
		
		################################################################
		gene_exp <- merge(gene_expression, eigenvectors, by=0)
		gene_exp$Barcode <- NA
				for(x in 1:length(gene_exp$sample)){
					gene_exp$Barcode[x] <- str_split_fixed (gene_exp$sample[x], "_", 2)[,1]
					if(gene_exp$sample[x] %like% "HV"){
					 gene_exp$Barcode[x] <- paste0(str_split_fixed(gene_exp$sample[x], "_", 3)[,1], "_", str_split_fixed(gene_exp$sample[x], "_", 3)[,2])
					}
				}
		
		gene_exp$DAY <- gsub("Day", "", gene_exp$DAY)
		gene_exp$DAY <-as.numeric(gene_exp$DAY)
		
		#####
		genes_subset <- read.delim(genes_use)
		genes_subset <- genes_subset$gene[genes_subset$Type !="AIR"]
		
		gene_new <- gene_exp[,c(genes_subset, "DAY", "DISEASE", "sample", "Barcode")]
		
		pdf(paste0(outputdir, "Allgenes_TSNE19.pdf"), width=15, height=5)
		x1 <- tsne(t(gene_new[,c(genes_subset)]), labels=as.factor(gene_new$DISEASE))+ggtitle("TSNE")
		x2 <- umap(t(gene_new[,c(genes_subset)]), labels=as.factor(gene_new$DISEASE))+ggtitle("UMAP")
		plot(ggarrange(x1, x2, ncol=2))
		dev.off()
		
		#### we want to add in clonal expansion etc!!!!!
		features <- read.delim(module_info)
		rownames(features) <- gsub("_productive", "", rownames(features))
		
		rownames(gene_new) <- gene_new$sample
		gene_new2 <- merge(gene_new, features, by=0)
		
		pdf(paste0(outputdir, "Allgenes_TSNE19.pdf"), width=15, height=5)
		x1 <- tsne(t(gene_new2[,c(genes_subset)]), labels=scale(gene_new2$BCR_READS_Cluster_Gini_Index__ALL), controlscale = TRUE, scale=1)+ggtitle("TSNE")
		x2 <- umap(t(gene_new2[,c(genes_subset)]), labels=scale(gene_new2$BCR_READS_Cluster_Gini_Index__ALL), controlscale = TRUE, scale=1)+ggtitle("UMAP")
		plot(ggarrange(x1, x2, ncol=2))
		x1 <- tsne(t(gene_new2[,c(genes_subset)]), labels=scale(gene_new2$BCR_READS_D50__ALL), controlscale = TRUE, scale=1)+ggtitle("TSNE")
		x2 <- umap(t(gene_new2[,c(genes_subset)]), labels=scale(gene_new2$BCR_READS_D50__ALL), controlscale = TRUE, scale=1)+ggtitle("UMAP")
		plot(ggarrange(x1, x2, ncol=2))
		dev.off()