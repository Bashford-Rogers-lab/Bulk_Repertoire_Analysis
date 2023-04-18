### Code to correlate gene and protein expression
## Lauren.overend@oriel.ox.ac.uk
#module load R/4.1.0-foss-2021a  
library(ggplotify)
library(grid)
library(gridExtra)
library(ggbiplot)
library(M3C)
## 
gex <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/REPERTOIRE_SAMPLES_GEX_GeneSymbol.txt', sep="\t")
protein <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/Yuxin_PROTEOMIC_PROCESSED/REPERTOIRE_SAMPLES.txt', sep="\t")
protein$sample <- rownames(protein)
### Look at immunoglobulin
genes_interest <- c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2")
protein <- protein[c(genes_interest)]
rownames(protein)<- gsub("GAUKJR010000", "UKJR010000", rownames(protein))

#############################################
####  FIRST WE CORRELATE PROTEINS AND RNASEQ
my_plots <- c()
for(i in 1:length(genes_interest)){
	 gene <- genes_interest[i]
	 gex_gene <- gex[, c("sample", gene)]
	 colnames(gex_gene) <- c("sample", paste0(gene, "_gex"))
	 prot_gene <- protein[, c(gene), drop=FALSE]
	 prot_gene$sample <- rownames(prot_gene)
	 colnames(prot_gene) <- c(paste0(gene, "_protein"), "sample")
	 datax <- merge(gex_gene, prot_gene, by="sample")

	 ## hmm whats missing 
	 unique_ids <- length(unique(c(gex_gene$sample,prot_gene$sample))) 
	 ## 127 in common!!!
	 ## Lets assign barcode 
	 datax$Barcode <- NA
		for(x in 1:length(datax$sample)){
			datax$Barcode[x] <- str_split_fixed (datax$sample[x], "_", 2)[,1]
			if(datax$sample[x] %like% "HV"){
			datax$Barcode[x] <- paste0(str_split_fixed(datax$sample[x], "_", 3)[,1], "_", str_split_fixed(datax$sample[x], "_", 3)[,2])
			}
		}		
	
	colnames(datax) <- c("sample", "gex", "protein", "Barcode")

	##RMCORR ANALYSIS::::   x then y 
	rmcorr_out <- rmcorr(Barcode, "gex", "protein", datax)
	pvalue <- rmcorr_out$p 
	plt_label <- paste0(gene, ": R = ", round(rmcorr_out$r, 3), ", p = ", round(pvalue, 3),"\n95% CI: [", paste0(round(rmcorr_out$CI, 3), collapse = ", "), "]")
	par(mar = rep(2, 4))
	p1 <<- as.ggplot(~plot(x=rmcorr_out, xlab="GEX log2(CPM+1)", ylab="log2(Protein Abundance)", main = plt_label,overall.col="black"))
	my_plots[[i]] <- p1 
}

pdf('/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/Comparisons/IsotypeGEXvsProtein.pdf', width=16, height=16)
do.call("grid.arrange", c(my_plots, ncol=ceiling(sqrt(length(genes_interest)[1]))))
dev.off()

### we can then model this!!!
protein$sample <- rownames(protein)
write.table(protein, '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/Yuxin_PROTEOMIC_PROCESSED/REPERTOIRE_SAMPLES_IMMUNOGLOBULIN.txt', sep="\t")

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
####  REPERTOIRE AND RNASEQ 
imp_data_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_BCR_PRODUCTIVE.txt'
imp_data <- read.delim(imp_data_file, header=TRUE)
eigenvectors <- read.delim('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt')

#### Repertoire DATA 
## First lets look at isotype usage 
iso_usage_cols <- colnames(imp_data)[colnames(imp_data) %like% "_vdjs_per_isotype"]
iso_usage <- imp_data[, c(iso_usage_cols)]
colnames(iso_usage) <- gsub("BCR_READS_percentage_", "", colnames(iso_usage))
colnames(iso_usage) <- gsub("_vdjs_per_isotype", "", colnames(iso_usage))
iso_usage$sample <- rownames(iso_usage)
iso_usage$sample <- gsub("_productive", "", iso_usage$sample )

csw_uniquex <- colnames(iso_usage)[colnames(iso_usage) %like% "unique"]
csw_unique <- csw_uniquex[!csw_uniquex %like% "IGHM" & !csw_uniquex %like% "IGHD" ]
nsw_unique <- csw_uniquex[csw_uniquex %like% "IGHM" |csw_uniquex %like% "IGHD" ]

## need to take out bad ids
bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
iso_usage <- iso_usage[!iso_usage$sample %like% "JR1795_1003_POSITIVE",]
iso_usage <- iso_usage[!iso_usage$sample %in% c(bad_ids),]

## RNA seq data
genes_interest <- c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE")
gex <- gex[c(genes_interest, "sample")]
gex_norm <- gex

## Lets inverse log first 
gex_norm$sample <- NULL
### lets inverse log (its reported as log2(CPM+1)
gex_norm2<- 2^gex_norm
gex_norm2 <- gex_norm2-1
gex_norm_nonlogged <- gex_norm2
gex_norm_nonlogged$sample <- rownames(gex_norm_nonlogged)
###------------------------------------------------
###------------------------------------------------
###------------------------------------------------
###------------------------------------------------
## First lets do a correlation showing that percentage measure and non-logged counts dont well correlate!!!!
genes_interest <- c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE")  #, "CSW", "NCSW"

my_plots <- c()
for(i in 1:length(genes_interest)){
	## from gene expressoion
	gene <- genes_interest[i]
	gex_gene <- gex_norm_nonlogged[, c("sample", gene)]
	colnames(gex_gene) <- c("sample", paste0(gene, "_gex"))
	
	## from repertoire
	rep_genex <- colnames(iso_usage)[colnames(iso_usage) %like% gene]
	rep_genex <-rep_genex[rep_genex %like% "total"]	
 	rep_gene <- iso_usage[, c(rep_genex), drop=FALSE]
	rep_gene$sample <- rownames(rep_gene)
	rep_gene$sample <- gsub("_productive", "", rep_gene$sample )
	
	## merge together (for BCR should be 180!!)
	datax <- merge(gex_gene, rep_gene, by="sample")
	print(dim(datax)[1])
	
	## Lets assign barcode 
	datax$Barcode <- NA
		for(x in 1:length(datax$sample)){
			datax$Barcode[x] <- str_split_fixed (datax$sample[x], "_", 2)[,1]
			if(datax$sample[x] %like% "HV"){
			datax$Barcode[x] <- paste0(str_split_fixed(datax$sample[x], "_", 3)[,1], "_", str_split_fixed(datax$sample[x], "_", 3)[,2])
			}
		}		
	colnames(datax) <- c("sample", "gex", "rep", "Barcode")

	##RMCORR ANALYSIS::::   x then y 
	rmcorr_out <- rmcorr(Barcode, "gex", "rep", datax)
	pvalue <- rmcorr_out$p 
	plt_label <- paste0(gene, ": r = ", round(rmcorr_out$r, 3), ", p = ", round(pvalue, 3),"\n95% CI: [", paste0(round(rmcorr_out$CI, 3), collapse = ", "), "]")
	par(mar = rep(2, 4))
	p1 <<- as.ggplot(~plot(x=rmcorr_out, xlab="CPM: Bulk RNAseq GEX", ylab="% of Repertoire (total): BCR Repertoire", main = plt_label,overall.col="black"))
	my_plots[[i]] <- p1 
}
pdf('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/ModalityComparisons/IsotypeGEX_nonnormalised_vs_REP.pdf', width=15, height=17)
do.call("grid.arrange", c(my_plots, ncol=ceiling(sqrt(length(genes_interest)[1]))))
dev.off() 

###------------------------------------------------
## PART 2 normalise gex as a proportion the same was as iso usage is done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
gex_norm2<-gex_norm2[, c(1:9)]

## normalise 
gex_norm2<-t(apply(gex_norm2,1, function(x) x/sum(x)))
gex_norm2 <- gex_norm2*100
gex_norm2 <- data.frame(gex_norm2)
gex_norm2$sample <- rownames(gex_norm2)

### Lets make the same comparison and show they correlate 
genes_interest <- c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE")  #, "CSW", "NCSW"
my_plots2 <- c()
for(i in 1:length(genes_interest)){
	## from gene expressoion
	gene <- genes_interest[i]
	gex_gene <- gex_norm2[, c("sample", gene)]
	colnames(gex_gene) <- c("sample", paste0(gene, "_gex"))
	
	## from repertoire
	rep_genex <- colnames(iso_usage)[colnames(iso_usage) %like% gene]
	rep_genex <-rep_genex[rep_genex %like% "total"]	
 	rep_gene <- iso_usage[, c(rep_genex), drop=FALSE]
	rep_gene$sample <- rownames(rep_gene)
	rep_gene$sample <- gsub("_productive", "", rep_gene$sample )
	
	## merge together (for BCR should be 180!!)
	datax <- merge(gex_gene, rep_gene, by="sample")
	print(dim(datax)[1])
	
	## Lets assign barcode 
	datax$Barcode <- NA
		for(x in 1:length(datax$sample)){
			datax$Barcode[x] <- str_split_fixed (datax$sample[x], "_", 2)[,1]
			if(datax$sample[x] %like% "HV"){
			datax$Barcode[x] <- paste0(str_split_fixed(datax$sample[x], "_", 3)[,1], "_", str_split_fixed(datax$sample[x], "_", 3)[,2])
			}
		}		
	colnames(datax) <- c("sample", "gex", "rep", "Barcode")

	##RMCORR ANALYSIS::::   x then y 
	rmcorr_out <- rmcorr(Barcode, "gex", "rep", datax)
	pvalue <- rmcorr_out$p 
	plt_label <- paste0(gene, ": r = ", round(rmcorr_out$r, 3), ", p = ", round(pvalue, 3),"\n95% CI: [", paste0(round(rmcorr_out$CI, 3), collapse = ", "), "]")
	par(mar = rep(2, 4))
	p1 <<- as.ggplot(~plot(x=rmcorr_out, xlab=" % CPM: Bulk RNAseq GEX", ylab="% of Repertoire (total): BCR Repertoire", main = plt_label,overall.col="black"))
	my_plots2[[i]] <- p1 
}

pdf('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/ModalityComparisons/IsotypeGEX_normalised_vs_REP.pdf', width=15, height=17)
do.call("grid.arrange", c(my_plots2, ncol=ceiling(sqrt(length(genes_interest)[1]))))
dev.off() 

### Look the correlation is so much better!!!!!!!
###------------------------------------------------
###------------------------------------------------
## Lets proceed with modeling the logCPM and how they change over time 
write.table(gex, paste0("/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/Isotype_lcpm.txt"), sep="\t")

## Lets plot how expression changes over time 
gex <- "/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/Isotype_lcpm.txt" 
eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/ModalityComparisons/'
metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'

source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Model_BULK_genecounts_hvs.R')
correlate_genes(eigenvectors, gex, outputdir, metadata, metahealth)


###------------------------------------------------
###------------------------------------------------




### Can we use the cyber sort code to model 

#################
gex_norm2 <- data.frame(gex_norm2)
gex_norm2$sample <-rownames(gex_norm2)
gex_norm2$DAY <- NA
gex_norm2$DAY[grep("_1", gex_norm2$sample)] <- "1"
gex_norm2$DAY[grep("_3", gex_norm2$sample)] <- "3"
gex_norm2$DAY[grep("_5", gex_norm2$sample)] <- "5"

### calculate proportion of reads 
gex_norm3 <- gather(gex_norm2, "isotype", "proportion_reads", -c(sample, DAY))
gex_norm3$DISEASE <- "SEPSIS"
gex_norm3$DISEASE[gex_norm2$sample %like% "HV"] <- "HEALTH"

## repertoire
keeps <- colnames(iso_usage)[colnames(iso_usage) %like% "unique"]
iso_2 <- iso_usage[, c(keeps, "sample")]
iso_2 <- gather(iso_2, "isotype", "proportion_repertoire", -c(sample))
iso_2$DAY <- NA
iso_2$DAY[grep("_1", iso_2$sample)] <- "1"
iso_2$DAY[grep("_3", iso_2$sample)] <- "3"
iso_2$DAY[grep("_5", iso_2$sample)] <- "5"
iso_2$DISEASE <- "SEPSIS"
iso_2$DISEASE[iso_2$sample %like% "HV"] <- "HEALTH"
iso_2$isotype <- gsub("unique__", "", iso_2$isotype)

## non normalised bulk 
gex2 <-  gather(gex, "isotype", "lcpm", -c(sample))
gex2 <- gex2[!gex2$isotype %like% "CSW",]
gex2$DAY <- NA
gex2$DAY[grep("_1", gex2$sample)] <- "1"
gex2$DAY[grep("_3", gex2$sample)] <- "3"
gex2$DAY[grep("_5", gex2$sample)] <- "5"
gex2$DISEASE <- "SEPSIS"
gex2$DISEASE[gex2$sample %like% "HV"] <- "HEALTH"
gex2$isotype <- gsub("unique__", "", gex2$isotype)


pdf('/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/Comparisons/GEXnormalised_isotype.pdf', width=13, height=10)
## First do a plot of gex normalised counts 
p1 <- ggline(gex_norm3, x = "DAY", y = "proportion_reads", add = "mean_ci", color = "isotype", ylab="Percentage of Repertoire", facet.by="DISEASE", main="Bulk RNAseq")
p2 <- ggline(iso_2, x = "DAY", y = "proportion_repertoire", add = "mean_ci", color = "isotype", ylab="Percentage of Repertoire", facet.by="DISEASE", main="BCR Repertoire")
p3 <- ggline(gex2, x = "DAY", y = "lcpm", add = "mean_ci", color = "isotype", ylab="log2(CPM)", facet.by="DISEASE", main="Bulk RNAseq: Expression")
plot_grid(p1, p2, p3, ncol=2)
dev.off()
pdf('/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/Comparisons/GEXnormalised_isotype_HvsD.pdf', width=8, height=8)
p3 <- ggline(gex2, x = "DAY", y = "lcpm", add = "mean_ci", color = "DISEASE", ylab="log2(CPM+1)", palette = c("#00AFBB", "#E7B800"),facet.by="isotype", main="Bulk RNAseq")+ stat_compare_means(aes(group = DISEASE), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
plot(p3)
dev.off()




