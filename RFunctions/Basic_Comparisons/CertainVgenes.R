imp_data_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_BCR_PRODUCTIVE.txt'
genes <- c("IGHV4-34", "IGHV4-28", "IGHV4-55", "IGHV4-31")
genes2 <- c("IGHV4_34", "IGHV4_28", "IGHV4_55", "IGHV4_31")


#### Gene expression
cohort1_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/logcpm_181_20600_cohort2.txt'
cohort2_file <- "/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/Logcpm_864_20416.txt"
cohort_assignment <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/RepertoireRNAseqCohorts.txt'
eigenvectors_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'



get_vgene_usage <- function(file=file, chain_vdj=chain_vdj, ids_all=ids_all, counts_used=counts_used, iso_type=iso_type, genes){
	###################
    datax <- read.delim(imp_data_file)
	keeps <- c()
	for(i in 1:length(genes2)){
		keep <- colnames(datax)[colnames(datax) %like% genes2[i]]
		keeps <- c(keep, keeps)
	}
	datax2 <- datax[c(keeps)]
	datax2$Sample <- rownames(datax2)
	datax2 <- gather(datax2, "V.gene", "percent_repertoire", -Sample)
	datax2$DISEASE="SEPSIS"
	datax2$DISEASE[datax2$Sample  %like% "HV"] <- "HEALTH"
	datax2$DAY <- NA
	datax2$DAY[grep("_1", datax2$Sample )] <- "1"
	datax2$DAY[grep("_3", datax2$Sample )] <- "3"
	datax2$DAY[grep("_5", datax2$Sample )] <- "5"
	datax2$Sample <- gsub("_productive", "", datax2$Sample)
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	datax2 <- datax2[!datax2$Sample %in% c(bad_ids),]
	datax2 <- datax2[!datax2$Sample %like% "JR1795_1003_POSITIVE",]
	
	pdf(paste0(outputdir, "/VGENE2.pdf"), width=10, height=10)
	ggplot(datax2, aes(x=DAY, y=percent_repertoire, fill=DISEASE)) + geom_boxplot() +facet_wrap(~V.gene, scales="free")
	p2 <- ggline(datax2, x = "DAY", y = "percent_repertoire", add = c("mean_ci", "point"), color = "DISEASE", palette = c("#00AFBB", "#E7B800"), ylab="Percentage", facet.by="V.gene", main="Unique Reads")+ stat_compare_means(aes(group = DISEASE), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p3 <- ggline(datax2, x = "DAY", y = "percent_repertoire", add = c("mean_ci", "point"), color = "DISEASE", palette = c("#00AFBB", "#E7B800"), ylab="Percentage", facet.by="V.gene", main="Total Reads")+ stat_compare_means(aes(group = DISEASE), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	plot_grid(p2, p3, ncol = 2, align="hv", axis="tblr")
	dev.off()
	
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
	
	gtf <- read.delim('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Homo_sapiens.GRCh38.94_gene_annotation_table.txt', sep="\t")
	gtf_lookup <- gtf[, c("gene_id", "GeneSymbol")]

	genes_keep <- gtf_lookup[gtf_lookup$GeneSymbol %in% genes,]
	exp_v <- 	gene_exp[, c(colnames(gene_exp)[colnames(gene_exp)%in% genes_keep$gene_id], "sample", "DAY", "DISEASE", "Barcode")]
	exp_v2 <- gather(exp_v, "V.gene", "lcp", -c(sample,DAY,DISEASE,Barcode))
	exp_v2 <- merge(exp_v2, genes_keep, by.x="V.gene", by.y="gene_id")
	
	## merge gene expression and repertoire counts
	datax2$V.gene <- gsub("BCR_READS_", "", datax2$V.gene)
	datax2$V.gene <- gsub("__ALL", "", datax2$V.gene)
	datax2$V.gene <- gsub("_", "-", datax2$V.gene)
	
	colnames(exp_v2)[7] <- "V.gene"
	colnames(exp_v2)[2] <- "Sample"
	exp_v2[,1] <- NULL
	final <- merge(datax2, exp_v2, by=c("Sample", "V.gene"))
	
	pdf(paste0(outputdir, "/VGENE21.pdf"), width=7, height=7)
	ggplot(final, aes(x=log2(percent_repertoire+1), y=lcp)) + geom_point(aes(colour=DISEASE.x, shape=DAY.x)) +facet_wrap(~V.gene, scales="free")+ geom_smooth(method='lm', formula= y~x)+
    stat_cor(method = "pearson") +theme_classic() +ylab("Gene Exp: log2(CPM=1)") +xlab("Repertoire: log2(Percent Repertoire+1)")+labs(colour="Disease", shape="Day")
	dev.off()
	
	
)

