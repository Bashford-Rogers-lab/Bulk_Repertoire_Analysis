assess_batcheffect <- function(cohort1_file, cohort2_file, outputdir, cohort_assignment){

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
	
	## Lets look for batch effect
	gene_expression2 <- t(gene_expression)
	
	pca_uncorrected_obj = prcomp(gene_expression2)
	pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)
	
	## Make sure order of Ids are the same so we can use this for plotting 
	compare_ids <- colnames(gene_expression2)==rownames(pca_uncorrected)
	print(paste0("Order of Ids is the same: ", all(compare_ids == TRUE)))
	###########################
	
	## Assign sources of batch and condition
	pca_uncorrected$condition = "S"
	pca_uncorrected$condition[rownames(pca_uncorrected) %like% "HV"] = "H"
	pca_uncorrected$Cohort[rownames(pca_uncorrected) %in% cohort1_samples$SampleID] <- "C1"
	pca_uncorrected$Cohort[rownames(pca_uncorrected) %in% cohort2_samples$SampleID] <- "C2"
	pca_uncorrected$DISEASE2 <- paste0(pca_uncorrected$condition, ".",  pca_uncorrected$Cohort)
	
	pdf(paste0(outputdir, "BATCHEFFECT_PCA.pdf"), width=20, height=20)
	p1x <- (M3C::pca(gene_expression2, pcx=1, pcy=2,  labels=as.factor(pca_uncorrected$DISEASE2))+ggtitle("PCA")) 
	p2x <- (M3C::pca(gene_expression2, pcx=2, pcy=3,  labels=as.factor(pca_uncorrected$DISEASE2))+ggtitle("PCA")) 
	p3x <- (M3C::pca(gene_expression2, pcx=3, pcy=4,  labels=as.factor(pca_uncorrected$DISEASE2))+ggtitle("PCA")) 
	p4x <- (M3C::pca(gene_expression2, pcx=4, pcy=5,  labels=as.factor(pca_uncorrected$DISEASE2))+ggtitle("PCA")) 
	plot(plot_grid(p1x, p2x, p3x, p4x, ncol=2))
	p1 = ggplot(data=pca_uncorrected, aes(x=PC1, y=PC2, color=DISEASE2, shape=Cohort))+geom_point()+theme_classic()+labs(color="DISEASE STATUS")+guides(color="none", shape="none")
	p2 = ggplot(data=pca_uncorrected, aes(x=PC2, y=PC3, color=DISEASE2, shape=Cohort))+geom_point()+theme_classic()+labs(color="DISEASE STATUS")+guides(color="none", shape="none")
	p3 = ggplot(data=pca_uncorrected, aes(x=PC3, y=PC4, color=DISEASE2, shape=Cohort))+geom_point()+theme_classic()+labs(color="DISEASE STATUS")+guides(color="none", shape="none")
	p4 = ggplot(data=pca_uncorrected, aes(x=PC4, y=PC5, color=DISEASE2, shape=Cohort))+geom_point()+theme_classic()+labs(color="DISEASE STATUS")
	plot(plot_grid(p1, p2, p3, p4, ncol=2))
	dev.off()
	
	pdf(paste0(outputdir, "BATCHEFFECT_ALL.pdf"), width=25, height=15)
	x1 <- tsne(gene_expression2, labels=as.factor(pca_uncorrected$DISEASE2))+ggtitle("TSNE")
	x2 <- umap(gene_expression2, labels=as.factor(pca_uncorrected$DISEASE2))+ggtitle("UMAP")
	plot(ggarrange(ggarrange(x1, x2, ncol=2), ggarrange(p1x, p2x, p3x, p4x, ncol=4), nrow=2))
	#plot(plot_grid(x1, x2, ncol=2))
	dev.off()
}
	