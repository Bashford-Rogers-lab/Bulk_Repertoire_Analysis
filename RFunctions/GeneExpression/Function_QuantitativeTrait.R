## Function to correlate all genes to a specified AIR Module using LME 
## Currently set for longitudinal analysis 
## Will generate a QTrait volcano plot and export tables for GO enrichment analysis
## Will then do a heatmap for those top interesting genes and see if it idenfities subgroups of patients
## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk


## Recquired Packages!
library(doParallel)
library(ggrepel)
library(gplots)
library(lme4)
library(lmerTest)
library(gghighlight)
library(XGR)
library(gprofiler2)
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

#module_use <- "Module_1"
#gene_expression_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/REPERTOIRE_SAMPLES_GEX_EnsembeID.txt'

get_best_genes <- function(gene_expression_file, eigenvectors_file, outputdir, module_use){
	
		## this will give us ensembl id and gene name 
		gtf <- read.delim('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Homo_sapiens.GRCh38.94_gene_annotation_table.txt', sep="\t")
	    gtf_lookup <- gtf[, c("gene_id", "GeneSymbol")]
		module_use <- module_use

		if(!dir.exists(outputdir)) {dir.create(outputdir)}
		plot_dir <- paste0(outputdir, "ALL_GENES/")
		if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
		###############################################################################
		### Read in files of interest (gene expression, sample cohort assignment etc 
		gene_expression <- (read.delim(gene_expression_file, header=TRUE))
		gene_expression$sample <- NULL
		gene_expression$Barcode<- NULL
		common_col_names <- colnames(gene_expression)
		###########################################################################################
		## Now lets get the module scores and subset to remove the bad ids 
		eigenvectors <- read.delim(eigenvectors_file, sep="\t", header=TRUE)
		bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
		eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
		
		if(module_use %in% colnames(eigenvectors)){
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
			
			print("Number of Overlapping Samples")
			print(dim(gene_exp)[1])
			#####################################################
			print("Getting LME correlations in parallel")
			cl <- makeCluster(30)		
			registerDoParallel(cl)
			##---------------------------------------------------
			## Lets do a lme for all the shared genes
			all_cors <- foreach(i=common_col_names, .combine=rbind, .packages=c("lme4", "lmerTest")) %dopar% {	
					gene <- i
					data_sub <- gene_exp[,c(module_use, gene, "Barcode", "DAY", "DISEASE")]
					
					## LME ANALYSIS: linear mixed effect model also add in batch 
					model_formula <- formula(paste0(module_use, "~", gene, "+ DISEASE + DAY +(1 | Barcode)"))
					lmm_out <<- lmerTest::lmer(model_formula, data = data_sub)
					p_vals <- summary(lmm_out)$coefficients[2,5]
					conf_est <- summary(lmm_out)$coef[2, 1]
					rows <- c(gene, p_vals, conf_est)
					return(rows)
			}
			stopCluster(cl)
			print("Got Effect Estimate for All genes")
		   
			##---------------------------------------------------
			
			all_cors <- data.frame(all_cors)
			colnames(all_cors) <- c("gene", "P", "Estimate")
			all_cors$P <- as.numeric(all_cors$P)
			all_cors$Estimate <- as.numeric(all_cors$Estimate)
			all_cors$BH <- p.adjust(as.numeric(all_cors$P),method="BH")
			
			all_cors <- merge(all_cors, gtf_lookup, by.x="gene", by.y="gene_id")
			## Lets add the gene name
			##-------------------------------------------
			## lets plot the correlation
			## we want to see ones that are adaptive immune reperroite 
			all_cors$Type <- "NON_AIR"
			all_cors$Type[all_cors$GeneSymbol %like% "^IGH" | all_cors$GeneSymbol %like% "^IGL" | all_cors$GeneSymbol %like% "^IGK"| all_cors$GeneSymbol %like% "TRAV"| all_cors$GeneSymbol %like% "TRBV"| all_cors$GeneSymbol %like% "TRGV"| all_cors$GeneSymbol %like% "TRDV"| all_cors$GeneSymbol %like% "TRAD"| all_cors$GeneSymbol %like% "TRBD"| all_cors$GeneSymbol %like% "TRGD"| all_cors$GeneSymbol %like% "TRDD"| all_cors$GeneSymbol %like% "TRAJ"| all_cors$GeneSymbol %like% "TRBJ"| all_cors$GeneSymbol %like% "TRGJ"| all_cors$GeneSymbol %like% "TRDJ"] <- "AIR"
			all_cors$Type <- factor(all_cors$Type, levels=c("NON_AIR", "AIR"))
			
			## Lets get stdeveiation and mean 
			est_std <- sd(all_cors$Estimate)
			est_mean <- mean(all_cors$Estimate)		
			threshold_up <- est_mean+est_std
			threshold_low <- est_mean-est_std
			##############
			
			## Which genes to plot 
			tolabel <- all_cors[all_cors$Estimate >= threshold_up | all_cors$Estimate <= threshold_low,]
			tolabel <- tolabel[tolabel$BH < 0.05,]
			tolabel <- tolabel[order(-tolabel$Estimate),]
			
			tolabel_airr <- tolabel[tolabel$Type=="AIR",]
			tolabel_gene <- tolabel[!tolabel$Type=="AIR",]
			
			############################################
			## Fist lets take those with the greates effect size
			## Do for AIRR genes
			if(dim(tolabel_airr)[1] >= 10){
				tolabel_airr_top <- tolabel_airr[1:5,]
				## if we have strongly negative genes
				if(any(tolabel_airr$Estimate <0)){
					tolabel_airr_low <- tolabel_airr[(dim(tolabel_airr)[1]-4):dim(tolabel_airr)[1],]
				} else {
					tolabel_airr_low <-tolabel_airr_top
				}
			} else {
				## If we have less than ten features in total we will take all of these 
				tolabel_airr_top <- tolabel_airr
				tolabel_airr_low <- tolabel_airr
			}
			
			###############################
			## Do for non AIRR genes 
			if(dim(tolabel_gene)[1] >= 20){
				tolabel_gene_top <- tolabel_gene[1:10,]
				if(any(tolabel_gene$Estimate <0)){
					tolabel_gene_low <- tolabel_gene[(dim(tolabel_gene)[1]-9):dim(tolabel_gene)[1],]
				} else {
				## if we have strongly negative genes
					tolabel_gene_low <-tolabel_gene_top
				}
			} else {
				## If we have less than 20 features in total we will take all of these 
				tolabel_gene_top <- tolabel_gene
				tolabel_gene_low <- tolabel_gene
			}
			
			all_labs <- rbind(tolabel_airr_top, tolabel_airr_low,tolabel_gene_top,tolabel_gene_low )
			all_labs <- unique(all_labs)
			
			#############################################
			### Now lets take those with the smallest p value 
			tolabel <- tolabel[order(tolabel$BH),]
			tolabel_airr <- tolabel[tolabel$Type=="AIR",]
			tolabel_gene <- tolabel[!tolabel$Type=="AIR",]
			tolabel_airr_top_bh <- tolabel_airr[1:5,]
			tolabel_gene_top_bh <- tolabel_gene[1:10,]
			all_labs <- rbind(all_labs, tolabel_airr_top_bh, tolabel_gene_top_bh)
			
			### make sure its all unique 
			all_labs <- unique(all_labs)
			all_labs <- all_labs[complete.cases(all_labs),]
			
			all_cors <- data.frame(all_cors)
			bb <- log10(0.05)
			
			### Lets do a volcano plot of the quantitive gene relationship 
			pdf(paste0(plot_dir, module_use, "Allgenes.pdf"), width=14, height=12)
			p <- ggplot(all_cors, aes(x=log10(BH), colour=Type, y=Estimate)) + geom_point(alpha=0.7) + theme_classic()+ggtitle(paste0("Quantitative Trait Analysis: ", module_use)) +geom_hline(yintercept=threshold_up, col="green")+geom_hline(yintercept=threshold_low, col="green")+geom_vline(xintercept=log10(0.05), col="blue")+geom_label_repel(data=all_labs, aes(label=as.character(GeneSymbol)), seed=10, fill = "white")+ylab(paste0("Effect Estimate: ", module_use))+
			gghighlight::gghighlight(log10(BH) < bb & ( Estimate>threshold_up | Estimate <threshold_low) , use_direct_label = FALSE)+labs(colour="Gene Type") +xlab("Log10(BH Adj P)")
			plot(p)
			dev.off()
			
			table_dir <- paste0(plot_dir, "GeneTables/")
			if (!dir.exists(paste0(table_dir))) {dir.create(table_dir)}

			### now we want to etract those to save 
			write.table(all_cors, paste0(table_dir, "/", module_use, "All_genes_LME.txt"), sep="\t")
			## Plot all those with an effect size greater than one std 
			write.table(tolabel, paste0(table_dir, "/", module_use, "All_genes_LME_significant.txt"), sep="\t")
			
			print("Lets do the Heat Map")	
			#return(all_labs)
			#---------------------------------------
			### Lets see if this pulls out the patient groups 
			#--------------------------------------------------------------------------------------------
			### Lets do a heatmap of the genes of interest 
			#################################################################################################
			#return(all_labs)
			if(any(!is.na(all_labs$gene))& length(unique(!is.na(all_labs$gene)))>10){
				print("Significant Genes Found")
				print("Genes Found for follow up analysis")
				all_labs <- all_labs[complete.cases(all_labs),]
				genes_to_investigate <- gene_exp[,c(all_labs$gene, colnames(eigenvectors), "Barcode")]
				rownames(genes_to_investigate) <- genes_to_investigate$sample
				## now lets do a heatmap for all the genes in health vs sepsis 
				color.map <- function(mol.biol) { if (mol.biol=="SEPSIS") "#FF0000" else "#0000FF" }
				patientcolors <- unlist(lapply(genes_to_investigate$DISEASE, color.map))
				
				## replace ensembl id with gene id 
				for(d in 1:length(colnames(genes_to_investigate))){
					 if(colnames(genes_to_investigate)[d] %in% all_labs$gene){
					 colnames(genes_to_investigate)[d] <- all_labs$GeneSymbol[all_labs$gene==colnames(genes_to_investigate)[d]]
					 }
				}
					
				pdf(paste0(plot_dir, module_use, "Allgenes_Heatmat.pdf"), width=15, height=20)
				heatmap.2(as.matrix(genes_to_investigate[,c(all_labs$GeneSymbol)]), col=topo.colors(100), RowSideColors=patientcolors, cexRow = 0.5, cexCol=0.5, key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="none")
				dev.off()
				
			#-----------------------------------------------------
				pdf(paste0(plot_dir, module_use, "Allgenes_TSNE.pdf"), width=15, height=5)
				x1 <- tsne(t(genes_to_investigate[,c(all_labs$GeneSymbol)]), labels=as.factor(genes_to_investigate$DISEASE))+ggtitle("TSNE")
				x2 <- umap(t(genes_to_investigate[,c(all_labs$GeneSymbol)]), labels=as.factor(genes_to_investigate$DISEASE))+ggtitle("UMAP")
				plot(ggarrange(x1, x2, ncol=2))
				dev.off()
			}
		}
}
