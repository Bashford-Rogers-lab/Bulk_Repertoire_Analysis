## Gene Ontology analsyis using XGR and gprofileR 
## Use genes annotated as non-AIR 
## Identify significant pathways and plot using REACTOME ANNOTATIONS

library(gghighlight)
library(XGR)
library(gprofiler2)
library(data.table)
library(stringr)

#gene_list <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/GeneExpPlots/ALL_GENES/GeneTables/Module_29All_genes_LME_significant.txt'
#outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/GeneExpPlots/'
 

gene_ontology_analysis <- function(gene_list, outputdir){
		
		## Set up plotting directory
		plot_dir <- paste0(outputdir, "ALL_GENES/")
		if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
		table_dir <- paste0(outputdir, "GO_tables/")
		if (!dir.exists(table_dir)) {dir.create(table_dir)}
		
		
		## read in the genelist 
		genes <- read.delim(gene_list, header=TRUE, sep="\t")
		
		## Get the background aka genes in our dataset 
		all_genes <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/GeneExpPlots/ALL_GENES/GeneTables/Module_31All_genes_LME.txt'
		all_genes <- read.delim(all_genes, header=TRUE, sep="\t")
		all_genesx <- all_genes$gene
		all_genesx2 <- all_genes$GeneSymbol
		
		### lets get rid of the non airr genes
		genes <- genes[!genes$Type =="AIR",]
		
		module_use <- unlist(str_split(basename(gene_list), "All_genes_LME_significant.txt"))[1]
		print(module_use)
		
		### Location of the gene look up 
		gtf <- read.delim('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Homo_sapiens.GRCh38.94_gene_annotation_table.txt', sep="\t")
	    gtf_lookup <- gtf[, c("gene_id", "GeneSymbol")]
		
		
		### Only run if we picked yup significant genes
		if(dim(genes)[1]>0){
	
			###--------------------------------------------------------------
			###--------------------------------------------------------------
			## PART 1: LETS DO IT WITH grpofileR
			gostres <- gost(query = genes$gene, 
					organism = "hsapiens", ordered_query = FALSE, 
					multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
					measure_underrepresentation = FALSE, evcodes = TRUE, 
					user_threshold = 0.05, correction_method = "g_SCS", 
					domain_scope = "annotated", custom_bg = NULL, 
					numeric_ns = "", sources = c("REAC"), as_short_link = FALSE)
					
			
			
			## results 
			if(!is.null(gostres)){
				gpro_results <- gostres$result
				pdf(paste0(plot_dir, module_use, "_GO_ALLPathways.pdf"), width=7, height=5)
				plot(p <- gostplot(gostres, capped = TRUE, interactive = FALSE))
				dev.off()
				
				
				## Lets get the ids and match it to the genes!!!
				gpro_results$AIRR_intercept <- NA
				gpro_results$GENE_intercept <- NA
				gpro_results$AIRR_intercept_count <- NA
				gpro_results$GENE_intercept_count <- NA
				
				for(i in 1:dim(gpro_results)[1]){
					genesx <- gpro_results$intersection[i]
					genesx <- str_split(genesx, ",")
					genesx <- unlist(genesx)
					genesx <- gtf_lookup[gtf_lookup$gene_id %in% genesx,]
					## Lets classify by type 
					airr <- genesx[genesx$GeneSymbol %like% "^IGH" | genesx$GeneSymbol %like% "^IGL" | genesx$GeneSymbol %like% "^IGK",]
					non_airr <- genesx[!genesx$GeneSymbol %like% "^IGH" & !genesx$GeneSymbol %like% "^IGL" & !genesx$GeneSymbol %like% "^IGK",]
					gpro_results$AIRR_intercept_count[i] <- length(airr[,1])
					gpro_results$GENE_intercept_count[i] <- length(non_airr[,1])
					### Concatenate 
					airr <- paste(airr$GeneSymbol, collapse = ",")
					non_airr <- paste(non_airr$GeneSymbol, collapse = ",")
					## Add to dataframe 
					gpro_results$AIRR_intercept[i] <- airr
					gpro_results$GENE_intercept[i] <- non_airr
			
				}
				gpro_results$intersection <- NULL
				
				## Lets pull out terms with genes!!!
				gene_terms <- gpro_results$term_id[gpro_results$GENE_intercept_count>=5]
				gene_terms2 <- gpro_results$term_id[!gpro_results$term_id %in% gene_terms]

				pdf(paste0(plot_dir, module_use, "_GO_SPLIT_Pathways.pdf"), width=12, height=17)
				if(length(gene_terms)>=1){
					plot(publish_gostplot(p, highlight_terms = c(gene_terms))+ggtitle(paste0(module_use, " Pathways including >5 non-AIRR Genes")))
				}
				if(length(gene_terms2)>=1){
					plot(publish_gostplot(p, highlight_terms = c(gene_terms2))+ggtitle(paste0(module_use, " Pathways including <5 non-AIRR Genes")))
				}
				dev.off()
				
				gpro_results[,15] <- NULL
				gpro_results[,10] <- NULL
				## Flatten this column which is a list
				gpro_results <- gpro_results %>% dplyr::rowwise() %>% dplyr::mutate(parents = paste(parents, collapse=',')) %>% dplyr::ungroup()
				write.table(gpro_results, paste0(table_dir, module_use, "_GO.txt"), sep="\t", row.names=FALSE)				
				print("GRPROFILER PATHWAYS FOUND")

			} else {
				print("GPROFILER No Pathway Found")
			}
			
			
			###--------------------------------------------------------------
			###--------------------------------------------------------------
			## PART2: LETS DO IT WITH XGR
			#"GOBP", "GOCC", "GOMF"

			myxgr <- xEnricherGenes(data = genes$GeneSymbol, ontology = c("REACTOME"), background=all_genesx2, p.adjust.method="BH", verbose=FALSE)
			#myxgr1 <- xEnricherGenes(data = genes$GeneSymbol, ontology = c("GOBP"), background=all_genesx2, p.adjust.method="BH", verbose=FALSE)
			#myxgr2 <- xEnricherGenes(data = genes$GeneSymbol, ontology = c("GOMF"), background=all_genesx2, p.adjust.method="BH", verbose=FALSE)
			#myxgr3 <- xEnricherGenes(data = genes$GeneSymbol, ontology = c("GOCC"), background=all_genesx2, p.adjust.method="BH", verbose=FALSE)

			if(!is.null(myxgr)){
				res <- xEnrichViewer(myxgr, top_num=length(myxgr$adjp), sortBy="adjp", details=TRUE)
				ressig <- res[res$pvalue <= 0.05,]
				if(dim(ressig)[1] != 0 ){
					ressig$adjp <- as.numeric(ressig$adjp)
					ressig <- ressig[order(ressig$adjp),]
					ressig$name <- str_wrap(ressig$name, width = 45)
					ressig$name <- factor(ressig$name, levels=c(ressig$name))

					if(dim(ressig)[1]>30){
						ressig2 <- ressig[1:30,]
					} else {
						ressig2 <- ressig
					}	
					
					ressig2 <- ressig2[order(-ressig2$adjp),]
					ressig2$name <- factor(ressig2$name, levels=c(ressig2$name))

					if(dim(ressig2)[1]>10){
						dimx = 8
						dimy = 7
					} else {
						dimx=7
						dim7 = 17
					}
					
					## Plot results
					pdf(paste0(plot_dir, module_use, "_XGR_Pathways.pdf"), width=dimx, height=dimx)
					## The inbuilt function doesnt work very well
					#plot(xEnrichBarplot(myxgr, top_num="auto", displayBy="adjp"))
					plot(ggplot(ressig2, aes(x=(-1*log10(adjp)), y=name, fill=log10(adjp))) +geom_col() + theme_classic()+xlab("-log10(Adjp)")+guides(fill="none")+ggtitle(paste0(module_use, ", p<0.05"))+geom_vline(xintercept=-log10(0.05), col="green") +ylab("REACTOME Term"))
					dev.off()
						
					## Look at genes it couldnt use 
					genes_used <- myxgr$data
					##notused - seems to be non cording rnas/pseudogenes and transcripts 
					genes$GeneSymbol[!genes$GeneSymbol %in% genes_used]

					## save xgr results 
					write.table(res, paste0(table_dir, module_use, "_XGR.txt"), sep="\t", row.names=FALSE)
					print("XGR PATHWAYS FOUND")
				}
			}  else {
				print("XGR No Pathway Found")
			}
		} else {
			print("No Signficant genes")
		}
	}

			   