
correlate_target_genes <- function(cohort1_file, cohort2_file, eigenvectors_file, gene_df, outputdir, cohort_assignment, module_list, list_name){
	
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

		#######################
		### Genes to use 
		genes_to_use <- gene_df 
		list_name <- list_name

		#######################################################
		### Lets subset the gene expression by the genes of interest 
		## genes in outputdir
		### Make sure we are only using genes that are present in the dataframe 
		genes_keep <- colnames(gene_expression)[colnames(gene_expression) %in% genes_to_use[,1]]
		## subset genes we are using 
		genes_to_use <- genes_to_use[genes_to_use$V1 %in% genes_keep, ]
		
		### Lets subset for genes of interest 
		genex <- gene_expression[,c(genes_to_use$V1)]
		genex <- data.frame(genex)
		
		## Lets use the other gene name (make sure in right order)
		colnames(genex) <- c(genes_to_use$V2)
		genex$sample <- rownames(genex)
		### convert to long format 
		genex_long <-  gather(genex, gene, expression, -sample, factor_key=TRUE)

		#### Merge with the eigenvectors..... 
		genex_long <- merge(genex_long, eigenvectors, by="sample")
		### There will be some samples missing due to missing rnaseq/repertoire data 
		genex_long$Cohort[genex_long$sample %in% cohort1_samples$SampleID] <- "Cohort1"
		genex_long$Cohort[genex_long$sample %in% cohort2_samples$SampleID] <- "Cohort2"
		## Assign a barcode ID 
		genex_long$Barcode <- NA
				for(x in 1:length(genex_long$sample)){
					genex_long$Barcode[x] <- str_split_fixed (genex_long$sample[x], "_", 2)[,1]
					if(genex_long$sample[x] %like% "HV"){
					 genex_long$Barcode[x] <- paste0(str_split_fixed(genex_long$sample[x], "_", 3)[,1], "_", str_split_fixed(genex_long$sample[x], "_", 3)[,2])
					}
				}
		genex_long$expression <- as.numeric(genex_long$expression)
		genex_long$DISEASE2 <- paste0(genex_long$DISEASE, "_",  genex_long$Cohort)
		
		## So that the plots are plotted the same way 
		order_genes <- levels(genex_long$gene)
		
		for(k in 1:length(module_list)){
			module_use <- module_list[k]
			
			## set dimensions 
			dimens <- dim(genes_to_use)[1]*0.6
			if(dimens < 10) {
				dimens <- 10
			}
			### First lets visualise 
			pdf(paste0(outputdir, list_name, "_", module_use, ".pdf"), width=dimens, height=dimens)
			plot(ggplot(genex_long, aes(x=get(module_use), y=expression))  +geom_line(aes(group=Barcode),col="lightgrey")+geom_point(aes(colour=DISEASE2, shape=DAY))+theme_classic() +facet_wrap(~gene,  scales = "free_y") + geom_smooth(method='lm', formula= y~x)+ggpubr::stat_cor(method = "pearson")+xlab(module_use) +ylab("lcm"))
			dev.off()
			
			my_comparisons <- list( c("SEPSIS_Cohort1", "SEPSIS_Cohort2") )
			pdf(paste0(outputdir, list_name, "_", module_use, "_BOXPLOT.pdf"), width=dimens, height=dimens)
			plot(ggplot(genex_long, aes(x=DISEASE2, y=expression)) +geom_boxplot(aes(fill=DAY))+facet_wrap(~gene, scales="free")+ theme_classic() + scale_x_discrete(guide = guide_axis(n.dodge=2)) + stat_compare_means(comparisons = my_comparisons)+
			scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))
			plot(ggplot(genex_long, aes(x=DISEASE2, y=expression)) +geom_boxplot(aes())+facet_wrap(~gene, scales="free")+ theme_classic() + scale_x_discrete(guide = guide_axis(n.dodge=2)) + stat_compare_means(comparisons = my_comparisons)+
			scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))
			dev.off()
			
			#facet_grid(rows=vars(gene), cols=vars(DAY), scales="free")
			
			## Lets get the correlations
			rm_corr_list <- list()
			lme_list <- list()
			results_stats <- c()
			##https://martakarass.github.io/post/2022-04-27-rmcorr_vs_lmm/ sourced from here to visualise mixed effects vs rmcorr 
			
			plot_dir <- paste0(outputdir, list_name)
			if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
			
			##-------------------------------------------------------------------------------
			## Now to test relationship between module and genes 
			for(i in 1:length(order_genes)){
				data_sub <- genex_long[genex_long$gene ==order_genes[i],]
				data_sub <- data_sub[,c(module_use, "Barcode", "DAY", "sample", "DISEASE", "Cohort", "gene", "expression")]
				gene_use <<- order_genes[i]
				colnames(data_sub)[1] <- "Module_Score"
				
				pdf(paste0(plot_dir, "/", list_name, "_", module_use, "_", gene_use, "_LME.pdf"), width=10, height=7)
				##RMCORR ANALYSIS::::   x then y 
				rmcorr_out <<- rmcorr(Barcode, "Module_Score", "expression", data_sub)
				pvalue <<- rmcorr_out$p 
				plt_label <<- paste0("r = ", round(rmcorr_out$r, 3), ", p = ", round(pvalue, 3),"\n95% CI: [", paste0(round(rmcorr_out$CI, 3), collapse = ", "), "]")
				par(mar = rep(2, 4))
				p1 <<- as.ggplot(~plot(x=rmcorr_out, overall=TRUE, xlab=module_use, ylab=paste0(gene_use, " lcm"), main = plt_label,overall.col="black"))
				rm_corr_list[[i]] <- p1 
				
				## Standard correlation
				std_corr <- cor.test(data_sub$Module_Score, data_sub$expression, method=c("pearson"))
			
				## LME ANALYSIS: linear mixed effect model also add in batch 
				colnames(data_sub)[1] <- "x"
				colnames(data_sub)[8] <- "y"
				lmm_out <<- lmerTest::lmer(y ~ x  + (1 | Barcode), data = data_sub)
				p_vals <- summary(lmm_out)$coefficients[,5]
				p_vals <- p_vals[2]
				conf_est <- summary(lmm_out)$coef[2, 1]
				conf_out <- as.vector(confint(lmm_out, parm = "x"))
				plt_label <- paste0(gene_use, "\nbeta = ", round(conf_est, 3), ", p = ", round(p_vals, 3),"\n95% CI: [", paste0(round(conf_out, 3), collapse = ", "), "]")
				data_sub$mu <- getME(lmm_out, "mu")[1 : nrow(data_sub)]
				plt <- ggplot(data_sub, aes(x = x, y = y, color = Barcode, group = Barcode)) +  geom_line(aes(x = x, y = mu), size = 0.3) + geom_point() + theme_classic(base_size = 12) + theme(legend.position = "none",)  + xlab(module_use) + ylab(paste0("lcm"))+labs(title = plt_label)
				lme_list[[i]] <- plt
			
				##Plot rmcorr and simple effects next to each other 
				plot(plot_grid(p1, plt, ncol=2, labels = 'AUTO', hjust = 0, vjust = 1, align="h", axis = "bt"))
				dev.off()
				
				## need to save values to a table 
				results <- c(module_use,as.character(gene_use), rmcorr_out$r, pvalue, conf_est, p_vals, std_corr$estimate, std_corr$p.value) 
				results <- data.frame(t(results))
				colnames(results) <- c("Module", "Gene", "RMCORR_R", "RMCORR_P", "LME_EST", "LME_P", "CORR_R", "CORR_P")
				results_stats <- rbind(results, results_stats)
				
				}
				
				##--------------------------------------------------------------------
				## lets save the stats 
				write.table(results_stats, paste0(outputdir,"/", list_name, "_", module_use, "_STATS.txt"), sep="\t")
				
				dimens <- dim(genes_to_use)[1]*0.65
				dimens <- dim(genes_to_use)[1]*0.6
				if(dimens < 10) {
					dimens <- 10
				}
			
				## Lets plot all the plots on one page!!!
				pdf(paste0(outputdir,"/", list_name, "_", module_use, "_LME_ALL.pdf"), width=dimens, height=dimens)
				do.call("grid.arrange", c(lme_list, ncol=ceiling(sqrt(dim(genes_to_use)[1]))))
				dev.off()	
				
				pdf(paste0(outputdir,"/", list_name, "_", module_use, "_RMCORR_ALL.pdf"), width=(dimens+6), height=(dimens+6))
				do.call("grid.arrange", c(rm_corr_list, ncol=ceiling(sqrt(dim(genes_to_use)[1]))))
				dev.off()
				
			}

}
