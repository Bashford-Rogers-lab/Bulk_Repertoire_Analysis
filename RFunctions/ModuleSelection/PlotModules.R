## Function to test imputation methods on repertoire data 
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## Jan 2022

library(VIM)
library(Hmisc)
library(missForest) 
library(data.table)
library(mice)
library(missCompare) 
library(xcms)
library(Amelia)
library(ggpubr)
library(missMDA)
library(reshape2)
library(ggplot2)
library(corrplot)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(foreach)
library(doParallel) 
library(ggforce)
library(plot3D)
library(umap)
library(psych)
library(lavaan)
library(cowplot)

plot_eigenvectors <- function(eigenvectors, meta_data1, productivity, outputdir, type_use, subsampled_deptha){
	
	mat_eigenvectors <- data.frame(eigenvectors)
	iso_type <- productivity
	plot_outputdir <- outputdir
	
	# Adding Day column 
	mat_eigenvectors$sample <-  rownames(mat_eigenvectors)
	mat_eigenvectors$DAY <- NA
	mat_eigenvectors$DAY[grep("_1", mat_eigenvectors$sample)] <- "Day1"
	mat_eigenvectors$DAY[grep("_3", mat_eigenvectors$sample)] <- "Day3"
	mat_eigenvectors$DAY[grep("_5", mat_eigenvectors$sample)] <- "Day5"
	mat_eigenvectors$DAY[grep("T1", mat_eigenvectors$sample)] <- "Day1"
	mat_eigenvectors$DAY[grep("T3", mat_eigenvectors$sample)] <- "Day3"
	mat_eigenvectors$DAY[grep("T5", mat_eigenvectors$sample)] <- "Day5"
	mat_eigenvectors$DAY[!mat_eigenvectors$DAY=="Day1" & !mat_eigenvectors$DAY=="Day3" & !mat_eigenvectors$DAY=="Day5"] <- "CONTROL"
	mat_eigenvectors$DISEASE <- "SEPSIS"
	mat_eigenvectors$DISEASE[mat_eigenvectors$sample %like% "HV"] <- "HEALTH"
	mat_eigenvectors$DISEASE[mat_eigenvectors$sample %like% "JR1795_1003"] <- "TECHNICAL"
	mat_eigenvectors$DAY[mat_eigenvectors$sample %like% "JR1795_1003_"] <- "TECHNICAL"
	
	technicals <- mat_eigenvectors[mat_eigenvectors$DISEASE=="TECHNICAL", ]
	mat_eigenvectors <- mat_eigenvectors[!mat_eigenvectors$DISEASE=="TECHNICAL", ]
	
	## Health vs disease onset 
	mat_eigenvectors_hd <- mat_eigenvectors[mat_eigenvectors$DAY %like% "1",]
	mat_eigenvectors_hd$State <- NA
	mat_eigenvectors_hd$State[mat_eigenvectors_hd$DAY %like% "Day1" & mat_eigenvectors_hd$DISEASE=="SEPSIS"] <- "Sepsis Day 1"
	mat_eigenvectors_hd$State[mat_eigenvectors_hd$DAY %like% "Day1" & mat_eigenvectors_hd$DISEASE=="HEALTH"] <- "Health Day 1"
	
	## Plotting eigenvectors by day 
	print("Plotting Eigenvectors by Day") 
	pdf(paste0(outputdir, "Plots/TRY_EIGEN_scaleTRUE_byDAY_", type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=7.5, width=10)
	for(i in 1:(length(colnames(mat_eigenvectors))-3)){
		coluse <- colnames(mat_eigenvectors)[i]
		f <- as.formula(paste0(coluse, "~DAY"))
		a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors)))
		my_comparisons <- list()
		for(z in 1:length(a[,1])){
			my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
		}
		p <- ggplot(mat_eigenvectors, aes_string(x="DAY", y=sym(coluse))) +geom_boxplot(aes(fill=DISEASE)) + geom_point(aes( group=DISEASE), position=position_dodge(width=0.75),)+theme_bw(base_size = 6)+facet_wrap(~DISEASE, scales="free_x") + stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE)
		pt <- ggplot(technicals, aes_string(x="DAY", y=sym(coluse))) +geom_boxplot(fill="yellow") + geom_point(aes( group=DISEASE), position=position_dodge(width=0.75),)+theme_bw(base_size = 6)+facet_wrap(~DISEASE, scales="free_x")
		
		coluse <- colnames(mat_eigenvectors)[i]
		f <- as.formula(paste0(coluse, "~DISEASE"))
		a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors)))
		my_comparisons <- list()
		for(z in 1:length(a[,1])){
			my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
		}
		q <- ggplot(mat_eigenvectors, aes_string(x="DISEASE", y=sym(coluse))) +geom_boxplot(aes(fill=DISEASE)) + geom_point(aes( group=DAY), position=position_dodge(width=0.75),)+theme_bw(base_size = 6)+facet_wrap(~DAY, scales="free_x") + stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE) 
		qt <- ggplot(technicals, aes_string(x="DISEASE", y=sym(coluse))) +geom_boxplot(fill="yellow") + geom_point(aes( group=DAY), position=position_dodge(width=0.75),)+theme_bw(base_size = 6)+facet_wrap(~DAY, scales="free_x") 

		plot(plot_grid(p, pt, q, qt, ncol=2, rel_widths=c(3,1)))
	}
	dev.off()
	
	print("Plotting Module Correllogram") 
	########### plot correlation matrix of eigenvectors
	mat_eigenvectors_numeric <- mat_eigenvectors[,c(1:(length(colnames(mat_eigenvectors))-3))]
	mat_eigenvectors_numeric <- data.frame(mat_eigenvectors_numeric)
	cortest = corr.test(mat_eigenvectors_numeric, adjust="holm")
	pval = cortest$p
	rval = cortest$r
	pdf(paste0(outputdir, "Plots/CorrModules",  type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=10, width=10)
	par(mfrow= c(1,1), mar = c(5,5,5,5))
	corrplot(rval, type="upper", p.mat=pval, insig="label_sig", tl.pos="td", sig.level=0.05, title ="Correlation of eigenvectors", method = "ellipse", diag = F, order = "hclust",  tl.cex = 0.7, mar=c(0,0,2,0))
	dev.off()
	
	### Save eigenvectors
	write.table(mat_eigenvectors, paste0(outputdir, "Eigenvectors_No_Technical_", type_use, "_", iso_type, ".txt"), sep='\t')
	print("DONE")
	
	if(file.exists(meta_data1) & (outputdir %like% "SEPSIS")){
		###############################################################################################################
		## Remove th
		## read in meta data 
		meta_data <- read.delim(meta_data1, sep='\t', header=TRUE)
		rownames(meta_data) <- meta_data$SampleID
		meta_data$SRS_New <- as.factor(meta_data$SRS_New)
		mat_eigenvectors_new <- merge(mat_eigenvectors, meta_data, by = 0, all.x=TRUE)
		rownames(mat_eigenvectors_new) <- mat_eigenvectors_new$Row.names
		mat_eigenvectors_new$Row.names <- NULL
		mat_eigenvectors_new$X7DAY_MORTALITY[mat_eigenvectors_new$DISEASE=="TECHNICAL"] <- "T.CONTROL"
		mat_eigenvectors_new <- mat_eigenvectors_new[mat_eigenvectors_new$DISEASE!="TECHNICAL",]
		mat_eigenvectors_new$SRS <- as.character(mat_eigenvectors_new$SRS)
		mat_eigenvectors_new$SRS[mat_eigenvectors_new$DISEASE=="HEALTH" ] <- "HEALTH"	
		mat_eigenvectors_new$outcome <- NA 
		mat_eigenvectors_new$outcome[mat_eigenvectors_new$DISEASE=="HEALTH"] <- "HEALTHY CONTROL"
		mat_eigenvectors_new$outcome[mat_eigenvectors_new$Days_death_from_ICU=="Alive"] <- "ALIVE"
		mat_eigenvectors_new$outcome[mat_eigenvectors_new$Days_death_from_ICU!="Alive"] <- "DEAD"
		mat_eigenvectors_new$X7DAY_MORTALITY[mat_eigenvectors_new$DISEASE=="HEALTH"] <- "HEALTH"
		mat_eigenvectors_new$X7DAY_MORTALITY[mat_eigenvectors_new$X7DAY_MORTALITY==0] <- "7 Day Mortality"
		mat_eigenvectors_new$X7DAY_MORTALITY[mat_eigenvectors_new$X7DAY_MORTALITY==1] <- "Post 7 Day Mortality"
		mat_eigenvectors_new$X7DAY_MORTALITY[mat_eigenvectors_new$X7DAY_MORTALITY==2] <- "Alive"
		
		mat_eigenvectors_new$X7DAY_MORTALITY <- factor(mat_eigenvectors_new$X7DAY_MORTALITY, levels = c("7 Day Mortality", "Post 7 Day Mortality", "Alive", "HEALTH"))
			
		pdf(paste0(outputdir, "Plots/TRY_EIGEN_scaleTRUE_byVariables_", type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=10, width=10)
		for(i in 1:(length(colnames(mat_eigenvectors))-3)){
			coluse <- colnames(mat_eigenvectors_new)[i]
			##############################
			f <- as.formula(paste0(coluse, "~DAY"))
			a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors_new)))
			my_comparisons <- list()
			for(z in 1:length(a[,1])){
				my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
			}
			p <- ggplot(mat_eigenvectors_new, aes_string(x="DAY", y=sym(coluse))) +geom_boxplot(alpha=0.5, aes(fill=DISEASE)) +theme_bw(base_size = 6)+ facet_wrap(~DISEASE, scales="free_x") + stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE) +guides(fill=FALSE)
			##############################
			f <- as.formula(paste0(coluse, "~X7DAY_MORTALITY"))
			a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors_new)))
			my_comparisons <- list()
			for(z in 1:length(a[,1])){
				my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
			}
			q <- ggplot(mat_eigenvectors_new, aes_string(x="X7DAY_MORTALITY", y=sym(coluse))) +geom_boxplot(alpha=0.5, aes(fill=X7DAY_MORTALITY)) +theme_bw(base_size = 6)+ facet_wrap(~DAY, scales="free_x") + stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+guides(fill=FALSE)
			##############################
			f <- as.formula(paste0(coluse, "~outcome"))
			a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors_new)))
			my_comparisons <- list()
			for(z in 1:length(a[,1])){
				my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
			}
			r <- ggplot(mat_eigenvectors_new, aes_string(x="outcome", y=sym(coluse))) +geom_boxplot(alpha=0.5, aes(fill=outcome)) +theme_bw(base_size = 6)+ facet_wrap(~DAY, scales="free_x") + stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+guides(fill=FALSE)
			##############################
			f <- as.formula(paste0(coluse, "~SRS"))
			a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors_new)))
			my_comparisons <- list()
			for(z in 1:length(a[,1])){
				my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
			}
			s <- ggplot(mat_eigenvectors_new, aes_string(x="SRS", y=sym(coluse))) +geom_boxplot(alpha=0.5, aes(fill=SRS)) +theme_bw(base_size = 6)+ facet_wrap(~DAY, scales="free_x") + stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+guides(fill=FALSE)
			##############################
			f <- as.formula(paste0(coluse, "~Sex"))
			a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors_new)))
			my_comparisons <- list()
			for(z in 1:length(a[,1])){
				my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
			}
			t <- ggplot(mat_eigenvectors_new, aes_string(x="Sex", y=sym(coluse))) +geom_boxplot(alpha=0.5, aes(fill=Sex)) +theme_bw(base_size = 6)+ facet_wrap(~DAY, scales="free_x") + stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+guides(fill=FALSE)
			##############################de
			plot(plot_grid(p, q, r, s, ncol=2))	
		}
	dev.off()
	print("Plotting Done")
}
}
