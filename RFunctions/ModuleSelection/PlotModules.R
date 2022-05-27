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
	
	mat_eigenvectors <- eigenvectors
	iso_type <- productivity
	plot_outputdir <- outputdir
	
	# Adding Day column 
	DAY <- rownames(mat_eigenvectors)
	DAY <- str_split_fixed(DAY, "_", 2)
	mat_eigenvectors <- data.frame(mat_eigenvectors)
	mat_eigenvectors$DAY <- DAY[,2]
	mat_eigenvectors$DAY[!mat_eigenvectors$DAY=="1" & !mat_eigenvectors$DAY=="3" & !mat_eigenvectors$DAY=="5"] <- "CONTROL"
	mat_eigenvectors$sample = rownames(mat_eigenvectors)
	mat_eigenvectors$sample <- gsub("_1", "", mat_eigenvectors$sample)
	mat_eigenvectors$sample <- gsub("_3", "", mat_eigenvectors$sample)
	mat_eigenvectors$sample <- gsub("_5", "", mat_eigenvectors$sample)
	mat_eigenvectors$sample[mat_eigenvectors$sample %like% "JR1795003_POSITIVE_"] <- "JR1795003_POSITIVE"
	mat_eigenvectors$sample[mat_eigenvectors$sample %like% "JR_HD_T0_POSITIVE"] <- "JR_HD_T0_POSITIVE"
	
	## Health vs disease onset 
	mat_eigenvectors_hd <- mat_eigenvectors[mat_eigenvectors$DAY %like% "1" | mat_eigenvectors$DAY %like% "CONTROL",]
	mat_eigenvectors_hd$State <- NA
	mat_eigenvectors_hd$State[mat_eigenvectors_hd$DAY %like% "1"] <- "Sepsis Day 1"
	mat_eigenvectors_hd$State[mat_eigenvectors_hd$DAY %like% "CONTROL"] <- "Health"
	
	## Plotting eigenvectors by day 
	print("Plotting Eigenvectors by Day") 
	pdf(paste0(outputdir, "Plots/TRY_EIGEN_scaleTRUE_byDAY_", type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=5, width=5)
	for(i in 1:(length(colnames(mat_eigenvectors))-2)){
		coluse <- colnames(mat_eigenvectors)[i]
		f <- as.formula(paste0(coluse, "~DAY"))
		a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors)))
		my_comparisons <- list()
		for(z in 1:length(a[,1])){
			my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
		}
		p <- ggplot(mat_eigenvectors, aes_string(x="DAY", y=sym(coluse), fill="DAY")) +geom_boxplot() + geom_point()+theme_bw(base_size = 6)+stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE)
		plot(p) #stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test"))
	}
	dev.off()
	
	print("Plotting Eigenvectors by HealthvsDisease") 
	## Plotting eigenvectors by Health Vs Disease
	pdf(paste0(outputdir, "Plots/TRY_EIGEN_scaleTRUE_Health_vs_Disease_", type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=5, width=5)
	for(i in 1:(length(colnames(mat_eigenvectors))-2)){
		coluse <- colnames(mat_eigenvectors_hd)[i]
		f <- as.formula(paste0(coluse, "~State"))
		a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors_hd)))
		my_comparisons <- list()
		for(z in 1:length(a[,1])){
			my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
		}
		p <- ggplot(mat_eigenvectors_hd, aes_string(x="State", y=sym(coluse), fill="State")) +geom_boxplot() + geom_point()+theme_bw(base_size = 6)+stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE)
		plot(p) #stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test"))
	}
	dev.off()
	
	print("Plotting Module Correllogram") 
	########### plot correlation matrix of eigenvectors
	mat_eigenvectors_numeric <- mat_eigenvectors[,c(1:(length(colnames(mat_eigenvectors))-2))]
	mat_eigenvectors_numeric <- data.frame(mat_eigenvectors_numeric)
	cortest = corr.test(mat_eigenvectors_numeric, adjust="holm")
	pval = cortest$p
	rval = cortest$r
	pdf(paste0(outputdir, "Plots/CorrModules",  type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=10, width=10)
	par(mfrow= c(1,1), mar = c(5,5,5,5))
	corrplot(rval, type="upper", p.mat=pval, insig="label_sig", tl.pos="td", sig.level=0.05, title ="Correlation of eigenvectors", method = "ellipse", diag = F, order = "hclust",  tl.cex = 0.7, mar=c(0,0,2,0))
	dev.off()
	
	### Save eigenvectors
	write.table(mat_eigenvectors, paste0(outputdir, "Eigenvectors_", type_use, "_", iso_type, ".txt"), sep='\t')
	print("DONE")
	
	if(file.exists(meta_data1)){
		###############################################################################################################
		## Remove th
		## read in meta data 
		meta_data <- read.delim(meta_data1, sep='\t', header=TRUE)
		rownames(meta_data) <- meta_data$SampleID
		meta_data$SRS_New <- as.factor(meta_data$SRS_New)
		mat_eigenvectors_new <- merge(mat_eigenvectors, meta_data, by = 0, all.x=TRUE)
		rownames(mat_eigenvectors_new) <- mat_eigenvectors_new$Row.names
		mat_eigenvectors_new$Row.names <- NULL
		mat_eigenvectors_new$X7DAY_MORTALITY[is.na(mat_eigenvectors_new$X7DAY_MORTALITY)] <- "Control"
		
		# Pre death 
		mat_eigenvectors_new_predeath <- mat_eigenvectors_new[mat_eigenvectors_new$X7DAY_MORTALITY %like% "0" | mat_eigenvectors_new$X7DAY_MORTALITY %like% "1" | mat_eigenvectors_new$DAY %like% "CONTROL",]
		
		df_new <- data.frame()
		for(i in 1:length(mat_eigenvectors_new_predeath$SampleID)){
			if(mat_eigenvectors_new_predeath$DAY[i] %like% "CONTROL"){
				df_new <- rbind(df_new, mat_eigenvectors_new_predeath[i,])
			} else {
				point <- mat_eigenvectors_new_predeath$DAY[i]
				sample <- mat_eigenvectors_new_predeath$sample[i]
				all_points <- mat_eigenvectors_new_predeath$DAY[mat_eigenvectors_new_predeath$sample ==sample]
				last_point <- max(all_points)
				if(point ==last_point){
					df_new <- rbind(df_new, mat_eigenvectors_new_predeath[i,])
					#print("Last point before death")
				} else {
					#print("Not last timepoint")
				} 
			}
		}
		df_new$Diseasepoint <- "Health"
		df_new$Diseasepoint[!df_new$DAY=="CONTROL"] <- "Sepsis Sample Pre-Death"
		mat_eigenvectors_hd_new <- mat_eigenvectors_hd[mat_eigenvectors_hd$sample %in% df_new$sample,]
		
		pdf(paste0(outputdir, "Plots/TRY_EIGEN_scaleTRUE_Health_vs_PreDeath_", type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=4, width=6)
		for(i in 1:(length(colnames(mat_eigenvectors))-2)){
			coluse <- colnames(df_new)[i]
				f <- as.formula(paste0(coluse, "~Diseasepoint"))
				a <- data.frame(do.call(compare_means, list(f, data=df_new)))
				my_comparisons <- list()
				for(z in 1:length(a[,1])){
					my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
				}
				p <- ggplot(df_new, aes_string(x="Diseasepoint", y=sym(coluse), fill="Diseasepoint")) +geom_boxplot() + geom_point()+theme_bw(base_size = 6)+stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE)+ guides(fill="none")+ggtitle("Non-Survivors")+ylim(-2, 10)
				###
				coluse <- colnames(mat_eigenvectors_hd_new)[i]
				f <- as.formula(paste0(coluse, "~State"))
				a <- data.frame(do.call(compare_means, list(f, data=mat_eigenvectors_hd_new)))
				## occasionally there will be 0 in one group causing issues 
				if(dim(a)[1]>0){
					my_comparisons <- list()
					for(z in 1:length(a[,1])){
						my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
					}
					p1<- ggplot(mat_eigenvectors_hd_new, aes_string(x="State", y=sym(coluse), fill="State")) +geom_boxplot() + geom_point()+theme_bw(base_size = 6)+stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE) + guides(fill="none")+ggtitle("Non-Survivors")+ylim(-2, 10)
					####
					plot(plot_grid(p, p1, ncol=2)) #stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test"))
				} else {
				print(coluse)
				}
		}
		dev.off()
		
		pdf(paste0(outputdir, "Plots/TRY_EIGEN_scaleTRUE_byVariables_", type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=5, width=7)
		for(i in 1:(length(colnames(mat_eigenvectors))-2)){
			coluse <- colnames(mat_eigenvectors_new)[i]
			p <- ggplot(mat_eigenvectors_new, aes_string(x="DAY", y=sym(coluse), fill="X7DAY_MORTALITY")) +geom_boxplot(alpha=0.5) +theme_bw(base_size = 6)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")
			plot(p) 
		}
		for(i in 1:(length(colnames(mat_eigenvectors))-2)){
			coluse <- colnames(mat_eigenvectors_new)[i]
			p <- ggplot(mat_eigenvectors_new, aes_string(x="DAY", y=sym(coluse), fill="SRS_New")) +geom_boxplot(alpha=0.5)  +theme_bw(base_size = 6)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")
			plot(p) 
		}
		for(i in 1:(length(colnames(mat_eigenvectors))-2)){
			coluse <- colnames(mat_eigenvectors_new)[i]
			p <- ggplot(mat_eigenvectors_new, aes_string(x="DAY", y=sym(coluse), fill="Sex")) +geom_boxplot(alpha=0.5)+theme_bw(base_size = 6)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")
			plot(p) 
		}
		for(i in 1:(length(colnames(mat_eigenvectors))-2)){
			coluse <- colnames(mat_eigenvectors_new)[i]
			p <- ggplot(mat_eigenvectors_new, aes_string(x="DAY", y=sym(coluse), fill="Shock_sepsis2")) +geom_boxplot(alpha=0.5)+theme_bw(base_size = 6)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")
			plot(p) 
		}
		dev.off()
	}
	print("Plotting Done")
}
