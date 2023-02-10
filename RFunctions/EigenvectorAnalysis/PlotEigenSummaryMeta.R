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

#eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
#outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'
#meta_data1 <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
#type_receptor <- "BCR"
#productivity <- "PRODUCTIVE"

plot_eigenvectors <- function(eigenvectors, meta_data1, outputdir ){
	mat_eigenvectors <- data.frame(read.delim(eigenvectors))
	plot_outputdir <- paste0(outputdir, "/MORTALITYPLOTS")
	if (!dir.exists(plot_outputdir)) {dir.create(plot_outputdir)}
	
	meta_data <- read.delim(meta_data1, sep='\t', header=TRUE)

	modules <- colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% "Module"]
	mortality <- colnames(meta_data)[colnames(meta_data) %like% "Mortality"]
	

	for(i in 1:length(modules)){
	  
	  modulesusex <- modules[i]
	  pdf(paste0(plot_outputdir, "/MortalityBoxPlotSummary_", modulesusex, ".pdf"), height=5, width=10)
	  for(c in 1:length(mortality)){
			datax <- mat_eigenvectors[, c(modules[i], "DAY", "DISEASE")]
			datay <- meta_data[, c("SampleID_alternative", mortality[c])]
			rownames(datay) <- datay$SampleID_alternative
			datay$SampleID_alternative <- NULL
			new <- merge(datay, datax, by=0, all.y=TRUE)
			new[,2][new$DISEASE=="HEALTH"] <- "HEALTH"

			moduleuse <- colnames(new)[3]
			variableuse <- colnames(new)[2]
			
			moduleuse <- colnames(new)[3]
			variableuse <- colnames(new)[2]
			##############################
			f <- as.formula(paste0(moduleuse, "~",variableuse ))
			a <- data.frame(do.call(compare_means, list(f, data=new)))
			my_comparisons <- list()
			for(z in 1:length(a[,1])){
				my_comparisons[[z]] <- c(a$group1[z], a$group2[z])
			}
			p <- ggplot(new, aes_string(x=sym(variableuse), y=sym(moduleuse))) +geom_boxplot(alpha=0.5, aes_string(fill=sym(variableuse))) +theme_classic(base_size = 6)+ xlab("Mortality")+ facet_wrap(~DAY, scales="free_x") + stat_compare_means(comparisons =my_comparisons, label = "p.signif", exact=FALSE) +guides(fill=FALSE)
			##############################
			plot(p)	
		}
	dev.off()
	print("Plotting Done")
	}
}
