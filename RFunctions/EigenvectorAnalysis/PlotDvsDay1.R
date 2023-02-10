## Functions to correlate eigenvectors to clinical variables!!
## Comparing Sepsis Day 1 to Health Day 1 function
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## September 2022
library(broom)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)

correlate_eigenvectors_t1 <- function(eigenvectors, metadata, outputdir,type_receptor, modules_keep=NA){
	
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	if(!is.na(modules_keep)){
		eigenvectors <- eigenvectors[, c( modules_keep,  "sample", "DAY", "DISEASE")]
	}	
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	## Include for just samples that are in the dataset (after filtering etc)
	metadata <- metadata[metadata$SampleID_alternative %in% rownames(eigenvectors),]
	
	## We just want day 1 
	eigenvectors_day1 <- eigenvectors[eigenvectors$DAY=="Day1",]
	exclude_cols <- c("dICU_discharge", "dhosp_discharge", "dDeath", "Death.unrelated_cause_description", "dICU_admission", "dBirth", "Diagnosis", "Mortality.Classification", "Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index", "Death.Days_from_ICU_admission") 
	comobs <- colnames(metadata)[colnames(metadata) %like% "Comorbidities"]
	mort <- colnames(metadata)[colnames(metadata)%like% "Mortality"]

	exclude_cols <- unique(c(exclude_cols, comobs, mort))
	cols_uses <- c()
	p_values <- c()
	
	for(i in 1:(length(colnames(eigenvectors))-3)){
		module <- colnames(eigenvectors_day1)[i]
		module <- eigenvectors_day1[,c(module, "sample")]
		
		### Just for sepsis data 
		for(x in 6:length(colnames(metadata))){
			discrete <- length(unique(metadata[,x][!is.na(metadata[,x])]))
			col_use <- colnames(metadata)[x] 
			if((discrete<4&discrete>=2 | col_use %like% "Recruitment"| col_use=="SRS" )& !col_use %in% c(exclude_cols) & col_use!="Hypertension"  & !col_use %like% "SOFA" ){
				## This is for discrete levels 
				#print(x)
				data_subset <- metadata[,c(4,x, 44, 13, 61, 70)]
				s <- table(data_subset[,2])
				if(all(s>1)|col_use=="SRS"){ 
					data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
					data_all[,3] <- as.factor(data_all[,3])
					### Should only be health not included...
					#module$sample[!module$sample %in%  data_subset$SampleID]
					model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2], "+Sex+Age+Comorbidities.Charlson_Index"))
					xmdl = lm(model_formula, data_all)
					## We really want the p value from anova here (as we are studying a categorical variable!) 
					xmdl.av <- aov(xmdl)
					p <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
					r <- NA
					typev <- "Discrete"
					row_data <- c(colnames(module)[1], colnames(data_subset)[2], p, r,typev)
					p_values <- rbind(row_data, p_values)
				}
			} else if ((discrete>=4& !col_use %in% c(exclude_cols)) | col_use=="Hypertension" | col_use %like% "SOFA"){
				## This is for continuous levels 
				print(x)
				data_subset <- metadata[,c(4,x, 44, 13, 61, 70)]
				data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
				data_all[,3] <- as.numeric(data_all[,3])
				model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2], "+Sex+Age+Comorbidities.Charlson_Index"))
				xmdl = lm(model_formula, data_all)
				p <- summary(xmdl)$coefficients[,4]
				p <- p[2]
				r <- NA
				typev <- "Continuous"
				row_data <- c(colnames(module)[1], colnames(data_subset)[2], p, r,typev)
				p_values <- rbind(row_data, p_values)
			} else {
					next
			}
		}
		}
		p_values<- data.frame(p_values)
		colnames(p_values) <- c("Module", "clinical_variable", "P", "Adj_R2", "VariableType")
		p_values$clinical_variable <- gsub("\\.", "\\: ", p_values$clinical_variable)
		p_values$clinical_variable <- gsub("_", " ", p_values$clinical_variable)
		p_values$clinical_variable <- as.factor(p_values$clinical_variable)
		if(all(p_values$Module %like% "Module")=="TRUE"){
			p_values$Module <- factor(p_values$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25","Module_26","Module_27", "Module_28","Module_29","Module_30", "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40"))
		}
		p_values <- with(p_values,  p_values[order(clinical_variable) , ])
		orders <- distinct(p_values[, c("clinical_variable", "VariableType")])
		a <- ifelse(orders$VariableType == "Continuous", "red", "black")
		
		## Very important to make numeric or it freaks out!
		p_values$P <- as.numeric(p_values$P)
		
		## Lets investigate the distribution of P values 
		if(!is.na(modules_keep)){
			pdf(paste0(outputdir,"/ClinicalData_PDistribution_INTERESTINGMODULES.pdf"), width=30, height=30)
		}else {
			pdf(paste0(outputdir,"/ClinicalData_PDistribution_AllModules.pdf"), width=30, height=30)
		}
		plot(ggplot(p_values) +
		geom_histogram(aes(x = P, fill=clinical_variable), breaks = seq(0, 1, 0.05),
					 color = "black") +
		scale_y_continuous(expand = expansion(c(0, 0.05))) +
		facet_wrap(vars(clinical_variable)) + # separate plots
		theme_bw(base_size = 12)	+ geom_vline(xintercept=0.05, col="red") + xlab("Raw p value") +ylab("Count")+ggtitle(paste0("Distribution of P values for LR Models ", type_receptor)) +labs(fill="Model") +guides(fill="none"))
		dev.off()  
	  
		p_values$BH <- p.adjust(as.numeric(p_values$P),method="BH")

		if(!is.na(modules_keep)){
			pdf(paste0(outputdir,"/ClinicalData_Day1_LR_INTERESTINGMODULES.pdf"), width=10, height=11)
		}else {
			pdf(paste0(outputdir,"/ClinicalData_Day1_LR_AllModules.pdf"), width=10, height=11)
		}
		plot(ggplot(p_values) +  geom_tile(aes(x=Module, y=clinical_variable, alpha=ifelse(P < 0.05, 1, 0.9), fill = as.numeric(P)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=a))+ggtitle(paste0(type_receptor, " Module.Day1~Variables, p<0.05")) + labs(fill="Model p Value") + guides(alpha="none")+xlab("Module") +ylab("Clinical Variables"))
		plot(ggplot(p_values) +  geom_tile(aes(x=Module, y=clinical_variable, fill = as.numeric(P)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=a))+ggtitle(paste0(type_receptor, " Module.Day1~Variables, p<0.05")) + labs(fill="Model p Value") + guides(alpha="none")+xlab("Module") +ylab("Clinical Variables"))
		plot(ggplot(p_values) +  geom_tile(aes(x=Module, y=clinical_variable, alpha=ifelse(BH < 0.05, 1, 0.9), fill = as.numeric(BH)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=a))+ggtitle(paste0( type_receptor, " Module.Day1~Variables, p<0.05")) + labs(fill="BH Asj p Value") + guides(alpha="none")+xlab("Module") +ylab("Clinical Variables"))
		plot(ggplot(p_values) +  geom_tile(aes(x=Module, y=clinical_variable, fill = as.numeric(BH)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=a))+ggtitle(paste0(type_receptor, " Module.Day1~Variables, p<0.05")) + labs(fill="BH Adj p Value") + guides(alpha="none")+xlab("Module") +ylab("Clinical Variables"))
		dev.off()
		
		##Lets save test values: 
		write.table(p_values, paste0(outputdir, "/Summary/Clinical_variables_lm_model.txt"), sep="\t")
}
##################################################################################################################
