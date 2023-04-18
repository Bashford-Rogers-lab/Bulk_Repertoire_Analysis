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
		
	### We need to exclude the two bad samples 
	### These samples had ID mix ups and need to be removed!!!
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
	
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	## Include for just samples that are in the dataset (after filtering etc)
	metadata <- metadata[metadata$SampleID_alternative %in% rownames(eigenvectors),]
	
	## We just want day 1 
	eigenvectors_day1 <- eigenvectors[eigenvectors$DAY=="Day1",]
	exclude_cols <- c("dICU_discharge", "dhosp_discharge", "dDeath", "Death.unrelated_cause_description", "dICU_admission", "dBirth", "Diagnosis", "Mortality.Classification", "Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index", "Death.Days_from_ICU_admission", "SRS", "SRSq") 
	comobs <- colnames(metadata)[colnames(metadata) %like% "Comorbidities"]
	mort <- colnames(metadata)[colnames(metadata)%like% "Mortality"]
	
	### we dont want to include these cols twice
	longitudinal_cols <- read.delim(paste0("/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/longitudinal_clinical_variables.txt"), header=FALSE)
	longitudinal_cols <- longitudinal_cols[,1]
	exclude_cols <- unique(c(exclude_cols, comobs, mort,longitudinal_cols))
	
	## On a previous run through I extracted all discrete columns and saved them as a txt file 
	#write.table(discrete_cols2, paste0("/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/discrtete_nonlongitunidal_clinical_variables.txt"), row.names=FALSE, col.names=FALSE)
	## Check format of discrete cols
	## Make No 0 and 1 Yes like the other discrete variables
	metadata$Cardiogenic_Pulmonary_Oedema[metadata$Cardiogenic_Pulmonary_Oedema=="No"] <- 0 
	metadata$Cardiogenic_Pulmonary_Oedema[metadata$Cardiogenic_Pulmonary_Oedema=="Yes"] <- 1
	metadata$Therapy.Mechantical_Ventilation_CPAP[metadata$Therapy.Mechantical_Ventilation_CPAP=="No"] <- 0 
	metadata$Therapy.Mechantical_Ventilation_CPAP[metadata$Therapy.Mechantical_Ventilation_CPAP=="Yes"] <- 1
	metadata$Renal_Failure[metadata$Renal_Failure=="No"] <- 0 
	metadata$Renal_Failure[metadata$Renal_Failure=="Yes"] <- 1
	
	### Lets do this!!!!!!!!
	cols_uses <- c()
	p_values <- c()
	behaviours_all <- c()
#
	for(i in 1:(length(colnames(eigenvectors))-3)){
		module <- colnames(eigenvectors_day1)[i]
		module <- eigenvectors_day1[,c(module, "sample")]
		
		### Just for sepsis data 
		for(x in 6:length(colnames(metadata))){
			discrete <- length(unique(metadata[,x][!is.na(metadata[,x])]))
			col_use <- colnames(metadata)[x] 
			if((discrete<4&discrete>=2 | col_use %like% "Recruitment"| col_use=="SRS" )& !col_use %in% c(exclude_cols) & col_use!="Hypertension"  & !col_use %like% "SOFA" & !col_use %like% "Infection.Number_of_ICU_Acquired_Infections" ){
				## This is for discrete levels 
				#print(x)
				data_subset <- metadata[,c(4,x, 44, 13, 61, 70)]
				s <- table(data_subset[,2])
				data_subset$Sex <- factor(data_subset$Sex, levels=c("Female", "Male"))
				
				if(all(s>1)|col_use=="SRS"){ 
					data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
					data_all[,3] <- as.factor(data_all[,3])
					### Should only be health not included...
					#module$sample[!module$sample %in%  data_subset$SampleID]
					model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2], "+ Age + Sex + Comorbidities.Charlson_Index"))
					xmdl = lm(model_formula, data_all)
					## Print if there are only two classes we take direct 
					if(length(unique(data_all[,3][!is.na(data_all[,3])]))==2){
						print("Only 2 factor levels present compare 0 vs 1")
						p <- summary(xmdl)$coefficients[,4]
						p <- p[2]
						r <- NA
						typev <- "Discrete"
						row_data <- c(colnames(module)[1], colnames(data_subset)[2], p, r,typev)
						p_values <- rbind(row_data, p_values)
						behaviour <- summary(xmdl)$coefficients[,1]
						behaviour <- behaviour[[2]]
						if(behaviour>0){
							direction = "increase"
						} else if (behaviour<0){
							direction = "decrease"
						} else {
							direction = "nochange"
						}
						behaviour_row <- c(colnames(module)[1], col_use, behaviour, direction)
						behaviours_all <- rbind(behaviour_row, behaviours_all)
					} else {
						print("WARNING More than  2 factor levels present compare using ANOVA ")
						## We really want the p value from anova here (as we are studying a categorical variable!) 
						xmdl.av <- aov(xmdl)
						p <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
						r <- NA
						typev <- "Discrete"
						row_data <- c(colnames(module)[1], colnames(data_subset)[2], p, r,typev)
						p_values <- rbind(row_data, p_values)
						#### Beacause these are all done via anova we cant just show direction of effect 
						behaviour <- summary(xmdl)$coefficients[,1]
						behaviour <- behaviour[[2]]
						if(behaviour>0){
							direction = "anova"
						} else if (behaviour<0){
							direction = "anova"
						} else {
							direction = "anova"
						}
						behaviour_row <- c(colnames(module)[1], col_use, behaviour, direction)
						behaviours_all <- rbind(behaviour_row, behaviours_all)
					}
				}
			} else if ((discrete>=4 & !col_use %in% c(exclude_cols) | col_use %like% "Infection.Number_of_ICU_Acquired_Infections")){
				## This is for continuous levels 
				#print(x)
				data_subset <- metadata[,c(4,x, 44, 13, 61, 70)]
				data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
				data_all[,3] <- as.numeric(data_all[,3])
				data_all$Sex <- factor(data_all$Sex, levels=c("Female", "Male"))
				
				model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2], "+Sex+Age+Comorbidities.Charlson_Index"))
				xmdl = lm(model_formula, data_all)
				p <- summary(xmdl)$coefficients[,4]
				p <- p[2]
				r <- NA
				typev <- "Continuous"
				row_data <- c(colnames(module)[1], colnames(data_subset)[2], p, r,typev)
				p_values <- rbind(row_data, p_values)
				### Lets get behaviour
				behaviour <- summary(xmdl)$coefficients[,1]
				behaviour <- behaviour[[2]]
				if(behaviour>0){
					direction = "increase"
				} else if (behaviour<0){
					direction = "decrease"
				} else {
					direction = "nochange"
				}
				behaviour_row <- c(colnames(module)[1], col_use, behaviour, direction)
				behaviours_all <- rbind(behaviour_row, behaviours_all)
				
			} else {
					next
			}
		}
		}
		p_values<- data.frame(p_values)
		colnames(p_values) <- c("Module", "clinical_variable", "P", "Adj_R2", "VariableType")
		
		## What is the direction of the effects 
		spare <- behaviours_all 
		behaviours_all <- data.frame(behaviours_all)
		colnames(behaviours_all) <- c("Module", "clinical_variable", "Effect", "Direction")
		p_values <- merge(p_values, behaviours_all, by=c("Module", "clinical_variable"))

		### Prepare for plotting
		p_values$clinical_variable <- gsub("\\.", "\\:\n", p_values$clinical_variable)
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
		p_values$psig <- "ns"
		p_values$psig[p_values$P < 0.05] <- "sig"
		p_values$psig <- factor(p_values$psig, levels=c("sig", "ns"))

		p_values$Direction[p_values$psig=="ns"] <- NA
		p_values$Direction <- factor(p_values$Direction, levels=c("increase", "decrease", "anova", NA))
		
		p_values$BH <- p.adjust(as.numeric(p_values$P),method="BH")
		p_values$BHsig <- "ns"
		p_values$BHsig[p_values$BH < 0.05] <- "sig"
		
		
		################
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
		
		## Lets plot this 
		if(!is.na(modules_keep)){
			pdf(paste0(outputdir,"/ClinicalData_Day1_LR_INTERESTINGMODULES.pdf"), width=8, height=8)
		}else {
			pdf(paste0(outputdir,"/ClinicalData_Day1_LR_AllModules.pdf"), width=8, height=8)
		}
		
		
		plot(ggplot(p_values, aes(x=Module, y=clinical_variable)) +  geom_tile(aes(alpha=psig, fill = as.numeric(P)),colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Analysis:  Clinical Variables"))+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj \nP < 0.05", values = c("red", '#00000000'), drop = FALSE) + labs(fill="P") +xlab("Module") +ylab(" Clinical Variable")+ geom_point(aes(x=Module, y=clinical_variable, shape = Direction), fill="white", col="black", size=2.5)+ scale_shape_manual(values=c(24,25, 4), labels = c('increase','decrease', 'anova'), na.translate=FALSE, drop=FALSE)+labs(shape="Direction\nof Effect", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE))
		plot(ggplot(p_values, aes(x=Module, y=clinical_variable)) +  geom_tile(aes(alpha=psig, fill = as.numeric(P)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Analysis:  Clinical Variables")) + labs(fill="P") +xlab("Module") +ylab(" Clinical Variable")+ geom_point(aes(x=Module, y=clinical_variable, shape = Direction), fill="white", col="black", size=2.5)+ scale_shape_manual(values=c(24,25, 4), labels = c('increase','decrease', 'anova'), na.translate=FALSE, drop=FALSE)+labs(shape="Direction\nof Effect", alpha="P <0.05")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj \nP < 0.05", values = c("red", '#00000000'), drop = FALSE)+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE))
		dev.off()
		
		##Lets save test values: 
		write.table(p_values, paste0(outputdir, "/Summary/Clinical_variables_lm_model.txt"), sep="\t")
}
##################################################################################################################
