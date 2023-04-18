## Functions to correlate eigenvectors to LONGITUDINAL clinical variables!!
## Uses linear mixed effect modeling
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## September 2022
## None have more than 2 factor levels so used LME for all of them and just took estimated effect for 0vs1 
library(broom)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)

	
correlate_eigenvectors_continuous <- function(eigenvectors, metadata, outputdir, type_receptor){
	
	### considering sample as well!!
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
		
	### We need to exclude the two bad samples 
	### These samples had ID mix ups and need to be removed!!!
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	
	#all_cols <- colnames(metadata)
	#keep_cols <- all_cols[all_cols %like% "Temperature" | all_cols %like% "WCC" | all_cols %like% "Heart_Rate"| all_cols %like% "Mean_Arteriol_Pressure"| all_cols %like% "Bilirubin.Highest" | all_cols %like% "Systolic_Blood_Pressure"| all_cols %like% "PC02"| all_cols %like% "MPa02"| all_cols %like% "Respiratory_Rate"| all_cols %like% "Platlets"| all_cols %like% "Renal_Failure" | all_cols %like% "Renal_Failure"| all_cols %like% "Creatine"| all_cols %like% "Urea"| all_cols %like% "Hypertension"| all_cols %like% "Bicarbonate"| all_cols %like% "Total.SOFA"| all_cols %like% "GCS."| all_cols %like% "SOFA" | all_cols %like% "Leukocytes"] 
	#keep_cols <- unique(keep_cols)
	#write.table(keep_cols, paste0("/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/longitudinal_clinical_variables.txt"), row.names=FALSE, col.names=FALSE)
	
	## We Want just those variables that vary over time 
	keep_cols <- read.delim("/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/longitudinal_clinical_variables.txt", header=FALSE)
	keep_cols <- keep_cols[,1]
	
	## Lets keep only those columns where we have data for all timepoints!!
	meta_continuous <- metadata[, c("SampleID_alternative", keep_cols, c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index"))]
	meta_continuous$Sex <- factor(meta_continuous$Sex, levels=c("Female", "Male"))

	## Lets add the barcode and the day 
	meta_continuous$Barcode <- NA
	meta_continuous$DAY <- NA
	meta_continuous$Renal_Failure[meta_continuous$Renal_Failure=="No"] <- 0 
	meta_continuous$Renal_Failure[meta_continuous$Renal_Failure=="Yes"] <- 1
	
	for(x in 1:length(meta_continuous$SampleID_alternative)){
			meta_continuous$Barcode[x] <- str_split_fixed (meta_continuous$SampleID_alternative[x], "_", 2)[,1]
			meta_continuous$DAY[x] <- str_split_fixed (meta_continuous$SampleID_alternative[x], "_", 2)[,2]
			if(meta_continuous$SampleID_alternative[x] %like% "HV"){
			 meta_continuous$Barcode[x] <- paste0(str_split_fixed(meta_continuous$SampleID_alternative[x], "_", 3)[,1], "_", str_split_fixed(meta_continuous$SampleID_alternative[x], "_", 3)[,2])
			 meta_continuous$DAY[x] <- paste0(str_split_fixed(meta_continuous$SampleID_alternative[x], "_", 3)[,3],)
			}
		}

	p_values_all <- c()
	behaviour_all <- c()
	for(i in 1:(length(colnames(eigenvectors))-3)){
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY")]
		module$DAY <- gsub("Day", "", module$DAY)
		### Just for sepsis data 
		for(x in 1:length(keep_cols)){
			col_use <- keep_cols[x] 
			discrete <- length(unique(meta_continuous[,col_use][!is.na(meta_continuous[,col_use])]))
				## is it continuous or categorical data 
				## Only have renal failure and that has just two levels!
				if((discrete<4&discrete>=2 | col_use %like% "Recruitment") & col_use!="Hypertension" & col_use!="Infection.Number_of_ICU_Acquired_Infections" & !col_use %like% "SOFA"){
					## This is for discrete levels 		
					print(x)
					print(discrete)
					print(col_use)
					data_subset <- meta_continuous[,c("SampleID_alternative",col_use, "Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")]
					data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
					data_all$Barcode <- str_split_fixed(data_all$sample, "_", 2)[,1]
					data_all[,3] <- as.numeric(data_all[,3])
					data_all[,4] <- as.factor(data_all[,4])
					
					## Run LME model 
					model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2], "+ DAY + Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))
					xmdl =  lmerTest::lmer(model_formula, data_all, REML=F)
					
					### Get P values
					p_vals <- summary(xmdl)$coefficients[,5]
					p_vals <- data.frame(t(p_vals[2:3]))
					no_levels <- length(p_vals)
					colnames(p_vals)[2] <- paste0(colnames(data_subset)[2], ".DAY")
					p_vals$Module <- colnames(module)[1]
					p_vals$typev <- "Continuous"
					
					## Make long format
					cola <- colnames(p_vals)[1]
					colb <- colnames(p_vals)[no_levels]
					p_vals <- gather(p_vals, "Effect", "p_value", all_of(cola):all_of(colb), factor_key=TRUE)
					p_vals$Effect <- gsub("Yes", "", p_vals$Effect)
					## Bind all together
					p_values_all <- rbind(p_vals, p_values_all)
					
					## Lets also get direction of effect 
					behaviour <- summary(xmdl)$coefficients[,1]
					behaviour <- data.frame(behaviour[2:3])
					behaviour$Module <- colnames(module)[1]
					behaviour$typev <- "Continuous"
					behaviour$direction <- "increase" 
					behaviour$direction[behaviour[,1] <0] <- "decrease"
					behaviour$direction[behaviour[,1]==0] <- "nochange"
					rownames(behaviour) <- c(colnames(data_subset)[2],  paste0(colnames(data_subset)[2], ".DAY"))
					behaviour$Effect <- rownames(behaviour)
					colnames(behaviour)[1] <- "EstimatedEffect"
					behaviour_all <- rbind(behaviour, behaviour_all)
				} else {
					## This is for continuous levels 
					print(x)
					data_subset <- meta_continuous[,c("SampleID_alternative",col_use, "Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index", "Barcode")]
					data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
					data_all$Barcode <- str_split_fixed(data_all$sample, "_", 2)[,1]
					data_all[,3] <- as.numeric(data_all[,3])
					data_all[,4] <- as.numeric(data_all[,4])
					
					## Run LME model 
					model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2], "+ DAY + Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))
					xmdl =  lmerTest::lmer(model_formula, data_all, REML=F)
					
					### Get P values
					p_vals <- summary(xmdl)$coefficients[,5]
					p_vals <- data.frame(t(p_vals[2:3]))
					no_levels <- length(p_vals)
					colnames(p_vals)[2] <- paste0(colnames(data_subset)[2], ".DAY")
					p_vals$Module <- colnames(module)[1]
					p_vals$typev <- "Continuous"
					
					## Make long format
					cola <- colnames(p_vals)[1]
					colb <- colnames(p_vals)[no_levels]
					p_vals <- gather(p_vals, "Effect", "p_value", all_of(cola):all_of(colb), factor_key=TRUE)
				
					## Bind all together
					p_values_all <- rbind(p_vals, p_values_all)
					
					## Lets also get direction of effect 
					behaviour <- summary(xmdl)$coefficients[,1]
					behaviour <- data.frame(behaviour[2:3])
					behaviour$Module <- colnames(module)[1]
					behaviour$typev <- "Continuous"
					behaviour$direction <- "increase" 
					behaviour$direction[behaviour[,1] <0] <- "decrease"
					behaviour$direction[behaviour[,1]==0] <- "nochange"
					rownames(behaviour) <- c(colnames(data_subset)[2],  paste0(colnames(data_subset)[2], ".DAY"))
					behaviour$Effect <- rownames(behaviour)
					colnames(behaviour)[1] <- "EstimatedEffect"
					behaviour_all <- rbind(behaviour, behaviour_all)
					
					} 
				}
			}
		#################################################
		### First we want to merge them together 
		p_values_all <- merge(p_values_all, behaviour_all, by=c("Effect", "Module", "typev"))
		
		### Get ready to plot 
		p_values_all$Effect <- gsub("\\.", "\\:\n", p_values_all$Effect)
		p_values_all$Effect <- gsub("_", " ", p_values_all$Effect)
		p_values_all$Effect <- as.factor(p_values_all$Effect)
		p_values_all$Effect <- as.factor(p_values_all$Effect)
		p_values_all$Module <- factor(p_values_all$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25","Module_26","Module_27", "Module_28","Module_29","Module_30", "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40"))
		p_values_all$p_value <- as.numeric(p_values_all$p_value)
		
		## Colour vector for ordering 
		p_values_all <- with(p_values_all,  p_values_all[order(Effect) , ])
		orders <- distinct(p_values_all[, c("Effect", "typev")])
		a <- ifelse(orders$typev == "Continuous", "red", "black")
		
		## with day 
		p_values_all$psig <- "ns"
		p_values_all$psig[p_values_all$p_value < 0.05] <- "sig"
		p_values_all$psig <- factor(p_values_all$psig, levels=c("sig", "ns"))
		
		
		p_values_all$direction[p_values_all$psig=="ns"] <- NA
		p_values_all$direction <- factor(p_values_all$direction, levels=c("increase", "decrease", NA))
		p_values_all$BH <- p.adjust(as.numeric(p_values_all$p_value),method="BH")
		p_values_all$BHsig <- "ns"
		p_values_all$BHsig[p_values_all$BH < 0.05] <- "sig"
		##------------------------------------------------------	
		## without day 
		p_values_all_plot <- p_values_all[!p_values_all$Effect %like% "DAY",]

		### Now lets actually plot!!!!!
		pdf(paste0(outputdir,"/ClinicalData_Continuouswithday.pdf"), width=10, height=20)
		plot(ggplot(p_values_all, aes(x=Module, y=Effect)) +  geom_tile(aes(alpha=psig, fill = as.numeric(p_value)),colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Analysis: Longitudinal Clinical Variables"))+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj \nP < 0.05", values = c("red", '#00000000'), drop = FALSE) + labs(fill="P") +xlab("Module") +ylab("Longitudinal Clinical Variable")+ geom_point(aes(x=Module, y=Effect, shape = direction), fill="white", col="black", size=2.5)+ scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE)+labs(shape="Direction\nof Effect", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE))
		plot(ggplot(p_values_all, aes(x=Module, y=Effect)) +  geom_tile(aes(alpha=psig, fill = as.numeric(p_value)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Analysis: Longitudinal Clinical Variables")) + labs(fill="P") +xlab("Module") +ylab("Longitudinal Clinical Variable")+ geom_point(aes(x=Module, y=Effect, shape = direction), fill="white", col="black", size=2.5)+ scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE)+labs(shape="Direction\nof Effect", alpha="P <0.05")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj \nP < 0.05", values = c("red", '#00000000'), drop = FALSE)+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE))
		dev.off()
		
		pdf(paste0(outputdir,"/ClinicalData_Continuous_no_day.pdf"), width=8, height=8)
		plot(ggplot(p_values_all_plot, aes(x=Module, y=Effect)) +  geom_tile(aes(alpha=psig, fill = as.numeric(p_value)),colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Analysis: Longitudinal Clinical Variables"))+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj \nP < 0.05", values = c("red", '#00000000'), drop = FALSE) + labs(fill="P") +xlab("Module") +ylab("Longitudinal Clinical Variable")+ geom_point(aes(x=Module, y=Effect, shape = direction), fill="white", col="black", size=2.5)+ scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE)+labs(shape="Direction\nof Effect", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE))
		plot(ggplot(p_values_all_plot, aes(x=Module, y=Effect)) +  geom_tile(aes(alpha=psig, fill = as.numeric(p_value)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Analysis: Longitudinal Clinical Variables")) + labs(fill="P") +xlab("Module") +ylab("Longitudinal Clinical Variable")+ geom_point(aes(x=Module, y=Effect, shape = direction), fill="white", col="black", size=2.5)+ scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE)+labs(shape="Direction\nof Effect", alpha="P <0.05")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj \nP < 0.05", values = c("red", '#00000000'), drop = FALSE)+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE))
		dev.off()

		##Lets save test values: 
		write.table(p_values_all, paste0(outputdir, "/Summary/Clinical_variables_lm_model_continous.txt"), sep="\t")
	
}
