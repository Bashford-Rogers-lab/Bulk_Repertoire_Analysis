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
library(cowplot)
library(AICcmodavg)

correlate_covariates <- function(eigenvectors, metadata, outputdir,type_receptor, metahealth){
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	
	#####################################
	### These samples had ID mix ups and need to be removed!!!
	print("Removing Samples where RNAseq suggests a sample mixup")
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
	############################################

	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	## Include for just samples that are in the dataset (after filtering etc)
	#metadata <- metadata[metadata$SampleID_alternative %in% rownames(eigenvectors),]
	#covariates <- c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")
	#metadata <- metadata[, c("SampleID_alternative", covariates)]
	#health <- read.delim(metahealth, sep="\t", header=TRUE)
	#health <- health[, c("SampleID", "Barcode", "Age", "Sex")]
	#colnames(health)[3:4] <- paste0(colnames(health)[3:4], "_Health")
	#metadata <- metadata[, c("SampleID_alternative", covariates)]
	
	## Preparing Health Metadata
	health <- read.delim(metahealth, sep="\t", header=TRUE)
	health <- health[, c("SampleID", "Barcode", "Age", "Sex")]
	colnames(health) <- c("SampleID_alternative", "Barcode", "Age", "Sex")
	health$Comorbidities.Charlson_Age <- 0 
	health$Comorbidities.Charlson_Index <- 0 
	## Add barcode to sepsis metadata 
	metadata$Barcode <- NA
		for(x in 1:length(metadata$SampleID_alternative)){
			metadata$Barcode[x] <- str_split_fixed (metadata$SampleID_alternative[x], "_", 2)[,1]
			if(metadata$SampleID_alternative[x] %like% "HV"){
			 metadata$Barcode[x] <- paste0(str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,1], "_", str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,2])
			}
		}
	
	## merge meta 
	meta_all <- plyr::rbind.fill(metadata, health)
	covariates <- c("Age", "Sex",  "Comorbidities.Charlson_Index")
	metadata <- meta_all[, c("SampleID_alternative", covariates)]
	## Set the Effect Direction for Sex to make sure its constant!!!!!
	metadata$Sex <- factor(metadata$Sex, levels=c("Female", "Male"))
	
	########################################################################################
	#################################################
	## Lets do sepsis first 
	## We just want day 1 
	eigenvectors_day1 <- eigenvectors[eigenvectors$DAY=="Day1",]
	
	##  P value for single term is the same as a P value for F statistic so we dont need to also report this!
	exclude_cols <- c()
	cols_uses <- c()
	p_values <- c()
	behaviour <- c()
	
	for(i in 1:(length(colnames(eigenvectors))-3)){
		module <- colnames(eigenvectors_day1)[i]
		module <- eigenvectors_day1[,c(module, "sample")]
		
		### Just for sepsis data 
		for(x in 2:length(colnames(metadata))){
			discrete <- length(unique(metadata[,x][!is.na(metadata[,x])]))
			col_use <- colnames(metadata)[x] 
			if((discrete<4&discrete>=2 | col_use %like% "Recruitment"| col_use=="SRS" )& !col_use %in% c(exclude_cols) & col_use!="Hypertension" & col_use!="Infection.Number_of_ICU_Acquired_Infections" & !col_use %like% "SOFA" & !col_use %like% "Comorbidities.Charlson_Age"){
				## This is for discrete levels 
				#print(x)
				data_subset <- metadata[,c(1,x)]
				s <- table(data_subset[,2])
				if(all(s>1)|col_use=="SRS"){ 
					data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
					#data_all[,3] <- as.factor(data_all[,3])
					### Should only be health not included...
					#module$sample[!module$sample %in%  data_subset$SampleID]
					model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2]))
					xmdl = lm(model_formula, data_all)
					xmdl_1 <- data.frame(glance(xmdl))
					p <- xmdl_1$p.value
					r <- xmdl_1$adj.r.squared
					typev <- "Discrete"
					row_data <- c(colnames(module)[1], paste0(colnames(data_subset)[2]), p, r,typev)
					p_values <- rbind(row_data, p_values)
					
					## whats direcition of effect 
					s <-  xmdl$coefficients[2]
					if(s > 0 ) {
						change <- "increase"
					} else if (s < 0 ){
						change <- "decrease"
					} else if (s==0){ 
						change <- "no.change"
					}
					behaviour_row <- c(colnames(module)[1], paste0(colnames(data_subset)[2]), s, change, "discrete")
					behaviour <- rbind(behaviour, behaviour_row)			
					###
					
				}
			} else if ((discrete>=4& !col_use %in% c(exclude_cols)) | col_use=="Hypertension" | col_use=="Infection.Number_of_ICU_Acquired_Infections" | col_use %like% "SOFA"| col_use %like% "Comorbidities.Charlson_Age"){
				## This is for continuous levels 
				print(x)
				data_subset <- metadata[,c(1,x)]
				data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
				data_all[,3] <- as.numeric(data_all[,3])
				model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2]))
				xmdl = lm(model_formula, data_all)
				xmdl_1 <- data.frame(glance(xmdl))
				p <- xmdl_1$p.value
				r <- xmdl_1$adj.r.squared
				typev <- "Continuous"
				row_data <- c(colnames(module)[1], paste0(colnames(data_subset)[2]), p, r,typev)
				p_values <- rbind(row_data, p_values)
				
				## whats direcition of effect 
				s <-  xmdl$coefficients[2]
				if(s > 0 ) {
					change <- "increase"
				} else if (s < 0 ){
					change <- "decrease"
				} else if (s==0){ 
					change <- "no.change"
				}
				behaviour_row <- c(colnames(module)[1], paste0(colnames(data_subset)[2]), s, change, "continuous")
				behaviour <- rbind(behaviour, behaviour_row)	
					
			} else {
					next
			}
		}
		}
		
		##############
		behaviour <- data.frame(behaviour)
		colnames(behaviour) <- c("Module", "clinical_variable", "Estimate", "Direction", "CovType")
		behaviour$clinical_variable <- gsub("\\.", "\\: ", behaviour$clinical_variable)
		behaviour$clinical_variable <- gsub("_", " ", behaviour$clinical_variable)
		behaviour$clinical_variable <- gsub("Comorbidities: Charlson Index", "COMORBIDITIES:CI", behaviour$clinical_variable)
		behaviour$clinical_variable <- gsub("Sex", "SEX.MALE", behaviour$clinical_variable)
		behaviour$clinical_variable <- gsub("Age", "AGE", behaviour$clinical_variable)

		############################################################
		### On to Plotting etc.
		p_values<- data.frame(p_values)
		colnames(p_values) <- c("Module", "clinical_variable", "P", "Adj_R2", "VariableType")
		p_values$clinical_variable <- gsub("\\.", "\\: ", p_values$clinical_variable)
		p_values$clinical_variable <- gsub("_", " ", p_values$clinical_variable)
		p_values$clinical_variable <- gsub("Comorbidities: Charlson Index", "COMORBIDITIES:CI", p_values$clinical_variable)
		p_values$clinical_variable <- gsub("Sex", "SEX.MALE", p_values$clinical_variable)
		p_values$clinical_variable <- gsub("Age", "AGE", p_values$clinical_variable)
		p_values$clinical_variable <- as.factor(p_values$clinical_variable)
		if(all(p_values$Module %like% "Module")=="TRUE"){
			p_values$Module <- factor(p_values$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25","Module_26","Module_27", "Module_28","Module_29","Module_30", "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40"))
		}
		
		## Very important to make numeric or it freaks out!
		p_values$P <- as.numeric(p_values$P)
		
		## Lets investigate the distribution of P values 
		pdf(paste0(outputdir,"/Covariates_Model_PDistribution.pdf"), width=7, height=5)
		plot(ggplot(p_values) +
		geom_histogram(aes(x = P, fill=clinical_variable), breaks = seq(0, 1, 0.05),
					 color = "black") +
		scale_y_continuous(expand = expansion(c(0, 0.05))) +
		facet_wrap(vars(clinical_variable)) + # separate plots
		theme_bw(base_size = 12)	+ geom_vline(xintercept=0.05, col="red") + xlab("Raw p value") +ylab("Count")+ggtitle(paste0("Distribution of P values for LR Models ", type_receptor)) +labs(fill="Model") +guides(fill="none"))
		dev.off()  
	  
		## Lets correct the p value using BH
		p_values$BH <- p.adjust(as.numeric(p_values$P),method="BH")

	 	p_values$psig <- "ns"
		p_values$psig[p_values$P < 0.05] <- "sig"
		p_values$BHsig <- "ns"
		p_values$BHsig[p_values$BH < 0.05] <- "sig"
		p_values_sub <- p_values[p_values$BHsig =="sig",]
		
		p_values$BHsig <- factor(p_values$BHsig, levels=c("sig", "ns"))
		p_values$psig <- factor(p_values$psig, levels=c("sig", "ns"))
		
		## Lets add in the effect direction so that we can add an arrow on the plot indicating direction
		p_values <- merge(p_values, behaviour, by=c("Module", "clinical_variable"))
		## Lets remove the direction for ns effect
		p_values$Direction[p_values$psig=="ns"] <- NA
		p_values$Direction <- factor(p_values$Direction, levels=c("increase", "decrease", NA))
		
		p_values$Mode <-  as.character(p_values$clinical_variable)
		p_values$Mode[p_values$clinical_variable=="AGE"] <- "Model 1"
		p_values$Mode[p_values$clinical_variable=="COMORBIDITIES:CI"] <- "Model 3"
		p_values$Mode[p_values$clinical_variable=="SEX.MALE"] <- "Model 2"
		p_values$Mode <- factor(p_values$Mode, levels=c("Model 3", "Model 2", "Model 1"))
		#######################################################
		## Do a plot showing BHP vs P
		px <- ggplot(p_values) +  geom_point(aes(x=P, y=BH, alpha=psig, colour = BHsig)) +theme_classic() +
		labs(alpha="P<0.05", colour="BH P<0.05")+xlab("P") +ylab("Benjamin Hochburg Adjusted P") + scale_alpha_manual(values = c(1, 0.1)) + scale_color_manual(values = c('red', 'lightblue'))+geom_vline(xintercept=0.05, col="blue") + geom_hline(yintercept=0.05, col="blue")+ ggtitle("Multiple Testing \n Correction")

		########## Lets plot 
		pdf(paste0(outputdir,"/Covariate_Models.pdf"), width=11, height=5.5)		
		p1 <- ggplot(p_values, aes(x=Module, y=clinical_variable)) +  geom_tile(aes(alpha=psig, fill = as.numeric(P)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Day1~Covariate, p<0.05")) + labs(fill="P", alpha="P<0.05") +xlab("Module") +ylab("Covariate")+
		geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop=FALSE)+ geom_point(aes(x=Module, y=clinical_variable, shape = Direction), fill="white", col="black", size=1.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)+facet_grid(rows=vars(Mode),  scales="free_y", space='free_y') 
		plot(plot_grid(p1, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))	
		
		p2 <- ggplot(p_values, aes(x=Module, y=clinical_variable)) +  geom_tile(aes(fill = as.numeric(P), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Day1~Covariate, p<0.05")) + labs(fill="P",alpha="P<0.05")+xlab("Module") +ylab("Covariate")+
		geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop=FALSE)+ geom_point(aes(x=Module, y=clinical_variable, shape = Direction), fill="white", col="black", size=1.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect")+facet_grid(rows=vars(Mode),  scales="free_y", space='free_y')+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		plot(plot_grid(p2, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))
		
		p3 <- ggplot(p_values, aes(x=Module, y=clinical_variable)) +  geom_tile(aes(alpha=psig, fill = as.numeric(BH)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Day1~Covariate, p<0.05")) + labs(fill="BH Adj P", alpha="P<0.05") +xlab("Module") +ylab("Covariate")+
		geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop=FALSE)+ geom_point(aes(x=Module, y=clinical_variable, shape = Direction), fill="white", col="black", size=1.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect")+facet_grid(rows=vars(Mode),  scales="free_y", space='free_y')+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		plot(plot_grid(p3, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))
		
		p4 <- ggplot(p_values, aes(x=Module, y=clinical_variable)) +  geom_tile(aes(fill = as.numeric(BH), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Module Day1~Covariate, p<0.05")) + labs(fill="BH Adj P", alpha="P<0.05") +xlab("Module") +ylab("Covariate")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop=FALSE)+ geom_point(aes(x=Module, y=clinical_variable, shape = Direction), fill="white", col="black", size=1.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect")+facet_grid(rows=vars(Mode),  scales="free_y", space='free_y')+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		plot(plot_grid(p4, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))
		dev.off()
		
		pdf(paste0(outputdir,"/Covariate_Model_BtoBH.pdf"), width=8, height=5)		
		plot(px)	
		dev.off()
		
		##Lets save test values: 
		write.table(p_values, paste0(outputdir, "/Summary/LM_MODELS_Covariates.txt"), sep="\t")
		return(p_values)
}
##################################################################################################################
