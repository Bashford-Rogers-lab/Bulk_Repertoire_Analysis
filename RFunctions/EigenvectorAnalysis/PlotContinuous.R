
	
correlate_eigenvectors_continuous <- function(eigenvectors, metadata, outputdir, type_receptor){
	
	### considering sample as well!!
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	
	## We Want just those variables that vary over time 
	all_cols <- colnames(metadata)
	keep_cols <- all_cols[all_cols %like% "Temperature" | all_cols %like% "WCC" | all_cols %like% "Heart_Rate"| all_cols %like% "Mean_Arteriol_Pressure"| all_cols %like% "Bilirubin.Highest" | all_cols %like% "Systolic_Blood_Pressure"| all_cols %like% "PC02"| all_cols %like% "MPa02"| all_cols %like% "Respiratory_Rate"| all_cols %like% "Platlets"| all_cols %like% "Renal_Failure" | all_cols %like% "Renal_Failure"| all_cols %like% "Creatine"| all_cols %like% "Urea"| all_cols %like% "Hypertension"| all_cols %like% "Bicarbonate"| all_cols %like% "Total.SOFA"| all_cols %like% "GCS."| all_cols %like% "SOFA"| all_cols %like% "SRS"| all_cols %like% "Leukocytes"] 
	keep_cols <- unique(keep_cols)
	
	## Lets keep only those columns where we have data for all timepoints!!
	meta_continuous <- metadata[, c("SampleID_alternative", keep_cols, c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index"))]

	p_values_all <- c()
	for(i in 1:(length(colnames(eigenvectors))-3)){
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY")]
		module$DAY <- gsub("Day", "", module$DAY)
		### Just for sepsis data 
		for(x in 1:length(keep_cols)){
			col_use <- keep_cols[x] 
			discrete <- length(unique(meta_continuous[,col_use][!is.na(meta_continuous[,col_use])]))
				## is it continuous or categorical data 
				if((discrete<4&discrete>=2 | col_use %like% "Recruitment") & col_use!="Hypertension" & col_use!="Infection.Number_of_ICU_Acquired_Infections" & !col_use %like% "SOFA"){
					## This is for discrete levels 		
					print(x)
					data_subset <- meta_continuous[,c("SampleID_alternative",col_use, "Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")]
					data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
					data_all$Barcode <- str_split_fixed(data_all$sample, "_", 2)[,1]
					data_all[,3] <- as.numeric(data_all[,3])
					data_all[,4] <- as.factor(data_all[,4])
					
					model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2], "+DAY+Age+Sex+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index+Barcode"))
					xmdl = lm(model_formula, data_all)
					## Should be anova as it is a categorical variable
					xmdl.av <- aov(xmdl)
					p <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
					d <- summary(xmdl.av)[[1]][["Pr(>F)"]][2]
					p_vals <- data.frame(t(c(p, d)))
					colnames(p_vals) <- c(col_use, paste0(col_use, "_DAY"))
					p_vals$Module <- colnames(module)[1]
					p_vals$typev <- "Discrete"
					cola <- colnames(p_vals)[1]
					no_levels <- length(p_vals)-2
					colb <- colnames(p_vals)[no_levels]
					p_vals <- gather(p_vals, "Effect", "p_value", all_of(cola):all_of(colb), factor_key=TRUE)				
					## Bind all together
					p_values_all <- rbind(p_vals, p_values_all)
					
				} else {
					## This is for continuous levels 
					print(x)
					data_subset <- meta_continuous[,c("SampleID_alternative",col_use, "Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")]
					data_all <- merge(module, data_subset, by.x="sample", by.y="SampleID_alternative")
					data_all$Barcode <- str_split_fixed(data_all$sample, "_", 2)[,1]
					data_all[,3] <- as.numeric(data_all[,3])
					data_all[,4] <- as.numeric(data_all[,4])
					
					model_formula <- formula(paste0(colnames(module)[1], "~", colnames(data_subset)[2], "+DAY+Age+Sex+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index+Barcode"))
					xmdl = lm(model_formula, data_all)
				
					p_vals <- summary(xmdl)$coefficients[,4]
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
					} 
				}
			}

		### Get ready to plot 
		p_values_all$Effect <- gsub("_", " ", p_values_all$Effect)
		p_values_all$Effect <- as.factor(p_values_all$Effect)
		p_values_all$Module <- factor(p_values_all$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25"))
		p_values_all$p_value <- as.numeric(p_values_all$p_value)
		p_values_all <- with(p_values_all,  p_values_all[order(Effect) , ])
		orders <- distinct(p_values_all[, c("Effect", "typev")])
		a <- ifelse(orders$typev == "Continuous", "red", "black")
		
		## with day 
		p_values_all$correct_bp <- as.numeric(p_values_all$p_value)
		p_values_all$correct_bp[p_values_all$correct_bp >= 0.05] <- NA
		
		## without day 
		p_values_all_plot <- p_values_all[!p_values_all$Effect %like% "DAY",]
		orders1 <- distinct(p_values_all_plot[, c("Effect", "typev")])
		b <- ifelse(orders1$typev == "Continuous", "red", "black")
		p_values_all_plot$BH <- p.adjust(as.numeric(p_values_all_plot$p_value),method="BH")
			
		pdf(paste0(outputdir,"/ClinicalData_Continuouswithday.pdf"), width=10, height=20)
		plot(ggplot(p_values_all, aes()) +  geom_tile(aes(x=Module, y=Effect, alpha=ifelse(p_value < 0.05, 1, 0.9), fill = as.numeric(p_value)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y=element_text(colour=a))+ggtitle(paste0(type_receptor, " Module Analysis: Module~Variable+Day+Individual, p<0.05")) + labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("Coefficient"))
		plot(ggplot(p_values_all) +  geom_tile(aes(x=Module, y=Effect, fill = as.numeric(correct_bp)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis",na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=a))+ggtitle(paste0(type_receptor, " Module Analysis: Module~Variable+Day+Individual, p<0.05")) + labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("Coefficient"))
		dev.off()
		
		pdf(paste0(outputdir,"/ClinicalData_Continuous_no_day.pdf"), width=10, height=10)
		plot(ggplot(p_values_all_plot)+  geom_tile(aes(x=Module, y=Effect, alpha=ifelse(p_value < 0.05, 1, 0.9), fill = as.numeric(p_value)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=b))+ggtitle(paste0(type_receptor, " Module Analysis: Module~Variable+Day, p<0.05")) + labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("Coefficient"))
		plot(ggplot(p_values_all_plot)+  geom_tile(aes(x=Module, y=Effect, fill = as.numeric(correct_bp)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis",na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=b))+ggtitle(paste0(type_receptor, " Module Analysis: Module~Variable+Day, p<0.05")) + labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("Coefficient"))
		plot(ggplot(p_values_all_plot)+  geom_tile(aes(x=Module, y=Effect, alpha=ifelse(BH < 0.05, 1, 0.9), fill = as.numeric(BH)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=b))+ggtitle(paste0(type_receptor, " Module Analysis: Module~Variable+Day, p<0.05")) + labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("Coefficient"))
		plot(ggplot(p_values_all_plot)+  geom_tile(aes(x=Module, y=Effect, fill = as.numeric(BH)), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis",na.value="white",limits=c(0,0.05))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=b))+ggtitle(paste0(type_receptor, " Module Analysis: Module~Variable+Day, p<0.05")) + labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("Coefficient"))
		dev.off()
		
		##Lets save test values: 
		write.table(p_values_all, paste0(outputdir, "/Summary/Clinical_variables_lm_model_continous.txt"), sep="\t")

}
