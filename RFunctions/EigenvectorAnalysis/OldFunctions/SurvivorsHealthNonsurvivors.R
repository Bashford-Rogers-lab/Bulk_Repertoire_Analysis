## Functions to correlathealth vs disease and timepoint!
## Comparing Sepsis Day 1 to Health Day 1 function
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## September 2022
library(broom)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)


correlate_eigenvectors_svsn <- function(eigenvectors, metadata, outputdir, type_receptor){
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	p_values_hd <- c()
	for(i in 1:(length(colnames(eigenvectors))-3)){
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		
		## Model with Interaction
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY", "*","DISEASE"))
		xmdl = lm(model_formula, module)
		p_vals <- summary(xmdl)$coefficients[,4]
		xmdl_1 <- data.frame(glance(xmdl))
		p_overall <- xmdl_1$p.value
		row_data <- c(colnames(module)[1], p_vals, p_overall)
		row_data <- data.frame(t(row_data))
		row_data[,2:6] <-as.numeric(row_data[,2:6])
		colnames(row_data) <- c("Module", "IM.Model_p_value", paste0("IM.", names(p_vals))) 
		
		## Model without interaction#
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY", "+","DISEASE"))
		xmdl = lm(model_formula, module)
		p_vals <- summary(xmdl)$coefficients[,4]
		xmdl_1 <- data.frame(glance(xmdl))
		p_overall <- xmdl_1$p.value
		row_data2 <- c(colnames(module)[1], p_vals, p_overall)
		row_data2 <- data.frame(t(row_data2))
		row_data2[,2:5] <-as.numeric(row_data2[,2:5])
		colnames(row_data2) <- c("Module", "Model_p_value", names(p_vals))
		final_row <- cbind(row_data2, row_data[,2:6])
		p_values_hd <- rbind(final_row, p_values_hd)
	}
	 
	colnames(p_values_hd) <- c("Module", "MODEL", "INTERCEPT", "DAY", "DISEASE",  "INT.MODEL", "INT.INTERCEPT", "INT.DAY", "INT.DISEASE","INT.INTERACTION") 
	p_hd <- gather(p_values_hd, "Effect", "p_value", MODEL:INT.INTERACTION, factor_key=TRUE)
	p_hd$Module <- factor(p_hd$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25"))
	orders <- unique(p_hd[, c("Effect")])
	a <- ifelse(orders %like% "INT\\.", "red", "black")
	pdf(paste0(outputdir, "/BCR_ClinicalData_HvsD_Time.pdf"), width=5, height=4)
	plot(ggplot(p_hd, aes(x=Module, y=Effect, alpha=ifelse(p_value <= 0.05, 1, 0.9))) +  geom_tile(aes(fill = as.numeric(p_value)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=a))+ggtitle(paste0(type_receptor, " Interaction vs Additive Model, p<0.05")) + labs(fill="P Value") + guides(alpha="none")+xlab("Module") +ylab("Clinical Variables"))
	dev.off()
	
	write.table(p_hd, paste0(outputdir, "/Summary/Clinical_variables_hvsd_interaction.txt"), sep="\t")
	### Because of the interaction we need to assess at levels of Time / Disease state rather than across
	
	p_hd_extended <- c()
	for(i in 1:(length(colnames(eigenvectors))-3)){
		####################################################
		## Lets do Health first
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		## Lets subset for health only 
		module <- module[module$DISEASE=="HEALTH",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY"))
		xmdl = lm(model_formula, module)
		p_vals <- summary(xmdl)$coefficients[,4]
		## get p val just for day (p overall will be the same!)
		p_vals <- p_vals[2]
		
		####################################################
		## Lets do Sepsis 
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		## Lets subset for health only 
		module <- module[module$DISEASE=="SEPSIS",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY"))
		xmdl = lm(model_formula, module)
		p_vals_s <- summary(xmdl)$coefficients[,4]
		## get p val just for day (p overall will be the same!)
		p_vals_s <- p_vals_s[2]
		
		###########################################
		## Day 1 HvsD
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		## Lets subset for health only 
		module <- module[module$DAY=="Day1",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE"))
		xmdl = lm(model_formula, module)
		p_vals_day1 <- summary(xmdl)$coefficients[,4]
		## get p val just for day (p overall will be the same!)
		p_vals_day1 <- p_vals_day1[2]
		
		###########################################
		## Day 3 HvsD
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		## Lets subset for health only 
		module <- module[module$DAY=="Day3",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE"))
		xmdl = lm(model_formula, module)
		p_vals_day3 <- summary(xmdl)$coefficients[,4]
		## get p val just for day (p overall will be the same!)
		p_vals_day3 <- p_vals_day3[2]
		
		###########################################
		## Day 3 HvsD
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		## Lets subset for health only 
		module <- module[module$DAY=="Day5",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE"))
		xmdl = lm(model_formula, module)
		p_vals_day5 <- summary(xmdl)$coefficients[,4]
		## get p val just for day (p overall will be the same!)
		p_vals_day5 <- p_vals_day5[2]
		
		###################Lets put it all together
		row_use <- c(colnames(eigenvectors)[i], p_vals, p_vals_s, p_vals_day1, p_vals_day3, p_vals_day5)
		row_use[2:5] <- as.numeric(row_use[2:5])
		p_hd_extended <- rbind(row_use, p_hd_extended)
	}
	p_hd_extended <- data.frame(p_hd_extended)
	colnames(p_hd_extended) <- c("Module", "Health.Timepoint", "Sepsis.Timepoint", "HvsS.Day1", "HvsS.Day3", "HvsS.Day5")
	p_hd_extended <- gather(p_hd_extended, "Effect", "p_value", Health.Timepoint:HvsS.Day5, factor_key=TRUE)
	p_hd_extended$p_value <- as.numeric(p_hd_extended$p_value)
	p_hd_extended$Module <- factor(p_hd_extended$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25"))
	
	p_hd_extended$Difference <- NA
	p_hd_extended$Normalised <- NA
	p_hd_extended_2 <- c()
	## Now we want to say wether its higher in sepsis or health 
	for(i in 1:length(p_hd_extended[,1])){
		row_used <- p_hd_extended[i,]
		if(row_used$Effect %like% "HvsS"){
			module_use <- as.character(row_used[,1])
			data_use <- eigenvectors[,c(module_use, "sample", "DAY", "DISEASE")]
			if(row_used$Effect %like% "Day1"){
				DAY <- "Day1"
			}
			if(row_used$Effect %like% "Day3"){
				DAY <- "Day3"
			}
			if(row_used$Effect %like% "Day5"){
				DAY <- "Day5"
			}
			data_use <- data_use[data_use$DAY ==DAY,]
			data_use$DAY <- gsub("Day", "", data_use$DAY)
			data_use$DAY <- as.numeric(data_use$DAY)
			my_sym <- sym(module_use) 
			means <- data_use %>% group_by(DISEASE) %>% summarize(Mean = mean(!!my_sym, na.rm=TRUE))
			means <- data.frame(means)
			diff_means <- means$Mean[means$DISEASE=="SEPSIS"] -  means$Mean[means$DISEASE=="HEALTH"]
			diff_means_normalised <- diff_means/means$Mean[means$DISEASE=="HEALTH"]
			row_used$Difference <- diff_means
			row_used$Normalised <- diff_means_normalised 
		}
		
		if(row_used$Effect %like% "Timepoint"){
			print(i)
			module_use <- as.character(row_used[,1])
			data_use <- eigenvectors[,c(module_use, "sample", "DAY", "DISEASE")]
			data_use$DAY <- gsub("Day", "", data_use$DAY)
			data_use$DAY <- as.numeric(data_use$DAY)
			my_sym <- sym(module_use) 
			means <- data_use %>% group_by(DISEASE, DAY) %>% summarize(Mean = mean(!!my_sym, na.rm=TRUE))
			means <- data.frame(means)
			if(row_used$Effect %like% "Health"){
				diff_means <- means$Mean[means$DISEASE=="HEALTH" & means$DAY==5] -  means$Mean[means$DISEASE=="HEALTH" & means$DAY==1]
				diff_means_normalised <- diff_means/means$Mean[means$DISEASE=="HEALTH" & means$DAY==5]
			} else {
				diff_means <- means$Mean[means$DISEASE=="SEPSIS" & means$DAY==5] -  means$Mean[means$DISEASE=="SEPSIS" & means$DAY==1]
				diff_means_normalised <- diff_means/means$Mean[means$DISEASE=="SEPSIS" & means$DAY==5]
			}
			row_used$Difference <- diff_means
			row_used$Normalised <- diff_means_normalised 
		}
		p_hd_extended_2 <- rbind(row_used, p_hd_extended_2)		
	}
	
	p_hd_extended_2$Difference <- as.numeric(p_hd_extended_2$Difference)
	p_hd_extended_2$Normalised <- as.numeric(p_hd_extended_2$Normalised)
	p_hd_extended_2$High <- NA
	p_hd_extended_2$High[p_hd_extended_2$Difference> 0 & p_hd_extended_2$Effect %like% "HvsS"] <- "High"
	p_hd_extended_2$High[p_hd_extended_2$Difference< 0& p_hd_extended_2$Effect %like% "HvsS"] <- "Low"
	p_hd_extended_2$High[p_hd_extended_2$Difference> 0 & !p_hd_extended_2$Effect %like% "HvsS"] <- "High"
	p_hd_extended_2$High[p_hd_extended_2$Difference< 0& !p_hd_extended_2$Effect %like% "HvsS"] <- "Low"
	p_hd_extended_2$Type <- NA
	p_hd_extended_2$Type[p_hd_extended_2$Effect %like% "HvsS"] <- "HvsS"
	p_hd_extended_2$Type[!p_hd_extended_2$Effect %like% "HvsS"] <- "Timepoint"
	q <- p_hd_extended_2[p_hd_extended_2$p_value <0.05,]
	
	pdf(paste0(outputdir,"/BCR_ClinicalData_Extended.pdf"), width=7.5, height=5)
	 plot(ggplot( p_hd_extended_2 ) +
	  geom_tile(aes(x=Module, y=Effect, fill = as.numeric(p_value), alpha=ifelse(p_value <= 0.05, 1, 0.9))) + 
	  scale_fill_continuous(type = "viridis")+
	  geom_point(data = q, aes(x=Module, y=Effect, shape = High, colour=Type), fill="White") + 
	  scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease')) +
	  scale_colour_discrete(labels = c("Health to Sepsis", "Timepoint 1 to 5")) +
	  theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	  labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("LR Model"))
	dev.off()
	write.table(p_hd_extended, paste0(outputdir, "/Summary/Clinical_variables_hvsd.txt"), sep="\t")
}