## Code to plot SHM in a neat format for thesis
## Explore CSR and non-silent mutations for downstream analysis

library(gpubr)
library(dplyr)
library(tidyr) 

#imp_data_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_BCR_PRODUCTIVE.txt'
#outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'

plot_summaries <- function(imp_data_file, outputdir){
	
	imp_data <- read.delim(imp_data_file, header=TRUE)
	
	## First lets look at isotype usage 
	iso_usage_cols <- colnames(imp_data)[colnames(imp_data) %like% "_vdjs_per_isotype"]
	iso_usage <- imp_data[, c(iso_usage_cols)]
	iso_usage$sample <- rownames(iso_usage)
	iso_usagex <- gather(iso_usage, "Metric", "percentage", -sample)
	iso_usagex$Type[iso_usagex$Metric %like% "total"] <- "TOTAL"
	iso_usagex$Type[iso_usagex$Metric %like% "unique"] <- "UNIQUE"
	iso_usagex$Metric <- gsub("BCR_READS_percentage_total_vdjs_per_isotype__", "", iso_usagex$Metric)
	iso_usagex$Metric <- gsub("BCR_READS_percentage_unique_vdjs_per_isotype__", "", iso_usagex$Metric)
	iso_usagex$sample <- gsub("_productive", "", iso_usagex$sample)
	iso_usagex$disease <- "SEPSIS"
	iso_usagex$disease[iso_usagex$sample %like% "HV"] <- "HEALTH"

	## need to take out bad ids
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	iso_usagex <- iso_usagex[!iso_usagex$sample %in% c(bad_ids),]
	iso_usagex <- iso_usagex[!iso_usagex$sample %like% "JR1795_1003_POSITIVE",]

	## Good to continue 
	iso_usagex$DAY <- NA
	iso_usagex$DAY[grep("_1", iso_usagex$sample)] <- "1"
	iso_usagex$DAY[grep("_3", iso_usagex$sample)] <- "3"
	iso_usagex$DAY[grep("_5", iso_usagex$sample)] <- "5"


	pdf(paste0(outputdir, "/ISOUSAGE.pdf"), width=10, height=5)
	ggplot(iso_usagex, aes(x=factor(DAY), y=percentage, fill=disease)) +geom_boxplot()+facet_grid(rows=vars(Type), cols=vars(Metric))+theme_classic()+ylab("% of Repertoire")+ggtitle("Isotype Usage")
	iso_usagex$DAY <- as.numeric(iso_usagex$DAY)
	ggplot(iso_usagex, aes(x=DAY, y=percentage, colour=disease)) +geom_point()+geom_smooth(method="lm")+facet_grid(rows=vars(Type), cols=vars(Metric))+theme_classic()+ylab("% of Repertoire")+ggtitle("Isotype Usage")
	dev.off()
	
	pdf(paste0(outputdir, "/ISOUSAGE2.pdf"), width=10, height=10)
	p2 <- ggline(iso_usagex[iso_usagex$Type=="UNIQUE",], x = "DAY", y = "percentage", add = "mean_se", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percentage", facet.by="Metric", main="Unique Reads")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p3 <- ggline(iso_usagex[iso_usagex$Type=="TOTAL",], x = "DAY", y = "percentage", add = "mean_se", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percentage", facet.by="Metric", main="Total Reads")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	plot_grid(p2, p3, ncol = 2, align="hv", axis="tblr")
	dev.off()


	##------------------------------------
	##------------------------------------
	
	
	## Now SHM
	shm_cols <- colnames(imp_data)[colnames(imp_data) %like% "mean_mutations"]
	
	shm_usage <- imp_data[, c(shm_cols)]
	shm_usage$sample <- rownames(shm_usage)
	shm_usagex <- gather(shm_usage, "Metric", "mean_mutations", -sample)

	shm_usagex$Metric <- gsub("BCR_READS_mean_mutations_per_vdj__", "", shm_usagex$Metric)
	shm_usagex$sample <- gsub("_productive", "", shm_usagex$sample)
	shm_usagex$disease <- "SEPSIS"
	shm_usagex$disease[shm_usagex$sample %like% "HV"] <- "HEALTH"

	## need to take out bad ids
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	shm_usagex <- shm_usagex[!shm_usagex$sample %in% c(bad_ids),]
	shm_usagex <- shm_usagex[!shm_usagex$sample %like% "JR1795_1003_POSITIVE",]

	## Good to continue 
	shm_usagex$DAY <- NA
	shm_usagex$DAY[grep("_1", shm_usagex$sample)] <- "1"
	shm_usagex$DAY[grep("_3", shm_usagex$sample)] <- "3"
	shm_usagex$DAY[grep("_5", shm_usagex$sample)] <- "5"
	
	## Lets remove those with different class switch tyoes 
	shm_usagex <- shm_usagex[!shm_usagex$Metric == "expanded" & !shm_usagex$Metric == "unexpanded" & !shm_usagex$Metric == "Class_switched",]
	shm_usagex1 <- shm_usagex
	shm_usagex1$Type <- "Mean Mutations per VDJ"

	pdf(paste0(outputdir, "/SHM.pdf"), width=13, height=6)
	ggplot(shm_usagex, aes(x=factor(DAY), y=mean_mutations, fill=disease)) +geom_boxplot()+facet_grid( cols=vars(Metric))+theme_classic()+ylab("Mean Mutations per VDJ")+ggtitle("Isotype Usage")
	shm_usagex$DAY <- as.numeric(shm_usagex$DAY)
	ggplot(shm_usagex, aes(x=DAY, y=mean_mutations, colour=disease)) +geom_point()+geom_smooth(method="lm")+facet_grid(cols=vars(Metric) )+theme_classic()+ylab("Mean Mutations per VDJ")+ggtitle("Isotype Usage")
	dev.off()
	
	######################################################
	
	shm_cols <- colnames(imp_data)[colnames(imp_data) %like% "nonsilent_silent"]
	
	shm_usage <- imp_data[, c(shm_cols)]
	shm_usage$sample <- rownames(shm_usage)
	shm_usagex <- gather(shm_usage, "Metric", "mean_mutations", -sample)

	shm_usagex$Metric <- gsub("BCR_READS_mean_nonsilent_silent_mutation_ratio_per_VDJ__", "", shm_usagex$Metric)
	shm_usagex$sample <- gsub("_productive", "", shm_usagex$sample)
	shm_usagex$disease <- "SEPSIS"
	shm_usagex$disease[shm_usagex$sample %like% "HV"] <- "HEALTH"

	## need to take out bad ids
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	shm_usagex <- shm_usagex[!shm_usagex$sample %in% c(bad_ids),]
	shm_usagex <- shm_usagex[!shm_usagex$sample %like% "JR1795_1003_POSITIVE",]

	## Good to continue 
	shm_usagex$DAY <- NA
	shm_usagex$DAY[grep("_1", shm_usagex$sample)] <- "1"
	shm_usagex$DAY[grep("_3", shm_usagex$sample)] <- "3"
	shm_usagex$DAY[grep("_5", shm_usagex$sample)] <- "5"
	shm_usagex <- shm_usagex[!shm_usagex$Metric == "IGHD_IGHM_mutated",]
	shm_usagex$Type <- "Non-silent:silent mutation ratio"
	
	### Specify levels
	shm_usagex$Metric <- factor(shm_usagex$Metric, levels=c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE"))
	shm_usagex1$Metric <- factor(shm_usagex1$Metric, levels=c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE"))
	
	pdf(paste0(outputdir, "/SHMRatio.pdf"), width=13, height=6)
	ggplot(shm_usagex, aes(x=factor(DAY), y=mean_mutations, fill=disease)) +geom_boxplot()+facet_grid( cols=vars(Metric))+theme_classic()+ylab("Non-silent:silent mutations ratio")+ggtitle("SHM")
	shm_usagex$DAY <- as.numeric(shm_usagex$DAY)
	ggplot(shm_usagex, aes(x=DAY, y=mean_mutations, colour=disease)) +geom_point()+geom_smooth(method="lm")+facet_grid(cols=vars(Metric) )+theme_classic()+ylab("Non-silent:silent mutations ratio")+ggtitle("SHM")
	dev.off()
	
	### Lets bind together for plotting 
	full <- rbind(shm_usagex, shm_usagex1)
	
	shm_usagex1$mean_mutations <- as.numeric(shm_usagex1$mean_mutations)
	shm_usagex$mean_mutations <- as.numeric(shm_usagex$mean_mutations)
	colnames(shm_usagex1)[2] <- "ISOTYPE"
	colnames(shm_usagex)[2] <- "ISOTYPE"
				  
	pdf(paste0(outputdir, "/SHM_Final.pdf"), width=10, height=7)
	p1 <- ggline(shm_usagex1, x = "DAY", y = "mean_mutations", add = "mean_se", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Mean Mutations per VDJ", facet.by="ISOTYPE")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p2 <- ggline(shm_usagex, x = "DAY", y = "mean_mutations", add = "mean_se", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Non-silent:silent mutations ratio", facet.by="ISOTYPE")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	plot_grid(p1, p2, ncol = 2, align="hv", axis="tblr")
	
	p3 <- ggline(shm_usagex1, x = "DAY", y = "mean_mutations", add = "mean_se", color = "ISOTYPE", ylab="Mean Mutations per VDJ", facet.by="disease")+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p3.1 <- ggline(shm_usagex1[shm_usagex1$ISOTYPE=="IGHM" |shm_usagex1$ISOTYPE=="IGHD",], x = "DAY", y = "mean_mutations", add = "mean_se", color = "ISOTYPE", ylab="Mean Mutations per VDJ", facet.by="disease")+ stat_compare_means(aes(group = ISOTYPE), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p4 <- ggline(shm_usagex, x = "DAY", y = "mean_mutations", add = "mean_se", color = "ISOTYPE", ylab="Non-silent:silent mutations ratio", facet.by="disease")+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p4.1 <- ggline(shm_usagex[shm_usagex$ISOTYPE=="IGHM" |shm_usagex$ISOTYPE=="IGHD",], x = "DAY", y = "mean_mutations", add = "mean_se", color = "ISOTYPE", ylab="Non-silent:silent mutations ratio", facet.by="disease")+ stat_compare_means(aes(group = ISOTYPE), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	plot_grid(p3, p3.1, p4, p4.1, ncol = 2, align="hv", axis="tblr")
	dev.off()
	
	pdf(paste0(outputdir, "/SHM_Final2.pdf"), width=11, height=14)
	plot_grid(p1, p2, p3, p3.1, p4, p4.1, ncol = 2, align="hv", axis="tblr" ,rel_heights=c(2,1,1), labels="AUTO")
	dev.off()
	
}