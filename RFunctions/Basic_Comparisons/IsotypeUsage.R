
library(dplyr)
library(tidyr) 

imp_data_file <- read.delim('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_BCR_PRODUCTIVE.txt', header=TRUE)
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'

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


	pdf(paste0(outputdir, "/SHMRatio.pdf"), width=13, height=6)
	ggplot(shm_usagex, aes(x=factor(DAY), y=mean_mutations, fill=disease)) +geom_boxplot()+facet_grid( cols=vars(Metric))+theme_classic()+ylab("Non-silent:silent mutations ratio")+ggtitle("SHM")
	shm_usagex$DAY <- as.numeric(shm_usagex$DAY)
	ggplot(shm_usagex, aes(x=DAY, y=mean_mutations, colour=disease)) +geom_point()+geom_smooth(method="lm")+facet_grid(cols=vars(Metric) )+theme_classic()+ylab("Non-silent:silent mutations ratio")+ggtitle("SHM")
	dev.off()

}