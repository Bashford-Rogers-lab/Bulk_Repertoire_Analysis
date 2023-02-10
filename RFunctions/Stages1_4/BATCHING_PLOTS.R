## Functions to plot the results of the Jaccard Analysis (batch effect)
## Lauren Overend

library(ggpubr)
library(rstatix) 
library(BBmisc)
library(gtools) 
library(coin)

## Multiplot Function
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)

  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                 ncol = cols, nrow = ceiling(numPlots/cols))
}

if (numPlots == 1) {
print(plots[[1]])

} else {
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

for (i in 1:numPlots) {
  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
 }
}
 }
 
 
create_jaccard_plots <- function(path_to_output, jaccard_matrix, layouts1){
	name <- unlist(str_split(i, paste0(path_to_output, "Summary/")))[2]
	name <- gsub(".txt", "", name) 
	layouts <- read.delim(layouts1, sep='\t')
	JACCARD_MATRIX_2 <- read.delim(jaccard_matrix, sep='\t')
	
	JACCARD_MATRIX_2$Sample1 <- gsub("BCR_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("BCR_", "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub("TCRA_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("TCRA_", "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub("TCRB_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("TCRB_", "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub("TCRG_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("TCRG_", "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub("TCRD_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("TCRD_", "", JACCARD_MATRIX_2$Sample2)
	
	## Merge Jaccard Matrix by Layouts 
	layouts_1 <- layouts
	colnames(layouts_1) <- c("SampleID", "Barcode1",  "Lane_S1", "Plate_S1", "Library1", "Position1", "PCRBarcode1")
	JACCARD_MATRIX_2_cohort1 <- merge(JACCARD_MATRIX_2, layouts_1, by.x="Sample1", by.y="SampleID")

	# IDS for sample S2:
	layouts_2 <- layouts
	colnames(layouts_2) <- c("SampleID", "Barcode2",  "Lane_S2", "Plate_S2", "Library2", "Position2", "PCRBarcode2")
	JACCARD_MATRIX_2_cohort1 <- merge(JACCARD_MATRIX_2_cohort1, layouts_2, by.x="Sample2", by.y="SampleID")
	
	## Columns of Sharing PCR Barcode: 
	JACCARD_MATRIX_2_cohort1$SharedInternalBarcode <- "NA" 
	JACCARD_MATRIX_2_cohort1$SharedInternalBarcode[JACCARD_MATRIX_2_cohort1$Sample1 == JACCARD_MATRIX_2_cohort1$Sample2 ] <- "SharedPCRBarcode.SameSample"
	JACCARD_MATRIX_2_cohort1$SharedInternalBarcode[JACCARD_MATRIX_2_cohort1$PCRBarcode1 == JACCARD_MATRIX_2_cohort1$PCRBarcode2 & JACCARD_MATRIX_2_cohort1$Barcode1 == JACCARD_MATRIX_2_cohort1$Barcode2 & JACCARD_MATRIX_2_cohort1$Sample1 != JACCARD_MATRIX_2_cohort1$Sample2 ] <- "SharedPCRBarcode.SameIndividual"
	JACCARD_MATRIX_2_cohort1$SharedInternalBarcode[JACCARD_MATRIX_2_cohort1$PCRBarcode1 == JACCARD_MATRIX_2_cohort1$PCRBarcode2 & JACCARD_MATRIX_2_cohort1$Barcode1 != JACCARD_MATRIX_2_cohort1$Barcode2] <- "SharedPCRBarcode.DifferentIndividual"
	JACCARD_MATRIX_2_cohort1$SharedInternalBarcode[JACCARD_MATRIX_2_cohort1$PCRBarcode1 != JACCARD_MATRIX_2_cohort1$PCRBarcode2 & JACCARD_MATRIX_2_cohort1$Barcode1 != JACCARD_MATRIX_2_cohort1$Barcode2] <- "DifferentPCRBarcode.DifferentIndividual"
	JACCARD_MATRIX_2_cohort1$SharedInternalBarcode[JACCARD_MATRIX_2_cohort1$PCRBarcode1 != JACCARD_MATRIX_2_cohort1$PCRBarcode2 & JACCARD_MATRIX_2_cohort1$Barcode1 == JACCARD_MATRIX_2_cohort1$Barcode2] <- "DifferentPCRBarcode.SameIndividual"
	
	## Columns of Sharing Sequencing Lane: 
	JACCARD_MATRIX_2_cohort1$SharedLane = "NA" 
	JACCARD_MATRIX_2_cohort1$SharedLane[JACCARD_MATRIX_2_cohort1$Lane_S1 == JACCARD_MATRIX_2_cohort1$Lane_S2] <- "SameSeqLane"
	JACCARD_MATRIX_2_cohort1$SharedLane[JACCARD_MATRIX_2_cohort1$Lane_S1 != JACCARD_MATRIX_2_cohort1$Lane_S2] <- "DifferentSeqLane"
	
	## Make a summary Column incoperating the different categories which will be used for plotting
	JACCARD_MATRIX_2_cohort1$Summary <- paste0(JACCARD_MATRIX_2_cohort1$SharedInternalBarcode, ".", JACCARD_MATRIX_2_cohort1$SharedLane)
	
	###
	JACCARD_MATRIX_2_cohort1$ANYNAS <- as.factor(JACCARD_MATRIX_2_cohort1$ANYNAS)
	
	# Make plots 
	pdf(paste0(path_to_output, "/Plots/BoxPlotComparisons_", name, ".pdf"), width=20, height=20)
	if ("LibCorrected" %in% colnames(JACCARD_MATRIX_2)){
		colNames <- names(JACCARD_MATRIX_2_cohort1)[c(3, 5, 6, 8, 9, 11, 12, 14, 16, 18, 19, 21, 22, 23, 26, 27)]
	} else {
		colNames <- names(JACCARD_MATRIX_2_cohort1)[c(3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26)]
	}
	for(s in colNames){
		print(s)
		e <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="SharedInternalBarcode", y=s, fill="SharedInternalBarcode")) + geom_boxplot(alpha = 0.50) + geom_point(aes(color = SharedInternalBarcode), position = position_jitterdodge())  + theme_classic() + facet_wrap(~SharedLane) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
		p <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(s)) + geom_density(aes(fill=SharedLane), alpha=0.8) + facet_wrap(~SharedInternalBarcode, scales = "free") + theme_classic() 
		qs <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(s)) + geom_density(aes(fill=SharedInternalBarcode), alpha=0.8) + facet_wrap(~SharedLane, scales = "free") + theme_classic() +geom_rug()
		tryCatch({ 
		plot(ggarrange(ggarrange(e, p, ncol = 2), qs, nrow = 2))
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		}
	dev.off()
	
	# Make Batch comparisons
	batch_combinations <- combinations(length(unique(JACCARD_MATRIX_2_cohort1$Plate_S2)), 2, unique(JACCARD_MATRIX_2_cohort1$Plate_S2), repeats.allowed=TRUE)
	pdf(paste0(path_to_output, "/Plots/BatchingComparisons_", name, ".pdf"), width=25, height=25)
	if ("LibCorrected" %in% colnames(JACCARD_MATRIX_2)){
		colNames <- names(JACCARD_MATRIX_2_cohort1)[c(3, 5, 6, 8, 9, 11, 12, 14, 16, 18, 19, 21, 22, 24, 25, 27)]
	} else {
		colNames <- names(JACCARD_MATRIX_2_cohort1)[c(3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26)]
	}
	for(s in colNames){
		print(s)
		plots <- list()
		for(i in 1:dim(batch_combinations)[1]){
			subset_1 <- batch_combinations[i,1]
			xlabel <- paste0("Batch ", subset_1)
			subset_2 <- batch_combinations[i,2]
			ylabel <- paste0("Batch ", subset_2)
			minimums <- JACCARD_MATRIX_2_cohort1[,s][!is.na(JACCARD_MATRIX_2_cohort1[,s])]
			minimums_2 <- min(minimums[minimums!="-Inf"])
			maxima <- JACCARD_MATRIX_2_cohort1[,s][!is.na(JACCARD_MATRIX_2_cohort1[,s])]
			maxima <- max(maxima)
			data_use <- JACCARD_MATRIX_2_cohort1[JACCARD_MATRIX_2_cohort1$Plate_S1==subset_1 & JACCARD_MATRIX_2_cohort1$Plate_S2==subset_2,]
			if ("LibCorrected" %in% colnames(JACCARD_MATRIX_2)){
				e <- ggplot(data_use, aes_string("Sample1", "Sample2"))  + geom_tile(aes_string(fill=s)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( limits = c(minimums_2, maxima), colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = data_use[data_use$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'), drop=FALSE) + xlab(paste0("Sample ID: Batch ", subset_1)) + ylab(paste0("Sample ID: Batch ", subset_2))
			} else {
				e <- ggplot(data_use, aes_string("Sample1", "Sample2", fill=s)) + geom_tile(aes_string(fill=s)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn(limits = c(minimums_2, maxima),  colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = data_use[data_use$ANYNAS=="YES",], fill = "white", color='#00000000')+ xlab(paste0("Sample ID: Batch ", subset_1)) + ylab(paste0("Sample ID: Batch ", subset_2))
			} 
			plots[[i]] <- e
			}
		tryCatch({ 
			multiplot(plotlist = plots, cols = 5)
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	dev.off()
	
	## Plot the statistics 
	calculate_statistics_plot(path_to_output, JACCARD_MATRIX_2_cohort1, name)
}



calculate_statistics_plot <- function(path_to_output, jaccard_matrix, name){
	name <- name
	JACCARD_MATRIX_2_cohort1 <- jaccard_matrix
	
	## Plotting comparisons between groups using GPUBR t test with unequal variance. 
	pdf(paste0(path_to_output, "/Plots/SummaryStatistics_", name, ".pdf"), width=10, height=10)
	
	#---------------------------------------------------------------------
	tryCatch({
	stat.test <- JACCARD_MATRIX_2_cohort1 %>% t_test(SharedSeq.MeanSubsample ~ Summary, var.equal=FALSE) %>%
	  adjust_pvalue(method = "bonferroni") %>%
	  add_significance()
	stat.test 
	stat.test <- stat.test %>% add_xy_position(x="Summary")
	e <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="Summary", y="SharedSeq.MeanSubsample")) + theme_classic() + geom_boxplot(alpha = 0.50, aes(fill=Summary)) + geom_point(aes(color = Summary), position = position_jitterdodge()) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
	x <- e + stat_pvalue_manual(stat.test, hide.ns = TRUE)  
	plot(x)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	#---------------------------------------------------------------------
	tryCatch({
	stat.test <- JACCARD_MATRIX_2_cohort1 %>% t_test(Jaccard.MeanSubsample ~ Summary, var.equal=FALSE) %>%
	  adjust_pvalue(method = "bonferroni") %>%
	  add_significance()
	stat.test 
	stat.test <- stat.test %>% add_xy_position(x="Summary")
	e <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="Summary", y="Jaccard.MeanSubsample")) + theme_classic() + geom_boxplot(alpha = 0.50, aes(fill=Summary)) + geom_point(aes(color = Summary), position = position_jitterdodge()) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
	x <- e + stat_pvalue_manual(stat.test, hide.ns = TRUE)  
	plot(x)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	#---------------------------------------------------------------------
	tryCatch({
	stat.test <- JACCARD_MATRIX_2_cohort1 %>% t_test(Jaccard.Full ~ Summary, var.equal=FALSE) %>%
	  adjust_pvalue(method = "bonferroni") %>%
	  add_significance()
	stat.test 
	stat.test <- stat.test %>% add_xy_position(x="Summary")
	e <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="Summary", y="Jaccard.Full")) + theme_classic() + geom_boxplot(alpha = 0.50, aes(fill=Summary)) + geom_point(aes(color = Summary), position = position_jitterdodge()) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
	x <- e + stat_pvalue_manual(stat.test, hide.ns = TRUE)  
	plot(x)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	#---------------------------------------------------------------------
	#---------------------------------------------------------------------
	tryCatch({
	stat.test <- JACCARD_MATRIX_2_cohort1 %>% t_test(SharedSeq.MeanSubsample.Weighted ~ Summary, var.equal=FALSE) %>%
	  adjust_pvalue(method = "bonferroni") %>%
	  add_significance()
	stat.test 
	stat.test <- stat.test %>% add_xy_position(x="Summary")
	e <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="Summary", y="SharedSeq.MeanSubsample.Weighted")) + theme_classic() + geom_boxplot(alpha = 0.50, aes(fill=Summary)) + geom_point(aes(color = Summary), position = position_jitterdodge()) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
	x <- e + stat_pvalue_manual(stat.test, hide.ns = TRUE)  
	plot(x)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	#---------------------------------------------------------------------
	tryCatch({
	stat.test <- JACCARD_MATRIX_2_cohort1 %>% t_test(Jaccard.MeanSubsample.Weighted ~ Summary, var.equal=FALSE) %>%
	  adjust_pvalue(method = "bonferroni") %>%
	  add_significance()
	stat.test 
	stat.test <- stat.test %>% add_xy_position(x="Summary")
	e <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="Summary", y="Jaccard.MeanSubsample.Weighted")) + theme_classic() + geom_boxplot(alpha = 0.50, aes(fill=Summary)) + geom_point(aes(color = Summary), position = position_jitterdodge()) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
	x <- e + stat_pvalue_manual(stat.test, hide.ns = TRUE)  
	plot(x)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	#---------------------------------------------------------------------
	tryCatch({
	stat.test <- JACCARD_MATRIX_2_cohort1 %>% t_test(Jaccard.Full.Weighted ~ Summary, var.equal=FALSE) %>%
	  adjust_pvalue(method = "bonferroni") %>%
	  add_significance()
	stat.test 
	stat.test <- stat.test %>% add_xy_position(x="Summary")
	e <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="Summary", y="Jaccard.Full.Weighted")) + theme_classic() + geom_boxplot(alpha = 0.50, aes(fill=Summary)) + geom_point(aes(color = Summary), position = position_jitterdodge()) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
	x <- e + stat_pvalue_manual(stat.test, hide.ns = TRUE)  
	plot(x)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	#---------------------------------------------------------------------
	dev.off() 
}










