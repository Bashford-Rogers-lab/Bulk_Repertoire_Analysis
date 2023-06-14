## Functions to plot the results of the Jaccard Analysis (batch effect)
## Lauren Overend

suppressMessages(library(ggpubr))
suppressMessages(library(rstatix) )
suppressMessages(library(BBmisc))
suppressMessages(library(gtools) )
suppressMessages(library(coin))
suppressMessages(library(stringr))
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

#jaccard_matrix <- "/well/immune-rep/shared/MISEQ/TEST_PIPE/Summary/JACCARDMATRIX_2ormore_AllSAMPLES_LIBCONTAM_CORRECTED_filter_TEST.txt"
#path_to_output <-"/well/immune-rep/shared/MISEQ/TEST_PIPE/"
#layouts1 <- "LEO_SEPSIS_BCR_ALL_layouts.txt"

create_jaccard_plots <- function(path_to_output, jaccard_matrix, layouts1){
	name <- unlist(str_split(jaccard_matrix, paste0(path_to_output, "Summary/")))[2]
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
	
	### Are they the same sample
	JACCARD_MATRIX_2_cohort1$SharedSample <- NA
	JACCARD_MATRIX_2_cohort1$SharedSample[JACCARD_MATRIX_2_cohort1$Sample1 == JACCARD_MATRIX_2_cohort1$Sample2 ] <- "Same RNA Sample"
	JACCARD_MATRIX_2_cohort1$SharedSample[JACCARD_MATRIX_2_cohort1$Sample1 != JACCARD_MATRIX_2_cohort1$Sample2 ] <- "Different RNA Sample"

	### Are they the same Individual 
	JACCARD_MATRIX_2_cohort1$SharedIndividual <- NA
	JACCARD_MATRIX_2_cohort1$SharedIndividual[JACCARD_MATRIX_2_cohort1$Barcode1 == JACCARD_MATRIX_2_cohort1$Barcode2 ] <- "Y"
	JACCARD_MATRIX_2_cohort1$SharedIndividual[JACCARD_MATRIX_2_cohort1$Barcode1 != JACCARD_MATRIX_2_cohort1$Barcode2 ] <- "N"

	### Are they the same PCR Barcode 
	JACCARD_MATRIX_2_cohort1$SharedPCRBarcode <- NA
	JACCARD_MATRIX_2_cohort1$SharedPCRBarcode[JACCARD_MATRIX_2_cohort1$PCRBarcode1 == JACCARD_MATRIX_2_cohort1$PCRBarcode2 ] <- "Y"
	JACCARD_MATRIX_2_cohort1$SharedPCRBarcode[JACCARD_MATRIX_2_cohort1$PCRBarcode1 != JACCARD_MATRIX_2_cohort1$PCRBarcode2 ] <- "N"
	
	### Are they the same Lane  
	JACCARD_MATRIX_2_cohort1$SharedLane <- NA
	JACCARD_MATRIX_2_cohort1$SharedLane[JACCARD_MATRIX_2_cohort1$Lane_S1 == JACCARD_MATRIX_2_cohort1$Lane_S2 ] <- "Same\nSequencing Lane"
	JACCARD_MATRIX_2_cohort1$SharedLane[JACCARD_MATRIX_2_cohort1$Lane_S1 != JACCARD_MATRIX_2_cohort1$Lane_S2 ] <- "Different\nSequencing Lane"
	
		
	## Make a summary Column incoperating the different categories which will be used for plotting
	#JACCARD_MATRIX_2_cohort1$Summary <- paste0(JACCARD_MATRIX_2_cohort1$SharedInternalBarcode, ".", JACCARD_MATRIX_2_cohort1$SharedLane)
	JACCARD_MATRIX_2_cohort1$ANYNAS <- as.factor(JACCARD_MATRIX_2_cohort1$ANYNAS)
	
	give.n <- function(x){
	  return(c(y = median(x)*1.05, label = length(x))) 
	  # experiment with the multiplier to find the perfect position
	}

	## lets FILTER NAS  
	#JACCARD_MATRIX_2_cohort_X <- JACCARD_MATRIX_2_cohort1[!is.na(JACCARD_MATRIX_2_cohort1$Comparison),]

	# Make plots 
	pdf(paste0(path_to_output, "/Plots/JACCARD_MATRIX_BOXPLOTS_BoxPlotComparisons_", name, ".pdf"), width=9, height=5)
	if ("LibCorrected" %in% colnames(JACCARD_MATRIX_2)){
		colNames <- names(JACCARD_MATRIX_2_cohort1)[c(3, 5, 6, 8, 9, 11, 12, 14, 16, 18, 19, 21, 22, 23, 26, 27)]
	} else {
		colNames <- names(JACCARD_MATRIX_2_cohort1)[c(3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26)]
	}
	for(s in colNames){
		print(s)
		labsx <- gsub("\\.", " ", s)
		labsx <- gsub("log", "log2(",labsx)
		labsx <- gsub("SharedSeq", "Sequence Overlap",labsx)
		labsx <- gsub("Jaccard", "Jaccard Index",labsx)
		labsx <- gsub("MeanSubsample", "Subsampled",labsx)
		labsx <- gsub("Full", "Full Data",labsx)
		if(length(grep("log2", labsx, value=TRUE)>1)){
			labsx <- paste0(labsx, " +c)")
		}
		
		#e <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="SharedInternalBarcode", y=s, fill="SharedInternalBarcode", color="SharedInternalBarcode")) + geom_boxplot(alpha = 0.50) + geom_point(aes(color = SharedInternalBarcode), position = position_jitterdodge())  + theme_classic() + facet_wrap(~SharedLane) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
		p1 <- ggplot(JACCARD_MATRIX_2_cohort1, aes_string(x="SharedPCRBarcode", y=s, fill="SharedIndividual")) +geom_boxplot(aes(fill=SharedIndividual, colour=SharedIndividual), alpha=0.5) + geom_point(aes(color = SharedIndividual), position = position_jitterdodge())+facet_grid(~SharedSample+SharedLane)+theme_classic()+xlab("Shared Internal PCR Barcode") +ylab(labsx)+
		scale_color_manual(values=c("#F8766D", "#619CFF"))+scale_fill_manual(values=c("#F8766D", "#619CFF"))  + stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75), vjust = -1, colour="black") +labs(colour="Same\nIndividual", fill="Same\nIndividual") +
		scale_y_continuous(expand=c(0.5, 0.5))
		#plot(e)
		plot(p1)
		}
	dev.off()
	
	# Make Batch comparisons
	batch_combinations <- combinations(length(unique(JACCARD_MATRIX_2_cohort1$Plate_S2)), 2, unique(JACCARD_MATRIX_2_cohort1$Plate_S2), repeats.allowed=TRUE)
	
	widthx <- 7+(length(unique(JACCARD_MATRIX_2$Sample1))*0.05)
	heightx <- 7+(length(unique(JACCARD_MATRIX_2$Sample1))*0.05)
	
	if(widthx<20){
		widthx=20
	}
	
	if(heightx<10){
		heightx=10
	}
	
	pdf(paste0(path_to_output, "/Plots/JACCARD_MATRIX_BATCH_COMPARISONS_", name, ".pdf"), width=widthx, height=heightx)
	if ("LibCorrected" %in% colnames(JACCARD_MATRIX_2)){
		colNames <- names(JACCARD_MATRIX_2_cohort1)[c(3, 5, 6, 8, 9, 11, 12, 14, 16, 18, 19, 21, 22, 24, 25, 27)]
	} else {
		colNames <- names(JACCARD_MATRIX_2_cohort1)[c(3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26)]
	}
	for(s in colNames){
		print(s)
		namex <- gsub("\\.", "\n", s)
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
				e <- ggplot(data_use, aes_string("Sample1", "Sample2"))  + geom_tile(aes_string(fill=s)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( limits = c(minimums_2, maxima), colours=c("navyblue", "darkorange1"), name=namex)  + geom_tile(aes(width = 1, height = 1), data = data_use[data_use$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'), drop=FALSE) + xlab(paste0("Sample ID: Batch ", subset_1)) + ylab(paste0("Sample ID: Batch ", subset_2))
			} else {
				e <- ggplot(data_use, aes_string("Sample1", "Sample2", fill=s)) + geom_tile(aes_string(fill=s)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn(limits = c(minimums_2, maxima),  colours=c("navyblue", "darkorange1"), name=namex)  + geom_tile(aes(width = 1, height = 1), data = data_use[data_use$ANYNAS=="YES",], fill = "white", color='#00000000')+ xlab(paste0("Sample ID: Batch ", subset_1)) + ylab(paste0("Sample ID: Batch ", subset_2))
			} 
			plots[[i]] <- e
			}
		lengthx <- length(plots)
		ncol = 5
		nrow = ceiling(lengthx/ncol)
		tryCatch({ 
			#multiplot(plotlist = plots, cols = 5)
			#grid_arrange_shared_legend(plotlist = plots, ncol = 5, position="right")
			#wrap_plots(plotlist = plots) + plot_layout(guides = "collect")
			eval(parse( text = paste0("grid_arrange_shared_legend(", paste0("plots", "[[", c(1:lengthx), "]]", sep = '', collapse = ','), ",ncol =", ncol, ",nrow =", nrow, ", position = 'right')", sep = '')))
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	dev.off()
	
	## Plot the statistics 
	#calculate_statistics_plot(path_to_output, JACCARD_MATRIX_2_cohort1, name)
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










