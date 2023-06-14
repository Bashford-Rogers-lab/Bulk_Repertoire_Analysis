library(corrplot)
library(psych)
library(rmcorr)
library(matrixStats)
library(RColorBrewer)
plot_correlation_across <- function(eigenvector_list, outputdir, type_receptor){
	
	## Function to merge files 
	my_merge <- function(df1, df2){                                # Create own merging function
		name <- basename(df1)
		name <- gsub("Eigenvectors_No_Technical_", "", name)
		name <- gsub("_PRODUCTIVE.txt", "", name)
		df1 <- read.delim(df1)
		removes <- colnames(df1)[!colnames(df1) %like% "sample" & !colnames(df1) %like% "DAY" & !colnames(df1) %like% "DISEASE"]
		df1 <- df1[, c(removes)]
		colnames(df1) <- paste0(name, "_", colnames(df1))
		#########
		name <- basename(df2)
		name <- gsub("Eigenvectors_No_Technical_", "", name)
		name <- gsub("_PRODUCTIVE.txt", "", name)
		df2 <- read.delim(df2)
		removes <- colnames(df2)[!colnames(df2) %like% "sample" & !colnames(df2) %like% "DAY" & !colnames(df2) %like% "DISEASE"]
		df2 <- df2[, c(removes)]
		colnames(df2) <- paste0(name, "_", colnames(df2))
		full <- merge(df1, df2, by = 0)
	}
	### Lets merge eigenvectors by rowname 
	full <- Reduce(my_merge, eigenvector_list) 
	colnames(full)[1] <- "Sample"
	
	##################################################################################################
	two_types <- str_split_fixed(type_receptor, "_", 2)
	type_1 <- two_types[1]
	type_2 <- two_types[2]
	
	###### Lets do Test!
	mat_eigenvectors <- full
	rownames(mat_eigenvectors) <- mat_eigenvectors$Sample
	mat_eigenvectors$Sample <- NULL
	mat1 <- mat_eigenvectors[,colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% type_1]]
	mat2 <- mat_eigenvectors[,colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% type_2]]
	cortest = corr.test(mat1, mat2, adjust="BH")
	pval = cortest$p
	rval = cortest$r
	## we only want tcr vs bcr 
	pdf(paste0(outputdir, "/AcrossModuleCorrelation_", type_receptor, ".pdf"), height=7, width=7)
	par(mfrow= c(1,1), mar = c(5,5,5,5))
	if(min(pval)<0.05){
		corrplot(rval, p.mat=pval, insig="label_sig", sig.level=0.05, title =paste0("Across Module Correlations BH correction\nAll Days ", type_receptor), method = "circle", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black")
	} else {
		corrplot(rval,  title =paste0("Across Module Correlations BH correction\nAll Days ", type_receptor), method = "circle", order = 'original', tl.cex = 0.7,diag = F, mar=c(0,0,2,0), tl.col = "black")
	}
	dev.off()
	
	#####################################################
	#####################################################################
	## Also want to do by day 
	mat_eigenvectors <- full
	for(x in 1:length(mat_eigenvectors$Sample)){
			mat_eigenvectors$Day[x] <- str_split_fixed(mat_eigenvectors$Sample[x], "_", 2)[,2]
			if(mat_eigenvectors$Sample[x] %like% "HV"){
			 mat_eigenvectors$Day[x] <- str_split_fixed(mat_eigenvectors$Sample[x], "_", 3)[,3]
			}
		}
	mat_eigenvectors_new <- mat_eigenvectors
	
	#######################################################################	
	mat_eigenvectors <- mat_eigenvectors_new[mat_eigenvectors_new$Day=="1",]
	mat_eigenvectors$Day <- NULL
	mat_eigenvectors$Sample <- NULL
	mat1 <- mat_eigenvectors[,colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% type_1]]
	mat2 <- mat_eigenvectors[,colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% type_2]]
	cortest = corr.test(mat1, mat2, adjust="BH")
	pval = cortest$p
	rval = cortest$r
	############ get rs
	rval1 <- data.frame(rval)
	rval1$Chain <- rownames(rval1)
	rval1long <- rval1 %>% gather(chain2, value, -c(Chain))
	colnames(rval1long)[3] <- "Day1"
	###
	colnames(pval) <- colnames(rval)
	rownames(pval) <- rownames(rval)
	#..........................................
	# SET UP PDF 
	pdf(paste0(outputdir, "/AcrossModuleCorrelation_SplitDay_", type_receptor, ".pdf"), height=10, width=15)
	layout(matrix(c(1,2,3,1,2,3,4,5,6), nrow = 3, ncol = 3, byrow = TRUE))
	#par( mfrow= c(2,3) )
	
	if(min(pval)<0.05){
		corrplot(rval, p.mat=pval, insig="label_sig", sig.level=0.05, title =paste0("Receptor Module Correlations (BH Adj p)\nDay1 ", type_receptor), method = "circle", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black")
	} else {
		corrplot(rval,  title =paste0("Receptor Module Correlations (BH Adj p)\nDay1 ", type_receptor), method = "circle", order = 'original', tl.cex = 0.7,diag = F, mar=c(0,0,2,0), tl.col = "black")
	}
	
	################################################
	mat_eigenvectors <- mat_eigenvectors_new[mat_eigenvectors_new$Day=="3",]
	mat_eigenvectors$Day <- NULL
	mat_eigenvectors$Sample <- NULL
	mat1 <- mat_eigenvectors[,colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% type_1]]
	mat2 <- mat_eigenvectors[,colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% type_2]]
	cortest = corr.test(mat1, mat2, adjust="BH")
	pval = cortest$p
	rval = cortest$r
	########
	rval2 <- data.frame(rval)
	rval2$Chain <- rownames(rval2)
	rval2long <- rval2 %>% gather(chain2, value, -c(Chain))
	colnames(rval2long)[3] <- "Day3"
	#########
	if(min(pval)<0.05){
		corrplot(rval, p.mat=pval, insig="label_sig", sig.level=0.05, title =paste0("Receptor Module Correlations (BH Adj p)\nDay3 ", type_receptor), method = "circle", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black")
	} else {
		corrplot(rval,  title =paste0("Receptor Module Correlations (BH Adj p)\nDay3 ", type_receptor), method = "circle", order = 'original', tl.cex = 0.7,diag = F, mar=c(0,0,2,0), tl.col = "black")
	}
	
	###############################
	mat_eigenvectors <- mat_eigenvectors_new[mat_eigenvectors_new$Day=="5",]
	mat_eigenvectors$Day <- NULL
	mat_eigenvectors$Sample <- NULL
	mat1 <- mat_eigenvectors[,colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% type_1]]
	mat2 <- mat_eigenvectors[,colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% type_2]]
	cortest = corr.test(mat1, mat2, adjust="BH")
	pval = cortest$p
	rval = cortest$r
	########
	rval3 <- data.frame(rval)
	rval3$Chain <- rownames(rval3)
	rval3long <- rval3 %>% gather(chain2, value, -c(Chain))
	colnames(rval3long)[3] <- "Day5"
	##############
	if(min(pval)<0.05){
		corrplot(rval, p.mat=pval, insig="label_sig", sig.level=0.05, title =paste0("Receptor Module Correlations (BH Adj p)\nDay5 ", type_receptor), method = "circle", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black")
	} else {
		corrplot(rval,  title =paste0("Receptor Module Correlations (BH Adj p)\nDay5 ", type_receptor), method = "circle", order = 'original', tl.cex = 0.7,diag = F, mar=c(0,0,2,0), tl.col = "black")
	}
	
	####### Look at how r values change 
	#############################################
	## We want to see how the correlations change across time!
	newer <- merge(rval1long, rval2long, by=c("Chain", "chain2"))
	newer <- merge(newer, rval3long, by=c("Chain", "chain2"))
	newer$mean <- rowMeans(newer[,3:5])
	newer$sd <- as.numeric(rowSds(as.matrix(newer[,3:5])))
	cols = brewer.pal(9, "Blues")
	pal = colorRampPalette(cols)
	# Use the following line with RColorBrewer
	newer$order = findInterval(newer$sd, sort(newer$sd))
	plot(newer$Day1, newer$Day3, main="Pairwise Correlation between Variables", xlab="Day 1 R", ylab="Day 3 R", pch=19, col=pal(nrow(newer))[newer$order])
	plot(newer$Day1, newer$Day5, main="Pairwise Correlation between Variables", xlab="Day 1 R", ylab="Day 5 R", pch=19, col=pal(nrow(newer))[newer$order])
	plot(newer$Day3, newer$Day5, main="Pairwise Correlation between Variables", xlab="Day 3 R", ylab="Day 3 R", pch=19, col=pal(nrow(newer))[newer$order])
	dev.off()

	#########################################################
	#### Statistically correct way with repeated measures correlation 
	####rmcorr
	mat_eigenvectors <- full

	## Add barcode
	for(x in 1:length(mat_eigenvectors$Sample)){
			mat_eigenvectors$Barcode[x] <- str_split_fixed (mat_eigenvectors$Sample[x], "_", 2)[,1]
			if(mat_eigenvectors$Sample[x] %like% "HV"){
			 mat_eigenvectors$Barcode[x] <- paste0(str_split_fixed(mat_eigenvectors$Sample[x], "_", 3)[,1], "_", str_split_fixed(mat_eigenvectors$sample[x], "_", 3)[,2])
			}
		}
	users <- colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% "Module"]
	rmcorrbigmatt <- rmcorr_mat(Barcode, users, mat_eigenvectors)
	to_plot <- rmcorrbigmatt$matrix 
	to_plot  <- to_plot[, colnames(to_plot)[colnames(to_plot) %like% type_2]]
	to_plot  <- to_plot[rownames(to_plot)[rownames(to_plot) %like% type_1],]
	
	#####################
	## get p valus 
	x <- rmcorrbigmatt$summary
	x <- x[, c("measure1", "measure2", "p.vals")]
    x <- x[x$measure1 %like% type_1 & x$measure2 %like% type_2 | x$measure1 %like% type_2 & x$measure2 %like% type_1,]
	x$p.vals <- as.numeric(x$p.vals)
	## Adjust p values
	x$p_adjust <- p.adjust(x$p.vals, method = "BH")
	
	## unajusted plot
	m <- matrix(NA, ncol = length(colnames(to_plot)), nrow = length(rownames(to_plot)))
	colnames(m) <- colnames(to_plot)
	rownames(m) <- rownames(to_plot)
	## fill in matrix will p values 
	for(i in 1:length(colnames(m))){
		measure2x <- colnames(m)[i]
			for(j in 1:length(rownames(m))){
				measure1x <- rownames(m)[j] 
				if(measure1x==measure2x){
					m[j,i] <- 1
				} else {
					pfil <- x$p.vals[x$measure1 == measure1x & x$measure2==measure2x]
					if(length(pfil)==0){
						pfil <- x$p.vals[x$measure1 == measure2x & x$measure2==measure1x]
					}
					m[j,i] <- pfil
				}
			}
		}
	
	## adjusted p values 
	madj <- matrix(NA, ncol = length(colnames(to_plot)), nrow = length(rownames(to_plot)))
	colnames(madj) <-  colnames(to_plot)
	rownames(madj) <- rownames(to_plot)
	## fill in matrix will p values 
	for(i in 1:length(colnames(madj))){
		measure2x <- colnames(madj)[i]
			for(j in 1:length(rownames(madj))){
				measure1x <- rownames(madj)[j] 
				if(measure1x==measure2x){
					madj[j,i] <- 1
				} else {
					pfil <- x$p_adjust[x$measure1 == measure1x & x$measure2==measure2x]
					if(length(pfil)==0){
						pfil <- x$p_adjust[x$measure1 == measure2x & x$measure2==measure1x]
					}
					madj[j,i] <- pfil
				}
			}
		}	
	pdf(paste0(outputdir, "/AcrossModuleCorrelation_RMCORR_", type_receptor, ".pdf"), height=7, width=7)
	corrplot(to_plot, p.mat=m, insig="label_sig",  sig.level=0.05, title =paste0("Receptor Module Correlations \nRMCORR: ", type_receptor), method = "circle",  order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black", pch.cex = 2)
	if(min(madj)<0.05){
		corrplot(to_plot, p.mat=madj, insig="label_sig",  sig.level=0.05, title =paste0("Receptor Module Correlations \nRMCORR (BH Adj p): ", type_receptor), method = "circle",  order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black", pch.cex = 2)
	}
	dev.off()
}
