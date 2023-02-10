library(corrplot)
library(psych)
library(rmcorr)
library(corrr)
library(tidyverse)
library(igraph)
library(ggraph)


plot_correlation_within <- function(eigenvectors, outputdir, type_receptor){
	########### plot correlation matrix of eigenvectors
	mat_eigenvectors <- data.frame(read.delim(eigenvectors))
	mat_eigenvectors_all <- mat_eigenvectors
	## Remove columns that are not eigenvectors
	removes <- colnames(mat_eigenvectors)[!colnames(mat_eigenvectors) %like% "sample" & !colnames(mat_eigenvectors) %like% "DAY" & !colnames(mat_eigenvectors) %like% "DISEASE"]
	mat_eigenvectors <- mat_eigenvectors[, c(removes)]
	cortest = corr.test(mat_eigenvectors, adjust="BH")
	pval = cortest$p
	rval = cortest$r
	pdf(paste0(outputdir, "/WithinModuleCorrelation_AllDays.pdf"), height=6, width=6)
	par(mfrow= c(1,1), mar = c(5,5,5,5))
	corrplot(rval, type="upper", p.mat=pval, insig="label_sig", tl.pos="td", sig.level=0.05, title =paste0("Within Module Correlations BH correction: ", type_receptor), method = "ellipse", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0))
	dev.off()
	
	pdf(paste0(outputdir, "/WithinModuleCorrelation_AllDays_Network.pdf"), height=7, width=10)
	rval2 <- data.frame(rval)
	rval2$measure1 <- rownames(rval2)
	rval2 <- rval2 %>% gather(measure2, r, -c(measure1)) 
	pval2 = data.frame(pval)
	pval2$measure1 <- rownames(pval2)
	pval2 <- pval2 %>% gather(measure2, p, -c(measure1)) 
	pval2$p <- as.numeric(pval2$p)
	all <- merge(rval2, pval2, by=c("measure1", "measure2"))
	all <- all %>%filter(abs(p) <0.05)
	rval2 <- all[, c("measure1", "measure2", "r")]
	rval2 <- rval2[!rval2$measure1==rval2$measure2,]
	rval2 <- rval2 %>%filter(abs(r) > 0.3)
	graph_cors <- rval2 %>%graph_from_data_frame(directed = FALSE)
	plot(ggraph(graph_cors) +geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +guides(edge_alpha = "none", edge_width = "none") +scale_edge_colour_viridis() +geom_node_point(color = "black", size = 3) +geom_node_text(aes(label = name), repel = TRUE) +theme_graph(base_family="sans") +labs(title = paste0("Correlation Network between ", type_receptor, " modules (BH adj p<0.05 & R>0.3)")))
	dev.off() 

	################################################################################
	## Also want to do by day 
	mat_eigenvectors <- mat_eigenvectors_all[mat_eigenvectors_all$DAY=="Day1",]
	removes <- colnames(mat_eigenvectors)[!colnames(mat_eigenvectors) %like% "sample" & !colnames(mat_eigenvectors) %like% "DAY" & !colnames(mat_eigenvectors) %like% "DISEASE"]
	mat_eigenvectors <- mat_eigenvectors[, c(removes)]
	cortest = corr.test(mat_eigenvectors, adjust="BH")
	pval = cortest$p
	rval = cortest$r
	rval1 <- rval
	pval1 <- pval
	pdf(paste0(outputdir, "/WithinModuleCorrelation_SplitDay.pdf"), height=6, width=15)
	par( mfrow= c(1,3) )
	corrplot(rval, type="upper", p.mat=pval, insig="label_sig", tl.pos="td", sig.level=0.05, title =paste0("Within Module Correlations BH correction Day 1: ", type_receptor), method = "ellipse", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0))
	##
	mat_eigenvectors <- mat_eigenvectors_all[mat_eigenvectors_all$DAY=="Day3",]
	removes <- colnames(mat_eigenvectors)[!colnames(mat_eigenvectors) %like% "sample" & !colnames(mat_eigenvectors) %like% "DAY" & !colnames(mat_eigenvectors) %like% "DISEASE"]
	mat_eigenvectors <- mat_eigenvectors[, c(removes)]
	cortest = corr.test(mat_eigenvectors, adjust="BH")
	pval = cortest$p
	rval = cortest$r
	rval3 <- rval
	pval3 <- pval
	corrplot(rval, type="upper", p.mat=pval, insig="label_sig", tl.pos="td", sig.level=0.05, title =paste0("Within Module Correlations BH correction Day 3: ", type_receptor), method = "ellipse", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0))
	##
	mat_eigenvectors <- mat_eigenvectors_all[mat_eigenvectors_all$DAY=="Day5",]
	removes <- colnames(mat_eigenvectors)[!colnames(mat_eigenvectors) %like% "sample" & !colnames(mat_eigenvectors) %like% "DAY" & !colnames(mat_eigenvectors) %like% "DISEASE"]
	mat_eigenvectors <- mat_eigenvectors[, c(removes)]
	cortest = corr.test(mat_eigenvectors, adjust="BH")
	pval = cortest$p
	rval = cortest$r
	rval5 <- rval
	pval5 <- pval
	corrplot(rval, type="upper", p.mat=pval, insig="label_sig", tl.pos="td", sig.level=0.05, title =paste0("Within Module Correlations BH correction Day 5: ", type_receptor), method = "ellipse", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0))
	dev.off()
	
	pdf(paste0(outputdir, "/WithinModuleCorrelation_SplitDay_Network.pdf"), height=15, width=15)
	### Day 1 
	rval2 <- data.frame(rval1)
	rval2$measure1 <- rownames(rval2)
	rval2 <- rval2 %>% gather(measure2, r, -c(measure1)) 
	pval2 = data.frame(pval1)
	pval2$measure1 <- rownames(pval2)
	pval2 <- pval2 %>% gather(measure2, p, -c(measure1)) 
	pval2$p <- as.numeric(pval2$p)
	all <- merge(rval2, pval2, by=c("measure1", "measure2"))
	all <- all %>%filter(abs(p) <0.05)
	rval2 <- all[, c("measure1", "measure2", "r")]
	rval2 <- rval2[!rval2$measure1==rval2$measure2,]
	rval2 <- rval2 %>%filter(abs(r) > 0.3)
	graph_cors <- rval2 %>%graph_from_data_frame(directed = FALSE)
	p1 <- ggraph(graph_cors) + geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +guides(edge_alpha = "none", edge_width = "none") +scale_edge_colour_viridis() + geom_node_point(color = "black", size = 3) +geom_node_text(aes(label = name), repel = TRUE) + theme_graph(base_family="sans") +labs(title = paste0(type_receptor, " modules Day 1 (BH adj p<0.05 & R>0.3)"))
	
	##
	rval2 <- data.frame(rval3)
	rval2$measure1 <- rownames(rval2)
	rval2 <- rval2 %>% gather(measure2, r, -c(measure1)) 
	pval2 = data.frame(pval3)
	pval2$measure1 <- rownames(pval2)
	pval2 <- pval2 %>% gather(measure2, p, -c(measure1)) 
	pval2$p <- as.numeric(pval2$p)
	all <- merge(rval2, pval2, by=c("measure1", "measure2"))
	all <- all %>%filter(abs(p) <0.05)
	rval2 <- all[, c("measure1", "measure2", "r")]
	rval2 <- rval2[!rval2$measure1==rval2$measure2,]
	rval2 <- rval2 %>%filter(abs(r) > 0.3)
	graph_cors <- rval2 %>%graph_from_data_frame(directed = FALSE)
	p2 <- ggraph(graph_cors) + geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +guides(edge_alpha = "none", edge_width = "none") +scale_edge_colour_viridis() + geom_node_point(color = "black", size = 3) +geom_node_text(aes(label = name), repel = TRUE) + theme_graph(base_family="sans") +labs(title = paste0(type_receptor, " modules Day 3 (BH adj p<0.05 & R>0.3)"))
	##
	rval2 <- data.frame(rval5)
	rval2$measure1 <- rownames(rval2)
	rval2 <- rval2 %>% gather(measure2, r, -c(measure1)) 
	pval2 = data.frame(pval5)
	pval2$measure1 <- rownames(pval2)
	pval2 <- pval2 %>% gather(measure2, p, -c(measure1)) 
	pval2$p <- as.numeric(pval2$p)
	all <- merge(rval2, pval2, by=c("measure1", "measure2"))
	all <- all %>%filter(abs(p) <0.05)
	rval2 <- all[, c("measure1", "measure2", "r")]
	rval2 <- rval2[!rval2$measure1==rval2$measure2,]
	rval2 <- rval2 %>%filter(abs(r) > 0.3)
	graph_cors <- rval2 %>%graph_from_data_frame(directed = FALSE)
	p3 <- ggraph(graph_cors) + geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +guides(edge_alpha = "none", edge_width = "none") +scale_edge_colour_viridis() + geom_node_point(color = "black", size = 3) +geom_node_text(aes(label = name), repel = TRUE) + theme_graph(base_family="sans") +labs(title = paste0(type_receptor, " modules Day 5 (BH adj p<0.05 & R>0.3)"))
	plot(plot_grid(p1, p2, p3, ncol=2))
	dev.off() 
	
	
	###################################################
	#######################################
	#### Statistically correct way with repeated measures correlation 
	####rmcorr
	mat_eigenvectors <- mat_eigenvectors_all
	
	## Add barcode
	for(x in 1:length(mat_eigenvectors$sample)){
			mat_eigenvectors$Barcode[x] <- str_split_fixed (mat_eigenvectors$sample[x], "_", 2)[,1]
			if(mat_eigenvectors$sample[x] %like% "HV"){
			 mat_eigenvectors$Barcode[x] <- paste0(str_split_fixed(mat_eigenvectors$sample[x], "_", 3)[,1], "_", str_split_fixed(mat_eigenvectors$sample[x], "_", 3)[,2])
			}
		}
	
	rmcorrbigmatt <- rmcorr_mat(Barcode, removes, mat_eigenvectors)
	to_plot <- rmcorrbigmatt$matrix 
	## get p valus 
	x <- rmcorrbigmatt$summary
	x <- x[, c("measure1", "measure2", "p.vals")]
	## Adjust p values
	x$p_adjust <- p.adjust(x$p.vals, method = "BH")
	
	## unajusted plot
	m <- matrix(NA, ncol = length(removes), nrow = length(removes))
	colnames(m) <- removes
	rownames(m) <- removes
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
	madj <- matrix(NA, ncol = length(removes), nrow = length(removes))
	colnames(madj) <- removes
	rownames(madj) <- removes
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
	pdf(paste0(outputdir, "/WithinModuleCorrelation_RMCORR.pdf"), height=6, width=6)
	corrplot(to_plot, type="upper", p.mat=m, insig="label_sig", tl.pos="td", sig.level=0.05, title =paste0("Within Module Correlations RMCORR: ", type_receptor), method = "ellipse", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0))
	corrplot(to_plot, type="upper", p.mat=madj, insig="label_sig", tl.pos="td", sig.level=0.05, title =paste0("Within Module Correlations RMCORR BH adj: ", type_receptor), method = "ellipse", diag = F, order = 'original', tl.cex = 0.7, mar=c(0,0,2,0))
	dev.off()
	
	## Lets do a network using rmcorr 
	pdf(paste0(outputdir, "/WithinModuleCorrelation_RMCORR_Network.pdf"), height=7, width=10)
	rval2 <- data.frame(to_plot)
	rval2$measure1 <- rownames(rval2)
	rval2 <- rval2 %>% gather(measure2, r, -c(measure1))
	pval2 = data.frame(x)
	pval2$p <- as.numeric(pval2$p.vals)
	pval2$p <- as.numeric(pval2$p_adjust)
	all <- merge(rval2, pval2, by=c("measure1", "measure2"))
	all <- all %>%filter(abs(p_adjust) <0.05)
	rval2 <- all[, c("measure1", "measure2", "r")]
	rval2 <- rval2[!rval2$measure1==rval2$measure2,]
	rval2 <- rval2 %>%filter(abs(r) > 0.3)
	rval2 <- rval2[!rval2$measure1==rval2$measure2,]
	graph_cors <- rval2 %>%graph_from_data_frame(directed = FALSE)
	p4 <- ggraph(graph_cors) +geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +guides(edge_alpha = "none", edge_width = "none") +scale_edge_colour_viridis() +geom_node_point(color = "black", size = 3) +geom_node_text(aes(label = name), repel = TRUE) +theme_graph(base_family="sans") +labs(title = paste0("RMCORR: ", type_receptor, " modules (BH adj p < 0.05 & R>0.3)"))
	plot(p4)
	dev.off() 
	
	pdf(paste0(outputdir, "/WithinModuleCorrelation_RMCORR_MethodCompare.pdf"), height=15, width=15)
	plot(plot_grid(p1, p2, p3, p4, ncol=2))
	dev.off()
}
