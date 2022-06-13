make_matrices9 <- function(file=file,ids_all=ids_all){
	info = file.info(file)
	if(info$size != 0) {
	
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		p=p[which(as.numeric(p[,"total_all"]) >10),]
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)[1]=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[,  c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		############################
		## Calculate the meaures 
		p$V4_34_AVY_NHS_total_unmut_count <- as.numeric(p$V4_34_AVY_NHS_total_unmut)
		p$V4_34_AVY_NHS_total_unmut <- as.numeric(p$V4_34_AVY_NHS_total_unmut_count)*100/as.numeric(p$total_all)
		p$V4_34_AVY_total_unmut <- as.numeric(p$V4_34_AVY_total_unmut)*100/as.numeric(p$total_all)
		p$V4_34_NHS_total_unmut <- as.numeric(p$V4_34_NHS_total_unmut)*100/as.numeric(p$total_all)
		
		d <- p[, c("X.sample", "isotype", "V4_34_AVY_total_unmut", "V4_34_NHS_total_unmut")]
		colnames(d) <- c("X.sample", "isotype", "Percentage_V4_34_AVY_unmutated", "Percentage_V4_34_NHS_unmutated")
	    e <- p[, c("X.sample", "isotype", "V4_34_AVY_NHS_total_unmut")]
		
		## classify 
		class_switched <- c("IGHA1", "IGHA2" , "IGHE",  "IGHG1" ,"IGHG2" ,"IGHG3" ,"IGHG4")
		IGHMD <- c("IGHD", "IGHM")
		## Calculate proportions for classified isotypes!! 
		p_switched <- p[p$isotype %in% class_switched,]
		df_switched <- data.frame()
		for(i in 1:length(unique(p_switched$X.sample))){
			samples <- unique(p_switched$X.sample)[i]
			total_reads <- sum(p_switched$total_all[p_switched$X.sample==samples])
			total_incidence <- sum(p_switched$V4_34_AVY_NHS_total_unmut_count[p_switched$X.sample==samples])
			percentage <- (total_incidence*100)/total_reads
			row <- c(samples, "class_switched", percentage)
			df_switched <- rbind(df_switched, row)
			}
		colnames(df_switched) <- colnames(e)	
		p_IGHMD <- p[p$isotype %in% IGHMD,]
		df_IGHMD <- data.frame()
		for(i in 1:length(unique(p_IGHMD$X.sample))){
			samples <- unique(p_IGHMD$X.sample)[i]
			total_reads <- sum(p_IGHMD$total_all[p_IGHMD$X.sample==samples])
			total_incidence <- sum(p_IGHMD$V4_34_AVY_NHS_total_unmut_count[p_IGHMD$X.sample==samples])
			percentage <- (total_incidence*100)/total_reads
			row <- c(samples, "IGHM/D", percentage)
			df_IGHMD <- rbind(df_IGHMD, row)
			}
		colnames(df_IGHMD) <- colnames(e)
		e <- rbind(e, df_switched, df_IGHMD)
		colnames(e) <- c("X.sample", "isotype", "Percentage_V4_34_AVY_NHS_unmutated")
		#################################	
		
		## reshape d
		colnames(d) <- gsub("\\.", "_", colnames(d))
		l <- reshape(d, idvar = "X_sample", timevar = "isotype", direction = "wide")
		colnames(l) <- gsub("\\.", "__", colnames(l))
		rownames(l) <- l$X_sample
		l$X_sample <- NULL
		
		## reshape e 
		colnames(e) <- gsub("\\.", "_", colnames(e))
		l2 <- reshape(e, idvar = "X_sample", timevar = "isotype", direction = "wide")
		colnames(l2) <- gsub("\\.", "__", colnames(l2))
		rownames(l2) <- l2$X_sample
		l2$X_sample <- NULL
		
		## Merge them both together
		new <- merge(l, l2, by=0)
		rownames(new) <- new$Row.names
		new$Row.names <- NULL
		
		q <- list(new)
		names(q) <- "auto_immunity"
		analysis_matrices = q	
		## Check for -1 values which need to be replaced with NA 
		for(i in 1:length(analysis_matrices)){
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		}
		## if there is no file
		} else { 
				print(paste0("File: ", file, " IS EMPTY"))
				analysis_matrices <- vector(mode = "list", length = 0)
		} 
		print("Done 9")
		return(analysis_matrices)
} 		