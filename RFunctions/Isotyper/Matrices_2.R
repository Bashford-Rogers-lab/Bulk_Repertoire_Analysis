make_matrices2 <- function(file=file, ids_all=ids_all){ 
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.delim(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.Id"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
			
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		###		
		
		c1 = c("Class_switched","IGHD,IGHM_mutated","IGHD,IGHM_unmutated")
		c2 = c( "IGHA1","IGHA2","IGHD","IGHE","IGHG1","IGHG2","IGHG3","IGHG4","IGHM"  )
		
		## Divide into smaller matrix 
		p_unique <- p[p$Isotype %in% c2,]
		p_group <- p[p$Isotype %in% c1,]
			
		## Per isotype 
		p_unique$total_reads_per_isotype <-  NA
		for(i in 1:length(p_unique$X.Id)){
			sample <- p_unique$X.Id[i]
			p_unique$total_reads_per_isotype[i] <- sum(p_unique$N.reads[p_unique$X.Id==p_unique$X.Id[i] & p_unique$Isotype %in% c2])
		}
		p_unique$total_Unique_reads_per_isotype <-  NA
		for(i in 1:length(p_unique$X.Id)){
			sample <- p_unique$X.Id[i]
			p_unique$total_Unique_reads_per_isotype[i] <- sum(p_unique$N.vertices[p_unique$X.Id==p_unique$X.Id[i] & p_unique$Isotype %in% c2])
		}
		
		## Per Class 
		p_group$total_reads_per_class <-  NA
		for(i in 1:length(p_group$X.Id)){
			sample <- p_group$X.Id[i]
			p_group$total_reads_per_class[i] <- sum(p_group$N.reads[p_group$X.Id==p_group$X.Id[i] & p_group$Isotype %in% c1])
		}
		p_group$total_Unique_reads_per_class <-  NA
		for(i in 1:length(p_group$X.Id)){
			sample <- p_group$X.Id[i]
			p_group$total_Unique_reads_per_class[i] <- sum(p_group$N.vertices[p_group$X.Id==p_group$X.Id[i] & p_group$Isotype %in% c1])
		}
		
		#####
		## Get percentages
		p_unique$percentage_total_vdjs_per_isotype <- (p_unique$N.reads / p_unique$total_reads_per_isotype)*100
		p_unique$percentage_unique_vdjs_per_isotype <- (p_unique$N.vertices / p_unique$total_Unique_reads_per_isotype)*100
		p_group$percentage_total_vdjs_per_class <- (p_group$N.reads / p_group$total_reads_per_class)*100
		p_group$percentage_unique_vdjs_per_class  <- (p_group$N.vertices / p_group$total_Unique_reads_per_class)*100
	 
		## Make into wide format 
		p_unique <- p_unique[, c("X.Id", "Isotype", "percentage_total_vdjs_per_isotype", "percentage_unique_vdjs_per_isotype")]
		l <- reshape(p_unique, idvar = "X.Id", timevar = "Isotype", direction = "wide")
		colnames(l) <- gsub("\\.", "__", colnames(l))
		rownames(l) <- l$X__Id
		l$X__Id <- NULL
		l[is.na(l)] <- 0 
		
		p_group <- p_group[, c("X.Id", "Isotype", "percentage_total_vdjs_per_class", "percentage_unique_vdjs_per_class")]
		l2 <- reshape(p_group, idvar = "X.Id", timevar = "Isotype", direction = "wide")
		colnames(l2) <- gsub("\\.", "__", colnames(l2))
		rownames(l2) <- l2$X__Id
		l2$X__Id <- NULL
		l2[is.na(l2)] <- 0
		
		## Bind them together 
		new <- merge(l, l2, by=0)
		rownames(new) <- new$Row.names
		new$Row.names <- NULL
		
		## make into matrix!
		q <- list(new)
		names(q) <- "Percentages"
		
		analysis_matrices1 = q	
		## Check for -1 values which need to be replaced with NA 
		for(i in 1:length(analysis_matrices1)){
			analysis_matrices1[[i]][analysis_matrices1[[i]]==-1 | analysis_matrices1[[i]]=="-1" | analysis_matrices1[[i]]=="-1.0" | analysis_matrices1[[i]]==-1.0] <- NA
		}
		## if there is no file
		} else { 
				print(paste0("File: ", file, " IS EMPTY"))
				analysis_matrices1 <- vector(mode = "list", length = 0)
		} 
	print("DONE 2")
	return(analysis_matrices1)
}	
