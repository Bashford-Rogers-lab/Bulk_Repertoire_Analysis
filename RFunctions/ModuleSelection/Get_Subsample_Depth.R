## Function to get the subsample depth for calculating correlation matrices 
## Lauren Overend and Rachael Bashford-Rogers
## lauren.overend@oriel.ox.ac.uk
## Jan 2022 
#filtered_matrix <- mat_filtered
#features <- features
#day <- 5

get_subsample_depth <- function(filtered_matrix, features){
			p <- filtered_matrix
			features <- features
			min_sample_input = 10^10
			for(i1 in c(1:length(features))){
			  for(i2 in c(i1:length(features))){
				if(i1<i2){
				  x = p[,features[i1]]
				  y = p[,features[i2]]
				  w = intersect(which(is.na(x)==F), which(is.na(y)==F))
				  if(length(w)>5){
					if(length(w)<min_sample_input){
					  min_sample_input = length(w)
					}
			}}}}
			## Multiply by 0.8 to ensure every sample isn't identical 
			min_sample_input = floor(0.8*min_sample_input)
			print(paste0("Sample Depth for Correlation Calculation: ", min_sample_input)) 
			return(min_sample_input)
}


get_subsample_depth_multiplesamples <- function(filtered_matrix, features, day){
			p <- filtered_matrix
			p$Sample <- rownames(p)
			for(x in 1:length(p$Sample)){
				p$Day[x] <- str_split_fixed(p$Sample[x], "_", 3)[,2]
				if(p$Sample[x] %like% "HV"){
				 p$Day[x] <- str_split_fixed(p$Sample[x], "_", 4)[,3]
				 }
			}
			###  We just want to use the first timpoint 
			p <- p[p$Day==day,]
			p$Sample <- NULL
			p$Day <- NULL
			
			########################
			features <- features
			min_sample_input = 10^10
			for(i1 in c(1:length(features))){
			  for(i2 in c(i1:length(features))){
				if(i1<i2){
				  x = p[,features[i1]]
				  y = p[,features[i2]]
				  w = intersect(which(is.na(x)==F), which(is.na(y)==F))
				  if(length(w)>5){
					if(length(w)<min_sample_input){
					  min_sample_input = length(w)
					  badcols <- c(features[i1],features[i2])
					}
			}}}}
			## Multiply by 0.8 to ensure every sample isn't identical 
			min_sample_input = floor(0.8*min_sample_input)
			print(paste0("Sample Depth for Correlation Calculation: ", min_sample_input)) 
			return(min_sample_input)
}