## Function to actually get the correlation 
foo <- function(x1, y1, w1, min_sample_input1 ){
    rand <- sample(w1, min_sample_input1, replace = FALSE)
    corx = cor(x1[rand],y1[rand])
	corx
}

## This is a new fast version of the same function!
Get_subsample_corr_matrix_fast <-function(features, p, min_sample_input, outputdir, iso_type, type, no_repeats=10, method=NA){
  # get subsampled corr matrix
  mat_feature_corr = matrix(data = 0, nrow = length(features), ncol = length(features), dimnames = c(list(features), list(features)))
  replicates <- no_repeats
  print(paste0("Running ", replicates, " Repeats"))
  ## New faster version
  for(i1 in c(1:length(features))){
    print (i1)
    for(i2 in c(i1:length(features))){
        x = p[,features[i1]]
        y = p[,features[i2]]
        w = intersect(which(is.na(p[,features[i1]])==F), which(is.na(y)==F))
        sample_corrs <- replicate(replicates, foo(x,y,w,min_sample_input))
		sample_corrs = sample_corrs[which(is.na(sample_corrs)==F)]
		sample_corrs <- mean(sample_corrs)    
		mat_feature_corr[i1,i2] = sample_corrs
        mat_feature_corr[i2,i1] = sample_corrs
      }
    }
 
  ## Sorting and saving 
  sort(apply(mat_feature_corr, 1, function(x){length(which(is.na(x)))}))
  if(!is.na(method)){
		if(method=="Day1"){
		saveRDS(file = paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_DAY1ONLY_", type, "_", iso_type, ".rds"), mat_feature_corr)
		}
		if(method=="Day3"){
		saveRDS(file = paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_DAY3ONLY_", type, "_", iso_type, ".rds"), mat_feature_corr)
		}
		if(method=="Day5"){
		saveRDS(file = paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_DAY5ONLY_", type, "_", iso_type, ".rds"), mat_feature_corr)
		}
		
  } else {
		saveRDS(file = paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_", type, "_", iso_type, ".rds"), mat_feature_corr)
  }
  
  return(mat_feature_corr)
}

