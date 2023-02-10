## Function to generate a pairwise correlation matrix 
## Lauren Overend and Rachael Bashford-Rogers
## lauren.overend@oriel.ox.ac.uk
## Jan 2022 

############ Where we correlate on the whole dataset! 
Get_subsample_corr_matrix <-function(features, p, min_sample_input, outputdir, iso_type, type, method=NA){
  # get subsampled corr matrix
  mat_feature_corr = matrix(data = 0, nrow = length(features), ncol = length(features), dimnames = c(list(features), list(features)))
  repeats = 10000
  for(i1 in c(1:length(features))){
    print (i1)
    for(i2 in c(i1:length(features))){
      if(i1<i2){
        x = p[,features[i1]]
        y = p[,features[i2]]
        w = intersect(which(is.na(x)==F), which(is.na(y)==F))
        if(length(w)>min_sample_input){
          sample_corrs = NULL
          for (r in c(1:repeats)){
            rand = sample(w, min_sample_input)
            cor = cor(x[rand],y[rand])
            sample_corrs = c(sample_corrs, cor)
          }
          sample_corrs = sample_corrs[which(is.na(sample_corrs)==F)]
          mat_feature_corr[i1,i2] = mean(sample_corrs)
          mat_feature_corr[i2,i1] = mean(sample_corrs)
        }else{
          if(length(w)==min_sample_input){
            cor = cor(x[w],y[w])
            mat_feature_corr[i1,i2] = cor
            mat_feature_corr[i2,i1] = cor
          }
        }
      }
    }
  }
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