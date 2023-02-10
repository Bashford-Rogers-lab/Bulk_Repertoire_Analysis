## Auxillary Script Functions
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## Jan 2022


add.alpha <- function(col, alpha=1){
	  if(missing(col))
		stop("Please provide a vector of colours.")
	  apply(sapply(col, col2rgb)/255, 2, 
			function(x) 
			  rgb(x[1], x[2], x[3], alpha=alpha)) }
concat = function(v) {
	  res = ""
	  for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	  res
	}