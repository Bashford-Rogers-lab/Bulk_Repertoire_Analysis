library(data.table) # for fread
library(fields)

read_output <- function(folder) {
  # reads all estimates 
  # e.g folder <- "my_output_folder/it2000"
  # if you are in the folder containing the estimates, folder <- "../it2000"
  out <- list()
  folder_exists <- system(paste0("[ -d ", folder, " ] && echo true || echo false"), intern=TRUE)

  stopifnot(folder_exists=='true')

  # get list of all files in folder
  files <- list.files(folder)
  # read all files
  for (file in files) {
    out[[file]] <- as.matrix(fread(paste0(folder, "/", file)))
  }

  out1 <- reformat_data(out) # from the estimates, get N, L, P, ncomps and reformat data
  for (i in names(out1)) {
    out[[i]] <- out1[[i]]
  }
  out$X1 <- NULL; out$S1 <- NULL; out$B1 <- NULL 

  return(out)
}

reformat_data <- function(out) {
  nam <- names(out)

  est <- list()
  est$A <- out$A

  est$N <- nrow(est$A) # number individuals
  est$C <- ncol(est$A) # number components

  num_X_mats <- length(nam[grep("X[0-9]?", nam)])
  stopifnot(length(nam[grep("S[0-9]?", nam)])==num_X_mats)
  stopifnot(num_X_mats != 0)
  num_B_mats <- length(nam[grep("B[0-9]?", nam)])

  if (num_B_mats == 0 && num_X_mats == 1) { # 2D matrix factorisation
    est$X <- out$X1
    est$S <- out$S1

  } else if (num_B_mats == 1 && num_X_mats == 1) { # 3D matrix factorisation
    est$X <- out$X1; est$S <- out$S1
    est$B <- out$B1
    est$P <- nrow(est$B)
  
  } else {
    print("problem with estimates")
  }

  return(est)
}

plot_scores <- function(data, out) {
  par(mfrow=c(1,3), mgp=c(2,1,0), mar=c(4,4,4,4), oma=c(1,1,1,1))
  image.plot(t(data$A), 
      main="Individual scores (truth)", ylab="Individuals", xlab="Components")
  image.plot(t(out$A), 
      main="Individual scores (estimates)", ylab="Individuals", xlab="Components")
  image.plot(abs(cor(data$A, out$A)), 
      main="Correlation between\ntrue and estimated\nindividual scores", 
      xlab="Truth", ylab="Estimates")
}
