#' split the training data into M groups
#' use k-means to find k subcentroids for each group and
#' assign cluster numbers to each member of the dataset

index_cs <- function(dat, M=100, k=NULL) {
  if (is.null(k)) {
    k <- floor(sqrt(dim(dat)[2]))
  }
  new <- scmap::getFeatures(dat, n_features = 1000)
  indices <- which(rowData(new)$scmap_features)
  rownames(dat) <- rowData(dat)$feature_symbol
  dat <- logcounts(dat)[indices,]
  rows <- as.integer(dim(dat)[1])
  cols <- as.integer(dim(dat)[2])
  
  #normalize dataset to perform k-means by cosine similarity
  norm_dat <- normalise(dat)
  chunksize <- floor(rows/M)
  
  chunks <- list()
  for (i in 1:M) {
    chunks[[i]] <- norm_dat[((i-1)*chunksize+1):(i*chunksize),]
  }
  
  subcentroids <- list()
  subclusters <- list()
  km <- list()
  error <- 0
  for (m in 1:M) {
    tryCatch({
      # invisible(capture.output(km[[m]] <- kmeanscpp(chunks[[m]], k)))
      # subcentroids[[m]] <- as.matrix(km[[m]]$centers)
      # subclusters[[m]] <- as.vector(km[[m]]$result)
      km[[m]] <- kmeans(x = t(chunks[[m]]), centers = k, iter.max = 50)
      subclusters[[m]] <- km[[m]]$cluster
      subcentroids[[m]] <- t(km[[m]]$centers)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  subcentroids <- lapply(subcentroids, normalise)
  for (m in 1:M) {
    rownames(subcentroids[[m]]) <- rownames(dat[((m-1)*chunksize+1):(m*chunksize),])
  }
  subclusters <- do.call(rbind, subclusters)
  l <- list()
  l[[1]] <- subcentroids
  l[[2]] <- subclusters
  l[[3]] <- k
  l[[4]] <- M
  return(l)
}