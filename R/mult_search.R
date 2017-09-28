mult_search <- function(list_dat, query_dat, w) {
  features_query <- rowData(query_dat)$feature_symbol
  exprs_query <- logcounts(query_dat)
  num_cells <- dim(exprs_query)[2]
  rownames(exprs_query) <- rowData(query_dat)$feature_symbol
  print(sprintf("Searching reference dataset %i...", 1))
  # evaluate first dataset on the reference list
  ref <- list_dat[[1]]
  subcentroids <- ref[[1]]
  subclusters <- ref[[2]]
  query_chunks <- list()
  k <- ref[[3]]
  M <- ref[[4]]
  SqNorm <- numeric(dim(exprs_query)[2])

  for (m in 1:M) {
    subcentroids_chunk <- subcentroids[[m]]
    features_chunk <- rownames(subcentroids_chunk)
    common_features <- intersect(features_chunk, features_query)
    if (length(common_features)==0) {
      # change to a more memory-efficient method later?
      query_chunks[[m]] <- matrix(numeric(num_cells), 1, num_cells)
      subcentroids[[m]] <- matrix(numeric(k), 1, k)
    } else {
      # reduce both the query and the reference to just the rows corr. to common features
      common_features <- sort(common_features)
      subcentroids[[m]] <- subcentroids_chunk[rownames(subcentroids_chunk) %in% common_features, , drop = FALSE]
      query_chunks[[m]] <- exprs_query[rownames(exprs_query) %in% common_features, , drop = FALSE]
      if (length(common_features) > 1) {
        subcentroids[[m]] <- subcentroids[[m]][order(rownames(subcentroids[[m]])), ]
        query_chunks[[m]] <- query_chunks[[m]][order(rownames(query_chunks[[m]])), ]
      }
      # find the squared Euclidean norm of every query after selecting features
      SqNorm <- SqNorm + EuclSqNorm(query_chunks[[m]])
    }
  }
  #print("Feature preparation done!")
  res = NNfirst(w, k, subcentroids, subclusters, query_chunks, M, SqNorm)

  # evaluate the other datasets in the reference list
  if (length(list_dat)>1) {
    for (dat_num in 2:length(list_dat)) {
      print(sprintf("Searching reference dataset %i...", dat_num))
      ref <- list_dat[[dat_num]]
      subcentroids <- ref[[1]]
      subclusters <- ref[[2]]
      query_chunks <- list()
      k <- ref[[3]]
      M <- ref[[4]]
      SqNorm <- numeric(dim(exprs_query)[2])
      for (m in 1:M) {
        subcentroids_chunk <- subcentroids[[m]]
        features_chunk <- rownames(subcentroids_chunk)
        common_features <- intersect(features_chunk, features_query)

        if (length(common_features)==0) {
          # change to a more memory-efficient method later?
          query_chunks[[m]] <- matrix(numeric(num_cells), 1, num_cells)
          subcentroids[[m]] <- matrix(numeric(k), 1, k)
        } else {
          common_features <- sort(common_features)
          subcentroids[[m]] <- subcentroids_chunk[rownames(subcentroids_chunk) %in% common_features, , drop = FALSE]
          query_chunks[[m]] <- exprs_query[rownames(exprs_query) %in% common_features, , drop = FALSE]
          if (length(common_features) > 1) {
            subcentroids[[m]] <- subcentroids[[m]][order(rownames(subcentroids[[m]])), ]
            query_chunks[[m]] <- query_chunks[[m]][order(rownames(query_chunks[[m]])), ]
          }
          # find the squared Euclidean norm of every query after selecting features
          SqNorm <- SqNorm + EuclSqNorm(query_chunks[[m]])
        }
      }
      res = NNmult(w, k, subcentroids, subclusters, query_chunks, M, SqNorm, res$cells, res$distances, res$dataset_inds, dat_num)
    }
  }
  for (i in 1:3) {
    colnames(res[[i]]) <- colnames(exprs_query)
  }
  return(res)
}