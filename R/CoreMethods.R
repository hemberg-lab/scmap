#' Find the most informative features (genes/transcripts) for projection
#'
#' This is a modification of the M3Drop method. Instead of fitting a 
#' Michaelis-Menten model to the log expression-dropout relation, we fit a 
#' linear model. Namely, the linear model is build on the log(expression) versus 
#' log(dropout) distribution. After fitting a linear model important features are
#'  selected as the top N residuals of the linear model.
#' 
#' Please note that \code{feature_symbol} column of \code{rowData(object)} must be 
#' present in the input object and should not contain any duplicated feature names. 
#' This column defines feature names used during projection. Feature symbols 
#' in the reference dataset must correpond to the feature symbols
#' in the projection dataset, otherwise the mapping will not work!
#'
#' @param object an object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class
#' @param n_features number of the features to be selected
#' @param suppress_plot boolean parameter, which defines whether to plot 
#' log(expression) versus log(dropout) distribution for all genes.
#' Selected features are highlighted with the red colour.
#'
#' @return an object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class with a new column in 
#' \code{rowData(object)} slot which is called \code{scmap_features}. It can be accessed
#' by using \code{as.data.frame(rowData(object))$scmap_features}.
#'
#' @name selectFeatures
#'
#' @importFrom SummarizedExperiment rowData rowData<-
selectFeatures.SingleCellExperiment <- function(object, n_features, suppress_plot) {
    if (is.null(as.data.frame(rowData(object))$feature_symbol)) {
        stop("There is no feature_symbol column in the rowData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
        return(object)
    }
    
    r_data <- as.data.frame(rowData(object))
    tmp <- linearModel(object, n_features)
    r_data$scmap_features <- tmp$scmap_features
    r_data$scmap_scores <- tmp$scmap_scores
    rowData(object) <- r_data
    
    if (!suppress_plot) {
        p <- ggplot_features(tmp$for_plotting, tmp$fit)
        print(p)
    }
    
    return(object)
}

#' @rdname selectFeatures
#' @aliases selectFeatures
setMethod("selectFeatures", "SingleCellExperiment", selectFeatures.SingleCellExperiment)

#' Set the most important features (genes/transcripts) for projection
#' 
#' This method manually sets the features to be used for projection.
#' 
#' Please note that \code{feature_symbol} column of \code{rowData(object)} must be 
#' present in the input object and should not contain any duplicated feature names. 
#' This column defines feature names used during projection. Feature symbols 
#' in the reference dataset must correpond to the feature symbols
#' in the projection dataset, otherwise the mapping will not work!
#'
#' @param object an object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class
#' @param features a character vector of feature names
#' 
#' @return an object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class with a new column in 
#' \code{rowData(object)} slot which is called \code{scmap_features}. It can be accessed
#' by using \code{as.data.frame(rowData(object))$scmap_features}.
#'
#' @name setFeatures
#'
#' @importFrom SummarizedExperiment rowData rowData<-
setFeatures.SingleCellExperiment <- function(object, features) {
    if (is.null(features)) {
        stop("Please provide a list of feature names using 'features' argument!")
        return(object)
    }
    if (is.null(as.data.frame(rowData(object))$feature_symbol)) {
        stop("There is no feature_symbol column in the rowData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
        return(object)
    }
    
    inds <- match(features, as.data.frame(rowData(object))$feature_symbol)
    
    if (!all(!is.na(inds))) {
        warning("Features ", paste(features[which(is.na(inds))], collapse = ", "), " are not present in the 'SCESet' object and therefore were not set.")
    }
    
    r_data <- as.data.frame(rowData(object))
    r_data$scmap_features <- FALSE
    r_data$scmap_features[inds[!is.na(inds)]] <- TRUE
    rowData(object) <- r_data
    
    return(object)
}

#' @rdname setFeatures
#' @aliases setFeatures
setMethod("setFeatures", "SingleCellExperiment", setFeatures.SingleCellExperiment)

#' Create a precomputed Reference
#' 
#' Calculates centroids of each cell type and merge them into a single table.
#'
#' @param object SingleCellExperiment object
#' @param cluster_col column name in the `colData` slot of the SingleCellExperiment object 
#' containing the cell classification information
#' 
#' @name indexCluster
#'
#' @return a `data.frame` containing calculated centroids of the cell types of
#' the Reference dataset
#'
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom dplyr group_by summarise %>%
#' @importFrom reshape2 melt dcast
#' @importFrom stats median
indexCluster.SingleCellExperiment <- function(object, cluster_col) {
    if(!checks_for_index(object)) {
      return(object)
    }
    if (!cluster_col %in% colnames(colData(object))) {
        stop("Please define an existing cluster column of the `colData` slot of the input object using the `cluster_col` parameter!")
        return(object)
    }
    tmp <- object[rowData(object)$scmap_features, ]
    gene <- cell_class <- exprs <- NULL
    exprs_mat <- logcounts(tmp)
    rownames(exprs_mat) <- as.data.frame(rowData(tmp))$feature_symbol
    colnames(exprs_mat) <- as.data.frame(colData(tmp))[[cluster_col]]
    
    # calculate median feature expression in every cell class of object
    exprs_mat <- reshape2::melt(exprs_mat)
    colnames(exprs_mat) <- c("gene", "cell_class", "exprs")
    exprs_mat <- exprs_mat %>% group_by(gene, cell_class) %>% summarise(med_exprs = median(exprs))
    exprs_mat <- reshape2::dcast(exprs_mat, gene ~ cell_class, value.var = "med_exprs")
    rownames(exprs_mat) <- exprs_mat$gene
    index <- exprs_mat[, 2:ncol(exprs_mat), drop = FALSE]
    index <- index[order(rownames(index)), , drop = FALSE]
    index <- index[, colSums(index) > 0, drop = FALSE]
    if (ncol(index) == 0) {
      stop("scmap index is empty because the median expression in the selected features is 0 in every cell cluster! Try to increase the number of selected features!")
      return(object)
    }
    metadata(object)$scmap_cluster_index <- index
    return(object)
}

#' @rdname indexCluster
#' @aliases indexCluster
setMethod("indexCluster", "SingleCellExperiment", indexCluster.SingleCellExperiment)

#' Create an index for a dataset to enable fast approximate nearest neighbour search
#' 
#' The method is based on product quantization for the cosine distance.
#' Split the training data into M identically sized chunks by genes.
#' Use k-means to find k subcentroids for each group.
#' Assign cluster numbers to each member of the dataset.
#' 
#' @param object an object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class
#' @param M number of chunks into which the expr matrix is split
#' @param k number of clusters per group for k-means clustering
#' 
#' @return a list of four objects: 1) a list of matrices containing the subcentroids of each group
#' 2) a matrix containing the subclusters for each cell for each group
#' 3) the value of M
#' 4) the value of k
#' 
#' @name indexCell
#' 
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom stats kmeans
#' 
#' @useDynLib scmap
#' @importFrom Rcpp sourceCpp
indexCell.SingleCellExperiment <- function(object, M, k) {
  if(!checks_for_index(object)) {
    return(object)
  }
  n_features <- length(which(rowData(object)$scmap_features))
  if (is.null(M)) {
    message("Parameter M was not provided, will use M = n_features / 10 (if n_features <= 1000), where n_features is the number of selected features, and M = 100 otherwise.")
    if (n_features > 1000) {
      M <- 100
    } else {
      M <- n_features / 10
    }
  }
  # if k is unspecified, we assign it to be the sqrt of the number of cells in the dataset
  if (is.null(k)) {
    message("Parameter k was not provided, will use k = sqrt(number_of_cells)")
    k <- floor(sqrt(ncol(object)))
  }
  
  tmp <- object[rowData(object)$scmap_features, ]
  exprs_mat <- logcounts(tmp)
  rownames(exprs_mat) <- as.data.frame(rowData(tmp))$feature_symbol
  features <- rownames(exprs_mat)
  
  # normalize dataset to perform k-means by cosine similarity
  norm_dat <- normalise(exprs_mat)
  colnames(norm_dat) <- colnames(exprs_mat)
  rownames(norm_dat) <- rownames(exprs_mat)
  chunksize <- floor(nrow(norm_dat)/M)
  
  # make chunks of the data
  chunks <- list()
  for (i in seq_len(M)) {
    chunks[[i]] <- norm_dat[((i - 1) * chunksize + 1):(i * chunksize), ]
  }
  
  subcentroids <- list()
  subclusters <- list()
  km <- list()
  
  # perform k-means for each chunk
  for (m in seq_len(M)) {
    tryCatch({
      km[[m]] <- kmeans(x = t(chunks[[m]]), centers = k, iter.max = 50)
      subclusters[[m]] <- km[[m]]$cluster
      subcentroids[[m]] <- t(km[[m]]$centers)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  # find chunks where k-means failed
  inds <- which(!unlist(lapply(subcentroids, is.null)))
  subcentroids <- subcentroids[inds]
  subclusters <- do.call(rbind, subclusters)
  
  # normalise chunks
  subcentroids <- lapply(
    subcentroids, 
    function(x) {
      y <- normalise(x)
      rownames(y) <- rownames(x)
      colnames(y) <- colnames(x)
      return(y)
    }
  )

  index <- list()
  index$subcentroids <- subcentroids
  index$subclusters <- subclusters
  metadata(object)$scmap_cell_index <- index
  return(object)
}

#' @rdname indexCell
#' @aliases indexCell
setMethod("indexCell", "SingleCellExperiment", indexCell.SingleCellExperiment)

#' scmap main function
#' 
#' Projection of one dataset to another
#' 
#' @param projection `SingleCellExperiment` object to project
#' @param index_list list of index objects each coming from the output of `indexCluster`
#' @param threshold threshold on similarity (or probability for SVM and RF)
#' 
#' @return The projection object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class with labels calculated by `scmap` and stored in 
#' the \code{scmap_labels} column of the \code{rowData(object)} slot.
#' 
#' @name scmapCluster
#' 
#' @importFrom SummarizedExperiment rowData colData colData<-
#' @importFrom S4Vectors DataFrame
#' @importFrom proxy simil
#' @importFrom stats cor
#' @importFrom matrixStats colMaxs rowMaxs
#' @importFrom methods is
scmapCluster.SingleCellExperiment <- function(projection, index_list, threshold) {
  if(!checks_for_projection(projection, index_list)) {
    return(projection)
  }
  labels <- list()
  simls <- list()
  for (n in seq_len(length(index_list))) {
    index <- index_list[[n]]
    # find and select only common features, then subset both datasets
    tmp <- setFeatures(projection, rownames(index))
    index <- index[rownames(index) %in% rowData(tmp)$feature_symbol[rowData(tmp)$scmap_features], , drop = FALSE]
    tmp <- tmp[rowData(tmp)$scmap_features, ]
    
    if (nrow(index) < 10) {
      warning("There are less than ten features in common between the `reference` and `projection` datasets. Most probably they come from different organisms! Please redefine your query!")
      return(projection)
    }  
    
    # get expression values of the projection dataset
    proj_exprs <- logcounts(tmp)
    rownames(proj_exprs) <- rowData(tmp)$feature_symbol
    
    # prepare projection dataset
    proj_exprs <- proj_exprs[order(rownames(proj_exprs)), ]
    
    # calculate similarities and correlations
    tmp <- t(index)
    res <- proxy::simil(tmp, t(proj_exprs), method = "cosine")
    res <- matrix(res, ncol = nrow(tmp), byrow = TRUE)
    max_inds1 <- max.col(res)
    maxs1 <- rowMaxs(res)
    
    res <- cor(index, proj_exprs, method = "pearson")
    max_inds2 <- max.col(t(res))
    maxs2 <- colMaxs(res)
    
    res <- cor(index, proj_exprs, method = "spearman")
    max_inds3 <- max.col(t(res))
    maxs3 <- colMaxs(res)
    
    cons <- cbind(colnames(index)[max_inds1], colnames(index)[max_inds2], 
                  colnames(index)[max_inds3])
    
    maximums <- cbind(maxs1, maxs2, maxs3)
    
    # cells with at least one NA correlation value become unassigned
    non_na_inds <- !is.na(max_inds1) & !is.na(max_inds2) & !is.na(max_inds3)
    
    # create labels
    maxs <- rep(NA, nrow(cons))
    labs <- rep("unassigned", nrow(cons))
    unique_labs <- unlist(apply(cons, 1, function(x) {
      length(unique(x))
    }))
    
    ## all similarities agree
    if (length(which(unique_labs == 1 & non_na_inds)) > 0) {
      labs[unique_labs == 1 & non_na_inds] <- cons[unique_labs == 1 & non_na_inds, 1]
      maxs_tmp <- rowMaxs(maximums[unique_labs == 1 & non_na_inds, , drop = FALSE])
      maxs[unique_labs == 1 & non_na_inds] <- maxs_tmp
    }
    
    ## only two similarities agree
    if (length(which(unique_labs == 2 & non_na_inds)) > 0) {
      tmp <- cons[unique_labs == 2 & non_na_inds, , drop = FALSE]
      inds <- unlist(apply(tmp, 1, function(x) {
        which(duplicated(x))
      }))
      labs[unique_labs == 2 & non_na_inds] <- tmp[cbind(seq_along(inds), inds)]
      
      ## calculate maximum similarity in case of two agreeing similarities
      tmp1 <- matrix(apply(tmp, 2, `==`, labs[unique_labs == 2 & non_na_inds]), ncol = 3)
      inds <- t(apply(tmp1, 1, which))
      maxs_tmp <- cbind(maximums[unique_labs == 2 & non_na_inds, , drop = FALSE][cbind(seq_along(inds[, 
                                                                                                      1]), inds[, 1])], maximums[unique_labs == 2 & non_na_inds, , drop = FALSE][cbind(seq_along(inds[, 
                                                                                                                                                                                                      1]), inds[, 2])])
      maxs_tmp <- rowMaxs(maxs_tmp)
      maxs[unique_labs == 2 & non_na_inds] <- maxs_tmp
    }
    
    ## check the similarity threshold
    labs[!is.na(maxs) & maxs < threshold] <- "unassigned"
    labels[[n]] <- labs
    simls[[n]] <- maxs
  }
  names(labels) <- names(index_list)
  names(simls) <- names(index_list)
  unassigned_rate_order <- order(
    unlist(
      lapply(labels, function(x) {
        length(x[x == "unassigned"])/length(x)
      })
    )
  )
  labels <- labels[unassigned_rate_order]
  simls <- simls[unassigned_rate_order]
  labels <- do.call(cbind, labels)
  simls <- do.call(cbind, simls)
  max_simls_inds <- apply(simls, 1, which.max)
  cons_labels <- labels[cbind(seq_along(max_simls_inds), max_simls_inds)]
  
  return(list(scmap_cluster_labs = labels, scmap_cluster_siml = simls, combined_labs = cons_labels))
}

#' @rdname scmapCluster
#' @aliases scmapCluster
setMethod("scmapCluster", "SingleCellExperiment", scmapCluster.SingleCellExperiment)

#' For each cell in a query dataset, we search for the nearest neighbours by cosine distance
#' within a collection of reference datasets.
#' 
#' @param projection an object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class
#' @param index_list list of index objects each coming from the output of `indexCell`
#' @param w a positive integer specifying the number of nearest neighbours to find
#' 
#' @return a list of 3 objects: 
#' 1) a matrix with the closest w neighbours by cell number of each query cell stored by column
#' 2) a matrix of integers giving the reference datasets from which the above cells came from
#' 3) a matrix with the cosine similarities corresponding to each of the nearest neighbours
#' 
#' @name scmapCell
#' 
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' 
#' @useDynLib scmap
#' @importFrom Rcpp sourceCpp
scmapCell.SingleCellExperiment <- function(projection, index_list, w) {
  if(!checks_for_projection(projection, index_list)) {
    return(projection)
  }
  res_all <- list()
  for (n in seq_len(length(index_list))) {
    index <- index_list[[n]]
    subcentroids <- index$subcentroids
    subclusters <- index$subclusters
    
    proj_exprs <- logcounts(projection)
    rownames(proj_exprs) <- rowData(projection)$feature_symbol
    
    dists <- dists_subcentroids(proj_exprs, subcentroids)
    
    # compute the w nearest neighbours and their similarities to the queries
    res <- NN(w, ncol(subcentroids[[1]]), dists$subcentroids, subclusters, dists$query_chunks, length(subcentroids), dists$SqNorm)
    for (i in 1:2) {
      colnames(res[[i]]) <- colnames(proj_exprs)
    }
    res_all[[n]] <- res
  }
  names(res_all) <- names(index_list)
  return(res_all)
}

#' @rdname scmapCell
#' @aliases scmapCell
setMethod("scmapCell", "SingleCellExperiment", scmapCell.SingleCellExperiment)

#' Approximate k-NN cell-type classification using scfinemap
#' 
#' Each cell in the query dataset is assigned a cell-type if the similarity between its
#' nearest neighbour exceeds a threshold AND its w nearest neighbours have the 
#' same cell-type. 
#' 
#' @param projection an object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class
#' @param scmapCell_results the output of `scmapCell()` with `projection` as its input.
#' @param cluster_list list of cell cluster labels correspondint to each index against which the `projection` has been projected 
#' @param w an integer specifying the number of nearest neighbours to find
#' @param threshold the threshold which the maximum similarity between the query and a reference cell must exceed
#' for the cell-type to be assigned
#' 
#' @return The query dataset with the predicted labels attached to colData(query_dat)$cell_type1
#' 
#' @name scmapCell2Cluster
#' 
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
#' 
#' @useDynLib scmap
#' @importFrom Rcpp sourceCpp
scmapCell2Cluster.SingleCellExperiment <- function(projection, scmapCell_results, cluster_list, w, threshold) {
  res <- list()
  for (i in seq_len(length(scmapCell_results))) {
    cells <- scmapCell_results[[i]]$cells
    distances <- scmapCell_results[[i]]$distances
    labs_new <- character(ncol(cells))
    for (j in seq_len(ncol(cells))) {
      celltypes <- character(w)
      for (k in 1:w) {
        celltypes[k] <- cluster_list[[i]][cells[k, j]]
      }
      # only assign the cell-type if the similarity exceeds the threshold and the w nearest neighbours
      # are the same cell-type
      if (max(distances[, j]) > threshold & length(unique(celltypes)) == 1) {
        labs_new[j] <- unique(celltypes)
      } else {
        labs_new[j] <- "unassigned"
      }
    }
    res[[i]] <- labs_new
  }
  names(res) <- names(scmapCell_results)
  return(res)
}

#' @rdname scmapCell2Cluster
#' @aliases scmapCell2Cluster
setMethod("scmapCell2Cluster", "SingleCellExperiment", scmapCell2Cluster.SingleCellExperiment)
