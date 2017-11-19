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

#' scmap main function
#' 
#' Projection of one dataset to another
#' 
#' @param projection SingleCellExperiment object to project
#' @param reference reference SingleCellExperiment object
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
scmapCluster.SingleCellExperiment <- function(projection, reference, threshold) {
    if (is.null(projection)) {
        stop("Please provide a `SingleCellExperiment` object for the `projection` parameter!")
        return(projection)
    }
    if (!"SingleCellExperiment" %in% is(projection)) {
        stop("`projection` dataset has to be of the `SingleCellExperiment` class!")
        return(projection)
    }
    if (is.null(rowData(projection)$feature_symbol)) {
        stop("There is no `feature_symbol` column in the `rowData` slot of the `projection` dataset! Please write your gene/transcript names to this column!")
        return(projection)
    }
    if (is.null(reference)) {
        stop("Please provide either a `SingleCellExperiment` object or precomputed `data.frame` reference for the `reference` parameter!")
        return(projection)
    }
    if ("SingleCellExperiment" %in% is(reference)) {
        if (is.null(metadata(reference)$scmap_cluster_index)) {
            stop("scmap index has not been calculated yet! Please run `indexCluster()` first!")
            return(projection)
        }
        index <- metadata(reference)$scmap_cluster_index
    } else {
        if (!"data.frame" %in% is(reference)) {
            stop("Your reference is neither of `SingleCellExperiment` class nor of `data.frame` class! Please define a correct reference!")
            return(projection)
        }
        index <- reference
    }
  
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
    
    # if (method == "svm") {
    #     # create the reference
    #     index <- logcounts(reference)
    #     rownames(index) <- as.data.frame(rowData(reference))$feature_symbol
    #     colnames(index) <- as.data.frame(colData(reference))[[cluster_col]]
    #     
    #     # prepare the datasets for projection
    #     tmp <- prepareData(index, dat)
    #     index <- tmp$reference
    #     dat <- tmp$dat
    #     
    #     # run support vector machine
    #     res <- support_vector_machines(index, dat)
    #     probs <- attr(res, "probabilities")
    #     max_inds <- max.col(probs)
    #     maxs <- rowMaxs(probs)
    #     
    #     # create labels
    #     labs <- rep("unassigned", length(max_inds))
    #     labs[maxs > threshold] <- colnames(probs)[max_inds][maxs > threshold]
    # }
    # 
    # if (method == "rf") {
    #     # create the reference
    #     index <- logcounts(reference)
    #     rownames(index) <- as.data.frame(rowData(reference))$feature_symbol
    #     colnames(index) <- as.data.frame(colData(reference))[[cluster_col]]
    #     
    #     # prepare the datasets for projection
    #     tmp <- prepareData(index, dat)
    #     index <- tmp$reference
    #     dat <- tmp$dat
    #     
    #     # run random forest
    #     res <- random_forest(index, dat)
    #     max_inds <- max.col(res)
    #     maxs <- rowMaxs(res)
    #     
    #     # create labels
    #     labs <- rep("unassigned", length(max_inds))
    #     labs[maxs > threshold] <- colnames(res)[max_inds][maxs > threshold]
    # }
    
    c_data <- as.data.frame(colData(projection))
    c_data$scmap_labs <- labs
    c_data$scmap_siml <- maxs
    colData(projection) <- DataFrame(c_data)
    
    return(projection)
}

#' @rdname scmapCluster
#' @aliases scmapCluster
setMethod("scmapCluster", "SingleCellExperiment", scmapCluster.SingleCellExperiment)

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
    if (is.null(object)) {
        stop("Please provide a `SingleCellExperiment` object using the `object` parameter!")
        return(object)
    }
    if (!"SingleCellExperiment" %in% is(object)) {
        stop("Input object is not of `SingleCellExperiment` class! Please provide an object of the correct class!")
        return(object)
    }
    if (is.null(rowData(object)$scmap_features)) {
        stop("Features are not selected! Please run `selectFeatures()` or `setFeatures()` first!")
        return(object)
    }
    if (is.null(rowData(object)$feature_symbol)) {
        stop("There is no `feature_symbol` column in the `rowData` slot of the `reference` dataset! Please write your gene/transcript names to this column!")
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
#' @param dat an object of \code{\link[SingleCellExperiment]{SingleCellExperiment}} class
#' @param k number of clusters per group for k-means clustering
#' @param M number of chunks into which the expr matrix is split
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
indexCell.SingleCellExperiment <- function(reference, M, k) {
  
  # if k is unspecified, we assign it to be the sqrt of the number of cells in the dataset
  if (is.null(k)) {
    message("Parameter k was not provided, will use k = sqrt(number_of_cells)")
    k <- floor(sqrt(ncol(reference)))
  }
  if(is.null(rowData(reference)$feature_symbol)) {
    stop("Please provide unique feature names in the `feature_symbol` column of the `rowData` slot!")
    return(NULL)
  }
  if(is.null(rowData(reference)$scmap_features)) {
    stop("Please select features first by running getFeatures()!")
    return(NULL)
  }
  rownames(reference) <- rowData(reference)$feature_symbol
  exprs_matrix <- logcounts(reference)[rowData(reference)$scmap_features, ]
  features <- rownames(exprs_matrix)
  
  # normalize dataset to perform k-means by cosine similarity
  norm_dat <- normalise(exprs_matrix)
  chunksize <- floor(nrow(exprs_matrix)/M)
  
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
      # invisible(capture.output(km[[m]] <- kmeanscpp(chunks[[m]], k))) subcentroids[[m]] <-
      # as.matrix(km[[m]]$centers) subclusters[[m]] <- as.vector(km[[m]]$result)
      km[[m]] <- kmeans(x = t(chunks[[m]]), centers = k, iter.max = 50)
      subclusters[[m]] <- km[[m]]$cluster
      subcentroids[[m]] <- t(km[[m]]$centers)
    }, error = function(e) {
      return(NULL)
    })
  }
  inds <- which(!unlist(lapply(subcentroids, is.null)))
  subcentroids <- lapply(subcentroids, function(x) {
    if (!is.null(x)) {
      return(normalise(x))
    } else {
      return(x)
    }
  }
  )
  for (m in seq_len(M)[inds]) {
    rownames(subcentroids[[m]]) <- features[((m - 1) * chunksize + 1):(m * chunksize)]
  }
  subclusters <- do.call(rbind, subclusters[inds])
  l <- list()
  l[[1]] <- subcentroids[inds]
  l[[2]] <- subclusters
  l[[3]] <- k
  l[[4]] <- length(inds)
  return(l)
}

#' @rdname indexCell
#' @aliases indexCell
setMethod("indexCell", "SingleCellExperiment", indexCell.SingleCellExperiment)
