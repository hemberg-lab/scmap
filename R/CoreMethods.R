#' Select the most informative features (genes/transcripts) for projection
#'
#' This is a modification of the M3Drop method. Instead of fitting a 
#' Michaelis-Menten model to the log expression-dropout relation, we fit a 
#' linear model. Namely, the linear model is build on the log(expression) versus 
#' log(dropout) distribution. After fitting a linear model important features are
#'  selected as the top N residuals of the linear model.
#' 
#' Please note that \code{object@featureData@data$feature_symbol} column must be 
#' present in the input object and should not contain any duplicated feature names. 
#' This column defines feature names used during projection. Feature symbols 
#' in the reference dataset must correpond to the feature symbols
#' in the projection dataset, otherwise the mapping will not work!
#'
#' @param object an object of \code{\link[scater]{SCESet}} class
#' @param n_features number of the features to be selected
#' @param suppress_plot boolean parameter, which defines whether to plot 
#' log(expression) versus log(dropout) distribution for all genes
#'
#' @return an object of \code{\link[scater]{SCESet}} class with a new column in 
#' \code{featureData} slot which is called \code{scmap_features}. It can be accessed
#' by using \code{fData(object)$scmap_features}.
#'
#' @name getFeatures
#'
#' @importFrom scater fData<-
#' @importFrom methods new
#'
#' @export
getFeatures.SCESet <- function(object, n_features, suppress_plot) {
    if (is.null(object@featureData@data$feature_symbol)) {
        message("There is no feature_symbol column in the featureData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
        return(object)
    }
    
    f_data <- object@featureData@data
    tmp <- linearModel(f_data, n_features)
    f_data$scmap_features <- tmp$scmap_features
    fData(object) <- new("AnnotatedDataFrame", data = f_data)
    
    if (!suppress_plot) {
        p <- ggplot_features(tmp$for_plotting, tmp$fit)
        print(p)
    }
    
    return(object)
}

#' @rdname getFeatures
#' @aliases getFeatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("getFeatures", signature(object = "SCESet"), function(object, n_features, suppress_plot) {
    getFeatures.SCESet(object, n_features, suppress_plot)
})

#' Set the most important features (genes/transcripts) for mapping
#' 
#' This method set the features to be used for mapping.
#' 
#' Please note that \code{object@featureData@data$feature_symbol} column must be 
#' present in the input object. This column defines feature names used during mapping
#' Feature symbols in the reference dataset must correpond to the feature symbols
#' in the mapping dataset, otherwise the mapping will not work!
#'
#' @param object an object of \code{\link[scater]{SCESet}} class
#' @param features a character vector of feature names
#' 
#' @return an object of \code{\link[scater]{SCESet}} class with a new column in 
#' \code{featureData} slot which is called \code{scmap_features}. It can be accessible
#' either by \code{fData(object)$scmap_features} or by \code{object@featureData@data$scmap_features}
#'
#' @name setFeatures
#'
#' @importFrom scater fData<-
#' @importFrom methods new
#'
#' @export
setFeatures.SCESet <- function(object, features) {
    if (is.null(features)) {
        message("Please provide a list of feature names using 'features' argument!")
        return(object)
    }
    if (is.null(object@featureData@data$feature_symbol)) {
        message("There is no feature_symbol column in the featureData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
        return(object)
    }
    
    inds <- match(features, object@featureData@data$feature_symbol)
    
    if (!all(!is.na(inds))) {
        message(paste0("Features ", paste(features[which(is.na(inds))], collapse = ", "), 
            " are not present in the 'SCESet' object and therefore were not set."))
    }
    f_data <- object@featureData@data
    f_data$scmap_features <- FALSE
    f_data$scmap_features[inds[!is.na(inds)]] <- TRUE
    fData(object) <- new("AnnotatedDataFrame", data = f_data)
    
    return(object)
}

#' @rdname setFeatures
#' @aliases setFeatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("setFeatures", signature(object = "SCESet"), function(object, features) {
    setFeatures.SCESet(object, features)
})

#' scmap main function
#' 
#' Mapping of one dataset to another
#' 
#' @param object_map SCESet to map
#' @param object_ref reference SCESet set
#' @param class_col column name in the pData slot of the reference SCESet containing the cell classification information
#' @param class_ref reference cell buckets
#' @param method which method to use
#' @param threshold threshold on similarity (or probability for SVM and RF)
#' 
#' @name projectData
#' 
#' @importFrom scater pData<- get_exprs
#' @importFrom proxy simil
#' @importFrom graphics abline hist plot points
#' @importFrom utils head
#' @importFrom stats cor
#' @export
projectData.SCESet <- function(object_map, object_ref, class_col, class_ref, method, threshold) {
    if (is.null(object_map)) {
        warning(paste0("Please define a scater object to map using the `object_map` parameter!"))
        return(object_map)
    }
    if (is.null(object_ref)) {
        if (is.null(class_ref)) {
            warning(paste0("Please define either a reference scater object using the `object_ref` parameter or scmap reference buckets using the `class_ref`!"))
            return(object_map)
        }
    }
    
    if (is.null(object_ref@featureData@data$scmap_features)) {
        warning("There are no features selected in the Reference dataset! Please run `getFeatures()` first!")
        return(object_map)
    }
    
    if (length(which(object_ref@featureData@data$scmap_features)) < 2) {
        warning("There are no common features between the Reference and Projection datasets! Either check that both datasets are from the same organism or increase the number of selected features (>100).")
        return(object_map)
    }
    
    if (!class_col %in% colnames(object_ref@phenoData@data)) {
        warning(paste0("Please define a correct class column of the reference scater object pData slot using the `class_col` parameter!"))
        return(object_map)
    }
    
    # get expression values of the projection dataset
    dat <- get_exprs(object_map, "exprs")
    rownames(dat) <- object_map@featureData@data$feature_symbol
    
    if (method == "scmap") {
        # create the reference
        if (is.null(class_ref)) {
            class_ref <- createReference(object_ref, class_col)
        }
        
        # find feature overlap between the datasets
        ovrlp <- overlapData(class_ref, dat)
        class_ref <- ovrlp$class_ref
        dat <- ovrlp$dat
        
        if (ncol(class_ref) == 0) {
            warning(paste0("Median expression in the selected features is 0 in every cell, please redefine your features!"))
            return(object_map)
        }
        
        # run scmap
        tmp <- t(class_ref)
        res <- proxy::simil(tmp, t(dat), method = "cosine")
        res <- matrix(res, ncol = nrow(tmp), byrow = T)
        max_inds1 <- max.col(res)
        maxs1 <- unlist(apply(res, 1, max))
        
        res <- cor(class_ref, dat, method = "pearson")
        max_inds2 <- max.col(t(res))
        maxs2 <- unlist(apply(res, 2, max))
        
        res <- cor(class_ref, dat, method = "spearman")
        max_inds3 <- max.col(t(res))
        maxs3 <- unlist(apply(res, 2, max))
        
        cons <- cbind(
            colnames(class_ref)[max_inds1],
            colnames(class_ref)[max_inds2],
            colnames(class_ref)[max_inds3]
        )
        
        maximums <- cbind(
            maxs1,
            maxs2,
            maxs3
        )
        
        # create labels
        maxs <- rep(NA, nrow(cons))
        labs <- rep("unassigned", nrow(cons))
        unique_labs <- unlist(apply(cons, 1, function(x) {length(unique(x))}))
        
        ## all similarities agree
        labs[unique_labs == 1] <- cons[unique_labs == 1, 1]
        maxs_tmp <- unlist(apply(maximums[unique_labs == 1, ], 1, max))
        maxs[unique_labs == 1] <- maxs_tmp
        
        ## only two similarities agree
        tmp <- cons[unique_labs == 2, ]
        inds <- unlist(apply(tmp, 1, function(x) {which(duplicated(x))}))
        labs[unique_labs == 2] <- tmp[cbind(seq_along(inds), inds)]
        
        ## calculate maximum similarity in case of two agreeing similarities
        inds <- t(apply(apply(tmp, 2, `==`, labs[unique_labs == 2]), 1, which))
        maxs_tmp <- cbind(
            maximums[unique_labs == 2, ][cbind(seq_along(inds[,1]), inds[,1])],
            maximums[unique_labs == 2, ][cbind(seq_along(inds[,1]), inds[,2])]
        )
        maxs_tmp <- apply(maxs_tmp, 1, max)
        maxs[unique_labs == 2] <- maxs_tmp
        
        ## check the similarity threshold
        labs[!is.na(maxs) & maxs < threshold] <- "unassigned"
    }
    
    if (method == "svm") {
        # create the reference
        class_ref <- get_exprs(object_ref[object_ref@featureData@data$scmap_features, ], "exprs")
        rownames(class_ref) <- object_ref[object_ref@featureData@data$scmap_features, ]@featureData@data$feature_symbol
        colnames(class_ref) <- object_ref@phenoData@data[[class_col]]
        
        # find feature overlap between the datasets
        ovrlp <- overlapData(class_ref, dat)
        class_ref <- ovrlp$class_ref
        dat <- ovrlp$dat
        
        # run support vector machine
        res <- support_vector_machines(class_ref, dat)
        probs <- attr(res, "probabilities")
        max_inds <- max.col(probs)
        maxs <- unlist(apply(probs, 1, max))
        
        # create labels
        labs <- rep("unassigned", length(max_inds))
        labs[maxs > threshold] <- colnames(probs)[max_inds][maxs > threshold]
    }
    
    if (method == "rf") {
        # create the reference
        class_ref <- get_exprs(object_ref[object_ref@featureData@data$scmap_features, ], "exprs")
        rownames(class_ref) <- object_ref[object_ref@featureData@data$scmap_features, ]@featureData@data$feature_symbol
        colnames(class_ref) <- object_ref@phenoData@data[[class_col]]
        
        # find feature overlap between the datasets
        ovrlp <- overlapData(class_ref, dat)
        class_ref <- ovrlp$class_ref
        dat <- ovrlp$dat
        
        # run random forest
        res <- random_forest(class_ref, dat)
        max_inds <- max.col(res)
        maxs <- unlist(apply(res, 1, max))
        
        # create labels
        labs <- rep("unassigned", length(max_inds))
        labs[maxs > threshold] <- colnames(res)[max_inds][maxs > threshold]
    }
    
    p_data <- object_map@phenoData@data
    p_data$scmap_labs <- labs
    p_data$scmap_siml <- maxs
    pData(object_map) <- new("AnnotatedDataFrame", data = p_data)
    
    return(object_map)
}

#' @rdname projectData
#' @aliases projectData
#' @importClassesFrom scater SCESet
#' @export
setMethod("projectData", signature(object_map = "SCESet"), function(object_map, object_ref, 
    class_col, class_ref, method, threshold) {
    projectData.SCESet(object_map, object_ref, class_col, class_ref, method, threshold)
})

