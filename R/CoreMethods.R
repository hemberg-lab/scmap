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
#' log(expression) versus log(dropout) distribution for all genes.
#' Selected features are highlighted with the red colour.
#'
#' @return an object of \code{\link[scater]{SCESet}} class with a new column in 
#' \code{featureData} slot which is called \code{scmap_features}. It can be accessed
#' by using \code{fData(object)$scmap_features}.
#'
#' @name getFeatures
#'
#' @importFrom Biobase fData fData<- AnnotatedDataFrame
getFeatures.SCESet <- function(object, n_features, suppress_plot) {
    if (is.null(fData(object)$feature_symbol)) {
        stop("There is no feature_symbol column in the featureData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
        return(object)
    }
    
    f_data <- fData(object)
    tmp <- linearModel(f_data, n_features)
    f_data$scmap_features <- tmp$scmap_features
    f_data$scmap_scores <- tmp$scmap_scores
    fData(object) <- AnnotatedDataFrame(data = f_data)
    
    if (!suppress_plot) {
        p <- ggplot_features(tmp$for_plotting, tmp$fit)
        print(p)
    }
    
    return(object)
}

#' @rdname getFeatures
#' @aliases getFeatures
#' @importClassesFrom scater SCESet
setMethod("getFeatures", "SCESet", getFeatures.SCESet)

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
#' @importFrom Biobase fData fData<- AnnotatedDataFrame
setFeatures.SCESet <- function(object, features) {
    if (is.null(features)) {
        stop("Please provide a list of feature names using 'features' argument!")
        return(object)
    }
    if (is.null(fData(object)$feature_symbol)) {
        stop("There is no feature_symbol column in the featureData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
        return(object)
    }
    
    inds <- match(features, fData(object)$feature_symbol)
    
    if (!all(!is.na(inds))) {
        warning("Features ", paste(features[which(is.na(inds))], collapse = ", "), " are not present in the 'SCESet' object and therefore were not set.")
    }
    f_data <- fData(object)
    f_data$scmap_features <- FALSE
    f_data$scmap_features[inds[!is.na(inds)]] <- TRUE
    fData(object) <- AnnotatedDataFrame(f_data)
    
    return(object)
}

#' @rdname setFeatures
#' @aliases setFeatures
#' @importClassesFrom scater SCESet
setMethod("setFeatures", "SCESet", setFeatures.SCESet)

#' scmap main function
#' 
#' Projection of one dataset to another
#' 
#' @param projection SCESet to project
#' @param reference reference SCESet set
#' @param cell_type_column column name in the pData slot of the reference SCESet containing the cell classification information
#' @param method which method to use
#' @param threshold threshold on similarity (or probability for SVM and RF)
#' 
#' @return Projection SCESet object with labels calculated by `scmap` stored in 
#' the `scmap_labels` column of the `phenoData` slot.
#' 
#' @name projectData
#' 
#' @importFrom Biobase fData pData pData<- AnnotatedDataFrame
#' @importFrom scater get_exprs
#' @importFrom proxy simil
#' @importFrom stats cor
#' @importFrom matrixStats colMaxs rowMaxs
#' @importFrom methods is
projectData.SCESet <- function(projection, reference, cell_type_column, method, threshold) {
    if (is.null(projection)) {
        stop("Please define a scater object to map using the `projection` parameter!")
        return(projection)
    }
    if (is.null(reference)) {
        stop("Please define either a reference scater object or precomputed scmap reference using the `reference` parameter!")
        return(projection)
    } else {
        if ("SCESet" %in% is(reference)) {
            if (is.null(fData(reference)$scmap_features)) {
                stop("There are no features selected in the Reference dataset! Please run `getFeatures()` first!")
                return(projection)
            }
            if (!cell_type_column %in% colnames(pData(reference))) {
                stop("Please define a correct class column of the reference scater object pData slot using the `cell_type_column` parameter!")
                return(projection)
            }
        } else {
            
        }
    }
    if (!("SCESet" %in% is(reference)) & (method == "svm" | method == "rf")) {
        stop("SVM/RF do not work with the precomputed reference provided for scmap. Please provide the full reference dataset.")
        return(projection)
    }
    
    projection_local <- projection
    
    # find and select only common features, then subset both datasets
    if ("SCESet" %in% is(reference)) {
        projection_local <- setFeatures(projection_local, fData(reference)$feature_symbol[fData(reference)$scmap_features])
        reference <- setFeatures(reference, fData(projection_local)$feature_symbol[fData(projection_local)$scmap_features])
        projection_local <- projection_local[fData(projection_local)$scmap_features, ]
        reference <- reference[fData(reference)$scmap_features, ]
    } else {
        projection_local <- setFeatures(projection_local, rownames(reference))
        reference <- reference[rownames(reference) %in% fData(projection_local)$feature_symbol[fData(projection_local)$scmap_features], 
            , drop = FALSE]
        projection_local <- projection_local[fData(projection_local)$scmap_features, ]
    }
    
    if (is.null(reference)) {
        if (nrow(reference) < 10) {
            warning("There are less than ten features in common between the Reference and Projection datasets. Most probably they come from different organisms!")
            return(projection)
        }
    } else {
        if (nrow(reference) < 10) {
            warning("There are less than ten features in common between the Reference and Projection datasets. Most probably they come from different organisms!")
            return(projection)
        }
    }
    
    # get expression values of the projection dataset
    dat <- get_exprs(projection_local, "exprs")
    rownames(dat) <- fData(projection_local)$feature_symbol
    
    if (method == "scmap") {
        # create the reference
        if ("SCESet" %in% is(reference)) {
            if (is.null(pData(reference)[[cell_type_column]])) {
                stop("Please define a correct class column of the reference scater object using the `cell_type_column` parameter!")
                return(projection)
            }
            reference_local <- createReference(reference, cell_type_column)
        } else {
            reference_local <- reference
        }
        
        # prepare the datasets for projection
        tmp <- prepareData(reference_local, dat)
        reference_local <- tmp$reference
        dat <- tmp$dat
        
        if (ncol(reference_local) == 0) {
            stop("Median expression in the selected features is 0 in every cell, please redefine your features!")
            return(projection)
        }
        
        # run scmap
        tmp <- t(reference_local)
        res <- proxy::simil(tmp, t(dat), method = "cosine")
        res <- matrix(res, ncol = nrow(tmp), byrow = TRUE)
        max_inds1 <- max.col(res)
        maxs1 <- rowMaxs(res)
        
        res <- cor(reference_local, dat, method = "pearson")
        max_inds2 <- max.col(t(res))
        maxs2 <- colMaxs(res)
        
        res <- cor(reference_local, dat, method = "spearman")
        max_inds3 <- max.col(t(res))
        maxs3 <- colMaxs(res)
        
        cons <- cbind(colnames(reference_local)[max_inds1], colnames(reference_local)[max_inds2], colnames(reference_local)[max_inds3])
        
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
            maxs_tmp <- cbind(maximums[unique_labs == 2 & non_na_inds, , drop = FALSE][cbind(seq_along(inds[, 1]), inds[, 
                1])], maximums[unique_labs == 2 & non_na_inds, , drop = FALSE][cbind(seq_along(inds[, 1]), inds[, 2])])
            maxs_tmp <- rowMaxs(maxs_tmp)
            maxs[unique_labs == 2 & non_na_inds] <- maxs_tmp
        }
        
        ## check the similarity threshold
        labs[!is.na(maxs) & maxs < threshold] <- "unassigned"
    }
    
    if (method == "svm") {
        # create the reference
        reference_local <- get_exprs(reference, "exprs")
        rownames(reference_local) <- fData(reference)$feature_symbol
        colnames(reference_local) <- pData(reference)[[cell_type_column]]
        
        # prepare the datasets for projection
        tmp <- prepareData(reference_local, dat)
        reference_local <- tmp$reference
        dat <- tmp$dat
        
        # run support vector machine
        res <- support_vector_machines(reference_local, dat)
        probs <- attr(res, "probabilities")
        max_inds <- max.col(probs)
        maxs <- rowMaxs(probs)
        
        # create labels
        labs <- rep("unassigned", length(max_inds))
        labs[maxs > threshold] <- colnames(probs)[max_inds][maxs > threshold]
    }
    
    if (method == "rf") {
        # create the reference
        reference_local <- get_exprs(reference, "exprs")
        rownames(reference_local) <- fData(reference)$feature_symbol
        colnames(reference_local) <- pData(reference)[[cell_type_column]]
        
        # prepare the datasets for projection
        tmp <- prepareData(reference_local, dat)
        reference_local <- tmp$reference
        dat <- tmp$dat
        
        # run random forest
        res <- random_forest(reference_local, dat)
        max_inds <- max.col(res)
        maxs <- rowMaxs(res)
        
        # create labels
        labs <- rep("unassigned", length(max_inds))
        labs[maxs > threshold] <- colnames(res)[max_inds][maxs > threshold]
    }
    
    p_data <- pData(projection)
    p_data$scmap_labs <- labs
    p_data$scmap_siml <- maxs
    pData(projection) <- AnnotatedDataFrame(p_data)
    
    return(projection)
}

#' @rdname projectData
#' @aliases projectData
#' @importClassesFrom scater SCESet
setMethod("projectData", "SCESet", projectData.SCESet)

#' Create a precomputed Reference
#' 
#' Calculates centroids of each cell type and merge them into a single table.
#'
#' @param reference reference SCESet set
#' @param cell_type_column column name in the pData slot of the reference SCESet 
#' containing the cell classification information
#' 
#' @name createReference
#'
#' @return a `data.frame` containing calculated centroids of the cell types of
#' the Reference dataset
#'
#' @importFrom Biobase fData pData
#' @importFrom scater get_exprs
#' @importFrom dplyr group_by summarise %>%
#' @importFrom reshape2 melt dcast
#' @importFrom stats median
createReference.SCESet <- function(reference, cell_type_column) {
    if (is.null(reference)) {
        stop("Please define a reference scater object using the `reference` parameter!")
    }
    gene <- cell_class <- exprs <- NULL
    reference_local <- get_exprs(reference, "exprs")
    rownames(reference_local) <- fData(reference)$feature_symbol
    colnames(reference_local) <- pData(reference)[[cell_type_column]]
    
    # calculate median feature expression in every cell class of reference
    reference_local <- reshape2::melt(reference_local)
    colnames(reference_local) <- c("gene", "cell_class", "exprs")
    reference_local <- reference_local %>% group_by(gene, cell_class) %>% summarise(med_exprs = median(exprs))
    reference_local <- reshape2::dcast(reference_local, gene ~ cell_class, value.var = "med_exprs")
    rownames(reference_local) <- reference_local$gene
    reference_local <- reference_local[, 2:ncol(reference_local), drop = FALSE]
    return(reference_local)
}

#' @rdname createReference
#' @aliases createReference
#' @importClassesFrom scater SCESet
setMethod("createReference", "SCESet", createReference.SCESet)
