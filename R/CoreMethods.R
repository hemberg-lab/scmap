#' Select the most important features (genes/transcripts) for mapping
#'
#' The features are selected by find positive residuals from the linear model of
#' dropouts versus expression dependency (both in a log scale). This method
#' is similar to the one used in M3Drop package, where the right side feautres
#' are selected as the differentially expressed.
#' 
#' Please note that \code{object@featureData@data$feature_symbol} column must be 
#' present in the input object. This column defines feature names used during mapping
#' Feature symbols in the reference dataset must correpond to the feature symbols
#' in the mapping dataset, otherwise the mapping will not work!
#'
#' @param object an object of \code{\link[scater]{SCESet}} class
#' @param n_features number of the features to be selected. Note that the number
#' of the output selected features can be lower than the requested 'n_features'.
#' @param pct_dropout_min lower threshold for dropout rate. Features with the
#' dropout rate lower than \code{pct_dropout_min} will not be selected.
#' @param pct_dropout_max upper threshold for dropout rate. Features with the
#' dropout rate lower than \code{pct_dropout_min} will not be selected.
#' @param suppress_plot defines whether to plot the linear model fit to the data
#' together with the highlighted residuals.
#'
#' @return an object of \code{\link[scater]{SCESet}} class with a new column in 
#' \code{featureData} slot which is called \code{scmap_features}. It can be accessible
#' either by \code{fData(object)$scmap_features} or by \code{object@featureData@data$scmap_features}
#'
#' @name getFeatures
#'
#' @importFrom scater fData<-
#' @importFrom methods new
#' @importFrom stats lm
#'
#' @export
getFeatures.SCESet <- function(object, n_features, pct_dropout_min, pct_dropout_max, 
    suppress_plot) {
    if (is.null(object@featureData@data$feature_symbol)) {
        message("There is no feature_symbol column in the featureData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
        return(object)
    }
    
    f_data <- object@featureData@data
    
    # do not consider ERCC spike-ins and genes with 0 dropout rate
    dropouts_filter <- which(f_data$pct_dropout != 0 & !grepl("ERCC-", object@featureData@data$feature_symbol))
    dropouts <- log2(f_data$pct_dropout[dropouts_filter])
    expression <- f_data$mean_exprs[dropouts_filter]
    
    fit <- lm(dropouts ~ expression)
    gene_inds <- as.numeric(names(head(sort(fit$residuals[fit$residuals > 0 & f_data$pct_dropout[dropouts_filter] > 
        pct_dropout_min & f_data$pct_dropout[dropouts_filter] < pct_dropout_max], decreasing = T), 
        n_features)))
    
    f_data$scmap_features <- FALSE
    f_data$scmap_features[dropouts_filter[gene_inds]] <- TRUE
    fData(object) <- new("AnnotatedDataFrame", data = f_data)
    
    if (!suppress_plot) {
        plot(expression, dropouts, xlab = "log2(Expression)", ylab = "log2(% of dropouts)", pch = 16, cex = 0.2)
        points(expression[gene_inds], dropouts[gene_inds], col = "red", pch = 16, cex = 0.7)
        abline(fit, col = "red")
    }
    
    return(object)
}

#' @rdname getFeatures
#' @aliases getFeatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("getFeatures", signature(object = "SCESet"), function(object, n_features, pct_dropout_min, 
    pct_dropout_max, suppress_plot) {
    getFeatures.SCESet(object, n_features, pct_dropout_min, pct_dropout_max, suppress_plot)
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
#' @importFrom dplyr group_by summarise %>%
#' @importFrom reshape2 melt dcast
#' @importFrom scater pData<- get_exprs
#' @importFrom proxy simil
#' @importFrom graphics abline hist plot points
#' @importFrom stats median
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
    
    if (length(which(object_ref@featureData@data$scmap_features)) < 2) {
        warning("There are no common features between the Reference and Projection datasets! Either check that both datasets are from the same organism or increase the number of selected features (>100).")
        return(object_map)
    }
    
    if (!class_col %in% colnames(object_ref@phenoData@data)) {
        warning(paste0("Please define a correct class column of the reference scater object pData slot using the `class_col` parameter!"))
        return(object_map)
    }
    
    if (method == "scmap") {
        gene <- cell_class <- exprs <- NULL
        
        if (is.null(class_ref)) {
            # calculate median feature expression in every cell class of object_ref
            classes <- object_ref@phenoData@data[[class_col]]
            if (is.null(classes)) {
                warning(paste0("Please define a correct class column of the reference scater object using the `class_col` parameter!"))
                return(object_map)
            }
            class_ref <- get_exprs(object_ref[object_ref@featureData@data$scmap_features, ], "exprs")
            rownames(class_ref) <- object_ref[object_ref@featureData@data$scmap_features, ]@featureData@data$feature_symbol
            colnames(class_ref) <- classes
            class_ref <- reshape2::melt(class_ref)
            colnames(class_ref) <- c("gene", "cell_class", "exprs")
            class_ref <- class_ref %>% group_by(gene, cell_class) %>% summarise(med_exprs = median(exprs))
            class_ref <- reshape2::dcast(class_ref, gene ~ cell_class, value.var = "med_exprs")
            rownames(class_ref) <- class_ref$gene
            class_ref <- class_ref[, 2:ncol(class_ref)]
        }
        
        dat <- get_exprs(object_map, "exprs")
        rownames(dat) <- object_map@featureData@data$feature_symbol
        dat <- dat[rownames(dat) %in% rownames(class_ref), ]
        dat <- dat[order(rownames(dat)), ]
        
        class_ref <- class_ref[rownames(class_ref) %in% rownames(dat), ]
        class_ref <- class_ref[order(rownames(class_ref)), ]
        class_ref <- class_ref[, colSums(class_ref) > 0]
        
        if (ncol(class_ref) == 0) {
            warning(paste0("Median expression in the selected features is 0 in every cell, please redefine your features!"))
            return(object_map)
        }
        
        original_classes <- colnames(class_ref)
        
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
        
        tmp <- cbind(
            original_classes[max_inds1],
            original_classes[max_inds2],
            original_classes[max_inds3]
        )
        
        class_assigned <- tmp[,1]
        class_assigned[apply(tmp, 1, function(x) (length(unique(x)) > 2))] <- "unassigned"
        
        maxs <- cbind(
            maxs1,
            maxs2,
            maxs3
        )
        
        maxs <- apply(maxs, 1, max)
        
        class_assigned[maxs < threshold] <- "unassigned"
        class_assigned[is.na(class_assigned)] <- "unassigned"
        p_data <- object_map@phenoData@data
        p_data$scmap_labs <- class_assigned
        p_data$scmap_siml <- maxs
        pData(object_map) <- new("AnnotatedDataFrame", data = p_data)
    }
    
    if (method == "svm") {
        
        classes <- object_ref@phenoData@data[[class_col]]
        if (is.null(classes)) {
            warning(paste0("Please define a correct class column of the reference scater object using the `class_col` parameter!"))
            return(object_map)
        }
        class_ref <- get_exprs(object_ref[object_ref@featureData@data$scmap_features, ], "exprs")
        rownames(class_ref) <- object_ref[object_ref@featureData@data$scmap_features, ]@featureData@data$feature_symbol
        colnames(class_ref) <- classes
        
        dat <- get_exprs(object_map, "exprs")
        rownames(dat) <- object_map@featureData@data$feature_symbol
        dat <- dat[rownames(dat) %in% rownames(class_ref), ]
        dat <- dat[order(rownames(dat)), ]
        
        class_ref <- class_ref[rownames(class_ref) %in% rownames(dat), ]
        class_ref <- class_ref[order(rownames(class_ref)), ]
        
        res <- support_vector_machines(class_ref, dat)
        probs <- attr(res, "probabilities")
        max_inds <- max.col(probs)
        maxs <- unlist(apply(probs, 1, max))
        
        labs <- rep("unassigned", length(max_inds))
        labs[maxs > threshold] <- colnames(probs)[max_inds][maxs > threshold]

        p_data <- object_map@phenoData@data
        p_data$scmap_labs <- labs
        p_data$scmap_siml <- maxs
        pData(object_map) <- new("AnnotatedDataFrame", data = p_data)
    }
    
    if (method == "rf") {
        
        if (is.null(object_ref@phenoData@data[[class_col]])) {
            warning(paste0("Please define a correct class column of the reference scater object using the `class_col` parameter!"))
            return(object_map)
        }
        class_ref <- get_exprs(object_ref[object_ref@featureData@data$scmap_features, ], "exprs")
        rownames(class_ref) <- object_ref[object_ref@featureData@data$scmap_features, ]@featureData@data$feature_symbol
        colnames(class_ref) <- object_ref@phenoData@data[[class_col]]
        
        dat <- get_exprs(object_map, "exprs")
        rownames(dat) <- object_map@featureData@data$feature_symbol
        dat <- dat[rownames(dat) %in% rownames(class_ref), ]
        dat <- dat[order(rownames(dat)), ]
        
        class_ref <- class_ref[rownames(class_ref) %in% rownames(dat), ]
        class_ref <- class_ref[order(rownames(class_ref)), ]
        
        res <- random_forest(class_ref, dat)
        max_inds <- max.col(res)
        maxs <- unlist(apply(res, 1, max))
        
        labs <- rep("unassigned", length(max_inds))
        labs[maxs > threshold] <- colnames(res)[max_inds][maxs > threshold]
        
        p_data <- object_map@phenoData@data
        p_data$scmap_labs <- labs
        p_data$scmap_siml <- maxs
        pData(object_map) <- new("AnnotatedDataFrame", data = p_data)
    }
    
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

