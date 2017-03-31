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
#' @importFrom scater fData<-
#' @importFrom methods new
#' @importFrom stats lm
#'
#' @export
getFeatures.scater <- function(object, n_features = 100, pct_dropout_min = 20, 
                               pct_dropout_max = 80, suppress_plot = TRUE) {
    if(is.null(object@featureData@data$feature_symbol)) {
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
    gene_inds <- as.numeric(
        names(
            head(
                sort(
                    fit$residuals[
                        fit$residuals > 0 &
                        f_data$pct_dropout[dropouts_filter] > pct_dropout_min &
                        f_data$pct_dropout[dropouts_filter] < pct_dropout_max
                    ],
                    decreasing = T
                ),
                n_features
            )
        )
    )

    f_data$scmap_features <- FALSE
    f_data$scmap_features[dropouts_filter[gene_inds]] <- TRUE
    fData(object) <- new("AnnotatedDataFrame", data = f_data)

    if(!suppress_plot) {
        plot(expression, dropouts, xlab = "log2(Expression)", ylab = "log2(% of dropouts)")
        points(expression[gene_inds], dropouts[gene_inds], col = "red")
        abline(fit, col = "red")
    }

    return(object)
}

#' @rdname getFeatures.scater
#' @aliases getFeatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("getFeatures", signature(object = "SCESet"), function(object, n_features = 100, pct_dropout_min = 20, pct_dropout_max = 80, suppress_plot = T) {
    getFeatures.scater(object, n_features, pct_dropout_min, pct_dropout_max, suppress_plot)
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
#' @param features a vector of feature names
#' 
#' @return an object of \code{\link[scater]{SCESet}} class with a new column in 
#' \code{featureData} slot which is called \code{scmap_features}. It can be accessible
#' either by \code{fData(object)$scmap_features} or by \code{object@featureData@data$scmap_features}
#'
#' @importFrom scater fData<-
#' @importFrom methods new
#'
#' @export
setFeatures.scater <- function(object, features) {
    if(is.null(object@featureData@data$feature_symbol)) {
        message("There is no feature_symbol column in the featureData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
        return(object)
    }
    
    inds <- match(features, object@featureData@data$feature_symbol)

    if(!all(!is.na(inds))) {
        message(paste0("Features ", paste(features[which(is.na(inds))], collapse = ", "), " are not present in the 'SCESet' object and therefore were not set."))
    }
    f_data <- object@featureData@data
    f_data$scmap_features <- FALSE
    f_data$scmap_features[inds[!is.na(inds)]] <- TRUE
    fData(object) <- new("AnnotatedDataFrame", data = f_data)

    return(object)
}

#' @rdname setFeatures.scater
#' @aliases setFeatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("setFeatures", signature(object = "SCESet"), function(object, features) {
    setFeatures.scater(object, features)
})


