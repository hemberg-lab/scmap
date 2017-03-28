
#' Define a set of genes used for creating cell binary signatures
#'
#' The important genes (hyperplanes) are defined as differentially expressed
#' genes identified using M3Drop (Michaelis-Menten Modelling of Dropouts for scRNASeq,
#' http://bioconductor.org/packages/M3Drop). A 'data.frame' with the genes
#' and corresponding p and q values is written to the 'hyperplanes' item of the
#' 'sc3' slot of the input object.
#'
#' @param object an object of 'SCESet' class
#' @param n_features number of the genes to be returned
#' @param pct_dropout_min lower threshold for dropout rate
#' @param pct_dropout_max upper threshold for dropout rate
#' @param suppress_plot defines whether to plot the fit
#'
#' @return an object of 'SCESet' class
#'
#' @importFrom scater fData<-
#' @importFrom methods new
#' @importFrom stats lm
#'
#' @export
getFeatures.SCESet <- function(object, n_features = 100, pct_dropout_min = 20, pct_dropout_max = 80, suppress_plot = T) {
    f_data <- object@featureData@data

    # do not consider ERCC spike-ins and genes with 0 dropout rate
    dropouts_filter <- which(f_data$pct_dropout != 0 & !grepl("ERCC-", featureNames(object)))
    dropouts <- log10(f_data$pct_dropout[dropouts_filter])
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
        plot(expression, dropouts, xlab = "log2(Expression)", ylab = "log10(% of dropouts)")
        points(expression[gene_inds], dropouts[gene_inds], col = "red")
        abline(fit, col = "red")
    }

    return(object)
}

#' @rdname getFeatures.SCESet
#' @aliases getFeatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("getFeatures", signature(object = "SCESet"), function(object, n_features = 100, pct_dropout_min = 20, pct_dropout_max = 80, suppress_plot = T) {
    getFeatures.SCESet(object, n_features, pct_dropout_min, pct_dropout_max, suppress_plot)
})

#' Title
#'
#' @param object an object of 'SCESet' class
#' @param features a vector of feature names
#'
#' @importFrom scater fData<-
#' @importFrom methods new
#' @importFrom Biobase featureNames
#'
#' @export
setFeatures.SCESet <- function(object, features) {
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

#' @rdname setFeatures.SCESet
#' @aliases setFeatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("setFeatures", signature(object = "SCESet"), function(object, features) {
    setFeatures.SCESet(object, features)
})


