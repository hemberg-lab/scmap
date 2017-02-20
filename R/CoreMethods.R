
#' Define a set of genes used for creating cell binary signatures
#'
#' The important genes (hyperplanes) are defined as differentially expressed
#' genes identified using M3Drop (Michaelis-Menten Modelling of Dropouts for scRNASeq,
#' http://bioconductor.org/packages/M3Drop). A 'data.frame' with the genes
#' and corresponding p and q values is written to the 'hyperplanes' item of the
#' 'sc3' slot of the input object.
#'
#' @param object an object of 'SCESet' class
#' @param n_genes number of the genes to be returned
#'
#' @return an object of 'SCESet' class
#'
#' @importFrom scater fData<-
#' @importFrom SC3 get_processed_dataset
#'
#' @export
getFeatures.SCESet <- function(object, n_genes = 100, pct_dropout_min = 20, pct_dropout_max = 80, suppress_plot = T) {
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
                n_genes
            )
        )
    )

    f_data$fsc3_features <- FALSE
    f_data$fsc3_features[dropouts_filter[gene_inds]] <- TRUE
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
setMethod("getFeatures", signature(object = "SCESet"), function(object, n_genes = 100, pct_dropout_min = 20, pct_dropout_max = 80, suppress_plot = T) {
    getFeatures.SCESet(object, n_genes, pct_dropout_min, pct_dropout_max, suppress_plot)
})

#' Title
#'
#'
#' @param features a vector of feature names
#'
#' @importFrom scater fData<-
#' @importFrom Biobase featureNames
#'
#' @export
setFeatures.SCESet <- function(object, features) {
    inds <- match(features, featureNames(object))

    if(!all(!is.na(inds))) {
        message(paste0("Features ", paste(features[which(is.na(inds))], collapse = ", "), " are not present in the 'SCESet' object and therefore were not set."))
    }
    f_data <- object@featureData@data
    f_data$fsc3_features <- FALSE
    f_data$fsc3_features[inds[!is.na(inds)]] <- TRUE
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

#' Create a binary signature for each cell
#'
#'
#' @param object an object of 'SCESet' class
#' @param exprs_values
#'
#' @importFrom SC3 get_processed_dataset
#' @importFrom scater pData<-
#'
#' @return an object of 'SCESet' class
#'
#' @export
getSignatures.SCESet <- function(object, exprs_values = "exprs") {
    if (is.null(object@featureData@data$fsc3_features)) {
        warning(paste0("Please run fsc3_get_features() first!"))
        return(object)
    }

    data <- object@assayData[[exprs_values]]
    data <- data[object@featureData@data$fsc3_features, ]
    data <- data[order(rownames(data)), ]

    p_data <- object@phenoData@data
    p_data$fsc3_signatures <- signature_mapper(data)
    pData(object) <- new("AnnotatedDataFrame", data = p_data)

    return(object)
}

#' @rdname getSignatures.SCESet
#' @aliases getSignatures
#' @importClassesFrom scater SCESet
#' @export
setMethod("getSignatures", signature(object = "SCESet"), function(object) {
    getSignatures.SCESet(object)
})
