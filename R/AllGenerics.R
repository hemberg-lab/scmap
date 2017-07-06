#' @export
#' 
#' @examples
#' library(scater)
#' pd <- AnnotatedDataFrame(ann)
#' sceset <- newSCESet(fpkmData = yan, phenoData = pd, logExprsOffset = 1)
#' sceset <- calculateQCMetrics(sceset)
#' # use gene names as feature symbols
#' fData(sceset)$feature_symbol <- featureNames(sceset)
#' # remove features with duplicated names
#' sceset <- sceset[!duplicated(fData(sceset)$feature_symbol), ]
#' sceset <- getFeatures(sceset)
#' 
setGeneric("getFeatures", signature = "object", function(object, n_features = 500, suppress_plot = TRUE) {
    standardGeneric("getFeatures")
})

#' @export
#' 
#' @examples
#' library(scater)
#' pd <- AnnotatedDataFrame(ann)
#' sceset <- newSCESet(fpkmData = yan, phenoData = pd, logExprsOffset = 1)
#' sceset <- calculateQCMetrics(sceset)
#' # use gene names as feature symbols
#' fData(sceset)$feature_symbol <- featureNames(sceset)
#' # remove features with duplicated names
#' sceset <- sceset[!duplicated(fData(sceset)$feature_symbol), ]
#' sceset <- setFeatures(sceset, c('MMP2', 'ZHX3'))
#' 
setGeneric("setFeatures", signature = "object", function(object, features = NULL) {
    standardGeneric("setFeatures")
})

#' @export
#' 
#' @examples
#' library(scater)
#' pd <- AnnotatedDataFrame(ann)
#' sceset <- newSCESet(fpkmData = yan, phenoData = pd, logExprsOffset = 1)
#' sceset <- calculateQCMetrics(sceset)
#' # use gene names as feature symbols
#' fData(sceset)$feature_symbol <- featureNames(sceset)
#' # remove features with duplicated names
#' sceset <- sceset[!duplicated(fData(sceset)$feature_symbol), ]
#' sceset <- getFeatures(sceset)
#' sceset <- projectData(projection = sceset, reference = sceset)
#' 
setGeneric("projectData", signature = "projection", function(projection = NULL, reference = NULL, cell_type_column = "cell_type1", 
    method = "scmap", threshold = 0.7) {
    standardGeneric("projectData")
})

#' @export
#' 
#' @examples
#' library(scater)
#' pd <- AnnotatedDataFrame(ann)
#' sceset <- newSCESet(fpkmData = yan, phenoData = pd, logExprsOffset = 1)
#' sceset <- calculateQCMetrics(sceset)
#' # use gene names as feature symbols
#' fData(sceset)$feature_symbol <- featureNames(sceset)
#' # remove features with duplicated names
#' sceset <- sceset[!duplicated(fData(sceset)$feature_symbol), ]
#' sceset <- getFeatures(sceset)
#' reference <- createReference(sceset[fData(sceset)$scmap_features, ])
#' 
setGeneric("createReference", signature = "reference", function(reference = NULL, cell_type_column = "cell_type1") {
    standardGeneric("createReference")
})

