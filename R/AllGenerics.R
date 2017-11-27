#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' sce <- selectFeatures(sce)
#' 
setGeneric("selectFeatures", signature = "object", function(object, n_features = 500, suppress_plot = TRUE) {
    standardGeneric("selectFeatures")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' sce <- setFeatures(sce, c('MMP2', 'ZHX3'))
#' 
setGeneric("setFeatures", signature = "object", function(object, features = NULL) {
    standardGeneric("setFeatures")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' sce <- getFeatures(sce)
#' reference <- indexCluster(sce[rowData(sce)$scmap_features, ])
#' 
setGeneric("indexCluster", signature = "object", function(object = NULL, cluster_col = "cell_type1") {
    standardGeneric("indexCluster")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' sce <- getFeatures(sce)
#' sce <- indexCell(sce)
#' 
setGeneric("indexCell", signature = "object", function(object = NULL, M = NULL, k = NULL) {
  standardGeneric("indexCell")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' sce <- getFeatures(sce)
#' sce <- indexCluster(sce)
#' sce <- scmapCluster(sce, sce)
#' 
setGeneric("scmapCluster", signature = "projection", function(projection = NULL, index_list = NULL, 
                                                              threshold = 0.7) {
  standardGeneric("scmapCluster")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' sce <- getFeatures(sce)
#' sce <- indexCell(sce)
#' scmapCell_results <- scmapCell(sce, sce)
#' 
setGeneric("scmapCell", signature = "projection", function(projection = NULL, index_list = NULL, w = 10) {
  standardGeneric("scmapCell")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' sce <- getFeatures(sce)
#' sce <- indexCell(sce)
#' scmapCell_results <- scmapCell(sce, sce)
#' sce <- scmapCell2Cluster(sce, scmapCell_results, cluster_list = list(colData(sce)$cell_type1))
#' 
setGeneric("scmapCell2Cluster", signature = "projection", function(projection = NULL, scmapCell_results = NULL, cluster_list = NULL, w = 3, threshold = 0.5) {
  standardGeneric("scmapCell2Cluster")
})

