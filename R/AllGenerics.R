#' @export
setGeneric("getFeatures", function(object, n_features = 100, pct_dropout_min = 20, pct_dropout_max = 80, 
    suppress_plot = TRUE) {
    standardGeneric("getFeatures")
})

#' @export
setGeneric("setFeatures", function(object, features = NULL) {
    standardGeneric("setFeatures")
})

#' @export
setGeneric("mapData", function(object_map = NULL, object_ref = NULL, class_col = "cell_type1", 
    class_ref = NULL, threshold = 0.7, suppress_plot = TRUE) {
    standardGeneric("mapData")
})
