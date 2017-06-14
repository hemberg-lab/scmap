#' @export
setGeneric("getFeatures", function(object, n_features = 100, pct_dropout_min = 0, pct_dropout_max = 100, 
    suppress_plot = TRUE) {
    standardGeneric("getFeatures")
})

#' @export
setGeneric("setFeatures", function(object, features = NULL) {
    standardGeneric("setFeatures")
})

#' @export
setGeneric("projectData", function(object_map = NULL, object_ref = NULL, class_col = "cell_type1", 
    class_ref = NULL, method = "scmap", threshold = 0.7) {
    standardGeneric("projectData")
})
