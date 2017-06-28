#' @export
setGeneric("getFeatures", function(object, n_features = 500, suppress_plot = TRUE) {
    standardGeneric("getFeatures")
})

#' @export
setGeneric("setFeatures", function(object, features = NULL) {
    standardGeneric("setFeatures")
})

#' @export
setGeneric("projectData", function(projection = NULL, reference = NULL, cell_type_column = "cell_type1", 
    method = "scmap", threshold = 0.7) {
    standardGeneric("projectData")
})
