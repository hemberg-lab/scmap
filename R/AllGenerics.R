#' @export
setGeneric("getFeatures", function(object, n_features, pct_dropout_min, pct_dropout_max, 
    suppress_plot) {
    standardGeneric("getFeatures")
})

#' @export
setGeneric("setFeatures", function(object, features) {
    standardGeneric("setFeatures")
})

