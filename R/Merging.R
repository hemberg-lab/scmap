#' @importFrom scater newSCESet is_exprs<- calculateQCMetrics
fsc3_merge_datasets <- function(object_reference, object_to_map, exprs_values = "exprs", cell_type_name = "stage") {
    if(class(object_reference) != "SCESet" | class(object_to_map) != "SCESet") {
        warning("Your arguments are not of `SCESet` class!")
        return()
    }
    features_reference <- featureNames(object_reference)
    features_to_map <- featureNames(object_to_map)
    common_features <- intersect(features_reference, features_to_map)
    common_features <- sort(common_features)
    f_data <- object_reference@featureData@data
    scmap_features <- common_features %in% features_reference[f_data$fsc3_features]
    
    dat_reference <- object_reference@assayData[[exprs_values]]
    dat_reference <- dat_reference[rownames(dat_reference) %in% common_features, ]
    dat_reference <- dat_reference[order(rownames(dat_reference)), ]

    dat_to_map <- object_to_map@assayData[[exprs_values]]
    dat_to_map <- dat_to_map[rownames(dat_to_map) %in% common_features, ]
    dat_to_map <- dat_to_map[order(rownames(dat_to_map)), ]

    res <- cbind(dat_reference, dat_to_map)
    colnames(res) <- 1:ncol(res)
    
    res_sceset <- newSCESet(exprsData = res)
    
    
    cell_type <- c(as.character(object_reference@phenoData@data[[cell_type_name]]),
        as.character(dat_to_map@phenoData@data$fsc3_buckets_assigned))
    
    p_data <- res_sceset@phenoData@data
    p_data <- cbind(p_data, cell_type)
    pData(res_sceset) <- new("AnnotatedDataFrame", data = p_data)
    
    res_sceset <- res_sceset[, p_data$cell_type != "unassigned"]
    is_exprs(res_sceset) <- dataset <- res_sceset@assayData[["exprs"]] > 0
    res_sceset <- calculateQCMetrics(res_sceset)
    
    f_data <- res_sceset@featureData@data
    f_data$fsc3_features <- scmap_features
    fData(res_sceset) <- new("AnnotatedDataFrame", data = f_data)
    
    return(res_sceset)
}
