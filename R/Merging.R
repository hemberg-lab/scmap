#' @importFrom scater newSCESet is_exprs<- calculateQCMetrics
#' @export
mergeData <- function(object_reference, object_to_map) {
    if (class(object_reference) != "SCESet" | class(object_to_map) != "SCESet") {
        warning("Your arguments are not of `SCESet` class!")
        return()
    }
    features_reference <- object_reference@featureData@data$feature_symbol
    features_to_map <- object_to_map@featureData@data$feature_symbol
    common_features <- intersect(features_reference, features_to_map)
    common_features <- sort(common_features)
    
    dat_reference <- object_reference@assayData[["exprs"]]
    rownames(dat_reference) <- object_reference@featureData@data$feature_symbol
    dat_reference <- dat_reference[rownames(dat_reference) %in% common_features, ]
    dat_reference <- dat_reference[order(rownames(dat_reference)), ]
    dat_reference <- dat_reference[!duplicated(rownames(dat_reference)), ]
    
    dat_to_map <- object_to_map@assayData[["exprs"]]
    rownames(dat_to_map) <- object_to_map@featureData@data$feature_symbol
    dat_to_map <- dat_to_map[rownames(dat_to_map) %in% common_features, ]
    dat_to_map <- dat_to_map[order(rownames(dat_to_map)), ]
    dat_to_map <- dat_to_map[!duplicated(rownames(dat_to_map)), ]
    
    res <- cbind(dat_reference, dat_to_map)
    colnames(res) <- 1:ncol(res)
    
    res_sceset <- newSCESet(exprsData = res, logExprsOffset = 1)
    
    if(is.null(object_reference@phenoData@data$scmap_labs)) {
        ref_cell_type <- as.character(object_reference@phenoData@data$cell_type1)
    } else {
        ref_cell_type <- as.character(object_reference@phenoData@data$scmap_labs)
    }
    
    if(is.null(object_to_map@phenoData@data$scmap_labs)) {
        map_cell_type <- as.character(object_to_map@phenoData@data$cell_type1)
    } else {
        map_cell_type <- as.character(object_to_map@phenoData@data$scmap_labs)
    }
    
    cell_type1 <- c(ref_cell_type, map_cell_type)
    
    p_data <- res_sceset@phenoData@data
    p_data <- cbind(p_data, cell_type1)
    pData(res_sceset) <- new("AnnotatedDataFrame", data = p_data)
    
    f_data <- res_sceset@featureData@data
    f_data$feature_symbol <- rownames(res)
    fData(res_sceset) <- new("AnnotatedDataFrame", data = f_data)
    
    # res_sceset <- res_sceset[, p_data$cell_type != "unassigned"]
    res_sceset <- calculateQCMetrics(res_sceset)
    
    return(res_sceset)
}
