
get_common_fsc3_features <- function(objects) {
    common_features <- NULL
    for(object in objects) {
        f_data <- object@featureData@data
        if(is.null(f_data$feature_symbol)) {
            message("Every input object has to contain 'feature_symbol' column in the featureData slot! Please check your input objects!")
            return(NULL)
        }
        object_features <- as.character(f_data$feature_symbol)[f_data$fsc3_features]
        common_features <- c(common_features, object_features)
    }
    common_features <- unique(common_features)
    common_inds <- rep(TRUE, length(common_features))
    for(object in objects) {
        f_data <- object@featureData@data
        inds <- match(common_features, f_data$feature_symbol)
        common_inds[is.na(inds)] <- FALSE
    }
    return(common_features[common_inds])
}

#' @importFrom scater pData<-
fsc3_assign_signatures <- function(object_to_assign, object_ref, threshold = 0.7) {
    # reference object
    buckets_ref <- unique(object_ref@phenoData@data[, c("fsc3_buckets", "fsc3_buckets_signatures")])
    # convert factors to strings
    if(is.factor(buckets_ref$fsc3_buckets)) {
        buckets_ref$fsc3_buckets <- levels(buckets_ref$fsc3_buckets)[buckets_ref$fsc3_buckets]
    }
    # object to assign
    sigs_to_assign <- object_to_assign@phenoData@data$fsc3_signatures

    buckets_assigned <- NULL
    for(sig in sigs_to_assign) {
        sig <- unlist(strsplit(sig, ""))
        common_bits <- c()
        for(sig_ref in buckets_ref$fsc3_buckets_signatures) {
            tmp <- unlist(strsplit(sig_ref, ""))
            common_bits <- c(common_bits, length(which(tmp == sig))/length(which(tmp != "_")))
        }
        if(max(common_bits) >= threshold) {
            tmp <- nnet::which.is.max(common_bits)
            buckets_assigned <- c(buckets_assigned, buckets_ref$fsc3_buckets[tmp])
        } else {
            buckets_assigned <- c(buckets_assigned, "unassigned")
        }
    }

    p_data <- object_to_assign@phenoData@data
    p_data$fsc3_buckets_assigned <- buckets_assigned
    pData(object_to_assign) <- new("AnnotatedDataFrame", data = p_data)

    return(object_to_assign)
}

fsc3_plot_sig_expression <- function(object, fColumn = NULL, n_cells = 10, hc = NULL) {
    if(is.null(fColumn)) {
        warning("Please provide the fColumn argument!")
        return(object)
    }
    p_data <- object@phenoData@data
    sigs <- p_data$fsc3_signatures[order(p_data[[fColumn]])]
    mat <- do.call(cbind, lapply(strsplit(sigs, ""), as.numeric))
    colnames(mat) <- p_data[[fColumn]][order(p_data[[fColumn]])]
    to_plot <- NULL
    gaps <- NULL
    gaps_it <- 0
    for(i in unique(colnames(mat))) {
        tmp <- mat[,colnames(mat) == i]
        if(ncol(tmp) > n_cells) {
            tmp <- tmp[, sample(1:ncol(tmp), n_cells)]
        }
        colnames(tmp) <- rep(i, ncol(tmp))
        to_plot <- cbind(to_plot, tmp)
        gaps <- c(gaps, gaps_it + ncol(tmp))
        gaps_it <- gaps_it + ncol(tmp)
    }
    if(is.null(hc)) {
        phm <- pheatmap::pheatmap(to_plot, cluster_cols = F, gaps_col = gaps)
    } else {
        phm <- pheatmap::pheatmap(to_plot, cluster_cols = F, cluster_rows = hc, gaps_col = gaps)
    }
}

