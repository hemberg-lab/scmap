#' Bimodality is calculated using dip test: https://en.wikipedia.org/wiki/Multimodal_distribution#General_tests
#' 
#' @param object_reference reference SCESet set
#' @param object_to_map SCESet to map
#' @param buckets_ref reference cell buckets
#' @param exprs_values expression values to use
#' @param similarity similarity measure
#' @param threshold threshold on similarity
#' @param suppress_plot defines whether to suppress an output plot
#' 
#' @importFrom dplyr group_by summarise %>%
#' @importFrom reshape2 melt dcast
#' @importFrom scater pData<-
#' @importFrom proxy dist
#' @importFrom diptest dip.test
#' @importFrom graphics abline hist plot points
#' @importFrom stats median
#' @importFrom utils head
#' @importFrom mixtools normalmixEM
#' @importFrom nnet which.is.max
#' @export
scmap <- function(object_reference = NULL, object_to_map, buckets_ref = NULL, exprs_values = "exprs", similarity = "cosine", threshold = 1, suppress_plot = TRUE) {
    
    gene <- bucket <- exprs <- NULL
    
    if(is.null(buckets_ref)) {
        # calculate median feature expression in every bucket of object_reference
        buckets <- object_reference@phenoData@data$scmap_buckets
        if (is.null(buckets)) {
            warning(paste0("Please run fsc3_get_buckets() on your reference dataset first!"))
            return(object_to_map)
        }
        buckets_ref <- object_reference@assayData[[exprs_values]]
        f_data <- object_reference@featureData@data
        buckets_ref <- buckets_ref[f_data$fsc3_features, ]
        colnames(buckets_ref) <- buckets
        buckets_ref <- reshape2::melt(buckets_ref)
        colnames(buckets_ref) <- c("gene", "bucket", "exprs")
        buckets_ref <- buckets_ref %>%
            group_by(gene, bucket) %>%
            summarise(
                med_exprs = median(exprs)
            )
        buckets_ref <- reshape2::dcast(buckets_ref, gene ~ bucket, value.var = "med_exprs")
        rownames(buckets_ref) <- buckets_ref$gene
        buckets_ref <- buckets_ref[ , 2:ncol(buckets_ref)]
    }
    
    dat <- object_to_map@assayData[[exprs_values]]
    dat <- dat[rownames(dat) %in% rownames(buckets_ref), ]
    dat <- dat[order(rownames(dat)), ]
    
    buckets_ref <- buckets_ref[rownames(buckets_ref) %in% rownames(dat), ]
    buckets_ref <- buckets_ref[order(rownames(buckets_ref)), ]
    buckets_ref <- buckets_ref[,colSums(buckets_ref) > 0]
    
    dat <- scale(dat, center = TRUE, scale = TRUE)
    buckets_ref <- scale(buckets_ref, center = TRUE, scale = TRUE)
    dat <- t(dat)
    buckets_ref <- t(buckets_ref)
    if(similarity == "euclidean") {
        res <- proxy::dist(buckets_ref, dat, method = "euclidean")
        res <- matrix(res, ncol = nrow(buckets_ref), byrow = T)
        res <- res / max(res, na.rm = T) * 2
    }
    if(similarity == "cosine") {
        res <- proxy::dist(buckets_ref, dat, method = "cosine")
        res <- matrix(res, ncol = nrow(buckets_ref), byrow = T)
    }
    min_inds <- unlist(apply(-res, 1, nnet::which.is.max))
    mins <- unlist(apply(res, 1, min))
    if(!suppress_plot) {
        unimod <- diptest::dip.test(mins)$p.value >= 0.05
        if(unimod) {
            hist(mins, xlim = c(0, 2), freq = FALSE, xlab = "Normalised distance", ylab = "Density", main = "Distribution of normalised distances")
        } else {
            mixmdl = mixtools::normalmixEM(mins)
            plot(mixmdl, which = 2, xlim = c(0, 2), xlab2 = "Normalised distance", main2 = "Distribution of normalised distances")
        }
    }
    buckets_assigned <- rownames(buckets_ref)[min_inds]
    buckets_assigned[mins > threshold] <- "unassigned"
    buckets_assigned[is.na(buckets_assigned)] <- "unassigned"
    p_data <- object_to_map@phenoData@data
    p_data$scmap_buckets_assigned <- buckets_assigned
    pData(object_to_map) <- new("AnnotatedDataFrame", data = p_data)
    return(object_to_map)
}
