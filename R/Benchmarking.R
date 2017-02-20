#' Run support vector machines (\code{SVM}) prediction
#'
#' Train an \code{SVM} classifier on a training dataset (\code{train}) and then
#' classify a study dataset (\code{study}) using the classifier.
#'
#' @param train training data.frame with colnames, corresponding to training labels
#' @param study study data.frame
#' @param kern kernel to be used with SVM
#' @return classification of the study dataset
#'
#' @importFrom e1071 svm
#' @importFrom stats predict
#' @export
support_vector_machines <- function(train, study, kern = "linear") {
    train <- t(train)
    labs <- factor(rownames(train))
    rownames(train) <- NULL
    model <- tryCatch(e1071::svm(train, labs, kernel = kern), error = function(cond) return(NA))
    pred <- stats::predict(model, t(study))
    return(pred = pred)
}

#' Run random forest prediction
#'
#' Create a random forest calssifier based on a training dataset (\code{train}) and then
#' classify a study dataset (\code{study}) using the classifier.
#'
#' @param train training data.frame with colnames, corresponding to training labels
#' @param study study data.frame
#' @param ntree number of trees to be used in random forest
#' @return classification of the study dataset
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @export
random_forest <- function(train, study, ntree = 50) {
    train <- t(train)
    train <- as.data.frame(train)
    y <- as.factor(rownames(train))
    study <- t(study)
    study <- as.data.frame(study)
    rownames(train) <- NULL
    rownames(study) <- NULL
    train_rf <- randomForest::randomForest(x = train, y = y, ntree = ntree)
    Prediction <- stats::predict(train_rf, study)
    return(Prediction)
}

#' Bimodality is calculated using dip test: https://en.wikipedia.org/wiki/Multimodal_distribution#General_tests
#' @importFrom dplyr group_by summarise %>%
#' @importFrom reshape2 melt dcast
#' @importFrom scater pData<-
#' @importFrom proxy dist
#' @importFrom diptest dip.test
#' @importFrom fitdistrplus fitdist
#' @export
fsc3_map_by_similarity <- function(object_reference = NULL, object_to_map, buckets_ref = NULL, exprs_values = "exprs", similarity = "cosine", threshold = 1, suppress_plot = TRUE) {
    if(is.null(buckets_ref)) {
        # calculate median feature expression in every bucket of object_reference
        buckets <- object_reference@phenoData@data$fsc3_buckets
        if (is.null(buckets)) {
            warning(paste0("Please run fsc3_get_buckets() on your reference dataset first!"))
            return(object_to_map)
        }
        buckets_ref <- object_reference@assayData[[exprs_values]]
        buckets_ref <- buckets_ref[fData(object_reference)$fsc3_features, ]
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
            # norm_distr <- fitdistrplus::fitdist(mins, "norm")
            # x <- seq(0, 2, length = 200)
            # y <- dnorm(x, mean = norm_distr$estimate[["mean"]], sd = norm_distr$estimate[["sd"]])
            hist(mins, xlim = c(0, 2), freq = FALSE, xlab = "Normalised distance", ylab = "Density", main = "Distribution of normalised distances")
            # lines(x, y, type = "l", lwd = 3, col = "red")
        } else {
            mixmdl = mixtools::normalmixEM(mins)
            plot(mixmdl, which = 2, xlim = c(0, 2), xlab2 = "Normalised distance", main2 = "Distribution of normalised distances")
        }
    }
    buckets_assigned <- rownames(buckets_ref)[min_inds]
    buckets_assigned[mins > threshold] <- "unassigned"
    buckets_assigned[is.na(buckets_assigned)] <- "unassigned"
    p_data <- object_to_map@phenoData@data
    p_data$fsc3_buckets_assigned <- buckets_assigned
    pData(object_to_map) <- new("AnnotatedDataFrame", data = p_data)
    return(object_to_map)
}

