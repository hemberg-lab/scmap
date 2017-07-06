#' Single cell RNA-Seq data extracted from a publication by Yan et al.
#'
#' @source \url{http://dx.doi.org/10.1038/nsmb.2660}
#'
#' Columns represent cells, rows represent genes expression values.
#'
"yan"

#' Cell type annotations for data extracted from a publication by Yan et al.
#'
#' @source \url{http://dx.doi.org/10.1038/nsmb.2660}
#'
#' Each row corresponds to a single cell from `yan` dataset
#'
"ann"

#' Plot Sankey diagram comparing two clusterings
#' 
#' Sometimes it is useful to see how the clusters in two different clustering
#' solutions correspond to each other. Sankey diagram is a good way to visualize
#' them. This function takes as input two clustering solutions and visualizes them
#' using a Sankey diagram. The order of the reference clusters is defined by their
#' labels in increasing order.
#' 
#' @param reference reference clustering labels
#' @param clusters clustering labels under investigations
#' @param plot_width width of the output plot in pixels
#' @param plot_height height of the output plot in pixels
#' @param colors colors of the links between two clusterings. If defined please
#' note that each cluster in the reference clustering has to have its own color.
#' This should be a normal text vector, e.g. c('#FF0000', '#FFA500', '#008000')
#' 
#' @return an object returned by `gvisSankey`
#' 
#' @importFrom dplyr %>% summarise group_by
#' @importFrom reshape2 melt
#' @importFrom googleVis gvisSankey
#' 
#' @examples
#' plot(getSankey(ann[ , 1], ann[ , 1]))
#' 
#' @export
getSankey <- function(reference, clusters, plot_width = 400, plot_height = 600, colors = NULL) {
    Var1 <- value <- NULL
    res.all <- NULL
    for (j in names(table(reference))) {
        res <- NULL
        for (i in names(table(clusters))) {
            tmp <- length(intersect(which(clusters == i), which(reference == j)))
            res <- c(res, tmp)
        }
        res.all <- rbind(res.all, res)
    }
    colnames(res.all) <- names(table(clusters))
    rownames(res.all) <- names(table(reference))
    
    if (ncol(res.all) > 1) {
        res.all <- res.all[order(as.numeric(table(reference)), decreasing = TRUE), order(as.numeric(table(clusters)), 
            decreasing = TRUE), drop = FALSE]
    }
    
    res <- reshape2::melt(res.all)
    res <- res[res$value != 0, ]
    
    if (ncol(res.all) > 1) {
        maxs <- res %>% dplyr::group_by(Var1) %>% dplyr::summarise(max = max(value))
        
        res <- merge(res, maxs)
        maxs <- res[res$value == res$max, ]
        maxs <- maxs[order(maxs$value, decreasing = TRUE), ]
        res <- res[res$value != res$max, ]
        res <- rbind(maxs, res)
        res <- res[, 1:3]
    }
    
    # remove cycles from the data
    res[, 1] <- paste0(res[, 1], " ")
    res[, 2] <- paste0(" ", res[, 2])
    
    colnames(res) <- c("From", "To", "# of cells")
    
    if (!is.null(colors)) {
        colors <- paste(colors, collapse = "', '")
        colors <- paste0("['", colors, "']")
    }
    
    Sankey <- gvisSankey(res, from = "From", to = "To", weight = "# of cells", options = list(width = plot_width, height = plot_height, 
        sankey = paste0("{
                node:{
                    label:{
                        fontName:'Arial',
                        fontSize:11,color:
                        '#000000',
                        bold:true,
                        italic:false
                    },
                    colors:'#FFFFFF',
                    nodePadding:12
                },", 
            if (!is.null(colors)) {
                paste0("link:{
                    colorMode: 'source',
                    colors: ", colors, "
                },")
            }, "iterations:0
            }")))
    
    return(Sankey)
}

#' @importFrom e1071 svm
#' @importFrom stats predict
support_vector_machines <- function(train, study, kern = "linear") {
    train <- t(train)
    labs <- factor(rownames(train))
    rownames(train) <- NULL
    model <- tryCatch(e1071::svm(train, labs, kernel = kern, probability = TRUE), error = function(cond) return(NA))
    pred <- stats::predict(model, t(study), probability = TRUE)
    return(pred = pred)
}

#' @importFrom randomForest randomForest
#' @importFrom stats predict
random_forest <- function(train, study, ntree = 50) {
    train <- t(train)
    train <- as.data.frame(train)
    y <- as.factor(rownames(train))
    study <- t(study)
    study <- as.data.frame(study)
    rownames(train) <- NULL
    rownames(study) <- NULL
    train_rf <- randomForest::randomForest(x = train, y = y, ntree = ntree, keep.forest = TRUE)
    Prediction <- stats::predict(train_rf, study, type = "prob")
    return(Prediction)
}

#' @importFrom utils head
#' @importFrom stats lm
linearModel <- function(f_data, n_features) {
    # do not consider ERCC spike-ins and genes with 0 dropout rate
    dropouts_filter <- which(f_data$pct_dropout != 0 & f_data$pct_dropout != 100 & !grepl("ERCC-", f_data$feature_symbol))
    dropouts <- log2(f_data$pct_dropout[dropouts_filter])
    expression <- f_data$mean_exprs[dropouts_filter]
    feature_symbols <- f_data$feature_symbol[dropouts_filter]
    
    fit <- lm(dropouts ~ expression)
    gene_inds <- as.numeric(names(head(sort(fit$residuals, decreasing = TRUE), n_features)))
    
    scmap_features <- rep(FALSE, nrow(f_data))
    scmap_features[dropouts_filter[gene_inds]] <- TRUE
    
    scmap_scores <- rep(NA, nrow(f_data))
    scmap_scores[dropouts_filter] <- fit$residuals
    
    d <- as.data.frame(cbind(expression, dropouts))
    d$Gene <- f_data$feature_symbol[dropouts_filter]
    d$Features <- "Other"
    d$Features[gene_inds] <- "Selected"
    d$Features <- factor(d$Features, levels = c("Selected", "Other"))
    
    return(list(scmap_features = scmap_features, scmap_scores = scmap_scores, for_plotting = d, fit = fit))
}

#' Plot feature selection plot
#' 
#' @param object an object of \code{\link[scater]{SCESet}} class
#' @param n_features number of features to select
#' 
#' @return a ggplot object to plot the feature selection plot
#' 
#' @importFrom Biobase fData
plotFeatures <- function(object, n_features = 500) {
    f_data <- fData(object)
    tmp <- linearModel(f_data, n_features)
    return(ggplot_features(tmp$for_plotting, tmp$fit))
}

#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual labs geom_abline theme_classic
ggplot_features <- function(d, fit) {
    dropouts <- Features <- NULL
    cols <- c("#d73027", "#4575b4")
    p <- ggplot(d, aes(x = expression, y = dropouts, colour = Features)) + geom_point(size = 0.7) + scale_colour_manual(values = cols) + 
        labs(x = "log2(Expression)", y = "log2(% of dropouts)") + geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1]) + 
        theme_classic(base_size = 12)
    return(p)
}

prepareData <- function(reference, dat) {
    dat <- dat[order(rownames(dat)), ]
    reference <- reference[order(rownames(reference)), , drop = FALSE]
    reference <- reference[, colSums(reference) > 0, drop = FALSE]
    
    return(list(reference = reference, dat = dat))
}

#' @importFrom Biobase fData fData<- pData pData<- AnnotatedDataFrame
#' @importFrom scater get_exprs newSCESet calculateQCMetrics
mergeData <- function(object_reference, object_to_map) {
    if (!("SCESet" %in% is(object_reference)) | !("SCESet" %in% is(object_to_map))) {
        stop("Your arguments are not of `SCESet` class!")
    }
    features_reference <- fData(object_reference)$feature_symbol
    features_to_map <- fData(object_to_map)$feature_symbol
    common_features <- intersect(features_reference, features_to_map)
    common_features <- sort(common_features)
    
    dat_reference <- get_exprs(object_reference, "exprs")
    rownames(dat_reference) <- fData(object_reference)$feature_symbol
    dat_reference <- dat_reference[rownames(dat_reference) %in% common_features, ]
    dat_reference <- dat_reference[order(rownames(dat_reference)), ]
    
    dat_to_map <- get_exprs(object_to_map, "exprs")
    rownames(dat_to_map) <- fData(object_to_map)$feature_symbol
    dat_to_map <- dat_to_map[rownames(dat_to_map) %in% common_features, ]
    dat_to_map <- dat_to_map[order(rownames(dat_to_map)), ]
    
    res <- cbind(dat_reference, dat_to_map)
    colnames(res) <- 1:ncol(res)
    
    res_sceset <- newSCESet(exprsData = res, logExprsOffset = 1)
    
    if (is.null(pData(object_reference)$scmap_labs)) {
        ref_cell_type <- as.character(pData(object_reference)$cell_type1)
    } else {
        ref_cell_type <- as.character(pData(object_reference)$scmap_labs)
    }
    
    if (is.null(pData(object_to_map)$scmap_labs)) {
        map_cell_type <- as.character(pData(object_to_map)$cell_type1)
    } else {
        map_cell_type <- as.character(pData(object_to_map)$scmap_labs)
    }
    
    cell_type1 <- c(ref_cell_type, map_cell_type)
    
    p_data <- pData(res_sceset)
    p_data <- cbind(p_data, cell_type1)
    pData(res_sceset) <- AnnotatedDataFrame(p_data)
    
    f_data <- fData(res_sceset)
    f_data$feature_symbol <- rownames(res)
    fData(res_sceset) <- AnnotatedDataFrame(f_data)
    
    res_sceset <- calculateQCMetrics(res_sceset)
    
    return(res_sceset)
}

