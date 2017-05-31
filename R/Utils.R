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
#' @importFrom dplyr %>% summarise group_by
#' @importFrom reshape2 melt
#' @importFrom googleVis gvisSankey
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
    
    res.all <- res.all[order(as.numeric(table(reference)), decreasing = T), order(as.numeric(table(clusters)), 
        decreasing = T)]
    
    res <- reshape2::melt(res.all)
    res <- res[res$value != 0, ]
    
    maxs <- res %>% dplyr::group_by(Var1) %>% dplyr::summarise(max = max(value))
    
    res <- merge(res, maxs)
    maxs <- res[res$value == res$max, ]
    maxs <- maxs[order(maxs$value, decreasing = T), ]
    res <- res[res$value != res$max, ]
    res <- rbind(maxs, res)
    res <- res[, 1:3]
    
    # remove cycles from the data
    res[, 1] <- paste0(res[, 1], " ")
    res[, 2] <- paste0(" ", res[, 2])
    
    colnames(res) <- c("From", "To", "Weight")
    
    if (!is.null(colors)) {
        colors <- paste(colors, collapse = "', '")
        colors <- paste0("['", colors, "']")
    }
    
    Sankey <- gvisSankey(res, from = "From", to = "To", weight = "Weight", options = list(width = plot_width, 
        height = plot_height, sankey = paste0("{
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
                    colors: ", 
                  colors, "
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
    pred <- stats::predict(model, t(study), probability=TRUE)
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

