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
    
    Sankey <- gvisSankey(res, from = "From", to = "To", weight = "# of cells", options = list(width = plot_width, 
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
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assayNames
#' @importFrom BiocGenerics counts
linearModel <- function(object, n_features) {
    log_count <- as.matrix(logcounts(object))
    cols <- ncol(log_count)
    if (!"counts" %in% assayNames(object)) {
        warning("Your object does not contain counts() slot. Dropouts were calculated using logcounts() slot...")
        dropouts <- rowSums(log_count == 0)/cols * 100
    } else {
        count <- as.matrix(counts(object))
        dropouts <- rowSums(count == 0)/cols * 100
    }
    # do not consider genes with 0 and 100 dropout rate
    dropouts_filter <- dropouts != 0 & dropouts != 100
    dropouts_filter <- which(dropouts_filter)
    dropouts <- log2(dropouts[dropouts_filter])
    expression <- rowSums(log_count[dropouts_filter, ])/cols
    
    fit <- lm(dropouts ~ expression)
    gene_inds <- fit$residuals
    names(gene_inds) <- 1:length(gene_inds)
    gene_inds <- as.numeric(names(head(sort(gene_inds, decreasing = TRUE), n_features)))
    
    scmap_features <- rep(FALSE, nrow(object))
    scmap_features[dropouts_filter[gene_inds]] <- TRUE
    
    scmap_scores <- rep(NA, nrow(object))
    scmap_scores[dropouts_filter] <- fit$residuals
    
    d <- as.data.frame(cbind(expression, dropouts))
    d$Gene <- rownames(object)[dropouts_filter]
    d$Features <- "Other"
    d$Features[gene_inds] <- "Selected"
    d$Features <- factor(d$Features, levels = c("Selected", "Other"))
    
    return(list(scmap_features = scmap_features, scmap_scores = scmap_scores, for_plotting = d, fit = fit))
}

#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual labs geom_abline theme_classic
ggplot_features <- function(d, fit) {
    dropouts <- Features <- NULL
    cols <- c("#d73027", "#4575b4")
    p <- ggplot(d, aes(x = expression, y = dropouts, colour = Features)) + geom_point(size = 0.7) + 
        scale_colour_manual(values = cols) + labs(x = "log2(Expression)", y = "log2(% of dropouts)") + 
        geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1]) + theme_classic(base_size = 12)
    return(p)
}

checks_for_index <- function(object) {
  if (is.null(object)) {
    stop("Please provide a `SingleCellExperiment` object using the `object` parameter!")
    return(FALSE)
  }
  if (!"SingleCellExperiment" %in% is(object)) {
    stop("Input object is not of `SingleCellExperiment` class! Please provide an object of the correct class!")
    return(FALSE)
  }
  if (is.null(rowData(object)$scmap_features)) {
    stop("Features are not selected! Please run `selectFeatures()` or `setFeatures()` first!")
    return(FALSE)
  }
  if (is.null(rowData(object)$feature_symbol)) {
    stop("There is no `feature_symbol` column in the `rowData` slot of the `reference` dataset! Please write your gene/transcript names to this column!")
    return(FALSE)
  }
  return(TRUE)
}

checks_for_projection <- function(projection, index_list) {
  if (is.null(projection)) {
    stop("Please provide a `SingleCellExperiment` object for the `projection` parameter!")
    return(FALSE)
  }
  if (!"SingleCellExperiment" %in% is(projection)) {
    stop("`projection` dataset has to be of the `SingleCellExperiment` class!")
    return(FALSE)
  }
  if (is.null(rowData(projection)$feature_symbol)) {
    stop("There is no `feature_symbol` column in the `rowData` slot of the `projection` dataset! Please write your gene/transcript names to this column!")
    return(FALSE)
  }
  if (is.null(index_list)) {
    stop("Please provide a list of precomputed scmap indexes as the `reference` parameter!")
    return(FALSE)
  }
  if (!"list" %in% is(index_list)) {
    stop("Please provide a list of precomputed scmap indexes as the `reference` parameter!")
    return(FALSE)
  }
  return(TRUE)
}

dists_subcentroids <- function(proj_exprs, subcentroids) {
  features_query <- rownames(proj_exprs)
  num_cells <- ncol(proj_exprs)
  SqNorm <- numeric(ncol(proj_exprs))
  query_chunks <- list()
  for (m in seq_len(length(subcentroids))) {
    subcentroids_chunk <- subcentroids[[m]]
    features_chunk <- rownames(subcentroids_chunk)
    common_features <- intersect(features_chunk, features_query)
    if (length(common_features) == 0) {
      # change to a more memory-efficient method later?
      query_chunks[[m]] <- matrix(num_cells, 1, num_cells)
      subcentroids[[m]] <- matrix(ncol(subcentroids_chunk), 1, ncol(subcentroids_chunk))
    } else {
      common_features <- sort(common_features)
      subcentroids[[m]] <- subcentroids_chunk[rownames(subcentroids_chunk) %in% common_features, 
                                              , drop = FALSE]
      query_chunks[[m]] <- proj_exprs[rownames(proj_exprs) %in% common_features, , drop = FALSE]
      if (length(common_features) > 1) {
        subcentroids[[m]] <- subcentroids[[m]][order(rownames(subcentroids[[m]])), ]
        query_chunks[[m]] <- query_chunks[[m]][order(rownames(query_chunks[[m]])), ]
      }
      # find the squared Euclidean norm of every query after selecting features
      SqNorm <- SqNorm + EuclSqNorm(query_chunks[[m]])
    }
  }
  return(list(subcentroids = subcentroids, query_chunks = query_chunks, SqNorm = SqNorm))
}

order_and_combine_labels <- function(labels, simls) {
  unassigned_rate_order <- order(
    unlist(
      lapply(labels, function(x) {
        length(x[x == "unassigned"])/length(x)
      })
    )
  )
  labels <- labels[unassigned_rate_order]
  simls <- simls[unassigned_rate_order]
  labels <- do.call(cbind, labels)
  simls <- do.call(cbind, simls)
  max_simls_inds <- apply(
    simls, 1, function(x) {
      if(!all(is.na(x))) {
        return(which.max(x))
      } else {
        return(NA)
      }
    }
  )
  inds <- which(!is.na(max_simls_inds))
  cons_labels <- rep("unassigned", nrow(labels))
  cons_labels[inds] <- labels[cbind(inds, max_simls_inds[inds])]
  return(list(scmap_cluster_labs = labels, scmap_cluster_siml = simls, combined_labs = cons_labels))
}