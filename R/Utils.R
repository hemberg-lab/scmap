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
#' @param colors colors of the links between two clusterings. If defined please
#' note that each cluster in the reference clustering has to have its own color.
#' This should be a normal text vector, e.g. c("#FF0000", "#FFA500", "#008000")
#' 
#' @importFrom dplyr %>% summarise group_by
#' @importFrom reshape2 melt
#' @importFrom googleVis gvisSankey
#' 
#' @export
plot_sankey <- function(reference, clusters, plot_width = 400, plot_height = 600, colors = NULL) {
    res.all <- NULL
    for(j in names(table(reference))){
        res <- NULL
        for(i in names(table(clusters))) {
            tmp <- length(intersect(which(clusters == i), which(reference == j)))
            res <- c(res, tmp)
        }
        res.all <- rbind(res.all, res)
    }
    colnames(res.all) <- names(table(clusters))
    rownames(res.all) <- names(table(reference))
    
    res.all <- res.all[order(as.numeric(table(reference)), decreasing = T), 
                       order(as.numeric(table(clusters)), decreasing = T)]
    
    res <- reshape2::melt(res.all)
    res <- res[res$value != 0, ]
    
    maxs <- res %>%
        dplyr::group_by(Var1) %>%
        dplyr::summarise(max = max(value))
    
    res <- merge(res, maxs)
    maxs <- res[res$value == res$max, ]
    maxs <- maxs[order(maxs$value, decreasing = T), ]
    res <- res[res$value != res$max, ]
    res <- rbind(maxs, res)
    res <- res[,1:3]
    
    # remove cycles from the data
    res[, 1] <- paste0(res[, 1], " ")
    res[, 2] <- paste0(" ", res[, 2])
    
    colnames(res) <- c("From", "To", "Weight")
    
    if(!is.null(colors)) {
        colors <- paste(colors, collapse = "', '")
        colors <- paste0("['", colors, "']")
    }
    
    Sankey <- gvisSankey(
        res,
        from="From",
        to="To",
        weight="Weight",
        options = list(
            width = plot_width,
            height = plot_height,
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
                },", if(!is.null(colors)) {
                    paste0("link:{
                    colorMode: 'source',
                    colors: ", colors, "
                },")},
                    "iterations:0
            }"
            ))
    )
    
    plot(Sankey)
}
