library(dplyr)
library(ggplot2)

Plot.Interac.gain.loss <- function(final.plot.data, size_n) { 
  # Define Nature journal-style color palette for axis text
  nature_colors <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B2DF8A"
  )
  
  # Original color palette for points (unchanged)
  color.use <- c('#3288BD', '#66C2A5', '#ABDDA4', '#E6F598', '#FEE08B', 
                 '#FDAE61', '#F46D43', '#D53E4F')
  
  # Prepare data
  final.plot.data <- final.plot.data %>% arrange(desc(dataset))
  final.plot.data$source.target <- as.character(final.plot.data$source.target)
  final.plot.data$group.names <- factor(final.plot.data$group.names, levels = sort(unique(final.plot.data$group.names)))
  
  # Map colors to unique 'source' values for axis text
  unique_sources <- unique(final.plot.data$source)
  if (length(unique_sources) > length(nature_colors)) {
    warning("More unique sources than colors available; recycling colors.")
    source_colors <- setNames(rep(nature_colors, length.out = length(unique_sources)), unique_sources)
  } else {
    source_colors <- setNames(nature_colors[1:length(unique_sources)], unique_sources)
  }
      label_colors <- source_colors[as.character(final.plot.data$source)[match(levels(final.plot.data$group.names), final.plot.data$group.names)]]
  # Create the base plot
  g <- ggplot(final.plot.data, aes(x = group.names, y = interaction_name_2, 
                                   color = prob, size = pval)) + 
    geom_point(pch = 16) + 
    theme_linedraw() + 
    theme(panel.grid.major = element_blank(),
          # axis.text.x = element_text(angle = 0),
          
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,
                                     color = label_colors),  # Apply Nature-style colors to x-axis text
          axis.title.x = element_blank(),
          
          axis.title.y = element_blank()) + 
    scale_x_discrete(position = "bottom") +
    scale_radius(range = c(min(final.plot.data$pval) + size_n, max(final.plot.data$pval) + size_n), 
                 breaks = sort(unique(final.plot.data$pval)), 
                 labels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")[sort(unique(final.plot.data$pval))], 
                 name = "p-value") +
    scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                           na.value = "white", 
                           limits = c(quantile(final.plot.data$prob, 0, na.rm = TRUE), 
                                      quantile(final.plot.data$prob, 1, na.rm = TRUE)), 
                           breaks = c(quantile(final.plot.data$prob, 0, na.rm = TRUE), 
                                      quantile(final.plot.data$prob, 1, na.rm = TRUE))
                           # , labels = c(round(quantile(final.plot.data$prob, 0, na.rm = TRUE), 3), round(quantile(final.plot.data$prob, 1, na.rm = TRUE),3))
                           , labels = c("min", "max")
                           
                          ) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
    theme(text = element_text(size = 12), 
          plot.title = element_text(size = 12)) +
    theme(legend.title = element_text(size = 8), 
          legend.text = element_text(size = 6), 
          strip.text = element_text(face = "bold", size = 12, color = "white")) +
    facet_grid(.~dataset)
  
  # Add vertical lines if multiple source.target values
  if (length(unique(final.plot.data$source.target)) > 1) {
    g <- g + geom_vline(xintercept = seq(1.5, length(unique(final.plot.data$source.target)) - 0.5, 1), 
                        lwd = 0.1, colour = "#c8c8c8")
  }
  
  # Add horizontal lines if multiple interaction_name_2 values
  if (length(unique(final.plot.data$interaction_name_2)) > 1) {
    g <- g + geom_hline(yintercept = seq(1.5, length(unique(final.plot.data$interaction_name_2)) - 0.5, 1), 
                        lwd = 0.1, colour = "#c8c8c8")
  }
  
  return(g)
}


getdata <- function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
    pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, 
    sort.by.source.priority = TRUE, thresh = 0.05, comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, 
    min.dataset = NULL, min.quantile = 0, max.quantile = 1, return.data = FALSE) 
{
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net", 
        sources.use = sources.use, targets.use = targets.use, 
        signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
        cells.level <- levels(object@idents[[comparison[ii]]])
        df.net <- df.net.all[[comparison[ii]]]
        df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
        df.net$source.target <- paste(df.net$source, df.net$target, 
            sep = " -> ")
        source.target <- paste(rep(sources.use, each = length(targets.use)), 
            targets.use, sep = " -> ")
        source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
        if (length(source.target.isolate) > 0) {
            df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                ncol = ncol(df.net)))
            colnames(df.net.isolate) <- colnames(df.net)
            df.net.isolate$source.target <- source.target.isolate
            df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
            df.net.isolate$pval <- 1
            a <- stringr::str_split(df.net.isolate$source.target, 
                " -> ", simplify = T)
            df.net.isolate$source <- as.character(a[, 1])
            df.net.isolate$target <- as.character(a[, 2])
            df.net <- rbind(df.net, df.net.isolate)
        }
        df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
            unique(df.net$source)])
        df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
            unique(df.net$target)])
        group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
            levels(df.net$target), sep = " -> ")
        group.names0 <- group.names
        group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
            ")")
        df.net$p.value = df.net$pval
        if (nrow(df.net) > 0) {
            df.net$pval[df.net$pval > 0.05] = 1
            df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
            df.net$pval[df.net$pval <= 0.01] = 3
            df.net$prob[df.net$prob == 0] <- NA
            df.net$prob.original <- df.net$prob
            df.net$prob <- -1/log(df.net$prob) ############
        }
        else {
            df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                ncol = 5))
            colnames(df.net) <- c("interaction_name_2", "source.target", 
                "prob", "pval", "prob.original")
            df.net$source.target <- group.names0
        }
        df.net$group.names <- as.character(df.net$source.target)
        df.net$source.target <- paste0(df.net$source.target, 
            " (", dataset.name[comparison[ii]], ")")
        df.net$dataset <- dataset.name[comparison[ii]]
        df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
        stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 0)
    if (sum(idx1) > 0) {
        values.assign <- seq(max(df.all$prob, na.rm = T) * 1.1, 
            max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
        position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
        df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
            position)]
    }
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2), ])
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
        for (j in 1:length(comparison)) {
            cells.order <- c(cells.order, paste0(group.names0[i], 
                " (", dataset.name[comparison[j]], ")"))
            dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
        }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
    min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
    max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
    df$prob[df$prob < min.cutoff] <- min.cutoff
    df$prob[df$prob > max.cutoff] <- max.cutoff
    if (remove.isolate) {
        df <- df[!is.na(df$prob), ]
    }
    if (!is.null(max.dataset)) {
        signaling <- as.character(unique(df$interaction_name_2))
        for (i in signaling) {
            df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
            cell <- as.character(unique(df.i$group.names))
            for (j in cell) {
                df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
                values <- df.i.j$prob
                idx.max <- which(values == max(values, na.rm = T))
                idx.min <- which(values == min(values, na.rm = T))
                dataset.na <- c(df.i.j$dataset[is.na(values)], 
                  setdiff(dataset.name[comparison], df.i.j$dataset))
                if (length(idx.max) > 0) {
                  if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
                    df.i.j$prob <- NA
                  }
                  else if ((idx.max != idx.min) & !is.null(min.dataset)) {
                    if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
                      df.i.j$prob <- NA
                    }
                    else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
                      dataset.na)) > 0) {
                      df.i.j$prob <- NA
                    }
                  }
                }
                df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
            }
            df[df$interaction_name_2 == i, "prob"] <- df.i$prob
        }
    }
    if (remove.isolate) {
        df <- df[!is.na(df$prob), ]
    }
    if (nrow(df) == 0) {
        stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    if (!is.null(pairLR.use)) {
        interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name, 
            ]$interaction_name_2, unique(df$interaction_name_2))
        df$interaction_name_2 <- factor(df$interaction_name_2, 
            levels = interaction_name_2.order)
    }
    df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
        unique(df$source.target)))
    if (sort.by.target & !sort.by.source) {
        if (!is.null(targets.use)) {
            df$target <- factor(df$target, levels = intersect(targets.use, 
                df$target))
            df <- with(df, df[order(target, source), ])
            source.target.order <- unique(as.character(df$source.target))
            df$source.target <- factor(df$source.target, levels = source.target.order)
        }
    }
    if (sort.by.source & !sort.by.target) {
        if (!is.null(sources.use)) {
            df$source <- factor(df$source, levels = intersect(sources.use, 
                df$source))
            df <- with(df, df[order(source, target), ])
            source.target.order <- unique(as.character(df$source.target))
            df$source.target <- factor(df$source.target, levels = source.target.order)
        }
    }
    if (sort.by.source & sort.by.target) {
        if (!is.null(sources.use)) {
            df$source <- factor(df$source, levels = intersect(sources.use, 
                df$source))
            if (!is.null(targets.use)) {
                df$target <- factor(df$target, levels = intersect(targets.use, 
                  df$target))
            }
            if (sort.by.source.priority) {
                df <- with(df, df[order(source, target), ])
            }
            else {
                df <- with(df, df[order(target, source), ])
            }
            source.target.order <- unique(as.character(df$source.target))
            df$source.target <- factor(df$source.target, levels = source.target.order)
        }
    }
    return(df)
}

netAnalysis_signalingRole_heatmap_new <- function (object=NULL, signaling = NULL, pattern = "all", slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", 
    title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, centr = NULL,
    cluster.rows = FALSE, cluster.cols = FALSE) 
{
    # pattern <- match.arg(pattern)
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
        outgoing[, i] <- centr[[i]]$outdeg
        incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
        mat <- t(outgoing)
        legend.name <- "Outgoing"
    }
    else if (pattern == "incoming") {
        mat <- t(incoming)
        legend.name <- "Incoming"
    }
    else if (pattern == "all") {
        mat <- t(outgoing + incoming)
        legend.name <- "Overall"
    }
    if (is.null(title)) {
        title <- paste0(legend.name, " signaling patterns")
    }
    else {
        title <- paste0(paste0(legend.name, " signaling patterns"), 
            " - ", title)
    }
    if (!is.null(signaling)) {
        mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
        mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
        idx <- match(rownames(mat1), signaling)
        mat[idx[!is.na(idx)], ] <- mat1
        dimnames(mat) <- list(signaling, colnames(mat1))
    }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    mat[mat == 0] <- NA
    
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    color.use <- color.use[colnames(mat)]
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
        which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
        simple_anno_size = grid::unit(0.2, "cm"))
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
        border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
        show_annotation_name = FALSE)
    pSum <- rowSums(mat.ori)
    pSum.original <- pSum
    pSum <- -1/log(pSum)
    pSum[is.na(pSum)] <- 0
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    if (length(idx1) > 0) {
        values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
            length.out = length(idx1))
        position <- sort(pSum.original[idx1], index.return = TRUE)$ix
        pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    # ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
    #     show_annotation_name = FALSE)
    if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
        legend.break <- max(mat, na.rm = T)
    }
    else {
        legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
            round(max(mat, na.rm = T), digits = 1))
    }
    
    library(ComplexHeatmap)
    library(circlize)  
    col_fun = colorRamp2(c(1, 0), c("#1F77B4", "white"))
    mat = replace(mat, is.na(mat) | is.nan(mat), 0)

    mat <- mat[sort(rownames(mat)),sort(colnames(mat))]
 
    
    ht1 = Heatmap(mat, 
                  col = colorRamp2(c(1, 0), c("#1F77B4", "white")),
                  na_col = "white", 
                  name = "Relative strength", 
                  bottom_annotation = col_annotation,
                  top_annotation = ha2,
                  # right_annotation = ha1,
                  cluster_rows = cluster.rows,
                  cluster_columns = cluster.rows,
                  row_names_side = "left",
                  row_names_rot = 0, 
                  row_names_gp = gpar(fontsize = font.size), 
                  column_names_gp = gpar(fontsize = font.size), 
                  width = unit(width, "cm"), 
                  height = unit(height, "cm"), 
                  column_title = title, 
                  row_order = sort(rownames(mat)),
                  column_title_gp = gpar(fontsize = font.size.title), 
                  column_names_rot = 45,  # x 轴标签 45 度
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 10,
                                                              fontface = "plain"),
                                              title_position = "leftcenter-rot",
                                              border = NA, 
                                              at = c(0, 0.5, 1),  # 图例刻度对应 0 到 1
                                              legend_height = unit(20, "mm"), 
                                              labels_gp = gpar(fontsize = 8), 
                                              grid_width = unit(2, "mm")))
    return(ht1)
}

centr_data <- function (object = NULL, slot.name = "netP", net = NULL, net.name = NULL, 
    thresh = 0.05, final.plot.data = NULL) 
{
    if (is.null(net)) {
        prob <- methods::slot(object, slot.name)$prob
        pval <- methods::slot(object, slot.name)$pval
        pval[prob == 0] <- 1
        prob[pval >= thresh] <- 0
        net = prob
    
        new.net <- net  
        new.net[] <- 0

        for(r in 1:nrow(final.plot.data)){      
        new.net[final.plot.data[r,"source"], 
                final.plot.data[r,"target"], 
                final.plot.data[r,"pathway_name"]] = net[final.plot.data[r,"source"],  final.plot.data[r,"target"], final.plot.data[r,"pathway_name"]] 
        
        }
    }
    if (is.null(net.name)) {
        net.name <- dimnames(new.net)[[3]]
    }
    if (length(dim(new.net)) == 3) {
        nrun <- dim(new.net)[3]
        my.sapply <- ifelse(test = future::nbrOfWorkers() == 
            1, yes = pbapply::pbsapply, no = future.apply::future_sapply)
        centr.all = my.sapply(X = 1:nrun, FUN = function(x) {
            net0 <- new.net[, , x]
            return(computeCentralityLocal(net0))
        }, simplify = FALSE)
    }
    else {
        centr.all <- as.list(computeCentralityLocal(new.net))
    }
    names(centr.all) <- net.name

    return(centr.all)
}

computeCentralityLocal <- function(net) {
  centr <- vector("list")
  G <- igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  centr$outdeg_unweighted <- rowSums(net > 0)
  centr$indeg_unweighted <- colSums(net > 0)
  centr$outdeg <- igraph::strength(G, mode="out")
  centr$indeg <- igraph::strength(G, mode="in")
  centr$hub <- igraph::hub_score(G)$vector
  centr$authority <- igraph::authority_score(G)$vector # A node has high authority when it is linked by many other nodes that are linking many other nodes.
  centr$eigen <- igraph::eigen_centrality(G)$vector # A measure of influence in the network that takes into account second-order connections
  centr$page_rank <- igraph::page_rank(G)$vector
  igraph::E(G)$weight <- 1/igraph::E(G)$weight
  centr$betweenness <- igraph::betweenness(G)
  #centr$flowbet <- try(sna::flowbet(net)) # a measure of its role as a gatekeeper for the flow of communication between any two cells; the total maximum flow (aggregated across all pairs of third parties) mediated by v.
  #centr$info <- try(sna::infocent(net)) # actors with higher information centrality are predicted to have greater control over the flow of information within a network; highly information-central individuals tend to have a large number of short paths to many others within the social structure.
  centr$flowbet <- tryCatch({
    sna::flowbet(net)
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  centr$info <- tryCatch({
    sna::infocent(net, diag = T, rescale = T, cmode = "lower")
    # sna::infocent(net, diag = T, rescale = T, cmode = "weak")
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  return(centr)
}

netAnalysis_signalingRole_heatmap <- function (object=object.list[[2]], signaling = NULL, pattern = "all", slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", 
    title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
    cluster.rows = FALSE, cluster.cols = FALSE) 
{
    pattern <- match.arg(pattern)
    if (length(slot(object, slot.name)$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
        outgoing[, i] <- centr[[i]]$outdeg
        incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
        mat <- t(outgoing)
        legend.name <- "Outgoing"
    }
    else if (pattern == "incoming") {
        mat <- t(incoming)
        legend.name <- "Incoming"
    }
    else if (pattern == "all") {
        mat <- t(outgoing + incoming)
        legend.name <- "Overall"
    }
    if (is.null(title)) {
        title <- paste0(legend.name, " signaling patterns")
    }
    else {
        title <- paste0(paste0(legend.name, " signaling patterns"), 
            " - ", title)
    }
    if (!is.null(signaling)) {
        mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
        mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
        idx <- match(rownames(mat1), signaling)
        mat[idx[!is.na(idx)], ] <- mat1
        dimnames(mat) <- list(signaling, colnames(mat1))
    }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    mat[mat == 0] <- NA
    if (is.null(color.use)) {
        color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
        name = color.heatmap))))(100)
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    color.use <- color.use[colnames(mat)]
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
        which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
        simple_anno_size = grid::unit(0.2, "cm"))
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
        border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
        show_annotation_name = FALSE)
    pSum <- rowSums(mat.ori)
    pSum.original <- pSum
    pSum <- -1/log(pSum)
    pSum[is.na(pSum)] <- 0
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    if (length(idx1) > 0) {
        values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
            length.out = length(idx1))
        position <- sort(pSum.original[idx1], index.return = TRUE)$ix
        pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    # ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
    #     show_annotation_name = FALSE)
    if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
        legend.break <- max(mat, na.rm = T)
    }
    else {
        legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
            round(max(mat, na.rm = T), digits = 1))
    }
    library(ComplexHeatmap)
    library(circlize)  

    col_fun = colorRamp2(c(1, 0), c("#1F77B4", "white"))
 
    ht1 = Heatmap(mat, 
                  col = col_fun,
                  na_col = "white", 
                  name = "Relative strength", 
                  bottom_annotation = col_annotation,
                  top_annotation = ha2,
                  # right_annotation = ha1,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  row_names_side = "left",
                  row_names_rot = 0, 
                  row_names_gp = gpar(fontsize = font.size), 
                  column_names_gp = gpar(fontsize = font.size), 
                  width = unit(width, "cm"), 
                  height = unit(height, "cm"), 
                  column_title = title, 
                  row_order = sort(rownames(mat)),
                  column_title_gp = gpar(fontsize = font.size.title), 
                  column_names_rot = 45,  # x 轴标签 45 度
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 10,
                                                              fontface = "plain"),
                                              title_position = "leftcenter-rot",
                                              border = NA, 
                                              at = c(0, 0.5, 1),  # 图例刻度对应 0 到 1
                                              legend_height = unit(20, "mm"), 
                                              labels_gp = gpar(fontsize = 8), 
                                              grid_width = unit(2, "mm")))
    return(ht1)
}