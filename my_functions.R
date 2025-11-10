query_ref_modulesgene <- function(query.obj, ref.obj, 
                                  ref.celltype = "maintype", 
                                  same.species = TRUE, 
                                  ngenes = 100,## gene number used for the signature calculation
                                  with_tie = TRUE, 
                                  rm.ribos.mt = TRUE
                                 ){ 
    library(R.utils)
    library(dplyr)
    library(Seurat)
    
    # ref.obj <- ref.obj %>% 
    #     NormalizeData() %>% 
    #     FindVariableFeatures() %>% 
    #     ScaleData() %>% 
    #     RunPCA()

    ###############################
    #return the top markers identified in the reference data -> for signiture calculation

    if(rm.ribos.mt){ 
        if(same.species){ 
            genes <- intersect(rownames(query.obj), rownames(ref.obj))
            gene.rm <- c(grep(genes, pattern = "^MT-", ignore.case = TRUE, value = TRUE), 
                        grep(genes, pattern = "^RPL", ignore.case = TRUE, value = TRUE), 
                        grep(genes, pattern = "^RPS", ignore.case = TRUE, value = TRUE))
            genes <- genes[!genes %in% gene.rm]
            allmarkers <- FindAllMarkers(ref.obj, 
                                         group.by = ref.celltype, 
                                         only.pos = TRUE, 
                                         min.pct = 0.25)
            markers <- allmarkers %>% 
                filter(p_val_adj < 0.05 & gene %in% genes) %>% 
                group_by(cluster) %>% 
                slice_min(order_by = p_val_adj, n = ngenes, with_ties = with_tie)
            return(markers)

    }else{
            genes <- intersect(rownames(query.obj) %>% tolower()  %>% capitalize(), 
                               rownames(ref.obj) %>% tolower()  %>% capitalize())
            
            gene.rm <- c(grep(genes, pattern = "^MT-", ignore.case = TRUE, value = TRUE), 
                        grep(genes, pattern = "^RPL", ignore.case = TRUE, value = TRUE), 
                        grep(genes, pattern = "^RPS", ignore.case = TRUE, value = TRUE))
            genes <- genes[!genes %in% gene.rm]
            
            allmarkers <- FindAllMarkers(ref.obj, 
                                         group.by = ref.celltype, 
                                         only.pos = TRUE, 
                                         min.pct = 0.25)
            allmarkers$gene <- allmarkers$gene %>% tolower()  %>% capitalize()
            
            markers <- allmarkers %>% 
                filter(p_val_adj < 0.05 & gene %in% genes)%>%
                group_by(cluster) %>% 
                slice_min(order_by = p_val_adj, n = ngenes, with_ties = with_tie)
            return(markers) 
        }
        
    }else{ 

        if(same.species){ 
            genes <- intersect(rownames(query.obj), rownames(ref.obj))
            allmarkers <- FindAllMarkers(ref.obj, 
                                         group.by = ref.celltype, 
                                         only.pos = TRUE, 
                                         min.pct = 0.25)
            markers <- allmarkers %>% 
                filter(p_val_adj < 0.05 & gene %in% genes) %>% 
                group_by(cluster) %>% 
                slice_min(order_by = p_val_adj, n = ngenes, with_ties = with_tie)
            return(markers)

    }else{
            genes <- intersect(rownames(query.obj) %>% tolower()  %>% capitalize(), 
                               rownames(ref.obj) %>% tolower()  %>% capitalize())
            
            allmarkers <- FindAllMarkers(ref.obj, 
                                         group.by = ref.celltype, 
                                         only.pos = TRUE, 
                                         min.pct = 0.25)
            allmarkers$gene <- allmarkers$gene %>% tolower()  %>% capitalize()
            
            markers <- allmarkers %>% 
                filter(p_val_adj < 0.05 & gene %in% genes)%>%
                group_by(cluster) %>% 
                slice_min(order_by = p_val_adj, n = ngenes, with_ties = with_tie)
            return(markers) 
        }
    }
}



query_ref_score <- function(
    query.obj,
    geneset = geneset,
    same.species = TRUE, 
    assay = "RNA",
    with_tie = TRUE){ 
    
   for(clu in unique(geneset$cluster)){ 
    marker <- subset(geneset, cluster == clu) %>% pull(gene)
    query.obj <- AddModuleScore(query.obj, features = list(clu = marker), name = clu)
}
    return(query.obj)
}
newseurat <- function(obj){  
    temp <- CreateSeuratObject(counts = GetAssayData(obj,assay = "RNA", layer = "counts"), 
                         meta.data = obj@meta.data)
    return(temp)
}

my_cortest <- function(mat = NULL,
                       gene = NULL, 
                       method = "spearman", # "kendall","pearson" 
                       p.adjust.method = "BH",
                       n_cores = 5) {

    options(future.globals.maxSize = 6 * 1024^3)
    
    ###### 1. filter data
    if (is.null(mat) || is.null(gene)) stop("mat and gene required")
    if (!gene %in% rownames(mat)) stop("gene not in rownames(mat)")
    
    cells_keep <- colnames(mat)[colSums(mat > 0) > 0]
    genes_keep <- rownames(mat)[rowSums(mat > 0) > 0]
    
    if (!length(cells_keep)) return(data.frame())
    mat <- mat[genes_keep, cells_keep, drop = FALSE]
    
    target <- mat[gene, ]
    other_genes <- setdiff(rownames(mat), gene)
    n <- length(other_genes)
    
    ###### 2. set cores
    if (is.null(n_cores)) n_cores <- 5
    future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(future::sequential), add = TRUE)

    ###### 3. calculated
    compute_one <- function(g) {
        ct <- cor.test(target, mat[g, ], method = method)
        data.frame(gene1 = gene, 
                   gene2 = g, 
                   corr = ct$estimate,
                   pval = ct$p.value, 
                   method = method,
                   stringsAsFactors = FALSE)
    }
    
    result_list <- future.apply::future_lapply(other_genes, compute_one)
    result <- do.call(rbind, result_list)
    result$pval_adj <- p.adjust(result$pval, 
                              method = p.adjust.method)
    
    result <- result[order(result$pval_adj), ]
    return(result)
}

my_volcan.plt <- function(meta, 
                   cutoff_lfc, 
                   cutoff_qval, 
                   x, y, label, 
                    max.overlaps = 20,
                    title = "DEGs",
                    label.size = 2.5, 
                    color
                   ){ 
    library(ggrepel)
    library(ggplot2)
    
    p <- ggplot(meta, aes(x = .data[[x]], y = -log10(.data[[y]]))) + 
        geom_hline(yintercept = -log10(cutoff_qval), linetype = "dashed", color = "#999999") +
        geom_vline(xintercept = c(-cutoff_lfc,cutoff_lfc), linetype = "dashed", color = "#999999") + 
        geom_point(size = 1, color = ifelse(meta[,color] < -cutoff_lfc, "red", ifelse(meta[,color]  > cutoff_lfc, "blue", "gray"))) + 
        geom_text_repel(data = meta, 
                        label = meta[,label], 
                        max.overlaps = max.overlaps, 
                        size = label.size, 
                        box.padding = 0, 
                        nudge_x = 0.05, 
                        nudge_y = 0.1) +
        theme_bw(base_size = 12)+
        ggsci::scale_color_jama() +
        theme(panel.grid = element_blank(),
              legend.position = 'right') + 
        labs(x = "Log2FC", y = "-Log10(q-value)", title = title)                                              
        return(p)
}




my_bar_percent.plt <- function(
    meta,
    var.x1, 
    var.x2, 
    levels.x2,
    var.fill,
    titlename = "percentage",
    label.size = 2.5,
    colors){
    
    library(ggplot2)
    library(dplyr)
    plotdata <- prop.table(table(meta[,c(var.x1, var.x2)]), margin = 1) %>% as.data.frame()
    plotdata[,var.x2]<- factor(plotdata[,var.x2], levels = levels.x2)
    plotdata$label <- paste0(round(plotdata$Freq, 3) * 100, "%")
    
    p <- ggplot(plotdata, aes(fill = .data[[var.fill]], y = Freq, x = .data[[var.x1]])) + 
        geom_bar(position = "fill", stat="identity", alpha = 0.8) + 
        theme_classic() + 
        geom_text(aes(label=label), position = position_stack(vjust = 0.5), size = label.size) + 
        scale_fill_manual(values = colors) + 
        theme(axis.text = element_text(size = 12, color = "black"), 
              axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
              axis.title = element_text(size = 12, face = "bold", color = "black") ) & 
        ggtitle(titlename) & labs(x = "", y = "Percentage")
    return(p)
}


my_density.plt <- function(meta, 
                           x.var = "pseudotime", 
                           fill.var = "time", 
                           color.var = "time", 
                           colors = colors, 
                           alpha = 0.5, 
                           legend.pos = "right", 
                           facet.var = "time",
                           facet = TRUE, 
                            scale =TRUE
                          ){ 
   library(ggplot2)
   if(facet){ 
       if(scale){
           p <- ggplot(meta, aes(x = .data[[x.var]], fill = .data[[fill.var]], color = .data[[color.var]])) + 
        geom_density(alpha = alpha) + 
        scale_fill_manual(values = colors) + 
        scale_color_manual(values = colors) + 
        theme_bw() + 
        facet_grid(.data[[facet.var]]~.) +
        theme(axis.text = element_text(size = 8),  
              axis.title = element_text(size = 12), 
              strip.text = element_text(size = 8), 
             legend.position = legend.pos) 
        return(p) 
       }else{

        p <- ggplot(meta, aes(x = .data[[x.var]], fill = .data[[fill.var]], color = .data[[color.var]])) + 
            geom_density(alpha = alpha) + 
            scale_fill_manual(values = colors) + 
            scale_color_manual(values = colors) + 
            theme_bw() +
            facet_grid(.data[[facet.var]]~., scales = "free") +
            theme(axis.text = element_text(size = 8),  
                  axis.title = element_text(size = 12), 
                  strip.text = element_text(size = 8), 
                 legend.position = legend.pos) 
        return(p) 
       }
       
   
   }else{ 
          p <- ggplot(meta, aes(x = .data[[x.var]], fill = .data[[fill.var]], color = .data[[color.var]])) + 
        geom_density(alpha = alpha) + 
        scale_fill_manual(values = colors) + 
        scale_color_manual(values = colors) + 
        theme_bw() + 
        theme(axis.text = element_text(size = 8),  
              axis.title = element_text(size = 12), 
              strip.text = element_text(size = 8), 
             legend.position = legend.pos) 
        return(p) 
   }
}

my_viol_compare.plt <- function(
    meta = meta,
    x = x, 
    y = y, 
    color = color, 
    title = "cytoTRACE", 
    ylab = "Potency score", 
    facet = TRUE, 
    facet.var = NULL, 
    colors = NULL,scale = TRUE,
    label.size = 6,
    label_y_positions = c(0.2,0.3, 0.4),
    label_x_positions = c(1.5, 2.5),
    
    comparison = NULL, legend.pos = "none"
){ 
    library(ggplot2)
    library(ggpubr)

    # if(length(comparison)>0){ 
    #     label_y_positions = c()
    #     for(i in 1:length(comparison)){ 
    #     label_y_positions = c(label_y_positions, max(meta[,y])+0.1)
    #     }
  
    # }

    if(facet){ 

        if(scale){ 
        p <- ggplot(meta, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) + 
            geom_violin() + 
            stat_summary(fun = "mean", geom = "crossbar", width = 0.3, size = 0.5, aes(color = .data[[color]]))+ 
            facet_grid(.~.data[[facet.var]], scales = "free") + 
            scale_color_manual(values = colors) +
            stat_compare_means(method = "wilcox.test", size = label.size, comparison = comparison, label.y =label_y_positions, label.x = label_x_positions) +
            theme_bw(base_size = 12) + 
            theme(plot.title = element_text(hjust = 0.5), 
        legend.position = legend.pos, 
        axis.text.x = element_text(angle = 45, 
                                         size = 12, 
                                         hjust = 1)) + 
            ggtitle(title) + labs(x = NULL, y = ylab)
        return(p)  
        
        }else{
    p <- ggplot(meta, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) + 
        geom_violin() + 
        stat_summary(fun = "mean", geom = "crossbar", width = 0.3, size = 0.5, aes(color = .data[[color]]))+ 
        facet_grid(.~.data[[facet.var]]) + 
        scale_color_manual(values = colors) +
        stat_compare_means(method = "wilcox.test", size = label.size, comparison = comparison, label.y =label_y_positions, label.x = label_x_positions) +
        theme_bw(base_size = 12) + 
        theme(plot.title = element_text(hjust = 0.5), 
        legend.position = legend.pos, 
        axis.text.x = element_text(angle = 45, 
                                         size = 12, 
                                         hjust = 1)) + 
        ggtitle(title) + labs(x = NULL, y = ylab)
     if(is.null(colors)){ 
        return(p)  
    }else{
        p <- p + scale_color_manual(values = colors) 
        return(p)  
    } 
        }
    

        
    }else{ 
        p <- ggplot(meta, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) + 
        geom_violin() + 
        stat_summary(fun = "mean", geom = "crossbar", width = 0.3, size = 0.5, aes(color = .data[[color]]))+ 
        # facet_grid(.~.data[[facet.var]]) + 
        # scale_color_manual(values = colors) +
        stat_compare_means(method = "wilcox.test", size = 4, comparison = comparison, label.y =label_y_positions, label.x = label_x_positions) +
        theme_bw(base_size = 12) + 
        theme(plot.title = element_text(hjust = 0.5), 
              legend.position = legend.pos, 
              axis.text.x = element_text(angle = 45, 
                                         size = 12, 
                                         hjust = 1)) + 
        ggtitle(title) + labs(x = NULL, y = ylab)

        if(is.null(colors)){ 
            return(p)  
        }else{
            p <- p + scale_color_manual(values = colors) 
            return(p)  
        }
        
    }
}
ensemble2symbol <- function(species, mat){ 
    library(dplyr)
    genes <- rownames(mat)
    if(species == "mouse"){ 

        library(org.Mm.eg.db)
        k <- keys(org.Mm.eg.db,keytype = "ENSEMBL")
        info <- AnnotationDbi::select(org.Mm.eg.db, 
                                       keys = k,
                                       columns = c("ENTREZID", "SYMBOL"), 
                                       keytype = "ENSEMBL")
            
        info.sub <- info[match(genes,info[,"ENSEMBL"]),] %>% na.omit()
        mat.sub <- mat[info.sub$ENSEMBL,]

        idx <- c(1:nrow(info.sub))[!duplicated(info.sub$SYMBOL)]
        mat.sub <- mat.sub[idx,]
        rownames(mat.sub) <- info.sub$SYMBOL[idx]
            
    }else if(species == "human"){ 
        library(org.Hs.eg.db, "/home/lfliuhku/conda/pkgs/jupyterlab/lib/R//library")
        info <- AnnotationDbi::select(org.Hs.eg.db, 
                                       keys = k,
                                       columns = c("ENTREZID", "SYMBOL"), 
                                       keytype = "ENSEMBL")
        
        info.sub <- info[match(genes,info[,"ENSEMBL"]),] %>% na.omit()
        mat.sub <- mat[info.sub$ENSEMBL,]

        idx <- c(1:nrow(info.sub))[!duplicated(info.sub$SYMBOL)]
        mat.sub <- mat.sub[idx,]
        rownames(mat.sub) <- info.sub$SYMBOL[idx]
    }
    return(mat.sub)
}

loom2seurat <- function(path_Data, subname, return.obj = TRUE, save.obj = TRUE){ 
    library(loomR)
    library(SeuratDisk)
    library(SeuratData)
    library(Seurat)
    data <- Connect(path_data, mode = "r") 
    obj  <- as.Seurat(data)

   if(save.obj){saveRDS(obj, file = paste0(subname, ".rds"))}
   if(return.obj){return(obj)}
    
}


loom2h5ad <- function(path_Data, subname){ 
    library(loomR)
    library(SeuratDisk)
    library(SeuratData)
    library(Seurat)
    data <- Connect(path_data, mode = "r") 
    obj  <- as.Seurat(data)
    
    obj[["RNA3"]] <- as(object = obj[["RNA"]], Class = "Assay")
    DefaultAssay(obj) <- "RNA3"
    obj[["RNA"]] <- NULL
    obj <- RenameAssays(object = obj, RNA3 = 'RNA')
    temp.file <- paste0(subname, ".h5Seurat")
    SaveH5Seurat(obj, file = temp.file, overwrite = TRUE)
    Convert(paste0(subname, ".h5Seurat"), 
            dest = "h5ad", 
            overwrite = TRUE) 

    if (file.exists(temp.file)) {
      file.remove(temp.file)
      cat(paste0("deleted", temp.file))
    } else {
      cat("No file found")
    }
}


seurat2h5ad <- function(obj, subname) {
    library(SeuratData)
    library(SeuratDisk)
    library(Seurat)
    # temp <- obj
    # obj <- CreateSeuratObject(counts = GetAssayData(temp, assay = "RNA", layer = "counts"), 
    #                           meta.data = temp@meta.data, min.cells = 10)
    
    obj[["RNA3"]] <- as(object = obj[["RNA"]], Class = "Assay")
    DefaultAssay(obj) <- "RNA3"
    obj[["RNA"]] <- NULL
    obj <- RenameAssays(object = obj, RNA3 = "RNA")
    temp.file <- paste0(subname, ".h5Seurat")
    SaveH5Seurat(obj, file = temp.file, overwrite = TRUE)
    Convert(paste0(subname, ".h5Seurat"), dest = "h5ad", overwrite = TRUE, 
        features = rownames(obj))
    if (file.exists(temp.file)) {
        file.remove(temp.file)
        cat(paste0("deleted", temp.file))
    }
    else {
        cat("No file found")
    }
}


plot_gene_trend <- function(plot.data, pseudotime = "pseudotime", group_by = "group_by", genes = NULL, title = NULL){ 
    library(ggplot2)
    
    plot.data$pseudotime <- (plot.data$pseudotime - min(plot.data$pseudotime))/(max(plot.data$pseudotime)- min(plot.data$pseudotime))
        
    if(length(genes) == 1){ 
        ggplot(plot.data, aes(x = .data[[pseudotime]], y = .data[[genes]], color = .data[[group_by]])) + 
        geom_point(size = 0.5, alpha = 0.2) + 
        geom_smooth(method = "loess", se = TRUE) + 
        scale_color_manual(values = c("Vcl cKO" = "#f7819f", "Control" = "#29b8e1")) + 
        theme_bw() + 
        theme(axis.text = element_text(size = 14), 
              axis.title = element_text(size = 16), 
              plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
              legend.text = element_text(size = 12), 
              legend.title = element_text(size = 14), 
              strip.text = element_text(size = 14, face = "bold")
    ) + 
        ylab(paste0("Scaled Expression")) & 
        scale_x_continuous(limits = c(0, 1), breaks = c(0.3, 0.6, 0.9)) & ggtitle(title)
    }else{ 

        plot.data2 <- data.frame("pseudotime" = plot.data[,pseudotime], 
                                 "group_by" = plot.data[,group_by],
                                 "avg_expr" = rowMeans(plot.data[,genes], na.rm = TRUE))
                                 
        rownames(plot.data2) <- rownames(plot.data)

        ggplot(plot.data2, aes(x = .data[[pseudotime]], y = .data[["avg_expr"]], color = .data[[group_by]])) + 
        geom_point(size = 0.5, alpha = 0.2) + 
        geom_smooth(method = "loess", se = TRUE) + 
        scale_color_manual(values = c("Vcl cKO" = "#f7819f", "Control" = "#29b8e1")) + 
        theme_bw() + 
        theme(axis.text = element_text(size = 14), 
              axis.title = element_text(size = 16), 
              plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
              legend.text = element_text(size = 12), 
              legend.title = element_text(size = 14), 
              strip.text = element_text(size = 14, face = "bold")
    ) + 
        ylab(paste0("Average Scaled Expression")) & 
        scale_x_continuous(limits = c(0, 1), breaks = c(0.3, 0.6, 0.9))  & ggtitle(title)
    
    } 

}



enrich.go <- function(genes, species, ont, name, convertID = TRUE){ 
    if(species == "mouse"){ 
        library(clusterProfiler)
        library(org.Mm.eg.db)
        library(ggplot2)
        # 
        if(convertID){ 
         gene.enriched <- AnnotationDbi::select(org.Mm.eg.db, 
                keys = genes, 
                columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                keytype = "SYMBOL")$ENTREZID 
         
        gene.enriched <- na.omit(gene.enriched)
             go <- enrichGO(gene = gene.enriched, # EntrezID
                       OrgDb = org.Mm.eg.db,  
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       ont = ont)
            go <- setReadable(go, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
        }else{
            gene.enriched = genes
            print(paste0(name, ":",length(gene.enriched)))
              go <- enrichGO(gene = gene.enriched, # EntrezID
                       OrgDb = org.Mm.eg.db,  
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       ont = ont, keyType = "SYMBOL")
       }
        
        
    
        ############## GO
       
        if(length(go$ID) > 0){ 
            data.go <- data.frame(
                "ID" = go$ID,
                "Description" = go$Description,
                "p.adjust" = go$p.adjust,
                "GeneRatio" = go$GeneRatio, 
                "Count" = go$Count,
                "group" = name, 
                "geneID" = go$geneID)
        }
        data.go <- dplyr::arrange(data.go,p.adjust)
        return(data.go)
        
    }else if(species == "human"){
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(ggplot2)
        library(dplyr)
        if(convertID){ 
         gene.enriched <- AnnotationDbi::select(org.Hs.eg.db, 
                keys = genes, 
                columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                keytype = "SYMBOL")$ENTREZID

         gene.enriched <- na.omit(gene.enriched)
        }else{gene.enriched = genes}
        
        print(paste0(name, ":",length(gene.enriched)))
    
        ############## GO
        go <- enrichGO(gene = gene.enriched, # EntrezID
                       OrgDb = org.Hs.eg.db,  
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       ont = "BP")
        go <- setReadable(go, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
        if(length(go$ID) > 0){ 
            data.go <- data.frame(
                "ID" = go$ID,
                "Description" = go$Description,
                "p.adjust" = go$p.adjust,
                "GeneRatio" = go$GeneRatio, 
                "Count" = go$Count,
                "group" = name, 
                "geneID" = go$geneID)
        }
        data.go <- dplyr::arrange(data.go, p.adjust)
        return(data.go)
    
    }
}




enrich.kegg <- function(genes, species, ont, name){ 
    if(species == "mouse"){ 
        library(clusterProfiler)
        library(org.Mm.eg.db)
        library(ggplot2)
        library(dplyr)
         gene.enriched <- AnnotationDbi::select(org.Mm.eg.db, 
                keys = genes, 
                columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                keytype = "SYMBOL")$ENTREZID ; 
        gene.enriched <- na.omit(gene.enriched)
        
        print(paste0(name, ":",length(gene.enriched)))
    
        ############## kegg
        kegg <- enrichKEGG(gene = gene.enriched, # EntrezID
                       organism = "mmu",  
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)
        kegg <- setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
        if(length(kegg$ID) > 0){ 
            data.kegg <- data.frame(
                "ID" = kegg$ID,
                "Description" =  sapply(kegg$Description, function(x){gsub(x, pattern = " - Mus musculus \\(house mouse\\)", replacement = "")}),
                "p.adjust" = kegg$p.adjust,
                "GeneRatio" = kegg$GeneRatio, 
                "Count" = kegg$Count,
                "group" = name, 
                "geneID" = kegg$geneID)
        }
        data.kegg <- dplyr::arrange(data.kegg, p.adjust)
        return(data.kegg)
        
    }else if(species == "human"){
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(ggplot2)
        library(dplyr)
        
        gene.enriched <- AnnotationDbi::select(org.Hs.eg.db, 
        keys = genes, 
        columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
        keytype = "SYMBOL")$ENTREZID; 
        gene.enriched <- na.omit(gene.enriched)
        print(paste0(name, ":",length(gene.enriched)))
    
        ############## kegg
        kegg <- enrichKEGG(gene = gene.enriched, # EntrezID
                      organism = "mmu",  
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)
        kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
        if(length(kegg$ID) > 0){ 
            data.kegg <- data.frame(
                "ID" = kegg$ID,
                "Description" =  sapply(kegg$Description, function(x){gsub(x, pattern = " - Mus musculus \\(house mouse\\)", replacement = "")}),
                "p.adjust" = kegg$p.adjust,
                "GeneRatio" = kegg$GeneRatio, 
                "Count" = kegg$Count,
                "group" = name, 
                "geneID" = kegg$geneID)
        }
        
        data.kegg <- dplyr::arrange(data.kegg, p.adjust)
        
        return(data.kegg)
    
    }
}

dotplt.enriched <- function(data = data, x = "group", y = "Description", size = "Count", color = "p.adjust", title){ 
    library(ggplot2)
    p <- ggplot(data = data, aes(x = .data[[x]], y = .data[[y]])) + 
  geom_point(aes(size = .data[[size]], color = .data[[color]])) +
  scale_color_gradient(low='#bc4749', high='#a7c957')+  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black")) + ggtitle(title) + xlab(NULL) + ylab("Pathways") + scale_size(range = c(4,7)) + 
    theme(axis.text = element_text(size = 14), 
          axis.text.x = element_text(angle = 30, hjust = 1), 
          plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 14)
          ) + 
#        scale_color_gradient(low = "#DF3A01", high = "#F6D8CE", breaks = c(0.01, 0.04),limits = c(0.00, 0.05));
        scale_color_gradient(low = "#DF3A01", high = "#F6D8CE");
    
    return(p)
}
