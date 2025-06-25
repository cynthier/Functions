
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
    p <- ggplot(data = data, aes(x = .data[[x]], y = .data[[y]])) + 
  geom_point(aes(size = .data[[size]], color = .data[[color]])) +
  scale_color_gradient(low='#bc4749', high='#a7c957')+  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black")) + ggtitle(title) + xlab(NULL) + ylab("Pathways") + scale_size(range = c(3,8)) + 
    theme(axis.text = element_text(size = 14), 
          axis.text.x = element_text(angle = 30, hjust = 1), 
          plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 14)
          ) 
    + 
        scale_color_gradient(low = "#DF3A01", high = "#F6D8CE", breaks = c(min(data$qvalue), 0.05), limits = c(0,0.05));
    return(p)
}
