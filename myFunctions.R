dpt.exp.cell.pseBin <- function(
    obj, 
    var.pseudotime,
    var.y = "bin",
    meta, #### meta data contained "Pseudotime" and cell barcodes
    plot.genes, 
    nbin = 20,
    rev.direction = FALSE, 
    order.gene = TRUE,      # Added: missing parameter
    col = "#a2042c",        # Added: missing color variable
    label.gene = NULL, 
    return.gene = FALSE, 
    ylab = "lineage",
    subtitle = "--->"
){

    library(dplyr)
    library(ggplot2)

    ##### 1. get the expressed cells
    mat <- Seurat::GetAssayData(obj, slot = "data", assay = "RNA")
    genes <- intersect(plot.genes, rownames(mat))

    if(length(intersect(meta[,"cell"], colnames(mat))) == 0){
        print("no cells found")
        return(NULL)
    }else{
        result <- lapply(genes, function(gene){
            # Get cells with expression > 0
            cell.exp <- colnames(mat)[which(mat[gene,] > 0)]
            # Match with meta cells
            valid.cells <- intersect(cell.exp, meta$cell)
            
            temp <- data.frame(
                "cell" = valid.cells,
                "Pseudotime" = meta[match(valid.cells, meta$cell), var.pseudotime],
                "exp" = as.numeric(mat[gene, valid.cells])
            )
            temp$gene <- gene
            return(temp)
        }) %>% do.call(rbind, .)
    }

    ##### 2. divided into bins accroding to the pseudotime
    bindata <- meta %>% 
        arrange(.data[[var.pseudotime]]) %>% 
        mutate(bin = ntile(.data[[var.pseudotime]], nbin))
    
    # Use the 'bin' column from bindata to assign bins to the result
    result$bin <- bindata$bin[match(result$cell, bindata$cell)]

    ##### 3. order the genes
    if(order.gene){
        gene_order <- result %>%
            group_by(gene) %>%
            summarise(mean_pseudo = mean(.data[[var.pseudotime]], na.rm = TRUE)) %>%
            arrange(desc(mean_pseudo)) %>%  
            pull(gene) 
        
        result$gene <- factor(result$gene, levels = gene_order)
    } else {
        gene_order <- genes
        result$gene <- factor(result$gene, levels = genes)
    }
 
    ##### 4. plot the data
    # Prepare the y-axis variable
    if(rev.direction){
        result[[var.y]] <- -result[[var.y]]
    }

    p1 <- ggplot(result, aes(x = gene, y = .data[[var.y]])) + 
        geom_jitter(size = 0.4, color = col, alpha = 0.5) +
        labs(x = "", y = ylab, title = "", subtitle = subtitle) +
        theme_bw() + coord_flip() +
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.title = element_text(size = 13), 
            legend.text = element_text(size = 10),
            panel.grid.major.y = element_line(color = "grey90"),
            panel.grid.minor.y = element_line(color = "grey90"),

            axis.text = element_text(angle = 0, size = 12, color = "black", face = "plain"),
            axis.text.y = element_text(
                angle = 0, size = 10, face = "italic", 
                color = ifelse(levels(result$gene) %in% label.gene, "#a2042c", "black")
            ),
            axis.title = element_text(size = 12, face = "plain")
        ) 

    if(return.gene){
        return(gene_order)
    }else{
        return(p1)
    }
}


getgradients <- function(n){
	gradients <- list(
		c(colorRampPalette(c("#5b51a3","#79c9a4","#f2faac","#fdb465","#a4104d"))(90))
	)
	if(n < length(gradients) + 1){
		return(gradients[[n]])
	}else{
		print(paste0("only",length(gradients), "gradient available"))
	}
}



aucell_score <- function(obj, 
                         geneset, 
                         scorename, 
                         comparison, 
                         return.plot = FALSE, 
                         x.var,
                         color.var,
                         fill.var, 
                         facet = F,
                         facet.var,
                         colors,
                         sub.title = NULL,
						 ylabel
                        ){ 
    library(AUCell)
    library(ggplot2)
    library(ggpubr)
    library(gridExtra)
	count.mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
	meta <- obj@meta.data

	orig.length.geneset <- length(geneset)
	geneset <- geneset[geneset %in% rownames(count.mat)]
	new.length.geneset <- length(geneset)

	cells.rankings <- AUCell_buildRankings(count.mat, plotStats = FALSE)
	cells_AUC <- AUCell_calcAUC(list("geneset" = geneset), cells.rankings)
	auc.mat <- getAUC(cells_AUC)
	meta[scorename] <- auc.mat[1,rownames(meta), drop = T]
    
    if(return.plot == FALSE){ 
        return(meta)
    }else{ 
        num <- length(intersect(geneset, rownames(count.mat)))
    	plot.data <- meta[meta[,x.var] %in% unlist(comparison),]
        p <- my_viol_compare.plt(
        	  meta = plot.data,
        	  x = x.var,
        	  y = scorename,
        	  color = color.var,
        	  fill = fill.var,
        	  color.fill = colors,
        	  color.color = colors,
        	  facet = facet,
        	  facet.var = facet.var,
        	  show.mean = TRUE,
        	  comparison = comparison,
        	  label.size = 4,
        	  mean.value.size = 4, 
        	  label_y_positions = max(plot.data[scorename]),
        	  title = paste0(scorename, "(",  new.length.geneset, "/",orig.length.geneset, " genes)")) + 
            ylim(min(plot.data[scorename])*1.2, max(plot.data[scorename])*1.5) +
            theme(axis.text.x = element_text(angle = 60)) + labs(subtitle = sub.title,y = ylabel)
      return(p)
    }	
}




monocle_deg_func <- function(obj, 
							cells.used, 
							genes.used = NULL, 
							number.normalization = FALSE, 
							group.var = "group", 
							assay.name = "RNA", 
							name.ctrl = "Control", 
							name.mutant = "Mutant"){

	library(monocle)
	library(Seurat)
	library(dplyr)
	library(foreach)
	library(doParallel)
	library(monocle)

	num_cores <- 5  
	registerDoParallel(cores = num_cores)


	  ##### subset the object
	  obj <- obj[,cells.used]
	  metadata <- obj@meta.data
	  if(is.null(genes.used)){
	  matrix <- GetAssayData(obj, layer = "counts", assay = assay.name)[, cells.used]
	  }else{
	  matrix <- GetAssayData(obj, layer = "counts", assay = assay.name)[genes.used, cells.used]
	  }
	  
	  ##### gene annotation
	  gene_anno <- data.frame("gene_short_name" = rownames(matrix))
	  rownames(gene_anno) <- gene_anno$gene_short_name
	  
	  ##### subset the cells 
	  cells_ko.all <- rownames(subset(metadata, group == name.mutant))
	  cells_con.all <- rownames(subset(metadata, group == name.ctrl))
	  
	  if(number.normalization){
		cell.ko.selected <- sample(cells_ko.all, size = min(length(cells_ko.all), length(cells_con.all)))
		cell.con.selected <- sample(cells_con.all, size = min(length(cells_ko.all), length(cells_con.all)))
	  }else{
		cell.ko.selected <- cells_ko.all
		cell.con.selected <- cells_con.all
	  }
	  
	  matrix <- matrix[,c(cell.ko.selected, cell.con.selected)]
	  metadata <- metadata[c(cell.ko.selected, cell.con.selected),]
	
	  ##### create the obj
	  pd <- new("AnnotatedDataFrame", data = metadata)
	  fd <- new("AnnotatedDataFrame", data = gene_anno)
	  
	  cds <- newCellDataSet(
	    matrix,
	    phenoData = pd,
	    featureData = fd,
	    expressionFamily = negbinomial.size()
	  )
	  
	  ##### estimate the size factors and dispersions
	  cds <- estimateSizeFactors(cds)
	  cds <- estimateDispersions(cds)
	  
	  ##### differential gene expression test
	  res <- differentialGeneTest(cds, fullModelFormulaStr = paste0("~", group.var))
	  res <- res %>% mutate(qval = p.adjust(pval, method = "fdr"))
	  
	  ##### add log2 fold change
	  a <- 1 + apply(matrix[rownames(res), cell.ko.selected], 1, mean)
	  b <- 1 + apply(matrix[rownames(res), cell.con.selected], 1, mean)
	  res$log2FC_mean <- log2(a / b)[res$gene]
	  
	  ##### add percentage expressed
	  res$bcg.pct <- (rowSums(matrix[, cell.con.selected] > 0) / length(cell.con.selected))[res$gene]
	  res$obj.pct <- (rowSums(matrix[, cell.ko.selected] > 0) / length(cell.ko.selected))[res$gene] 
	  
	  res$number_control_cells <- length(cell.con.selected)
	  res$number_mutant_cells <- length(cell.ko.selected)
	  
	  if(number.normalization){
		return(list(res, cell.ko.selected, cell.con.selected))
	  }else{
		return(res)
	  }

}



my_positive_proportion <- function(obj, gene, xvar, xvar_order = NULL, substring = NULL, assay = "RNA"){ 
    library(dplyr)
    library(Seurat)
    library(ggplot2)
    
    plotdata <- obj@meta.data
    values <- GetAssayData(obj, layer = "counts",assay = assay)[gene,rownames(plotdata)]
    if(all(values %% 1 == 0)){ 
        plotdata[,gene] <- values
        plotdata[,"exp"] <- ifelse(plotdata[,gene] > 0, "detected", "undetected") 
        summ <- table(plotdata[,c(xvar, "exp")]) %>% 
            prop.table(margin = 1) %>% 
            round(digits = 3) %>% 
            as.matrix() %>% 
            as.data.frame()
        
        if(length(xvar_order) > 0){ 
            summ$xvar <- factor(summ$xvar, levels = xvar_order)
        }

        summ$exp <- factor(summ$exp, levels = c("undetected", "detected"))
       p <-  ggplot(summ, aes(x = .data[[xvar]],
                         y = .data[["Freq"]], 
                         fill =  .data[["exp"]])) + 
            geom_bar(position = "stack", stat = "identity", width = 0.6) + 
            theme_bw(base_size = 12) + 
            theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
                  axis.text.y = element_text(color = "black"),
                  plot.title = element_text(hjust = 0.5),
                  legend.title = element_blank()) + 
            scale_fill_manual(values = c("undetected" = "#a9bcf5", "detected" = "#f5a9a9")) + 
            labs(x = "", 
                 y = "Percentage", 
                 title = paste0("Proportion cells expressed ",gene, " ", substring))   
        return(p)
    }else{echo("they are not the count value")}
}




my_cortest <- function(mat = NULL,
                       gene = NULL, 
                       method = "spearman", # "kendall","pearson" 
                       p.adjust.method = "BH",
                       n_cores = 5, 
                       rm.zero = TRUE) {

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
        
        data <- data.frame("g1" = mat[gene, ], 
                           "g2" = mat[g, ])
        if(rm.zero == TRUE){
            data <- data[-which(data[,"g1"]  == 0 & data[,"g2"] == 0),]
        }else{data <- data}

        if(nrow(data) >= 2){ 
            ct <- cor.test(data[,"g1"], data[,"g2"], method = method)
            data.frame(gene1 = gene, 
                   gene2 = g, 
                   corr = ct$estimate,
                   pval = ct$p.value, 
                   method = method,
                   ncell = nrow(data),
                   stringsAsFactors = FALSE)
        }else{ 
            data.frame(gene1 = gene, 
                   gene2 = g, 
                   corr = NA,
                   pval = NA, 
                   method = NA,
                   ncell = NA,
                   stringsAsFactors = FALSE)
        
        }
        
        
    }
    
    result_list <- future.apply::future_lapply(other_genes, compute_one)
    result <- do.call(rbind, result_list)
    result$pval_adj <- p.adjust(result$pval, 
                              method = p.adjust.method)
    
    result <- result[order(result$pval_adj), ]
    return(result)
}



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

getcolors = function(){ 


	colors = c(
        '10X scRNA'="pink",   # 深蓝
    'Smart-seq2'="#117A65",     # 深绿
    'Smart-Seq2'="#117A65",     # 深绿
    '10X snRNA'="#b9c2dd",       # 深橙棕
    
    "Whole_Embryo"    ="#FF6F61",  # 鲜橙红
    "Anterior"        ="#88B04B",  # 明亮绿
    "Trunk"           ="#92A8D1",  # 明亮蓝
    "Post otic"       ="#6B5B95",  # 明亮紫
    "Posterior"       ="#F7CAC9",  # 浅粉
    "whole bowel"     ="#C94C4C",  # 鲜红棕
    "Stomach"         ="#F7786B",  # 鲜红
    "Small intestine" ="#034F84",  # 明亮蓝绿
    "Colon"           ="#FFB347",  # 鲜橙黄
    "Hindgut"         ="#B5E7A0",  # 浅绿
    "Ileum"           ="#46B3E6",  # 明亮天蓝
    "Duodenum"        ="#F4E285",   # 明亮黄
    "Duodenum"        ="#F4E285",   # 明亮黄
    
    "https=//marionilab.cruk.cam.ac.uk/MouseGastrulation2018/"="#e41a1c", # 鲜红
    "GSE129114" ="#377eb8",  # 鲜蓝
    "CRA011744" ="#4daf4a",  # 鲜绿
    "in-house"  ="#984ea3",  # 鲜紫
    "GSE154007" ="#ff7f00",  # 橙色
    "GSE149524" ="#ffff33",  # 亮黄
    "GSE156905" ="#a65628",  # 棕色
    "GSE184981" ="#f781bf",  # 粉色
    "GSE263422" ="#999999",  # 灰色
    "http=//mousebrain.org/"="#00CED1", # 亮蓝绿
    "SCP1038"   ="#d95f02",  # 橙红
    "GSE102962" ="#1b9e77",  # 青绿
    "GSE153192" ="#A9E2F3",   # 紫色
    "GSE156146"="#0099FF",
    "GSE213604"="pink",
    

    "Progenitor"="#FF3333",
    "Neuronal lineage"="#6A3D9A",
    "Glial lineage"="#0174df",
    "Neural progenitor"="#FF66CC",
    "BP"="#9933FF",
    "GP"="#FF9900",
    "Glia"="#00E5EE",
    "Glia_Neuron"="#CCCC00", 
    "E-MTAB-5553"= "black",

    #### neuronal lineage
    "nNOS-cre;R26R-TdTomato;ChAT-EGFP-L10A"="#9A66C4",  # 紫色
    "Retrograde tracer"="#C199E0",              # 中灰色
    "Wnt1-Cre;Rosa26-YPF+"="#E6C1F2",           # 青色
    
    
    #### Glial lineage
    "Plp1==GFP"="#66AFFF",                      # 春绿色
    
    
    ### Glia, neuron, and progenitor
    "Baf53b-Cre;R26R-Tomato"="#FF7F0E",           # 深橙色
    "Wnt1-Cre;R26R-Tomato"="#F4D35E",           # 深红色
    "Wnt1-Cre;YPF"="#FFC107",                   # 明黄色
    "Wnt1-Cre;Sox10-Cre;Uchl1-H2B mCherry"="#FFD700",  # 棕色
    
    ### glia and neuron
        "Phox2b-CFP" = "#d8f6ce",                     # 浅绿色
        "Wnt1-Cre(cre/wt);H2B-mCherry(+)" = "#82947c", # 绿色
    "Tyr-Cre+_(G+N)"="#04a256",
    ### neuron and progenitor
    "Sox10+"="#2CA02C",                         # 鲜绿色
    "Whole_Embryo"="#66C266",                    # 亮蓝色
    "unsorting"="#BDBDBD",

    "Baf53b-Cre_(N)" = "#de8ef6", 
    "Tyr_(G+N)"="#04a256",
    "Baf53b_(N)" = "#de8ef6", 
    'nNOS-cre;R26R-TdTomato;ChAT-EGFP-L10A_(N)' = "#E49CEA",
    'Plp1==GFP_(G)' = "#43ade5",
    'Plp1_(G)' = "#43ade5",
    'Wnt1-Cre;Sox10-Cre;Uchl1-H2B mCherry_(P+N)' = "#9be14f",
    'Uchl1_(N)'= "#642EFE",
    'Wnt1_(P+G+N)' = "#facc2e",
    'Sox10_(P+G)' = "#fa5882",
    "Nos_(N)"="#a88afc", 
    "Chat_(N)"="#a9a9f5", 
    "Phox2b_(N)"="#fbecff", 
    


    ############# 
    "E9.5"= "#f7fcf0",  
    "E10.5"= "#e0f3db", 
    "E12.5"= "#ccebc5", 
    "E13.5"= "#7bccc4", 
    "E14.5"= "#4eb3d3", 
    "E15.5"= "#2b8cbe", 
    "E17.5"= "#0868ac", 
    "E18.5"= "#084081",

    "P5"= "#fff7f3",
    "P7"= "#fde0dd",
    "P10"= "#fcc5c0",
    "P14"= "#fa9fb5",
    "P16-18"= "#f768a1",
    "P19_P21"= "#dd3497",
    "P20"= "#ae017e",
    "P21"= "#7a0177",
    "P21_P23"= "#49006a",
    "P24"= "#fed976",        
    "P60"= "#fecc5c",
    "W7"= "#fd8d3c",
    "P180"= "#f03b20",
    "Adult"= "#bd0026",

    ############# 
    "ENC1"= "#b2ddeb",
    "ENC2"= "#9db5cb",
    "ENC3"= "#65a6d9",
    "ENC4"= "#5ac4c3",
    "ENC5"= "#d94e91",
    "ENC6_Nog-"= "#7CCD7C",
    "ENC6_Nog+"= "#B3EE3A",
    
    "ENC7_Kcna5-"= "#FFD700",
    "ENC7_Kcna5+"= "#F4FA58",
    
    "ENC8"= "#e8957a",
    "ENC9"= "#f27d53",
    "ENC10"= "#cbb894",
    "ENC11"= "#c980b6",
    "ENC12"= "#58bc7e", 
    "Glia"="#1E90FF",
    "Mesoderm-derived neuron"="#3c813c",
    "Neuroblast"="#8B658B",

    ##########
'GSE263422'= '#B71C1C',         # Dark Red
    'GSE149524'= '#1B5E20',         # Dark Green
    'SCP1038'= '#F57F17',           # Dark Yellow
    'CRA011744'= '#0D47A1',         # Dark Blue
    'GSM7747465'= '#4A148C',        # Dark Purple
    'in-house'= '#FF1744',          # Bright Red
    'GSE213604'= '#00E676',         # Bright Green
    'GSE156905'= '#FFFF00',         # Bright Yellow
    'http=//mousebrain.org/'= '#03A9F4',  # Bright Blue
    'GSE102962'= '#E040FB',         # Bright Purple
    'GSE154007'= '#FF8A80',         # Light Red
    'GSE156146'= '#B9F6CA',         # Light Green
    'GSE184981'= '#FFF9C4',         # Light Yellow
    'GSE129114'= '#B3E5FC',         # Light Blue
    'E-MTAB-5553'= '#E1BEE7',        # Light Purple 

    "Embryonic_mid"="#A9E2F3", 
    "Embryonic_late"="#5882FA", 
    "Postnatal"="#BE81F7", 
    "Adult"="#FA5882",
	# "neuronal progenitor" = "#1f78b4", 
	"neuroblast" = "#fb9a99",
	"neurons" = "#ff7f00",
	"Doublet" = "#BE81F7",
	"Singlet" = "#aaaaaa",
	"Progenitor" = "#B3DE69", 
	"SMC" = "#E6E6FA", 
	"Neuronal" = "#FCCDE5",

	"transient stage" = "#ceecf5",
	"committed stage" = "#68b8d0",
      # "cycling neuronal progenitor"='#1f77b4',
    "neuronal progenitor 5"='#e1b829',

	"Mutant"="#f82401", 
	"Vcl cKO"="#f82401", 
	"Control" = "#0f7dbd",
     "OTHERS"='grey',

	"hENCC" = "#ace9fb",
	"hENCC1" = "#ace9fb", 

	"Day9" = "#57a6bd", 
	"Day20" = "#476d78",

	"Vagal NCC" = "#385450",
	"hENCC_2" = "#8dd3c7", 
	"hENCC_1" = "#385450", 

	# "hNP_day9" = "#bebada", 
	# "hNP_day20" = "#fed976", 
	# "HDAC1_KD_Day9" = "#766eb1", 
	# "HDAC1_KD_Day20" = "#fd8d3c", 
	# "HDAC1_2_KD_Day20" = "#b10026", 

	"hNP_day9"= "#bebada", 
	"hNP_day20"="#fcd90c",
	"hNP_Day20"="#fcd90c", 

	"HDAC1_KD_Day9"= "#4000FF", 
	"HDAC1_KD_day9"= "#4000FF",  
	"HDAC1_KD_day20"="#f75716", 
	"HDAC1_KD_Day20"="#f75716", 
	"HDAC1_2_KD_Day20"="#8a1518",
    "HDAC1_2_KD_day20"="#8a1518",
	"ctrl_ENCC" = "#F5A9E1", 
	"ctrl_N_D20" = "#9A2EFE", 
	"ctrl_N_D40" = "#00BFFF", 
	"ctrl_iPSC" = "#23715d", 

	"IMR_ENCC" = "#21963e", 
	"IMR_N_Diff_D20" = "#9A2EFE", 
	"IMR_N_Diff_D40" = "#F5A9A9", 
	"IMR_iPSC" = "#00FF80", 

	"0" = "#8ea9f8",  # Blue
	"1" = "#1F77B4",  # Blue
	"2" = "#FF7F0E",  # Orange
	"3" = "#2CA02C",  # Green
	"4" = "#D62728",  # Red
	"5" = "#9467BD",  # Purple
	"6" = "#8C564B",  # Brown
	"7" = "#E377C2",  # Pink
	"8" = "#7F7F7F",  # Gray
	"9" = "#BCBD22",  # Olive
	"10" = "#17BECF", # Cyan
	"11" = "#FFBB78", # Light Orange
	"12" = "#98DF8A", # Light Green
	"13" = "#FF9896",  
   "14" = "#C5B0D5", # 浅淡紫 (Light Purple - 呼应 C5 但更亮)
  "15" = "#C49C94", # 浅褐 (Light Brown - 呼应 C6)
  "16" = "#F7B6D2", # 浅粉 (Light Pink - 呼应 C7)
  "17" = "#DBDB8D", # 浅橄榄绿 (Light Olive - 呼应 C9)
  "18" = "#9EDAE5", # 浅天蓝 (Light Cyan - 呼应 C10)
  "19" = "#E7BA52", # 金色 (Golden - 新色调，区分于 Orange)
  "20" = "#AD494A",  # 砖红 (Brick Red - 区分于 C4)
  "Vagal" = "pink", 
  "Cranial" = "#A9E2F3", 
  "NC" = "#85b0a5",

  "cycling neuronal progenitor" = "#1f77b4", # 亮红色 (代表活跃分裂)
  "neuronal progenitor"         = "#c0dffb", # 橙色
  "cycling neuroblast"          = "#34634c", # 金黄色
	"neuroblast2"="#98df8a",
	"neuroblast3"='#a9dd98',
	"neuroblast1"="#3e9621",
  "ENC(KCNJ6_PHOX2B)"           = "#a9bcf5", # 天蓝色
  "IN?(NTNG1_PAX5)"             = "#c49c94", # 宝蓝色
  "IPAN/IN(LHX1_GALNTL6)"       = "#9966FF", # 薰衣草紫
  "ENC(ST18_high)"              = "#dce14f", # 亮紫色
  "ENC(GAD2_SLC6A5)"            = "#ff9896", # 桃粉色
  "ENC(ADCY8_HES6)"             = "#F7BE81", # 玫红色
  "ENC(ISL1_SST)"               = "#00BFFF", # 深天蓝
  "ENC(C12orf57)"               = "#e2a9f3",  # 棕红色 (特别标注你关注的 C12orf57)



	"C0" = "#8ea9f8",  # Blue
	"C1" = "#1F77B4",  # Blue
	"C2" = "#FF7F0E",  # Orange
	"C3" = "#2CA02C",  # Green
	"C4" = "#D62728",  # Red
	"C5" = "#9467BD",  # Purple
	"C6" = "#8C564B",  # Brown
	"C7" = "#E377C2",  # Pink
	"C8" = "#7F7F7F",  # Gray
	"C9" = "#BCBD22",  # Olive
	"C10" = "#17BECF", # Cyan
	"C11" = "#FFBB78", # Light Orange
	"C12" = "#98DF8A", # Light Green
	"C13" = "#FF9896",  
  "C14" = "#C5B0D5", # 浅淡紫 (Light Purple - 呼应 C5 但更亮)
  "C15" = "#C49C94", # 浅褐 (Light Brown - 呼应 C6)
  "C16" = "#F7B6D2", # 浅粉 (Light Pink - 呼应 C7)
  "C17" = "#DBDB8D", # 浅橄榄绿 (Light Olive - 呼应 C9)
  "C18" = "#9EDAE5", # 浅天蓝 (Light Cyan - 呼应 C10)
  "C19" = "#E7BA52", # 金色 (Golden - 新色调，区分于 Orange)
  "C20" = "#AD494A",  # 砖红 (Brick Red - 区分于 C4)

	"cycling neuronal progenitor 1"='#1f77b4',
	"cycling neuronal progenitor 2"='#aec7e8',

	"neuronal progenitor 1"='#ff7f0e', 
	"neuronal progenitor 2"='#FFDE59',
	"neuronal progenitor 4"='#a0a156',
	"neuronal progenitor 3"='#e1b829',

	"neuroblast 3"="#98df8a",
	"neuroblast 1"='#34634c',
	"neuroblast 2"="#3e9621",

	"ENC (SST+)"='#ff9896',
	"IPAN/IN"='#8c564b',
	"IN?"='#c49c94',
	"ENC (GAD2+)"='#e377c2',
	# "ENC (PHOX2B+_KCNJ6+)"='#f7b6d2',
	"ENC (PHOX2B+_KCNJ6+)"='#D20103',
	# "ENC (TH+_PHOX2B+)"='#7f7f7f',
	# "ENC (TH+_PHOX2B+)"='#060270',
	"ENC (NOVA1_high)"='#c7c7c7',

	"BP" = "#6b853e", 
	"BPearly" = "#94c243", 
	"BPlate" = "#607838", 
	"GP" = "#feb462", 
	"GPlate" = "#FFDE59",
	"NPearly" = "#fdcee6", 
	"NPlate" = "#ff2f9a", 
	"ENMFB" = "#d9d9d9",
	"Neuroblast" = "#fdcee6",
	"Branch A" = "#bca9f5", 
	"Neuroblast" = "#E41A1C",   
	"Branch B" = "#f180ba",
	"GPearly" = "#D7DF01",
	"GPlate" = "#fbb469",
    "ENMFB" = "#D9D9D9",
    "Others" = "#E7D040",
    "GP" = "#FD9D52",
    "Differentiating GP" = "#F9D151",
    "BP" = "#647E3B",
    "Cycling BP" = "#1F8A1F",
    "Cycling Neuroblast" = "#FBEDF3",
    "Neuroblast" = "#F2AFD1",
    "Branch B" = "#E661AC",
    "Branch A" = "#B4B1D9",
    "Differantiating GP" = "#F9D151"
	)
	return(colors)
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

my_volcan.plt <- function(
				meta, 
				cutoff_lfc, 
				cutoff_qval, 
				x, 
                rev.result = FALSE,
				y, 
				label = TRUE, 
				label.var = "label",
				max.overlaps = 20,
				title = "DEGs",
				label.size = 2.5, 
				color.var = "avg_log2FC"
                   ){ 
    library(ggrepel)
    library(ggplot2)

    if(rev.result){
        meta[,x] <- -meta[,x]
    }

	down.temp <- meta %>% dplyr::filter(.data[[x]] < -cutoff_lfc & .data[[y]] < cutoff_qval)
	up.temp <- meta %>% dplyr::filter(.data[[x]] > cutoff_lfc & .data[[y]] < cutoff_qval)

	label.data <- data.frame(
		"num" = c(nrow(down.temp), nrow(up.temp)),
		"x" = c(min(down.temp[,x]) * 0.8, max(up.temp[,x]) * 0.8),
		"y" = rep(max(c(max(-log10(down.temp[,y])) * 0.8, max(-log10(up.temp[,y])) * 0.8),2)),
		"col" = c("blue","red")

	)

	if(label){
		p <- ggplot(meta, aes(x = .data[[x]], y = -log10(.data[[y]]))) + 
	        geom_hline(yintercept = -log10(cutoff_qval), linetype = "dashed", color = "#999999") +
	        geom_vline(xintercept = c(-cutoff_lfc,cutoff_lfc), linetype = "dashed", color = "#999999") + 
	        geom_point(size = 1, color = ifelse(meta[,color.var] < -cutoff_lfc, "blue", ifelse(meta[,color.var]  > cutoff_lfc, "red", "gray"))) + 
			geom_text_repel(data = meta, 
                        label = meta[,label.var], 
                        max.overlaps = max.overlaps, 
                        size = label.size, 
                        box.padding = 0, 
                        nudge_x = 0.05, 
                        nudge_y = 0.1)+
	        theme_bw(base_size = 14)+
	        ggsci::scale_color_jama() + 
			geom_text(data = label.data, aes(x = x, y = y, label = num), color = ifelse(label.data$col == "red", "red", "blue"),size = 6) +
	        theme(panel.grid = element_blank(),
	              legend.position = 'right',
				 axis.text = element_text(size = 12, color = "black"),
				 axis.title = element_text(size = 12, color = "black"),
				  plot.title = element_text(hjust = 0.5)) + 
	        labs(x = "Log2FC", y = "-Log10(q-value)", title = title)                                              
        return(p)
		
	}else{
		p <- ggplot(meta, aes(x = .data[[x]], y = -log10(.data[[y]]))) + 
	        geom_hline(yintercept = -log10(cutoff_qval), linetype = "dashed", color = "#999999") +
	        geom_vline(xintercept = c(-cutoff_lfc,cutoff_lfc), linetype = "dashed", color = "#999999") + 
	        geom_point(size = 1, color = ifelse(meta[,color.var] < -cutoff_lfc, "blue", ifelse(meta[,color.var]  > cutoff_lfc, "red", "gray"))) + 
	        theme_bw(base_size = 14)+
	        ggsci::scale_color_jama() +
			geom_text(data = label.data, aes(x = x, y = y, label = num), color = ifelse(label.data$col == "red", "red", "blue"),size = 6) +
	        theme(panel.grid = element_blank(),
				 axis.text = element_text(size = 12, color = "black"),
				 axis.title = element_text(size = 12, color = "black"),
	              legend.position = 'right', 
				  plot.title = element_text(hjust = 0.5)) + 
	        labs(x = "Log2FC", y = "-Log10(q-value)", title = title)                                              
        return(p)
	}
}    

#####

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
  meta, x, y, color = NULL,fill = NULL,
  title = "cytoTRACE",
  subtitle = NULL,
  base_size = 14,
  ylab = "Potency score",
  facet = TRUE,
  facet.var = NULL,
  color.fill = NULL,
  color.color = NULL,
  mean.value.size = 4,
  mean.digit = 2,
  scale = TRUE,
  label.size = 5,
  label_y_positions = rep(max(meta[[y]], na.rm = TRUE) * 1.05, length(comparisons)),
  label_x_positions = c(1.5, 2.5),
  show.mean = TRUE,
  comparison = NULL,  # list(c("Group1", "Group2"), c("Group2", "Group3"))
  legend.pos = "none",
  p.format = TRUE,     # 是否格式化 p 值（如 0.00123 → 1.23e-03）
  p.digits = 4        # 小数位数或科学计数法精度
) {
  library(ggplot2)
  library(ggsignif)
  library(dplyr)

  # 判断 color 参数，并填充color映射
  p <- ggplot(meta, aes(x = .data[[x]], y = .data[[y]])) +
    {
      if (!is.null(color)) {
        geom_violin(aes(color = .data[[color]], fill = .data[[fill]]), alpha = 0.8)
      } else {
        geom_violin(fill = "grey", color = "black", alpha = 0.8)
      }
    } +
    # 十字表示均值，也考虑 color 是否为空
    stat_summary(fun = mean, geom = "crossbar", width = 0.3, size = 0.1,color = "black"
				 ) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = legend.pos,
      axis.text.y = element_text( size = base_size),
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size)
    ) +
    ggtitle(title) + labs(x = NULL, y = ylab, subtitle = subtitle) +  ylim(min(meta[[y]])*1.2, max(meta[[y]])*1.5) 

  # 分面
  if (facet && !is.null(facet.var)) {
    p <- p + facet_grid(. ~ .data[[facet.var]], scales = if (scale) "free_x" else "fixed")
  }

  # 配色
  if (!is.null(colors) & !is.null(color)) {
    p <- p + scale_color_manual(values = color.color) +
             scale_fill_manual(values = color.fill)
  }

  # 均值标签
  if (show.mean) {
    p <- p + stat_summary(
      fun = mean, geom = "text",
      aes(label = round(after_stat(y), 4)),
      color = "black", size = mean.value.size, vjust = -0.5
    )
  }

  # 添加 p 值标注
  if (!is.null(comparison) && length(comparison) > 0) {
    label_y_positions_full <- rep(label_y_positions, length.out = length(comparison))
    for (i in seq_along(comparison)) {
      grp_pair <- comparison[[i]]
      g1 <- grp_pair[1]
      g2 <- grp_pair[2]

      # 避免因分组不存在而报错
      idx1 <- which(meta[[x]] == g1)
      idx2 <- which(meta[[x]] == g2)
      if (length(idx1) == 0 || length(idx2) == 0) next

      mean1 <- mean(meta[[y]][idx1], na.rm = TRUE)
      mean2 <- mean(meta[[y]][idx2], na.rm = TRUE)

      # Wilcoxon 检验
      test_res <- tryCatch({
        wilcox.test(meta[[y]][idx1], meta[[y]][idx2], alternative = "two.sided", paired = FALSE)
      }, error = function(e) return(list(p.value = NA)))
      p_val <- test_res$p.value

      # 格式化P值
      if (is.na(p_val)) {
        p_text <- "p = NA"
      } else if (p.format) {
        if (p_val < 1e-3) {
          p_text <- paste0(formatC(p_val, format = "e", digits = p.digits))
        } else {
          p_text <- paste0(signif(p_val, p.digits))
        }
      } else {
        p_text <- paste0(p_val)
      }

      # 颜色判断
      text_color <- ifelse((mean1 - mean2) > 0, "blue", "red")

      # y 轴位置
      y_pos <- if (i <= length(label_y_positions_full)) {
        label_y_positions_full[i]
      } else {
        max(meta[[y]], na.rm = TRUE) * 1.05 + 0.08 * (i - length(label_y_positions_full))
      }

      # 实际添加标注
      p <- p + geom_signif(
        comparisons = list(grp_pair),
        y_position = y_pos,
        tip_length = 0.01,
        vjust = 0,
        textsize = label.size,
        map_signif_level = FALSE,
        annotations = p_text,
        color = text_color
      ) &  ylim(min(meta[[y]])*1.5, max(meta[[y]])*1.5) 
    }
  }

  return(p)
}



my_box_compare.plt <- function(
  meta, x, y, color = NULL,fill = NULL,
  title = "cytoTRACE",
  subtitle = NULL,
  base_size = 14,
  ylab = "Potency score",
  facet = TRUE,
  facet.var = NULL,
  color.fill = NULL,
  color.color = NULL,
  mean.value.size = 4,
  mean.digit = 2,
  scale = TRUE,
  label.size = 5,
  label_y_positions = rep(max(meta[[y]], na.rm = TRUE) * 1.05, length(comparisons)),
  label_x_positions = c(1.5, 2.5),
  show.mean = TRUE,
  comparison = NULL,  # list(c("Group1", "Group2"), c("Group2", "Group3"))
  legend.pos = "none",
  p.format = TRUE,     # 是否格式化 p 值（如 0.00123 → 1.23e-03）
  p.digits = 3         # 小数位数或科学计数法精度
) {
  library(ggplot2)
  library(ggsignif)
  library(dplyr)

  # 判断 color 参数，并填充color映射
  p <- ggplot(meta, aes(x = .data[[x]], y = .data[[y]])) +
    {
      if (!is.null(color)) {
        geom_boxplot(aes(color = .data[[color]]), alpha = 0.6) 
      } else {
        geom_boxplot(aes(color = .data[[color]]), alpha = 0.6) 
      }
    } + geom_jitter(aes(color = .data[[color]]),alpha = 0.1, size = 0.2) +
    # 十字表示均值，也考虑 color 是否为空
    stat_summary(fun = mean, geom = "point", size = 1, alpha = 0.5,color = "black"
				 ) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = legend.pos,
      axis.text.y = element_text( size = base_size),
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size)
    ) +
    ggtitle(title) + labs(x = NULL, y = ylab, subtitle = subtitle) +  ylim(min(meta[[y]])*1.2, max(meta[[y]])*1.5) 

  # 分面
  if (facet && !is.null(facet.var)) {
    p <- p + facet_grid(. ~ .data[[facet.var]], scales = if (scale) "free_x" else "fixed")
  }

  # 配色
  if (!is.null(colors) & !is.null(color)) {
    p <- p + scale_color_manual(values = color.color) +
             scale_fill_manual(values = color.fill)
  }

  # 均值标签
  if (show.mean) {
    p <- p + stat_summary(
      fun = mean, geom = "text",
      aes(label = round(after_stat(y), 4)),
      color = "black", size = mean.value.size, vjust = -0.5
    )
  }

  # 添加 p 值标注
  if (!is.null(comparison) && length(comparison) > 0) {
    label_y_positions_full <- rep(label_y_positions, length.out = length(comparison))
    for (i in seq_along(comparison)) {
      grp_pair <- comparison[[i]]
      g1 <- grp_pair[1]
      g2 <- grp_pair[2]

      # 避免因分组不存在而报错
      idx1 <- which(meta[[x]] == g1)
      idx2 <- which(meta[[x]] == g2)
      if (length(idx1) == 0 || length(idx2) == 0) next

      mean1 <- mean(meta[[y]][idx1], na.rm = TRUE)
      mean2 <- mean(meta[[y]][idx2], na.rm = TRUE)

      # Wilcoxon 检验
      test_res <- tryCatch({
        wilcox.test(meta[[y]][idx1], meta[[y]][idx2])
      }, error = function(e) return(list(p.value = NA)))
      p_val <- test_res$p.value

      # 格式化P值
      if (is.na(p_val)) {
        p_text <- "p = NA"
      } else if (p.format) {
        if (p_val < 1e-3) {
          p_text <- paste0(formatC(p_val, format = "e", digits = p.digits))
        } else {
          p_text <- paste0(signif(p_val, p.digits))
        }
      } else {
        p_text <- paste0(p_val)
      }

      # 颜色判断
      text_color <- ifelse((mean1 - mean2) > 0, "blue", "red")

      # y 轴位置
      y_pos <- if (i <= length(label_y_positions_full)) {
        label_y_positions_full[i]
      } else {
        max(meta[[y]], na.rm = TRUE) * 1.05 + 0.08 * (i - length(label_y_positions_full))
      }

      # 实际添加标注
      p <- p + geom_signif(
        comparisons = list(grp_pair),
        y_position = y_pos,
        tip_length = 0.01,
        vjust = 0,
        textsize = label.size,
        map_signif_level = FALSE,
        annotations = p_text,
        color = text_color
      ) &  ylim(min(meta[[y]])*1.5, max(meta[[y]])*1.5) 
    }
  }

  return(p)
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
        # library(ggplot2)
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
        # library(ggplot2)
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
          ) + scale_color_gradient(low = "#DF3A01", high = "#F6D8CE");
    
    return(p)
}
