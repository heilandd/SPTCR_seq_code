# Functions for the SPTCR Analysis Script (by D. H. Heiland) ---------------------------------------------------------------
preprocessTCR <- function(object, TCRs, minTCRsize=5, maxTCRsize=50, min_UMI=2){
  
  TCRs <- TCRs %>% dplyr::select(index, barcodes, expression, Locus, V, D, J, CDR3_aa, v_identity,d_identity,j_identity,tcrbert_leiden)
  bcs <- SPATA2::getBarcodes(object) %>% unlist() %>% as.character()
  
  without_filter <- nrow(TCRs)
  
  # filter non spot barcodes
  org <- nrow(TCRs)
  TCRs <- TCRs %>% filter(barcodes %in% bcs)
  after <- nrow(TCRs)
  message(paste0("Percentage of valid barcodes:", (after/org)*100, "%"))
  
  
  # TCR Filter
  TCRs[TCRs=="NaN"] <- 0
  TCRs$nchar <- map(TCRs$CDR3_aa, ~nchar(.x)) %>% unlist()
  TCRs <- 
    TCRs %>% 
    filter(expression>=min_UMI) %>% 
    filter(nchar>minTCRsize) %>% 
    filter(nchar<maxTCRsize) %>% 
    filter(V!="nan") %>% 
    filter(D!="nan") %>% 
    filter(J!="nan") %>% 
    filter(CDR3_aa!="nan") %>% 
    filter(Locus=="TRB")
  #mutate(valid=(as.numeric(v_identity)+as.numeric(d_identity)+as.numeric(d_identity))/3) %>% 
  #filter(valid==100)
  
  with_filter <- nrow(TCRs)
  
  message(paste0("Percentage of TCRs after filtering: ", round((with_filter/without_filter)*100, digits = 2), "%"))
  
  return(TCRs)
  
}
preprocessLD <- function(TCRs, min_dist=3){
  
  message(paste0("Number of : ", length(unique(TCRs$CDR3_aa)), " TCRs will be analyzed"))
  message(paste0("Run CDR3 distance ... "))
  new <- map_dfr(unique(TCRs$CDR3_aa), .f=function(CDR3){
    
    dist <- 
      TCRs %>% 
      mutate(dist= stringdist::stringdist(CDR3, TCRs$CDR3_aa)) %>% 
      filter(dist<=min_dist)
    
    CDR3aa <- dist$CDR3_aa[which.max(dist$nchar)]
    
    group <- 
      dist %>% 
      group_by(barcodes) %>% 
      summarise(expression=sum(expression),
                v_identity =mean(v_identity),
                d_identity = mean(d_identity),
                j_identity = mean(j_identity)) %>% 
      ungroup()
    group$CDR3_aa <- CDR3aa
    group <- cbind(group, dist[1:nrow(group),c("Locus","V","D","J","tcrbert_leiden")])
    group$nchar <- nchar(CDR3aa)
    group$valid <- rowMeans(group[,c("v_identity","d_identity","j_identity")])
    
    group <- group[,names(dist)[c(!names(dist) %in% c("dist", "index")) ] ]
    
    return(group)
    
  }, .progress = T)
  
  
  new_TCRs <- new %>% group_by(CDR3_aa) %>% count(CDR3_aa)
  
  message(paste0("Mean number of spots per TCRs: ", mean(new_TCRs$n)))
  message(paste0("Number of : ", length(unique(new_TCRs$CDR3_aa)), " clones after summary"))
  message(paste0("------------------------------------------------"))
  message(paste0("Percentage of TCRs after filtering: ", round((nrow(new_TCRs)/nrow(TCRs))*100, digits = 2), "%"))
  
  return(new)
  
}
normalizeTCR <- function(object, spot_TCR){
  new <- spot_TCR
  # Analyse the clones that were found
  Total_UMI <- 
    SPATA2::getFeatureDf(object) %>% 
    dplyr::select(barcodes,nCount_Spatial,Nr_of_cells)
  
  Total_UMI <- mutate(Total_UMI, nCount_Spatial_per_cell=nCount_Spatial/Nr_of_cells)
  
  new_filter <- new %>% filter(expression>2)
  new_filter <- left_join(new_filter, Total_UMI, by="barcodes")
  
  message(paste0("Percentage of TCRs after filtering: ", round((nrow(new_filter)/nrow(new))*100, digits = 2), "%"))
  
  clones <- 
    new_filter %>% 
    group_by(CDR3_aa) %>% 
    summarize(reads=sum(expression/nCount_Spatial_per_cell),
              bert=(tcrbert_leiden %>% unique())[1],
              length_bert=length(tcrbert_leiden %>% unique()),
              nr_spots=length(barcodes %>% unique())) %>% 
    arrange(desc(nr_spots))
  
  
  return(clones)
}
normalizedTCRMatrix <- function(object, spot_TCR){
  new <- spot_TCR
  # Analyse the clones that were found
  Total_UMI <- 
    SPATA2::getFeatureDf(object) %>% 
    dplyr::select(barcodes,nCount_Spatial,Nr_of_cells)
  
  Total_UMI <- mutate(Total_UMI, nCount_Spatial_per_cell=nCount_Spatial/Nr_of_cells)
  
  new_filter <- 
    left_join(new, Total_UMI, by="barcodes") %>% 
    mutate(reads=log((expression/nCount_Spatial_per_cell)) %>% scales::rescale(., c(0,10)))
  
  new_mat <- reshape2::acast(CDR3_aa~barcodes, data=new_filter, value.var = "reads", fun.aggregate= mean, fill=0)
  
  
  
  
  return(new_mat)
}
Dimplot_dhh <- function(object,color,pt.b=6,pt.w=5, ...){
  
  plot <- DimPlot(object,...)
  data <- plot$data
  dim_names <- names(data)
  
  p <- 
    ggplot(data=data)+
    scattermore::geom_scattermore(data=data, mapping = aes(x=!!sym(dim_names[1]),y=!!sym(dim_names[2])),pointsize = pt.b, color="black")+
    scattermore::geom_scattermore(data=data, mapping = aes(x=!!sym(dim_names[1]),y=!!sym(dim_names[2])),pointsize = pt.w, color="white")
  
  p <- p+plot$layers[[1]]
  
  p <- p+
    scale_color_manual(values=color)+
    #scale_colour_gradientn(colours = col(50), oob = scales::squish, na.value = "white")+
    ylab("UMAP2")+xlab("UMAP1")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))+
    coord_fixed()
  
  
  return(p)
  
  
  
  
}
plotSurface_DHH <- function(object,pt.b=6,pt.w=5, ...){
  
  plot <- SPATA2::plotSurface(object, ...)
  
  data <- plot$data
  dim_names <- names(data)
  
  p <- 
    ggplot(data=data)+
    scattermore::geom_scattermore(data=data, mapping = aes(x=!!sym(dim_names[2]),y=!!sym(dim_names[3])),pointsize = pt.b, color="black")+
    scattermore::geom_scattermore(data=data, mapping = aes(x=!!sym(dim_names[2]),y=!!sym(dim_names[3])),pointsize = pt.w, color="white")
  
  p <- 
    p+plot$layers[[1]]+
    SPATA2::ggpLayerAxesSI(object)
  
  
  
  
}
add_hull <- function(object,plot,pt.b=6,pt.w=5){
  
  data <- plot$data
  dim_names <- names(data)
  
  p <- 
    ggplot(data)+
    scattermore::geom_scattermore(data=SPATA2::getCoordsDf(object), mapping = aes(x,y),pointsize = pt.b, color="black")+
    scattermore::geom_scattermore(data=SPATA2::getCoordsDf(object), mapping = aes(x,y),pointsize = pt.w, color="white")
  
  p <- 
    p+plot$layers[[1]]+
    SPATA2::ggpLayerAxesSI(object)
  
  return(p)
  
}
CytoSpace2SPATA2 <- function(object,
                             cytospace.out,
                             reference){
  
  st <- object
  
  message(paste0(Sys.time(), "--- Read in Cell Type fractions -----"))
  setwd(cytospace.out)
  sample <- st %>% SPATA2::getSampleName()
  
  #Add Cells to org SPATA2 object
  cell_types_spot <- read.csv(paste0(cytospace.out,"/cell_type_assignments_by_spot.csv"))
  names(cell_types_spot)[1] <- "barcodes"
  names(cell_types_spot)[2:ncol(cell_types_spot)] <- paste0(names(cell_types_spot)[2:ncol(cell_types_spot)], "_CS")
  
  cell_types_spot <- left_join(data.frame(barcodes=SPATA2::getBarcodes(st)[[1]]), cell_types_spot, by="barcodes")
  cell_types_spot[is.na(cell_types_spot)]=0
  st <- SPATA2::addFeatures(st, cell_types_spot,overwrite = T)
  
  message(paste0(Sys.time(), "--- Create single cell SPATA2 object -----"))
  sc_loc <- read.csv(paste0(cytospace.out,"/assigned_locations.csv"))
  
  # jitter position
  sc_loc <- sc_loc %>% left_join(., sc_loc %>% group_by(SpotID) %>% count())
  sc_loc <- map_dfr(.x=sc_loc$SpotID %>% unique(), .f=function(x){
    a <- sc_loc %>% filter(SpotID==x)
    if(a$n[1]==1){a$row.new=a$row; a$col.new=a$col}else{
      a$row.new <- runif(nrow(a), unique(a$row)-0.6, unique(a$row)+0.6) %>% sample()
      a$col.new <- runif(nrow(a), unique(a$col)-0.6, unique(a$col)+0.6) %>% sample()
    }
    return(a)
  }, .progress=T)
  
  
  #Get count mat
  mat <- Seurat::GetAssayData(reference, "counts")
  mat <- mat[,sc_loc$OriginalCID]
  mat <- mat %>% as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("OriginalCID")
  mat <- sc_loc %>% dplyr::select(UniqueCID,OriginalCID ) %>% left_join(.,mat)
  
  #UMAP
  umap <- reference@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column("OriginalCID")
  umap <- 
    sc_loc %>% 
    dplyr::select(UniqueCID,OriginalCID ) %>% 
    left_join(.,umap) %>% dplyr::select(-2) %>% 
    mutate(sample=sample) %>% 
    rename("umap1":=UMAP_1) %>% 
    rename("umap2":=UMAP_2) %>% 
    rename("barcodes":=UniqueCID) %>% 
    dplyr::select(barcodes, sample, umap1, umap2)
  
  meta_info <- left_join(sc_loc[,1:2],reference@meta.data %>% as.data.frame() %>% rownames_to_column("OriginalCID"))
  meta_info$barcodes <- meta_info$UniqueCID
  
  counts <- mat[,3:ncol(mat)]
  counts[is.na(counts)]=0
  counts[1:10, 1:10]
  
  counts <- counts %>% t() %>% Matrix::Matrix(sparse = T)
  colnames(counts) <- mat$UniqueCID
  
  coords <- data.frame(barcodes=sc_loc$UniqueCID, x=sc_loc$col.new, y=sc_loc$row.new)
  
  
  message(paste0(Sys.time(), "--- Reshape space of single cells -----"))
  coords_org <- getCoordsDf(st)
  coords$x <- scales::rescale(coords$x, c(min(coords_org$x), max(coords_org$x)))
  coords$y <- scales::rescale(coords$y, c(min(coords_org$y), max(coords_org$y)))
  coords_df <- coords
  
  
  
  message(paste0(Sys.time(), "--- Initiate object -----"))
  sc_obj <- SPATA2::initiateSpataObject_CountMtr(coords_df = coords_df, count_mtr = counts, sample_name = sample,
                                                 RunPCA=T, FindNeighbors=F, FindClusters=F,
                                                 RunTSNE=F, RunUMAP=F)
  
  sc_obj <- SPATA2::addFeatures(sc_obj, data.frame(barcodes=sc_loc$UniqueCID, Cell_lables=sc_loc$CellType, SpotID=sc_loc$SpotID))
  sc_obj@used_genesets <- st@used_genesets
  img <- st@images
  img[[sample]]@coordinates <- SPATA2::getCoordsDf(sc_obj)
  sc_obj@images <-img
  sc_obj <- SPATA2::flipCoordinates(sc_obj, axis="x")
  
  #Adapt the umap coords
  sc_obj@dim_red[[sample]]$umap <- umap
  names(meta_info) <- paste0(names(meta_info), "_meta")
  meta_info$barcodes <- meta_info$UniqueCID_meta
  
  sc_obj <- SPATA2::addFeatures(sc_obj, meta_info)
  
  #Add method info 
  sc_obj@information$pxl_scale_fct <- st@information$pxl_scale_fct
  sc_obj@information$method <- st@information$method
  
  return(list(st, sc_obj))
  
}
runHistologyAnnotation <- function(object, geneSet){
  
  # Get IVY expression clusters
  object@used_genesets <- rbind(object@used_genesets, geneSet)
  gs <- geneSets$ont %>% unique()
  ls <- joinWithGeneSets(object, gene_sets = gs) %>% select(barcodes, {{gs}})
  ls <- ls %>% mutate(Histology_IVY= map(1:nrow(ls), ~which.max(ls[.x,gs])) %>% unlist() %>% names() )
  
  # Remove the Spots without Tumor or change to Infiltartive 
  CNV <- joinWith(object,  features = c("Chr7", "Chr10")) %>% filter(Chr7<1.02 & Chr10>0.98) %>% pull(barcodes)
  ls[ls$barcodes %in% CNV, ]$Histology_IVY <- "IVY_IT"
  
  # Remove Necrosis 
  nec <- joinWith(object,  features = "nCount_Spatial") %>% filter(nCount_Spatial < quantile(nCount_Spatial, 0.02)) %>% pull(barcodes)
  ls[ls$barcodes %in% nec, ]$Histology_IVY <- "IVY_Necrosis"
  
  object <- SPATA2::addFeatures(object,ls %>% select(barcodes, Histology_IVY), overwrite = T)
  
  return(object)
  
}
getClonaladjacency <- function(object, spot_TCR){
  
  ## Get adjacent matrix by Triangulation
  DI <- RCDT::delaunay(SPATA2::getCoordsDf(object)[,c("x", "y")] %>% as.matrix())
  nn=data.frame(from=SPATA2::getCoordsDf(object)$barcodes[DI$edges[,1]],
                to=SPATA2::getCoordsDf(object)$barcodes[DI$edges[,2]])
  
  clone_NN <- map_dfr(unique(spot_TCR$CDR3_aa), .f=function(CDR3){
    TCR_unique <- 
      spot_TCR %>% 
      filter(CDR3_aa==CDR3) %>% 
      distinct(barcodes, .keep_all = T) 
    
    joinNN <- 
      TCR_unique %>% 
      mutate(from=barcodes) %>% 
      left_join(., nn, by="from") %>% 
      group_by(barcodes) %>% 
      summarise(NN=list(to))
    
    NN <- map(joinNN$barcodes, .f=function(bc){
      matches <- intersect(bc, joinNN$NN %>% unlist() %>% unique())
      if(is_empty(matches)){matches=0}else{matches=length(matches)}
    })
    
    TCR_unique$NN <- NN %>% unlist()
    
    return(TCR_unique)
    
    
  }, .progress = T)
  
  
  return(clone_NN)
  
  
  
}
getSpotwiseDiversity <- function(object, spot_TCR){
  
  feature_level <- map_dfr(unique(spot_TCR$barcodes), .f=function(bc){
    
    spot_level <- spot_TCR %>% filter(barcodes==bc)
    mean_string_dist <- map(unique(spot_level$CDR3_aa), function(CDR3){
      dist= mean(stringdist::stringdist(CDR3, spot_level$CDR3_aa))
    }) %>% unlist() %>% mean()
    
    out <- data.frame(barcodes=bc, Diversity_index=mean_string_dist)
    return(out)
    
  },.progress=T)
  return(feature_level)
  
}
createCellGraph <- function(object, features){
  
  Axis <- SPATA2::joinWith(object, features = features)
  DI <- RCDT::delaunay(Axis[,c("x", "y")] %>% as.matrix())
  nn=data.frame(from=Axis$barcodes[DI$edges[,1]],
                to=Axis$barcodes[DI$edges[,2]])
  nn$barcodes <- nn$to
  nn <- left_join(nn, Axis)
  nn$to <- nn[,features]
  nn <- nn[,1:2]
  nn$barcodes <- nn$from
  nn <- left_join(nn, Axis)
  nn$from <- nn[,features]
  nn <- nn[,c("from", "to")]
  return(nn)
}
plotColorOverlap_RL <- function (object,R_mat, L_mat, R, L, pt_size = 3, get.map = F, R_enh=0, L_enh=0,
                                 as.layer = F, two.colors = c("purple", "lightgreen"), negative.color = "darkgrey", 
                                 col.threshold = 0.8, smooth = F, smooth_span = 0.2, normalize = T) 
{
  
  
  #R
  object@information$active_mtr=R_mat
  data_1<- SPATA2::hlpr_join_with_aes(object, df = SPATA2::getCoordsDf(object), 
                                      color_by = R, normalize = normalize, smooth = smooth, 
                                      smooth_span = smooth_span)
  
  #L
  object@information$active_mtr=L_mat
  data_2 <- SPATA2::hlpr_join_with_aes(object, df = SPATA2::getCoordsDf(object), 
                                       color_by = L, normalize = normalize, smooth = smooth, 
                                       smooth_span = smooth_span)
  
  features <- c(R,L)
  
  data <- data_1 %>% mutate(L=data_2 %>% pull(!!sym(L)))
  names(data)[8] <- L
  
  data[,R] <- data[,R]+R_enh
  data[,L] <- data[,L]+L_enh
  
  
  library(reshape2)
  BlendExpression <- function(data) {
    if (ncol(x = data) != 2) {
      stop("'BlendExpression' only blends two features")
    }
    features <- colnames(x = data)
    data <- as.data.frame(x = apply(X = data, MARGIN = 2, 
                                    FUN = function(x) {
                                      return(round(x = 9 * (x - min(x))/(max(x) - min(x))))
                                    }))
    data[, 3] <- data[, 1] + data[, 2] * 10
    colnames(x = data) <- c(features, paste(features, collapse = "_"))
    for (i in 1:ncol(x = data)) {
      data[, i] <- factor(x = data[, i])
    }
    return(data)
  }
  BlendMatrix <- function(n = 10, col.threshold = 0.5, two.colors = c("#ff0000", 
                                                                      "#00ff00"), negative.color = "black") {
    if (0 > col.threshold || col.threshold > 1) {
      stop("col.threshold must be between 0 and 1")
    }
    C0 <- as.vector(col2rgb(negative.color, alpha = TRUE))
    C1 <- as.vector(col2rgb(two.colors[1], alpha = TRUE))
    C2 <- as.vector(col2rgb(two.colors[2], alpha = TRUE))
    blend_alpha <- (C1[4] + C2[4])/2
    C0 <- C0[-4]
    C1 <- C1[-4]
    C2 <- C2[-4]
    merge.weight <- min(255/(C1 + C2 + C0 + 0.01))
    sigmoid <- function(x) {
      return(1/(1 + exp(-x)))
    }
    blend_color <- function(i, j, col.threshold, n, C0, C1, 
                            C2, alpha, merge.weight) {
      c.min <- sigmoid(5 * (1/n - col.threshold))
      c.max <- sigmoid(5 * (1 - col.threshold))
      c1_weight <- sigmoid(5 * (i/n - col.threshold))
      c2_weight <- sigmoid(5 * (j/n - col.threshold))
      c0_weight <- sigmoid(5 * ((i + j)/(2 * n) - col.threshold))
      c1_weight <- (c1_weight - c.min)/(c.max - c.min)
      c2_weight <- (c2_weight - c.min)/(c.max - c.min)
      c0_weight <- (c0_weight - c.min)/(c.max - c.min)
      C1_length <- sqrt(sum((C1 - C0)^2))
      C2_length <- sqrt(sum((C2 - C0)^2))
      C1_unit <- (C1 - C0)/C1_length
      C2_unit <- (C2 - C0)/C2_length
      C1_weight <- C1_unit * c1_weight
      C2_weight <- C2_unit * c2_weight
      C_blend <- C1_weight * (i - 1) * C1_length/(n - 1) + 
        C2_weight * (j - 1) * C2_length/(n - 1) + (i - 
                                                     1) * (j - 1) * c0_weight * C0/(n - 1)^2 + C0
      C_blend[C_blend > 255] <- 255
      C_blend[C_blend < 0] <- 0
      return(rgb(red = C_blend[1], green = C_blend[2], 
                 blue = C_blend[3], alpha = alpha, maxColorValue = 255))
    }
    blend_matrix <- matrix(nrow = n, ncol = n)
    for (i in 1:n) {
      for (j in 1:n) {
        blend_matrix[i, j] <- blend_color(i = i, j = j, 
                                          col.threshold = col.threshold, n = n, C0 = C0, 
                                          C1 = C1, C2 = C2, alpha = blend_alpha, merge.weight = merge.weight)
      }
    }
    return(blend_matrix)
  }
  BlendMap <- function(color.matrix) {
    color.heat <- matrix(data = 1:prod(dim(x = color.matrix)) - 
                           1, nrow = nrow(x = color.matrix), ncol = ncol(x = color.matrix), 
                         dimnames = list(1:nrow(x = color.matrix), 1:ncol(x = color.matrix)))
    xbreaks <- seq.int(from = 0, to = nrow(x = color.matrix), 
                       by = 2)
    ybreaks <- seq.int(from = 0, to = ncol(x = color.matrix), 
                       by = 2)
    color.heat <- melt(color.heat)
    names(color.heat) <- c("rows", "cols", "vals")
    color.heat$rows <- as.numeric(x = as.character(x = color.heat$rows))
    color.heat$cols <- as.numeric(x = as.character(x = color.heat$cols))
    color.heat$vals <- factor(x = color.heat$vals)
    plot <- ggplot(data = color.heat, mapping = aes_string(x = "rows", 
                                                           y = "cols", fill = "vals")) + geom_raster(show.legend = FALSE) + 
      theme(plot.margin = unit(x = rep.int(x = 0, times = 4), 
                               units = "cm")) + scale_x_continuous(breaks = xbreaks, 
                                                                   expand = c(0, 0), labels = xbreaks) + scale_y_continuous(breaks = ybreaks, 
                                                                                                                            expand = c(0, 0), labels = ybreaks) + scale_fill_manual(values = as.vector(x = color.matrix)) + 
      theme_classic()
    return(plot)
  }
  BlendExpression(data[, features])
  color.matrix <- BlendMatrix(n = 50, two.colors = two.colors, 
                              col.threshold = col.threshold, negative.color = negative.color)
  if (get.map == T) {
    p = BlendMap(color.matrix)
  }
  else {
    a <- scales::rescale(data %>% pull(R), to = c(1, 
                                                  ncol(color.matrix))) %>% as.integer()
    b <- scales::rescale(data %>% pull(L), to = c(1, 
                                                  nrow(color.matrix))) %>% as.integer()
    color = map_chr(.x = 1:nrow(data), .f = function(i) {
      color.matrix[a[i], b[i]]
    })
    if (as.layer == T) {
      
      p = geom_point(data = data, mapping = aes(x = x, 
                                                y = y), size = pt_size, color = color)
    }
    else {
      p = ggplot(data, mapping = aes(x = x, y = y)) + 
        scattermore::geom_scattermore(aes(x = x, y = y),pointsize = 7, color="black")+
        scattermore::geom_scattermore(aes(x = x, y = y),pointsize = 6, color="white")+
        geom_point(size = pt_size,  color = color) + coord_fixed() + theme_void()
    }
  }
  return(p)
}
