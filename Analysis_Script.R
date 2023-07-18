## Analysis Script of the SPTCR-seq Data analysis
## By D. H. Heiland, Microenvironment and Immunology Research Laboratory
## The ST data can be downloaded via SPATAData: https://github.com/theMILOlab/SPATAData
## The SPTCR-seq raw data are placed in the repository: 



# Load packages -----------------------------------------------------------
library(Seurat)
library(SPATA2)
library(tidyverse)
library(igraph)
library(ggrepel)
library(ggforce)
library(ggraph)
library(SPATAwrappers)



# Load the add-on Functions
source("Functions.R")

# Load single cell reference data (GBmap)
seurat <- readRDS("path/to/GBMap")
seurat$annotation_level_4 <- 
  str_replace_all(seurat$annotation_level_4, "[ ]", "_") %>% 
  str_replace_all(., "-", "_") %>% 
  str_replace_all(., "/", "_")


# These are the samples used: "248_T" "259_T" "260_T" "269_T" "275_T" "296_T" "304_T" "313_T" "334_T"
# Please download the SPATA2 object in SPATAData

# SPTCR-seq Raw data are a rds file of all samples named by the ID ("248_T",...)
TCRs_list <- readRDS("path/to/RAW_SPTCR_seq_.rds")

# Run Analysis for a single sample (Preprocess and Figure 1): ---------------------------------------

# Load Data
TCRs <- TCRs_list[["275"]]
object <- readRDS("path/to/275_T.rds")


TCRs <- preprocessTCR(object, TCRs)
spot_TCR <- preprocessLD(TCRs)
clones <- normalizeTCR(object, spot_TCR)


### Plot UMI per Spot and Mean UMI per cell
Total_UMI <- 
  SPATA2::getFeatureDf(object) %>% 
  dplyr::select(barcodes,nCount_Spatial,Nr_of_cells) %>% 
  mutate(nCount_Spatial_per_cell=nCount_Spatial/Nr_of_cells) %>% 
  arrange(desc(nCount_Spatial))

ggplot(Total_UMI)+
  geom_point(mapping=aes(x=1:nrow(Total_UMI), y=nCount_Spatial))+
  geom_point(mapping=aes(x=1:nrow(Total_UMI), y=nCount_Spatial_per_cell))+
  theme_bw()


clone_mat <- normalizedTCRMatrix(object, spot_TCR)

## PCA Analysis of Clone expression / Spot
pca <- 
  Seurat::CreateSeuratObject(clone_mat) %>% 
  FindVariableFeatures() %>% 
  Seurat::ScaleData() %>% 
  Seurat::RunPCA()

ggplot(pca@reductions$pca@feature.loadings %>% as.data.frame())+
  geom_point(mapping=aes(x=PC_1, y=PC_2))+
  theme_classic()


#### Plot Clonal / Variance

ggplot(clones) + 
  geom_point(mapping=aes(x=1:nrow(clones), y=nr_spots/reads, size= reads), size=2)+
  theme_classic()

ggplot(clones) + 
  geom_point(mapping=aes(x=reads %>% log(), y=nr_spots,alpha=nr_spots, color= as.character(bert)), size=5)+
  geom_text_repel(data=clones %>% head(10), mapping=aes(x=reads %>% log(), y=nr_spots,label = CDR3_aa),size = 5)+
  geom_text_repel(data=clones %>% filter(nr_spots<100) %>% arrange(desc(reads)) %>% head(10), mapping=aes(x=reads %>% log(), y=nr_spots,label = CDR3_aa),color="darkgreen",size = 5)+
  xlab("Log UMI")+
  ylab("Number of Spots")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black"))



#### Measurement of the clone NN (if a clone is expressed in the adj. spot)
clone_NN <- getClonaladjacency(object, spot_TCR)
clone_NN_index <- 
  clone_NN %>% 
  group_by(CDR3_aa) %>% 
  summarise(C_NN=sum(NN))
clones <- 
  left_join(clones, clone_NN_index) %>% 
  mutate(C_NN_I =nr_spots/(C_NN+1))


ggplot(clones) + 
  geom_point(mapping=aes(x=reads %>% log(), y=C_NN_I,alpha=C_NN_I, color= as.character(bert)), size=5)+
  geom_text_repel(data=clones %>% arrange(desc(C_NN_I))%>% head(10), mapping=aes(x=reads %>% log(), y=C_NN_I, label = CDR3_aa),color="darkgreen",size = 5)+
  xlab("Log UMI")+
  ylab("Neighborhood Index")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black"))

### Evaluation of a clone shows local clonal expansion or is distributed across the sample by merging the
### clonal / non clonal index and the clonal NN estimation as a  - spatial clonality index -

# Merge the clonal CDR3 and diverse CDR3 to clonal / non clonal index
clonal_index <- clones %>% mutate(clonal=ifelse(log(reads)>-1 & C_NN_I>5, 1, 0))


Clonality_Index <- 
  clone_NN %>% 
  left_join(., clonal_index %>% dplyr::select(CDR3_aa,clonal), by="CDR3_aa") %>% 
  group_by(barcodes) %>% 
  summarise(Clonality_Index=sum(clonal),
            TCR_UMI = sum(expression),
            number_clones = length(CDR3_aa))

Clonality_Index <- left_join(data.frame(barcodes=SPATA2::getBarcodes(object)[[1]]), Clonality_Index)
Clonality_Index[is.na(Clonality_Index)] <- 0



object <- addFeatures(object, Clonality_Index, overwrite = T)
col=colorRampPalette(RColorBrewer::brewer.pal(9, "Greys") )
p <- plotSurface(object, color_by = "TCR_UMI", display_image = F)
add_hull(object,p,pt.b = 7)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,100),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()

p <- plotSurface(object, color_by = "Clonality_Index", display_image = F)
add_hull(object,p,pt.b = 7)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,10),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()

p <- plotSurface(object, color_by = "number_clones", display_image = F)
add_hull(object,p,pt.b = 7)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,30),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()


### Isolate spots with high clonal expansion, high TCR variance or low TCR abundance
Clonality_Index <- mutate(Clonality_Index,TCR_prams=case_when(
  number_clones>15 ~ "Diverse",
  Clonality_Index>2 & number_clones<15 ~ "Clonal",
  TCR_UMI==0 ~ "NoTcell",
  TRUE ~ "undetermined"
)) 
object <- addFeatures(object, Clonality_Index, overwrite = T)
plotSurface(object, color_by = "TCR_prams", display_image = F)















# Figure 2-3 ----------------------------------------------------------------

SPTCR_Index_all <- readRDS("path/to/Processed_SPTCR_seq_Histology.rds")
df <- SPTCR_Index_all %>% group_by(sample, Histology_IVY) %>% summarise(mean=mean(TCR_UMI) ) 

level <-  SPTCR_Index_all %>% group_by(sample) %>% summarise(mean=mean(TCR_UMI) ) %>% arrange((mean)) %>% pull(sample)


ggplot(data = df, aes(x = factor(sample, levels = level), y = mean, fill = Histology_IVY)) + 
  geom_bar(position="dodge",stat="identity") + 
  theme_classic()+
  scale_fill_manual(values=colors_ivy)

SPATAImmune::plotBarCompare(SPTCR_Index_all, "Histology_IVY", "TCR_prams")+scale_fill_manual(values=colors_ivy)
SPTCR_Index_all[SPTCR_Index_all$Clonality_Index>0, ]$TCR_prams <- "Clonal"
SPATAImmune::plotBarCompare(SPTCR_Index_all, "TCR_prams", "sample")+scale_fill_manual(values=color_class)


# Figure 4 ----------------------------------------------------------------
## Figure 4b

allSample_TCR_list <- readRDS("path/to/Processed_SPTCR_seq.rds")

clones <- map_dfr(1:9, ~allSample_TCR_list[[.x]][[2]] %>% mutate(sample=samples_all[.x]))
color_samples <- RColorBrewer::brewer.pal(9, "Set3") %>% scales::muted(l = 70, c = 60)
names(color_samples) <- samples_all

ggplot(clones) + 
  geom_point(mapping=aes(x=reads %>% log(), 
                         y=nr_spots,
                         alpha=0.8, 
                         color= as.character(sample)), size=4)+
  #geom_text_repel(data=clones %>% head(10), mapping=aes(x=reads %>% log(), y=nr_spots,label = CDR3_aa),size = 5)+
  #geom_text_repel(data=clones %>% filter(nr_spots<100) %>% arrange(desc(reads)) %>% head(10), mapping=aes(x=reads %>% log(), y=nr_spots,label = CDR3_aa),color="darkgreen",size = 5)+
  xlab("Log UMI")+
  ylab("Number of Spots")+
  ylim(0,150)+
  xlim(-5,5)+
  scale_color_manual(values = color_samples)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black"))

## Figure 4c
ggplot(clones) + 
  geom_point(mapping=aes(x=reads, y=jitter(C_NN_I, 500), alpha=0.8, color= as.character(sample)), size=2)+
  #geom_text_repel(data=clones %>% arrange(desc(C_NN_I))%>% head(10), mapping=aes(x=reads %>% log(), y=C_NN_I, label = CDR3_aa),color="darkgreen",size = 5)+
  xlab("Log UMI")+
  ylab("Neighborhood Index")+
  theme_bw() +
  xlim(0,10)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black"))


### Representative Examples of the sample: 259_T


object <- readRDS("path/to/259_T")
Clonality_Index <- 
  allSample_TCR_list[[2]]$Clonality_Index %>% 
  left_join(data.frame(barcodes=SPATA2::getBarcodes(object)[[1]]), Clonality_Index)
Clonality_Index[is.na(Clonality_Index)] <- 0
diversity <- getSpotwiseDiversity(object, allSample_TCR_list[[2]]$spot_TCR)
Clonality_Index <- left_join(Clonality_Index, diversity)
Clonality_Index[is.na(Clonality_Index$Diversity_index), "Diversity_index"] <- 0
object <- addFeatures(object, Clonality_Index, overwrite = T)

# Add Histology
object <- runHistologyAnnotation(object, geneSet)


## H&E Figure 4 d-e
## Surface plots
col=colorRampPalette(RColorBrewer::brewer.pal(9, "Greys") )
p <- 
  plotSurface(object, color_by = "TCR_UMI", display_image = F, pt_alpha = 0)
add_hull(object,p,pt.b = 7)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,100),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()


## Histo
colors_ivy <- c(colors_ivy, IVY_Necrosis="#000000")
p <- 
  plotSurface(object, color_by = "Histology_IVY",  pt_alpha = 0.8, display_image = F)
add_hull(object,p, pt.b = 7)+
  scale_color_manual(values = colors_ivy)+coord_fixed()


df <- SPATA2::joinWithFeatures(list_obj_259[[2]], features = "Cell_lables")
colors <- readRDS("~/colors_cell_deconv.RDS")
cc <- colors$colors ;names(cc) <- colors$annotation_level_4
lymphoid <- colors[colors$annotation_level_2=="Myeloid", ]$annotation_level_4
df$alpha=0.5
df[df$Cell_lables %in% lymphoid, ]$alpha=1
ggplot(df, aes(x, y, group = -1L)) +
  geom_voronoi_tile(aes(fill = Cell_lables, alpha=alpha), max.radius = 10)+
  scale_fill_manual(values=cc)+
  theme_void()+
  coord_fixed()+
  Seurat::NoLegend()


col=colorRampPalette(RColorBrewer::brewer.pal(9, "Greens") )
p <- plotSurface(object, color_by = "TCR_UMI", display_image = F)
add_hull(object,p,pt.b = 7)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,100),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()


p <- plotSurface(object, color_by = "number_clones", display_image = F)
add_hull(object,p,pt.b = 7)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,80),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()


p <- plotSurface(object, color_by = "Clonality_Index", display_image = F)
add_hull(object,p,pt.b = 7)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,5),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()

p <- plotSurface(object, color_by = "Diversity_index", display_image = F)
add_hull(object,p,pt.b = 7)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(3,10),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()


color_class <- c(RColorBrewer::brewer.pal(3, "Set1") %>% scales::muted(l = 60, c = 70), "#E5E5E5")
names(color_class) <- c("Clonal", "Diverse", "NoTcell", "undetermined")

p=plotSurface(object, color_by = "TCR_prams",  pt_alpha = 1, display_image = F)
add_hull(object,p, pt.b = 7)+
  scale_color_manual(values = color_class)+coord_fixed()


# Figure 5 Histology and Barplots ----------------------------------------------------------------
Index <- map_dfr(1:9, ~allSample_TCR_list[[.x]][[3]] %>% mutate(sample=samples_all[.x]))
Index <- mutate(Index,TCR_prams=case_when(
  number_clones>5 ~ "Diverse",
  Clonality_Index>0  ~ "Clonal",
  TCR_UMI==0 ~ "NoTcell",
  TRUE ~ "undetermined"
))


## Figure 5b
SPATAImmune::plotBarCompare(Index, "TCR_prams", "sample")

## Figure 5c
filter_index <- Index %>% filter(TCR_prams=="Clonal") 
barplot(table(filter_index$TCR_prams, filter_index$sample))





# Figure 5 Spatial visualization of clones --------------------------------


object <- readRDS("259_T.rds")
SPTCR_Index <- SPTCR_Index_all[[2]]

# Clone level Data
clones <- map_dfr(1:9, ~allSample_TCR_list[[.x]][[2]] %>% mutate(sample=samples_all[.x]))
clones_select <- clones %>% filter(sample == sample_use)


target_clones <- clones_select %>% arrange(desc(nr_spots)) %>% head(10) %>% pull(CDR3_aa) %>% sample(10)

clones_plot <- 
  allSample_TCR_list[[sample_use_i]]$spot_TCR %>% 
  filter(CDR3_aa %in% target_clones) %>% 
  group_by(barcodes, CDR3_aa) %>% 
  summarise(expression=mean(expression))


clones_plot <- left_join(clones_plot, SPATA2::getCoordsDf(object), by="barcodes")
clones_plot$x=jitter(clones_plot$x, factor = 100);clones_plot$y=jitter(clones_plot$y, factor = 100)
clones_plot <- clones_plot %>% filter(expression>3)


background <- SPATA2::getCoordsDf(object)
col_clones <- RColorBrewer::brewer.pal(10, "Set3")

ggplot()+
  scattermore::geom_scattermore(data=background, mapping = aes(x,y),pointsize = 6, color="black")+
  scattermore::geom_scattermore(data=background, mapping = aes(x,y),pointsize = 5, color="white")+
  geom_point(data=clones_plot, mapping = aes(x,y,color=CDR3_aa,size=expression))+
  scale_color_manual(values=col_clones)+
  coord_fixed()+
  SPATA2::ggpLayerAxesSI(object)


prod <- c("PRF1", "GZMB", "IFNG", "FASLG", "TNF")
exh <- c("EOMES", "TOX", "TIGIT",  "HAVCR2", "LAG3", "CTLA4", "PDCD1")

object@information$active_expr_mtr <- "denoised" 
object@information$active_mtr <- "denoised" 
score <- joinWith(object, genes = prod, average_genes=T) %>% select(barcodes, mean_genes)
names(score)[2] <- "Cyto"
object <- addFeatures(object, score, overwrite = T)


score <- joinWith(object, genes = exh, average_genes=T) %>% select(barcodes, mean_genes)
names(score)[2] <- "exh"
object <- addFeatures(object, score, overwrite = T)

p <- plotSurface(object, color_by = "exh", display_image = F, pt_alpha = 1)
add_hull(object,p, pt.b = 7)+
  coord_fixed()+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0.3,0.7),
                         oob = scales::squish, na.value = "white")

p <- plotSurface(object, color_by = "Cyto", display_image = F, pt_alpha = 1)
add_hull(object,p, pt.b = 7)+
  coord_fixed()+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0.2,1),
                         oob = scales::squish, na.value = "white")



# Spatial cell-cell interaction model (Figure 6) ------------------------------

### Extract the T cells from this regions (single cell level) and perform single cell trajectory analysis
### Characterize the cell NN of clonal expansion, high TCR variance or low TCR abundance spots


seurat <- readRDS("~/path/to/Seurat_GBMap.rds")
nn_clonal <- readRDS("nn_clonal.rds")
nn_diverse <- readRDS("nn_diverse.rds")


#Graph from clones:
nn <- nn_clonal %>% group_by(from, to) %>% count()
#Graph from T cell diversity:
nn <- map_dfr(1:8, .f=function(i){
  sample=samples_all[i]
  out <- nn_diverse %>% filter(sample==sample)
  out <- out[runif(700, 1, nrow(out)) %>% round(), ]
})
nn <- nn %>% group_by(from, to) %>% count()

## Figure 6b

#Create Graph
size <- nn %>% group_by(from) %>% count
Lymph <- colors[colors$annotation_level_2=="Lymphoid",]$annotation_level_4
nn$n[nn$from %in% Lymph]=nn$n[nn$from %in% Lymph]*5

#Filter
nn <- nn %>% filter(n>5)
size$n[size$from %in% Lymph]=size$n[size$from %in% Lymph]*5
size <- left_join(data.frame(from=c(nn$from, nn$to) %>% unique()), size)
size[is.na(size)]=0

# igraph object
net <- igraph::graph_from_data_frame(d=nn[,1:2], vertices=c(nn$from, nn$to) %>% unique(), directed=T)
E(net)$weight <- nn$n
V(net)$type <- c(nn$from, nn$to) %>% unique()
V(net)$size <- size$n
g <- delete.vertices(simplify(net), degree(net)==0)


ggraph(g, 'igraph',  algorithm = 'fr')+
  geom_edge_link0(aes(edge_alpha = weight, width=weight), color="lightgrey") +
  geom_node_point(aes(color=type, size=size))+
  theme_classic()+
  geom_text_repel(aes(x=x, y=y,label=type %>% str_replace_all(., "_", "-")), size=3)+
  scale_color_manual(values=col_n)+
  coord_fixed()+
  Seurat::NoLegend()


ggraph(g, 'tree')+
  geom_edge_link(color="lightgrey", alpha=0.5, width=1)+
  geom_node_point(aes(color=type, size=size))+
  geom_text_repel(aes(x=x, y=y,label=type %>% str_replace_all(., "_", "-")), size=3)


graph <- tbl_graph(igraph::as_data_frame(g, "vertices"), igraph::as_data_frame(g, "edges"))


set.seed(1)
p <- ggraph(graph, 'circlepack',weight = size) + 
  geom_edge_link(color="lightgrey", alpha=0.5, width=1) + 
  geom_node_point(aes(colour = type),  size=5) +
  theme_classic()+
  scale_color_manual(values=col_n)+
  #coord_fixed()+
  Seurat::NoLegend()

p+geom_text_repel(data=p$data %>% filter(type %in% Lymph), 
                  aes(x=x, y=y,label=type %>% str_replace_all(., "_", "")), size=3)


## Figure 6 Plots of the NCEM

#load files

path= "/path/to/NCEM/"  
files <- path %>% dir()
files_ex <-files %>% str_detect(., pattern = "Exhausted") 


load.data.ex <- map(files[files_ex], function(i){
  
  print(i)
  data <- read.csv(paste0(path, i))
  
  #covert into graph
  g <- data.frame(from=data$sender,
                  to=data$receiver,
                  magnitude=data$magnitude,
                  de = data$de_genes)
  
  #g <- g %>% filter(magnitude>20)
  
  library(igraph)
  
  g$from <- g$from %>% str_replace_all(., " ", "_") %>% str_replace_all(., "-", "_") %>% toupper()
  g$to <- g$to %>% str_replace_all(., " ", "_") %>% str_replace_all(., "-", "_") %>% toupper()
  
  graph <- igraph::graph_from_data_frame(g[1:2], directed=T)
  E(graph)$magnitude <- g$magnitude
  E(graph)$de <- g$de
  E(graph)$color <- g$from
  type <- c(g$from, g$to) %>% unique()
  
  names(col_n)[4] <- "MES_LIKE_HYPOXIA"
  names(col_n)[3] <- "MES_LIKE_HYPOXIA_INDEPENDENT"
  names(col_n)[names(col_n)=="TAM_BDM_HYPOXIA_MES"] <- "TAM_BDM_HYPOXIA"
  
  V(graph)$type <- type 
  #graph <- delete.vertices(simplify(graph), degree(graph)==0)
  
  names(col_n) <- names(col_n) %>% str_replace_all(., " ", "_") %>% str_replace_all(., "-", "_") %>% toupper()
  
  library(tidygraph)  
  library(ggraph)
  library(ggrepel)
  
  mat <- igraph::as_adjacency_matrix(graph, attr="magnitude") %>% as.matrix()
  
  return(list(g, mat))
  
})


g <- map_dfr(1:3, ~load.data.ex[[.x]][[1]] %>% mutate(sample=.x))
g <- g %>% group_by(from, to) %>% summarise(magnitude=mean(magnitude), de=mean(de))
g <- g %>% filter(magnitude>15)

library(igraph)
graph <- igraph::graph_from_data_frame(g[1:2], directed=T)
E(graph)$magnitude <- g$magnitude
E(graph)$de <- g$de
E(graph)$color <- g$from
type <- c(g$from, g$to) %>% unique()
V(graph)$type <- type 

V(graph)$type
order=c("CDC2", "PDC", "B_CELL", "PLASMA_B", "NK","CD8_NK_SIG", "CD8_CYTOTOXIC","CD4_INF","REG_T",
        "MES_LIKE_HYPOXIA", "MES_LIKE_HYPOXIA_INDEPENDENT","AC_LIKE", "NPC_LIKE_NEURAL","NPC_LIKE_OPC" ,
        "TAM_BDM_MHC","TAM_BDM_HYPOXIA", 
        "PERIVASCULAR_FIBROBLAST", "SMC" )

ggraph(graph, layout = 'linear', circular = TRUE, sort.by = factor(type, levels = order))+
  geom_edge_arc2(aes(width=magnitude, alpha=magnitude, edge_color = color))+
  geom_node_point(aes(color=type), size=8)+
  coord_fixed()+
  theme_classic()+
  scale_color_manual(values=col_n)+
  scale_edge_color_manual(values=col_n)+
  geom_text_repel(aes(x=x, y=y,label=type %>% str_replace_all(., "_", "-")), size=3)


g <- map_dfr(1:3, ~load.data.ex[[.x]][[1]] %>% mutate(sample=.x))
g <- g %>% group_by(from, to) %>% summarise(magnitude=mean(magnitude), de=mean(de))

library(igraph)
graph <- igraph::graph_from_data_frame(g[1:2], directed=T)
E(graph)$magnitude <- g$magnitude
E(graph)$de <- g$de
E(graph)$color <- g$from
type <- c(g$from, g$to) %>% unique()
V(graph)$type <- type 

V(graph)$type
order=c("CDC2", "PDC", "B_CELL", "PLASMA_B", "NK","CD8_NK_SIG", "CD8_CYTOTOXIC","CD4_INF","REG_T",
        "MES_LIKE_HYPOXIA", "MES_LIKE_HYPOXIA_INDEPENDENT","AC_LIKE", "NPC_LIKE_NEURAL","NPC_LIKE_OPC" ,
        "TAM_BDM_MHC","TAM_BDM_HYPOXIA", 
        "PERIVASCULAR_FIBROBLAST", "SMC" )


mat <- igraph::as_adjacency_matrix(graph, attr="magnitude") %>% as.matrix()

pheatmap::pheatmap(mat[order, order], 
                   cluster_rows = F, 
                   cluster_cols = F, 
                   col = colorRampPalette(c("#FFFFFF",RColorBrewer::brewer.pal(9, "Reds")) )(50) ) 



load.data.ex <- map(files[!files_ex], function(i){
  
  print(i)
  data <- read.csv(paste0(path, i))
  
  #covert into graph
  g <- data.frame(from=data$sender,
                  to=data$receiver,
                  magnitude=data$magnitude,
                  de = data$de_genes)
  
  #g <- g %>% filter(magnitude>20)
  
  library(igraph)
  
  g$from <- g$from %>% str_replace_all(., " ", "_") %>% str_replace_all(., "-", "_") %>% toupper()
  g$to <- g$to %>% str_replace_all(., " ", "_") %>% str_replace_all(., "-", "_") %>% toupper()
  
  graph <- igraph::graph_from_data_frame(g[1:2], directed=T)
  E(graph)$magnitude <- g$magnitude
  E(graph)$de <- g$de
  E(graph)$color <- g$from
  type <- c(g$from, g$to) %>% unique()
  
  names(col_n)[4] <- "MES_LIKE_HYPOXIA"
  names(col_n)[3] <- "MES_LIKE_HYPOXIA_INDEPENDENT"
  names(col_n)[names(col_n)=="TAM_BDM_HYPOXIA_MES"] <- "TAM_BDM_HYPOXIA"
  
  V(graph)$type <- type 
  #graph <- delete.vertices(simplify(graph), degree(graph)==0)
  
  names(col_n) <- names(col_n) %>% str_replace_all(., " ", "_") %>% str_replace_all(., "-", "_") %>% toupper()
  
  library(tidygraph)  
  library(ggraph)
  library(ggrepel)
  
  mat <- igraph::as_adjacency_matrix(graph, attr="magnitude") %>% as.matrix()
  
  return(list(g, mat))
  
})
g <- map_dfr(1:7, ~load.data.ex[[.x]][[1]] %>% mutate(sample=.x))
g <- g %>% group_by(from, to) %>% summarise(magnitude=mean(magnitude), de=mean(de))
g <- g %>% filter(magnitude>15)
library(igraph)
graph <- igraph::graph_from_data_frame(g[1:2], directed=T)
E(graph)$magnitude <- g$magnitude
E(graph)$de <- g$de
E(graph)$color <- g$from
type <- c(g$from, g$to) %>% unique()
V(graph)$type <- type 

V(graph)$type
order=c("CDC2", "PDC", "B_CELL", "PLASMA_B", "NK","CD8_NK_SIG", "CD8_CYTOTOXIC","CD4_INF","REG_T",
        "MES_LIKE_HYPOXIA", "MES_LIKE_HYPOXIA_INDEPENDENT","AC_LIKE", "NPC_LIKE_NEURAL","NPC_LIKE_OPC" ,
        "TAM_BDM_MHC","TAM_BDM_HYPOXIA", 
        "PERIVASCULAR_FIBROBLAST", "SMC" )

ggraph(graph, layout = 'linear', circular = TRUE, sort.by = factor(type, levels = order))+
  geom_edge_arc2(aes(width=magnitude, alpha=magnitude, edge_color = color))+
  geom_node_point(aes(color=type), size=8)+
  coord_fixed()+
  theme_classic()+
  scale_color_manual(values=col_n)+
  scale_edge_color_manual(values=col_n)+
  geom_text_repel(aes(x=x, y=y,label=type %>% str_replace_all(., "_", "-")), size=3)


g <- map_dfr(1:3, ~load.data.ex[[.x]][[1]] %>% mutate(sample=.x))
g <- g %>% group_by(from, to) %>% summarise(magnitude=mean(magnitude), de=mean(de))

library(igraph)


graph <- igraph::graph_from_data_frame(g[1:2], directed=T)
E(graph)$magnitude <- g$magnitude
E(graph)$de <- g$de
E(graph)$color <- g$from
type <- c(g$from, g$to) %>% unique()
V(graph)$type <- type 

V(graph)$type
order=c("CDC2", "PDC","DC2", "B_CELL", "PLASMA_B", "NK","CD8_NK_SIG", "CD8_CYTOTOXIC","CD4_INF","REG_T",
        "MES_LIKE_HYPOXIA", "MES_LIKE_HYPOXIA_INDEPENDENT","AC_LIKE", "NPC_LIKE_NEURAL","NPC_LIKE_OPC" ,
        "TAM_BDM_MHC","TAM_BDM_HYPOXIA","TAM_BDM_INF", 
        "PERIVASCULAR_FIBROBLAST", "SMC" )


mat <- igraph::as_adjacency_matrix(graph, attr="magnitude") %>% as.matrix()


pheatmap::pheatmap(mat[order[order %in% rownames(mat)], order[order %in% colnames(mat)] ], 
                   cluster_rows = F, 
                   cluster_cols = F, 
                   col = colorRampPalette(c("#FFFFFF",RColorBrewer::brewer.pal(9, "Reds")) )(50) ) 






## Figure 7:

# Cell-Cell communication  ------------------------------------------------

#Cell Chat library
library(CellChat)
library(patchwork)

CellChatDB <- CellChatDB.human 
DB <- CellChatDB$interaction
pathways <- DB$interaction_name_2
inter <- map_dfr(.x= 1:length(pathways), .f=function(i){
  
  pw <- pathways[i]
  pw_l1 <- pw %>% str_split(., " - ") %>% unlist()
  Ligand <- pw_l1[1]
  Receptor <- pw_l1[2] %>% str_remove_all(., "[)]") %>% str_remove_all(., "[(]") %>% str_split(., "[+]") %>% unlist()
  
  out <- data.frame(R = Receptor, L = Ligand) 
  
},.progress=T)


#B-cell on NK
object <- getCelltypeSpecificGeneExpression(object, 
                                            celltype = "C2L_B.Cell",
                                            mtr_name="B_cell",
                                            active_mtr="denoised", 
                                            enhance_factor=10)
object %>% getFeatureNames()
object <- getCelltypeSpecificGeneExpression(object, 
                                            celltype = "C2L_Nk",
                                            mtr_name="NK",
                                            active_mtr="denoised", 
                                            enhance_factor=10)

object %>% getFeatureNames()
object <- getCelltypeSpecificGeneExpression(object, 
                                            celltype = "C2L_Cd4.Inf",
                                            mtr_name="CD4_IFN",
                                            active_mtr="denoised", 
                                            enhance_factor=10)

#subset inter

inter <- 
  inter %>% 
  filter(R %in% SPATA2::getGenes(object)) %>% 
  filter(L %in% SPATA2::getGenes(object))


#B cell sender genes
active_mtr <- "B_cell"
object@information$active_mtr=active_mtr
object@information$active_expr_mtr= active_mtr
B_mat <- getExpressionMatrix(object, mrt_name=active_mtr)


#adapt
ligand_unique <- inter$L %>% unique()
B_mat <- B_mat[ligand_unique, ]
B_mat %>% dim()


#NK receiver genes
active_mtr <- "NK"
object@information$active_mtr=active_mtr
object@information$active_expr_mtr= active_mtr
NK_mat <- getExpressionMatrix(object, mrt_name=active_mtr)

#adapt
receptor_unique <- inter$R %>% unique()
NK_mat <- NK_mat[receptor_unique, ]
NK_mat %>% dim()


# Only on spot level
eval <- map_dfr(1:nrow(inter), .f=function(i){
  
  signal <- inter[i, ]
  
  R <- NK_mat[signal$R, ]
  L <- B_mat[signal$L,] 
  
  cor <- cor(R,L)
  KLD <- LaplacesDemon::KLD(R,L)$mean.sum.KLD
  cos <- lsa::cosine(R,L)
  
  data.frame(R=signal$R, L=signal$L, spatial_weight=cor, coesine=cos, KLD=KLD)
  
  
}, .progress=T)
plot_genes <- eval %>% arrange(KLD) %>% head(100)

## plot res


col=colorRampPalette(c("#FFFFFF",RColorBrewer::brewer.pal(9, "Greys") ))
object@information$active_mtr="B_cell"
p1 <- plotSurface(object, color_by = plot_genes$L[1], display_image = F, pt_alpha = 1)
p1 <- add_hull(object,p1, pt.b = 7)+
  coord_fixed()+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,0.1),
                         oob = scales::squish, na.value = "white")
object@information$active_mtr="NK"
p2 <- plotSurface(object, color_by = plot_genes$R[1], display_image = F, pt_alpha = 1)
p2 <- add_hull(object,p2, pt.b = 7)+
  coord_fixed()+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,0.8),
                         oob = scales::squish, na.value = "white")

p1+p2

object@information$active_mtr="B_cell"
a <- joinWith(object, genes=plot_genes$L[1])
object@information$active_mtr="NK"
b <-joinWith(object, genes=plot_genes$R[1])
a <- mutate(a, b=b$SDC4) %>% mutate(c=COL6A1+b)


## Quantify for regions

plot_genes <- eval %>% arrange(KLD) %>% head(20)
quantify <- map_dfr(1:nrow(plot_genes), .f=function(i){
  
  object@information$active_mtr="B_cell"
  a <- joinWith(object, genes=plot_genes$L[i])
  object@information$active_mtr="NK"
  b <-joinWith(object, genes=plot_genes$R[i])
  
  df <- data.frame(barcodes=a$barcodes, sum=a %>% pull(!!sym(plot_genes$L[i]))+b %>% pull(!!sym(plot_genes$R[i])) )
  df <- left_join(df, joinWith(object, features="Tcell_type") %>% select(barcodes, Tcell_type))
  
  df %>% group_by(Tcell_type) %>% summarise(val=mean(sum)) %>% mutate(sig=paste0(plot_genes$R[i], "_",plot_genes$L[i]))
  
  
  
  
}, .progress=T)

mat <- quantify %>% reshape2::acast(Tcell_type~sig, value.var = "val", fun.aggregate = mean)

quantify <- quantify %>% arrange(desc(val)) %>% mutate(sig=factor(sig, level=unique(sig) ))

ggplot(quantify, mapping = aes(y=sig, x=val, fill=Tcell_type))+
  geom_bar(position = "dodge",stat="identity")+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(3, "Set1")[1:2], "white"))+
  theme_classic()


plotColorOverlap_RL(object,
                    R_mat="CD4_IFN",
                    L_mat = "NK",
                    L_enh=0,
                    two.colors = c("red", "green"),
                    negative.color="white",
                    col.threshold = 0.1,
                    R=plot_genes$R[i], L=plot_genes$L[i])+
  SPATA2::ggpLayerAxesSI(object)




#Graph
plot_genes <- eval %>% arrange(KLD) %>% head(10)
graph <- igraph::graph_from_data_frame(plot_genes[1:2], directed=T)
E(graph)$coesine <- plot_genes$coesine
E(graph)$KLD <- plot_genes$KLD
E(graph)$spatial_weight <- plot_genes$spatial_weight
V(graph)$type <- c(plot_genes$R, plot_genes$L) %>% unique()


df <- data.frame(genes= c(plot_genes$R, plot_genes$L), 
                 class=c(rep("R", nrow(plot_genes)), rep("L", nrow(plot_genes))))
df <- df[!duplicated(df$genes), ]
V(graph)$class <- df$class


library(ggraph)
library(ggrepel)
ggraph(graph, layout = 'linear', circular = TRUE)+
  geom_edge_arc2(aes(width=coesine, alpha=coesine), color="lightgrey")+
  geom_node_point(aes(col=class), size=9)+
  coord_fixed()+
  theme_classic()+
  geom_text_repel(aes(x=x, y=y,label=type %>% str_replace_all(., "_", "-")), size=5)




# Cell Type specific expression --------------------------------------------


## Add cell2Location
object <- readRDS("313_T") ## Make sure that the object contains the Cell2Location or other deconvolution in the fdata slot


# 1. Get weighted cell type specific gene signature

DownScaleSeurat <- function (seurat, maintain_var = "cell_states", max = 10000, 
                             min = 50, only_var = T, factor = 5, n_feature = 3000) {
  tab_quant <- seurat@meta.data %>% as.data.frame() %>% pull(!!sym(maintain_var)) %>% 
    table() %>% as.data.frame() %>% mutate(nr = Freq/min(Freq)) %>% 
    mutate(cells = round(nr * factor)) %>% mutate(cells = ifelse(cells > 
                                                                   max, max, cells)) %>% mutate(cells = ifelse(cells < min, 
                                                                                                               Freq, cells)) %>% dplyr::rename(`:=`("cluster", .)) %>% 
    dplyr::select(cluster, cells) %>% mutate(cluster = as.character(cluster))
  cells <- map(.x = 1:nrow(tab_quant), .f = function(i) {
    sample(rownames(seurat@meta.data[seurat@meta.data[, maintain_var] == 
                                       tab_quant$cluster[i], ]), tab_quant$cells[i])
  }) %>% unlist()
  seurat.out <- subset(seurat, cells = cells)
  
  seurat.out <- Seurat::FindVariableFeatures(seurat.out, return.only.var.genes = only_var, nfeatures = n_feature) %>% 
    Seurat::ScaleData()
  
  return(seurat.out)
}
down_seurat <- DownScaleSeurat(seurat, maintain_var="annotation_level_4", max = 5000)
Seurat::Idents(down_seurat) <- "annotation_level_4"
diff <- FindAllMarkers(down_seurat, logfc.threshold=0.5, only.pos=T)


# 2. Define cell type specific suppressor and enhancer

## Cell type naming
translate <- data.frame(
  C2L=names(Cell2loc)[2:length(names(Cell2loc))][!names(Cell2loc)[2:length(names(Cell2loc))]=="C2L_Neuron"] ,
  annotation_level_4=colors$annotation_level_4[order(colors$annotation_level_4)]
)

diff <- diff[diff$cluster %in% translate$annotation_level_4, ]
translate <- translate[translate$annotation_level_4 %in% diff$cluster, ]
diff$cluster <- map(diff$cluster, ~translate[translate$annotation_level_4==.x, "C2L"]) %>% unlist()

getCelltypeSpecificGeneExpression <- function(object, 
                                              celltype, 
                                              enhance_factor=10, 
                                              mtr_name="T_cell", 
                                              active_mtr="scaled"){
  
  if(length(celltype)==1){
    message(paste0("Create expr Matrix from a single cell type"))
    
    sig <- 
      diff %>% 
      filter(cluster==celltype) %>% 
      mutate(enhancer=avg_log2FC+1) %>% 
      dplyr::select(gene, enhancer) 
    
    object@information$active_mtr=active_mtr
    object@information$active_expr_mtr= active_mtr
    
    mat <- SPATA2::getExpressionMatrix(object)
    include=rownames(mat)[rownames(mat) %in% sig$gene]
    sig <- sig %>% filter(gene %in% include)
    
    message(paste0("Run cell type expression enhancer"))
    factor <-   abs(t(mat[sig$gene, ] )) %*% diag(sig$enhancer)
    mat[sig$gene, ] <- mat[sig$gene, ]+t(factor)
    
    
    message(paste0("Run spot supressor"))
    sup <- joinWith(object, features = celltype) %>% mutate(f=exp(!!sym(celltype)))
    
    f <- sup %>% pull(f)
    f <- f+c(c(f-1)*enhance_factor)
    
    factor <- t(t(abs(mat)) * f)
    mat <- mat+factor
    
    message(paste0("Return matrix under the name: ", mtr_name ))
    object <- addExpressionMatrix(object, expr_mtr = mat,  mtr_name=mtr_name)
    
  }else{
    message(paste0("Create expr Matrix from a multiple cell types"))
    
    sig <- 
      diff %>% 
      filter(cluster %in% celltype) %>% 
      mutate(enhancer=avg_log2FC+1) %>% 
      dplyr::select(gene, enhancer) %>% 
      group_by(gene) %>% 
      summarize(enhancer=max(enhancer)) %>% 
      as.data.frame()
    
    object@information$active_mtr=active_mtr
    object@information$active_expr_mtr= active_mtr
    mat <- SPATA2::getExpressionMatrix(object)
    include=rownames(mat)[rownames(mat) %in% sig$gene]
    sig <- sig %>% filter(gene %in% include)
    
    
    message(paste0("Run cell type expression enhancer"))
    factor <-   abs(t(mat[sig$gene, ] )) %*% diag(sig$enhancer)
    mat[sig$gene, ] <- mat[sig$gene, ]+t(factor)
    
    
    message(paste0("Run spot supressor"))
    sup <- joinWith(object, features = celltype) %>% as.data.frame()
    sup_x <- sup[, celltype] %>% rowSums()/length(celltype)
    sup <- mutate(sup, f=exp(sup_x))
    
    f <- sup %>% pull(f)
    f <- f+c(c(f-1)*enhance_factor)
    
    factor <- t(t(abs(mat)) * f)
    mat <- mat+factor
    
    message(paste0("Return matrix under the name: ", mtr_name ))
    object <- addExpressionMatrix(object, expr_mtr = mat,  mtr_name=mtr_name)
    
    
  }
  
  return(object)
}



celltype <- c("C2L_Cd8.Cytotoxic", "C2L_Cd8.Em")

object <- getCelltypeSpecificGeneExpression(object, 
                                            celltype = celltype,
                                            active_mtr="denoised", 
                                            enhance_factor=10)

object@information$active_mtr="scaled"
object@information$active_mtr="T_cell"

## Plot 
col=colorRampPalette(c("#FFFFFF",RColorBrewer::brewer.pal(9, "Greys") ))
genes <- c("CD8A", "CD3D", "GZMK","GZMB", "PDCD1", "HAVCR2")
ggpubr::ggarrange(plotlist=map(genes, .f=function(g){
  p1 <- plotSurface(object, color_by = g, display_image = F, pt_alpha = 1)
  p1 <- add_hull(object,p1, pt.b = 7)+
    coord_fixed()+
    scale_colour_gradientn(colours = col(50), 
                           limits= c(0,0.3),
                           oob = scales::squish, na.value = "white")
}), ncol = 3, nrow = 2)




## Figure 8:

# Metabolic Data ----------------------------------------------------------



list_matabolomic <-  readRDS("ListMetabolicData")# This list of files can be created based on the repository: https://osf.io/8qbdz/

for(s in intersect(samples_all, list_matabolomic$sample)){
  
  
  sample_use <- s
  sample_use_i <- which(samples_all==sample_use)
  ## preprocessing
  list_data <- readRDS("~/Desktop/NIO/ALA/Metabolomic/list_data.RDS")
  object <- readRDS(list_data[list_data$samples==sample_use, ]$dir) %>% SPATA2::updateSpataObject()
  #SPATA2::getExpressionMatrix(object, "MALDI") %>% rownames()
  your_path_to_anno <- "~/Desktop/SpatialTranscriptomics/Visium/Visium/MALDI_RAW/Images and metaspace annotation/Metaspace_annotations_943peaks.csv"
  metabolic_pw <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/MALDI_RAW/Images and metaspace annotation/meta_pw.RDS")
  
  # We first we annotate each moleculeIds to the corresponding m/z 
  annotations <- read.csv(your_path_to_anno, sep=";")
  anno_df <- 
    map_dfr(1:nrow(annotations),  .f=function(i){data.frame(molID=
                                                              str_split(annotations$moleculeIds[i], pattern=", ") %>% 
                                                              unlist(),
                                                            mz=annotations$mz[i])}, .progress=T) %>% 
    mutate(match=MALDIquant::match.closest(anno_df$mz, as.numeric(rownames(SPATA2::getExpressionMatrix(object)))) )
  
  
  meta_mx <- SPATA2::getExpressionMatrix(object)[anno_df$match, ]
  meta_mx <- 
    cbind(anno_df, meta_mx %>% as.data.frame()) %>% 
    dplyr::select(-mz, -match) %>% 
    dplyr::group_by(molID) %>% 
    summarise_all(.,.funs = mean)
  row_n <- meta_mx$molID
  meta_mx <- meta_mx[,2:ncol(meta_mx)] %>% as.matrix()
  rownames(meta_mx) <- row_n
  object <- SPATA2::addExpressionMatrix(object, expr_mtr = meta_mx, mtr_name = "Metabolites")
  
  sample_use_i <- which(samples_all==sample_use)
  path_set <- paste0("UKF", SPTCR[[sample_use_i]][1], "_", SPTCR[[sample_use_i]][2])
  Cell2loc <- read.csv(paste0("~/Desktop/ImmunoSpatial/New_data_04_05/Tutorial_Output/",path_set,"/",path_set,"_labels.csv"))
  Cell2loc$barcodes <- Cell2loc$spot_id
  cell_types <- names(Cell2loc)[11:64]
  Cell2loc <- Cell2loc %>% select(barcodes, {{cell_types}})
  names(Cell2loc)[2:ncol(Cell2loc)] <- paste0("C2L_", names(Cell2loc)[2:ncol(Cell2loc)])
  object <- addFeatures(object, Cell2loc)
  names(Cell2loc)
  ## Add T cell prams
  SPTCR_Index_all <- readRDS("~/Desktop/ImmunoSpatial/New_data_04_05/Tutorial_Output/SPTCR_Index_all.RDS")
  SPTCR_Index <- SPTCR_Index_all %>% filter(sample==sample_use) %>% dplyr::select(barcodes, TCR_UMI, number_clones,Diversity_index,Histology_IVY)
  SPTCR_Index <- left_join(data.frame(barcodes=SPATA2::getBarcodes(object)[[1]]), SPTCR_Index)
  
  object <- addFeatures(object, SPTCR_Index)
  
  ## Get the cytotoxic and exhausted clusters
  celltype <- c("C2L_Cd8.Cytotoxic",  "C2L_Cd8.Nk.Sig")
  
  object@information$active_mtr="scaled"
  object@information$active_expr_mtr="scaled"
  
  object <- SPATA2::runAutoencoderDenoising(object, bottleneck = 32, activation = "relu")
  
  object <- getCelltypeSpecificGeneExpression(object, 
                                              celltype = celltype,
                                              mtr_name="T_cell",
                                              active_mtr="denoised", 
                                              enhance_factor=10)
  
  cytotoxic <- c("GZMA", "GZMB", "GZMK", "PRF1", "IL2", "IFNG", "FASLG", "CCL5", "NKG7", "CD28")
  exhausted <- c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "EOMES", "BATF")
  
  object@information$active_mtr="T_cell"
  a <- joinWith(object, genes=cytotoxic, average_genes = T)
  b <- joinWith(object, genes=exhausted, average_genes = T)
  c <- joinWith(object, features = "TCR_UMI")
  a <- mutate(a, exhausted=b$mean_genes, TCR_UMI=c$TCR_UMI) %>% rename("cytotoxic":=mean_genes)
  
  
  
  a <- 
    a %>% 
    mutate(type_x=cytotoxic-exhausted) %>%
    mutate(type_y=c(cytotoxic+exhausted)/2) %>% 
    mutate(Tcell_type=case_when(
      TCR_UMI<5~"No",
      TCR_UMI>4 & type_x>0 & abs(type_x)>0.01 ~ "cytotoxic",
      TCR_UMI>5 & type_x<0 & abs(type_x)>0.06 ~ "exhausted",
      TRUE~"No"  ))
  
  
  a$Tcell_type <- factor(a$Tcell_type)
  
  
  ggplot(a, aes(x,y,color=Tcell_type))+
    scattermore::geom_scattermore(aes(x = x, y = y),pointsize = 7, color="black")+
    scattermore::geom_scattermore(aes(x = x, y = y),pointsize = 6, color="white")+
    geom_point()+theme_void()+
    scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set1")[1:2], "white"))+
    coord_fixed()+SPATA2::ggpLayerAxesSI(object)
  
  
  object <- addFeatures(object, a[,c("barcodes", "Tcell_type")], overwrite = T)
  
  saveRDS(object, paste0("~/Desktop/ImmunoSpatial/New_data_04_05/Tutorial_Output/Integrated_Metabolomics_",sample_use,".RDS"))
  
  
}


## Get T cell specific Metabolism
all_clones_metabolism <- map(.x=intersect(samples_all, list_matabolomic$sample), function(s){
  
  sample_use <- s
  sample_use_i <- which(samples_all==sample_use)
  
  object <- readRDS(paste0("~/Desktop/ImmunoSpatial/New_data_04_05/Tutorial_Output/Integrated_Metabolomics_",sample_use,".RDS"))
  
  object <- getCelltypeSpecificMetabolism(object, 
                                          celltype = celltype,
                                          mtr_name="T_cell_Metabolism",
                                          active_mtr="Metabolites", 
                                          enhance_factor=50)
  
  
  ## Clones
  # Clone level Data
  allSample_TCR_list <- map(1:9, ~readRDS(paste0("~/Desktop/ImmunoSpatial/New_data_04_05/Tutorial_Output/", samples_all[.x], "_TCR_output.RDS")))
  clones <- map_dfr(1:9, ~allSample_TCR_list[[.x]][[2]] %>% mutate(sample=samples_all[.x]))
  clones_select <- clones %>% filter(sample == sample_use)
  sample_use_i
  
  
  cytotoxic <- c("GZMA", "GZMB", "GZMK", "PRF1", "IL2", "IFNG", "FASLG", "CCL5", "NKG7", "CD28")
  exhausted <- c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "EOMES", "BATF")
  
  object@information$active_mtr="T_cell"
  a <- joinWith(object, genes=cytotoxic, average_genes = T)
  b <- joinWith(object, genes=exhausted, average_genes = T)
  c <- joinWith(object, features = "TCR_UMI")
  a <- mutate(a, exhausted=b$mean_genes, TCR_UMI=c$TCR_UMI) %>% rename("cytotoxic":=mean_genes)
  
  
  
  a <- 
    a %>% 
    mutate(type_x=cytotoxic-exhausted) %>%
    mutate(type_y=c(cytotoxic+exhausted)/2) %>% 
    mutate(Tcell_type=case_when(
      TCR_UMI<5~"No",
      TCR_UMI>4 & type_x>0 & abs(type_x)>0.01 ~ "cytotoxic",
      TCR_UMI>5 & type_x<0 & abs(type_x)>0.06 ~ "exhausted",
      TRUE~"No"  ))
  
  
  
  metabolic_pathways <- c("Glutathione_Metabolism", "Gluconeogenesis", 
                          "Pentose_Phosphate_Pathway", "Glycine_and_Serine_Metabolism","Amino_Sugar_Metabolism",
                          "Glycerolipid_Metabolism", "Glycolysis", "Tyrosine_Metabolism", "Fructose_and_Mannose_Degradation")
  
  metabolic_pw$ont[metabolic_pw$ont %>% str_detect(., pattern = "Amino")]
  
  SPATA2::getExpressionMatrixNames(object)
  object@information$active_mtr="T_cell_Metabolism"
  object@information$active_expr_mtr= "T_cell_Metabolism"
  object@used_genesets <- metabolic_pw
  
  
  df <- 
    dplyr::left_join(object %>%
                       SPATA2::joinWith(gene_sets = metabolic_pathways) %>% 
                       dplyr::select(-x,-y,-row,-col,-sample),
                     a %>%  
                       dplyr::select(barcodes, cytotoxic, exhausted),
                     by="barcodes") %>% 
    dplyr::select(-barcodes)
  
  
  
  target_clones <- clones_select %>% arrange(desc(nr_spots)) %>% head(50) %>% pull(CDR3_aa)
  
  clones_plot <- 
    allSample_TCR_list[[sample_use_i]]$spot_TCR %>% 
    filter(CDR3_aa %in% target_clones) %>% 
    group_by(barcodes, CDR3_aa) %>% 
    summarise(expression=mean(expression))
  
  
  clones_plot <- left_join(clones_plot, SPATA2::getCoordsDf(object), by="barcodes")
  
  anno_func <- dplyr::left_join(object %>%
                                  SPATA2::joinWith(gene_sets = metabolic_pathways) %>% 
                                  dplyr::select(-x,-y,-row,-col,-sample),
                                a %>%  
                                  dplyr::select(barcodes, cytotoxic, exhausted),
                                by="barcodes")
  
  clones_plot <- left_join(clones_plot, anno_func[,c("barcodes", "Glycolysis", "exhausted", "cytotoxic")])
  
  
  clones_plot <- clones_plot %>% group_by(CDR3_aa) %>% summarise(Glycolysis=mean(Glycolysis, na.rm=T), 
                                                                 exhausted=mean(exhausted,na.rm=T),
                                                                 cytotoxic=mean(cytotoxic,na.rm=T),
                                                                 expression=mean(expression,na.rm=T))
  
  
  return(list(df,clones_plot))
  
}, .progress = T)

clones_cor <- map_dfr(1:length(all_clones_metabolism), .f= 
                        ~all_clones_metabolism[[.x]][[2]] %>% 
                        mutate(sample=intersect(samples_all, list_matabolomic$sample)[.x]) %>% 
                        mutate(Glycolysis=scales::rescale(Glycolysis, c(0,1)),
                               exhausted=scales::rescale(exhausted, c(0,1)))
)

clones_cor$sample %>% unique()
ggplot(clones_cor)+# %>% filter(sample %in% c("275_T", "260_T", "269_T")))+
  geom_point(mapping = aes(x=exhausted, y=Glycolysis, size=expression, color=sample))+
  theme_classic()



genes_cor <- map(1:length(all_clones_metabolism), ~ cor(all_clones_metabolism[[.x]][[1]]))
cor.mat <- Reduce(`+`, genes_cor) / 5
corrplot::corrplot(cor.mat[colnames(cor.mat)[1:9], c("cytotoxic", "exhausted")], 
                   is.corr = F,
                   tl.col = "black",
                   tl.srt = 45,
                   pch.cex = 2,
                   sig.level=0.05,
                   insig = "label_sig",
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Greys")))(50))


clones_plot <- clones_cor
clones_plot[is.na(clones_plot)] <- 0
clones_plot$Glycolysis <- scales::rescale(clones_plot$Glycolysis, c(min(clones_plot$exhausted), max(clones_plot$exhausted)))


clones_plot[clones_plot$exhausted>0.5, ]$Glycolysis <- clones_plot[clones_plot$exhausted>0.5, ]$Glycolysis+0.3
clones_plot[clones_plot$exhausted<0.5, ]$Glycolysis <- clones_plot[clones_plot$exhausted<0.5, ]$Glycolysis-0.2
inter=0.5
ggplot(data=clones_plot %>% arrange(exhausted, desc(cytotoxic)))+
  geom_point(mapping=aes(x=1:nrow(clones_plot),y=exhausted, size=expression, color=sample))+
  geom_hline(yintercept = inter)+
  geom_segment(mapping=aes(x=1:nrow(clones_plot), xend=1:nrow(clones_plot), y=exhausted, yend=inter,color=sample, size=expression))+
  #geom_point(mapping=aes(x=1:nrow(clones_plot),y=Glycolysis),size=3, color="black")+
  #geom_segment(mapping=aes(x=1:nrow(clones_plot), xend=1:nrow(clones_plot), y=Glycolysis, yend=inter), color="black")+
  theme_classic()+
  scale_color_manual(values=c(RColorBrewer::brewer.pal(5, "Set1")))


ggplot(clones_plot)+
  geom_point(mapping = aes(x=exhausted, y=Glycolysis, size=expression, color=sample))+
  theme_classic()

sample_use <- "275_T"
sample_use_i <- which(samples_all==sample_use)

object <- readRDS(paste0("~/Desktop/ImmunoSpatial/New_data_04_05/Tutorial_Output/Integrated_Metabolomics_",sample_use,".RDS"))
object@used_genesets <- metabolic_pw

getCelltypeSpecificMetabolism <- function(object, 
                                          celltype, 
                                          enhance_factor=10, 
                                          mtr_name="T_cell_Metabolism", 
                                          active_mtr="Metabolites"){
  
  if(length(celltype)==1){
    message(paste0("Create expr Matrix from a single cell type"))
    
    
    
    object@information$active_mtr=active_mtr
    object@information$active_expr_mtr= active_mtr
    
    mat <- SPATA2::getExpressionMatrix(object)
    
    message(paste0("Run spot supressor"))
    sup <- joinWith(object, features = celltype) %>% mutate(f=exp(!!sym(celltype)))
    
    f <- sup %>% pull(f)
    f <- f+c(c(f-1)*enhance_factor)
    
    factor <- t(t(abs(mat)) * f)
    mat <- mat+factor
    
    message(paste0("Return matrix under the name: ", mtr_name ))
    object <- addExpressionMatrix(object, expr_mtr = mat,  mtr_name=mtr_name)
    
  }else{
    message(paste0("Create expr Matrix from a multiple cell types"))
    
    
    object@information$active_mtr=active_mtr
    object@information$active_expr_mtr= active_mtr
    mat <- SPATA2::getExpressionMatrix(object)
    
    
    message(paste0("Run spot supressor"))
    sup <- joinWith(object, features = celltype) %>% as.data.frame()
    sup_x <- sup[, celltype] %>% rowSums()/length(celltype)
    sup <- mutate(sup, f=exp(sup_x))
    
    f <- sup %>% pull(f)
    f <- f+c(c(f-1)*enhance_factor)
    
    factor <- t(t(abs(mat)) * f)
    mat <- mat+factor
    
    message(paste0("Return matrix under the name: ", mtr_name ))
    object <- addExpressionMatrix(object, expr_mtr = mat,  mtr_name=mtr_name)
    
    
  }
  
  return(object)
}

object <- getCelltypeSpecificMetabolism(object, 
                                        celltype = celltype,
                                        mtr_name="T_cell_Metabolism",
                                        active_mtr="Metabolites", 
                                        enhance_factor=50)



## Clones
# Clone level Data
allSample_TCR_list <- map(1:9, ~readRDS(paste0("~/Desktop/ImmunoSpatial/New_data_04_05/Tutorial_Output/", samples_all[.x], "_TCR_output.RDS")))
clones <- map_dfr(1:9, ~allSample_TCR_list[[.x]][[2]] %>% mutate(sample=samples_all[.x]))
clones_select <- clones %>% filter(sample == sample_use)
sample_use_i

target_clones <- clones_select %>% arrange(desc(nr_spots)) %>% head(10) %>% pull(CDR3_aa) %>% sample(10)

clones_plot <- 
  allSample_TCR_list[[sample_use_i]]$spot_TCR %>% 
  filter(CDR3_aa %in% target_clones) %>% 
  group_by(barcodes, CDR3_aa) %>% 
  summarise(expression=mean(expression))


clones_plot <- left_join(clones_plot, SPATA2::getCoordsDf(object), by="barcodes")
clones_plot$x=jitter(clones_plot$x, factor = 100);clones_plot$y=jitter(clones_plot$y, factor = 100)
clones_plot <- clones_plot %>% filter(expression>5)

cytotoxic <- c("GZMA", "GZMB", "GZMK", "PRF1", "IL2", "IFNG", "FASLG", "CCL5", "NKG7", "CD28")
exhausted <- c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "EOMES", "BATF")

object@information$active_mtr="T_cell"
a <- joinWith(object, genes=cytotoxic, average_genes = T)
b <- joinWith(object, genes=exhausted, average_genes = T)
c <- joinWith(object, features = "TCR_UMI")
a <- mutate(a, exhausted=b$mean_genes, TCR_UMI=c$TCR_UMI) %>% rename("cytotoxic":=mean_genes)



a <- 
  a %>% 
  mutate(type_x=cytotoxic-exhausted) %>%
  mutate(type_y=c(cytotoxic+exhausted)/2) %>% 
  mutate(Tcell_type=case_when(
    TCR_UMI<5~"No",
    TCR_UMI>4 & type_x>0 & abs(type_x)>0.01 ~ "cytotoxic",
    TCR_UMI>5 & type_x<0 & abs(type_x)>0.06 ~ "exhausted",
    TRUE~"No"  ))



getExpressionMatrixNames(object)
object@information$active_mtr="T_cell_Metabolism"
anno_func <- dplyr::left_join(object %>%
                                SPATA2::joinWith(gene_sets = metabolic_pathways) %>% 
                                dplyr::select(-x,-y,-row,-col,-sample),
                              a %>%  
                                dplyr::select(barcodes, cytotoxic, exhausted),
                              by="barcodes")

clones_plot <- left_join(clones_plot, anno_func[,c("barcodes", "Glycolysis", "exhausted", "cytotoxic")])
clones_plot <- na.omit(clones_plot)



background <- SPATA2::getCoordsDf(object)
col_clones <- RColorBrewer::brewer.pal(10, "Set3")


names(clones_plot)
col=colorRampPalette(c("#FFFFFF",RColorBrewer::brewer.pal(9, "Blues") ))
ggplot()+
  scattermore::geom_scattermore(data=background, mapping = aes(x,y),pointsize = 6, color="black")+
  scattermore::geom_scattermore(data=background, mapping = aes(x,y),pointsize = 5, color="white")+
  geom_point(data=clones_plot, mapping = aes(x,y,color=Glycolysis,size=expression))+
  #scale_color_manual(values=col_clones)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0.2,0.6),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()+
  SPATA2::ggpLayerAxesSI(object)

col=colorRampPalette(c("#FFFFFF",RColorBrewer::brewer.pal(9, "Reds") ))
col=colorRampPalette(c(RColorBrewer::brewer.pal(9, "RdBu") ))
ggplot()+
  scattermore::geom_scattermore(data=background, mapping = aes(x,y),pointsize = 6, color="black")+
  scattermore::geom_scattermore(data=background, mapping = aes(x,y),pointsize = 5, color="white")+
  geom_point(data=clones_plot, mapping = aes(x,y,color=exhausted,size=expression))+
  #scale_color_manual(values=col_clones)+
  scale_colour_gradientn(colours = col(50), 
                         limits= c(0,0.1),
                         oob = scales::squish, na.value = "white")+
  coord_fixed()+
  SPATA2::ggpLayerAxesSI(object)

## Clonal enrichment of metabolism and exhaustion


target_clones <- clones_select %>% arrange(desc(nr_spots)) %>% head(50) %>% pull(CDR3_aa) %>% sample(50)

clones_plot <- 
  allSample_TCR_list[[sample_use_i]]$spot_TCR %>% 
  filter(CDR3_aa %in% target_clones) %>% 
  group_by(barcodes, CDR3_aa) %>% 
  summarise(expression=mean(expression))


clones_plot <- left_join(clones_plot, SPATA2::getCoordsDf(object), by="barcodes")

anno_func <- dplyr::left_join(object %>%
                                SPATA2::joinWith(gene_sets = metabolic_pathways) %>% 
                                dplyr::select(-x,-y,-row,-col,-sample),
                              a %>%  
                                dplyr::select(barcodes, cytotoxic, exhausted),
                              by="barcodes")

clones_plot <- left_join(clones_plot, anno_func[,c("barcodes", "Glycolysis", "exhausted", "cytotoxic")])


clones_plot <- clones_plot %>% group_by(CDR3_aa) %>% summarise(Glycolysis=mean(Glycolysis, na.rm=T), 
                                                               exhausted=mean(exhausted,na.rm=T),
                                                               cytotoxic=mean(cytotoxic,na.rm=T),
                                                               expression=mean(expression,na.rm=T))

clones_plot$Glycolysis <- scales::rescale(clones_plot$Glycolysis, c(min(clones_plot$exhausted), max(clones_plot$exhausted)))


ggplot(data=clones_plot %>% arrange(exhausted, desc(cytotoxic)))+
  geom_point(mapping=aes(x=1:nrow(clones_plot),y=exhausted, size=expression))+
  geom_hline(yintercept = 0.057)+
  geom_segment(mapping=aes(x=1:nrow(clones_plot), xend=1:nrow(clones_plot), y=exhausted, yend=0.057, size=expression))+
  geom_point(mapping=aes(x=1:nrow(clones_plot),y=Glycolysis),size=3, color="red")+
  geom_segment(mapping=aes(x=1:nrow(clones_plot), xend=1:nrow(clones_plot), y=Glycolysis, yend=0.057), color="red")+
  theme_classic()