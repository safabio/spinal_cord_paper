### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 2 ####
## Fabio Sacher, 21.05.2024
### ### ### ### ### ### ### ### ### ### ###

### ### ### ### ### ### ### ### ### ### ###
#### scWGCNA module network plots ####
### ### ### ### ### ### ### ### ### ### ###

# code is taken and adapted from scWGCNA::scW.p.network

library(ggplot2)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(patchwork)

# scWGNA.data
scWGCNA.data = readRDS("~/spinal_cord_paper/output/Gg_devel_int_scWGCNA_250723.rds")
gnames = modplots::gnames
rownames(gnames) <- gnames[,1]
modules = c("darkred", "darkgreen", "lightgreen", "tan")


# get the module colors
my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))

nwrk <- list()

text_size <- list(5,5,5,2)
maxN_size <- list(23,23,23,15)

for (i in seq(modules)) {
  col <- modules[i]
  # get index if color is provided
  module = which(my.cols == col)
  # extract the module
  mynet = scWGCNA.data[["networks"]][[module]]
  
  gene.labs = gnames[network::network.vertex.names(mynet),2]
  set.seed(42)
  
  GGally::ggnet2(mynet,
                  mode = "fruchtermanreingold",
                  layout.par = list(repulse.rad = network::network.size(mynet)^1.1,
                                    area = network::network.size(mynet)^2.3),
                  node.size = network::get.vertex.attribute(mynet, "membership01"),
                  max_size = maxN_size[[i]],
                  node.color = col, 
                  edge.size = "weight02", 
                  edge.color = "black",
                  edge.alpha = network::get.edge.attribute(mynet, "weight01")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_text(ggplot2::aes(label = gene.labs),
                       alpha = 1,
                       color = "black",
                       fontface = "italic",
                       size = text_size[[i]]
  )
  
  ggsave(paste0("~/spinal_cord_paper/figures/Fig_2_", col,"_network.pdf"), width = 5, height = 5)
}

### ### ### ### ### ### ### ### ### ### ###
#### tan and darkgreen avg exp tsne ####
### ### ### ### ### ### ### ### ### ### ###

###
# 10.01.2025
# Ubuntu R-4.4.1
# code taken from spinal_cord_paper/markdown/Gg_devel_scWGCNA_module_analysis.html
# plotting tan tsne and violin plot
###

for (i in seq(my.ses)) {
  # prepare average expression table
  tmp <- avg.mod.eigengenes[,1:22] %>%
    tibble::rownames_to_column("cell_ID") %>%
    dplyr::filter(grepl(paste0("_", i, "$"), cell_ID)) %>%
    dplyr::mutate(cell_ID = stringr::str_remove_all(cell_ID, paste0("_", i))) %>%
    tibble::column_to_rownames("cell_ID")
  
  identical(rownames(tmp), colnames(my.ses[[i]]))
  # add meta data to the seurat objects
  my.ses[[i]] <- AddMetaData(my.ses[[i]], tmp)
}

#max and min expression per module (column max)
mod_max <- apply(avg.mod.eigengenes[,1:22], MARGIN = 2, FUN = max)[module_order]
mod_min <- apply(avg.mod.eigengenes[,1:22], MARGIN = 2, FUN = min)[module_order]

modplots <- list()
modplots[[1]] <- list()
modplots[[2]] <- list()
modplots[[3]] <- list()

modules_in_order <- colnames(tmp)[module_order]


# plot the modules split to the stages
for (i in seq(my.ses)) {
  for (j in seq(ncol(tmp))) {
    
    modplots[[i]][[j]]  <- FeaturePlot(
      my.ses[[i]], order = TRUE,
      features = modules_in_order[j],
      reduction = "tsne"
    ) +
      ggtitle(stringr::str_remove(modules_in_order[j],"^AE")) +
      scale_color_gradient(low="gray90", high=substring(modules_in_order[j], 3), #colors in the scale
                           limits=c(mod_min[j], mod_max[j])) #same limits for plots
    
    
  }
}

full_plot <- c(modplots[[1]], modplots[[2]], modplots[[3]])

# tSNE plot
pdf("~/spinal_cord_paper/figures/Supp_Fig_2_modules_darkgreen_tan_AE_plots.pdf", width = 13, height = 8)
(full_plot[[8]] + full_plot[[30]] + full_plot[[52]] +
    full_plot[[9]] + full_plot[[31]] + full_plot[[53]]) +
  plot_layout(ncol = 3, guides = "collect")
dev.off()

### ### ### ### ### ### ### ### ### ### ###
#### Fig 2 violin plots scWGCNA mod exp ####
### ### ### ### ### ### ### ### ### ### ###

# load seurat objects and add annotations
se_path <- c("Gg_D05_ctrl_seurat_070323",
             "Gg_D07_ctrl_seurat_070323",
             "Gg_ctrl_1_seurat_070323")

# broad cluster order
broad_order <- c("progenitors",
                 "FP",
                 "RP",
                 "FP/RP",
                 "neurons",
                 "OPC",
                 "MFOL",
                 "pericytes",
                 "microglia",
                 "blood",
                 "vasculature"
)

# custom colors
clust_col <- read.csv("~/spinal_cord_paper/annotations/broad_cluster_marker_colors.csv") %>% 
  rename(broad = broad_cluster) %>% 
  select(-marker)

my.ses <- list()
col_table <- list()
ord_levels <- list()

for (i in seq(se_path)) {
  # load the data sets
  my.se <- readRDS(paste0("~/spinal_cord_paper/data/", se_path[i], ".rds"))
  annot <- read.csv(list.files("~/spinal_cord_paper/annotations",
                               pattern = str_remove(se_path[i], "_seurat_\\d{6}"),
                               full.names = TRUE))

  if(length(table(annot$number)) != length(table(my.se$seurat_clusters))) {
    stop("Number of clusters must be identical!")
    
  }
  # rename for left join
  annot <- annot %>% 
    mutate(fine = paste(fine, number, sep = "_")) %>% 
    mutate(number = factor(number, levels = 1:nrow(annot))) %>% 
    rename(seurat_clusters = number) 
  # create index for color coding
  col_table[[i]] <- annot %>%
    left_join(clust_col, by = "broad") %>% 
    select(c("fine", "color"))
  # cluster order for vln plots
  ord_levels[[i]] <- annot$fine[order(match(annot$broad, broad_order))]
  
  my.se@meta.data <- my.se@meta.data %>% 
    rownames_to_column("rowname") %>% 
    left_join(annot, by = "seurat_clusters") %>% 
    mutate(fine = factor(fine, levels = annot$fine)) %>% 
    column_to_rownames("rowname")

  my.ses[[i]] <- my.se
}


names(my.ses) <- c("D05", "D07", "D10")
names(col_table) <- c("D05", "D07", "D10")
names(ord_levels) <- c("D05", "D07", "D10")

# rm(my.se, annot)

# scWGCNA data
WGCNA_data = list()
WGCNA_data[[1]] = readRDS("~/spinal_cord_paper/output/Gg_devel_int_scWGCNA_250723.rds")
my.wsub =list()
my.wsub[[1]]= c(1:22)

# This is just to add a little bit more sense to the modules, so that we don't get just a color. Corresponds to WGCNA_data
my.modulenames = list()
my.modulenames[[1]] = c(1:22)

# broad cluster color table
all_col <- do.call(rbind, col_table) %>%
  rownames_to_column("sample") %>%
  mutate(sample = substr(sample, 1, 3)) %>%
  mutate(sample_celltype = paste(sample, fine, sep = "_")) %>%
  select(c("color", "sample_celltype", "sample"))

#Get a dataframe with annotations for all the samples and colors we need.
my.metam <- list()

for (i in seq(my.ses)) {
  my.metam[[i]] <- my.ses[[i]][[]]
  rownames(my.metam[[i]]) <- paste0(rownames(my.metam[[i]]), "_", i)
}

my.metam <- do.call(rbind, my.metam)

my.metam$orig.ident <- str_replace_all(my.metam$orig.ident, pattern =  "Gg_D05_ctrl", "D05")
my.metam$orig.ident <- str_replace_all(my.metam$orig.ident, pattern =  "Gg_D07_ctrl", "D07")
my.metam$orig.ident <- str_replace_all(my.metam$orig.ident, pattern =  "Gg_ctrl_1", "D10")
# my.metam$sample_celltype = paste0(substr(my.metam$orig.ident,7,9),"_",my.metam$seurat_clusters)
my.metam$sample_celltype = paste0(my.metam$orig.ident, "_", my.metam$fine)

my.metam <- my.metam %>% 
  tibble::rownames_to_column(var = "cell_ID") %>%
  dplyr::left_join(all_col, by = "sample_celltype") %>%
  tibble::column_to_rownames(var = "cell_ID")

cell_table <- my.metam %>% 
  rownames_to_column("cell_ID") %>% 
  select("sample_celltype", "cell_ID")

avg.mod.eigengenes <- WGCNA_data[[1]]$sc.MEList$averageExpr %>% 
  rownames_to_column("cell_ID") %>% 
  left_join(cell_table, by = "cell_ID") %>% 
  column_to_rownames("cell_ID")

# mean eigengenes
avg.mod.eigengenes.mean <- avg.mod.eigengenes %>%
  group_by(sample_celltype) %>%
  summarise_all("mean") %>%
  column_to_rownames("sample_celltype")

# integrated data set of all devel samples.
my.sec = readRDS("~/spinal_cord_paper/data/Gg_devel_int_seurat_250723.rds")

identical(rownames(my.metam), colnames(my.sec))

my.metam$sample_celltype <- factor(my.metam$sample_celltype, levels = all_col$sample_celltype)
#Set the identities of the integrated data, to the annotated clusters
my.sec = SetIdent(my.sec, value = my.metam$sample_celltype)

identical(rownames(avg.mod.eigengenes), colnames(my.sec))

# add meta data to the int seurat object
my.sec <- AddMetaData(my.sec, avg.mod.eigengenes)
my.sec <- AddMetaData(my.sec, my.metam[c("fine", "sample_celltype")])

custom_order <- c(paste(names(ord_levels)[1], ord_levels[[1]], sep = '_'),
                  paste(names(ord_levels)[2], ord_levels[[2]], sep = '_'),
                  paste(names(ord_levels)[3], ord_levels[[3]], sep = '_'))

my.sec$sample_celltype <- factor(
  my.sec$sample_celltype,
  levels = custom_order
)

vln_ind <- VlnPlot(my.sec,
                   features = c("AEdarkred", "AElightgreen", "AEdarkgreen", "AEtan"),
                   group.by = "sample_celltype",
                   stack = TRUE,
                   flip = TRUE,
                   cols = c("darkred", "lightgreen", "darkgreen", "tan"),pt.size = 1
) +
  geom_boxplot(outlier.shape = NA) +
  NoLegend()

pdf("~/spinal_cord_paper/figures/Fig_2_AE_selected_mod.pdf", height = 25, width = 20)
vln_ind
dev.off()