### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 4 ####
## Fabio Sacher, 03.01.2024
### ### ### ### ### ### ### ### ### ### ###

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(magrittr)
library(tibble)
library(cowplot)
library(gridExtra)


### ### ### ### ### ### ### ### ###
#### cluster sizes bar plots ####
### ### ### ### ### ### ### ### ###


my.se <- list()
# load seurat objects
my.se[[1]] <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds")
my.se[[2]] <- readRDS("~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds")
my.se[[3]] <- readRDS("~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

names(my.se) <- c("ctrl_int", "lumb_int", "poly_int")

annot <- list()
annot[[1]] <- read.csv("~/spinal_cord_paper/annotations/Gg_ctrl_int_cluster_annotation.csv") %>% 
  mutate(fine = paste(number, fine, sep = "_"))
annot[[2]] <- read.csv("~/spinal_cord_paper/annotations/Gg_lumb_int_cluster_annotation.csv") %>% 
  mutate(fine = paste(number, fine, sep = "_"))
annot[[3]] <- read.csv("~/spinal_cord_paper/annotations/Gg_poly_int_cluster_annotation.csv") %>% 
  mutate(fine = paste(number, fine, sep = "_"))

plots <- list()

for (i in seq_along(my.se)) {
  # cluster ids and id summary (frequency table)
  my.se[[i]]$seurat_clusters <- factor(my.se[[i]]$seurat_clusters, labels = annot[[i]]$fine)
  cl_sizes <- data.frame(cluster = my.se[[i]]$seurat_clusters)
  # plot cluster size 
  plots[[i]] <- ggplot(cl_sizes, aes(x = cluster)) +
    geom_bar(fill = "lightgrey", color = "black",stat = "count") + 
    stat_count(geom = "text", colour = "black", size = 2.5,
               aes(label = after_stat(count)),position=position_stack(vjust=0.5)) +
    ggtitle(names(my.se)[i]) +
    theme_cowplot() +
    scale_x_discrete(labels = levels(cl_sizes$cluster)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # plot total number of cells
    annotate("text", x = 20, y = 500, label = paste("N cells =", dim(my.se[[i]])[2]))
  
}

grid.arrange(grobs = plots, ncol = 1)

pdf(file = "~/spinal_cord_paper/figures/Int_data_cluster_sizes.pdf", height = 13)
grid.arrange(grobs = plots, ncol = 1)
dev.off()

### ### ### ### ### ### ### ### 
#### Ctrl_lumb int dotplot ####
### ### ### ### ### ### ### ###

# seurat object
my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_lumb_int_seurat_250723.rds")
# cluster labels from B10int and L10int
ctrl_lumb_int_combined_labels <- readRDS("~/spinal_cord_paper/annotations/ctrl_lumb_int_combined_labels.rds")

identical(colnames(my.se), rownames(ctrl_lumb_int_combined_labels))
my.se$annot_sample  <- ctrl_lumb_int_combined_labels$annot_sample

my.se@active.assay <- "RNA"

my.se[[]]

my.se$sample <- str_extract(my.se$orig.ident, "ctrl|lumb")
my.se$split_clusters <- paste(my.se$seurat_clusters, my.se$sample, sep = "_")

Idents(my.se) <- "split_clusters" 
my.se@active.assay <- "RNA"

gnames = modplots::gnames

## select HOX genes
hox_select <- rev(c("HOXB5","HOXC6","HOXC9","HOXA10","HOXC10","HOXD10","HOXA11","HOXD11"))

cand <- modplots::gnames %>% 
  filter(Gene.name %in% hox_select)

cand <- cand[match(hox_select, cand$Gene.name),]

dpl_hox_select <- modplots::mDotPlot2(
  my.se, 
  features = cand$Gene.stable.ID,
  cols = c("lightgrey", "black"),
  cluster.idents = TRUE,
  gnames = gnames
)

dpl_hox_select[[1]] +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# hclust as dendrogram
dend.idents.hox <- as.dendrogram(dpl_hox_select[[2]])

pdf("~/spinal_cord_paper/figures/Fig_4_ctrl_lumb_hox_selected_dotplot_dendro.pdf")
plot(dend.idents.hox)
dev.off()  

ggsave(
  filename = "~/spinal_cord_paper/figures/Fig_4_ctrl_lumb_hox_selected_dotplot.pdf",
  width = 13, height = 3,
  plot = dpl_hox_select[[1]] +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
)

# Dotplot and clustering on the selected HOX genes from above 
# order the factors 
my.se@meta.data$split_clusters <- factor(
  my.se@meta.data$split_clusters,
  levels =  dpl_hox_select[[2]][["labels"]][dpl_hox_select[[2]][["order"]]])

# my.se@meta.data$split_clusters <- forcats::fct_rev(my.se@meta.data$split_clusters)
Idents(my.se) <- "split_clusters" 

hox <- modplots::gnames %>% 
  filter(grepl("^HOX", Gene.name)) %>% 
  filter(Gene.stable.ID %in% rownames(my.se)) %>% 
  mutate(group = str_extract(Gene.name, "HOX.")) %>% 
  mutate(number = as.numeric(str_extract(Gene.name, "\\d{1,2}$"))) %>% 
  arrange(group) %>% 
  arrange(number)

dpl_hox <- modplots::mDotPlot2(my.se,
                               features = rev(hox$Gene.stable.ID),
                               cols = c("lightgrey", "black"),
                               gnames = gnames)

dpl_hox +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(
  filename = "~/spinal_cord_paper/figures/Supp_fig_4_ctrl_lumb_hox_dotplot.pdf",
  width = 13, height = 8,
  plot = dpl_hox +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
)


## Marker gene dotplot
GOI <- rev(c("SOX9","SHH","RSPO1","TUBB3","NRXN3","TLX3","OLIG2","PLP1","IGFBP7","IFI30","HBBA","CDH5"))

cand <- modplots::gnames %>% 
  filter(Gene.name %in% GOI)

cand <- cand[match(GOI, cand$Gene.name),]


# Dotplot of candidate genes (GSI clust. from my.htmp)
dpl_cand <- modplots::mDotPlot2(
  my.se,
  features = cand$Gene.stable.ID,
  cols = c("lightgrey", "black"), 
  cluster.idents = TRUE,
  gnames = gnames
)

dpl_cand[[1]] +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# hclust as dendrogram
dend.idents.cand <- as.dendrogram(dpl_cand[[2]])

pdf("~/spinal_cord_paper/figures/Fig_4_ctrl_lumb_cand_dotplot_dendro.pdf")
plot(dend.idents.cand)
dev.off()  


ggsave(
  filename = "~/spinal_cord_paper/figures/Fig_4_ctrl_lumb_cand_dotplot.pdf",
  width = 13, height = 4,
  plot = dpl_cand[[1]] +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
)

### ### ### ### ### ### ### ### 
#### GB module network plot ###
### ### ### ### ### ### ### ###

# scWGNA.data
scWGCNA.data = readRDS("~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds")
gnames = modplots::gnames
rownames(gnames) <- gnames[,1]
modules = c("brown4","magenta")

# get the module colors
my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))

nwrk <- list()

text_size <- list(5,2)
maxN_size <- list(20,15)

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
  
  ggsave(paste0("~/spinal_cord_paper/figures/Fig_4_", col,"_network.pdf"), width = 5, height = 5)
}

