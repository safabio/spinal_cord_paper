### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 4 ####
## Fabio Sacher, 03.01.2024
### ### ### ### ### ### ### ### ### ### ###

library(Seurat)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(modplots)
library(miloR)
library(patchwork)

### ### ### ### ### ### ### ### ###
#### Milo Neighborhood DA plot ####
### ### ### ### ### ### ### ### ###

my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_lumb_int_seurat_250723.rds")
my.se$cond <- substr(my.se$orig.ident, 1, 4)

my.milo <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_lumb_int_milo_050225.rds")
da.results <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_lumb_int_milo_da_results050225.rds")

cust_pal <- c(
  "#923b37",
  "#b36d66",
  "#d0a099",
  "#edd4d1",
  "#F7F7F7",
  "#F7F7F7",
  "#F7F7F7",
  "#c8c2e0",
  "#9b92c7",
  "#6f64ad",
  "#3d3e99"
)

nhg <- plotNhoodGraphDA(my.milo, da.results, alpha=1, layout = "TSNE") +
  scale_fill_gradientn(colours = cust_pal, limits = c(-7.5,7.5))

pdf("~/spinal_cord_paper/figures/Fig_4_milo_network.pdf", width = 4, height = 4)
nhg +
  theme_void() +
  NoLegend()
nhg
dev.off()

pdf("~/spinal_cord_paper/figures/Fig_4_milo_volplot.pdf", height = 10, width = 5)
plotDAbeeswarm(da.results, group.by = "seurat_clusters", alpha=1.0) +
  ylim(c(-7.5,7.5)) +
  scale_color_gradientn(colours = cust_pal, limits = c(-7.5,7.5))
dev.off()

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
#### Glycogen Body dotplot #####
### ### ### ### ## ### ### ## ###

# Gg_brach_lumb_int
my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_lumb_int_seurat_250723.rds")
my.se@active.assay <- "RNA"
# indiidual labels
ctrl_lumb_int_combined_labels <- readRDS("~/spinal_cord_paper/annotations/ctrl_lumb_int_combined_labels.rds")

identical(colnames(my.se), rownames(ctrl_lumb_int_combined_labels))
my.se$annot_sample  <- ctrl_lumb_int_combined_labels$annot_sample

# split clusters by sample
my.se$sample <- str_extract(my.se$orig.ident, "ctrl|lumb")
my.se$split_clusters <- paste(my.se$seurat_clusters, my.se$sample, sep = "_")

Idents(my.se) <- "split_clusters" 
my.se@active.assay <- "RNA"

# Get the top 50 markers (cluster vs rest)
bl10int_markers <- read.table("~/spinal_cord_paper/data/Gg_ctrl_lumb_int_fullDE_250723.txt") %>% 
  filter(cluster %in% c(21, 22)) %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 50) 

# shared top 50 genes
dup_genes <- bl10int_markers$Gene.name[duplicated(bl10int_markers$Gene.name)]

markers_avihu <- c("GBE1","GYG2","MSX1","MSX2","PAX3","LMX1A","ENSGALG00000028041","ROR2","WNT7B","GDF7")

intersect(markers_avihu, dup_genes)
# [1] "GBE1" "MSX2"

# remove duplicates from shared markers and cand. list
GOI <- c(markers_avihu, dup_genes)
GOI <- GOI[!duplicated(GOI)]

gnames <- modplots::gnames

IOI <- gnames[gnames$Gene.name %in% GOI,] %>% 
  mutate(Gene.name = factor(Gene.name, levels = GOI))

IOI <- IOI[match(GOI, IOI$Gene.name), ]

my.sub <- subset(my.se, subset = sample == "lumb")

mdot_clust <- modplots::mDotPlot2(my.sub,
                        group.by = "split_clusters",
                        features = rev(IOI$Gene.stable.ID), 
                        gnames = modplots::gnames,
                        cluster.idents = TRUE,
                        cols = c("lightgrey", "black"))

# order factors on ordered dotplot
my.sub$split_clusters <- factor(my.sub$split_clusters,
                                levels = mdot_clust[[2]]$labels[mdot_clust[[2]]$order])

p1 <- modplots::mDotPlot2(my.sub,
                          group.by = "split_clusters",
                          features = rev(IOI$Gene.stable.ID), 
                          gnames = modplots::gnames,
                          cols = c("lightgrey", "black")) +
  coord_flip() +
  ggtitle(paste0(my.se@project.name, " marker Dotplot")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5))

# hclust as dendrogram
dend.idents.cand <- as.dendrogram(mdot_clust[[2]])

pdf("~/spinal_cord_paper/figures/Fig_4_GB_dotplot.pdf", height = 6.5, width = 10.5)
plot(dend.idents.cand)
p1
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

cl_ord <- c(2,5,6,10,11,13,19,23,25,15,12,21,22,3,7,17,20,27,28,30,1,4,8,9,24,26,29,33,16,18,34,14,31,32)

my.se$sample <- str_extract(my.se$orig.ident, "ctrl|lumb")
my.se$split_clusters <- paste(my.se$seurat_clusters, my.se$sample, sep = "_")
my.se$split_clusters <- factor(my.se$split_clusters, 
                               levels = c(paste0(c(cl_ord), "_ctrl"), paste0(c(cl_ord), "_lumb")))

Idents(my.se) <- "split_clusters" 
my.se@active.assay <- "RNA"

gnames = modplots::gnames

## select HOX genes
hox_select <- rev(c("HOXA11","HOXD3","HOXD4","HOXD8","HOXD9","HOXD10","HOXD11"))

cand <- modplots::gnames %>% 
  filter(Gene.name %in% hox_select)

cand <- cand[match(hox_select, cand$Gene.name),]

dpl_hox_select <- modplots::mDotPlot2(
  my.se, 
  features = cand$Gene.stable.ID,
  cols = c("lightgrey", "black"),
  gnames = gnames
)

lab_cols <- ifelse(grepl("lumb", levels(Idents((my.se)))), "#419c73", "black")


ggsave(
  filename = "~/spinal_cord_paper/figures/Fig_4_ctrl_lumb_hox_selected_dotplot.pdf",
  width = 15, height = 3,
  plot = dpl_hox_select +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = lab_cols))
)

## all HOX genes
hox_alphanum <- modplots::gnames %>% 
  filter(grepl("^HOX", Gene.name)) %>% 
  mutate(hoxnum = as.numeric(str_remove(Gene.name, "HOX[ABCD]"))) %>% 
  mutate(hoxgroup = substr(Gene.name, 4, 4)) %>% 
  arrange(hoxgroup) %>% 
  arrange(hoxnum)

dpl_hox_alphanum <- modplots::mDotPlot2(
  my.se, 
  features = rev(hox_alphanum$Gene.stable.ID),
  cols = c("lightgrey", "black"),
  gnames = gnames) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = lab_cols))

hox_groups <- modplots::gnames %>% 
  filter(grepl("^HOX", Gene.name)) %>% 
  mutate(hoxnum = as.numeric(str_remove(Gene.name, "HOX[ABCD]"))) %>% 
  mutate(hoxgroup = substr(Gene.name, 4, 4)) %>% 
  arrange(hoxnum) %>% 
  arrange(hoxgroup)

dpl_hox_groups <- modplots::mDotPlot2(
  my.se, 
  features = rev(hox_groups$Gene.stable.ID),
  cols = c("lightgrey", "black"),
  gnames = gnames)+
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = lab_cols))

pdf("~/spinal_cord_paper/figures/Fig_4_ctrl_lumb_hox_all_dotplots.pdf", width = 16, height = 9)
dpl_hox_alphanum
dpl_hox_groups
dev.off()

### ### ### ### ### ### ### ### 
#### GB module network plot ####
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

### ### ### ### ### ### ### ### ### ### 
#### GB module networks on B/L10int ####
### ### ### ### ### ### ### ### ### ### 

# devel_int WGCNA data
wgcna <- readRDS("~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds")
# MN modules magenta and brown4
modules <- c(26, 5)

dev_mods <- wgcna$module.genes[modules]
names(dev_mods) <- names(wgcna$modules)[modules]

# Lumb int seurat 
my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_lumb_int_seurat_250723.rds")
# calculate scores
my.se <- AddModuleScore(
  my.se,
  dev_mods,
  name = "module",
  assay = "RNA"
)

colnames(my.se[[]])[(ncol(my.se[[]])-length(dev_mods)+1):ncol(my.se[[]])] <- names(dev_mods)

head(my.se[[]])

mod_plots <- list()
void_plots <- list()

for (i in seq(modules)) {
  # get the module color
  curr_mod <- str_remove(names(dev_mods)[i], "^\\d{1,2}_")
  # plot
  mod_plots[[i]] <- FeaturePlot(
    my.se,
    names(dev_mods)[i],
    cols = c("gray90", curr_mod),
    reduction = "tsne",
    order = TRUE
  ) 
  void_plots[[i]] <- FeaturePlot(
    my.se,
    names(dev_mods)[i],
    cols = c("gray90", curr_mod),
    reduction = "tsne", 
    order = TRUE
  )  + theme_void()+ NoLegend()
  
}

pdf("~/spinal_cord_paper/figures/Fig_4_L10int_GB_modules_on_BL10int_addmodulescore.pdf")
gridExtra::grid.arrange(grobs = mod_plots[1])
gridExtra::grid.arrange(grobs = mod_plots[2])
gridExtra::grid.arrange(grobs = void_plots[1])
gridExtra::grid.arrange(grobs = void_plots[2])
dev.off()

# avg expression of module of interest
# the lumbar cells have the _3 and _4 added to cellIDs
MOI <- wgcna[["sc.MEList"]][["averageExpr"]] %>% 
  tibble::rownames_to_column("CellID") %>% 
  mutate(CellID = str_replace(CellID, "_1$", "_3")) %>% 
  mutate(CellID = str_replace(CellID, "_2$", "_4")) %>% 
  tibble::column_to_rownames("CellID") %>% 
  select(AEmagenta, AEbrown4)

table(rownames(MOI) %in% rownames(my.se[[]]))

my.lumb <- readRDS("~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds")

MOI_lumb <- wgcna[["sc.MEList"]][["averageExpr"]] %>%
  select(AEmagenta, AEbrown4)

table(rownames(MOI_lumb) %in% rownames(my.lumb[[]]))

my.se <- AddMetaData(my.se, MOI)
my.lumb <- AddMetaData(my.lumb, MOI_lumb)

my.se <- AddMetaData(my.se, Embeddings(my.se, reduction = "tsne"))
my.lumb <- AddMetaData(my.lumb, Embeddings(my.lumb, reduction = "tsne"))

pt_size <- 2

# brachial lumb integrated
BL10int_M <- my.se[[]] %>% 
  arrange(AEmagenta) %>%  
  drop_na(AEmagenta) %>%
  ggplot(aes(
    x = tsne_1,
    y = tsne_2)
  ) +
  geom_point(aes(color = AEmagenta), size = pt_size) +
  scale_colour_gradientn(colours = c("gray90","gray80","magenta")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"))

BL10int_B <- my.se[[]] %>% 
  arrange(AEbrown4) %>%  
  drop_na(AEbrown4) %>%
  ggplot(aes(
    x = tsne_1,
    y = tsne_2)
  ) +
  geom_point(aes(color = AEbrown4), size = pt_size) +
  scale_colour_gradientn(colours = c("gray90","gray80","brown4")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"))

# lumbar only
L10int_M <- my.lumb[[]] %>% 
  arrange(AEmagenta) %>%  
  drop_na(AEmagenta) %>%
  ggplot(aes(
    x = tsne_1,
    y = tsne_2)
  ) +
  geom_point(aes(color = AEmagenta), size = pt_size) +
  scale_colour_gradientn(colours = c("gray90","gray80","magenta")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"))

L10int_B <- my.lumb[[]] %>% 
  arrange(AEbrown4) %>%  
  drop_na(AEbrown4) %>%
  ggplot(aes(
    x = tsne_1,
    y = tsne_2)
  ) +
  geom_point(aes(color = AEbrown4), size = pt_size) +
  scale_colour_gradientn(colours = c("gray90","gray80","brown4")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"))


pdf("~/spinal_cord_paper/figures/Fig_4_L10int_GB_modules_on_B10int_and_BL10int_labeltransfer.pdf")

BL10int_M + ggtitle("Magenta BL10int")
BL10int_M + theme_void() + NoLegend()

BL10int_B + ggtitle("Brown4 BL10int")
BL10int_B + theme_void() + NoLegend()

L10int_M + ggtitle("Magenta L10int")
L10int_M + theme_void() + NoLegend()

L10int_B + ggtitle("Brown4 L10int")
L10int_B + theme_void() + NoLegend()

dev.off()

### ### ### ### ### ### ### ### 
#### B/L10int cl-3 vs cl-7 DE ####
### ### ### ### ### ### ### ### 

# seurat object
my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_lumb_int_seurat_250723.rds")
my.se@active.assay <- "RNA"

DimPlot(my.se, reduction = "tsne", label = TRUE)

gnames <- modplots::gnames
# markers cl3 vs cl7
tmp <- FindMarkers(my.se,
                   ident.1 = 3,
                   ident.2 = 7,
                   only.pos = FALSE,
                   min.pct = 0.25,
                   logfc.threshold =  0.2,
                   latent.vars = c("CC.Difference.seurat"),
                   test.use = "MAST",
                   assay = "RNA")
# filter for avg_log2FC > 0.5
cl3v7DE <- tmp %>%
  rownames_to_column("Gene.stable.ID") %>% 
  left_join(gnames) %>% 
  filter(abs(avg_log2FC) > 0.5)

# dorsal and ventral markers for cl3 and 7
GOI <- c("PBX3","MEIS1","SNCG","PROX1","PRDM8","BHLHE22","PAX2")

gnames <- modplots::gnames

IOI <- gnames[gnames$Gene.name %in% GOI,] %>% 
  mutate(Gene.name = factor(Gene.name, levels = GOI))

IOI <- IOI[match(GOI, IOI$Gene.name), ]
# subset for cl3 and cl7 only
my.sub <- subset(x = my.se, idents = c(3, 7))
# dotplot of selected markers
dot_sel <- mDotPlot2(my.sub,
                     group.by = "seurat_clusters",
                     features = rev(IOI$Gene.stable.ID), 
                     gnames = modplots::gnames,
                     cols = c("lightgrey", "black")) +
  coord_flip() +
  ggtitle(paste0(my.se@project.name, " dorsal ventral markers Dotplot")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5))
# top 25 markers per cluster by log2FC
top_markers <- cl3v7DE %>% 
  mutate(cluster = case_when(
    avg_log2FC > 0 ~ "cl3",
    avg_log2FC < 0 ~ "cl7"
  )) %>% 
  group_by(cluster) %>% 
  slice_max(abs(avg_log2FC), n = 25) 
# dotplot of top 25 markers
dot_top <- mDotPlot2(my.sub,
                     group.by = "seurat_clusters",
                     features = rev(top_markers$Gene.stable.ID), 
                     gnames = modplots::gnames,
                     cols = c("lightgrey", "black")) +
  coord_flip() +
  ggtitle(paste0(my.se@project.name, " cl3 v cl7 top DE marker dotplot")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5))

pdf("~/spinal_cord_paper/figures/Supp_fig_4_BL10_int_cl3_v_cl7_dotplots.pdf", height = 8, width = 4)
dot_top
dot_sel
dev.off()
