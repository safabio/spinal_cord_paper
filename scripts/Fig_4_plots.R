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

### ### ### ### ### ### ### ### ### ### ###
#### detailed venn diagram of RP and FP module overlap ####
### ### ### ### ### ### ### ### ### ### ###

wgcna_dat <- list(
  ctrl = readRDS("~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds"),
  lumb = readRDS("~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds"),
  poly = readRDS("~/spinal_cord_paper/output/Gg_poly_int_scWGCNA_250723.rds")
)

# get the module genes
ctrl_greenyellow  <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[1]]$module.genes[[7]]) 
ctrl_tan          <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[1]]$module.genes[[19]]) 
ctrl_purple       <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[1]]$module.genes[[15]]) 

lumb_darkred      <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[2]]$module.genes[[13]])
lumb_lightgreen   <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[2]]$module.genes[[23]])
lumb_yellow       <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[2]]$module.genes[[49]])

poly_purple       <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[3]]$module.genes[[14]])
poly_ligthgreen   <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[3]]$module.genes[[9]])
poly_black        <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% wgcna_dat[[3]]$module.genes[[1]])

venn_rp = list(ctrl_greenyellow = ctrl_greenyellow$Gene.name,
               lumb_darkred = lumb_darkred$Gene.name, 
               poly_purple = poly_purple$Gene.name)

venn_fp = list(ctrl_tan = ctrl_tan$Gene.name,
               lumb_lightgreen = lumb_lightgreen$Gene.name, 
               poly_ligthgreen = poly_ligthgreen$Gene.name)


venn_cilia = list(ctrl_purple = ctrl_purple$Gene.name,
                  lumb_yellow = lumb_yellow$Gene.name, 
                  poly_black = poly_black$Gene.name)

p_lab_rp <- ggvenn(venn_rp, fill_color = c("greenyellow", "darkred", "purple"), show_elements = TRUE, text_size = 2) +
  ggtitle("Roof plate modules D10 int")
p_rp <- ggvenn(venn_rp, fill_color = c("greenyellow", "darkred", "purple")) +
  ggtitle("Roof plate modules D10 int")

p_lab_fp <- ggvenn(venn_fp, fill_color = c("tan", "lightgreen", "lightgreen"), show_elements = TRUE, text_size = 2) +
  ggtitle("Floor plate modules D10 int")
p_fp <- ggvenn(venn_fp, fill_color = c("tan", "lightgreen", "lightgreen")) +
  ggtitle("Floor plate modules D10 int")

p_lab_cilia <- ggvenn(venn_cilia, fill_color = c("purple", "yellow", "black"), show_elements = TRUE, text_size = 2) +
  ggtitle("Cilia/Flagella modules D10 int")
p_cilia <- ggvenn(venn_cilia, fill_color = c("purple", "yellow", "black")) +
  ggtitle("Cilia/Flagella modules D10 int")

pdf("~/spinal_cord_paper/figures/RP_FP_cilia_modules_venn_D10.pdf", width = 7, height = 7)
p_lab_rp
p_rp
p_lab_fp
p_fp
p_lab_cilia
p_cilia
dev.off()

Reduce(intersect, venn_rp)
Reduce(intersect, venn_fp)
Reduce(intersect, venn_cilia)


rm(wgcna_dat)

### ### ### ### ### ### ### ### ### ### ### ### ###
#### ctrl_int m-greenyellow and m-tan overlay ####
### ### ### ### ### ### ### ### ### ### ### ### ###

my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds")

my.wgcna <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds")

avgExp <- my.wgcna[["sc.MEList"]]$averageExpr

if (identical(rownames(my.se[[]]), rownames(avgExp))) {
  my.se <- AddMetaData(my.se, avgExp)
}

cl17 <- subset(x = my.se, idents = 17)

p <-
  FeaturePlot(
    my.se,
    features = c("AEgreenyellow", "AEtan"),
    cols = c("grey90", "greenyellow", "tan"),
    blend = TRUE,
    reduction = "tsne",
    order = TRUE,
    combine = FALSE
  )

p_cl17 <-
  FeaturePlot(
    cl17,
    features = c("AEgreenyellow", "AEtan"),
    cols = c("grey90", "greenyellow", "tan"),
    blend = TRUE,
    reduction = "tsne",
    order = TRUE,
    combine = FALSE,
    pt.size = 3
  )


pdf("~/spinal_cord_paper/figures/Ctrl_int_RP_FP_module_overlap.pdf", width = 8, height = 7)
(p_cl17[[1]] + p_cl17[[2]])/
  (p_cl17[[4]] + p_cl17[[3]])

(p[[1]] + p[[2]])/
  (p[[4]] + p[[3]])
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ###
#### tSNE plots of RP/FP/cilia modules overlap ####
### ### ### ### ### ### ### ### ### ### ### ### ###

data_sets <- c("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

# RP overlap
RP_overlap <- c("GYG2","COL1A2","GDF10","ENSGALG00000052026","LMX1A","OLFML3",
             "WNT3A","ENSGALG00000033591","KLHDC8B","UTRN","TSPAN18")

gnames <- modplots::gnames

plots <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = RP_overlap, my.slot = "data", size = 0.5, return = TRUE)
  # set empty legend title
  for (j in seq(RP_overlap)) {
    plots[[i]][[j]][["labels"]][["colour"]] <- ""
  }
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/RP_module_overlap_tsne.pdf", width = 18, height = 14)
grid.arrange(grobs = plots[[1]])
grid.arrange(grobs = plots[[2]])
grid.arrange(grobs = plots[[3]])
dev.off()

# FP overlap
FP_overlap <- c("ASB9","SLC6A4","ECM2","CNP1","TLL1","NTN1",
                "SHH","HPGDS","PTGDS","EXFABP","KRT24","CMTM8","ADIPOQ")

plots <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = FP_overlap, my.slot = "data", size = 0.5, return = TRUE)
  # set empty legend title
  for (j in seq(FP_overlap)) {
    plots[[i]][[j]][["labels"]][["colour"]] <- ""
  }
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/FP_module_overlap_tsne.pdf", width = 18, height = 14)
grid.arrange(grobs = plots[[1]])
grid.arrange(grobs = plots[[2]])
grid.arrange(grobs = plots[[3]])
dev.off()

rm(data_sets, RP_overlap, FP_overlap, gnames, plots)

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
