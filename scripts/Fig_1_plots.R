### ### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 1 ###
## Fabio Sacher, 04.06.2024 ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ###

### ### ### ### ### ### ### ### ###
#### Gg_all_int pseudobulk pca ####
### ### ### ### ### ### ### ### ###

library(stringr)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)

se_path <- c(
    "Gg_ctrl_1_seurat_070323",
    "Gg_ctrl_2_seurat_070323",
    "Gg_lumb_1_seurat_070323",
    "Gg_lumb_2_seurat_070323",
    "Gg_poly_1_seurat_070323",
    "Gg_poly_2_seurat_070323",
    "Gg_D05_ctrl_seurat_070323",
    "Gg_D07_ctrl_seurat_070323"
  )

my.samples <- c(
    "Gg_ctrl_1",
    "Gg_ctrl_1",
    "Gg_lumb_1",
    "Gg_lumb_2",
    "Gg_poly_1",
    "Gg_poly_2",
    "Gg_ctrl_D05",
    "Gg_ctrl_D07"
  )

clust_col <- read.csv("~/spinal_cord_paper/annotations/broad_cluster_marker_colors.csv") %>% 
  dplyr::rename(broad = broad_cluster) %>% 
  select(-marker) %>% 
  filter(!grepl("FP/RP", broad))


broad_order <- c(
    "progenitors",
    "FP",
    "RP",
    "neurons",
    "OPC",
    "MFOL",
    "pericytes",
    "microglia",
    "blood",
    "vasculature"
  )

my.meta <- list()
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
    dplyr::rename(seurat_clusters = number) 
  
  # cluster order for vln plots
  ord_levels[[i]] <- annot$fine[order(match(annot$broad, broad_order))]
  
  # create index for color coding
  col_table[[i]] <- annot %>%
    left_join(clust_col, by = "broad") %>% 
    select(c("fine", "color"))
  
  # add cluster annotation to meta data
  my.meta[[i]] <- my.se@meta.data %>% 
    rownames_to_column("rowname") %>% 
    left_join(annot, by = "seurat_clusters") %>% 
    mutate(fine = factor(fine, levels = annot$fine)) %>% 
    column_to_rownames("rowname") %>% 
    mutate(label = paste(old.ident, fine, sep = "_")) %>%
    select(label, broad) %>% 
    rownames_to_column("cell_id") %>% 
    mutate(cell_id = paste(cell_id, i, sep = "_"))
  
  rm(my.se)
  
}

names(my.meta) <- my.samples
names(ord_levels) <- my.samples
names(col_table) <- my.samples

# select broad and sample level fine clustering
meta <- do.call(rbind, my.meta) %>% 
  remove_rownames() %>% 
  column_to_rownames("cell_id") %>% 
  mutate(broad = factor(broad, levels = clust_col$broad))

write.table(meta, sep = ",", file = "~/spinal_cord_paper/annotations/Gg_all_int_transferred_labels.csv")

my.int <- readRDS("~/spinal_cord_paper/data/Gg_all_int_seurat_270524.rds")

# add sample level fine clustering
my.int <- AddMetaData(my.int, meta)

Idents(my.int) <- "broad" 

clust_col <- clust_col %>% 
  filter(broad %in% broad_order)

DimPlot(
  object = my.int,
  reduction = "tsne",
  cols = clust_col$color,
  label = TRUE,
  label.size = 5
) +
  ggtitle(NULL) +
  theme_void() +
  NoLegend()

Idents(my.int) <- "label"

# Pseudobulk
pb <- AggregateExpression(my.int, assays = "RNA", slot = "count")
names(colnames(pb[["RNA"]])) <- str_remove(colnames(pb[["RNA"]]), "data[, 1]")


meta_idx <- meta %>% filter(!duplicated(label))

coldata <- data.frame(colnames(pb[[1]])) %>% 
  dplyr::rename(label = colnames.pb..1...) %>% 
  mutate(orig.ident = str_extract(label, "Gg_\\D{4}_[12]|Gg_D?\\d{1,2}_ctrl")) %>% 
  rownames_to_column("rownames") %>% 
  left_join(meta_idx, by = "label") %>% 
  column_to_rownames("rownames")

library(DESeq2)
# Create DDS object, vstransform and calculate pca
dds <- DESeqDataSetFromMatrix(countData = pb[[1]],
                              colData = coldata,
                              design = ~ orig.ident)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

ngenes <- 1000
point.size <- 2
pch.order <- c(0,5,2,6,3,4,1,8)
# PC1 and 2
my.pca12 <-
  modplots::mPCA(
    vsd,
    PCs = c("PC1", "PC2"),
    ntop = ngenes,
    group = "broad",
    colors = clust_col$color,
    shape = "orig.ident",
    pch = pch.order,
    pt.size = point.size
  )
# PC1 and 3
my.pca13 <-
  modplots::mPCA(
    vsd,
    PCs = c("PC1", "PC3"),
    ntop = ngenes,
    group = "broad",
    colors = clust_col$color,
    shape = "orig.ident",
    pch = pch.order,
    pt.size = point.size
  )
# PC4 and 2 (reversed so plot aligns)
my.pca42 <-
  modplots::mPCA(
    vsd,
    PCs = c("PC4", "PC2"),
    ntop = ngenes,
    group = "broad",
    colors = clust_col$color,
    shape = "orig.ident",
    pch = pch.order,
    pt.size = point.size
  )

# combine PC plots into one
pca_ptchw <- my.pca12[[2]] +
  my.pca42[[2]] +
  my.pca13[[2]] +
  guide_area() +
  plot_layout(guides = "collect") & theme(legend.direction = "vertical", legend.box = "horizontal")

ggsave(
  filename = "~/spinal_cord_paper/figures/Supp_Fig_1_Gg_all_int_pca.pdf",
  width = 8, height = 7,
  plot = pca_ptchw
)

### ### ### ### ### ### ###
#### HINTW violin plot ####
### ### ### ### ### ### ###

# combined samples seurat object
my.int <- readRDS("~/spinal_cord_paper/data/Gg_all_int_seurat_270524.rds")
DefaultAssay(my.int) <- "RNA"

my.int[[]]$orig.ident <- factor(my.int[[]]$orig.ident,
                                levels = c("Gg_D05_ctrl",
                                           "Gg_D07_ctrl",
                                           "Gg_ctrl_1",
                                           "Gg_ctrl_2",
                                           "Gg_lumb_1",
                                           "Gg_lumb_2",
                                           "Gg_poly_1",
                                           "Gg_poly_2"))

colnames(my.int[[]])

gnames <- modplots::gnames
# HINTW ID (marker for female cells/samples)
hintw <- gnames[grepl("HINTW", gnames$Gene.name),]

cols <- c(
  "#A4A4A4",
  "#515151",
  "#000000",
  "#000000",
  "#419c73",
  "#419c73",
  "goldenrod3",
  "goldenrod3"
)

svol <- VlnPlot(
  my.int,
  features = hintw$Gene.stable.ID,
  group.by = "orig.ident",
  cols = cols
) +
  ggtitle("HINTW expression by sample")
# remove beeswarm
svol_empty <- VlnPlot(
  my.int,
  features = hintw$Gene.stable.ID,
  group.by = "orig.ident",
  pt.size = 0,
  cols = cols
) +
  ggtitle("HINTW expression by sample")

# export plot
ggsave("~/spinal_cord_paper/figures/Violinplot_HINTW_by_sample.pdf",
       width = 7,
       height = 7,
       (svol / svol_empty) + plot_layout(guides = "collect"))


