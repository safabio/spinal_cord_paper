### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 5 ####
## Fabio Sacher, 21.12.23
### ### ### ### ### ### ### ### ### ### ###

library(Seurat)
library(tidyverse)
library(gridExtra)
library(magrittr)
library(ggbeeswarm)
library(miloR)
library(patchwork)

### ### ### ### ### ### ### ### ###
#### Milo Neighborhood DA plot ####
### ### ### ### ### ### ### ### ###

my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_seurat_250723.rds")
my.se$cond <- substr(my.se$orig.ident, 1, 4)

my.milo <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_poly_int_milo_050225.rds")
da.results <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_poly_int_milo_da_results050225.rds")

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

pdf("~/spinal_cord_paper/figures/Fig_5_milo_network.pdf", width = 4, height = 4)
nhg +
  theme_void() +
  NoLegend()
nhg
dev.off()

pdf("~/spinal_cord_paper/figures/Fig_5_milo_volplot.pdf", height = 10, width = 5)
plotDAbeeswarm(da.results, group.by = "seurat_clusters", alpha=1.0) +
  ylim(c(-7.5,7.5)) +
  scale_color_gradientn(colours = cust_pal, limits = c(-7.5,7.5))
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Compare ctrl and poly 1 & 2 FP/RP cluster size ####
### ### ### ### ### ### ### ### ### ### ### ### ### ###

data_sets <- c("~/spinal_cord_paper/data/Gg_ctrl_1_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_ctrl_2_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_poly_1_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_poly_2_seurat_070323.rds")
# cluster ids (ctlr 2 has only one FP/RP cluster)
RP_clust <- c(22,17,19,10)
FP_clust <- c(18,17,13,12)

gnames <- modplots::gnames

plots <- list()
# cluster size data frame
RP_size <-
  data.frame(
    "size" = c(0, 0, 0, 0),
    "sample" = c("ctrl_1", "ctrl_2", "poly_1", "poly_2"),
    "clust" = RP_clust,
    "type" = rep("RP", 4)
  )
FP_size <-
  data.frame(
    "size" = c(0, 0, 0, 0),
    "sample" = c("ctrl_1", "ctrl_2", "poly_1", "poly_2"),
    "clust" = FP_clust,
    "type" = rep("FP", 4)
  )

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  # plot RP, FP, and markers for both
  plots[[i]] <- mFeaturePlot(
    my.se, #          RP       Both     FP
    my.features = c("LMX1A", "SLIT2", "NTN1",
                    "RSPO1", "LECT1", "SHH"),
    my.slot = "scale.data",
    size = 1
  )
  # cluster size
  RP_size[i, 1] <- table(my.se@active.ident)[[RP_clust[i]]]
  FP_size[i, 1] <- table(my.se@active.ident)[[FP_clust[i]]]
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/RP_FP_ctrl_poly_1_2_tsne.pdf", width = 12, height = 7)
grid.arrange(plots[[1]])
grid.arrange(plots[[2]])
grid.arrange(plots[[3]])
dev.off()

toplot <- rbind(RP_size, FP_size) %>% 
  mutate(ID = str_extract(sample, "\\d{1}$")) %>% 
  mutate(sample = str_remove(sample, "_\\d{1}$")) %>% 
  mutate(type = factor(type)) %>% 
  mutate(sample = factor(sample))

p <- ggplot(toplot, aes(x = sample, y = type, color = size, size = size, label = size)) +
  geom_count() +
  scale_color_gradient(low = "grey", high = "forestgreen") +
  facet_wrap("ID") +
  geom_text(cex = 3, col = "black") +
  scale_size_area(max_size = 30) +
  geom_text(aes(label = paste0("cl-", clust)),
            col = "black", cex = 5, nudge_y = 0.2) +
  ggtitle("Cluster sizes of FP and RP clusters in ctrl and poly 1 & 2 samples")

pdf("~/spinal_cord_paper/figures/RP_FP_cluster_sizes.pdf", width = 7, height = 7)
p + theme_classic()
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ###
### BP10int DE genes flat volcanoplot (miloR insp) ####
### ### ### ### ### ### ### ### ### ### ### ### ### ###

# ### ### ### ### ### #
# Do not run again. DE calculation takes a long time!
# ### ### ### ### ### #

# my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_seurat_250723.rds")
# # cluster labels from B10int and L10int
# ctrl_poly_int_combined_labels <- readRDS("~/spinal_cord_paper/annotations/ctrl_poly_int_combined_labels.rds")
# 
# identical(colnames(my.se), rownames(ctrl_poly_int_combined_labels))
# 
# my.se$annot_sample  <- ctrl_poly_int_combined_labels$annot_sample
# 
# my.se@active.assay <- "RNA"
# 
# markers <- list()
# numbers <- list()
# composition <- list()
# 
# for (i in seq(levels(Idents(my.se)))) {
#   # subset for individual clusters
#   mn.se <- subset(x = my.se, idents = levels(Idents(my.se))[i])
#   mn.se$sample <- str_extract(mn.se$orig.ident, "ctrl|poly")
#   
#   composition[[i]] <- mn.se[[]] %>% 
#     select(sample, annot_sample)
#   
#   Idents(mn.se) <- "sample"
#   
#   tmp_markers <- FindMarkers(mn.se,
#                              ident.1 = "ctrl",
#                              only.pos = FALSE, 
#                              min.pct = 0.25, 
#                              logfc.threshold =  0.2,
#                              latent.vars = c("CC.Difference.seurat"), 
#                              test.use = "MAST", 
#                              assay = "RNA")
#   # cell numbers per sample
#   numbers[[i]] <- data.frame(table(mn.se$sample))
#   
#   tmp_markers <- tmp_markers %>%
#     rownames_to_column("Gene.stable.ID") %>% 
#     left_join(gnames)
#   
#   
#   markers[[i]] <- tmp_markers
# }
# 
# names(markers) <- paste0("cl-", levels(Idents(my.se)))
# names(numbers) <- paste0("cl-", levels(Idents(my.se)))
# names(composition) <- paste0("cl-", levels(Idents(my.se)))
# 
# # bind lists into data frames
# poly_markers <- bind_rows(markers, .id = "cluster") %>% 
#   mutate(cluster = factor(cluster, levels = paste0("cl-", levels(Idents(my.se)))))
# poly_numbers <- bind_rows(numbers, .id = "cluster") %>% 
#   mutate(cluster = factor(cluster, levels = paste0("cl-", levels(Idents(my.se)))))
# poly_composition <- bind_rows(composition, .id = "cluster") %>% 
#   mutate(cluster = factor(cluster, levels = paste0("cl-", levels(Idents(my.se)))))
# 
# saveRDS(poly_markers, "~/spinal_cord_paper/data/Gg_ctrl_poly_int_markers.rds")
# saveRDS(poly_numbers, "~/spinal_cord_paper/data/Gg_ctrl_poly_int_numbers.rds")
# saveRDS(poly_composition, "~/spinal_cord_paper/data/Gg_ctrl_poly_int_composition.rds")

poly_markers <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_markers.rds") %>% 
  mutate(clust_id = str_remove(cluster, "^cl-")) %>% 
  mutate(clust_id = factor(clust_id, levels = c(1:32))) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(sample = case_when(
    avg_log2FC > 0 ~ "B10int",
    avg_log2FC < 0 ~ "Poly10int"
  )) %>% 
  # transform p_val_adj for alpha plotting
  mutate(neg_log_10_p_val_adj = -log10(p_val_adj)) %>% 
  # for plotting we treshold the pvals it
  mutate(p_alpha = case_when(
    neg_log_10_p_val_adj > 20 ~ 20,
    TRUE ~ neg_log_10_p_val_adj
  )) 


# extract DA values from da.results
da_BL10int <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_poly_int_milo_da_results050225.rds")

miloR::plotDAbeeswarm(da_BL10int, group.by = "seurat_clusters", alpha=1.0)

# get the avg abs DA logFC per cluster
da_abs_logFC <- da_BL10int %>% 
  select(logFC, seurat_clusters) %>% 
  mutate(abs_logFC = abs(logFC)) %>% 
  group_by(seurat_clusters) %>% 
  summarise_all(mean) %>% 
  mutate(seurat_clusters_num = as.numeric(as.character(seurat_clusters))) %>% 
  mutate(cluster = factor(paste0("cl-", seurat_clusters)))

da_abs_logFC[da_abs_logFC$seurat_clusters=="Mixed", 
             "seurat_clusters_num"] <- 0

ggplot(poly_markers, aes(x = avg_log2FC, y = cluster)) +
  geom_point() +
  scale_y_discrete(drop = FALSE)

# testing wether the logFC calculations were correct
vol_logFC <- miloR::plotDAbeeswarm(da_BL10int, group.by = "seurat_clusters", alpha=1.0) + 
  geom_point(data = da_abs_logFC, 
             mapping = aes(x=rev(seurat_clusters_num+1),
                           y=logFC),
             pch = 3,
             inherit.aes = FALSE) +
  ggtitle("mean logFC change")

vol_abslogFC <- miloR::plotDAbeeswarm(da_BL10int, group.by = "seurat_clusters", alpha=1.0) + 
  geom_point(data = da_abs_logFC, 
             mapping = aes(x=rev(seurat_clusters_num+1),
                           y=abs_logFC),
             pch = 3,
             inherit.aes = FALSE) +
  ggtitle("mean abs(logFC) change")

#code taken and modified from miloR::plotDAbeeswarm()
beeswarm_pos <- ggplot_build(
  poly_markers %>% 
    arrange(clust_id) %>%
    ggplot(aes(clust_id, avg_log2FC)) +
    geom_quasirandom() +
    scale_x_discrete(drop = FALSE)
)

pos_x <- beeswarm_pos$data[[1]]$x
pos_y <- beeswarm_pos$data[[1]]$y
n_groups <- length(levels(poly_markers$clust_id))

abs_logFC_df <- da_abs_logFC %>% 
  select(cluster, abs_logFC) %>% 
  # rename to prevent confusion
  rename("abs_DA_logFC" = "abs_logFC")

alpha <- poly_markers %>% 
  arrange(clust_id) %>% 
  mutate(pos_x = pos_x, pos_y = pos_y) %>% 
  ggplot(aes(pos_x, pos_y, color = sample, alpha = p_alpha)) +
  scale_color_manual(values = c("black", "goldenrod3")) +
  xlab("clust_id") +
  ylab("Log Fold Change") + 
  scale_alpha_continuous(range = c(0,1)) +
  scale_x_continuous(breaks = seq(1, n_groups),
                     limits = c(0.5,n_groups+0.5),
                     labels = setNames(levels(poly_markers$clust_id), seq(1, n_groups))) +
  geom_point(pch = 16) + 
  ylim(c(-3,3)) +
  coord_flip() + 
  theme_bw(base_size = 22) +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(angle = 180))

alpha_size <- poly_markers %>% 
  arrange(clust_id) %>%
  left_join(abs_logFC_df, by = "cluster") %>%
  # min max (scale to 0 - 1) and invert (-1*x)+1
  mutate(scale_inv_abs_DA_logFC = (-1*((abs_DA_logFC-min(abs_DA_logFC))/(max(abs_DA_logFC)-min(abs_DA_logFC))))+1) %>%
  mutate(pos_x = pos_x, pos_y = pos_y) %>% 
  ggplot(aes(pos_x,
             pos_y, 
             color = sample, 
             alpha = p_alpha, 
             size = scale_inv_abs_DA_logFC)) +
  scale_size(range = c(0.1,2)) +
  scale_color_manual(values = c("black", "goldenrod3")) +
  xlab("clust_id") +
  ylab("Log Fold Change") + 
  scale_alpha_continuous(range = c(0.1,1)) +
  scale_x_continuous(breaks = seq(1, n_groups),
                     limits = c(0.5, n_groups+1), 
                     labels = setNames(levels(poly_markers$clust_id), seq(1, n_groups))) +
  geom_point(pch = 16) + 
  ylim(c(-3,3)) +
  coord_flip() + 
  theme_bw(base_size = 22) +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(angle = 180))

pdf("~/spinal_cord_paper/figures/Fig_5_BPoly10int_flat_DE_volplot.pdf", height = 10, width = 5)
vol_logFC
vol_abslogFC
alpha + 
  ggtitle("aes(alpha)")
alpha + 
  NoLegend()
alpha_size + 
  labs(alpha = "-log10(p_val_adj)\nthres = 20",
                  size = "scaled_inv_\nabs_DA_logFC") +
  ggtitle("aes(alpha, size)")
alpha_size + 
  NoLegend()
dev.off()

### ### ### ### ### ### ### 
###  GO terms of ctrl_poly DE ####
### ### ### ### ### ### ### 

library(tidyverse)
library(org.Gg.eg.db)
library(limma)

# load the marker genes
poly_markers <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_markers.rds") %>% 
  mutate(clust_id = str_remove(cluster, "^cl-")) %>% 
  mutate(clust_id = factor(clust_id, levels = c(1:32))) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(sample = case_when( # marker for B or poly?
    avg_log2FC > 0 ~ "B10int",
    avg_log2FC < 0 ~ "Poly10int"
  )) %>% 
  mutate(GO_groups = factor(paste0(sample, "_", clust_id))) %>% 
  droplevels()

# order factor groups
levs <- c(paste0(rep("B10int", 32), "_", 1:32), paste0(rep("Poly10int", 32), "_", 1:32))
levs <- levs[levs %in% levels(poly_markers$GO_groups)]

poly_markers <- poly_markers %>% 
  mutate(GO_groups = factor(GO_groups, levels = levs))

load("~/spinal_cord_paper/data/gga_KEGG_gene_pathway.Rda")
load("~/spinal_cord_paper/data/gga_KEGG_path_name.Rda")

my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_seurat_250723.rds")
# expressed universe
my.gouni = rownames(my.se[["RNA"]]@counts)[which(apply(my.se[["RNA"]]@counts, 1, function(x) length(which(x >0))) > 0)]

rm(my.se)

# archive export of entrez ids
entrez_id <- read.delim("~/crosstalk/data/mart_export_102_entrez_ids.tsv") 
colnames(entrez_id) <- c("id","entrez","name","entrez_name")

universe <- entrez_id %>% 
  na.omit() %>%
  filter(id %in% my.gouni) %>%
  pull(entrez) %>%
  unique() 

#Create lists
entrez_ids <- list()
gowg <- list()
kewg <- list()

my.sp <-  "Gg"

# A loop that goes through the DEG files we created earlier
for(i in seq(length(levs))){ #as many levels as we have
  gene_ids <- poly_markers %>%
    filter(GO_groups == levs[i]) %>%
    pull(Gene.stable.ID)
  
  entrez_ids[[i]] <- entrez_id %>%
    na.omit() %>%
    filter(id %in% gene_ids) %>%
    pull(entrez) %>%
    unique() 
  
  gowg[[i]] = goana(entrez_ids[[i]], species = my.sp, universe = universe) #The GO analysis
  kewg[[i]] = kegga(entrez_ids[[i]], species = my.sp, universe = universe, gene.pathway = GeneID.PathID, pathway.names = PathID.PathName)
}

names(gowg) <- levs
names(kewg) <- levs

keggpath=list()
goterms=list()

for (i in seq(length(levs))) {
  #### TopGO terms ####
  toplot=limma::topGO(gowg[[i]], n= 50, ontology = c("BP"))
  
  #Change to character and then back to factor, to keep the order from TopGO
  toplot$Term = as.character(toplot$Term)
  toplot$Term = factor(toplot$Term, levels = unique(toplot$Term))
  
  goterms[[i]]=toplot
  
  rm(toplot)
  
  #### TopKEGG pathways ####
  toplot=limma::topKEGG(kewg[[i]], n= 50)
  
  #Change to character and then back to factor, to keep the order from TopGO
  toplot$Pathway = as.character(toplot$Pathway)
  toplot$Pathway = factor(toplot$Pathway, levels = unique(toplot$Pathway))
  
  keggpath[[i]]=toplot

  rm(toplot)
}

names(goterms) <- levs
names(keggpath) <- levs

coltoplot <- data.frame(levs) %>% 
  mutate(color = case_when(
    grepl("B10", levs) ~ "black",
    TRUE ~ "goldenrod3"
  )) %>% 
  pull(color)

goplots <- list()

for (i in seq(goterms)) {
  goplots[[i]] <- ggplot(goterms[[i]], aes(x=Term, y=-log10(P.DE), shape=Ont)) +
    geom_point(size = 3, color = coltoplot[[i]]) +
    labs(title=paste0("Top GOTerm of ", names(goterms)[[i]], " DE genes")) +
    scale_x_discrete(labels=strtrim(goterms[[i]]$Term, 50)) +
    coord_flip()
}

names(goplots) <- levs

# export top 50 tables 
saveRDS(goterms, "~/spinal_cord_paper/data/Gg_ctrl_poly_int_markers_top50_GOterms.rds")
# export GO term enrich barplot
pdf("~/spinal_cord_paper/figures/Fig_5_Gg_ctrl_poly_int_markers_top50_GOterms.pdf")
goplots
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ###
# barplot of n DEG for B/poly10int ####
### ### ### ### ### ### ### ### ### ### ### ### ### ###

library(tidyverse)

poly_markers <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_markers.rds") %>% 
  mutate(clust_id = str_remove(cluster, "^cl-")) %>% 
  mutate(clust_id = factor(clust_id, levels = c(1:32))) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(sample = case_when( # marker for B or poly?
    avg_log2FC > 0 ~ "B10int",
    avg_log2FC < 0 ~ "Poly10int"
  )) %>% 
  mutate(DE_groups = factor(paste0(sample, "_", clust_id)))

DE_bar <- ggplot(poly_markers, aes(x = clust_id)) +
  geom_bar(position="dodge") +
  # scale_fill_manual(values = c("goldenrod3", "black")) +
  scale_x_discrete(drop = FALSE) +
  cowplot::theme_cowplot()

factor_order <- c(1,3,7,9,13,16,17,18,19,27,23,22,5,6,11,20,25,29,30,32,2,4,8,12,14,15,26,10,24,21,28,31)

DE_bar_ordered <- poly_markers %>% 
  mutate(clust_id = factor(clust_id, levels = factor_order)) %>% 
  ggplot(aes(x = clust_id)) +
  geom_bar(position="dodge") +
  # scale_fill_manual(values = c("goldenrod3", "black")) +
  scale_x_discrete(drop = FALSE) +
  cowplot::theme_cowplot()

pdf("~/spinal_cord_paper/figures/Supp_Fig_5_Gg_ctrl_poly_int_markers_barplot.pdf", width = 10, height = 7)
DE_bar/
DE_bar_ordered
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ###
# bulk heatmaps DEG B/Poly10int ####
### ### ### ### ### ### ### ### ### ### ### ### ### ###

library(tidyverse)
library(Seurat)
library(pheatmap)

gnames <- modplots::gnames

my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_seurat_250723.rds")

my.se[[]] <- my.se[[]] %>% 
  mutate(side = str_extract(orig.ident, "ctrl|poly")) %>% 
  mutate(clust_side = paste(seurat_clusters, side, sep = "-"))

NPC <- c(paste(c("g1","g3","g7","g9","g13","g16","g17","g18","g19","g22","g23","g27"), "ctrl", sep = "-"),
         paste(c("g1","g3","g7","g9","g13","g16","g17","g18","g19","g22","g23","g27"), "poly", sep = "-"))
OPC <- c(paste(c("g2","g4","g8","g10","g12","g14","g15","g24","g26"), "ctrl", sep = "-"),
         paste(c("g2","g4","g8","g10","g12","g14","g15","g24","g26"), "poly", sep = "-"))
neuron <-  c(paste(c("g5","g6","g11","g20","g25","g29","g30","g32"), "ctrl", sep = "-"),
             paste(c("g5","g6","g11","g20","g25","g29","g30","g32"), "poly", sep = "-"))

pb <- AggregateExpression(
  my.se,
  assays = "RNA",
  group.by = "clust_side")


top_markers_noMT <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_markers.rds") %>% 
  mutate(clust_id = str_remove(cluster, "^cl-")) %>% 
  mutate(clust_id = factor(clust_id, levels = c(1:32))) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(sample = case_when( # marker for B or poly?
    avg_log2FC > 0 ~ "ctrl",
    avg_log2FC < 0 ~ "poly"
  ))  %>% 
  mutate(DE_groups = factor(paste0("g", clust_id, "-", sample))) %>% 
  droplevels() %>% 
  mutate(abs_logFC = abs(avg_log2FC)) %>% 
  group_by(cluster) %>% 
  filter(!grepl("^ND\\d$|^MT-|^CYTB|^ATP|^COX\\d|^COII$", Gene.name)) %>% 
  slice_max(n = 10, order_by = abs_logFC) %>% 
  mutate(Gene.name = str_replace(Gene.name, "^ENSGALG0+", "EG"))

mark_NPC_noMT <- top_markers_noMT[top_markers_noMT$DE_groups %in% NPC,]
mark_OPC_noMT <- top_markers_noMT[top_markers_noMT$DE_groups %in% OPC,]
mark_neuron_noMT <- top_markers_noMT[top_markers_noMT$DE_groups %in% neuron,]

### ### ### ### ### ### ### ### ### ### ### ### ### ###
# DE genes row and col cluster noMT ####
### ### ### ### ### ### ### ### ### ### ### ### ### ###

pb_NPC <- pb[["RNA"]][unique(mark_NPC_noMT$Gene.stable.ID),NPC] %>% 
  data.frame() %>% 
  rownames_to_column("Gene.stable.ID") %>% 
  mutate(Gene.stable.ID = str_remove(Gene.stable.ID, "\\.\\d+$")) %>% 
  left_join(gnames) %>% 
  mutate(Gene.name = str_replace(Gene.name, "^ENSGALG0+", "EG")) %>% 
  column_to_rownames("Gene.name") %>% 
  select(-Gene.stable.ID)

pb_OPC <- pb[["RNA"]][unique(mark_OPC_noMT$Gene.stable.ID),OPC] %>% 
  data.frame() %>% 
  rownames_to_column("Gene.stable.ID") %>% 
  mutate(Gene.stable.ID = str_remove(Gene.stable.ID, "\\.\\d+$")) %>% 
  left_join(gnames) %>% 
  mutate(Gene.name = str_replace(Gene.name, "^ENSGALG0+", "EG")) %>% 
  column_to_rownames("Gene.name") %>% 
  select(-Gene.stable.ID)

pb_neuron <- pb[["RNA"]][unique(mark_neuron_noMT$Gene.stable.ID),neuron] %>% 
  data.frame() %>% 
  rownames_to_column("Gene.stable.ID") %>% 
  mutate(Gene.stable.ID = str_remove(Gene.stable.ID, "\\.\\d+$")) %>% 
  left_join(gnames) %>% 
  mutate(Gene.name = str_replace(Gene.name, "^ENSGALG0+", "EG")) %>% 
  column_to_rownames("Gene.name") %>% 
  select(-Gene.stable.ID)

ann_color <- list(
  `cl-1` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-2` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-3` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-4` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-5` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-6` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-7` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-8` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-9` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-10` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-11` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-12` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-13` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-14` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-15` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-16` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-17` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-18` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-19` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-20` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-21` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-22` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-23` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-24` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-25` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-26` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-27` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-28` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-29` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-30` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-31` = c(
    yes = "grey60",
    no = "white"
  ),
  `cl-32` = c(
    yes = "grey60",
    no = "white"
  ),
  origin = c(
    ctrl = "black",
    poly = "goldenrod3"
  )
) 

ann_col_NPC <- data.frame(
  cluster = c(NPC)
) %>% 
  mutate(cluster = str_replace(cluster, "-", "\\.")) %>% 
  mutate(origin = str_extract(cluster, "ctrl|poly")) %>% 
  column_to_rownames("cluster")

ann_col_OPC <- data.frame(
  cluster = c(OPC)
) %>% 
  mutate(cluster = str_replace(cluster, "-", "\\.")) %>% 
  mutate(origin = str_extract(cluster, "ctrl|poly")) %>% 
  column_to_rownames("cluster")

ann_col_neuron <- data.frame(
  cluster = c(neuron)
) %>% 
  mutate(cluster = str_replace(cluster, "-", "\\.")) %>% 
  mutate(origin = str_extract(cluster, "ctrl|poly")) %>% 
  column_to_rownames("cluster")

ann_row <- top_markers_noMT %>% 
  select(c(Gene.name,cluster)) %>% 
  mutate(DEG = "yes") %>% 
  spread(key= cluster, value = DEG, fill = "no") %>% 
  column_to_rownames("Gene.name")

NPC <- c(paste(c("g1","g3","g7","g9","g13","g16","g17","g18","g19","g22","g23","g27"), "ctrl", sep = "-"),
         paste(c("g1","g3","g7","g9","g13","g16","g17","g18","g19","g22","g23","g27"), "poly", sep = "-"))
OPC <- c(paste(c("g2","g4","g8","g10","g12","g14","g15","g24","g26"), "ctrl", sep = "-"),
         paste(c("g2","g4","g8","g10","g12","g14","g15","g24","g26"), "poly", sep = "-"))
neuron <-  c(paste(c("g5","g6","g11","g20","g25","g29","g30","g32"), "ctrl", sep = "-"),
             paste(c("g5","g6","g11","g20","g25","g29","g30","g32"), "poly", sep = "-"))

ann_row_NPC <- ann_row %>% 
  select(c(`cl-1`,`cl-3`,`cl-7`,`cl-9`,`cl-13`,`cl-16`,`cl-17`,`cl-18`,`cl-19`,`cl-22`,`cl-23`))

ann_row_OPC <- ann_row %>% 
  select(c(`cl-2`,`cl-4`,`cl-8`,`cl-10`,`cl-12`,`cl-15`,`cl-24`))

ann_row_neuron <- ann_row %>% 
  select(c(`cl-5`,`cl-6`,`cl-11`,`cl-20`,`cl-25`))

heat_col <- colorRampPalette(colors = c("darkblue","dodgerblue4", "white", "red", "darkred"))


## heatamaps of ctrl parts only to get col clustering
## NPC
NPC_log1p_zscore_ctrl <- pheatmap(
  log1p(pb_NPC[,grepl("ctrl", colnames(pb_NPC))]),
  annotation_col = ann_col_NPC,
  annotation_colors = ann_color,
  annotation_row = ann_row_NPC, 
  scale = "row",
  color = heat_col(700),
  border_color = NA,
  cellwidth = 9,
  cellheight = 9,
  main = "NPC DEG no MT log1p zscore"
)


NPC_ctrl <- NPC_log1p_zscore_ctrl$tree_col$labels[NPC_log1p_zscore_ctrl$tree_col$order]
NPC_poly <- str_replace(NPC_ctrl, "ctrl", "poly")
NPC_order <- c(rbind(NPC_ctrl, NPC_poly))

## OPC
OPC_log1p_zscore_ctrl <- pheatmap(
  log1p(pb_OPC[,grepl("ctrl", colnames(pb_OPC))]),
  annotation_col = ann_col_OPC,
  annotation_colors = ann_color,
  annotation_row = ann_row_OPC, 
  scale = "row",
  color = heat_col(700),
  border_color = NA,
  cellwidth = 9,
  cellheight = 9,
  main = "OPC DEG no MT log1p zscore"
)


OPC_ctrl <- OPC_log1p_zscore_ctrl$tree_col$labels[OPC_log1p_zscore_ctrl$tree_col$order]
OPC_poly <- str_replace(OPC_ctrl, "ctrl", "poly")
OPC_order <- c(rbind(OPC_ctrl, OPC_poly))

## neuron
neuron_log1p_zscore_ctrl <- pheatmap(
  log1p(pb_neuron[,grepl("ctrl", colnames(pb_neuron))]),
  annotation_col = ann_col_neuron,
  annotation_colors = ann_color,
  annotation_row = ann_row_neuron, 
  scale = "row",
  color = heat_col(700),
  border_color = NA,
  cellwidth = 9,
  cellheight = 9,
  main = "neuron DEG no MT log1p zscore"
)


neuron_ctrl <- neuron_log1p_zscore_ctrl$tree_col$labels[neuron_log1p_zscore_ctrl$tree_col$order]
neuron_poly <- str_replace(neuron_ctrl, "ctrl", "poly")
neuron_order <- c(rbind(neuron_ctrl, neuron_poly))


# z-score heatmap
NPC_log1p_zscore <- pheatmap(
  log1p(pb_NPC[,NPC_order]),
  annotation_col = ann_col_NPC,
  annotation_colors = ann_color,
  annotation_row = ann_row_NPC, 
  scale = "row",
  color = heat_col(700), 
  cluster_cols = FALSE,
  border_color = NA,
  cellwidth = 9,
  cellheight = 9,
  main = "NPC DEG no MT log1p zscore"
)
# z-score heatmap
OPC_log1p_zscore <- pheatmap(
  log1p(pb_OPC[,OPC_order]),
  annotation_col = ann_col_OPC, 
  annotation_colors = ann_color,
  annotation_row = ann_row_OPC,
  scale = "row",
  color = heat_col(700), 
  cluster_cols = FALSE,
  border_color = NA,
  cellwidth = 9,
  cellheight = 9,
  main = "OPC DEG no MT log1p zscore"
)
# z-score heatmap
neuron_log1p_zscore <- pheatmap(
  log1p(pb_neuron[,neuron_order]),
  annotation_col = ann_col_neuron, 
  annotation_colors = ann_color,
  annotation_row = ann_row_neuron,
  scale = "row",
  color = heat_col(700), 
  cluster_cols = FALSE,
  border_color = NA,
  cellwidth = 9,
  cellheight = 9,
  main = "neuron DEG no MT log1p zscore"
)

pdf("~/spinal_cord_paper/figures/Supp_Fig_5_DEG_heatmap_noMT.pdf", height = 10, width = 9)
gridExtra::grid.arrange(NPC_log1p_zscore_ctrl[[4]])
gridExtra::grid.arrange(NPC_log1p_zscore[[4]])
gridExtra::grid.arrange(OPC_log1p_zscore_ctrl[[4]])
gridExtra::grid.arrange(OPC_log1p_zscore[[4]])
gridExtra::grid.arrange(neuron_log1p_zscore_ctrl[[4]])
gridExtra::grid.arrange(neuron_log1p_zscore[[4]])
dev.off()

### ### ### ### ### ###
### Venn diagram DEG ####
### ### ### ### ### ###
library(tidyverse)
library(VennDiagram) 
library(eulerr)

top_markers_noMT <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_markers.rds") %>% 
  mutate(clust_id = str_remove(cluster, "^cl-")) %>% 
  mutate(clust_id = factor(clust_id, levels = c(1:32))) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(sample = case_when( # marker for B or poly?
    avg_log2FC > 0 ~ "ctrl",
    avg_log2FC < 0 ~ "poly"
  ))  %>% 
  mutate(DE_groups = factor(paste0("g", clust_id, "-", sample))) %>% 
  droplevels() %>% 
  mutate(abs_logFC = abs(avg_log2FC)) %>% 
  group_by(cluster) %>% 
  filter(!grepl("^ND\\d$|^MT-|^CYTB|^ATP|^COX\\d|^COII$", Gene.name)) %>% 
  mutate(Gene.name = str_replace(Gene.name, "^ENSGALG0+", "EG"))

mark_NPC_all_noMT <- top_markers_noMT[top_markers_noMT$DE_groups %in% NPC,] %>% 
  pull(Gene.name) %>% 
  unique()
mark_OPC_all_noMT <- top_markers_noMT[top_markers_noMT$DE_groups %in% OPC,] %>% 
  pull(Gene.name) %>% 
  unique()
mark_neuron_all_noMT <- top_markers_noMT[top_markers_noMT$DE_groups %in% neuron,] %>% 
  pull(Gene.name) %>% 
  unique()

venn.diagram(list(NPC = mark_NPC_all_noMT, 
                  OPC = mark_OPC_all_noMT,
                  neurons =mark_neuron_all_noMT),
             fill = c("#edc919",
                      "#008cb5",
                      "#cd2b91"), 
             alpha = c(0.5, 0.5, 0.5), lwd =0, 
             imagetype = "svg",
             "~/spinal_cord_paper/figures/Supp_fig_5_venn_diagram.svg")

venn <- list(NPC = mark_NPC_all_noMT, 
           OPC = mark_OPC_all_noMT,
           neurons =mark_neuron_all_noMT)

pdf("~/spinal_cord_paper/figures/Supp_fig_5_DEGvenn_diagram.pdf")
plot(euler(venn, shape = "ellipse"), quantities = TRUE)
dev.off()
