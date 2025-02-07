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

DimPlot(my.se,
        reduction = "tsne",
        group.by = "cond") + 
  plotNhoodGraphDA(my.milo, da.results, alpha=0.9, layout = "TSNE") +
  plot_layout(guides="collect")


nhg <- plotNhoodGraphDA(my.milo, da.results, alpha=0.9, layout = "TSNE") +
  theme_void() +
  NoLegend()

ggsave(filename = "~/spinal_cord_paper/figures/Fig_5_milo_network.pdf",
       width = 4,
       height = 4,
       nhg)

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


