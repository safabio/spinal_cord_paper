### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 5 ####
## Fabio Sacher, 21.12.23
### ### ### ### ### ### ### ### ### ### ###


library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)


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

### ### ### ### ### ### 
#### DE heatmaps ####
### ### ### ### ### ### 

gnames <- modplots::gnames

# load the ctrl poly integrated data and add the combined labels.
my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_seurat_250723.rds")

ctrl_poly_int_combined_labels <- readRDS("~/spinal_cord_paper/annotations/ctrl_poly_int_combined_labels.rds")

my.se <- AddMetaData(my.se, ctrl_poly_int_combined_labels)

# subset for two 
Idents(my.se) <- "annot_sample"

sub <- subset(my.se, idents = c("excitatory neurons_11_ctrl", "excitatory_neurons_15_poly") )

cand_c11_p15 <- c("SPOCK1", "RUNX1T1", "TAC1", "RELN")

gnames[gnames$Gene.name %in% cand_c11_p15,]

c11_p15 <- DoHeatmap(
  subset(
    my.se, 
    idents = c("excitatory neurons_11_ctrl",
               "excitatory_neurons_15_poly")),
  raster = FALSE,
  features = gnames[gnames$Gene.name %in% cand_c11_p15, 1]
)

cand_c16_p14 <- c("ST18", "EPB41", "ASCL1", "NOTCH1")

gnames[gnames$Gene.name %in% cand_c16_p14,]

c16_p14 <- DoHeatmap(
  subset(
    my.se, 
    idents = c("inhibitory_neurons_16_ctrl",
               "inhibitory_neurons_14_poly")),
  raster = FALSE,
  features = gnames[gnames$Gene.name %in% cand_c16_p14, 1]
)

pdf("~/spinal_cord_paper/figures/Fig_5_ctrl_poly_DE_htmp.pdf")
c11_p15
c16_p14
dev.off()

