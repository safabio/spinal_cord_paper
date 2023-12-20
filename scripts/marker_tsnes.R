### ### ### ### ### ### ### ### ### ### ###
## Additional plots and code for figure 1
## Fabio Sacher, 19.12.23
### ### ### ### ### ### ### ### ### ### ###

library(Seurat)
library(ggplot2)
library(modplots)
library(gridExtra)

### ### ### ### ### ### ### ### ### ### ###
#### tSNE plots of RP candidate markers ####
### ### ### ### ### ### ### ### ### ### ###


data_sets <- c("~/spinal_cord_paper/data/Gg_D05_ctrl_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_D07_ctrl_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_ctrl_1_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

# Many of the genes from the lightgreen module
probes <- c("CHODL","SLIT3","SEMA3C","ALDH1A2","MECOM","LMX1A","RSPO1","MSX1", "MSX2")

gnames <- modplots::gnames

plots <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = probes, size = 0.5)
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/RP_tsne.pdf", width = 12, height = 10)
grid.arrange(plots[[1]])
grid.arrange(plots[[2]])
grid.arrange(plots[[3]])
grid.arrange(plots[[4]])
grid.arrange(plots[[5]])
grid.arrange(plots[[6]])
dev.off()

rm(probes, gnames, plots)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### tSNE plots of SLC10A4, SLC18A2, and A3 candidate markers ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


data_sets <- c("~/spinal_cord_paper/data/Gg_D05_ctrl_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_D07_ctrl_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_ctrl_1_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

probes <- c("SLC18A3","SLC10A4","SLC18A2")

gnames <- modplots::gnames

plots <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = probes, my.slot = "scale.data", size = 0.5)
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/SLC_tsne.pdf", width = 12, height = 10)
grid.arrange(plots[[1]])
grid.arrange(plots[[2]])
grid.arrange(plots[[3]])
grid.arrange(plots[[4]])
grid.arrange(plots[[5]])
grid.arrange(plots[[6]])
dev.off()

rm(probes, gnames, plots)

### ### ### ### ### ### ### ### ### ### ###
#### tSNE plots of MN column markers ####
### ### ### ### ### ### ### ### ### ### ###

# The goal is to see the extent of Lumbar MNs

data_sets <- c("~/spinal_cord_paper/data/Gg_D05_ctrl_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_D07_ctrl_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_ctrl_1_seurat_070323.rds",
               "~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

probes <-
  c("SLC18A3", # MN
    "CHAT", # MN
    "SLC10A4", # MN
    "ALDH1A2", # LMC
    "FOXP1", # LMC
    "TAC1", # HMC https://www.nature.com/articles/s41467-022-35574-x
    "LHX3", # MMC
    "ISL1", # MN
    "MNX1") # MN

gnames <- modplots::gnames

plots <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = probes, my.slot = "scale.data", size = 0.5)
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/MN_tsne.pdf", width = 12, height = 10)
grid.arrange(plots[[1]])
grid.arrange(plots[[2]])
grid.arrange(plots[[3]])
grid.arrange(plots[[4]])
grid.arrange(plots[[5]])
grid.arrange(plots[[6]])
dev.off()

rm(probes, gnames, plots)




