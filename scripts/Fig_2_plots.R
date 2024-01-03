### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 2 ####
## Fabio Sacher, 03.01.2024
### ### ### ### ### ### ### ### ### ### ###


library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

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



