### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 3 ####
## Fabio Sacher, 06.12.23
### ### ### ### ### ### ### ### ### ### ###

library(Seurat)
library(ggplot2)
library(modplots)
library(dplyr)
library(tibble)
library(stringr)
library(ggvenn)
library(cowplot)
library(gridExtra)
library(pheatmap)
library(patchwork)
library(scWGCNA)

### ### ### ### ### ### ### ### ### ### ###
#### detailed venn diagram of CSF-module overlap ####
### ### ### ### ### ### ### ### ### ### ###

wgcna_dat <- list(
  ctrl = readRDS("~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds"),
  lumb = readRDS("~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds"),
  poly = readRDS("~/spinal_cord_paper/output/Gg_poly_int_scWGCNA_250723.rds")
)

# get the module genes
ctrl_royalb <- wgcna_dat[[1]]$module.genes[[17]]
lumb_darkma <- wgcna_dat[[2]]$module.genes[[9]]
poly_salmon <- wgcna_dat[[3]]$module.genes[[16]]

ctrl_names <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% ctrl_royalb) 

lumb_names <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% lumb_darkma) 

poly_names <- modplots::gnames %>% 
  filter(Gene.stable.ID %in% poly_salmon) 


venn = list(ctrl_royalb = ctrl_names$Gene.name,
            lumb_darkma = lumb_names$Gene.name, 
            poly_salmon = poly_names$Gene.name)

p_lab <- ggvenn(venn, fill_color = c("royalblue", "darkmagenta", "salmon"), show_elements = TRUE)
p <- ggvenn(venn, fill_color = c("royalblue", "darkmagenta", "salmon"))

pdf("~/spinal_cord_paper/figures/CSF-cNS_modules_venn.pdf", width = 10, height = 10)
p_lab
p
dev.off()

Reduce(intersect, venn)

rm(wgcna_dat)

### ### ### ### ### ### ### ### ### ### ###
#### tSNE plots of new 3 in situ probes ####
### ### ### ### ### ### ### ### ### ### ###

data_sets <- c("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

# PKD2L1 is already in figure 4
probes <- c("CACNA1G","SFRP5","CRTAC1")

gnames <- modplots::gnames

plots <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(
    my.se, 
    my.features = probes, 
    my.slot = "scale.data", 
    size = 0.5, 
    raster = TRUE,
    return = TRUE)
  
  rm(my.se)
}

### ### ### ### ### ### ### ### ### ###
#### CSF module membership values ####
### ### ### ### ### ### ### ### ### ### 

ctrl_wgcna <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds")

stars <- c("PKD2L1","GATA2","CACNA1G","SFRP5","CRTAC1")

text_size = 7

df <- ctrl_wgcna$modules[[17]] %>% 
  tibble::rownames_to_column("Gene.stable.ID") %>% 
  dplyr::left_join(gnames) %>% # add gene names 
  arrange(Membership) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when( # add column with '*' indicating overlap genes
    Gene.name %in% stars ~ "*",
    TRUE ~ ""
  ))

p1 <- ggplot(df, aes(x = Gene.name,
                     y = Membership,
                     label = signif(Membership,3))) +
  geom_col(fill = "royalblue") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  xlab("Module genes") +
  ylab("Module Membership") +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = text_size),
        axis.text.y = element_text(size = text_size),
        axis.title.x = element_text(size = text_size + 2, color = "white"),
        axis.title.y = element_text(size = text_size + 2),
        plot.title = element_text(
          size = text_size + 4,
          face = "plain",
          hjust = -0.5,
          color = "royalblue"
        )) +
  geom_text(hjust = -0.2, cex = 2) +
  ggtitle("royalblue") +
  geom_text(aes(label = highlight), hjust = 1.5, vjust = 0.75, cex = 8, color = "white") +
  coord_flip()

lumb_wgcna <- readRDS("~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds")

df <- lumb_wgcna$modules[[9]] %>% 
  tibble::rownames_to_column("Gene.stable.ID") %>% 
  dplyr::left_join(gnames) %>%  # add gene names
  arrange(Membership) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when( # add column with '*' indicating overlap genes
    Gene.name %in% stars ~ "*",
    TRUE ~ ""
  ))

p2 <- ggplot(df, aes(x = Gene.name,
                     y = Membership,
                     label = signif(Membership,3))) +
  geom_col(fill = "darkmagenta") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  xlab("Module genes") +
  ylab("Module Membership                             ") +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = text_size),
        axis.text.y = element_text(size = text_size),
        axis.title.x = element_text(size = text_size + 2),
        axis.title.y = element_text(size = text_size + 2, color = "white"),
        plot.title = element_text(
          size = text_size + 4,
          face = "plain",
          hjust = -2,
          color = "darkmagenta"
        )) +
  geom_text(hjust = -0.2, cex = 2) +
  ggtitle("darkmagenta") +
  geom_text(aes(label = highlight), hjust = 1.5, vjust = 0.75, cex = 5, color = "white") +
  coord_flip()

poly_wgcna <- readRDS("~/spinal_cord_paper/output/Gg_poly_int_scWGCNA_250723.rds")

df <- poly_wgcna$modules[[16]] %>% 
  tibble::rownames_to_column("Gene.stable.ID") %>% 
  dplyr::left_join(gnames) %>%  # add gene names
  arrange(Membership) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when( # add column with '*' indicating overlap genes
    Gene.name %in% stars ~ "*",
    TRUE ~ ""
  ))

p3 <- ggplot(df, aes(x = Gene.name,
                     y = Membership,
                     label = signif(Membership,3))) +
  geom_col(fill = "salmon") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  xlab("Module genes") +
  ylab("Module Membership") +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = text_size),
        axis.text.y = element_text(size = text_size),
        axis.title.x = element_text(size = text_size + 2, color = "white"),
        axis.title.y = element_text(size = text_size + 2, color = "white"),
        plot.title = element_text(
          size = text_size + 4,
          face = "plain",
          hjust = -0.3,
          color = "salmon"
        )) +
  geom_text(hjust = -0.2, cex = 2) +
  ggtitle("salmon") +
  geom_text(aes(label = highlight), hjust = 1.5, vjust = 0.75, cex = 4) +
  coord_flip()

### ### ### ### ### ###
#### Supp figure 3 ####
### ### ### ### ### ##

t11 <- plots[[1]][[1]] + 
  theme_void() + 
  theme(plot.title = element_blank()) +
  NoLegend() +
  annotate("text", label = "B10_int", x = 0, y = 110) +
  annotate("text", label = "A", fontface = "bold", x = -110, y = 110) + 
  annotate("text", label = probes[1], fontface = "italic", x = -110, y = 0, angle = 90)

t21 <- plots[[2]][[1]] + 
  theme_void() +
  theme(plot.title = element_blank()) + 
  NoLegend() + 
  annotate("text", label = "L10_int", x = 0, y = 110) + 
  annotate("text", label = probes[1], fontface = "italic", x = -110, y = 0, angle = 90, color = "white")

t31 <- plots[[3]][[1]] +
  theme_void() +
  theme(plot.title = element_blank())  +
  NoLegend() + 
  annotate("text", label = "P10_int", x = 0, y = 110) + 
  annotate("text", label = probes[1], fontface = "italic", x = -110, y = 0, angle = 90, color = "white")
# SFRP5
t12 <- plots[[1]][[2]] +
  theme_void() + 
  theme(plot.title = element_blank()) +
  NoLegend() + 
  annotate("text", label = "brachial", x = 0, y = 110, color = "white") +
  annotate("text", label = probes[2], fontface = "italic", x = -110, y = 0, angle = 90)

t22 <- plots[[2]][[2]] +
  theme_void() + 
  theme(plot.title = element_blank()) +
  annotate("text", label = "lumbar", x = 0, y = 110, color = "white") +
  annotate("text", label = probes[2], fontface = "italic", x = -110, y = 0, angle = 90, color = "white") +
  NoLegend()

t32 <- plots[[3]][[2]] +
  theme_void() + 
  theme(plot.title = element_blank()) +
  annotate("text", label = "polydactyl", x = 0, y = 110, color = "white") + 
  annotate("text", label = probes[2], fontface = "italic", x = -110, y = 0, angle = 90, color = "white") +
  NoLegend() 
# CRTAC1
t13 <- plots[[1]][[3]] +
  theme_void() + 
  theme(plot.title = element_blank()) +
  NoLegend() + 
  annotate("text", label = "brachial", x = 0, y = 110, color = "white") +
  annotate("text", label = probes[3], fontface = "italic", x = -110, y = 0, angle = 90)

t23 <- plots[[2]][[3]] +
  theme_void() + 
  annotate("text", label = "lumbar", x = 0, y = 110, color = "white") +
  theme(plot.title = element_blank()) + 
  annotate("text", label = probes[3], fontface = "italic", x = -110, y = 0, angle = 90, color = "white") +
  NoLegend()

t33 <- plots[[3]][[3]] +
  theme_void() + 
  annotate("text", label = "polydactyl", x = 0, y = 110, color = "white") +
  theme(plot.title = element_blank()) + 
  annotate("text", label = probes[3], fontface = "italic", x = -110, y = 0, angle = 90, color = "white") +
  NoLegend()

pdf("~/spinal_cord_paper/figures/Supp_Fig_3.pdf", paper = "a4", width = 8.5, height = 11)
grid.arrange(t11, t21, t31,
             t12, t22, t32,
             t13, t23, t32,
             p1, p2, p3, heights=c(1,1,1,2))
dev.off()

### ### ### ### ### ### ### ### ### ### ###
#### tSNE plots of all CSF-cNS modules overlap ####
### ### ### ### ### ### ### ### ### ### ###

data_sets <- c("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

overlap <- c("PKD2L1", "GATA2", "CACNA1G", "SFRP5", "CRTAC1", "ABCG2", "GATA3",
             "MYO3B", "ST3GAL1", "SST", "TAL1", "TAL2", "ENSGALG00000000645", "ENSGALG00000007596")
gnames <- modplots::gnames

plots <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = overlap, my.slot = "scale.data", size = 0.5, return = TRUE)
  # set empty legend title
  for (j in seq(overlap)) {
    plots[[i]][[j]][["labels"]][["colour"]] <- ""
  }
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/CSF_overlap_tsne.pdf", width = 18, height = 14)
grid.arrange(grobs = plots[[1]])
grid.arrange(grobs = plots[[2]])
grid.arrange(grobs = plots[[3]])
dev.off()

rm(data_sets, overlap, gnames, plots)


### ### ### ### ### ### ### ### ### ### ###
#### tSNE plots of NRP2 (neuropilin 2) ####
### ### ### ### ### ### ### ### ### ### ###


data_sets <- c("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

blend_names <- modplots::gnames %>% 
  filter(Gene.name %in% c("PKD2L1", "NRP2")) 

plots <- list()
blend <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = c("PKD2L1","NRP2"), my.slot = "scale.data", size = 0.5, return = TRUE)
  tmp <- FeaturePlot(
    my.se,
    reduction = "tsne",
    features = blend_names$Gene.stable.ID,
    blend = TRUE,
    cols = "white",
    order = TRUE, 
    combine = FALSE
  )
  # create patchwork
  blend[[i]] <- (tmp[[1]] + tmp[[3]])/
    (tmp[[2]] + tmp[[4]]) +
    plot_layout(guides = 'collect')
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/NRP2_tsne.pdf", width = 6, height = 8)
(plots[[1]][[1]]  + theme(plot.title = element_blank()) + theme_void() + 
    NoLegend() + annotate("text", label = paste(probes, collapse = '_'), x = 0, y = 0)) +
  (plots[[1]][[2]]  + theme(plot.title = element_blank()) + theme_void() + 
     NoLegend() + annotate("text", label = paste(probes, collapse = '_'), x = 0, y = 0)) +
  (plots[[2]][[1]]  + theme(plot.title = element_blank()) + theme_void() + 
     NoLegend() + annotate("text", label = paste(probes, collapse = '_'), x = 0, y = 0)) +
  (plots[[2]][[2]]  + theme(plot.title = element_blank()) + theme_void() + 
     NoLegend() + annotate("text", label = paste(probes, collapse = '_'), x = 0, y = 0)) +
  (plots[[3]][[1]]  + theme(plot.title = element_blank()) + theme_void() + 
     NoLegend() + annotate("text", label = paste(probes, collapse = '_'), x = 0, y = 0)) +
  (plots[[3]][[2]]  + theme(plot.title = element_blank()) + theme_void() + 
     NoLegend() + annotate("text", label = paste(probes, collapse = '_'), x = 0, y = 0)) +
  plot_layout(nrow = 3, byrow = TRUE) & theme(text = element_blank())
dev.off()

pdf("~/spinal_cord_paper/figures/NRP2_PKD2L1_blend_tsne.pdf", width = 7, height = 6)
blend[[1]]
blend[[2]]
blend[[3]]
dev.off()

rm(data_sets, plots)

### ### ### ### ### ### ### ### ### ### ###
#### CSF modules network plot  ####
### ### ### ### ### ### ### ### ### ### ###

# code is taken and adapted from scWGCNA::scW.p.network

# scWGNA.data
wgcna_dat <- list(
  ctrl = readRDS("~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds"),
  lumb = readRDS("~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds"),
  poly = readRDS("~/spinal_cord_paper/output/Gg_poly_int_scWGCNA_250723.rds")
)

modules <- c("royalblue", "darkmagenta", "salmon")

gnames <- modplots::gnames
rownames(gnames) <- gnames[,1]

for (i in seq(modules)) {
  # get the module colors
  my.cols <- levels(as.factor(wgcna_dat[[i]][["dynamicCols"]]))
  
  col <- modules[i]
  # get index if color is provided
  module <- which(my.cols == col)
  # extract the module
  mynet <- wgcna_dat[[i]][["networks"]][[module]]
  
  gene.labs = gnames[network::network.vertex.names(mynet),2]
  set.seed(42)
  
  GGally::ggnet2(mynet,
                 mode = "fruchtermanreingold",
                 layout.par = list(repulse.rad = network::network.size(mynet)^1.1,
                                   area = network::network.size(mynet)^2.3),
                 node.size = network::get.vertex.attribute(mynet, "membership01"),
                 max_size = 23,
                 node.color = col, 
                 edge.size = "weight02", 
                 edge.color = "black",
                 edge.alpha = network::get.edge.attribute(mynet, "weight01")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_text(ggplot2::aes(label = gene.labs),
                       alpha = 1,
                       color = "black",
                       fontface = "italic",
                       size = 5
    )
  
  ggsave(paste0("~/spinal_cord_paper/figures/Fig_3_",
                names(wgcna_dat)[i], 
                "_int_", 
                col,
                "_network.pdf"),
         width = 6, height = 6)
}

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


