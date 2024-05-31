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

# cell highlight cl-26 poly
my.se <- readRDS("~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")  

p <- DimPlot(my.se,
             reduction = "tsne",
             cells.highlight = WhichCells(my.se, idents = "26"),
             cols.highlight = "black") +
  theme_minimal() +
  NoLegend()

pdf("~/spinal_cord_paper/figures/Gg_int_poly_cl26_highlight_tsne.pdf", width = 10, height = 10)
p
dev.off()

rm(my.se, p)

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
#### tSNE plots of all 5 in situ probes ####
### ### ### ### ### ### ### ### ### ### ###


data_sets <- c("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

# PKD2L1 is already in figure 4
probes <- c("PKD2L1","GATA2","CACNA1G","SFRP5","CRTAC1")

gnames <- modplots::gnames

plots <- list()

for (i in seq(data_sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = probes, my.slot = "scale.data", size = 0.5, return = TRUE)
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/CSF_probes_tsne.pdf", width = 8.5, height = 5)

  # brachial 
(plots[[1]][[1]] + theme(plot.title = element_blank()) +
    geom_segment(aes(x=-5,y=28,xend=-6,yend=29),
                 size = 0.8, color = "black", 
                 linejoin = "mitre", lineend = "square",
                 arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
    theme_void() +
    NoLegend() + annotate("text", label = probes[1], fontface = "italic", x = 0, y = 110)) +
  annotate("text", label = "A", fontface = "bold", x = -100, y = 110)  + 
  annotate("text", label = "brachial", x = -110, y = 0, angle = 90)  + 
  (plots[[1]][[2]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=-5,y=28,xend=-6,yend=29),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[2], fontface = "italic", x = 0, y = 110)) +  
  annotate("text", label = "brachial", col = "white", x = -110, y = 0, angle = 90)  + 
  (plots[[1]][[3]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=-5,y=28,xend=-6,yend=29),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[3], fontface = "italic", x = 0, y = 110)) + 
  annotate("text", label = "brachial", col = "white", x = -110, y = 0, angle = 90)  + 
  (plots[[1]][[4]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=-5,y=28,xend=-6,yend=29),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[4], fontface = "italic", x = 0, y = 110)) + 
  annotate("text", label = "brachial", col = "white", x = -110, y = 0, angle = 90)  + 
  (plots[[1]][[5]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=-5,y=28,xend=-6,yend=29),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[5], fontface = "italic", x = 0, y = 110)) +  
  annotate("text", label = "brachial", col = "white", x = -110, y = 0, angle = 90)  + 
  # lumbar
  (plots[[2]][[1]] + theme(plot.title = element_blank()) + 
     geom_segment(aes(x=40,y=18,xend=41,yend=19),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     
     theme_void() +
     NoLegend() + annotate("text", label = probes[1], fontface = "italic", x = 0, y = 110,
                           color = "white")) +
  annotate("text", label = "lumbar", x = -100, y = 0, angle = 90)  + 
  (plots[[2]][[2]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=40,y=18,xend=41,yend=19),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     
     theme_void() + 
     NoLegend() + annotate("text", label = probes[2], fontface = "italic", x = 0, y = 110,
                           color = "white")) + 
  annotate("text", label = "lumbar", col = "white", x = -100, y = 0, angle = 90)  + 
  (plots[[2]][[3]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=40,y=18,xend=41,yend=19),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     
     theme_void() + 
     NoLegend() + annotate("text", label = probes[3], fontface = "italic", x = 0, y = 110,
                           color = "white")) + 
  annotate("text", label = "lumbar", col = "white", x = -100, y = 0, angle = 90)  + 
  (plots[[2]][[4]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=40,y=18,xend=41,yend=19),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     
     theme_void() + 
     NoLegend() + annotate("text", label = probes[4], fontface = "italic", x = 0, y = 110,
                           color = "white")) + 
  annotate("text", label = "lumbar", col = "white", x = -100, y = 0, angle = 90)  + 
  (plots[[2]][[5]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=40,y=18,xend=41,yend=19),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[5], fontface = "italic", x = 0, y = 110,
                           color = "white")) + 
  annotate("text", label = "lumbar", col = "white", x = -100, y = 0, angle = 90)  + 
  # polydactyl
  (plots[[3]][[1]] +  theme(plot.title = element_blank())  +
     geom_segment(aes(x=-21,y=54,xend=-22,yend=53),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() +
     NoLegend() + annotate("text", label = probes[1], fontface = "italic", x = 0, y = 110,
                           color = "white")) + 
  annotate("text", label = "polydactyl", x = -100, y = 0, angle = 90)  + 
  (plots[[3]][[2]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=-21,y=54,xend=-22,yend=53),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[2], fontface = "italic", x = 0, y = 110,
                           color = "white")) + 
  annotate("text", label = "polydactyl", col = "white", x = -100, y = 0, angle = 90)  + 
  (plots[[3]][[3]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=-21,y=54,xend=-22,yend=53),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[3], fontface = "italic", x = 0, y = 110,
                           color = "white")) + 
  annotate("text", label = "polydactyl", col = "white", x = -100, y = 0, angle = 90)  + 
  (plots[[3]][[4]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=-21,y=54,xend=-22,yend=53),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[4], fontface = "italic", x = 0, y = 110,
                           color = "white")) + 
  annotate("text", label = "polydactyl", col = "white", x = -100, y = 0, angle = 90)  + 
  (plots[[3]][[5]] + theme(plot.title = element_blank()) +
     geom_segment(aes(x=-21,y=54,xend=-22,yend=53),
                  size = 0.8, color = "black", 
                  linejoin = "mitre", lineend = "square",
                  arrow = arrow(angle = 20, length = unit(0.25, "cm"), type = "closed")) +
     theme_void() + 
     NoLegend() + annotate("text", label = probes[5], fontface = "italic", x = 0, y = 110,
                           color = "white")) +
  annotate("text", label = "polydactyl", col = "white", x = -100, y = 0, angle = 90)  + 
  plot_layout(nrow = 3, byrow = TRUE) & theme(text = element_blank())

dev.off()

rm(data_sets, probes, gnames, plots)

### ### ### ### ### ### ### ### ### ###
#### CSF module membership values ####
### ### ### ### ### ### ### ### ### ### 

ctrl_wgcna <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds")

probes <- c("PKD2L1","GATA2","CACNA1G","SFRP5","CRTAC1")

df <- ctrl_wgcna$modules[[17]] %>% 
  tibble::rownames_to_column("Gene.stable.ID") %>% 
  dplyr::left_join(gnames) %>% # add gene names 
  arrange(desc(Membership)) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when( # add column with '*' indicating overlap genes
    Gene.name %in% probes ~ "*",
    TRUE ~ ""
  ))

p1 <- ggplot(df, aes(x = Gene.name,
                     y = Membership,
                     label = signif(Membership,3))) +
  geom_col(fill = "royalblue", color = "black") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  xlab("Module genes") +
  ylab("Module Membership") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) +
  geom_text(hjust = -0.2, angle = 90, cex = 3) +
  ggtitle("Module 17_royalblue gene Module Membership (ctrl_int)") +
  geom_text(aes(label = highlight), vjust = 1.5, cex = 8)

rm(ctrl_wgcna)

lumb_wgcna <- readRDS("~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds")

df <- lumb_wgcna$modules[[9]] %>% 
  tibble::rownames_to_column("Gene.stable.ID") %>% 
  dplyr::left_join(gnames) %>%  # add gene names
  arrange(desc(Membership)) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when( # add column with '*' indicating overlap genes
    Gene.name %in% probes ~ "*",
    TRUE ~ ""
  ))

p2 <- ggplot(df, aes(x = Gene.name,
                     y = Membership,
                     label = signif(Membership,3))) +
  geom_col(fill = "darkmagenta", color = "black") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  xlab("Module genes") +
  ylab("Module Membership") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) +
  geom_text(hjust = -0.2, angle = 90, cex = 3) +
  ggtitle("Module 9_darkmagenta gene Module Membership (lumb_int)") +
  geom_text(aes(label = highlight), vjust = 1.5, cex = 8)

rm(lumb_wgcna)

poly_wgcna <- readRDS("~/spinal_cord_paper/output/Gg_poly_int_scWGCNA_250723.rds")

df <- poly_wgcna$modules[[16]] %>% 
  tibble::rownames_to_column("Gene.stable.ID") %>% 
  dplyr::left_join(gnames) %>%  # add gene names
  arrange(desc(Membership)) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when( # add column with '*' indicating overlap genes
    Gene.name %in% probes ~ "*",
    TRUE ~ ""
  ))

p3 <- ggplot(df, aes(x = Gene.name,
                     y = Membership,
                     label = signif(Membership,3))) +
  geom_col(fill = "salmon", color = "black") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  xlab("Module genes") +
  ylab("Module Membership") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) +
  geom_text(hjust = -0.2, angle = 90, cex = 3) +
  ggtitle("Module 16_salmon gene Module Membership (poly_int)") +
  geom_text(aes(label = highlight), vjust = 1.5, cex = 8)

rm(poly_wgcna)

pdf("~/spinal_cord_paper/figures/CSF_module_memberships.pdf", paper = "a4", width = 10, height = 12)
p1 / p2 / p3
dev.off()

rm(df, p1, p2, p3)

### ### ### ### ### ### ### ### ### ### ###
#### pairwise gene expression correlation in Poly_int ####
### ### ### ### ### ### ### ### ### ### ###

overlap <- c("PKD2L1", "GATA2", "CACNA1G", "SFRP5", "CRTAC1", "ABCG2", "GATA3",
             "MYO3B", "ST3GAL1", "SST", "TAL1", "TAL2", "ENSGALG00000000645", "ENSGALG00000007596")
gnames <- modplots::gnames
rownames(gnames) <- gnames$Gene.name

overlap_ids <- gnames[overlap, ]

my.se <- readRDS("~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")  

overlap_subset <- data.frame(t(as.matrix(my.se@assays[["RNA"]]@data[overlap_ids$Gene.stable.ID,])))

# pearson correlation 
cor_table_pear <- cor(overlap_subset, method = "pear")
# replace gene IDs with gene names
if (identical(rownames(cor_table_pear),overlap_ids$Gene.stable.ID) & identical(colnames(cor_table_pear),overlap_ids$Gene.stable.ID)) {
  colnames(cor_table_pear) <- overlap_ids$Gene.name
  rownames(cor_table_pear) <- overlap_ids$Gene.name
}

# spearman correlation 
cor_table_spea <- cor(overlap_subset, method = "spear")
# replace gene IDs with gene names
if (identical(rownames(cor_table_spea),overlap_ids$Gene.stable.ID) & identical(colnames(cor_table_spea),overlap_ids$Gene.stable.ID)) {
  colnames(cor_table_spea) <- overlap_ids$Gene.name
  rownames(cor_table_spea) <- overlap_ids$Gene.name
}

pdf("~/spinal_cord_paper/figures/CSF_probes_poly_int_corr_heatmap.pdf", paper = "a4", width = 10, height = 12)
# pearson 
pheatmap(
  cor_table_pear,
  cellwidth = 20,
  cellheight = 20,
  display_numbers = TRUE,
  border_color = NA,
  main = "Pearson correlation of CSF probe expression in Poly int"
)
# spearman
pheatmap(
  cor_table_spea,
  cellwidth = 20,
  cellheight = 20,
  display_numbers = TRUE,
  border_color = NA,
  main = "Spearman correlation of CSF probe expression in Poly int"
)
dev.off()

### ### ### ### ### ### ### ### ### ### ###
#### pairwise gene expression correlation in Poly-cl26 ####
### ### ### ### ### ### ### ### ### ### ###

overlap <- c("PKD2L1", "GATA2", "CACNA1G", "SFRP5", "CRTAC1", "ABCG2", "GATA3",
             "MYO3B", "ST3GAL1", "SST", "TAL1", "TAL2", "ENSGALG00000000645", "ENSGALG00000007596")
gnames <- modplots::gnames
rownames(gnames) <- gnames$Gene.name

overlap_ids <- gnames[overlap, ]

my.se <- readRDS("~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")  

cl26 <- subset(x = my.se, idents = 26)

overlap_subset <- data.frame(t(as.matrix(cl26@assays[["RNA"]]@data[overlap_ids$Gene.stable.ID,])))

# pearson correlation 
cor_table_pear <- cor(overlap_subset, method = "pear")
# replace gene IDs with gene names
if (identical(rownames(cor_table_pear),overlap_ids$Gene.stable.ID) & identical(colnames(cor_table_pear),overlap_ids$Gene.stable.ID)) {
  colnames(cor_table_pear) <- overlap_ids$Gene.name
  rownames(cor_table_pear) <- overlap_ids$Gene.name
}

# spearman correlation 
cor_table_spea <- cor(overlap_subset, method = "spear")
# replace gene IDs with gene names
if (identical(rownames(cor_table_spea),overlap_ids$Gene.stable.ID) & identical(colnames(cor_table_spea),overlap_ids$Gene.stable.ID)) {
  colnames(cor_table_spea) <- overlap_ids$Gene.name
  rownames(cor_table_spea) <- overlap_ids$Gene.name
}

pdf("~/spinal_cord_paper/figures/CSF_probes_poly_cl26_corr_heatmap.pdf", paper = "a4", width = 10, height = 12)
# pearson 
pheatmap(
  cor_table_pear,
  cellwidth = 20,
  cellheight = 20,
  display_numbers = TRUE,
  border_color = NA,
  main = "Pearson correlation of CSF probe expression in Poly cl-26"
)
# spearman
pheatmap(
  cor_table_spea,
  cellwidth = 20,
  cellheight = 20,
  display_numbers = TRUE,
  border_color = NA,
  main = "Spearman correlation of CSF probe expression in Poly cl-26"
)
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

### ### ### ### ### ###
#### NRP2 in WGCNA ####
### ### ### ### ### ###

wgcna_dat <- c(
  "~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds",
  "~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds",
  "~/spinal_cord_paper/output/Gg_poly_int_scWGCNA_250723.rds"
  )

nrp2_id <- gnames$Gene.stable.ID[gnames$Gene.name %in% "NRP2"]

output <- list()

for (i in seq(wgcna_dat)) {
  my.wgcna <- readRDS(wgcna_dat[3])
  
  output[[i]] <- nrp2_id %in% head(unlist(my.wgcna[["module.genes"]]))
}

print(output)
# [[1]]
# [1] FALSE
# 
# [[2]]
# [1] FALSE
# 
# [[3]]
# [1] FALSE

### ### ### ### ### ### ### ### ### ### ###
#### cl-26 poly module salmon expression ####
### ### ### ### ### ### ### ### ### ### ###

my.se <- readRDS("~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

my.wgcna <- readRDS("~/spinal_cord_paper/output/Gg_poly_int_scWGCNA_250723.rds")

avgExp <- my.wgcna[["sc.MEList"]]$averageExpr

if (identical(rownames(my.se[[]]), rownames(avgExp))) {
  my.se <- AddMetaData(my.se, avgExp)
}
# cluster color and marker table
clust_col <- read.csv("~/spinal_cord_paper/annotations/broad_cluster_marker_colors.csv")
# cluster annotation
annot <- read.csv("~/spinal_cord_paper/annotations/Gg_poly_int_cluster_annotation.csv")

broad_order <- c("progenitors",
                 "FP",
                 "RP",
                 "FP/RP",
                 "neurons",
                 "OPC",
                 "MFOL",
                 "pericytes",
                 "microglia",
                 "blood",
                 "vasculature"
)

# rename for left join
annot <- annot %>% 
  mutate(fine = paste(fine, number, sep = "_")) %>% 
  mutate(number = factor(number, levels = 1:nrow(annot))) %>% 
  rename(seurat_clusters = number)

ord_levels <- annot$fine[order(match(annot$broad, broad_order))]

# add cluster annotation to meta data
my.se@meta.data <- my.se@meta.data %>% 
  rownames_to_column("rowname") %>% 
  left_join(annot, by = "seurat_clusters") %>% 
  mutate(fine = factor(fine, levels = ord_levels)) %>% 
  mutate(seurat_clusters = factor(seurat_clusters, levels = str_extract(ord_levels, "\\d{1,2}$"))) %>% 
  column_to_rownames("rowname")
# color palette
col_plot <- c(rep("#edc919", 8),
              "#853335",
              rep("#99ca3c", 2),
              rep("#cd2b91", 8),
              rep("#008cb5", 6),
              "#333399","#f26a2d",
              "#000000","#bebebe","#996633")

# violin plot of module salmon expression by cluster
pdf("~/spinal_cord_paper/figures/poly_module_salmon_vln_plot.pdf", width = 10, height = 6)
  VlnPlot(my.se, "AEsalmon", group.by = "seurat_clusters", cols = col_plot) +
    ggtitle("AEsalmon avgExp by cell across poly_int clusters")
  
  VlnPlot(my.se, "AEsalmon", group.by = "seurat_clusters", cols = col_plot, pt.size = 0) +
    geom_boxplot() +
    ggtitle("AEsalmon avgExp by cell across poly_int clusters with boxplot")
  
  VlnPlot(my.se, "AEsalmon", group.by = "seurat_clusters", cols = col_plot, pt.size = 0.5) +
    scale_y_log10() +
    ggtitle("AEsalmon avgExp by cell across poly_int clusters (log scale, 0 removed")
dev.off()

# 3D plot of module salmon expression
library(plotly)
library(htmlwidgets)

x_tsne <- my.se@reductions[["tsne"]]@cell.embeddings[,"tsne_1"]
y_tsne <- my.se@reductions[["tsne"]]@cell.embeddings[,"tsne_2"]
salmon <- my.se[[]]$AEsalmon

saveWidget(
  plot_ly(
    x = x_tsne,
    y = y_tsne,
    z = salmon,
    type = "scatter3d",
    mode = "markers",
    color = salmon,
    size = 1,
    alpha_stroke = 1,
    colors = c("#E5E5E5", "#FA8072")
  ),
  file = "~/spinal_cord_paper/figures/poly_mod_salmon_3d.html")

### ### ### ### ### ### ### ### ### ### ###
#### lumb module darkmagent expression ####
### ### ### ### ### ### ### ### ### ### ###

my.se <- readRDS("~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds")

my.wgcna <- readRDS("~/spinal_cord_paper/output/Gg_lumb_int_scWGCNA_250723.rds")

avgExp <- my.wgcna[["sc.MEList"]]$averageExpr

if (identical(rownames(my.se[[]]), rownames(avgExp))) {
  my.se <- AddMetaData(my.se, avgExp)
}
# cluster color and marker table
clust_col <- read.csv("~/spinal_cord_paper/annotations/broad_cluster_marker_colors.csv")
# cluster annotation
annot <- read.csv("~/spinal_cord_paper/annotations/Gg_lumb_int_cluster_annotation.csv")

broad_order <- c("progenitors",
                 "FP",
                 "RP",
                 "FP/RP",
                 "neurons",
                 "OPC",
                 "MFOL",
                 "pericytes",
                 "microglia",
                 "blood",
                 "vasculature"
)

# rename for left join
annot <- annot %>% 
  mutate(fine = paste(fine, number, sep = "_")) %>% 
  mutate(number = factor(number, levels = 1:nrow(annot))) %>% 
  rename(seurat_clusters = number)

ord_levels <- annot$fine[order(match(annot$broad, broad_order))]

# add cluster annotation to meta data
my.se@meta.data <- my.se@meta.data %>% 
  rownames_to_column("rowname") %>% 
  left_join(annot, by = "seurat_clusters") %>% 
  mutate(fine = factor(fine, levels = ord_levels)) %>% 
  mutate(seurat_clusters = factor(seurat_clusters, levels = str_extract(ord_levels, "\\d{1,2}$"))) %>% 
  column_to_rownames("rowname")
# color palette
col_plot <- c(rep("#edc919", 9),
              "#853335",
              rep("#99ca3c", 2),
              rep("#cd2b91", 7),
              rep("#008cb5", 5),
              "#333399","#f26a2d",
              "#000000","#bebebe","#996633")

# violin plot of module salmon expression by cluster
pdf("~/spinal_cord_paper/figures/lumb_module_darkmagenta_vln_plot.pdf", width = 10, height = 6)
VlnPlot(my.se, "AEdarkmagenta", group.by = "seurat_clusters", cols = col_plot) +
  ggtitle("AEdarkmagenta avgExp by cell across lumb_int clusters")

VlnPlot(my.se, "AEdarkmagenta", group.by = "seurat_clusters", cols = col_plot, pt.size = 0) +
  geom_boxplot() +
  ggtitle("AEdarkmagenta avgExp by cell across lumb_int clusters with boxplot")

VlnPlot(my.se, "AEdarkmagenta", group.by = "seurat_clusters", cols = col_plot, pt.size = 0.5) +
  scale_y_log10() +
  ggtitle("AEdarkmagenta avgExp by cell across lumb_int clusters (log scale, 0 removed)")
dev.off()

### ### ### ### ### ### ### ### ### ### ###
#### ctrl module royalblue expression ####
### ### ### ### ### ### ### ### ### ### ###

my.se <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds")

my.wgcna <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds")

avgExp <- my.wgcna[["sc.MEList"]]$averageExpr

if (identical(rownames(my.se[[]]), rownames(avgExp))) {
  my.se <- AddMetaData(my.se, avgExp)
}
# cluster color and marker table
clust_col <- read.csv("~/spinal_cord_paper/annotations/broad_cluster_marker_colors.csv")
# cluster annotation
annot <- read.csv("~/spinal_cord_paper/annotations/Gg_ctrl_int_cluster_annotation.csv")

broad_order <- c("progenitors",
                 "FP",
                 "RP",
                 "FP/RP",
                 "neurons",
                 "OPC",
                 "MFOL",
                 "pericytes",
                 "microglia",
                 "blood",
                 "vasculature"
)

# rename for left join
annot <- annot %>% 
  mutate(fine = paste(fine, number, sep = "_")) %>% 
  mutate(number = factor(number, levels = 1:nrow(annot))) %>% 
  rename(seurat_clusters = number)

ord_levels <- annot$fine[order(match(annot$broad, broad_order))]

# add cluster annotation to meta data
my.se@meta.data <- my.se@meta.data %>% 
  rownames_to_column("rowname") %>% 
  left_join(annot, by = "seurat_clusters") %>% 
  mutate(fine = factor(fine, levels = ord_levels)) %>% 
  mutate(seurat_clusters = factor(seurat_clusters, levels = str_extract(ord_levels, "\\d{1,2}$"))) %>% 
  column_to_rownames("rowname")
# color palette
col_plot <- c(rep("#edc919", 6),
              "#00ff00",
              rep("#cd2b91", 5),
              rep("#008cb5", 6),
              rep("#333399", 2),
              "#f26a2d","#bebebe","#996633")

# violin plot of module salmon expression by cluster
pdf("~/spinal_cord_paper/figures/ctrl_module_royalblue_vln_plot.pdf", width = 10, height = 6)
VlnPlot(my.se, "AEroyalblue", group.by = "seurat_clusters", cols = col_plot) +
  ggtitle("AEroyalblue avgExp by cell across ctrl_int clusters")

VlnPlot(my.se, "AEroyalblue", group.by = "seurat_clusters", cols = col_plot, pt.size = 0) +
  geom_boxplot() +
  ggtitle("AEroyalblue avgExp by cell across ctrl_int clusters with boxplot")

VlnPlot(my.se, "AEroyalblue", group.by = "seurat_clusters", cols = col_plot, pt.size = 0.5) +
  scale_y_log10() +
  ggtitle("AEroyalblue avgExp by cell across ctrl_int clusters (log scale, 0 removed)")
dev.off()


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




