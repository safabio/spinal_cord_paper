##############
## Additional plots and code for figure 4
## Fabio Sacher, 06.12.23
##############

library(Seurat)
library(ggplot2)
library(modplots)
library(dplyr)
library(ggvenn)
library(cowplot)
library(gridExtra)
library(pheatmap)

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

##############
## detailed venn diagram of CSF-module overlap
##############

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

##############
## tSNE plots of all 5 in situ probes
############## 


data_sets <- c("~/spinal_cord_paper/data/Gg_ctrl_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_lumb_int_seurat_250723.rds",
               "~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")

probes <- c("PKD2L1","GATA2","CACNA1G","SFRP5","CRTAC1")

gnames <- modplots::gnames

plots <- list()

for (i in seq(data.sets)) {
  my.se <- readRDS(data_sets[i])
  
  plots[[i]] <- mFeaturePlot(my.se, my.features = probes, my.slot = "scale.data", size = 1)
  
  rm(my.se)
}

pdf("~/spinal_cord_paper/figures/CSF_probes_tsne.pdf", width = 12, height = 7)
  grid.arrange(plots[[1]])
  grid.arrange(plots[[2]])
  grid.arrange(plots[[3]])
dev.off()

rm(data_sets, probes, gnames, plots)

##############
## CSF module membership values 
############## 

ctrl_wgcna <- readRDS("~/spinal_cord_paper/output/Gg_ctrl_int_scWGCNA_250723.rds")

probes <- c("PKD2L1","GATA2","CACNA1G","SFRP5","CRTAC1")

df <- ctrl_wgcna$modules[[17]] %>% 
  tibble::rownames_to_column("Gene.stable.ID") %>% 
  dplyr::left_join(gnames) %>% 
  arrange(desc(Membership)) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when(
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
  dplyr::left_join(gnames) %>% 
  arrange(desc(Membership)) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when(
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
  dplyr::left_join(gnames) %>% 
  arrange(desc(Membership)) %>% 
  mutate(Gene.name = factor(Gene.name, levels = Gene.name)) %>% 
  mutate(highlight = case_when(
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

##############
## pairwise gene expression correlation in Poly-cl26
############## 

probes <- c("PKD2L1", "GATA2", "CACNA1G", "SFRP5", "CRTAC1")
gnames <- modplots::gnames
rownames(gnames) <- gnames$Gene.name

probe_ids <- gnames[probes, ]

my.se <- readRDS("~/spinal_cord_paper/data/Gg_poly_int_seurat_250723.rds")  

cl26 <- subset(x = my.se, idents = 26)

probe_subset <- data.frame(t(as.matrix(cl26@assays[["RNA"]]@data[probe_ids$Gene.stable.ID,])))

# pearson correlation 
cor_table_pear <- cor(probe_subset, method = "pear")

if (identical(rownames(cor_table_pear),probe_ids$Gene.stable.ID) & identical(colnames(cor_table_pear),probe_ids$Gene.stable.ID)) {
  colnames(cor_table_pear) <- probe_ids$Gene.name
  rownames(cor_table_pear) <- probe_ids$Gene.name
}

# spearman correlation 
cor_table_spea <- cor(probe_subset, method = "spear")

if (identical(rownames(cor_table_spea),probe_ids$Gene.stable.ID) & identical(colnames(cor_table_spea),probe_ids$Gene.stable.ID)) {
  colnames(cor_table_spea) <- probe_ids$Gene.name
  rownames(cor_table_spea) <- probe_ids$Gene.name
}

pdf("~/spinal_cord_paper/figures/CSF_probes_poly_cl26_corr_heatmap.pdf", paper = "a4", width = 10, height = 12)
# pearson 
pheatmap(
  cor_table_pear,
  cellwidth = 50,
  cellheight = 50,
  display_numbers = TRUE,
  border_color = NA,
  main = "Pearson correlation of CSF probe expression in Poly cl-26"
)
# spearman
pheatmap(
  cor_table_spea,
  cellwidth = 50,
  cellheight = 50,
  display_numbers = TRUE,
  border_color = NA,
  main = "Spearman correlation of CSF probe expression in Poly cl-26"
)
dev.off()

