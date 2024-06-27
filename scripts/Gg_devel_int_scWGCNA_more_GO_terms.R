
library(scWGCNA)
library(Seurat)
library(gridExtra)
library(ggplot2)
library(magrittr)

load("~/spinal_cord_paper/data/gga_KEGG_gene_pathway.Rda")
load("~/spinal_cord_paper/data/gga_KEGG_path_name.Rda")

# load seurat object
my.se <- readRDS("~/spinal_cord_paper/data/Gg_devel_int_seurat_250723.rds")
my.scWGCNA <- readRDS("~/spinal_cord_paper/output/Gg_devel_int_scWGCNA_250723.rds") 

DefaultAssay(my.se) <- "SCT"
my.se <- RunPCA(my.se,
                features = rownames(my.se[["SCT"]]@scale.data),
                verbose = FALSE)

my.se

# set variable features for SCT assay
my.se[["SCT"]]@var.features <- my.se[["integrated"]]@var.features

# set reduction for plots 
my.reduction <- "tsne"

my.sp <-  "Gg"

# var features for scWGCNA
var.feat <- my.se[["SCT"]]@var.features

# gnames (important: rownames == IDs)
my.gnames <- modplots::gnames
rownames(my.gnames) <- my.gnames$Gene.stable.ID

my.gouni = rownames(my.se[["RNA"]]@counts)[which(apply(my.se[["RNA"]]@counts, 1, function(x) length(which(x >0))) > 0)]

#Get the ENSEMBL to ENTREZ list from the BioMart
library(org.Gg.eg.db)
IDs=as.list(org.Gg.egENSEMBL2EG)

#Create two lists, one for the DEGs and one for the GOs
dewg = list()
gowg = list()
kewg = list()


# A loop that goes through the DEG files we created earlier
for(i in seq(length(names(my.scWGCNA$modules)))){ #as many clusters as we have
  dewg[[i]] = my.scWGCNA$module.genes[[i]] #Read the table
  ex = sapply(dewg[[i]], function(x) AnnotationDbi::exists(x, org.Gg.egENSEMBL2EG)) #Which genes have an ENTREZ?
  dewg[[i]] = dewg[[i]][ex] #Only those genes
  dewg[[i]] = unlist(IDs[dewg[[i]]]) #The ENTREZ IDs
  gowg[[i]] = limma::goana(dewg[[i]], species = my.sp, universe = unlist(IDs[my.gouni])) #The GO analysis
  kewg[[i]] = limma::kegga(dewg[[i]], species = my.sp, universe = unlist(IDs[my.gouni]), gene.pathway = GeneID.PathID, pathway.names = PathID.PathName)
}

keggpath=list()
goterms=list()

for (i in seq(length(names(my.scWGCNA$modules)))) {
  #### TopGO terms ####
  gowg_filt <- gowg[[i]] %>% dplyr::filter(P.DE < 0.05)
  toplot=limma::topGO(gowg_filt, n= Inf, ontology = c("BP"))
  
  #Change to character and then back to factor, to keep the order from TopGO
  toplot$Term = as.character(toplot$Term)
  toplot$Term = factor(toplot$Term, levels = unique(toplot$Term))
  
  goterms[[i]]=toplot
  
  rm(toplot)
  
  #### TopKEGG pathways ####
  kewg_filt <- kewg[[i]] %>% dplyr::filter(P.DE < 0.05)
  toplot=limma::topKEGG(kewg_filt, n= Inf)
  
  #Change to character and then back to factor, to keep the order from TopGO
  toplot$Pathway = as.character(toplot$Pathway)
  toplot$Pathway = factor(toplot$Pathway, levels = unique(toplot$Pathway))
  
  keggpath[[i]]=toplot
}

names(goterms) <- names(my.scWGCNA$modules)
names(keggpath) <- names(my.scWGCNA$modules)

saveRDS(goterms, paste0("~/spinal_cord_paper/output/", my.se@project.name, "_module_all_GOTerms_",format(Sys.Date(), "%d%m%y"),".rds"))
saveRDS(keggpath, paste0("~/spinal_cord_paper/output/", my.se@project.name, "_module_all_KEGGPath_",format(Sys.Date(), "%d%m%y"),".rds"))
