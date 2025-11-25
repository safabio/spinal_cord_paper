################

# Code: Patrick Tschopp
# Input produced in: ~/spinal_cord_paper/markdown/module_genes_TSS_tables.Rmd

#################
################
# TSS to BED [R]

TSS=readRDS("Gg_lumb_int_GB_module_TSS_table.rds")
colnames(TSS)
# [1] "seqid"        "type"         "start"        "end"          "strand"       "strand_start" "gene_id"     
# [8] "gene_name"    "module"       "Membership"   "p.val"
TSS_bed=cbind(TSS[,1],(as.numeric(TSS[,6])-500),(as.numeric(TSS[,6])+500),TSS[,8],TSS[,9])

# write BED file
write.table(TSS_bed,"TSS-Prom_BrownMagentaLumbar.bed",sep = "\t",col.names = F,row.names = F,quote = F)


################
# BED to FASTA [Unix shell]

module purge
module load BEDTools/2.30.0-GCC-10.3.0

dir=$(pwd)
for file in $(ls $dir/*.bed); do 
name=${file##*/}
  base=${name%.bed}
  bedtools getfasta -fo "${base}.fasta" -fi ../references/genomes/ENS_g6/10XGenome_ENSG6_extended_060420/fasta/genome.fa -bed ${file}
  done
  
  
  ################
  # SEA analysis [Unix shell]
  
  module purge
  module load MEME/5.5.7-gompi-2023b
  
  dir=$(pwd)
  for file in $(ls $dir/*.fasta); do 
  name=${file##*/}
    base=${name%.fasta}
    sea --verbosity 1 --o ../SEA/GB/${base} --thresh 100.0 --align center --p ${file} --m ../SEA/LNS_denovoMotifs_PWM_addLogOdds.pssm.meme --m ../SEA/JASPAR2024_vertebrate_non-red.meme
    done
    
    
    ################
    # motif-enrich vs TF-module-corr [R]
    
    setwd("")
    library(stringr)
    library(dplyr)
    library("tidyverse")
    library(ggplot2)
    library(ggrepel)
    
    # read in corresponding TF Pearson correlation file
    Pcorr=read.csv("../module_GB_TF_corr_pearson.csv", sep = ",",header = T)
    
    # read in and analyze SEA motif enrichment files
    TFmod=data.frame(matrix(ncol = 7, nrow = 0))
    colnames(TFmod)=c("gene.name","module.corr","pearson.corr","motif.name","motif.ID","motif.enrich","module")
    
    folders=list.files()
    mods=sapply(strsplit(folders,"_"), `[`, 2)
    
    for(i in 1:length(folders)){
      
      sea=read.table(paste(getwd(),folders[i],"sea.tsv",sep="/"), sep = "\t",header = T)
      SingleMotifs=toupper(unlist(strsplit(sea$ALT_ID,"::")))
      nmotifs=length(SingleMotifs)
      
      commonPcorr=Pcorr[Pcorr$Gene.name %in% intersect(Pcorr$Gene.name,SingleMotifs),]
      commonPcorr=commonPcorr[,c(5,6,4)]
      commonPcorrMod=commonPcorr %>%
        filter(str_detect(tf_member, mods[i]))
      rownames(commonPcorrMod)=commonPcorrMod[,1]
      commonPcorrMod=commonPcorrMod[ order(row.names(commonPcorrMod)), ]
      nTFcorr=nrow(commonPcorrMod)
      
      sea$single_ID=toupper(sapply(strsplit(sea$ALT_ID,"::"), `[`, 1))
      commonseaMod=sea[toupper(sea$single_ID) %in% intersect(Pcorr$Gene.name,SingleMotifs),]
      commonseaMod=commonseaMod[,c(4,18,3,10)]
      rownames(commonseaMod)=make.names(commonseaMod$single_ID, unique = TRUE)
      commonseaMod=commonseaMod[order(row.names(commonseaMod)),]
      commonseaMod$Gene.name=toupper(commonseaMod$single_ID)
      
      motifsTFs=dplyr::left_join(commonseaMod,commonPcorrMod,  by  ="Gene.name")
      rownames(motifsTFs)=rownames(commonseaMod)
      nmotifsTFs=nrow(motifsTFs)
      
      sp=ggplot(motifsTFs, aes(corr, ENR_RATIO, label = rownames(motifsTFs)))+
        geom_point(colour = mods[i])+
        labs(x="TF-module Pearson corr",y="motif enrich ratio")+
        ggtitle(paste0("module: ",mods[i], " [motifs enriched: ", nmotifs, ", TFs expr.: ", nTFcorr, ", motif/TF combos: ", nmotifsTFs, "]"))
      pdf_filename=paste0(mods[i], "_motifEnrich_TFcorr.pdf")
      pdf(pdf_filename, width = 10, height = 12)
      print(sp + ggrepel::geom_text_repel(color = "black", max.overlaps=50, size = 3))
      dev.off()
      
      temp=data.frame(motifsTFs$Gene.name,motifsTFs$tf_member,motifsTFs$corr,motifsTFs$ALT_ID,motifsTFs$ID,motifsTFs$ENR_RATIO)
      temp$module=strsplit(motifsTFs$tf_member,"_")[[1]][1]
      colnames(temp)=c("gene.name","module.corr","pearson.corr","motif.name","motif.ID","motif.enrich","module")
      TFmod=rbind(TFmod,temp)
      
    }
    
    saveRDS(TFmod,file="GB_TFmoduleCorr_motifEnrich.rds")
    
#########    

# Output is used in ~/spinal_cord_paper/scripts/Fig_2_plots.R and ~/spinal_cord_paper/scripts/Fig_4_plots.R
        