### ### ### ### ### ### ### ### ### ### ###
#### Additional plots and code for figure 2 ####
## Fabio Sacher, 21.05.2024
### ### ### ### ### ### ### ### ### ### ###

### ### ### ### ### ### ### ### ### ### ###
#### scWGCNA module network plots ####
### ### ### ### ### ### ### ### ### ### ###

# code is taken and adapted from scWGCNA::scW.p.network

# scWGNA.data
scWGCNA.data = readRDS("~/spinal_cord_paper/output/Gg_devel_int_scWGCNA_250723.rds")
gnames = modplots::gnames
rownames(gnames) <- gnames[,1]
modules = c("darkred", "darkgreen", "lightgreen")


# get the module colors
my.cols = levels(as.factor(scWGCNA.data[["dynamicCols"]]))

nwrk <- list()

for (i in seq(modules)) {
  col <- modules[i]
  # get index if color is provided
  module = which(my.cols == col)
  # extract the module
  mynet = scWGCNA.data[["networks"]][[module]]
  
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
  
  ggsave(paste0("~/spinal_cord_paper/figures/Fig_2_", col,"network.pdf"), width = 5, height = 5)
}


