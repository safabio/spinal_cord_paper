#' Plots a barplot with the fraction of the expressed and orthologous genes for each module of a scWGCNA comparison.
#' 
#' This function will plot a barplot showing what's the fraction of genes in each module that is (if its the case) present in the provided 1-2-1 orthologues table, and the genes that are expressed in all test samples, and therefore used for the comparison.
#' @param scWGCNA.comp.data scWGCNA comparative data as calculated by scWGNA.compare().
#' @return A ggplot barplot
#' @export
#' @examples
#' 
#' # A pre-calculated list scWGCNA comparative data, calculated with scWGCNA.compare
#' class(MmvGg.comparative)
#' 
#' # Plot the fraction of genes used for the comparison.
#' scW.p.modulefrac(MmvGg.comparative)

# The function, taking only the pre-calculated data
scW.p.modulefrac.mod = function(scWGCNA.comp.data){
  
  # Take the precalculated number of genes lost
  to.plot = scWGCNA.comp.data$misc$modulefrac
  
  # Make them into decimal fractions
  to.plot[,-1] = to.plot[,-1] / to.plot[,2]
  
  # If we had genes lost to orthology
  if (ncol(to.plot) == 4) {
    # Make the difference, to plot them stacked
    to.plot$diff = to.plot$ortho - to.plot$expr
    # Get rid of that column, no needed
    to.plot = to.plot[,-3]
  }
  
  # melt it, to have it ggplot friendly
  to.plot = reshape::melt(to.plot[,c(1,3:ncol(to.plot))])
  
  # Plot it!
  ggplot2::ggplot(to.plot, 
                  ggplot2::aes(x=module, 
                               y=value, 
                               fill = variable,
                               color=module)) +
    ggplot2::geom_bar(position=ggplot2::position_stack(reverse = TRUE), 
                      stat = "identity",
                      size=2,
                      width = 0.7) + 
    ggplot2::scale_color_manual(values = as.character(to.plot$module)) + 
    ggplot2::scale_fill_manual(values=c("gray","gray93"),
                               labels = c("expressed in\nall samples","orthologous\nnon-expressed")) + 
    ggplot2::theme_light() + 
    ggplot2::guides(color="none") +
    ggplot2::labs(y="fraction", 
                  fill="Genes in module") + 
    ggplot2::theme(legend.text = ggplot2::element_text(margin = ggplot2::margin(t = 10, b=10))) +
    coord_flip()
  
}

#' Plots dotplot showing the z-score values for preservation and other aspects of it
#' 
#' This function will plot a dotplot, with the zscore values for global preservation, as well as density and connectivity for each module and each test sample that was compared. Can also plot the median rank.
#' @param scWGCNA.comp.data scWGCNA comparative data as calculated by scWGNA.compare().
#' @param to.plot character or character vector. Which aspects of the preservation should be plotted? Options: "preservation", "median.rank", "density", "connectivity".
#' @param test.samples character or character vector. Which test samples to plot. Default is all
#' @param pt.style vector of integers which relate to the pch for each sample. 
#' @return Either a single ggplot dotoplot showing the desired aspect of preservation. If several were requested, a gridExtra of the different ggplot dotplots.
#' @export
#' @examples
#' 
#' # S pre-calculated list scWGCNA comparative data, calculated with scWGCNA.compare
#' class(MmvGg.comparative)
#' 
#' # Plot the overall preservation and median rank.
#' scW.p.preservation(scWGCNA.comp.data = MmvGg.comparative, to.plot=c("preservation", "median.rank"))

# The function takinf the comparative data and what will be plotted.
scW.p.preservation.mod = function(
  scWGCNA.comp.data,
  to.plot=c("preservation", "median.rank"),
  test.samples = NULL,
  pt.style = NULL) {
  
  if (!is.null(pt.style)) {
    if (!all(pt.style %in% 0:25)) {
      stop("elements of pt.style must be valid pch value (integers from 0 to 25)!")
    }
    if (!length(pt.style) == length(scWGCNA.comp.data$misc$testnames)) {
      stop("pt.style must be the same lenght as the test sample!")
    }
  }
  
  # leave grey and gold modules out
  modColors = rownames(scWGCNA.comp.data$preservation$observed[[1]][[2]])
  plotMods = !(modColors %in% c("grey", "gold"))
  
  # The variable where we keep the data to plot
  plotData = data.frame(Rank = numeric(), Zsum= numeric(), Density = numeric(), Connectivity = numeric(),
                        Size = numeric(), Cols = character(), Sample = character())
  
  # Fill in with data from the samples
  for (test in 1:(length(scWGCNA.comp.data$quality$observed[[1]])-1)) {
    plotData = rbind(plotData,data.frame(Rank = scWGCNA.comp.data$preservation$observed[[1]][[test+1]][, 2],
                                         Zsum= scWGCNA.comp.data$preservation$Z[[1]][[test+1]][, 2],
                                         Density = scWGCNA.comp.data$preservation$Z[[1]][[test+1]][, 3],
                                         Connectivity = scWGCNA.comp.data$preservation$Z[[1]][[test+1]][, 4],
                                         Size = rank(scWGCNA.comp.data$preservation$observed[[1]][[test+1]][, 1], ties.method = "first"),
                                         Cols = rownames(scWGCNA.comp.data$preservation$observed[[1]][[test+1]]),
                                         Sample = rep(scWGCNA.comp.data$misc$testnames[test],nrow(scWGCNA.comp.data$preservation$observed[[1]][[test+1]]))))
  }
  
  # We kick out the golden and grey modules
  plotData = plotData[rep(plotMods, length(scWGCNA.comp.data$quality$observed[[1]])-1),]
  
  # For the density and connectivity, we will set the same limits
  my.lim = (max(plotData[,3:4])-min(plotData[,3:4]))*0.025
  my.lim = c(min(plotData[,3:4])-my.lim, max(plotData[,3:4])+my.lim)
  
  #For the median rank, we will show transparency for modules without evidence of preservation
  my.alpha = rep(1,dim(plotData)[1])
  my.alpha[which(plotData$Zsum < 2)] = 0.75
  plotData$alpha = my.alpha
  
  # The plotting of the Zsummary and Median Rank, against the module size
  my.p=list()
  
  # The labels for the plots
  my.plotData = plotData
  my.plotData$Sample = factor(my.plotData$Sample, levels = scWGCNA.comp.data$misc$testnames)
  my.plotData$mylabels = my.plotData$Cols
  my.plotData$mylabels[which(plotData$Zsum < 2)] = NA
  
  plotData = my.plotData #?????????
  
  # Do we plot only certain test samples?
  if (!is.null(test.samples)) {
    plotData = plotData[which(plotData$Sample %in% test.samples),]
  }
  
  # A global theme, for all plots
  my.theme = list(ggplot2::scale_color_manual(values= as.character(unique(plotData$Cols))),
                  ggplot2::scale_fill_manual(values= rep("grey",length(unique(plotData$Cols)))),
                  ggplot2::theme_classic())
  
  if (!is.null(pt.style)) {
    my.theme = c(my.theme, ggplot2::scale_shape_manual(name=  "Test\nsample:", values= as.integer(pt.style),
                                                labels=scWGCNA.comp.data$misc$testnames))
  } else{
    my.theme = c(my.theme, ggplot2::scale_shape_manual(name= "Test\nsample:",values= c(0:25),
                                                       labels=scWGCNA.comp.data$misc$testnames))
  }
  
  # Some other shared theme items
  my.gp = list(ggplot2::geom_point(ggplot2::aes(color=Cols, fill=Cols, stroke=1.2, shape=Sample), size=3.5),
               ggplot2::geom_hline(yintercept = c(2,10), linetype="dashed", color=c("red", "limegreen")))
  my.gpa = ggplot2::geom_point(ggplot2::aes(color=Cols, fill=Cols, stroke=1.2, shape=Sample, alpha=alpha), size=3.5)
  my.box = annotate(
    "rect",
    fill = "white",
    xmin = 0,
    xmax = max(plotData$Size) + 1,
    ymin = -13,
    ymax = -3
    )
  my.label = annotate(
    'text',
    x=plotData$Size,
    y=-8,
    label = plotData$Cols,
    angle = 90
    )
  
  
  # Check which plots to make
  if (any(to.plot %in% "preservation")) {
    
    my.p[["preservation"]] = ggplot2::ggplot(plotData, ggplot2::aes(x=Size, y=Zsum, label=Cols)) +
      geom_vline(xintercept = 1:max(plotData$Size),
                 cex = 0.2,
                 col = "grey") +
      my.gp +
      my.theme +
      ggplot2::labs(x ="Module size rank", y = "Zsummary") +
      my.box +
      my.label +
      ylim(-13,70)
  }
  
  if (any(to.plot %in% "median.rank")) {
    
    my.p[["median.rank"]] = ggplot2::ggplot(plotData, ggplot2::aes(x=Size, y=Rank)) +
      ggplot2::scale_y_continuous(trans = "reverse", expand = c(0.1,0), limits = c(max(plotData$Rank) + 8, 0)) +
      geom_vline(xintercept = 1:max(plotData$Size),
                 cex = 0.2,
                 col = "grey") + 
      my.gpa +
      my.theme +
      ggplot2::labs(x ="Module size rank", y = "Median Rank") +
      annotate(
        "rect",
        fill = "white",
        xmin = 0,
        xmax = max(plotData$Size) + 1,
        ymin = max(plotData$Rank) + 1,
        ymax = max(plotData$Rank) + 8
      ) +
      annotate(
        'text',
        x=plotData$Size,
        y=max(plotData$Rank) + 4.5,
        label = plotData$Cols,
        angle = 90
      )
  }
  
  if (any(to.plot %in% "density")) {
    
    my.p[["density"]] = ggplot2::ggplot(plotData, ggplot2::aes(x=Size, y=Density, label=Cols)) +
      geom_vline(xintercept = 1:max(plotData$Size),
                 cex = 0.2,
                 col = "grey") +
      my.gp +
      my.theme +
      ggplot2::labs(x ="Module size rank", y = "Density") +
      my.box +
      my.label +
      ylim(-13,70)
    
  }
  
  if (any(to.plot %in% "connectivity")) {
    
    my.p[["connectivity"]] = ggplot2::ggplot(plotData, ggplot2::aes(x=Size, y=Connectivity, label=Cols)) +
      geom_vline(xintercept = 1:max(plotData$Size),
                 cex = 0.2,
                 col = "grey") +
      my.gp +
      my.theme +
      ggplot2::labs(x ="Module size rank", y = "Connectivity") +
      my.box +
      my.label +
      ylim(-13,70)
  }
  
  # Arrange the plots as needed
  my.p = my.p[to.plot]
  
  # If put the legend in the last plot
  my.p[[length(my.p)]] = my.p[[length(my.p)]] + 
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8), legend.title = ggplot2::element_text(size=10)) +
    ggplot2::guides(fill = "none", alpha="none", color="none")
  
  # If we have more than one, remove the legend from the non-last plots
  if (length(my.p) > 1) {
    my.p[-length(my.p)] = lapply(my.p[-length(my.p)], function(gp){
      gp + ggplot2::theme(legend.position="none")
    })
    gridExtra::grid.arrange(grobs=my.p, ncol=2)
  } else{my.p}
  
}