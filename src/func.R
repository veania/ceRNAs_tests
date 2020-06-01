library(tidyverse)
library(hrbrthemes)
library(viridis)

DrawTwoBoxplots <- function(data, main){
  ggplot(data, aes(x=ceRNA, y=log2FC, fill=ceRNA)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme_ipsum(axis_title_size = 15,axis_text_size = 15) +
    theme(
      legend.position="none",
      plot.title = element_text(size=13)
    ) +
    ggtitle(main) +
    xlab("")
}

DrawBaseBoxplot <- function(list.val, cex.lab, ylim, proportion, at.x, xlim, col){
  boxplot(list.val,
          pch = 20, 
          col = col,
          axes = FALSE, # Don't plot the axes
          frame.plot = FALSE,
          range = 0,
          ylim = ylim,
          width = proportion,
          at = at.x,
          boxwex = 0.4,
          xlim = xlim)
}
DrawThreeBoxplots <- function(list.val, ylab, main, ylim = c(-10,10), cex.lab.y = 1.2, cex.lab.x = 1){
  at.y = pretty(unlist(list.val), 7)
  proportion <- unlist(lapply(list.val, length))/sum(unlist(lapply(list.val, length)))
  at.x <- seq(1,by=.3,length.out=3)
  xlim = c(.8,1.8)
  DrawBaseBoxplot(list.val, cex.lab, ylim, proportion, at.x, xlim, col = F)
  abline(h = at.y, col = 'grey80')
  # grid(5,5)
  par(new=TRUE)
  DrawBaseBoxplot(list.val, cex.lab, ylim, proportion, at.x, xlim, col = viridis(3, alpha=0.6))
  for (i in 1:length(list.val)){
    mtext(side = 1, text = names(list.val)[i], at = at.x[i],
          col = "grey20", line = 1, cex = cex.lab.x)
  }
  mtext(side = 2, text = ylab, at = 4,
        col = "grey20", line = 2, cex = cex.lab.y)
  #at = at[-c(1,length(at))]
  mtext(side = 2, text = at.y, at = at.y, col = "grey20", line = 1, cex = 0.9)
  mtext(side = 3, text = main,
        col = "grey20", line = 1, cex = 1.5)
}
