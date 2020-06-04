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
DrawMultipleBoxplots <- function(list.val, ylab = '', main = '', 
                              xlim = c(.8,1.8), ylim = c(-10,10), cex.lab.y = 1.2, cex.lab.x = 1){
  at.y = pretty(unlist(list.val), 7)
  proportion <- unlist(lapply(list.val, length))/sum(unlist(lapply(list.val, length)))
  at.x <- seq(1,by=.3,length.out=length(list.val))
  DrawBaseBoxplot(list.val, cex.lab, ylim, proportion, at.x, xlim, col = F)
  abline(h = at.y, col = 'grey80')
  # grid(5,5)
  par(new=TRUE)
  DrawBaseBoxplot(list.val, cex.lab, ylim, proportion, at.x, xlim, col = viridis(length(list.val), alpha=0.6))
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

is.even <- function(x) x %% 2 == 0

PrettyScatter <- function(x, y, bg = "#EEAEEE96", 
                          main = '', panel.first.step = 1, cex.lab = 1.5, 
                          xlab = '', ylab = ''){
  plot(x, y, xlim = c(min(x), max(x)), ylim = c(min(y), max(y)),
       xlab = xlab, ylab = ylab, 
       cex.lab = cex.lab,
       bg = bg, # Fill colour
       col = "grey10",
       pch = 21,
       axes = FALSE, # Don't plot the axes
       frame.plot = FALSE, # Remove the frame 
       #panel.first = ,
       main = main)
  abline(h = seq(round(min(x)), round(max(y)), panel.first.step), col = 'grey60')
  
  at = pretty(x)
  at = at[-c(1,length(at))]
  mtext(side = 1, text = at, at = at, 
        col = "grey20", line = 1, cex = 0.9)
  at = pretty(y)
  at = at[-c(1,length(at))]
  mtext(side = 2, text = at, at = at, col = "grey20", line = 1, cex = 0.9, las = 2)
}