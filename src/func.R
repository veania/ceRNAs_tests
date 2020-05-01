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