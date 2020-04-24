library(data.table)
home.dir <- '/home/anna/medved/home/aelizarova/'

DE <- fread(file.path(home.dir, 'ceRNAs/oligo_DE_Summary_gene_filtered.tsv'))
lncRNAs_summary <- fread(file.path(home.dir, 'ceRNAs/lncRNAs_summary.tsv'))

lncRNAs_summary$lncRNA[lncRNAs_summary$lncRNA == "RP11-398K22"] <- "RP11-398K22.12"
lncRNAs_summary$lncRNA[lncRNAs_summary$lncRNA == "AC091729"] <- "AC091729.9"
lncRNAs_summary$lncRNA[lncRNAs_summary$lncRNA == "RP11-38L15"] <- "RP11-38L15.3"
lncRNAs_summary$lncRNA[lncRNAs_summary$lncRNA == "NEAT2(MALAT1)"] <- "MALAT1"
lncRNAs_summary$lncRNA[lncRNAs_summary$lncRNA == "ANRIL"] <- "CDKN2B-AS1"

ceRNA.logFC.dt <- merge(lncRNAs_summary[, c('lncRNA', 'ceRNA')],
                        DE[, c('geneSymbol', 'log2FC', 'KD.geneSymbol')],
                        by.x = 'lncRNA',
                        by.y= 'KD.geneSymbol')
ceRNA.logFC.dt[is.na(ceRNA)]$ceRNA <- 0
ceRNA.logFC.dt.cut <- ceRNA.logFC.dt[log2FC>1 | log2FC< -1]

contingency.table <- as.data.table(table(ceRNA.logFC.dt.cut[ceRNA == 1]$log2FC>1))
contingency.table[, expr.change:=ifelse(V1 == T, 'rise', 'fall')]
contingency.table$V1 <- NULL
setnames(contingency.table, 'N', 'sponge')
contingency.table[, not.sponge:= c(as.data.table(table(ceRNA.logFC.dt.cut[ceRNA == 0]$log2FC>1))[V1==F]$N,
                                   as.data.table(table(ceRNA.logFC.dt.cut[ceRNA == 0]$log2FC>1))[V1==T]$N)]
setcolorder(contingency.table, c('expr.change', 'sponge', 'not.sponge'))
# contingency.table <- as.data.table(t(contingency.table))
# contingency.table <- rbindlist(list(contingency.table,
#                                     as.data.table(t(colnames(contingency.table)))))
# contingency.table <- as.data.table(t(contingency.table))
# colnames(contingency.table) <- as.character(contingency.table[1,])
# contingency.table <- contingency.table[2:3,]
# setnames(contingency.table, 'expr.change', 'sponge.status')
# contingency.table[,-c('sponge.status')]$fall <- as.numeric(contingency.table[,-c('sponge.status')]$fall)
# 

test.res <- chisq.test(contingency.table[,-c('expr.change')])
df <- as.data.frame(contingency.table[,-c('expr.change')])
colnames(df)
rownames(df) <- c('log2FC < -1', 'logFC > 1')
library("graphics")

png('out/chi.squared.sponges.png', width = 700)
mosaicplot(df,
           main = "Gene expression falls more often in case of non-sponge lncRNAs", 
           sub = paste0('chi squared test p-value = ', round(test.res$p.value,33)),
           cex.axis = 1, 
           color = viridis(4, alpha = 0.5))
dev.off()

library(reshape2 )
library(ggmosaic)
data.mosaic <- melt(contingency.table)
ggplot(data = data.mosaic) +
  geom_mosaic(aes(x = product(expr.change, variable), fill=expr.change, weight = value)) +
  labs(tag = c('lala', 'lalalala'), x="Is it rude recline? ", title="Gene expression falls more often in case of non-sponge lncRNAs") +
  scale_fill_viridis(begin = 0.9, end = 0.6, discrete = TRUE, alpha=0.6) 


boxplot(ceRNA.logFC.dt[ceRNA==1]$log2FC, ceRNA.logFC.dt[ceRNA==0]$log2FC, 
        names = c('lncRNA is sponge', 'lncRNA is not a sponge'), col = viridis(2))


library(tidyverse)
library(hrbrthemes)
library(viridis)

data <- ceRNA.logFC.dt[, c('ceRNA', 'log2FC')]
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == '0', 'not a sponge', 'a sponge')


ylim1 = boxplot.stats(data$log2FC)$stats[c(1, 5)]

png('out/boxplot.sponges.png', width = 700)
# Plot
data %>%
  ggplot( aes(x=ceRNA, y=log2FC, fill=ceRNA)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("p-value < 2.2e-16") +
  xlab("")
dev.off()

wilcox.test(ceRNA.logFC.dt[ceRNA == 1]$log2FC, ceRNA.logFC.dt[ceRNA == 0]$log2FC, paired = F)
t.test(ceRNA.logFC.dt[ceRNA == 1]$log2FC, ceRNA.logFC.dt[ceRNA == 0]$log2FC)


# what if we looked at only significantly DE genes?
data <- ceRNA.logFC.dt.cut[, c('ceRNA', 'log2FC')]
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == '0', 'not a sponge', 'a sponge')


ylim1 = boxplot.stats(data$log2FC)$stats[c(1, 5)]

png('out/boxplot.logFC.cut.sponges.png', width = 700)
# Plot
data %>%
  ggplot( aes(x=ceRNA, y=log2FC, fill=ceRNA)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("p-value < 2.2e-16") +
  xlab("")
dev.off()

wilcox.test(ceRNA.logFC.dt.cut[ceRNA == 1]$log2FC, ceRNA.logFC.dt.cut[ceRNA == 0]$log2FC, paired = F)
t.test(ceRNA.logFC.dt.cut[ceRNA == 1]$log2FC, ceRNA.logFC.dt.cut[ceRNA == 0]$log2FC)