library(readxl)
library(data.table)
library(miRBaseConverter)
library(biomaRt)
library(lattice)

source('src/func.R')



lastal.res <- fread("grep -v '^#' out/lnc_mir_prediction/myalns.tab")

lncrnas.precursors <- fread('data/lncrna_mir_precursors.txt', header = F)
setnames(lncrnas.precursors, 
         colnames(lncrnas.precursors), 
         c('lncrna', 'mirna.precursor'))

fantom_DE <- fread('../oligo_DE_Summary_gene_filtered.tsv')
fantom_DE$geneSymbol <- stringr::str_remove(fantom_DE$geneSymbol, 'HG')


version=checkMiRNAVersion(lncrnas.precursors$mirna.precursor, verbose = TRUE)
mirnatable <- as.data.table(getMiRNATable(version = "v22", species = "hsa"))
mir.mature.precursor <- 
  mirnatable[Precursor %in% lncrnas.precursors$mirna.precursor, c('Precursor', 'Mature1', 'Mature2')]
mature1 <- mir.mature.precursor[, c('Precursor', 'Mature1')]
setnames(mature1, 'Mature1', 'mature')
mature2 <- mir.mature.precursor[, c('Precursor', 'Mature2')]
setnames(mature2, 'Mature2', 'mature')
mir.mature.precursor <- rbindlist(list(mature1, mature2))
mir.mature.precursor <- mir.mature.precursor[!is.na(mature)]


lncrnas.precursors.mature <- merge(x = lncrnas.precursors,
                                   y= mir.mature.precursor,
                                   by.x = 'mirna.precursor',
                                   by.y = 'Precursor',
                                   allow.cartesian = T)
lastal.res.filtered <- lastal.res[!(V7 %in% lncrnas.precursors$lncrna & V2 %in% lncrnas.precursors.mature$mature)]



ensembl = useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

attributes = listAttributes(ensembl)
filters = listFilters(ensembl)

lncRNAs <- getBM(attributes = c('ensembl_gene_id', "external_transcript_name"), 
                 filters = 'ensembl_gene_id',
                 values = unique(fantom_DE$KD.geneID), 
                 mart = ensembl)
lncRNAs <- as.data.table(lncRNAs)
unique(lncRNAs[external_transcript_name %in% lastal.res.filtered[,.N,by=V7]$V7]$ensembl_gene_id)

lastal.res.cut <- unique(lastal.res[, c('V2', 'V7')])
setnames(lastal.res.cut, colnames(lastal.res.cut), c('mature.mirna', 'lncrna.transcript'))
lastal.res.cut <- merge(x = lastal.res.cut,
                        y = lncRNAs,
                        by.x = 'lncrna.transcript',
                        by.y = 'external_transcript_name')
lastal.res.cut <- unique(lastal.res.cut[,c('mature.mirna', 'ensembl_gene_id')])
setnames(lastal.res.cut, 'ensembl_gene_id', 'lncrna')

# add info about targets of the mirnas from literature
mirtarbase_valid <- as.data.table(read_excel('data/hsa_MTI.xlsx'))[`Species (miRNA)`=='Homo sapiens']
lnc.mir.tar <- unique(merge(x = lastal.res.cut,
                           y = mirtarbase_valid[, c("miRNA", "Target Gene")],
                           by.x = 'mature.mirna',
                           by.y = "miRNA",
                           allow.cartesian = T))

fantom_DE <- fread('../oligo_DE_Summary_gene_filtered.tsv')
fantom_DE$geneSymbol <- stringr::str_remove(fantom_DE$geneSymbol, 'HG')

lnc.mir.tar.fantom <- merge(x = lnc.mir.tar,
                            y = fantom_DE,
                            by.x = c('lncrna', 'Target Gene'),
                            by.y = c('KD.geneID', 'geneSymbol'))
lnc.mir.tar.fantom <- 
  unique(lnc.mir.tar.fantom[,c('lncrna', 'mature.mirna', 'Target Gene', 
                               'log2FC', 'KD.geneSymbol')])
# write.table(lnc.mir.tar.fantom, 'out/lnc_mir_prediction/lnc.mir.tar.fantom.mir.included.tsv')

num_mirnas_per_target <- 
  lnc.mir.tar.fantom[, length(unique(.SD$mature.mirna)), by = 'Target Gene']
histogram(num_mirnas_per_target$V1, 
          type = 'count',
          nint = 100, 
          main = 'Distribution of miRNAs per target',
          xlab = 'Number of miRNAs that sponge the target',
          col = "#AEAEEE96")

la <- lnc.mir.tar.fantom[, length(unique(.SD$mature.mirna)), by = "KD.geneSymbol"]
setnames(la, 'V1', 'num_mirnas_per_lncrna')
setorder(la, num_mirnas_per_lncrna)
write.table(la, 'out/lnc_mir_prediction/mirnas_per_lncrna.tsv', 
			quote = F, row.names = F, sep = '\t')

histogram(la$V1, 
          type = 'count',
          nint = 100, 
          main = 'Distribution of predicted miRNA targets per lncRNA',
          xlab = 'Number of predicted miRNA targets per lncRNA',
          col = "#EEAEEE96")

la$KD.geneSymbol <- factor(la$KD.geneSymbol, 
                           levels = unique(la$KD.geneSymbol))
dev.off()


is.even <- function(x) x %% 2 == 0


# dotplot num of predicted mirnas per lncrna
png('out/lnc_mir_prediction/mirnas_per_lncrna.png', width = 2100)
mypanel<-function(x,y,...){
  panel.dotplot(x, y, ...)
  panel.text(x = c(1:137)[1:26],
             y = as.numeric(la$num_mirnas_per_lncrna)[1:26]+70,
             labels=la$num_mirnas_per_lncrna[1:26],
             cex = 0.85)
  panel.text(x = c(1:137)[27:137][is.even(c(27:137))],
             y = as.numeric(la$num_mirnas_per_lncrna)[27:137][is.even(c(27:137))]-60,
             labels=la$num_mirnas_per_lncrna[27:137][is.even(c(27:137))],
             cex = 0.8)
  panel.text(x = c(1:137)[27:137][!is.even(c(27:137))]-0.1,
             y = as.numeric(la$num_mirnas_per_lncrna)[27:137][!is.even(c(27:137))]+60,
             labels=la$num_mirnas_per_lncrna[27:137][!is.even(c(27:137))],
             cex = 0.8)
}

dotplot(num_mirnas_per_lncrna ~ KD.geneSymbol, la,
        main = list(label = 'Distribution of predicted miRNA targets per lncRNA',
                    cex = 1.2),
        ylab = list(label = 'Number of predicted miRNAs per lncRNA',
                    cex = 1),
        scales=list(y=list(rot=0, tick.number = 20, cex = 1), 
                    x=list(rot=45, cex = 1),
                    relation = 'free'),
        ylim = c(110,2500),
        panel=mypanel)
dev.off()

lnc.tar.fantom <- unique(lnc.mir.tar.fantom[,c('lncrna', 'Target Gene', 'log2FC', 'KD.geneSymbol')])
boxplot(lnc.tar.fantom$log2FC)
# write.table(lnc.tar.fantom, 'out/lnc_mir_prediction/lnc.tar.fantom.tsv', 
# 	quote = F, row.names = F, sep = '\t')

fantom_DE_not_sponge <- fantom_DE[!(KD.geneID %in% lnc.tar.fantom$lncrna &
                                      geneSymbol %in% lnc.tar.fantom$`Target Gene`)]

# scatterplot num of predicted mirnas vs 
# num of targets of the mirna in fatnom per lncrna 

merge(x = la,
      y = lnc.tar.fantom[, length(unique(.SD$`Target Gene`)), 
                             by = "KD.geneSymbol"],
      by = "KD.geneSymbol")


# beautiful pictures
lnc.tar.fantom <- fread('out/lnc_mir_prediction/lnc.tar.fantom.tsv')
# compare logFC of targets of mirna sponged by lncrna & all targets of lncrna
data <- data.table(ceRNA = T,
                   log2FC =lnc.tar.fantom$log2FC)
data <- rbindlist(list(data,
                       data.table(ceRNA = F,
                                  log2FC = fantom_DE$log2FC)))
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'all', 'a sponge')
wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'all']$log2FC, paired = F)

#ylim1 = boxplot.stats(data$log2FC)$stats[c(1, 5)]

png('out/mirtarbase+lastal/boxplot.mir.targets.vs.all.lnc.targets.png', width = 700)
DrawTwoBoxplots(data,
                main = "All targets of miRNAs sponged by lncRNAs vs all targets of lncRNAs\np-value < 2.2e-16")
dev.off()

# compare only DE targets of mirna sponged by lncrna & DE targets of non-sponge lncrna
data <- data.table(ceRNA = T,
                   log2FC = lnc.tar.fantom[log2FC > 1 | log2FC < -1]$log2FC)
data <- rbindlist(list(data,
                       data.table(ceRNA = F,
                                  log2FC = fantom_DE[log2FC > 1 | log2FC < -1]$log2FC)))
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'all', 'a sponge')
wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'all']$log2FC, paired = F)

png('out/mirtarbase+lastal/boxplot.DE.mir.targets.vs.DE.all.lnc.targets.png', width = 700)
DrawTwoBoxplots(data,
                main = "DE targets of miRNAs sponged by lncRNAs vs DE targets of non-sponge lncRNAs\np-value < 2.2e-16")
dev.off()

# all targets of fantom lncrnas except targets of sponged mirnas
data <- data.table(ceRNA = T,
                   log2FC =lnc.tar.fantom$log2FC)
data <- rbindlist(list(data,
                       data.table(ceRNA = F,
                                  log2FC = fantom_DE_not_sponge$log2FC)))
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'not a sponge', 'a sponge')
wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'not a sponge']$log2FC, paired = F)

#ylim1 = boxplot.stats(data$log2FC)$stats[c(1, 5)]

png('out/mirtarbase+lastal/boxplot.mir.targets.vs.all.targets.non-sponge.lnc.png', width = 700)
DrawTwoBoxplots(data,
                main = "All targets of miRNAs sponged by lncRNAs vs all targets of lncRNAs\np-value < 2.2e-16")
dev.off()


# all three in one
list.val = list(
  all = fantom_DE$log2FC,
  `not targeted by\nsponged miRNAs` = fantom_DE_not_sponge$log2FC,
  `targeted by\nsponged miRNAs` = lnc.tar.fantom$log2FC,
  `#(mir_per_lnc) > 20` = 
    lnc.tar.fantom[KD.geneSymbol %in% la[num_mirnas_per_lncrna > 20]$KD.geneSymbol]$log2FC)
png('out/mirtarbase+lastal/three.boxplots.signif.png', width = 700, height = 600)
DrawThreeBoxplots(list.val, ylab = 'logFC', main = '', ylim = c(min(unlist(list.val)), 11),
                  cex.lab.x = 1.1)
y <- max(list.val$all)+1
# set an offset for tick lengths
offset <- 0.2
# draw first horizontal line
lines(c(1.3, 1.6),c(y, y))
# draw ticks
lines(c(1.3, 1.3),c(y, y-offset))
lines(c(1.6, 1.6),c(y, y-offset))
text(1.45,y+0.5,"***")

y <- max(list.val$all)+3
# set an offset for tick lengths
offset <- 0.2
# draw first horizontal line
lines(c(1,1.3, 1.6),c(y, y, y))
# draw ticks
lines(c(1,1),c(y, y-offset))
lines(c(1.6,1.6),c(y, y-offset))
text(1.3,y+0.5,"***")
dev.off()

# chi squared
contingency.table.sponge <- data.table('logFC < 0' = sum(lnc.tar.fantom$log2FC < 0),
                                       'logFC > 0' = sum(lnc.tar.fantom$log2FC>0))
contingency.table.not.sponge <- data.table('logFC < 0' = sum(fantom_DE_not_sponge$log2FC < -0),
                                           'logFC > 0' = sum(fantom_DE_not_sponge$log2FC>0))
contingency.table <- rbindlist(list(contingency.table.not.sponge,
                                    contingency.table.sponge))
contingency.table
chisq.test(contingency.table)
df <- contingency.table
rownames(df) <- c('not a sponge', 'a sponge')
library("graphics")
library(viridis)

#png('out/chi.squared.sponges.png', width = 700)
mosaicplot(df,
           main = "Gene expression falls more often in case of non-sponge lncRNAs",
           sub = 'chi squared test p-value < 2.2e-16',
           cex.axis = 1,
           color = viridis(4, alpha = 0.5))
mtext(df[2, 'logFC < 0'], at = 0.75, line = -2)
mtext('lalala', at = 0.75, line = -5.5)
mtext('lalala', at = 0.25, line = -6.5)
mtext('lalala', at = 0.25, line = -3)
# dev.off()

library("vcd")
data <- data.table(`is the lncRNA a sponge?` = T,
                   log2FC =lnc.tar.fantom$log2FC)
data <- rbindlist(list(data,
                       data.table(`is the lncRNA a sponge?` = F,
                                  log2FC = fantom_DE_not_sponge$log2FC)))
data$`is the lncRNA a sponge?` <- as.factor(data$`is the lncRNA a sponge?`)
levels(data$`is the lncRNA a sponge?`) <- ifelse(levels(data$`is the lncRNA a sponge?`) == F, 'not a sponge', 'a sponge')
data$`sign of logFC` <- ifelse(data$log2FC>0, 'logFC > 0', 'logFC < 0')
data$`sign of logFC` <- as.factor(data$`sign of logFC`)

struct <- structable(`sign of logFC` ~ `is the lncRNA a sponge?`, data = data)
png('out/mirtarbase+lastal/mosaic.png', 700)
mosaic(struct, shade = TRUE, direction = "v", pop = FALSE)
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
dev.off()
