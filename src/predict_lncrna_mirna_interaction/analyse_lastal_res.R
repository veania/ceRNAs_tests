library(readxl)
library(data.table)
library(miRBaseConverter)
library(biomaRt)
library(lattice)

source('src/func.R')



lastal.res <- fread("grep -v '^#' out/lnc_mir_prediction/myalns.tab")

lncrnas.precursors <- fread('data/lncrna_mir_precursors_EXACT_ASO_TRANSCRIPTS_142.txt', 
                            sep = '\t', header = F)
setnames(lncrnas.precursors, 
         colnames(lncrnas.precursors), 
         c('lncrna', 'mirna.precursor'))

fantom_DE <- fread('../oligo_DE_Summary_gene_filtered.tsv')
fantom_DE$geneSymbol <- stringr::str_remove(fantom_DE$geneSymbol, 'HG')
fantom_DE$perturb_id <- stringr::str_remove(fantom_DE$perturb_id, 'ASO_')
fantom_DE$perturb_id_prefix <- unlist(str_split_fixed(fantom_DE$perturb_id, '\\_', 2)[,1])
fantom_DE$KD.geneID_ASO_fa_format = paste(fantom_DE$KD.geneID, fantom_DE$perturb_id_prefix, sep = '|')

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
#datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# attributes = listAttributes(ensembl)
# filters = listFilters(ensembl)

lncRNAs <- getBM(attributes = c('ensembl_gene_id', "external_transcript_name"), 
                 filters = 'ensembl_gene_id',
                 values = unique(fantom_DE$KD.geneID), 
                 mart = ensembl)
lncRNAs <- as.data.table(lncRNAs)
# unique(lncRNAs[external_transcript_name %in% lastal.res.filtered[,.N,by=V7]$V7]$ensembl_gene_id)

lastal.res.cut <- unique(lastal.res[, c('V2', 'V7')])
setnames(lastal.res.cut, colnames(lastal.res.cut), c('mature.mirna', 'lncrna.transcript'))
# lastal.res.cut <- merge(x = lastal.res.cut,
#                         y = lncRNAs,
#                         by.x = 'lncrna.transcript',
#                         by.y = 'external_transcript_name')
# lastal.res.cut <- unique(lastal.res.cut[,c('mature.mirna', 'ensembl_gene_id')])
# setnames(lastal.res.cut, 'ensembl_gene_id', 'lncrna')
lastal.res.cut[, `KD.gene.ID|ASO_id`:=paste(
  str_split_fixed(lastal.res.cut$lncrna.transcript,'\\|', 4)[,c(2)],
  str_split_fixed(lastal.res.cut$lncrna.transcript,'\\|', 4)[,c(3)], 
  sep='|')]
lastal.res.cut$lncrna <- fantom_DE[match(lastal.res.cut$`KD.gene.ID|ASO_id`,
                                         KD.geneID_ASO_fa_format)]$KD.geneSymbol
lastal.res.cut[, ASO_id:=
  str_split_fixed(lastal.res.cut$lncrna.transcript,'\\|', 4)[,c(3)]]
lastal.res.cut$lncrna <- fantom_DE[match(lastal.res.cut$ASO_id, perturb_id_prefix)]$KD.geneSymbol
unique(lastal.res.cut$lncrna)

# add miRNAs with high confidence even if there's not enough evidence for MTI
high_conf_mir = fread('data/mature_high_conf_v22_1.txt', header = FALSE)


# add info about targets of the mirnas from literature
mirtarbase_valid <- as.data.table(read_excel('data/hsa_MTI.xlsx'))[(`Species (miRNA)`=='Homo sapiens'
                                                                   &
                                                                     `Support Type` %in% c('Functional MTI',
                                                                                           'Functional MTI (Weak)')) |
                                                                     miRNA %in% high_conf_mir$V1]
lnc.mir.tar <- unique(merge(x = lastal.res.cut,
                           y = mirtarbase_valid[, c("miRNA", "Target Gene")],
                           by.x = 'mature.mirna',
                           by.y = "miRNA",
                           allow.cartesian = T))


lnc.mir.tar.fantom <- merge(x = lnc.mir.tar,
                            y = fantom_DE,
                            by.x = c('lncrna', 'Target Gene'),
                            by.y = c('KD.geneSymbol', 'geneSymbol'))
lnc.mir.tar.fantom <- 
  unique(lnc.mir.tar.fantom[,c('lncrna', 'mature.mirna', 'Target Gene', 
                               'log2FC', 'KD.geneID')])
# write.table(lnc.mir.tar.fantom, 'out/lnc_mir_prediction/lnc.mir.tar.fantom.mir.included.tsv')

num_mirnas_per_target <- 
  lnc.mir.tar.fantom[, length(unique(.SD$mature.mirna)), by = 'Target Gene']
histogram(num_mirnas_per_target$V1, 
          type = 'count',
          nint = 100, 
          main = 'Distribution of miRNAs per target',
          xlab = 'Number of miRNAs that sponge the target',
          ylab = 'Number of targets with this number of miRNAs',
          col = "#AEAEEE96",
          scales=list(tick.number = 20, cex = 1),
          xlim = c(-10, 70))

num_mirnas_per_lnc <- lnc.mir.tar.fantom[, length(unique(.SD$mature.mirna)), by = "lncrna"]
setnames(num_mirnas_per_lnc, 'V1', 'num_mirnas_per_lncrna')
setorder(num_mirnas_per_lnc, num_mirnas_per_lncrna)
write.table(num_mirnas_per_lnc, 'out/lnc_mir_prediction/mirnas_per_lncrna_EXACT_TRANSCRIPTS.tsv', 
			quote = F, row.names = F, sep = '\t')

histogram(num_mirnas_per_lnc$num_mirnas_per_lncrna, 
          type = 'count',
          nint = 100, 
          main = 'Distribution of predicted miRNA targets per lncRNA',
          xlab = 'Number of predicted miRNA targets per lncRNA',
          col = "#EEAEEE96")

num_mirnas_per_lnc$lncrna <- factor(num_mirnas_per_lnc$lncrna, 
                           levels = unique(num_mirnas_per_lnc$lncrna))


# dotplot num of predicted mirnas per lncrna
png('out/lnc_mir_prediction/mirnas_per_lncrna_EXACT_TRANSCRIPTS_SUPPORT_TYPE_FILTER.png', width = 2100)
mypanel<-function(x,y,...){
  panel.dotplot(x, y, ...)
  panel.text(x = c(1:78)[1:100],
             y = as.numeric(num_mirnas_per_lnc$num_mirnas_per_lncrna)[1:78]+10,
             labels=num_mirnas_per_lnc$num_mirnas_per_lncrna[1:78],
             cex = 0.85)
  panel.text(x = c(1:137)[79:133][is.even(c(79:133))],
             y = as.numeric(num_mirnas_per_lnc$num_mirnas_per_lncrna)[79:133][is.even(c(79:133))]-10,
             labels=num_mirnas_per_lnc$num_mirnas_per_lncrna[79:133][is.even(c(79:133))],
             cex = 0.8)
  panel.text(x = c(1:133)[79:133][!is.even(c(79:133))]-0.1,
             y = as.numeric(num_mirnas_per_lnc$num_mirnas_per_lncrna)[79:133][!is.even(c(79:133))]+10,
             labels=num_mirnas_per_lnc$num_mirnas_per_lncrna[79:133][!is.even(c(79:133))],
             cex = 0.8)
}

dotplot(num_mirnas_per_lncrna ~ lncrna, num_mirnas_per_lnc,
        main = list(label = 'Distribution of predicted miRNA targets per lncRNA',
                    cex = 1.2),
        ylab = list(label = 'Number of predicted miRNAs per lncRNA',
                    cex = 1),
        scales=list(y=list(rot=0, tick.number = 20, cex = 1), 
                    x=list(rot=45, cex = 1),
                    relation = 'free'),
        ylim = c(0,150),
        panel=mypanel)
dev.off()

lnc.tar.fantom <- unique(lnc.mir.tar.fantom[,c('lncrna', 'Target Gene', 'log2FC', 'KD.geneID')])
boxplot(lnc.tar.fantom$log2FC)
# write.table(lnc.tar.fantom, 'out/lnc_mir_prediction/lnc.tar.fantom.tsv', 
# 	quote = F, row.names = F, sep = '\t')

fantom_DE_not_sponge <- fantom_DE[!(KD.geneID %in% lnc.tar.fantom$KD.geneID &
                                      geneSymbol %in% lnc.tar.fantom$`Target Gene`)]

# scatterplot num of predicted mirnas vs 
# num of targets of the mirna in fatnom per lncrna 

num_mirnas_tar_per_lnc <- 
  merge(x = num_mirnas_per_lnc,
        y = lnc.tar.fantom[, length(unique(.SD$`Target Gene`)), 
                               by = "lncrna"],
        by = "lncrna")
setnames(num_mirnas_tar_per_lnc, 'V1', 'num_targets_per_lnc')

# scatterplot num of predicted mirnas vs 
# num of targets of the mirna in fatnom per lncrna 
PrettyScatter(x = num_mirnas_tar_per_lnc$num_targets_per_lnc,  
              y = num_mirnas_tar_per_lnc$num_mirnas_per_lncrna,
              panel.first.step = 500,
              bg = "#0000FF7D",
              cex.lab = 1,
              main = '',
              xlab = 'Number of FANTOM targets per lncRNA among targets of miRNA',
              ylab = 'Number of predicted miRNAs per lncRNA')

# beautiful pictures
#lnc.tar.fantom <- fread('out/lnc_mir_prediction/lnc.tar.fantom.tsv')
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

# png('out/mirtarbase+lastal/boxplot.mir.targets.vs.all.lnc.targets.png', width = 700)
# DrawTwoBoxplots(data,
#                 main = "All targets of miRNAs sponged by lncRNAs vs all targets of lncRNAs\np-value < 2.2e-16")
# dev.off()
# 
# # compare only DE targets of mirna sponged by lncrna & DE targets of non-sponge lncrna
# data <- data.table(ceRNA = T,
#                    log2FC = lnc.tar.fantom[log2FC > 1 | log2FC < -1]$log2FC)
# data <- rbindlist(list(data,
#                        data.table(ceRNA = F,
#                                   log2FC = fantom_DE[log2FC > 1 | log2FC < -1]$log2FC)))
# data$ceRNA <- as.factor(data$ceRNA)
# levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'all', 'a sponge')
# wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'all']$log2FC, paired = F)
# 
# png('out/mirtarbase+lastal/boxplot.DE.mir.targets.vs.DE.all.lnc.targets.png', width = 700)
# DrawTwoBoxplots(data,
#                 main = "DE targets of miRNAs sponged by lncRNAs vs DE targets of non-sponge lncRNAs\np-value < 2.2e-16")
# dev.off()
# 
# # all targets of fantom lncrnas except targets of sponged mirnas
data <- data.table(ceRNA = T,
                   log2FC =lnc.tar.fantom$log2FC)
data <- rbindlist(list(data,
                       data.table(ceRNA = F,
                                  log2FC = fantom_DE_not_sponge$log2FC)))
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'not a sponge', 'a sponge')
wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'not a sponge']$log2FC, paired = F)

# #ylim1 = boxplot.stats(data$log2FC)$stats[c(1, 5)]
# 
# png('out/mirtarbase+lastal/boxplot.mir.targets.vs.all.targets.non-sponge.lnc.png', width = 700)
# DrawTwoBoxplots(data,
#                 main = "All targets of miRNAs sponged by lncRNAs vs all targets of lncRNAs\np-value < 2.2e-16")
# dev.off()


# all three in one
list.val = list(
  all = fantom_DE$log2FC,
  `not targeted by\nsponged miRNAs` = fantom_DE_not_sponge$log2FC,
  `targeted by\nsponged miRNAs` = lnc.tar.fantom$log2FC,
  `#(mir_per_lnc) > 20` =
    lnc.tar.fantom[lncrna %in% num_mirnas_per_lnc[num_mirnas_per_lncrna > 20]$lncrna]$log2FC,
  `#(mir_per_lnc) > 50` =
    lnc.tar.fantom[lncrna %in% num_mirnas_per_lnc[num_mirnas_per_lncrna > 50]$lncrna]$log2FC,
  `#(mir_per_lnc) > 100` =
    lnc.tar.fantom[lncrna %in% num_mirnas_per_lnc[num_mirnas_per_lncrna > 100]$lncrna]$log2FC,
  `#(mir_per_lnc) > 20 &\n#(targets_per_lnc)>20` =
    lnc.tar.fantom[lncrna %in% num_mirnas_tar_per_lnc[num_mirnas_per_lncrna > 20 &
                                                               num_targets_per_lnc > 20]$lncrna]$log2FC,
  `#(mir_per_lnc) > 50 &\n#(targets_per_lnc)>50` =
    lnc.tar.fantom[lncrna %in% num_mirnas_tar_per_lnc[num_mirnas_per_lncrna > 50 &
                                                               num_targets_per_lnc > 50]$lncrna]$log2FC,
  `#(mir_per_lnc) > 100 &\n#(targets_per_lnc)>100` =
    lnc.tar.fantom[lncrna %in% num_mirnas_tar_per_lnc[num_mirnas_per_lncrna > 100 &
                                                               num_targets_per_lnc > 100]$lncrna]$log2FC)
png('out/mirtarbase+lastal/boxplots.change.num.targets.thresholds.SUPPORT_TYPE_FILTER.png', width = 1400, height = 600)
DrawMultipleBoxplots(list.val, ylab = 'logFC', main = '', ylim = c(min(unlist(list.val)), 11),
                     cex.lab.x = 0.9, xlim = c(0.8,3.9))
dev.off()


list.val = list(
  all = fantom_DE$log2FC,
  `not targeted by\nsponged miRNAs` = fantom_DE_not_sponge$log2FC,
  `targeted by\nsponged miRNAs` = lnc.tar.fantom$log2FC)
png('out/mirtarbase+lastal/three.boxplots.signif.SUPPORT_TYPE_FILTER.png', width = 700, height = 600)
at.x.by = 0.5
DrawMultipleBoxplots(list.val, ylab = 'logFC', main = '', ylim = c(min(unlist(list.val)), 11),
                  cex.lab.x = 1.1, at.x.by = at.x.by, xlim = c(0.8, 1.2+2*at.x.by), proportion.custom = T, proportion = rep(1/3,3))
y <- max(list.val$all)+1
# set an offset for tick lengths
offset <- 0.2
# draw first horizontal line
lines(c(1+at.x.by, 1 + 2*at.x.by),c(y, y))
# draw ticks
lines(c(1+at.x.by, 1+at.x.by),c(y, y-offset))
lines(c(1+2*at.x.by, 1+2*at.x.by),c(y, y-offset))
text(1+3*at.x.by/2,y+0.5,"***")

y <- max(list.val$all)+3
# set an offset for tick lengths
offset <- 0.2
# draw first horizontal line
lines(c(1,1+at.x.by, 1+2*at.x.by),c(y, y, y))
# draw ticks
lines(c(1,1),c(y, y-offset))
lines(c(1+2*at.x.by,1+2*at.x.by),c(y, y-offset))
text(1+at.x.by,y+0.5,"***")
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
png('out/mirtarbase+lastal/mosaic.SUPPORT_TYPE_FILTER.png', 700)
mosaic(struct, shade = TRUE, direction = "v", pop = FALSE)
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
dev.off()



################# group by ncRNA and find out which one has the largest effect
lnc.mir.tar.fantom.l <- split(lnc.mir.tar.fantom, 
                              lnc.mir.tar.fantom$lncrna)
la <- sapply(names(lnc.mir.tar.fantom.l), function(lncrna){
  print(lncrna)
  dt.i <- unique(lnc.mir.tar.fantom.l[[lncrna]][,c("Target Gene", "log2FC")])
  data <- data.table(ceRNA = T,
                     log2FC =dt.i$log2FC)
  dt.j <- fantom_DE[KD.geneSymbol == lncrna & 
                      !(geneSymbol %in% dt.i$`Target Gene`)][,c("geneSymbol", "log2FC")]
  if(dt.i[,.N]>3 & dt.j[,.N]>3){
    data <- rbindlist(list(data,
                           data.table(ceRNA = F,
                                      log2FC = dt.j$log2FC)))
    data$ceRNA <- as.factor(data$ceRNA)
    levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 
                                 'not targeted by\nsponged miRNAs', 
                                 'targeted by\nsponged miRNAs')
    wilcox.test(data[ceRNA == 'targeted by\nsponged miRNAs']$log2FC,
                data[ceRNA == 'not targeted by\nsponged miRNAs']$log2FC, 
                paired = F)$`p.value`
  }else
    NULL
}, simplify = F, USE.NAMES = T)
la <- Filter(Negate(is.null), la)
pvals.wilcox <- data.table(lncrna = names(la),
                           pval = unlist(la))
pvals.wilcox$FDR <- p.adjust(pvals.wilcox$pval, 'BH')

thresh = 0.05
pvals.wilcox <- pvals.wilcox[FDR<thresh]
signif.lncrnas <- pvals.wilcox$lncrna
# "CDKN2B-AS1" 1 value
#lncrna <- names(lnc.mir.tar.fantom.l)[1]

library(rcompanion)
# ATTENTION! With CIs works really slow, takes about 0.5h:
eff.size <- sapply(signif.lncrnas, function(lncrna){
  print(lncrna)
  dt.i <- unique(lnc.mir.tar.fantom.l[[lncrna]][,c("Target Gene", "log2FC")])
  data <- data.table(ceRNA = T,
                     log2FC =dt.i$log2FC)
  dt.j <- fantom_DE[KD.geneSymbol == lncrna & 
                      !(geneSymbol %in% dt.i$`Target Gene`)][,c("geneSymbol", "log2FC")]
  if(dt.i[,.N]>3 & dt.j[,.N]>3){
    data <- rbindlist(list(data,
                           data.table(ceRNA = F,
                                      log2FC = dt.j$log2FC)))
    data$ceRNA <- as.factor(data$ceRNA)
    levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 
                                 'not targeted by\nsponged miRNAs', 
                                 'targeted by\nsponged miRNAs')
    # data$ceRNA <- factor(data$ceRNA, levels=rev(levels(data$ceRNA)))
    eff.size.calc = wilcoxonR(data$log2FC, data$ceRNA, ci = TRUE)
    vda = vda(x = data[ceRNA == 'not targeted by\nsponged miRNAs']$log2FC, 
              y = data[ceRNA == 'targeted by\nsponged miRNAs']$log2FC,
              ci=TRUE)
    data.table(lncrna = lncrna,
              r = eff.size.calc[,1],
              r.lower.ci = eff.size.calc[,2],
              r.upper.ci = eff.size.calc[,3],
              vda = vda[,1],
              vda.lower.ci = vda[,2],
              vda.upper.ci = vda[,3],
              sample.size = data[,.N])
    
  }else
    NULL
}, simplify = F, USE.NAMES = T)
eff.size.dt <- rbindlist(eff.size)



histogram(eff.size.dt$r, 
          type = 'count',
          nint = 10, 
          xlab = 'Effect size',
          ylab = 'Number of lncRNAs',
          col = "#AEAEEE96")
PrettyScatter(x = eff.size.dt$r, 
              y = eff.size.dt$sample.size, 
              bg = "#0000FF7D", 
              main = '', 
              xlab = 'effect size', 
              ylab = 'sample size', 
              cex.lab = 1)

setorder(eff.size.dt, -vda)

eff.size.add <- sapply(signif.lncrnas, function(lncrna){
    print(lncrna)
    dt.i <- unique(lnc.mir.tar.fantom.l[[lncrna]][,c("Target Gene", "log2FC")])
    data <- data.table(ceRNA = T,
                       log2FC =dt.i$log2FC)
    dt.j <- fantom_DE[KD.geneSymbol == lncrna & 
                        !(geneSymbol %in% dt.i$`Target Gene`)][,c("geneSymbol", "log2FC")]
    if(dt.i[,.N]>3 & dt.j[,.N]>3){
      data <- rbindlist(list(data,
                             data.table(ceRNA = F,
                                        log2FC = dt.j$log2FC)))
      data$ceRNA <- as.factor(data$ceRNA)
      data.table(lncrna = lncrna,
                 sponged.sample.size = data[ceRNA == T, .N],
                 not.sponged.sample.size = data[ceRNA == F, .N])
      
    }else
      NULL
  }, simplify = F, USE.NAMES = T)
eff.size.add.dt <- rbindlist(eff.size.add)
eff.size.full <- merge(eff.size.dt,
                       eff.size.add.dt,
                       by = 'lncrna')
eff.size.full$FDR <- pvals.wilcox[match(eff.size.full$lncrna, lncrna)]$FDR
saveRDS(eff.size.full, file = 'out/eff.size.full.SUPPORT_TYPE_FILTERED_WEAK_SUPPORT_INCLUDED.RDS')

library(gridExtra)
library(grid)
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.8)),
  colhead = list(fg_params=list(cex = 0.8)),
  rowhead = list(fg_params=list(cex = 0.8)))

setorder(eff.size.full, -vda)
myt <- gridExtra::tableGrob(eff.size.full, theme = mytheme)

grid.draw(myt)



# plot all
png('out/142tr_MTI+MTI_weak.png', width = 1000)
par(mfrow=c(4,4))
for(lncrna.i in eff.size.full$lncrna[1:16]){
  dt.T <- unique(lnc.mir.tar.fantom.l[[lncrna.i]][,c("Target Gene", "log2FC")])
  dt.F <- fantom_DE[KD.geneSymbol == lncrna.i & 
                      !(geneSymbol %in% dt.T$`Target Gene`)][,c("geneSymbol", "log2FC")]
  dts.plot = list(not.sponged = dt.F$log2FC,
                  sponged = dt.T$log2FC)
  print(lncrna.i)
  Sys.sleep(0.5)
  DrawMultipleBoxplots(dts.plot, 
                       ylab = 'logFC', 
                       main = paste0(lncrna.i, '\nvda = ', eff.size.full[lncrna == lncrna.i]$vda,
                                     ', FDR = ', round(pvals.wilcox[lncrna == lncrna.i]$FDR, 5)), 
                       ylim = c(-6, 7),
                       cex.lab.x = 0.7, 
                       cex.lab.y = 0.7,
                       xlim = c(0.8,1.6),
                       main.cex = 0.8,
                       at.y.custom = T,
                       at.y = pretty(c(-6:7), 5),
                       at.y.cex = 0.7,
                       at.x.by = .35)
}
dev.off()




library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
require(gridExtra)

png('out/MTI_+weak+high_conf_mir_1.png', width = 1900, height = 1200)
plots = lapply(eff.size.full$lncrna[25:39], function(lncrna.i){
  dt.T <- unique(lnc.mir.tar.fantom.l[[lncrna.i]][,c("Target Gene", "log2FC")])
  dt.F <- fantom_DE[KD.geneSymbol == lncrna.i & 
                      !(geneSymbol %in% dt.T$`Target Gene`)][,c("geneSymbol", "log2FC")]
  
  data <- data.frame(value = c(dt.T$log2FC, dt.F$log2FC),
                     name = c(rep('sponged', length(dt.T$log2FC)), rep('not.sponged', length(dt.F$log2FC))))
  
  
  # sample size
  sample_size = data %>% group_by(name) %>% summarize(num=n())
  
  plot = data %>%
    left_join(sample_size) %>%
    mutate(myaxis = paste0(name, "\n", "n=", num)) %>%
    ggplot( aes(x=myaxis, y=value, fill=name)) +
    geom_violin(width=1) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2, position = position_dodge(width = 0.7)) +
    scale_fill_viridis(discrete = TRUE) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste0(lncrna.i, '\nvda = ', eff.size.full[lncrna == lncrna.i]$vda,
            ', FDR = ', round(pvals.wilcox[lncrna == lncrna.i]$FDR, 5))) +
    xlab("") +
    ylab('logFC')
      
  plot
})
do.call("grid.arrange", c(plots, ncol=8, nrow = 2))
dev.off()
