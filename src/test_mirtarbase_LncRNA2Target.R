library("readxl")
library(data.table)
source('src/func.R')

# mirtarbase_all <- as.data.table(read_excel("data/miRTarBase_SE_WR.xls"))[`Species (miRNA)`=='Homo sapiens' & 
#                                                                            `Species (Target Gene)`=='Homo sapiens']
mirtarbase_valid <- as.data.table(read_excel('data/hsa_MTI.xlsx'))[`Species (miRNA)`=='Homo sapiens']
#?? why so many validated but less in all?..

fantom_DE <- fread('../oligo_DE_Summary_gene_filtered.tsv')
fantom_DE$geneSymbol <- stringr::str_remove(fantom_DE$geneSymbol, 'HG')
#fantom_DE.mirna.target <- fantom_DE[grep('mir', geneSymbol,ignore.case = T)]


# experimentally confirmed from DB 2019  
lncRNA.miRNA.interaction <- fread('data/lncTarD.txt')[grep('mir', Target, ignore.case = T), c('Regulator', 'RegulatorEnsembleID', 'RegulatorAliases', 'Target')]

# leave only lncRNAs present in FANTOM DE
lncrna.mirna.int.fantom <- lncRNA.miRNA.interaction[RegulatorEnsembleID %in% fantom_DE$KD.geneID]
lncrna.mirna.int.fantom$Target <- paste0('hsa-', lncrna.mirna.int.fantom$Target)


# library("biomaRt")
# listMarts()
# ensembl <- useMart("ensembl")
# datasets <- listDatasets(ensembl)
# ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# attributes = listAttributes(ensembl)
# lnc2target.NOT.in.fantom.biomart <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
#                                filters = 'ensembl_gene_id',
#                                values = unique(lnc2target.NOT.in.fantom$RegulatorEnsembleID),
#                                mart = ensembl)
# lnc2target.NOT.in.fantom.biomart <- as.data.table(lnc2target.NOT.in.fantom.biomart)
# setnames(lnc2target.NOT.in.fantom.biomart, colnames(lnc2target.NOT.in.fantom.biomart),
#          c('RegulatorEnsembleID', 'biomart_hgnc'))
# lnc2target.NOT.in.fantom.biomart <- merge(lnc2target.NOT.in.fantom,
#                                           lnc2target.NOT.in.fantom.biomart,
#                                           by = 'RegulatorEnsembleID')
# 
# 
# fantomDE.biomart <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
#                                           filters = 'ensembl_gene_id',
#                                           values = unique(fantom_DE$KD.geneID),
#                                           mart = ensembl)
# fantomDE.biomart <- as.data.table(fantomDE.biomart)
# setnames(fantomDE.biomart, colnames(fantomDE.biomart), c('KD.geneID', 'biomart_hgnc'))
# fantomDE.biomart <- merge(fantomDE.biomart,
#                           unique(fantom_DE[, c('KD.geneID',"KD.geneSymbol")]))
# fantomDE.biomart[biomart_hgnc !=KD.geneSymbol & biomart_hgnc !='']$biomart_hgnc %in% lnc2target.NOT.in.fantom.biomart$biomart_hgnc
# 
# 
# 
# 
# 
# 
# 
# 
# 
# lncRNA.knownCE.in.fantom <- unique(lncRNA.miRNA.interaction[RegulatorEnsembleID %in% fantom_DE$KD.geneID]$Regulator)
# lncRNA.knownCE.not.in.fantom <- unique(lncRNA.miRNA.interaction[!RegulatorEnsembleID %in% fantom_DE$KD.geneID]$Regulator)
# 
# lncrna_aliases <- str_split_fixed(lncRNA.miRNA.interaction$RegulatorAliases[!is.na(lncRNA.miRNA.interaction$RegulatorAliases)], pattern = '\\|', 11)
# lncrna_aliases <- unique(as.vector(lncrna_aliases))
# lncrna_aliases <- lncrna_aliases[!lncrna_aliases == ""]
# lncRNA.knownCE.in.fantom2 <- lncrna_aliases[lncrna_aliases %in% fantom_DE$KD.geneSymbol]
# aliases_search <- sapply(lncrna_aliases, function(alias){
#   print(alias)
#   unique(fantom_DE[grep(alias, geneSymbol,ignore.case = T)]$geneSymbol)
# }, USE.NAMES = T, simplify = F)
# aliases_search <- lapply(aliases_search, function(x){if(length(x)>0){x}})
# aliases_search <- Filter(Negate(is.null), aliases_search)
# lncRNA.knownCE.in.fantom.to.add <- c("HOXD-AS1", "ZEB2-AS", "ILF3-AS1", "LINC00152", "CDKN2B-AS", "NME1")
# #
# lncRNA.knownCE.in.fantom <- c(lncRNA.knownCE.in.fantom,
#                               lncRNA.knownCE.in.fantom2,
#                               lncRNA.knownCE.in.fantom.to.add)
# 
# lncrna.mirna.int.fantom <-
#   rbindlist(sapply(lncRNA.knownCE.in.fantom, function(alias){
#   res <- lncRNA.miRNA.interaction[grep(alias, RegulatorAliases)]
#   if(res[,.N] == 0){
#     res <- lncRNA.miRNA.interaction[grep(alias, Regulator)]
#   }
#   data.table(lncrna = alias,
#              mirna = res$Target)
# },USE.NAMES = T, simplify = F))
# lncrna.mirna.int.fantom$Target <- paste0('hsa-', lncrna.mirna.int.fantom$mirna)


ensembl.mirna.corresp.mirbase_id <- getBM(attributes=c('ensembl_gene_id', 'mirbase_id', 'mirbase_accession', 'hgnc_symbol'),
                                           filters = 'mirbase_id',
                                           values = unique(lncrna.mirna.int.fantom$Target),
                                           mart = ensembl)

not.found.mirna <- lncrna.mirna.int.fantom$Target[
  !tolower(lncrna.mirna.int.fantom$Target) %in% ensembl.mirna.corresp.mirbase_id$mirbase_id]

not.found.mirna %in% mirtarbase_valid$miRNA
#
# 
mature.mirna.ids <- unique(lncrna.mirna.int.fantom$Target)[grep('-3p', unique(lncrna.mirna.int.fantom$Target))]
mature.mirna.ids2 <- unique(lncrna.mirna.int.fantom$Target)[grep('-5p', unique(lncrna.mirna.int.fantom$Target))]
mature.mirna.ids <- c(mature.mirna.ids, mature.mirna.ids2)
rm(mature.mirna.ids2)
to.find.mature.id <- unique(lncrna.mirna.int.fantom$Target)[!unique(lncrna.mirna.int.fantom$Target) %in% mature.mirna.ids]

version=checkMiRNAVersion(to.find.mature.id, verbose = TRUE)
mirnatable <- as.data.table(getMiRNATable(version = "v22", species = "hsa"))
mirnatable <- mirnatable[Precursor %in% tolower(to.find.mature.id), c('Precursor', 'Mature1', 'Mature2')]
mirnatable.mature1 <- mirnatable[,c('Precursor', 'Mature1')]
setnames(mirnatable.mature1, 'Mature1', 'mature')
mirnatable.mature2 <- mirnatable[,c('Precursor', 'Mature2')]
setnames(mirnatable.mature2, 'Mature2', 'mature')
mirnatable_long <- rbindlist(list(mirnatable.mature1, mirnatable.mature2))[!is.na(mature)]
mirnatable_long <- rbindlist(list(mirnatable_long,
                                  data.table(Precursor = mature.mirna.ids,
                                             mature = mature.mirna.ids)))
setnames(mirnatable_long, colnames(mirnatable_long), c("mirna.lnc.inter.db_id", "mature_mirna"))

to.check.mature.id <- tolower(to.find.mature.id)[!tolower(to.find.mature.id) %in% mirnatable$Precursor]
# # library(miRNAmeConverter)
# # MiRNANameConverter <- MiRNANameConverter()
# # translateMiRNAName(MiRNANameConverter, miRNAs = "hsa-mir-9", versions = '17')
# 

mirna_id.corresp <- sapply(to.check.mature.id, function(mirna){
  print(mirna)
  search.res <- unique(mirtarbase_valid[grep(mirna, miRNA, ignore.case = T)]$miRNA)
  if(length(search.res)>0){
    data.table(mirna.lnc.inter.db_id = mirna,
               mature_mirna = search.res)
  }else
    NULL
},USE.NAMES = T, simplify = F)
mirna_id.corresp <- rbindlist(Filter(Negate(is.null), mirna_id.corresp))

mirna_id.corresp.fixed <- data.table(
  mirna.lnc.inter.db_id = c('hsa-mir-1', 'hsa-mir-1', 'hsa-mir-153', 'hsa-mir-153', 
                            'hsa-mir-147', 'hsa-mir-147', 'hsa-mir-7', 'hsa-mir-7', 'hsa-mir-7'),
  mature_mirna = c('hsa-miR-1-3p', 'hsa-miR-1-5p', 'hsa-miR-153-3p', 'hsa-miR-153-5p',
                   'hsa-miR-147a', 'hsa-miR-147b', 'hsa-miR-7-1-3p', 'hsa-miR-7-2-3p', 'hsa-miR-7-5p'))
mirna.correct <- c('hsa-mir-29b', 'hsa-mir-199a', 'hsa-mir-125b', 'hsa-mir-218', 'hsa-mir-203')
mirna_id.corresp.fixed <- rbindlist(list(mirna_id.corresp.fixed,
                                         mirna_id.corresp[mirna.lnc.inter.db_id %in% mirna.correct]))
mirna_id.corresp.fixed <- rbindlist(list(mirnatable_long, mirna_id.corresp.fixed))
mirna_id.corresp <- mirna_id.corresp.fixed

lncrna.mirna.int.fantom <- merge(x = lncrna.mirna.int.fantom,
                                 y = mirna_id.corresp,
                                 by.x = 'Target',
                                 by.y = "mirna.lnc.inter.db_id")
lncrna.mirna.int.fantom <- lncrna.mirna.int.fantom[, c('Regulator', 'RegulatorEnsembleID', 'mature_mirna')]




# ensemble.mirna <- rbindlist(list(ensembl.mirna.corresp.gene.symbol,
#                                  ensembl.mirna.corresp.mirbase_id))
# rm(ensembl.mirna.corresp.gene.symbol, ensembl.mirna.corresp.mirbase_id)

# library(miRBaseConverter)
# version=checkMiRNAVersion(ensemble.mirna$mirbase_id, verbose = TRUE)
# mirnatable <- as.data.table(getMiRNATable(version = "v22", species = "hsa"))
# mirnatable <- mirnatable[Precursor %in% ensemble.mirna$mirbase_id, c('Precursor', 'Mature1', 'Mature2')]
# mirnatable.mature1 <- mirnatable[,c('Precursor', 'Mature1')]
# setnames(mirnatable.mature1, 'Mature1', 'mature')
# mirnatable.mature2 <- mirnatable[,c('Precursor', 'Mature2')]
# setnames(mirnatable.mature2, 'Mature2', 'mature')
# mirnatable_long <- rbindlist(list(mirnatable.mature1, mirnatable.mature2))[!is.na(mature)]
# mirtable_long_all_ids <- merge(x = mirnatable_long,
#                                y = ensemble.mirna,
#                                by.x = 'Precursor',
#                                by.y = 'mirbase_id')

lnc.mir.targets <- merge(x = mirtarbase_valid,
                     y = lncrna.mirna.int.fantom,
                     by.x = 'miRNA',
                     by.y = 'mature_mirna',
                     allow.cartesian = T)
lnc.mir.targets <- 
  unique(lnc.mir.targets[,c('Regulator', 'RegulatorEnsembleID', 'miRNA', 'Target Gene')])
setnames(lnc.mir.targets, 
         c('Regulator', 'RegulatorEnsembleID', 'miRNA'), 
         c('lncrna', 'lncrna.ensemble.id', 'mirna'))

lncRNA.mirna.target.logFC <- merge(x = lnc.mir.targets,
                                   y = fantom_DE[, c("KD.geneID", "KD.geneSymbol", 'geneSymbol', 'log2FC')],
                                   by.x = c('lncrna.ensemble.id', 'Target Gene'),
                                   by.y = c("KD.geneID", 'geneSymbol'))


fantom_DE_not_sponge.exp <- fantom_DE[!(KD.geneID %in% lncRNA.mirna.target.logFC $lncrna.ensemble.id &
                                      geneSymbol %in% lncRNA.mirna.target.logFC $mir.target)]

list.val = list(
  all = fantom_DE$log2FC,
  `not targeted by\nsponged miRNAs` = fantom_DE_not_sponge$log2FC,
  `targeted by\nsponged miRNAs` = lncRNA.mirna.target.logFC$log2FC)
wilcox.test(list.val$`not targeted by\nsponged miRNAs`, 
            list.val$`targeted by\nsponged miRNAs`, paired = F)
wilcox.test(list.val$all,
            list.val$`targeted by\nsponged miRNAs`, paired = F)

png('out/mirtarbase+lnctard/three.boxplots.signif.png', width = 700, height = 600)
at.x.by = 0.5
DrawMultipleBoxplots(list.val, ylab = 'logFC', main = '', ylim = c(min(unlist(list.val)), 11),
                     cex.lab.x = 1.1, at.x.by = at.x.by, xlim = c(0.8, 1.2+2*at.x.by),
                     proportion.custom = T, proportion = rep(1/3, 3))

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

################# group by ncRNA and find out which one has the largest effect
lnc.mir.tar.fantom.l <- split(lncRNA.mirna.target.logFC, 
                              lncRNA.mirna.target.logFC$lncrna)
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
#saveRDS(eff.size.full, file = 'out/eff.size.full.SUPPORT_TYPE_FILTERED.RDS')

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
par(mfrow=c(1,3))
for(lncrna.i in eff.size.full$lncrna[1:3]){
  dt.T <- unique(lnc.mir.tar.fantom.l[[lncrna.i]][,c("Target Gene", "log2FC")])
  dt.F <- fantom_DE[KD.geneSymbol == lncrna.i & 
                      !(geneSymbol %in% dt.i$`Target Gene`)][,c("geneSymbol", "log2FC")]
  dts.plot = list(not.sponged = dt.F$log2FC,
                  sponged = dt.T$log2FC)
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




