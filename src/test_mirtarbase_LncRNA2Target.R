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
lncRNA.miRNA.interaction <- as.data.table(read_excel('data/lnc2target.xlsx'))[
  `Species`=='Homo sapiens'][
    grep('mir', Target_official_symbol,ignore.case = T), c("lncRNA_name_from_paper",
                                                           "LncRNA_official_symbol",
                                                           "Ensembl_ID",
                                                           "Target_symbol_from_paper",
                                                           "Target_official_symbol",
                                                           "Target_entrez_gene_ID")]
setnames(lncRNA.miRNA.interaction, c("LncRNA_official_symbol", "Ensembl_ID"), 
                                   c('Regulator', 'RegulatorEnsembleID'))
lncRNA.miRNA.interaction$lncRNA_name_from_paper <- toupper(lncRNA.miRNA.interaction$lncRNA_name_from_paper)
lnc2target.int <- intersect(toupper(unique(lncRNA.miRNA.interaction$lncRNA_name_from_paper)), 
                            unique(lncRNA.miRNA.interaction$Regulator))
lnc2target.lncrnas <- c(unique(lncRNA.miRNA.interaction$lncRNA_name_from_paper)[
  !unique(lncRNA.miRNA.interaction$lncRNA_name_from_paper) %in% lnc2target.int],
  unique(lncRNA.miRNA.interaction$Regulator)[
    !unique(lncRNA.miRNA.interaction$Regulator) %in% lnc2target.int],
  lnc2target.int)
lnc2target.NOT.in.fantom <- unique(lncRNA.miRNA.interaction[!RegulatorEnsembleID %in% fantom_DE$KD.geneID,
                                                     c('Regulator', 'RegulatorEnsembleID')])[
                                                       Regulator != 'NA' & RegulatorEnsembleID != 'NA']
library("biomaRt")
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)
lnc2target.NOT.in.fantom.biomart <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                               filters = 'ensembl_gene_id',
                               values = unique(lnc2target.NOT.in.fantom$RegulatorEnsembleID),
                               mart = ensembl)
lnc2target.NOT.in.fantom.biomart <- as.data.table(lnc2target.NOT.in.fantom.biomart)
setnames(lnc2target.NOT.in.fantom.biomart, colnames(lnc2target.NOT.in.fantom.biomart), 
         c('RegulatorEnsembleID', 'biomart_hgnc'))
lnc2target.NOT.in.fantom.biomart <- merge(lnc2target.NOT.in.fantom,
                                          lnc2target.NOT.in.fantom.biomart,
                                          by = 'RegulatorEnsembleID')


fantomDE.biomart <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                                          filters = 'ensembl_gene_id',
                                          values = unique(fantom_DE$KD.geneID),
                                          mart = ensembl)
fantomDE.biomart <- as.data.table(fantomDE.biomart)
setnames(fantomDE.biomart, colnames(fantomDE.biomart), c('KD.geneID', 'biomart_hgnc'))
fantomDE.biomart <- merge(fantomDE.biomart,
                          unique(fantom_DE[, c('KD.geneID',"KD.geneSymbol")]))
fantomDE.biomart[biomart_hgnc !=KD.geneSymbol & biomart_hgnc !='']$biomart_hgnc %in% lnc2target.NOT.in.fantom.biomart$biomart_hgnc









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
# 
# lncRNA.knownCE.in.fantom <- c(lncRNA.knownCE.in.fantom,
#                               lncRNA.knownCE.in.fantom2,
#                               lncRNA.knownCE.in.fantom.to.add)

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


# ensembl.mirna.corresp.mirbase_id <- getBM(attributes=c('ensembl_gene_id', 'mirbase_id', 'mirbase_accession', 'hgnc_symbol'), 
#                                            filters = 'mirbase_id', 
#                                            values = unique(lncrna.mirna.int.fantom$Target), 
#                                            mart = ensembl)
# 
# 
# not.found.mirna <- lncrna.mirna.int.fantom$Target[
#   !tolower(lncrna.mirna.int.fantom$Target) %in% ensembl.mirna.corresp.mirbase_id$mirbase_id]
# 
# not.found.mirna %in% mirtarbase_valid$miRNA
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
lncrna.mirna.int.fantom <- lncrna.mirna.int.fantom[, c('lncrna', 'mature_mirna')]




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
  unique(lnc.mir.targets[,c('lncrna', 'miRNA', 'Target Gene')])
setnames(lnc.mir.targets, 
         colnames(lnc.mir.targets), 
         c('lncrna', 'mirna', 'mir.target'))

lncRNA.mirna.target.logFC <- merge(x = lnc.mir.targets,
                                   y = fantom_DE[, c("KD.geneSymbol", 'geneSymbol', 'log2FC')],
                                   by.x = c("lncrna", 'mir.target'),
                                   by.y = c("KD.geneSymbol", 'geneSymbol'))

setnames(lncRNA.mirna.target.logFC, 'log2FC', 'mir_target_log2FC')
boxplot(lncRNA.mirna.target.logFC$mir_target_log2FC)
boxplot(lncRNA.mirna.target.logFC[mir_target_log2FC>1 | mir_target_log2FC < -1]$mir_target_log2FC,
        fantom_DE$log2FC)
boxplot(lncRNA.mirna.target.logFC[mir_target_log2FC>1 | mir_target_log2FC < -1]$mir_target_log2FC,
        fantom_DE[log2FC>1 | log2FC < -1]$log2FC)


fantom_DE[, ceRNA:=ifelse(KD.geneSymbol %in% lncRNA.knownCE.in.fantom, T, F)]
fantom_DE.not.sponge <- fantom_DE[ceRNA == F]
fantom_DE.sponge <- fantom_DE[ceRNA == T]

boxplot(lncRNA.mirna.target.logFC$mir_target_log2FC,
        fantom_DE.not.sponge$log2FC)
boxplot(lncRNA.mirna.target.logFC[mir_target_log2FC>1 | mir_target_log2FC < -1]$mir_target_log2FC,
        fantom_DE.not.sponge[log2FC>1 | log2FC < -1]$log2FC)
wilcox.test(lncRNA.mirna.target.logFC[mir_target_log2FC>1 | mir_target_log2FC < -1]$mir_target_log2FC,
            fantom_DE.not.sponge[log2FC>1 | log2FC < -1]$log2FC,
            paired = F)

# chi squared
contingency.table.sponge <- data.table('rise' = sum(lncRNA.mirna.target.logFC$mir_target_log2FC>1),
                                       'fall' = sum(lncRNA.mirna.target.logFC$mir_target_log2FC < -1))
contingency.table.not.sponge <- data.table('rise' = sum(fantom_DE.not.sponge$log2FC>1),
                                           'fall' = sum(fantom_DE.not.sponge$log2FC < -1))
contingency.table <- rbindlist(list(contingency.table.not.sponge,
                                    contingency.table.sponge))
chisq.test(contingency.table)
df <- contingency.table

library("graphics")
library(viridis)

#png('out/chi.squared.sponges.png', width = 700)
mosaicplot(df,
           main = "Gene expression falls more often in case of non-sponge lncRNAs", 
           sub = paste0('chi squared test p-value = ', round(test.res$p.value,33)),
           cex.axis = 1, 
           color = viridis(4, alpha = 0.5))
dev.off()


# beautiful pictures

# compare logFC of targets of mirna sponged by lncrna & all targets of lncrna
data <- data.table(ceRNA = T,
                   log2FC = lncRNA.mirna.target.logFC$mir_target_log2FC)
data <- rbindlist(list(data,
                       data.table(ceRNA = F,
                                  log2FC = fantom_DE.not.sponge$log2FC)))
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'not a sponge', 'a sponge')
wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'not a sponge']$log2FC, paired = F)

#ylim1 = boxplot.stats(data$log2FC)$stats[c(1, 5)]

png('out/mirtarbase+lnctard/boxplot.mir.targets.vs.all.targets.non-sponge.lnc.png', width = 700)
DrawTwoBoxplots(data, 
                main = "All targets of miRNAs sponged by lncRNAs vs all targets of non-sponge lncRNAs\np-value = 0.0003")
dev.off()

# compare only DE targets of mirna sponged by lncrna & DE targets of non-sponge lncrna
data <- data.table(ceRNA = T,
                   log2FC = lncRNA.mirna.target.logFC[
                     mir_target_log2FC > 1 | mir_target_log2FC < -1]$mir_target_log2FC)
data <- rbindlist(list(data,
                       data.table(ceRNA = F,
                                  log2FC = fantom_DE.not.sponge[log2FC > 1 | log2FC < -1]$log2FC)))
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'not a sponge', 'a sponge')
wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'not a sponge']$log2FC, paired = F)

png('out/mirtarbase+lnctard/boxplot.DE.mir.targets.vs.DE.non-sponge.lnc.targets.png', width = 700)
DrawTwoBoxplots(data, 
                main = "DE targets of miRNAs sponged by lncRNAs vs DE targets of non-sponge lncRNAs\np-value = 0.037")
dev.off()

# compare only DE targets of mirna sponged by lncrna & DE targets all lncrna
data <- data.table(ceRNA = T,
                   log2FC = lncRNA.mirna.target.logFC[
                     mir_target_log2FC > 1 | mir_target_log2FC < -1]$mir_target_log2FC)
data <- rbindlist(list(data,
                       data.table(ceRNA = F,
                                  log2FC = fantom_DE[log2FC > 1 | log2FC < -1]$log2FC)))
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'all', 'a sponge')
wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'all']$log2FC, paired = F)

png('out/mirtarbase+lnctard/boxplot.DE.mir.targets.vs.DE.all.lnc.targets.png', width = 700)
DrawTwoBoxplots(data, 
                main = "DE targets of miRNAs sponged by lncRNAs vs DE targets of all lncRNAs\np-value = 0.044")
dev.off()

# compare all targets of mirna sponged by lncrna & all targets all lncrna
data <- data.table(ceRNA = T,
                   log2FC = lncRNA.mirna.target.logFC$mir_target_log2FC)
data <- rbindlist(list(data,
                       data.table(ceRNA = F,
                                  log2FC = fantom_DE$log2FC)))
data$ceRNA <- as.factor(data$ceRNA)
levels(data$ceRNA) <- ifelse(levels(data$ceRNA) == F, 'all', 'a sponge')
wilcox.test(data[ceRNA == 'a sponge']$log2FC, data[ceRNA == 'all']$log2FC, paired = F)

png('out/mirtarbase+lnctard/boxplot.all.mir.targets.vs.all.lnc.targets.png', width = 700)
DrawTwoBoxplots(data, 
                main = "All targets of miRNAs sponged by lncRNAs vs all targets of all lncRNAs\np-value = 0.0003")
dev.off()
