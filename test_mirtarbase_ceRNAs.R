library("readxl")
library(data.table)

mirtarbase_all <- as.data.table(read_excel("data/miRTarBase_SE_WR.xls"))[`Species (miRNA)`=='Homo sapiens']
mirtarbase_lit <- as.data.table(read_excel('data/MicroRNA_Target_Sites.xlsx'))[`Species (miRNA)`=='Homo sapiens' & 
                                                                                 `Species (Target Gene)`=='Homo sapiens']
# mirtarbase_valid <- as.data.table(read_excel('data/hsa_MTI.xlsx'))[`Species (miRNA)`=='Homo sapiens']
#?? why so many validated but less in all?..

fantom_DE <- fread('../oligo_DE_Summary_gene_filtered.tsv')
fantom_DE$geneSymbol <- stringr::str_remove(fantom_DE$geneSymbol, 'HG')
#fantom_DE.mirna.target <- fantom_DE[grep('mir', geneSymbol,ignore.case = T)]
  
lncRNA.miRNA.interaction <- unique(fantom_DE$KD.geneSymbol)[
  unique(fantom_DE$KD.geneSymbol) %in% mirtarbase_lit$`Target Gene`]











library("biomaRt")
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)
ensembl.mirna.corresp.gene.symbol <- getBM(attributes=c('ensembl_gene_id', 'mirbase_id', 'mirbase_accession', 'hgnc_symbol'), 
                               filters = 'hgnc_symbol', 
                               values = unique(fantom_DE.mirna.target$geneSymbol), 
                               mart = ensembl)
ensembl.mirna.corresp.mirbase_id <- getBM(attributes=c('ensembl_gene_id', 'mirbase_id', 'mirbase_accession', 'hgnc_symbol'), 
                                           filters = 'mirbase_id', 
                                           values = unique(fantom_DE.mirna.target$geneSymbol), 
                                           mart = ensembl)
ensemble.mirna <- rbindlist(list(ensembl.mirna.corresp.gene.symbol,
                                 ensembl.mirna.corresp.mirbase_id))
rm(ensembl.mirna.corresp.gene.symbol, ensembl.mirna.corresp.mirbase_id)

library(miRBaseConverter)
version=checkMiRNAVersion(ensemble.mirna$mirbase_id, verbose = TRUE)
mirnatable <- as.data.table(getMiRNATable(version = "v22", species = "hsa"))
mirnatable <- mirnatable[Precursor %in% ensemble.mirna$mirbase_id, c('Precursor', 'Mature1', 'Mature2')]
mirnatable.mature1 <- mirnatable[,c('Precursor', 'Mature1')]
setnames(mirnatable.mature1, 'Mature1', 'mature')
mirnatable.mature2 <- mirnatable[,c('Precursor', 'Mature2')]
setnames(mirnatable.mature2, 'Mature2', 'mature')
mirnatable_long <- rbindlist(list(mirnatable.mature1, mirnatable.mature2))[!is.na(mature)]
mirtable_long_all_ids <- merge(x = mirnatable_long,
                               y = ensemble.mirna,
                               by.x = 'Precursor',
                               by.y = 'mirbase_id')

mir.all_ids <- merge(x = mirtarbase_lit,
                     y = mirtable_long_all_ids,
                     by.x = 'miRNA',
                     by.y = 'mature')
mir.all_ids <- 
  unique(mir.all_ids[,c('Precursor', "ensembl_gene_id", "mirbase_accession", "hgnc_symbol", 'Target Gene')])
setnames(mir.all_ids, 
         colnames(mir.all_ids), 
         c('mirbase_id', "mir_ensembl", "mirbase_accession", "mir_hgnc", 'mir_target'))
lncRNA.mirna.target1 <- unique(
  merge(x = fantom_DE.mirna.target,
        y = mir.all_ids,
        by.x = c('geneSymbol'),
        by.y = c('mir_hgnc'),
        allow.cartesian = T)[, -c('log2FC', 'geneID', 'mir_ensembl', 'mirbase_accession')])
lncRNA.mirna.target2 <- unique(
  merge(x = fantom_DE.mirna.target,
        y = mir.all_ids,
        by.x = c('geneSymbol'),
        by.y = c('mirbase_id'),
        allow.cartesian = T)[, -c('log2FC', 'geneID', 'mir_ensembl', 'mirbase_accession')])
lncRNA.mirna.target <- rbindlist(list(lncRNA.mirna.target1, lncRNA.mirna.target2))
setnames(lncRNA.mirna.target, 'geneSymbol', 'mir_hgnc')
rm(lncRNA.mirna.target1, lncRNA.mirna.target2)
lncRNA.mirna.target.logFC <- merge(x = lncRNA.mirna.target,
                                   y = fantom_DE[, c("KD.geneID_ASO", 'geneSymbol', 'log2FC')],
                                   by.x = c("KD.geneID_ASO", 'mir_target'),
                                   by.y = c("KD.geneID_ASO", 'geneSymbol'))

setnames(lncRNA.mirna.target.logFC, 'log2FC', 'mir_target_log2FC')
boxplot(lncRNA.mirna.target.logFC$mir_target_log2FC)
boxplot(lncRNA.mirna.target.logFC[mir_target_log2FC>1 | mir_target_log2FC < -1]$mir_target_log2FC,
        fantom_DE$log2FC)
boxplot(lncRNA.mirna.target.logFC[mir_target_log2FC>1 | mir_target_log2FC < -1]$mir_target_log2FC,
        fantom_DE[log2FC>1 | log2FC < -1]$log2FC)


fantom_DE[, ceRNA:=ifelse(KD.geneSymbol %in% lncRNA.mirna.target$KD.geneSymbol, T, F)]
fantom_DE.not.sponge <- fantom_DE[ceRNA == F]

boxplot(lncRNA.mirna.target.logFC$mir_target_log2FC,
        fantom_DE.not.sponge$log2FC)
boxplot(lncRNA.mirna.target.logFC[mir_target_log2FC>1 | mir_target_log2FC < -1]$mir_target_log2FC,
        fantom_DE.not.sponge[log2FC>1 | log2FC < -1]$log2FC)
wilcox.test(lncRNA.mirna.target.logFC[mir_target_log2FC>1 | mir_target_log2FC < -1]$mir_target_log2FC,
            fantom_DE.not.sponge[log2FC>1 | log2FC < -1]$log2FC,
            paired = F)

# chi squared
contingency.table <- as.data.table(table(lncRNA.mirna.target.logFC$))
contingency.table[, expr.change:=ifelse(V1 == T, 'rise', 'fall')]
contingency.table$V1 <- NULL
setnames(contingency.table, 'N', 'sponge')
contingency.table[, not.sponge:= c(as.data.table(table(ceRNA.logFC.dt.cut[ceRNA == 0]$log2FC>1))[V1==F]$N,
                                   as.data.table(table(ceRNA.logFC.dt.cut[ceRNA == 0]$log2FC>1))[V1==T]$N)]
setcolorder(contingency.table, c('expr.change', 'sponge', 'not.sponge'))

# beautiful pictures