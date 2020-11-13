library(data.table)

DIANA_res = fread('data/DIANA1.csv')
mirtarbase_valid <- as.data.table(read_excel('data/hsa_MTI.xlsx'))[`Species (miRNA)`=='Homo sapiens']

lnc.mir.tar <- unique(merge(x = DIANA_res[,c('Gene Name', 'Gene Id', 'Mirna')],
                            y = mirtarbase_valid[, c("miRNA", "Target Gene")],
                            by.x = 'Mirna',
                            by.y = "miRNA",
                            allow.cartesian = T))
setnames(lnc.mir.tar, c('Gene Name', 'Gene Id'), c('lncrna', 'lncrna.id'))

fantom_DE <- fread('../oligo_DE_Summary_gene_filtered.tsv')
fantom_DE$geneSymbol <- stringr::str_remove(fantom_DE$geneSymbol, 'HG')

lnc.mir.tar.fantom <- merge(x = lnc.mir.tar,
                            y = fantom_DE,
                            by.x = c('lncrna', 'Target Gene'),
                            by.y = c('KD.geneSymbol', 'geneSymbol'))
lnc.mir.tar.fantom <- 
  unique(lnc.mir.tar.fantom[,c('lncrna', 'Mirna', 'Target Gene', 
                               'log2FC', 'KD.geneID')])
num_mirnas_per_target <- 
  lnc.mir.tar.fantom[, length(unique(.SD$Mirna)), by = 'Target Gene']

num_mirnas_per_lnc <- lnc.mir.tar.fantom[, length(unique(.SD$Mirna)), by = "lncrna"]


lnc.tar.fantom <- unique(lnc.mir.tar.fantom[,c('lncrna', 'Target Gene', 'log2FC', 'KD.geneID')])





DIANA_found_lnc = unlist(str_split_fixed(unique(lnc.mir.tar$lncrna.id), '\\|', 5)[,1])
unique(fantom_DE[!KD.geneID %in% DIANA_found_lnc]$KD.geneID)

write.table(unique(fantom_DE[!KD.geneID %in% DIANA_found_lnc]$KD.geneID), 'out/pass_to_DIANA.txt',
            quote = F, row.names = F
)
