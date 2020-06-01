library("biomaRt")
ensembl = useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

attributes = listAttributes(ensembl)
filters = listFilters(ensembl)

# degs.biotypes <- sapply(degs, function(dt.per.hour){
#   ids <- dt.per.hour$Symbol
#   getBM(attributes = c("mgi_symbol", "gene_biotype"), 
#         filters = "mgi_symbol",
#         values = ids, 
#         mart = ensembl)
# },USE.NAMES = T, simplify = F)
# 
# DE.lncRNAs <- unique(unlist(sapply(degs.biotypes, function(dt.per.hour){
#   as.data.table(dt.per.hour)[gene_biotype == 'lncRNA']$mgi_symbol
# },USE.NAMES = T, simplify = F)))

lncRNAs <- getBM(attributes = c('ensembl_gene_id', "external_transcript_name"), 
                    filters = 'ensembl_gene_id',
                    values = unique(fantom_DE$KD.geneID), 
                    mart = ensembl)
sequences <- getSequence(id = lncRNAs$external_transcript_name, 
                         type = "external_transcript_name", 
                         seqType = "transcript_exon_intron", 
                         mart = ensembl)
sequences <- sequences[order(sequences$external_transcript_name),]
exportFASTA(sequences, file = 'data/lncRNAs.fasta')

hist(str_count(sequences$transcript_exon_intron), breaks = 100)

# transcripts.pos <- getBM(attributes = c('external_transcript_name', 'chromosome_name',
#                                         'transcript_start', 'transcript_end'),
#                          filters = 'external_transcript_name', 
#                          values = unique(DE.lncRNAs$external_transcript_name), 
#                          mart = ensembl)
# colnames(transcripts.pos) <- c('transcript', 'chr', 'start', 'end')
# transcripts.pos$chr <- paste0('chr', transcripts.pos$chr)

# library(liftOver)
# ch = import.chain('data/mm10ToMm9.over.chain')
# gr <- GRanges(transcripts.pos)
# genome(gr) = "mm10"
# transcripts.pos.liftover <- as.data.table(as.data.frame(liftOver(gr, ch)))
# transcripts.pos.liftover <- merge(x = transcripts.pos.liftover,
#                                   y = DE.lncRNAs,
#                                   by.x = 'transcript',
#                                   by.y = 'external_transcript_name')
# transcripts.pos.liftover[,c('group', 'group_name'):=NULL]
