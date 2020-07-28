library("biomaRt")
library(data.table)
library(stringr)
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



# get ONLY needed RNAs from annotation
ann <- fread('/home/mazurovev/fantom6/hg19.cage_peak_phase1and2combined_ann.txt', header = T)
la <- str_split_fixed(ann$short_description, '\\@', 2)
ann$name <- unlist(str_split_fixed(ann$short_description, '\\@', 2)[,2])
ensemble.gene <- getBM(attributes = c('ensembl_gene_id', "hgnc_symbol"), 
                      filters = 'hgnc_symbol',
                      values = unique(ann$name)[!unique(ann$name) %in% unique(ann$name)[grep('chr', unique(ann$name))]], 
                      mart = ensembl)


ann.lnc <- ann[name %in% unique(fantom_DE$KD.geneSymbol)]



# get entrez ids for all lncrnas
lnc.entrez <- getBM(attributes = c('entrezgene_id', "hgnc_symbol", "ensembl_gene_id", "hgnc_id"), 
                    filters = 'ensembl_gene_id',
                    values = unique(fantom_DE$KD.geneID), 
                    mart = ensembl)








# get ids to retrieve ASOs to extract further fastq of the transcripts
fantom_DE <- fread('../oligo_DE_Summary_gene_filtered.tsv')
fantom_DE$geneSymbol <- stringr::str_remove(fantom_DE$geneSymbol, 'HG')
fantom_DE$perturb_id <- stringr::str_remove(fantom_DE$perturb_id, 'ASO_')
fantom_DE$perturb_id_prefix <- unlist(str_split_fixed(fantom_DE$perturb_id, '\\_', 2)[,1])
fantom_DE$KD.geneID_ASO_fa_format = paste(fantom_DE$KD.geneID, fantom_DE$perturb_id_prefix, sep = '|')
unique(fantom_DE$KD.geneID_ASO)
write.table(unique(paste(fantom_DE$KD.geneID, fantom_DE$perturb_id_prefix, sep = '|')),
            file = 'out/KD.gene.ID|ASO_id.txt', 
            quote = F,
            row.names = F,
            col.names = F)



# check names
lncRNAs$gene.name <- fantom_DE[match(lncRNAs$ensembl_gene_id, KD.geneID)]$KD.geneSymbol
