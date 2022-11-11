source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(doParallel)
registerDoParallel(4)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')

#### Read data ####
marks <- read_rds('data/all_marks_list_v3.4lines.rds')
gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')
links_df <- read_tsv('~/projects/early/data/RNA_ATAC/RNA_ATAC_linked_peaks_1mb.tsv')

gene_annot_use <- gene_annot[gene_annot$gene_biotype=='protein_coding']

#### Get gene activities for gene bodies ####
gene_annot_gb <- CollapseToLongestTranscript(gene_annot_use)
gene_act_gb <- map_par(marks, function(srt){
    gene_act <- FeatureMatrix(
        Fragments(srt), features=gene_annot_gb, process_n=10000, cells=colnames(srt)
    )
    ranges_names <- gene_annot_gb$gene_name[match(rownames(gene_act), GRangesToString(gene_annot_gb))]
    rownames(gene_act) <- ranges_names
    return(gene_act)
}, parallel=T)

gene_act_gb %>% write_rds('data/CT_gene_body_activity.rds')


#### Get gene activities for promoters ####
gene_annot_prom <- promoters(gene_annot_gb, upstream=5000, downstream=0)

gene_act_prom <- map_par(marks, function(srt){
    gene_act <- FeatureMatrix(
        Fragments(srt), features=gene_annot_prom, process_n=10000, cells=colnames(srt)
    )
    ranges_names <- gene_annot_prom$gene_name[match(rownames(gene_act), GRangesToString(gene_annot_prom))]
    gene_act <- Pando::aggregate_matrix(gene_act, ranges_names)
    return(gene_act)
}, parallel=T)

gene_act_prom %>% write_rds('data/CT_promoter_activity.rds')


#### Get gene activities for distal CREs / enhancers ####
gene_annot_prom <- promoters(gene_annot_gb, upstream=2000, downstream=2000)

links_ranges <- links_df$peak %>% StringToGRanges()
links_ranges$gene_name <- links_df$gene
prom_olaps <- queryHits(findOverlaps(links_ranges, gene_annot_prom))
gb_olaps <- queryHits(findOverlaps(links_ranges, gene_annot_gb))

links_ranges <- links_ranges[-c(prom_olaps, gb_olaps)]

gene_act_prom <- map_par(marks, function(srt){
    gene_act <- FeatureMatrix(
        Fragments(srt), features=links_ranges, process_n=10000, cells=colnames(srt)
    )
    ranges_names <- links_ranges$gene_name[match(rownames(gene_act), GRangesToString(links_ranges))]
    gene_act <- Pando::aggregate_matrix(gene_act, ranges_names)
    return(gene_act)
}, parallel=T)

gene_act_prom %>% write_rds('data/CT_distal_activity.rds')








