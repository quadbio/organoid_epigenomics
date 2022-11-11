source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')

library(Pando)
library(doParallel)

registerDoParallel(36)

filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist
Matrix <- Matrix::Matrix

setwd('~/projects/cutntag/')

#### Read stuff ####
rna <- read_rds('data/RNA/RNA_all_srt_v2.3lines.rds')
rna_pt_meta <- read_tsv('data/RNA/cellrank/RNA_full_cellrank_probs.tsv') %>% 
    dplyr::rename('cell'=1) %>% 
    select(cell, velocity_pseudotime, pseudotime_ranks) %>% 
    column_to_rownames('cell')

rna <- AddMetaData(rna, rna_pt_meta)

##### Get pseudotimes for all genes ####
rna_expr <- rna@assays$RNA@data
rna_cluster_expr <- rna@assays$RNA@misc$summary$clusters

weighted_mean_pt <- map_par(seq(nrow(rna_expr)), function(i){
    x <- rna_expr[i,]
    non0 <- which(x!=0)
    wex <- weighted.mean(rna$pseudotime_ranks[non0], x[non0])
    return(wex)
}, parallel=T) %>% unlist()

max_pt <- map_par(seq(nrow(rna_expr)), function(i){
    x <- rna_expr[i,]
    wex <- rna$pseudotime_ranks[which.max(x)]
    return(wex)
}, parallel=T) %>% unlist()

cluster_pt <- aggregate_matrix(as.matrix(rna$pseudotime_ranks), rna$clusters) %>% as.numeric()
max_cluster_pt <- map_par(seq(ncol(rna_cluster_expr)), function(i){
    x <- rna_cluster_expr[,i]
    wex <- cluster_pt[which.max(x)]
    return(wex)
}, parallel=T) %>% unlist()

weighted_sd_pt <- map_par(seq(nrow(rna_expr)), function(i){
    x <- rna_expr[i,]
    non0 <- which(x!=0)
    wvar <- Hmisc::wtd.var(rna$pseudotime_ranks[non0], x[non0])
    return(wvar)
}, parallel=T) %>% unlist()

weighted_sd_pt <- sqrt(weighted_sd_pt)


gene_pt <- tibble(
    gene = rownames(rna_expr),
    mean = weighted_mean_pt,
    sd = weighted_sd_pt,
    max = max_pt,
    cluster_max = max_cluster_pt
)

gene_pt %>% write_tsv('data/RNA/RNA_gene_pseudotime.tsv')




