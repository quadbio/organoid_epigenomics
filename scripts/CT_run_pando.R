library(tidyverse)
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(doParallel)
library(BSgenome.Hsapiens.UCSC.hg38)

registerDoParallel(24)

rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist
Matrix <- Matrix::Matrix

setwd('~/projects/cutntag/')

data(motifs)
data(motif2tf)

#### Read data ####
marks <- read_rds('data/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')
pc_genes <- read_tsv('~/resources/gene_sets/protein_coding_genes_goodlist.txt', col_names=F)

rna <- FindVariableFeatures(rna, nfeatures=4000)
genes_use <- rna %>% VariableFeatures() %>% intersect(pc_genes$X1)


# rnatac <- read_rds('~/projects/early/data/RNA_ATAC/integration/RNA_ATAC_pseudocells_v2.3pando_srt.rds')
# 
# marks$H3K27ac[['peaks_bin']] <- CreateAssayObject((marks$H3K27ac[['peaks']]@data > 0)*1)
# marks$H3K27me3[['peaks_bin']] <- CreateAssayObject((marks$H3K27me3[['peaks']]@data > 0)*1)
# marks$H3K4me3[['peaks_bin']] <- CreateAssayObject((marks$H3K4me3[['peaks']]@data > 0)*1)
# 
# marks$H3K27ac <- Pando::aggregate_assay(marks$H3K27ac, assay='peaks_bin', group_name='clusters')
# marks$H3K27me3 <- Pando::aggregate_assay(marks$H3K27me3, assay='peaks_bin', group_name='clusters')
# marks$H3K4me3 <- Pando::aggregate_assay(marks$H3K4me3, assay='peaks_bin', group_name='clusters')
# 
# H3K27ac_cluster_detect <- marks$H3K27ac@assays$peaks_bin@misc$summary$clusters
# H3K27me3_cluster_detect <- marks$H3K27me3@assays$peaks_bin@misc$summary$clusters
# H3K4me3_cluster_detect <- marks$H3K4me3@assays$peaks_bin@misc$summary$clusters
# 
# H3K27ac_maxdetect_feats <- colMaxs(H3K27ac_cluster_detect)
# H3K27me3_maxdetect_feats <- colMaxs(H3K27me3_cluster_detect)
# H3K4me3_maxdetect_feats <- colMaxs(H3K4me3_cluster_detect)
# 
# H3K27ac_regions <- colnames(H3K27ac_cluster_detect)[H3K27ac_maxdetect_feats > 0.05]
# H3K27me3_regions <- colnames(H3K27me3_cluster_detect)[H3K27me3_maxdetect_feats > 0.05]
# H3K4me3_regions <- colnames(H3K4me3_cluster_detect)[H3K4me3_maxdetect_feats > 0.05]
# 
# 
# all_regions <- c(H3K27ac_regions, H3K27me3_regions, H3K4me3_regions) %>% 
#     StringToGRanges()
# 
# 
# all_regions <- intersect(all_regions, all_regions)
# 
# rnatac[['peaks_bin']] <- CreateAssayObject((rnatac[['peaks']]@data > 0)*1)
# rnatac <- Pando::aggregate_assay(rnatac, assay='peaks_bin', group_name='highres_clusters')
# atac_cluster_detect <- rnatac@assays$peaks_bin@misc$summary$highres_clusters
# atac_maxdetect_feats <- colMaxs(atac_cluster_detect)
# atac_regions <- colnames(atac_cluster_detect)[atac_maxdetect_feats > 0.1] %>% 
#     StringToGRanges()
# 
# atac_intersects <- intersect(all_regions, atac_regions)
# 
# rnatac_links <- Links(rnatac[['peaks']])
# 
# links_intersects <- intersect(atac_intersects, rnatac_links)
# 
# 
# rnatac_links <- initiate_grn(
#     rnatac, 
#     region=links_intersects,
#     exclude_exons=FALSE
# )
# 
# rnatac_links <- find_motifs(
#     rnatac_links, 
#     genome=BSgenome.Hsapiens.UCSC.hg38,
#     pfm=motifs,
#     motif_tfs=motif2tf
# )



rnatac_links <- read_rds('data/GRN/RNA_ATAC_pseudocells_pando_links_srt.rds')

# rnatac_links <- infer_grn(
#     rnatac_links, 
#     genes=genes_use,
#     peak_to_gene_method='GREAT',
#     upstream=100000,
#     downstream=100000,
#     tf_cor=0.2,
#     aggregate_peaks_col='highres_clusters',
#     method='bagging_ridge',
#     network_name='bagging_ridge_network',
#     parallel=TRUE,
# )
# 
# rnatac_links %>% write_rds('data/GRN/RNA_ATAC_pseudocells_pando_links_srt.rds')

rnatac_links <- infer_grn(
    rnatac_links, 
    genes=genes_use,
    peak_to_gene_method='GREAT',
    upstream=100000,
    downstream=100000,
    tf_cor=0.2,
    aggregate_peaks_col='highres_clusters',
    method='glm',
    network_name='glm_network',
    parallel=TRUE, 
)

rnatac_links %>% write_rds('data/GRN/RNA_ATAC_pseudocells_pando_links_srt.rds')







rnatac_atac <- read_rds('data/GRN/RNA_ATAC_pseudocells_pando_atac_srt.rds')


# rnatac_atac <- initiate_grn(
#     rnatac, 
#     region=atac_intersects,
#     exclude_exons=FALSE
# )
# 
# rnatac_atac <- find_motifs(
#     rnatac_atac, 
#     genome=BSgenome.Hsapiens.UCSC.hg38,
#     pfm=motifs,
#     motif_tfs=motif2tf
# )


# rnatac_atac <- infer_grn(
#     rnatac_atac, 
#     genes=genes_use,
#     peak_to_gene_method='GREAT',
#     upstream=100000,
#     downstream=100000,
#     tf_cor=0.2,
#     aggregate_peaks_col='highres_clusters',
#     method='bagging_ridge',
#     network_name='bagging_ridge_network',
#     parallel=TRUE, 
# )
# 
# rnatac_atac %>% write_rds('data/GRN/RNA_ATAC_pseudocells_pando_atac_srt.rds')

rnatac_atac <- infer_grn(
    rnatac_atac, 
    genes=genes_use,
    peak_to_gene_method='GREAT',
    upstream=100000,
    downstream=100000,
    tf_cor=0.2,
    aggregate_peaks_col='highres_clusters',
    method='glm',
    network_name='glm_network',
    parallel=TRUE, 
)

rnatac_atac %>% write_rds('data/GRN/RNA_ATAC_pseudocells_pando_atac_srt.rds')



