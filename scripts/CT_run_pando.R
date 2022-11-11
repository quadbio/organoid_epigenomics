source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(doParallel)
library(BSgenome.Hsapiens.UCSC.hg38)

data('motif2tf')
data('motifs')
registerDoParallel(24)

rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist

setwd('~/projects/cutntag/')


#### Read data ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')
rnatac <- read_rds('/local2/USERS/jfleck/projects/early/data/RNA_ATAC_pseudocells_v2.2links_srt.rds')
grn_features <- read_rds('~/projects/early/data/gene_sets/RNA_feature_sets.rds')$grn_features

all_ct_ranges <- map(marks, function(srt){
        srt[['peaks_bin']] <- CreateAssayObject((srt[['peaks']]@data > 0)*1)
        srt <- Pando::aggregate_assay(srt, assay='peaks_bin', group_name='clusters')
        clust_summary <- srt@assays$peaks_bin@misc$summary$clusters
        peak_detection <- colMaxs(clust_summary)
        det_peaks <- colnames(clust_summary)[peak_detection>0.05]
        return(StringToGRanges(det_peaks))
    }) %>% 
    purrr::reduce(IRanges::union)


#### Run Pando with CT peaks + links ####
links_ranges <- Links(rnatac[['peaks']])$peak %>% StringToGRanges()
ranges_use <- IRanges::intersect(links_ranges, all_ct_ranges)

rnatac_pando <- initiate_grn(
    rnatac,
    regions = ranges_use,
    peak_assay = 'peaks',
    rna_assay = 'RNA',
    exclude_exons = F
)

rnatac_pando <- find_motifs(
    rnatac_pando,
    pfm = motifs,
    motif_tfs = motif2tf,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

rnatac_pando <- find_motifs(
    rnatac_pando,
    pfm = motifs,
    motif_tfs = motif2tf,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

rnatac_pando <- infer_grn(
    rnatac_pando,
    genes = grn_features,
    peak_to_gene_method = 'GREAT',
    upstream = 1e+5,
    downstream = 1e+5,
    aggregate_peaks_col = 'highres_clusters',
    method = 'glm',
    parallel = T
)

links_grn <- GetGRN(rnatac_pando)
links_grn %>% write_rds('data/grn/pando_CT_links_grn.rds')


#### Run Pando with CT peaks + 20% max detection ####
rnatac[['peaks_bin']] <- CreateAssayObject((rnatac[['peaks']]@data > 0)*1)
rnatac <- Pando::aggregate_assay(rnatac, assay='peaks_bin', group_name='highres_clusters')
clust_summary <- rnatac@assays$peaks_bin@misc$summary$highres_clusters
peak_detection <- colMaxs(clust_summary)
det_peaks <- colnames(clust_summary)[peak_detection>0.2]
peaks_ranges <- StringToGRanges(det_peaks)

ranges_use <- IRanges::intersect(peaks_ranges, all_ct_ranges)

rnatac_pando <- initiate_grn(
    rnatac, 
    regions = ranges_use,
    peak_assay = 'peaks', 
    rna_assay = 'RNA',
    exclude_exons = F
)

rnatac_pando <- find_motifs(
    rnatac_pando, 
    pfm = motifs, 
    motif_tfs = motif2tf, 
    genome = BSgenome.Hsapiens.UCSC.hg38
)

rnatac_pando <- find_motifs(
    rnatac_pando, 
    pfm = motifs, 
    motif_tfs = motif2tf, 
    genome = BSgenome.Hsapiens.UCSC.hg38
)

rnatac_pando <- infer_grn(
    rnatac_pando, 
    genes = grn_features, 
    peak_to_gene_method = 'GREAT',
    upstream = 1e+5,
    downstream = 1e+5,
    aggregate_peaks_col = 'highres_clusters',
    method = 'glm',
    parallel = T
)

pando_grn <- GetGRN(rnatac_pando)
pando_grn %>% write_rds('data/grn/pando_CT_det20_grn.rds')



#### Run Pando with CT peaks + 10% max detection ####
det_peaks <- colnames(clust_summary)[peak_detection>0.1]
peaks_ranges <- StringToGRanges(det_peaks)

ranges_use <- IRanges::intersect(peaks_ranges, all_ct_ranges)

rnatac_pando <- initiate_grn(
    rnatac, 
    regions = ranges_use,
    peak_assay = 'peaks', 
    rna_assay = 'RNA',
    exclude_exons = F
)

rnatac_pando <- find_motifs(
    rnatac_pando, 
    pfm = motifs, 
    motif_tfs = motif2tf, 
    genome = BSgenome.Hsapiens.UCSC.hg38
)

rnatac_pando <- find_motifs(
    rnatac_pando, 
    pfm = motifs, 
    motif_tfs = motif2tf, 
    genome = BSgenome.Hsapiens.UCSC.hg38
)

rnatac_pando <- infer_grn(
    rnatac_pando, 
    genes = grn_features, 
    peak_to_gene_method = 'GREAT',
    upstream = 1e+5,
    downstream = 1e+5,
    aggregate_peaks_col = 'highres_clusters',
    method = 'glm',
    parallel = T
)

pando_grn <- GetGRN(rnatac_pando)
pando_grn %>% write_rds('data/grn/pando_CT_det10_grn.rds')







