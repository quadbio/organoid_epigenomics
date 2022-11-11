source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')

library(Pando)
library(SeuratWrappers)
library(harmony)

filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist
Matrix <- Matrix::Matrix

setwd('~/projects/cutntag/')

#### Read data ####
drugs_d158 <- read_rds('data/drugs/drugs_d15_d18_v1_v1.2lines_srt.rds')

drugs_d158$inhibitor_target <- case_when(
    str_detect(drugs_d158$inhib_annotation, '485') ~ 'H3K27ac',
    str_detect(drugs_d158$inhib_annotation, '395') ~ 'H3K27me3',
    T ~ 'DMSO'
)

rna <- read_rds('data/RNA/RNA_all_srt_v2.3lines.rds')
rna_pt_meta <- read_tsv('data/RNA/cellrank/RNA_full_cellrank_probs.tsv') %>% 
    dplyr::rename('cell'=1) %>% 
    select(cell, velocity_pseudotime, pseudotime_ranks) %>% 
    column_to_rownames('cell')

rna <- AddMetaData(rna, rna_pt_meta)

drugs_d158_H3K27me3 %>% feature_plot(features='percent.mt', reduction='humap')


#### Get H3K27me3 inhib + DMSO and reintegrate ####
drugs_d158_H3K27me3 <- drugs_d158 %>% subset(inhibitor_target%in%c('H3K27me3', 'DMSO'))

drugs_d158_H3K27me3 <- RunHarmony(drugs_d158_H3K27me3, group.by.vars = 'orig.ident')
drugs_d158_H3K27me3 <- RunUMAP(
    drugs_d158_H3K27me3, 
    reduction='harmony', 
    dims=1:ncol(drugs_d158_H3K27me3[['harmony']]),
    reduction.name='humap'
)
p1 <- dim_plot(drugs_d158_H3K27me3, group.by=c('DMSO', 'inhib_annotation'), reduction='humap')
p2 <- feature_plot(drugs_d158_H3K27me3, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='humap')
p1 / p2 + plot_layout(heights=c(1,3))

ggsave('plots/drugs/drugs_d15_d18_H3K27me3_inhib_harmony_umap.png', width=7, height=10)


drugs_d158_H3K27me3 <- FindNeighbors(drugs_d158_H3K27me3, reduction='harmony')
drugs_d158_H3K27me3 <- FindClusters(drugs_d158_H3K27me3, reduction='harmony', resolution=0.2)

dim_plot(drugs_d158_H3K27me3, group.by=c('seurat_clusters', 'inhib_annotation', 'orig.ident'), reduction='humap', label=T)
ggsave('plots/drugs/drugs_d15_d18_H3K27me3_inhib_clusters_umap.png', width=10, height=5)


dim_plot(drugs_d158_H3K27me3, group.by=c('seurat_clusters'), reduction='humap', label=T, split.by='inhib_annotation')


feature_plot(drugs_d158_H3K27me3, group.by=c('seurat_clusters', 'inhibitor_target', 'orig.ident'), reduction='humap', label=T)


drugs_d158_H3K27me3 %>% write_rds('data/drugs/drugs_d15_d18_A395_v1_v1.2lines_srt.rds')

library(SeuratDisk)
SaveH5Seurat(drugs_d158_H3K27me3, 'data/drugs/drugs_d15_d18_A395_v1_v1.2lines.h5seurat', overwrite = T)
Convert('data/drugs/drugs_d15_d18_A395_v1_v1.2lines.h5seurat', dest='h5ad')





#### Get H3K27ac inhib + DMSO and reintegrate ####
drugs_d158_H3K27ac <- drugs_d158 %>% subset(inhibitor_target%in%c('H3K27ac', 'DMSO'))

drugs_d158_H3K27ac <- RunHarmony(drugs_d158_H3K27ac, group.by.vars = 'orig.ident')
drugs_d158_H3K27ac <- RunUMAP(
    drugs_d158_H3K27ac, 
    reduction='harmony', 
    dims=1:ncol(drugs_d158_H3K27ac[['harmony']]),
    reduction.name='humap'
)
p1 <- dim_plot(drugs_d158_H3K27ac, group.by=c('DMSO', 'inhib_annotation'), reduction='humap')
p2 <- feature_plot(drugs_d158_H3K27ac, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='humap')
p1 / p2 + plot_layout(heights=c(1,3))

ggsave('plots/drugs/drugs_d15_d18_H3K27ac_inhib_harmony_umap.png', width=7, height=10)


drugs_d158_H3K27ac <- FindNeighbors(drugs_d158_H3K27ac, reduction='harmony')
drugs_d158_H3K27ac <- FindClusters(drugs_d158_H3K27ac, reduction='harmony', resolution=0.2)

dim_plot(drugs_d158_H3K27ac, group.by=c('seurat_clusters', 'inhib_annotation', 'orig.ident'), reduction='humap', label=T)
ggsave('plots/drugs/drugs_d15_d18_H3K27ac_inhib_clusters_umap.png', width=10, height=5)


dim_plot(drugs_d158_H3K27ac, group.by=c('seurat_clusters'), reduction='humap', label=T, split.by='inhib_annotation')


drugs_d158_H3K27ac %>% write_rds('data/drugs/drugs_d15_d18_A485_v1_v1.2lines_srt.rds')

library(SeuratDisk)
SaveH5Seurat(drugs_d158_H3K27ac, 'data/drugs/drugs_d15_d18_A485_v1_v1.2lines.h5seurat', overwrite = T)
Convert('data/drugs/drugs_d15_d18_A485_v1_v1.2lines.h5seurat', dest='h5ad')






#### Rank DA peaks ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')

lineage_da <- read_tsv('data/results/diff_expression/all_marks_DA_lineage_coarse.tsv')
pt_clust <- read_tsv('data/trajectories/ctx/all_mod_pseudotime_genes_dtw_clust.tsv')

#### Get detected peaks ####
marks$H3K27ac[['peaks_bin']] <- CreateAssayObject((marks$H3K27ac[['peaks']]@data > 0)*1)
marks$H3K27ac <- Pando::aggregate_assay(marks$H3K27ac, assay='peaks_bin', group_name='clusters')
H3K27ac_peak_detect <- marks$H3K27ac@assays$peaks_bin@misc$summary$clusters
H3K27ac_detected <- colnames(H3K27ac_peak_detect)[colMaxs(H3K27ac_peak_detect)>0.05]

marks$H3K27me3[['peaks_bin']] <- CreateAssayObject((marks$H3K27me3[['peaks']]@data > 0)*1)
marks$H3K27me3 <- Pando::aggregate_assay(marks$H3K27me3, assay='peaks_bin', group_name='clusters')
H3K27me3_peak_detect <- marks$H3K27me3@assays$peaks_bin@misc$summary$clusters
H3K27me3_detected <- colnames(H3K27me3_peak_detect)[colMaxs(H3K27me3_peak_detect)>0.05]

marks$H3K4me3[['peaks_bin']] <- CreateAssayObject((marks$H3K4me3[['peaks']]@data > 0)*1)
marks$H3K4me3 <- Pando::aggregate_assay(marks$H3K4me3, assay='peaks_bin', group_name='clusters')
H3K4me3_peak_detect <- marks$H3K4me3@assays$peaks_bin@misc$summary$clusters
H3K4me3_detected <- colnames(H3K4me3_peak_detect)[colMaxs(H3K4me3_peak_detect)>0.05]


#### Get bigwigs and summarize over peaks ####
bw_path <- '/local2/USERS/jfleck/projects/cutntag/spike_norm/bamcoverage'
bw_names <- list.files(bw_path)
bw_rpkm_files <- list.files(bw_path, recursive=T, pattern='*_rpkm.bw', full.names=T)

names(bw_rpkm_files) <- bw_names

bw_ranges_list <- map(bw_rpkm_files, rtracklayer::import)

bw_H3K27ac_A485_list <- bw_ranges_list[str_detect(names(bw_ranges_list), 'H3K27ac.*A485.*(18|15)')]
bw_H3K27me3_A395_list <- bw_ranges_list[str_detect(names(bw_ranges_list), 'H3K27me3.*A395.*(18|15)')]

bw_H3K27ac_DMSO_list <- bw_ranges_list[str_detect(names(bw_ranges_list), 'H3K27ac.*DMSO.*(18|15)')]
bw_H3K27me3_DMSO_list <- bw_ranges_list[str_detect(names(bw_ranges_list), 'H3K27me3.*DMSO.*(18|15)')]


H3K27me3_A395_peak_score_mat <- bw_H3K27me3_A395_list %>% map(function(x){
    olaps <- findOverlaps(StringToGRanges(H3K27me3_detected), x)
    score_mat <- as.matrix(x$score[subjectHits(olaps)])
    peaks_scores <- Pando::aggregate_matrix(score_mat, groups = queryHits(olaps)) %>% 
        as.numeric()
    names(peaks_scores) <- H3K27me3_detected
    return(peaks_scores)
}) %>% purrr::reduce(rbind) %>% Matrix(sparse=T)

rownames(H3K27me3_A395_peak_score_mat) <- names(bw_H3K27me3_A395_list)


H3K27ac_A485_peak_score_mat <- bw_H3K27ac_A485_list %>% map(function(x){
    olaps <- findOverlaps(StringToGRanges(H3K27ac_detected), x)
    score_mat <- as.matrix(x$score[subjectHits(olaps)])
    peaks_scores <- Pando::aggregate_matrix(score_mat, groups = queryHits(olaps)) %>% 
        as.numeric()
    names(peaks_scores) <- H3K27ac_detected
    return(peaks_scores)
}) %>% purrr::reduce(rbind) %>% Matrix(sparse=T)

rownames(H3K27ac_A485_peak_score_mat) <- names(bw_H3K27ac_A485_list)


H3K27me3_DMSO_peak_score_mat <- bw_H3K27me3_DMSO_list %>% map(function(x){
    olaps <- findOverlaps(StringToGRanges(H3K27me3_detected), x)
    score_mat <- as.matrix(x$score[subjectHits(olaps)])
    peaks_scores <- Pando::aggregate_matrix(score_mat, groups = queryHits(olaps)) %>% 
        as.numeric()
    names(peaks_scores) <- H3K27me3_detected
    return(peaks_scores)
}) %>% purrr::reduce(rbind) %>% Matrix(sparse=T)

rownames(H3K27me3_DMSO_peak_score_mat) <- names(bw_H3K27me3_DMSO_list)


H3K27ac_DMSO_peak_score_mat <- bw_H3K27ac_DMSO_list %>% map(function(x){
    olaps <- findOverlaps(StringToGRanges(H3K27ac_detected), x)
    score_mat <- as.matrix(x$score[subjectHits(olaps)])
    peaks_scores <- Pando::aggregate_matrix(score_mat, groups = queryHits(olaps)) %>% 
        as.numeric()
    names(peaks_scores) <- H3K27ac_detected
    return(peaks_scores)
}) %>% purrr::reduce(rbind) %>% Matrix(sparse=T)

rownames(H3K27ac_DMSO_peak_score_mat) <- names(bw_H3K27ac_DMSO_list)


#### Calculate mean logFC to control ####
H3K27me3_A395_logmean <- (H3K27me3_A395_peak_score_mat+1) %>% log2() %>% colMeans() 
H3K27me3_DMSO_logmean <- (H3K27me3_DMSO_peak_score_mat+1) %>% log2() %>% colMeans() 

H3K27ac_A485_logmean <- (H3K27ac_A485_peak_score_mat+1) %>% log2() %>% colMeans() 
H3K27ac_DMSO_logmean <- (H3K27ac_DMSO_peak_score_mat+1) %>% log2() %>% colMeans() 


H3K27me3_peak_genes <- Signac::ClosestFeature(marks$H3K27ac, StringToGRanges(H3K27me3_detected)) %>% 
    as_tibble()
H3K27ac_peak_genes <- Signac::ClosestFeature(marks$H3K27ac, StringToGRanges(H3K27ac_detected)) %>% 
    as_tibble()


H3K27me3_fc <- H3K27me3_A395_logmean - H3K27me3_DMSO_logmean
hist(H3K27me3_fc)

H3K27ac_fc <- H3K27ac_A485_logmean - H3K27ac_DMSO_logmean
hist(H3K27ac_fc)

H3K27me3_fc_df <- H3K27me3_fc %>% 
    enframe('peak', 'fc') %>% inner_join(H3K27me3_peak_genes, by=c('peak'='query_region'))

H3K27ac_fc_df <- H3K27ac_fc %>% 
    enframe('peak', 'fc') %>% inner_join(H3K27ac_peak_genes, by=c('peak'='query_region'))

drug_ct_fc_df <- bind_rows('H3K27me3'=H3K27me3_fc_df, 'H3K27ac'=H3K27ac_fc_df, .id='mark') 

drug_ct_fc_df %>% write_tsv('data/drugs/drugs_CT_peaks_logfc.tsv')



#### Plot ranked peaks and logFC histograms ####
drug_ct_fc_df <- read_tsv('data/drugs/drugs_CT_peaks_logfc.tsv')

plot_df <- drug_ct_fc_df %>% 
    arrange(desc(fc)) %>% 
    mutate(peak=factor(peak, levels=unique(.$peak)))

p1 <- ggplot(filter(plot_df, mark=='H3K27me3'), aes(peak, fc)) +
    geom_hline(yintercept = 0) +
    geom_point(size=0.5) +
    no_x_text() +
    ggtitle('H3K27me3 peaks / A395') +
    labs(y='log2 fold change')

ph <- ggplot(filter(plot_df, mark=='H3K27me3'), aes(fc)) +
    geom_histogram(color='black', fill='grey') +
    coord_flip() +
    theme_void() 

pm <- p1 + ph + plot_layout(widths=c(3,1))



p1 <- ggplot(filter(plot_df, mark=='H3K27ac'), aes(peak, fc)) +
    geom_hline(yintercept = 0) +
    geom_point(size=0.5) +
    no_x_text() +
    ggtitle('H3K27ac peaks / A485') +
    labs(y='log2 fold change')

ph <- ggplot(filter(plot_df, mark=='H3K27ac'), aes(fc)) +
    geom_histogram(color='black', fill='grey') +
    coord_flip() +
    theme_void() 

pa <- p1 + ph + plot_layout(widths=c(3,1))

pm / pa
ggsave('plots/drugs/CT/CT_peaks_logfc_ranks_hist.png', width=8, height=8)




#### Join with DE results to see which genes/peaks are depleted ####
#### H3K27me3 inhibition ####

H3K27me3_ct_fc <- drug_ct_fc_df %>% filter(mark=='H3K27me3') %>% 
    arrange(desc(fc)) %>% 
    mutate(peak=factor(peak, levels=unique(.$peak)))

A395_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_global_inhib_de.tsv')
A395_nepi_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_nepi_inhib_de.tsv')
A395_cluster_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_cluster_de.tsv')

cl6_de <- A395_cluster_de %>% filter(cluster==6)
top_cl6_de <- A395_cluster_de %>% filter(cluster==6, p_val_adj<1e-4, avg_log2FC>0.5) %>% pull(feature) %>% unique()
top_cl7_de <- A395_cluster_de %>% filter(cluster==7, p_val_adj<1e-4, avg_log2FC>0.5) %>% pull(feature) %>% unique()
top_global_de <- A395_global_de %>% filter(p_val_adj<1e-2, avg_log2FC>0.1) %>% pull(feature) %>% unique()
top_nepi_de <- A395_nepi_de %>% filter(p_val_adj<1e-2, avg_log2FC>0.1) %>% pull(feature) %>% unique()

feats_show <- intersect(top_cl6_de, top_nepi_de) 


p1 <- ggplot(arrange(H3K27me3_ct_fc, gene_name%in%feats_show), aes(peak, fc, color=gene_name%in%feats_show, fill=gene_name%in%feats_show)) +
    geom_hline(yintercept = 0) +
    geom_point(size=0.5) +
    no_x_text() +
    ggtitle('H3K27me3 peaks / A395') +
    labs(y='log2 fold change') +
    no_legend()

p2 <- ggplot(arrange(H3K27me3_ct_fc, gene_name%in%feats_show), aes(gene_name%in%feats_show, fc, fill=gene_name%in%feats_show)) +
    geom_boxplot() +
    theme_void() +
    no_legend()

(p1 | p2) + plot_layout(widths=c(3,1))


plot_df <- inner_join(A395_global_de, H3K27me3_ct_fc, by=c('feature'='gene_name')) %>% 
    mutate(fc_clip=clip_abs(avg_log2FC, 1.5))
    
ggplot(arrange(plot_df, avg_log2FC), aes(peak, '1', fill=fc_clip)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(pals::brewer.rdbu(100)), limits=c(-1.5, 1.5))




#### Upset plot with intersections ####
library(ggupset)

upset_tbl <- tibble(
    feature = A395_cluster_de$feature,
    cluster_6 = A395_cluster_de$feature %in% top_cl6_de,
    cluster_7 = A395_cluster_de$feature %in% top_cl7_de,
    # global = A395_cluster_de$feature %in% top_global_de,
    nepi = A395_cluster_de$feature %in% top_nepi_de
) %>% pivot_longer(!feature, values_to='is_there', names_to='set') %>%
    filter(is_there) %>%
    distinct() %>% 
    select(-is_there) %>% 
    group_by(feature) %>%
    summarize(set = list(set))

ggplot(upset_tbl, aes(x = set)) +
    geom_bar() +
    scale_x_upset()


intersect(top_cl6_de, top_nepi_de)


#### Stratify peaks into bivalent, switching, repressed ####
# Bivalent: H3K4me3 and H3K27me3 togehter in any of psc, nepi, non_nect
# Switching: H3K27me3 and H3K27me3 intersecting peak, exclusive for one in any of psc, nepi, non_nect
# Repressed: H3K27me3-exclusive in all of psc, nepi, non_nect

K27_isect_matches <- read_tsv('data/intersect/H3K27_poised_marks_intersect_matches.tsv')
Kme_isect_matches <- read_tsv('data/intersect/H3Kme_bivalent_marks_intersect_matches.tsv')

Kme_isect_detect <- Kme_isect_matches %>% 
    filter(H3K27me3%in%H3K27me3_detected & H3K4me3%in%H3K4me3_detected)

K27_isect_detect <- K27_isect_matches %>% 
    filter(H3K27me3%in%H3K27me3_detected & H3K27ac%in%H3K27ac_detected)


# Tweek detection in early 
H3K27ac_early_clusters <- marks$H3K27ac@meta.data %>% 
    filter(celltype_jf%in%c('non_nect', 'psc', 'nect') & age%in%c('EB', '15d', '35d')) %>% 
    pull(clusters) %>% unique()
H3K27ac_early_detect <- marks$H3K27ac@assays$peaks_bin@misc$summary$clusters[H3K27ac_early_clusters, ]
H3K27ac_detected_early <- colnames(H3K27ac_early_detect)[colMaxs(H3K27ac_early_detect)>0.05]

H3K27me3_early_clusters <- marks$H3K27me3@meta.data %>% 
    filter(celltype_jf%in%c('non_nect', 'psc', 'nect') & age%in%c('EB', '15d', '35d')) %>% 
    pull(clusters) %>% unique()
H3K27me3_early_detect <- marks$H3K27me3@assays$peaks_bin@misc$summary$clusters[H3K27me3_early_clusters, ]
H3K27me3_detected_early <- colnames(H3K27me3_early_detect)[colMaxs(H3K27me3_early_detect)>0.05]

H3K4me3_early_clusters <- marks$H3K4me3@meta.data %>% 
    filter(celltype_jf%in%c('non_nect', 'psc', 'nect') & age%in%c('EB', '15d', '35d')) %>% 
    pull(clusters) %>% unique()
H3K4me3_early_detect <- marks$H3K4me3@assays$peaks_bin@misc$summary$clusters[H3K4me3_early_clusters, ]
H3K4me3_detected_early <- colnames(H3K4me3_early_detect)[colMaxs(H3K4me3_early_detect)>0.05]


switching_peaks <- K27_isect_detect %>% filter(H3K27me3%in%H3K27me3_detected_early & H3K27ac%in%H3K27ac_detected_early) %>% pull(H3K27me3) %>% unique()
bivalent_peaks <- Kme_isect_detect %>% filter(H3K27me3%in%H3K27me3_detected_early & H3K4me3%in%H3K4me3_detected_early & !H3K27me3%in%switching_peaks) %>% pull(H3K27me3) %>% unique()
repressed_peaks <- H3K27me3_detected %>% intersect(H3K27me3_detected_early) %>% unique() %>% setdiff(bivalent_peaks) %>% setdiff(switching_peaks)
all_peaks <- H3K27me3_detected %>% intersect(H3K27me3_detected_early) %>% unique()

bivalent_peaks_ranges <- bivalent_peaks %>% StringToGRanges() %>% {names(.) <- bivalent_peaks; .}
switching_peaks_ranges <- switching_peaks %>% StringToGRanges() %>% {names(.) <- switching_peaks; .}
repressed_peaks_ranges <- repressed_peaks %>% StringToGRanges() %>% {names(.) <- repressed_peaks; .}
all_peaks_ranges <- all_peaks %>% StringToGRanges() %>% {names(.) <- all_peaks; .}

# bivalent_peaks_ranges %>% export.bed('data/intersect/H3K27me3_early5_bivalent_peaks.bed')
# switching_peaks_ranges %>% export.bed('data/intersect/H3K27me3_early5_switching_peaks.bed')
# repressed_peaks_ranges %>% export.bed('data/intersect/H3K27me3_early5_repressed_peaks.bed')
# all_peaks_ranges %>% export.bed('data/intersect/H3K27me3_early5_all_peaks.bed')

H3K27me3_ct_fc <- H3K27me3_ct_fc %>% 
    mutate(
        bivalent=peak%in%bivalent_peaks,
        switching=peak%in%switching_peaks,
        repressed=peak%in%repressed_peaks
    ) 


H3K27me3_d15_detect <- marks$H3K27me3@assays$peaks_bin@misc$summary$clusters[H3K27me3_early_clusters[str_detect(H3K27me3_early_clusters, 'mid')], ] %>% 
    {colnames(.)[colMaxs(.)>0.05]}


H3K27me3_status_de <- H3K27me3_ct_fc %>% 
    pivot_longer(cols = c(bivalent, switching, repressed), names_to='status') %>% 
    filter(value) %>% select(-value) %>% 
    inner_join(A395_global_de, by=c('gene_name'='feature')) %>% 
    mutate(
        status=factor(status, levels=c('switching', 'bivalent', 'repressed', 'other')),
        d15_detect=peak%in%H3K27me3_d15_detect
    )

p1 <- ggplot(H3K27me3_status_de, aes(avg_log2FC, status, fill=d15_detect)) +
    geom_vline(xintercept = 0) +
    geom_density_ridges(alpha=0.5, scale=0.8) +
    coord_flip() +
    ggtitle('gene logFC')

p2 <- ggplot(H3K27me3_status_de, aes(status, fc, fill=d15_detect)) +
    geom_hline(yintercept = 0) +
    geom_quasirandom(position='dodge', dodge.width=0.8, shape=21) +
    geom_boxplot(alpha=0.5, outlier.shape=NA) +
    ggtitle('peak logFC')

p1 | p2


#### Compare FC between deeptools clusters ####
H3K27me3_395_deeptools_bed <- read_tsv('/local2/USERS/fzenk/drug_treatment/heatmaps/bed/H3K27me3_395_all_JF_classification_all_sort-ac-me_raw.bed')
colnames(H3K27me3_395_deeptools_bed)[1] <- 'chrom'

bed_ranges <- makeGRangesFromDataFrame(H3K27me3_395_deeptools_bed)


peak_order <- subjectHits(findOverlaps(bed_ranges, StringToGRanges(H3K27me3_detected)))
deeptools_clusters <- H3K27me3_395_deeptools_bed %>% select('deeptools_range'=peak, 'cluster'=deepTools_group)
deeptools_order <- H3K27me3_detected[peak_order]
deeptools_clusters$peak <- deeptools_order


H3K27me3_dt_cluster_de <- H3K27me3_ct_fc %>% 
    inner_join(global_de, by=c('gene_name'='feature')) %>% 
    inner_join(deeptools_clusters) %>% 
    mutate(peak=factor(peak, levels=deeptools_order))


p1 <- ggplot(H3K27me3_dt_cluster_de, aes(avg_log2FC, cluster)) +
    geom_vline(xintercept = 0) +
    geom_density_ridges(alpha=0.5, scale=0.8) +
    coord_flip() +
    ggtitle('gene logFC')

p2 <- ggplot(H3K27me3_dt_cluster_de, aes(cluster, fc)) +
    geom_hline(yintercept = 0) +
    geom_quasirandom(position='dodge', dodge.width=0.8, shape=21) +
    geom_boxplot(alpha=0.5, outlier.shape=NA) +
    ggtitle('peak logFC')

p1 | p2


#### Test enrichment of genes in pos FC ####
global_de_sig <- global_de %>% filter(abs(avg_log2FC)>0.1, p_val_adj<0.01)

gene_sets <- H3K27me3_dt_cluster_de %>% 
    # mutate(status=paste0(status, '_', d15_detect)) %>% 
    mutate(status=cluster) %>% 
    group_by(status) %>% group_split() %>% 
    {names(.) <- map_chr(., ~as.character(.x$status)[1]);.} %>% 
    map(~unique(.x$gene_name))

H3K27me3_de_enrich <- map_dfr(gene_sets, function(x){
    test_table <- rbind(
        table(sign(global_de_sig$avg_log2FC)),
        table(sign(filter(global_de_sig, feature%in%x)$avg_log2FC))
    )
    ftest <- fisher.test(test_table[,c('-1','1')])
    enrich_df <- tibble(
        pval = ftest$p.value,
        odds_ratio = ftest$estimate,
        logodds = log2(odds_ratio),
        ngenes_pos = test_table[2,2],
        ngenes_neg = test_table[2,1],
        ngenes = length(x)
    )
}, .id='status')


plot_df <- H3K27me3_de_enrich 

ggplot(plot_df, aes(status, logodds, size=ngenes)) +
    geom_bar(stat='identity', position=position_dodge(width=0.2), size=0, width=0.02, fill='black') +
    geom_point(position=position_dodge(width=0.2), shape=21, color='black', fill='grey') +
    scale_size_continuous(range=c(6,8)) +
    geom_hline(yintercept = 0) +
    labs(y='log2 fold enrichment\n(upregulated genes)', color='detected in day 15')







#### Check where DEG are usually expressed ####
rna[['RNA_bin']] <- CreateAssayObject((rna[['RNA']]@counts > 0) * 1)
rna <- aggregate_assay(rna, group_name='clusters', assay='RNA_bin')

marks$H3K27me3[['cRNA_bin']] <- CreateAssayObject((marks$H3K27me3[['cRNA']]@counts > 0) * 1)
marks$H3K27me3 <- aggregate_assay(marks$H3K27me3, group_name='clusters', assay='cRNA_bin')


A395_global_de_genes <- A395_global_de %>% 
    mutate(padj=p.adjust(p_val, method='fdr')) %>% 
    filter(abs(avg_log2FC)>0.1, padj<1e-4) 

genes_use <- A395_global_de_genes$feature

A395_deg_rna_expr_df <- rna@assays$RNA@misc$summary$clusters[, A395_global_de_genes$feature] %>% 
    colMaxs() %>% {names(.) <- A395_global_de_genes$feature;.} %>% enframe('feature', 'max_expr')
A395_deg_rna_prcexpr_df <- rna@assays$RNA_bin@misc$summary$clusters[, A395_global_de_genes$feature] %>% 
    colMaxs() %>% {names(.) <- A395_global_de_genes$feature;.} %>% enframe('feature', 'max_perc_expr')

plot_df <- A395_global_de_genes %>% 
    inner_join(A395_deg_rna_expr_df) %>% 
    inner_join(A395_deg_rna_prcexpr_df)

ggplot(plot_df, aes(max_perc_expr)) +
    geom_histogram(color='black', fill='grey') +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Max frac expressed in highres clusters', y='Count') 

ggsave('plots/drugs/diff_expression/A395_DEG_max_expr_hist.png', width=6, height=4, bg='white')



#### Pseudotime axis for DEG expression and repression in timecourse ####
rna_marks <- read_rds('data/RNA/all_RNA_marks_combined_clusters_srt.rds')
cluster_pt <- aggregate_matrix(as.matrix(rna$pseudotime_ranks), rna$clusters)[colnames(rna_marks),] %>% as.numeric()

rna_marks$pseudotime_ranks <- cluster_pt

A395_global_up_genes <- A395_global_de %>% 
    mutate(padj=p.adjust(p_val, method='fdr')) %>% 
    filter(avg_log2FC>0.1, padj<1e-4, feature%in%rownames(rna_marks@assays$H3K27me3_RNA)) 

A395_deg_cluster_expr <- rna_marks@assays$RNA[A395_global_up_genes$feature, ]
A395_deg_expr <- rna@assays$RNA@data[A395_global_up_genes$feature, ]

A395_deg_cluster_repr <- rna_marks@assays$H3K27me3_RNA[A395_global_up_genes$feature, ]

max_expr_pt <- apply(A395_deg_expr, 1, function(x){
    wex <- rna$pseudotime_ranks[which.max(x)]
    return(wex)
})

max_cluster_expr_pt <- apply(A395_deg_cluster_expr, 1, function(x){
    wex <- cluster_pt[which.max(x)]
    return(wex)
})

max_cluster_repr_pt <- apply(A395_deg_cluster_repr, 1, function(x){
    wex <- cluster_pt[which.max(x)]
    return(wex)
})

gene_pt <- tibble(
    feature = names(max_cluster_expr_pt),
    expr_max = max_expr_pt[names(max_cluster_expr_pt)],
    cluster_expr_max = max_cluster_expr_pt,
    cluster_repr_max = max_cluster_repr_pt
)

hist(max_expr_pt)
hist(max_cluster_expr_pt)
hist(max_cluster_repr_pt)


p1 <- ggplot(gene_pt, aes(expr_max)) +
    geom_histogram(breaks=seq(0,1,0.1), color='black', fill='grey') +
    ggtitle('Max expression over cells')

p2 <- ggplot(gene_pt, aes(cluster_expr_max)) +
    geom_histogram(breaks=seq(0,1,0.1), color='black', fill='grey') +
    ggtitle('Max expression over clusters')

p3 <- ggplot(gene_pt, aes(cluster_repr_max)) +
    geom_histogram(breaks=seq(0,1,0.1), color='black', fill='grey') +
    ggtitle('Max repression/H3K27me3 over clusters')

p1 / p2 / p3

ggsave('plots/drugs/diff_expression/A395_DEG_pt_expr_hist.png', width=6, height=8, bg='white')



genes_plot <- c('SFRP4', 'VGLL3', 'NRP2', 'S100A11', 'KRT8', 'KRT19', 'FTL', 'DSP', 'STMN2')

rna_marks@active.assay <- 'RNA'
p1 <- feature_plot(rna_marks, features=genes_plot, reduction='umap', pt.size=3) &
    scale_color_gradientn(colors=gyylorrd2())

rna_marks@active.assay <- 'H3K27me3_RNA'
p2 <- feature_plot(rna_marks, features=genes_plot, reduction='umap', pt.size=3, order=T) &
    scale_color_gradientn(colors=gybu2())

p1 / p2
ggsave('plots/drugs/diff_expression/A395_DEG_feature_cluster_umap.png', width=6, height=10, bg='white')


deg_mean_repr <- colMeans(t(scale(t(A395_deg_cluster_repr))))[colnames(rna_marks)]
deg_mean_expr <- colMeans(t(scale(t(A395_deg_cluster_expr))))[colnames(rna_marks)]

rna_marks$A395_deg_mean_repr <- deg_mean_repr
rna_marks$A395_deg_mean_expr <- deg_mean_expr

p1 <- feature_plot(rna_marks, features=c('A395_deg_mean_expr'), reduction='umap', pt.size=4) +
    scale_color_gradientn(colors=gyylorrd2()) +
    ggtitle('Mean expression of A395 DEG')
p2 <- feature_plot(rna_marks, features=c('A395_deg_mean_repr'), reduction='umap', pt.size=4) +
    scale_color_gradientn(colors=gybu2()) +
    ggtitle('Mean repression/H3K27me3 of A395 DEG')

p_umap <- p1 | p2

p1 <- feature_plot(rna_marks, features=c('A395_deg_mean_expr'), reduction='fr', pt.size=4) +
    scale_color_gradientn(colors=gyylorrd2()) +
    ggtitle('')
p2 <- feature_plot(rna_marks, features=c('A395_deg_mean_repr'), reduction='fr', pt.size=4) +
    scale_color_gradientn(colors=gybu2()) +
    ggtitle('')

p_fr <- p1 | p2

p_umap / p_fr
ggsave('plots/drugs/diff_expression/A395_DEG_mean_expr_repr_cluster_umap.png', width=9, height=8, bg='white')






#### H3K27ac inhibition ####
H3K27ac_ct_fc <- drug_ct_fc_df %>% filter(mark=='H3K27ac') %>% 
    arrange(desc(fc)) %>% 
    mutate(peak=factor(peak, levels=unique(.$peak)))

A485_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_global2_inhib_de.tsv')
A485_cluster_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_cluster_inhib_de.tsv')

top_global_de <- A485_global_de %>% filter(p_val_adj<1e-2, avg_log2FC>0.1) %>% pull(feature) %>% unique()

feats_show <- top_global_de

p1 <- ggplot(arrange(H3K27ac_ct_fc, gene_name%in%feats_show), aes(peak, fc, color=gene_name%in%feats_show, fill=gene_name%in%feats_show)) +
    geom_hline(yintercept = 0) +
    geom_point(size=0.5) +
    no_x_text() +
    ggtitle('H3K27me3 peaks / A395') +
    labs(y='log2 fold change') +
    no_legend()

p2 <- ggplot(arrange(H3K27ac_ct_fc, gene_name%in%feats_show), aes(gene_name%in%feats_show, fc, fill=gene_name%in%feats_show)) +
    geom_boxplot() +
    theme_void() +
    no_legend()

(p1 | p2) + plot_layout(widths=c(3,1))


plot_df <- inner_join(A485_global_de, H3K27ac_ct_fc, by=c('feature'='gene_name')) %>% 
    mutate(fc_clip=clip_abs(avg_log2FC, 1.5))

ggplot(arrange(plot_df, avg_log2FC), aes(peak, fc, color=avg_log2FC)) +
    geom_hline(yintercept = 0) +
    geom_point(size=0.5) +
    scale_color_gradientn(colors = rev(pals::brewer.rdbu(100)), limits=c(-4, 4)) +
    no_x_text() +
    ggtitle('H3K27ac peaks / A485') +
    labs(y='log2 fold change') 


ggplot(arrange(plot_df, avg_log2FC), aes(avg_log2FC, fc, color=avg_log2FC)) +
    geom_hline(yintercept = 0) +
    geom_point(size=0.5) +
    scale_color_gradientn(colors = rev(pals::brewer.rdbu(100)), limits=c(-4, 4)) +
    ggtitle('H3K27me3 peaks / A395') +
    labs(y='log2 fold change') 


ggplot(arrange(plot_df, fc_clip), aes(peak, '1', fill=avg_log2FC)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(pals::brewer.rdbu(100)), limits=c(-1.5, 1.5))



#### Stratify peaks into switching, active promoters, enhancers ####

K27_isect_matches <- read_tsv('data/intersect/H3K27_poised_marks_intersect_matches.tsv')
prom_isect_matches <- read_tsv('data/intersect/H3K27K4_promoters_marks_intersect_matches.tsv')

prom_isect_detect <- prom_isect_matches %>% 
    filter(H3K27ac%in%H3K27ac_detected & H3K4me3%in%H3K4me3_detected)

K27_isect_detect <- K27_isect_matches %>% 
    filter(H3K27me3%in%H3K27me3_detected & H3K27ac%in%H3K27ac_detected)


switching_peaks <- K27_isect_detect %>% filter(H3K27me3%in%H3K27me3_detected_early & H3K27ac%in%H3K27ac_detected_early) %>% pull(H3K27me3) %>% unique()
bivalent_peaks <- Kme_isect_detect %>% filter(H3K27me3%in%H3K27me3_detected_early & H3K4me3%in%H3K4me3_detected_early & !H3K27me3%in%switching_peaks) %>% pull(H3K27me3) %>% unique()
repressed_peaks <- H3K27me3_detected %>% intersect(H3K27me3_detected_early) %>% unique() %>% setdiff(bivalent_peaks) %>% setdiff(switching_peaks)
all_peaks <- H3K27me3_detected %>% intersect(H3K27me3_detected_early) %>% unique()


promoter_peaks <- prom_isect_detect %>% filter(H3K4me3%in%H3K4me3_detected_early & H3K27ac%in%H3K27ac_detected_early) %>% pull(H3K27ac) %>% unique()
switching_peaks <- K27_isect_detect %>% filter(H3K27me3%in%H3K27me3_detected_early & H3K27ac%in%H3K27ac_detected_early) %>% pull(H3K27ac) %>% unique()
enhancer_peaks <- H3K27ac_detected %>% unique() %>% setdiff(promoter_peaks) %>% setdiff(switching_peaks)
all_peaks <- H3K27ac_detected %>% unique()


promoter_peaks %>% StringToGRanges() %>% export.bed('data/intersect/H3K27ac_early5_promoter_peaks.bed')
switching_peaks %>% StringToGRanges() %>% export.bed('data/intersect/H3K27ac_early5_switching_peaks.bed')
enhancer_peaks %>% StringToGRanges() %>% export.bed('data/intersect/H3K27ac_early5_enhancer_peaks.bed')
all_peaks %>% StringToGRanges() %>% export.bed('data/intersect/H3K27ac_early5_all_peaks.bed')


H3K27ac_ct_fc <- H3K27ac_ct_fc %>% 
    mutate(
        promoter=peak%in%promoter_peaks,
        switching=peak%in%switching_peaks,
        enhancer=peak%in%enhancer_peaks
    ) 


d15_detect <- marks$H3K27ac@assays$peaks_bin@misc$summary$clusters[H3K27ac_early_clusters[str_detect(H3K27ac_early_clusters, 'mid')], ] %>% 
    {colnames(.)[colMaxs(.)>0.05]}


H3K27ac_status_de <- H3K27ac_ct_fc %>% 
    pivot_longer(cols = c(promoter, switching, enhancer), names_to='status') %>% 
    filter(value) %>% select(-value) %>% 
    inner_join(global_de, by=c('gene_name'='feature')) %>% 
    mutate(
        status=factor(status, levels=c('switching', 'promoter', 'enhancer', 'other')),
        d15_detect=peak%in%d15_detect
    )

p1 <- ggplot(H3K27ac_status_de, aes(avg_log2FC, status, fill=d15_detect)) +
    geom_vline(xintercept = 0) +
    geom_density_ridges(alpha=0.5, scale=0.8) +
    coord_flip() +
    ggtitle('gene logFC')

p2 <- ggplot(H3K27ac_status_de, aes(status, fc, fill=d15_detect)) +
    geom_hline(yintercept = 0) +
    geom_quasirandom(position='dodge', dodge.width=0.8, shape=21) +
    geom_boxplot(alpha=0.5, outlier.shape=NA) +
    ggtitle('peak logFC')

p1 | p2



#### Test enrichment of genes in neg FC ####
global_de_sig <- A485_global_de %>% filter(abs(avg_log2FC)>0.1, p_val_adj<1e-100)

gene_sets <- H3K27ac_status_de %>% 
    group_by(status) %>% group_split() %>% 
    {names(.) <- map_chr(., ~as.character(.x$status)[1]);.} %>% 
    map(~unique(.x$gene_name))

H3K27ac_de_enrich <- map_dfr(gene_sets, function(x){
    test_table <- rbind(
        table(sign(global_de_sig$avg_log2FC)),
        table(sign(filter(global_de_sig, feature%in%x)$avg_log2FC))
    )
    ftest <- fisher.test(test_table[,c('1','-1')])
    enrich_df <- tibble(
        pval = ftest$p.value,
        odds_ratio = ftest$estimate,
        logodds = log2(odds_ratio),
        ngenes_pos = test_table[2,2],
        ngenes_neg = test_table[2,1],
        ngenes = length(x)
    )
}, .id='status')


plot_df <- H3K27ac_de_enrich 
ggplot(plot_df, aes(status, logodds, size=ngenes)) +
    geom_bar(stat='identity', position=position_dodge(width=0.2), size=0, width=0.02, fill='black') +
    geom_point(position=position_dodge(width=0.2), shape=21, color='black', fill='grey') +
    scale_size_continuous(range=c(6,8)) +
    geom_hline(yintercept = 0) +
    labs(y='log2 fold enrichment\n(downregulated genes)', color='detected in day 15')






#### Check where DEG are usually expressed ####
rna[['RNA_bin']] <- CreateAssayObject((rna[['RNA']]@counts > 0) * 1)
rna <- aggregate_assay(rna, group_name='clusters', assay='RNA_bin')

marks$H3K27me3[['cRNA_bin']] <- CreateAssayObject((marks$H3K27me3[['cRNA']]@counts > 0) * 1)
marks$H3K27me3 <- aggregate_assay(marks$H3K27me3, group_name='clusters', assay='cRNA_bin')

A485_global_de_genes <- A485_global_de %>% 
    mutate(padj=p.adjust(p_val, method='fdr')) %>% 
    filter(abs(avg_log2FC)>0.1, padj<1e-4) 

genes_use <- A485_global_de_genes$feature

A485_deg_rna_expr_df <- rna@assays$RNA@misc$summary$clusters[, A485_global_de_genes$feature] %>% 
    colMaxs() %>% {names(.) <- A485_global_de_genes$feature;.} %>% enframe('feature', 'max_expr')
A485_deg_rna_prcexpr_df <- rna@assays$RNA_bin@misc$summary$clusters[, A485_global_de_genes$feature] %>% 
    colMaxs() %>% {names(.) <- A485_global_de_genes$feature;.} %>% enframe('feature', 'max_perc_expr')

plot_df <- A485_global_de_genes %>% 
    inner_join(A485_deg_rna_expr_df) %>% 
    inner_join(A485_deg_rna_prcexpr_df)

ggplot(plot_df, aes(max_perc_expr)) +
    geom_histogram(color='black', fill='grey') +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Max frac expressed in highres clusters', y='Count') 

ggsave('plots/drugs/diff_expression/A485_DEG_max_expr_hist.png', width=6, height=4, bg='white')


#### Pseudotime axis for DEG expression and repression in timecourse ####
rna_marks <- read_rds('data/RNA/all_RNA_marks_combined_clusters_srt.rds')
cluster_pt <- aggregate_matrix(as.matrix(rna$pseudotime_ranks), rna$clusters)[colnames(rna_marks),] %>% as.numeric()

rna_marks$pseudotime_ranks <- cluster_pt

A485_global_down_genes <- A485_global_de %>% 
    mutate(padj=p.adjust(p_val, method='fdr')) %>% 
    filter(avg_log2FC<(-0.1), padj<1e-4, feature%in%rownames(rna_marks@assays$H3K27me3_RNA)) 

A485_deg_cluster_expr <- rna_marks@assays$RNA[A485_global_up_genes$feature, ]
A485_deg_expr <- rna@assays$RNA@data[A485_global_up_genes$feature, ]

A485_deg_cluster_act <- rna_marks@assays$H3K27ac_RNA[A485_global_down_genes$feature, ]

max_expr_pt <- apply(A485_deg_expr, 1, function(x){
    wex <- rna$pseudotime_ranks[which.max(x)]
    return(wex)
})

max_cluster_expr_pt <- apply(A485_deg_cluster_expr, 1, function(x){
    wex <- cluster_pt[which.max(x)]
    return(wex)
})

max_cluster_act_pt <- apply(A485_deg_cluster_act, 1, function(x){
    wex <- cluster_pt[which.max(x)]
    return(wex)
})

gene_pt <- tibble(
    feature = names(max_cluster_expr_pt),
    expr_max = max_expr_pt[names(max_cluster_expr_pt)],
    cluster_expr_max = max_cluster_expr_pt,
    cluster_act_max = max_cluster_act_pt
)

hist(max_expr_pt)
hist(max_cluster_expr_pt)
hist(max_cluster_act_pt)


p1 <- ggplot(gene_pt, aes(expr_max)) +
    geom_histogram(breaks=seq(0,1,0.1), color='black', fill='grey') +
        ggtitle('Max expression over cells')

p2 <- ggplot(gene_pt, aes(cluster_expr_max)) +
    geom_histogram(breaks=seq(0,1,0.1), color='black', fill='grey') +
        ggtitle('Max expression over clusters')

p3 <- ggplot(gene_pt, aes(cluster_act_max)) +
    geom_histogram(breaks=seq(0,1,0.1), color='black', fill='grey') +
    ggtitle('Max activation/H3K27ac over clusters')

p1 / p2 / p3

ggsave('plots/drugs/diff_expression/A485_DEG_pt_expr_hist.png', width=6, height=8, bg='white')


genes_plot <- c('LDHA', 'SFRP1', 'SFRP2', 'NR2F1', 'PTN', 'PAX6', 'HES4', 'BNIP3', 'HES5')

rna_marks@active.assay <- 'RNA'
p1 <- feature_plot(rna_marks, features=genes_plot, reduction='umap', pt.size=3) &
    scale_color_gradientn(colors=gyylorrd2())

rna_marks@active.assay <- 'H3K27ac_RNA'
p2 <- feature_plot(rna_marks, features=genes_plot, reduction='umap', pt.size=3, order=T) &
    scale_color_gradientn(colors=gyylgn2())

p1 / p2
ggsave('plots/drugs/diff_expression/A485_DEG_feature_cluster_umap.png', width=6, height=10, bg='white')



deg_mean_act <- colMeans(t(scale(t(A485_deg_cluster_act))))[colnames(rna_marks)]
deg_mean_expr <- colMeans(t(scale(t(A485_deg_cluster_expr))))[colnames(rna_marks)]

rna_marks$A485_deg_mean_act <- deg_mean_act
rna_marks$A485_deg_mean_expr <- deg_mean_expr

p1 <- feature_plot(rna_marks, features=c('A485_deg_mean_expr'), reduction='umap', pt.size=4) +
    scale_color_gradientn(colors=gyylorrd2()) +
    ggtitle('Mean expression of A485 DEG')
p2 <- feature_plot(rna_marks, features=c('A485_deg_mean_act'), reduction='umap', pt.size=4) +
    scale_color_gradientn(colors=gyylgn2()) +
    ggtitle('Mean activation/H3K27ac of A485 DEG')

p_umap <- p1 | p2

p1 <- feature_plot(rna_marks, features=c('A485_deg_mean_expr'), reduction='fr', pt.size=4) +
    scale_color_gradientn(colors=gyylorrd2()) +
    ggtitle('')
p2 <- feature_plot(rna_marks, features=c('A485_deg_mean_act'), reduction='fr', pt.size=4) +
    scale_color_gradientn(colors=gyylgn2()) +
    ggtitle('')

p_fr <- p1 | p2

p_umap / p_fr
ggsave('plots/drugs/diff_expression/A485_DEG_mean_expr_repr_cluster_umap.png', width=9, height=8, bg='white')

