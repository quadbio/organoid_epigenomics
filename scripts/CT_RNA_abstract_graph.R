source('~/scripts/single_cell/graphs.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')

library(Pando)

setwd('~/projects/cutntag/')

select <- dplyr::select

#### Read data ###
marks <- read_rds('data/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')

cluster_corr_df <- read_tsv('data/RNA/integration/cluster_corr_df.tsv')

cluster_meta <- as_tibble(rna@meta.data, rownames='cell') %>%
    inner_join(as_tibble(rna[['cssumap']]@cell.embeddings, rownames='cell')) %>%
    group_by(clusters, lineage, state, celltype_jf) %>%
    summarize(CSSUMAP_1=mean(CSSUMAP_1), CSSUMAP_2=mean(CSSUMAP_2)) %>%
    filter(celltype_jf!='other')


#### Get matches ####
H3K4me3_corr_matches <- filter(cluster_corr_df, mark=='H3K4me3') %>%
    mutate(mark_stage=str_replace(H3K4me3_clusters, '(\\w+)_\\d+', '\\1'), rna_stage=str_replace(rna_clusters, '(\\w+)_\\d+', '\\1')) %>%
    filter(mark_stage==rna_stage) %>%
    group_by(rna_clusters, mark_stage) %>%
    filter(corr==max(corr)) %>% ungroup() %>% select(H3K4me3_clusters, rna_clusters)

H3K27me3_corr_matches <- filter(cluster_corr_df, mark=='H3K27me3') %>%
    mutate(mark_stage=str_replace(H3K27me3_clusters, '(\\w+)_\\d+', '\\1'), rna_stage=str_replace(rna_clusters, '(\\w+)_\\d+', '\\1')) %>%
    filter(mark_stage==rna_stage) %>%
    group_by(rna_clusters, mark_stage) %>%
    filter(corr==min(corr)) %>% ungroup() %>% select(H3K27me3_clusters, rna_clusters)

H3K27ac_corr_matches <- filter(cluster_corr_df, mark=='H3K27ac') %>%
    mutate(mark_stage=str_replace(H3K27ac_clusters, '(\\w+)_\\d+', '\\1'), rna_stage=str_replace(rna_clusters, '(\\w+)_\\d+', '\\1')) %>%
    filter(mark_stage==rna_stage) %>%
    group_by(rna_clusters, mark_stage) %>%
    filter(corr==max(corr)) %>% ungroup() %>% select(H3K27ac_clusters, rna_clusters)



#### Put together with embedding coords ####
clusters_matched <- cluster_meta %>%
    left_join(H3K4me3_corr_matches, by=c('clusters'='rna_clusters')) %>%
    left_join(H3K27me3_corr_matches, by=c('clusters'='rna_clusters')) %>%
    left_join(H3K27ac_corr_matches, by=c('clusters'='rna_clusters'))


ggplot(clusters_matched, aes(CSSUMAP_1, CSSUMAP_2, color=lineage)) +
    geom_point()



#### Combine in object ####
marks$H3K27ac <- Pando::aggregate_assay(marks$H3K27ac, assay = 'cRNA', group_name = 'clusters')
marks$H3K27me3 <- Pando::aggregate_assay(marks$H3K27me3, assay = 'cRNA', group_name = 'clusters')
marks$H3K4me3 <- Pando::aggregate_assay(marks$H3K4me3, assay = 'cRNA', group_name = 'clusters')

H3K27ac_RNA <- marks$H3K27ac@assays$cRNA@misc$summary$clusters[clusters_matched$H3K27ac_clusters, ]
rownames(H3K27ac_RNA) <- clusters_matched$clusters
H3K27me3_RNA <- marks$H3K27me3@assays$cRNA@misc$summary$clusters[clusters_matched$H3K27me3_clusters, ]
rownames(H3K27me3_RNA) <- clusters_matched$clusters
H3K4me3_RNA <- marks$H3K4me3@assays$cRNA@misc$summary$clusters[clusters_matched$H3K4me3_clusters, ]
rownames(H3K4me3_RNA) <- clusters_matched$clusters

H3K27ac_peaks <- marks$H3K27ac@assays$peaks@misc$summary$clusters[clusters_matched$H3K27ac_clusters, ]
rownames(H3K27ac_peaks) <- clusters_matched$clusters
H3K27me3_peaks <- marks$H3K27me3@assays$peaks@misc$summary$clusters[clusters_matched$H3K27me3_clusters, ]
rownames(H3K27me3_peaks) <- clusters_matched$clusters
H3K4me3_peaks <- marks$H3K4me3@assays$peaks@misc$summary$clusters[clusters_matched$H3K4me3_clusters, ]
rownames(H3K4me3_peaks) <- clusters_matched$clusters

clust_RNA <- rna@assays$RNA@misc$summary$clusters[clusters_matched$clusters, ]

cluster_srt <- CreateSeuratObject(t(clust_RNA), meta.data = column_to_rownames(clusters_matched, 'clusters'))
cluster_srt[['H3K27ac_RNA']] <- CreateAssayObject(t(H3K27ac_RNA))
cluster_srt[['H3K27me3_RNA']] <- CreateAssayObject(t(H3K27me3_RNA))
cluster_srt[['H3K4me3_RNA']] <- CreateAssayObject(t(H3K4me3_RNA))
cluster_srt[['H3K27ac_peaks']] <- CreateAssayObject(t(H3K27ac_peaks))
cluster_srt[['H3K27me3_peaks']] <- CreateAssayObject(t(H3K27me3_peaks))
cluster_srt[['H3K4me3_peaks']] <- CreateAssayObject(t(H3K4me3_peaks))

umap_mat <- as.matrix(column_to_rownames(clusters_matched, 'clusters')[,c('CSSUMAP_1', 'CSSUMAP_2')])
colnames(umap_mat) <- NULL
cluster_srt[['umap']] <- CreateDimReducObject(umap_mat, key='UMAP_')

genes_plot <- unique(c('C1orf61', 'NFIX', 'NFIB', 'EMX1', 'POU3F1', 'FEZF2', 'LHX2', 'ZIC5', 'ZIC2', 'ZIC3', 'NEUROD2', 'THTPA', 'NHLH1', 'FOXG1', 'KCNH2', 'UNC5A', 'FEZF2', 'SEZ6', 'SLA', 'OTX2', 'EBF3', 'POU5F1', 'DLX2'))

cluster_srt@active.assay <- 'H3K27ac_RNA'
p1 <- feature_plot(cluster_srt, features=genes_plot) &
    scale_color_gradientn(colors=pals::brewer.blues(100))

cluster_srt@active.assay <- 'RNA'
p2 <- feature_plot(cluster_srt, features=genes_plot)

p1 | p2

cluster_srt %>% write_rds('data/all_RNA_marks_combined_clusters_srt.rds')



#### Add FR embedding for neuronal only ####
cluster_srt <- read_rds('data/all_RNA_marks_combined_clusters_srt.rds')
cluster_meta <- read_tsv('data/noastro_RNA_cluster_meta.tsv')
cluster_fr <- cluster_meta %>%
    select(clusters, FR1, FR2) %>%
    column_to_rownames('clusters')

noastro_cluster_srt <- subset(cluster_srt, cells=cluster_meta$clusters)

noastro_cluster_srt <- AddMetaData(noastro_cluster_srt, cluster_fr)
noastro_cluster_srt$lineage <- case_when(
    noastro_cluster_srt$celltype_jf %in% c('mesen_ex', 'rhom_ex', 'nt_npc') ~ 'mesen_rhom',
    noastro_cluster_srt$celltype_jf %in% c('ctx_npc', 'ctx_ip', 'ctx_ex') ~ 'ctx',
    noastro_cluster_srt$celltype_jf %in% c('RPC', 'RGC') ~ 'retina',
    noastro_cluster_srt$celltype_jf %in% c('dien_npc', 'dien_ex') ~ 'dien',
    noastro_cluster_srt$celltype_jf %in% c('nect') ~ 'nect',
    noastro_cluster_srt$celltype_jf %in% c('psc') ~ 'EB',
    noastro_cluster_srt$celltype_jf %in% c('astrocytes') ~ 'astro',
    noastro_cluster_srt$celltype_jf %in% c('non_nect') ~ 'nn',
    noastro_cluster_srt$celltype_jf %in% c('choroid_plexus') ~ 'chp',
    noastro_cluster_srt$celltype_jf %in% c('OPC') ~ 'OPC'
)

noastro_cluster_srt[['fr']] <- CreateDimReducObject(as.matrix(cluster_fr), key = 'FR_')

noastro_cluster_srt@active.assay <- 'H3K27ac_RNA'

dim_plot(noastro_cluster_srt, group.by='celltype_jf', pt.size=4, reduction='fr') +
    scale_color_manual(values=pantone_celltype)

feature_plot(
    noastro_cluster_srt,
    reduction='fr',
    features=c('GLI3', 'VSX2', 'LHX5', 'RSPO3', 'POU5F1', 'NFIB'),
    pt.size=1, order=T
)

noastro_cluster_srt@active.assay <- 'H3K27me3_RNA'

feature_plot(
    noastro_cluster_srt,
    reduction='fr',
    features=c('GLI3', 'VSX2', 'LHX5', 'RSPO3', 'POU5F1', 'NFIB'),
    pt.size=1, order=T
)


noastro_cluster_srt %>% write_rds('data/noastro_RNA_marks_combined_clusters_srt.rds')




