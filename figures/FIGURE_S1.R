source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist

setwd('~/projects/cutntag/')


#### Read data ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')


#### Matching heatmaps ####
cluster_corr_list <- read_rds('data/RNA/integration/cluster_correlation_list.rds')
H3K27ac_final <- read_tsv('data/RNA/integration/cluster_matches_H3K27ac.tsv')
H3K27me3_final <- read_tsv('data/RNA/integration/cluster_matches_H3K27me3.tsv')
H3K4me3_final <- read_tsv('data/RNA/integration/cluster_matches_H3K4me3.tsv')

cluster_meta <- rna@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    distinct(clusters, stage, celltype_jf)

#### Plot heatmaps ####
ct_order <- c('psc', 'non_nect', 'nect', 'ctx_npc', 'ctx_ip', 'ctx_ex', 'dien_npc', 'dien_ex', 
              'nt_npc', 'mesen_ex', 'rhom_ex', 'RPC', 'RGC', 'astrocytes', 'OPC', 'choroid_plexus', 'other')

stage_order <- c('EB', 'mid', 'late', 'mo8', 'retina')

rna_H3K27ac_cor <- cluster_corr_list$H3K27ac
H3K27ac_clust <- rna_H3K27ac_cor %>% t() %>% dist() %>% hclust()
H3K27ac_order <- H3K27ac_clust$labels[H3K27ac_clust$order]
rna_clust <- rna_H3K27ac_cor %>% dist() %>% hclust()
rna_order <- rna_clust$labels[rna_clust$order]

rna_H3K27ac_cor_df <- rna_H3K27ac_cor %>% 
    as_tibble(rownames='rna_clusters') %>% 
    pivot_longer(!rna_clusters, names_to='H3K27ac_clusters', values_to='corr') %>% 
    inner_join(cluster_meta, by=c('rna_clusters'='clusters')) %>% 
    left_join(select(H3K27ac_final, rna_clusters, mark_clusters, corr), by='rna_clusters') %>% 
    mutate(
        H3K27ac_clusters = factor(H3K27ac_clusters, levels=H3K27ac_order),
        rna_clusters = factor(rna_clusters, levels=unique(.$rna_clusters)),
        H3K27ac_stage = factor(str_replace(H3K27ac_clusters, '(\\w+)_\\d+', '\\1'), levels=stage_order),
        rna_stage = factor(str_replace(rna_clusters, '(\\w+)_\\d+', '\\1'), levels=stage_order),
        celltype_jf=factor(celltype_jf, levels=ct_order)
    ) 


ggplot(rna_H3K27ac_cor_df, aes(rna_clusters, H3K27ac_clusters, fill=corr.x)) +
    geom_tile() +
    geom_point(data=filter(rna_H3K27ac_cor_df, corr.x==corr.y), size=0.1, shape=4, stroke=0.2) +
    scale_fill_gradientn(limits=c(-0.65, 0.65), colors=rev(pals::brewer.spectral(100))) +
    facet_grid(H3K27ac_stage~celltype_jf, scales = 'free', space='free') +
    article_text() +
    no_x_text() + no_y_text() +
    theme(
        panel.border = element_blank(),
        panel.spacing = unit(0.02, 'cm'),
        plot.margin = unit(rep(0, 4), 'cm')
    ) 
ggsave('plots/paper/sfig1/H3K27ac_matches_heatmap.pdf', width=8, height=10, units='cm')



rna_H3K27me3_cor <- cluster_corr_list$H3K27me3
H3K27me3_clust <- rna_H3K27me3_cor %>% t() %>% dist() %>% hclust()
H3K27me3_order <- H3K27me3_clust$labels[H3K27me3_clust$order]
rna_clust <- rna_H3K27me3_cor %>% dist() %>% hclust()
rna_order <- rna_clust$labels[rna_clust$order]

rna_H3K27me3_cor_df <- rna_H3K27me3_cor %>% 
    as_tibble(rownames='rna_clusters') %>% 
    pivot_longer(!rna_clusters, names_to='H3K27me3_clusters', values_to='corr') %>% 
    inner_join(cluster_meta, by=c('rna_clusters'='clusters')) %>% 
    left_join(select(H3K27me3_final, rna_clusters, mark_clusters, corr), by='rna_clusters') %>% 
    mutate(
        H3K27me3_clusters = factor(H3K27me3_clusters, levels=H3K27me3_order),
        rna_clusters = factor(rna_clusters, levels=unique(.$rna_clusters)),
        H3K27me3_stage = factor(str_replace(H3K27me3_clusters, '(\\w+)_\\d+', '\\1'), levels=stage_order),
        rna_stage = factor(str_replace(rna_clusters, '(\\w+)_\\d+', '\\1'), levels=stage_order),
        celltype_jf=factor(celltype_jf, levels=ct_order)
    ) 


ggplot(rna_H3K27me3_cor_df, aes(rna_clusters, H3K27me3_clusters, fill=corr.x)) +
    geom_tile() +
    geom_point(data=filter(rna_H3K27me3_cor_df, corr.x==corr.y), size=0.1, shape=4, stroke=0.2) +
    scale_fill_gradientn(limits=c(-0.45, 0.45), colors=rev(pals::brewer.spectral(100))) +
    facet_grid(H3K27me3_stage~celltype_jf, scales = 'free', space='free') +
    article_text() +
    no_x_text() + no_y_text() +
    theme(
        panel.border = element_blank(),
        panel.spacing = unit(0.02, 'cm'),
        plot.margin = unit(rep(0, 4), 'cm')
    ) 
ggsave('plots/paper/sfig1/H3K27me3_matches_heatmap2.pdf', width=8, height=10, units='cm')



rna_H3K4me3_cor <- cluster_corr_list$H3K4me3
H3K4me3_clust <- rna_H3K4me3_cor %>% t() %>% dist() %>% hclust()
H3K4me3_order <- H3K4me3_clust$labels[H3K4me3_clust$order]
rna_clust <- rna_H3K4me3_cor %>% dist() %>% hclust()
rna_order <- rna_clust$labels[rna_clust$order]

rna_H3K4me3_cor_df <- rna_H3K4me3_cor %>% 
    as_tibble(rownames='rna_clusters') %>% 
    pivot_longer(!rna_clusters, names_to='H3K4me3_clusters', values_to='corr') %>% 
    inner_join(cluster_meta, by=c('rna_clusters'='clusters')) %>% 
    left_join(select(H3K4me3_final, rna_clusters, mark_clusters, corr), by='rna_clusters') %>% 
    mutate(
        H3K4me3_clusters = factor(H3K4me3_clusters, levels=H3K4me3_order),
        rna_clusters = factor(rna_clusters, levels=unique(.$rna_clusters)),
        H3K4me3_stage = factor(str_replace(H3K4me3_clusters, '(\\w+)_\\d+', '\\1'), levels=stage_order),
        rna_stage = factor(str_replace(rna_clusters, '(\\w+)_\\d+', '\\1'), levels=stage_order),
        celltype_jf=factor(celltype_jf, levels=ct_order)
    ) 


ggplot(rna_H3K4me3_cor_df, aes(rna_clusters, H3K4me3_clusters, fill=corr.x)) +
    geom_tile() +
    geom_point(data=filter(rna_H3K4me3_cor_df, corr.x==corr.y), size=0.1, shape=4, stroke=0.2) +
    scale_fill_gradientn(limits=c(-0.7, 0.7), colors=rev(pals::brewer.spectral(100))) +
    facet_grid(H3K4me3_stage~celltype_jf, scales = 'free', space='free') +
    article_text() +
    no_x_text() + no_y_text() +
    theme(
        panel.border = element_blank(),
        panel.spacing = unit(0.02, 'cm'),
        plot.margin = unit(rep(0, 4), 'cm')
    ) 
ggsave('plots/paper/sfig1/H3K4me3_matches_heatmap.pdf', width=8, height=10, units='cm')



print_scale(rev(pals::brewer.spectral(100)))
ggsave('plots/paper/sfig1/matches_heatmap_scale.pdf', width=4, height=4, units='cm')





