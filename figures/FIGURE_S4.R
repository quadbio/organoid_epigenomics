source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist

setwd('~/projects/cutntag/')


#### Read data ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.3lines.rds')



#### Pseudotime dynamics for genes ####
#### CTX pseudotime clusters ####
H3K27ac_ctx_subs <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_subs_dpt_srt.rds')
H3K4me3_ctx_subs <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_subs_dpt_srt.rds')
H3K27me3_ctx_subs <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_subs_dpt_srt.rds')
rna_ctx_subs <- read_rds('data/trajectories/ctx/RNA_eb_ctx_subs_dpt_srt.rds')

mark_feats <- intersect(rownames(H3K27ac_ctx_subs[['cRNA']]), rownames(H3K27me3_ctx_subs[['cRNA']])) %>% intersect(rownames(H3K4me3_ctx_subs[['cRNA']]))


#### Plot age dist over pt bins #### 
H3K27ac_meta <- H3K27ac_ctx_subs@meta.data %>% as_tibble(rownames='cell') %>% filter(pt_bins>28)
H3K4me3_meta <- H3K4me3_ctx_subs@meta.data %>% as_tibble(rownames='cell') %>% filter(pt_bins>28)
H3K27me3_meta <- H3K27me3_ctx_subs@meta.data %>% as_tibble(rownames='cell') %>% filter(pt_bins>28)
rna_meta <- rna_ctx_subs@meta.data %>% as_tibble(rownames='cell') %>% filter(pt_bins>28)

p1 <- ggplot(H3K27ac_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27ac')

p2 <- ggplot(H3K4me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K4me3')

p3 <- ggplot(H3K27me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27me3')

p4 <- ggplot(rna_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('rna')

p1 / p2 / p3 / p4 & theme_void() & scale_fill_manual(values=colours_timescale)
ggsave('plots/paper/sfig4/sfig4_ctx_neurogen_age_dist_bar.pdf', width=8, height=5)


#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), ]
H3K27me3_clusters <- H3K27me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), ]
H3K4me3_clusters <- H3K4me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), ]
rna_clusters <- rna_ctx_subs@assays$RNA@misc$summary$pt_bins[as.character(28:50), ]

genes_plot <- c(
    'SOX2', 'VIM', 'GLI3', 'GRIA2', 'NEUROD1', 'DCX',
    'STMN2', 'BCL11A', 'NEUROD6', 'NEUROD2', 
    'UNCX', 'LHX5', 'LHX1', 'LHX2', 'OTX2', 'SCRT1', 'TUBB3', 'FOXG1', 'EMX1'
) %>% intersect(mark_feats)

H3K27ac_expr <- H3K27ac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

H3K27me3_expr <- H3K27me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

H3K4me3_expr <- H3K4me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr=scale01(expr))

p1 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    # geom_point(size=0.2) +
    scale_color_manual(values=modality_colors) +
    facet_grid(~gene, scales='free') +
    labs(x='Pseudotime bins', y='Expression', color='Modality')

p2 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=modality_colors) +
    facet_grid(modality~gene, scales='free') +
    labs(x='Pseudotime bins', y='Expression', color='Modality')

p1 / p2 + plot_layout(heights=c(1,3))
ggsave('plots/paper/sfig4/sfig4_ctx_neurogen_gene_expr_line.pdf', width=16, height=5)




#### DIEN pseudotime clusters ####
H3K27ac_dien_subs <- read_rds('data/trajectories/dien/H3K27ac_dien_subs_dpt_lsi_regress_srt.rds')
H3K4me3_dien_subs <- read_rds('data/trajectories/dien/H3K4me3_dien_subs_dpt_lsi_regress_srt.rds')
H3K27me3_dien_subs <- read_rds('data/trajectories/dien/H3K27me3_dien_subs_dpt_lsi_regress_srt.rds')
rna_dien_subs <- read_rds('data/trajectories/dien/RNA_dien_subs_dpt_srt.rds')

mark_feats <- intersect(rownames(H3K27ac_dien_subs[['cRNA']]), rownames(H3K27me3_dien_subs[['cRNA']])) %>% intersect(rownames(H3K4me3_dien_subs[['cRNA']]))

H3K27ac_dien_subs <- Pando::aggregate_assay(H3K27ac_dien_subs, assay='cRNA', group_name='pt_bins')
H3K27me3_dien_subs <- Pando::aggregate_assay(H3K27me3_dien_subs, assay='cRNA', group_name='pt_bins')
H3K4me3_dien_subs <- Pando::aggregate_assay(H3K4me3_dien_subs, assay='cRNA', group_name='pt_bins')
rna_dien_subs <- Pando::aggregate_assay(rna_dien_subs, assay='RNA', group_name='pt_bins')


#### Plot age dist over pt bins #### 
H3K27ac_meta <- H3K27ac_dien_subs@meta.data %>% as_tibble(rownames='cell') 
H3K4me3_meta <- H3K4me3_dien_subs@meta.data %>% as_tibble(rownames='cell') 
H3K27me3_meta <- H3K27me3_dien_subs@meta.data %>% as_tibble(rownames='cell') 
rna_meta <- rna_dien_subs@meta.data %>% as_tibble(rownames='cell') 

p1 <- ggplot(H3K27ac_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27ac')

p2 <- ggplot(H3K4me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K4me3')

p3 <- ggplot(H3K27me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27me3')

p4 <- ggplot(rna_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('rna')

p1 / p2 / p3 / p4 & theme_void() & scale_fill_manual(values=colours_timescale)
ggsave('plots/paper/sfig4/sfig4_dien_subs_age_dist_bar.pdf', width=8, height=5)


#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_clusters <- H3K27me3_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_clusters <- H3K4me3_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
rna_clusters <- rna_dien_subs@assays$RNA@misc$summary$pt_bins[as.character(1:20), ]


H3K27ac_expr <- H3K27ac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

H3K27me3_expr <- H3K27me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

H3K4me3_expr <- H3K4me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr=scale01(expr))

p1 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    # geom_point(size=0.2) +
    scale_color_manual(values=modality_colors) +
    facet_grid(~gene, scales='free') +
    labs(x='Pseudotime bins', y='Expression', color='Modality')

p2 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=modality_colors) +
    facet_grid(modality~gene, scales='free') +
    labs(x='Pseudotime bins', y='Expression', color='Modality')

p1 / p2 + plot_layout(heights=c(1,3))
ggsave('plots/paper/sfig4/sfig4_dien_neurogen_gene_expr_line.pdf', width=16, height=5)



#### NT pseudotime clusters ####
H3K27ac_nt_subs <- read_rds('data/trajectories/nt/H3K27ac_nt_subs_dpt_lsi_regress_srt.rds')
H3K4me3_nt_subs <- read_rds('data/trajectories/nt/H3K4me3_nt_subs_dpt_lsi_regress_srt.rds')
H3K27me3_nt_subs <- read_rds('data/trajectories/nt/H3K27me3_nt_subs_dpt_lsi_regress_srt.rds')
rna_nt_subs <- read_rds('data/trajectories/nt/RNA_nt_subs_dpt_srt.rds')

mark_feats <- intersect(rownames(H3K27ac_nt_subs[['cRNA']]), rownames(H3K27me3_nt_subs[['cRNA']])) %>% intersect(rownames(H3K4me3_nt_subs[['cRNA']]))

H3K27ac_nt_subs <- Pando::aggregate_assay(H3K27ac_nt_subs, assay='cRNA', group_name='pt_bins')
H3K27me3_nt_subs <- Pando::aggregate_assay(H3K27me3_nt_subs, assay='cRNA', group_name='pt_bins')
H3K4me3_nt_subs <- Pando::aggregate_assay(H3K4me3_nt_subs, assay='cRNA', group_name='pt_bins')
rna_nt_subs <- Pando::aggregate_assay(rna_nt_subs, assay='RNA', group_name='pt_bins')

#### Plot age dist over pt bins #### 
H3K27ac_meta <- H3K27ac_nt_subs@meta.data %>% as_tibble(rownames='cell') 
H3K4me3_meta <- H3K4me3_nt_subs@meta.data %>% as_tibble(rownames='cell') 
H3K27me3_meta <- H3K27me3_nt_subs@meta.data %>% as_tibble(rownames='cell') 
rna_meta <- rna_nt_subs@meta.data %>% as_tibble(rownames='cell') 

p1 <- ggplot(H3K27ac_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27ac')

p2 <- ggplot(H3K4me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K4me3')

p3 <- ggplot(H3K27me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27me3')

p4 <- ggplot(rna_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('rna')

p1 / p2 / p3 / p4 & theme_void() & scale_fill_manual(values=colours_timescale)
ggsave('plots/paper/sfig4/sfig4_nt_subs_age_dist_bar.pdf', width=8, height=5)


#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_clusters <- H3K27me3_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_clusters <- H3K4me3_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
rna_clusters <- rna_nt_subs@assays$RNA@misc$summary$pt_bins[as.character(1:20), ]


H3K27ac_expr <- H3K27ac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

H3K27me3_expr <- H3K27me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

H3K4me3_expr <- H3K4me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr=scale01(expr))

p1 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    # geom_point(size=0.2) +
    scale_color_manual(values=modality_colors) +
    facet_grid(~gene, scales='free') +
    labs(x='Pseudotime bins', y='Expression', color='Modality')

p2 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=modality_colors) +
    facet_grid(modality~gene, scales='free') +
    labs(x='Pseudotime bins', y='Expression', color='Modality')

p1 / p2 + plot_layout(heights=c(1,3))
ggsave('plots/paper/sfig4/sfig4_nt_neurogen_gene_expr_line.pdf', width=16, height=5)






#### Cluster neurogen peaks ####
H3K27ac_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K27ac_ctx_subs[['cRNA']]@data > 0)*1)
H3K27me3_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K27me3_ctx_subs[['cRNA']]@data > 0)*1)
H3K4me3_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K4me3_ctx_subs[['cRNA']]@data > 0)*1)

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='cRNA_bin', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='cRNA_bin', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='cRNA_bin', group_name='pt_bins')

mark_feats <- intersect(rownames(H3K27ac_ctx_subs[['cRNA']]), rownames(H3K27me3_ctx_subs[['cRNA']])) %>% intersect(rownames(H3K4me3_ctx_subs[['cRNA']]))

H3K27ac_cluster_detect <- H3K27ac_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), mark_feats]
H3K27me3_cluster_detect <- H3K27me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), mark_feats]
H3K4me3_cluster_detect <- H3K4me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), mark_feats]

H3K27ac_maxdetect_feats <- colMaxs(H3K27ac_cluster_detect)
H3K27me3_maxdetect_feats <- colMaxs(H3K27me3_cluster_detect)
H3K4me3_maxdetect_feats <- colMaxs(H3K4me3_cluster_detect)

mark_detect_feats <- mark_feats[
    (H3K27ac_maxdetect_feats > 0.02) & (H3K27me3_maxdetect_feats > 0.02) & (H3K4me3_maxdetect_feats > 0.02)
]

rna_ctx_subs <- FindVariableFeatures(rna_ctx_subs, nfeatures=6000)
H3K27ac_ctx_subs <- FindVariableFeatures(H3K27ac_ctx_subs, assay='cRNA', nfeatures=6000)
H3K27me3_ctx_subs <- FindVariableFeatures(H3K27me3_ctx_subs, assay='cRNA', nfeatures=6000)
H3K4me3_ctx_subs <- FindVariableFeatures(H3K4me3_ctx_subs, assay='cRNA', nfeatures=6000)

all_var_feats <- bind_rows(
    'RNA' = as_tibble(rna_ctx_subs[['RNA']]@meta.features, rownames='feature'),
    'H3K27ac' = as_tibble(H3K27ac_ctx_subs[['cRNA']]@meta.features, rownames='feature'),
    'H3K27me3' = as_tibble(H3K27me3_ctx_subs[['cRNA']]@meta.features, rownames='feature'),
    'H3K4me3' = as_tibble(H3K4me3_ctx_subs[['cRNA']]@meta.features, rownames='feature'),
    .id = 'modality'
)

rna_var_feats <- all_var_feats %>% 
    filter(modality=='RNA') %>% 
    top_n(10000, vst.variance) %>% 
    pull(feature) %>% 
    intersect(mark_detect_feats) -> feats_use



#### Smooth with GAMs ####
H3K27ac_clusters <- H3K27ac_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), ]
H3K27me3_clusters <- H3K27me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), ]
H3K4me3_clusters <- H3K4me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(28:50), ]
rna_clusters <- rna_ctx_subs@assays$RNA@misc$summary$pt_bins[as.character(28:50), ]

x <- 28:50
H3K27ac_gams <- H3K27ac_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
rna_smooths <- rna_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27ac_smooths_scale <- H3K27ac_smooths %>% scale()
H3K27me3_smooths_scale <- H3K27me3_smooths %>% scale()
H3K4me3_smooths_scale <- H3K4me3_smooths %>% scale()
rna_smooths_scale <- rna_smooths %>% scale()

all_smooths_list <- map(set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths_scale[,n],
            H3K27me3=H3K27me3_smooths_scale[,n],
            H3K4me3=H3K4me3_smooths_scale[,n],
            rna=rna_smooths_scale[,n]
        )
    )
})

all_smooths_raw_list <- map(set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths[,n],
            H3K27me3=H3K27me3_smooths[,n],
            H3K4me3=H3K4me3_smooths[,n],
            rna=rna_smooths[,n]
        )
    )
})

#### Distances with DTW and L2 (euclidean) ####
library(dtw)
library(doParallel)
registerDoParallel(36)

dtw_smooth_dist <- map_par(names(all_smooths_list), function(i){
    map_dfr(names(all_smooths_list), function(j){
        ij_align <- dtw(all_smooths_list[[i]], all_smooths_list[[j]])
        ij_dist <- ij_align$distance
        return(tibble(
            i=i, j=j, dist=ij_dist
        ))
    })
}, parallel=T)

dtw_smooth_dist_df <- bind_rows(dtw_smooth_dist)
dtw_smooth_dist_mat <- dtw_smooth_dist_df %>%
    pivot_wider(names_from=j, values_from=dist) %>%
    column_to_rownames('i') %>% as.matrix()

pheatmap::pheatmap(dtw_smooth_dist_mat)


#### K-means clustering ####
library(FCPS)

dtw_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 10
)

dtw_5_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 5
)

dtw_kmeans_df <- dtw_kmeans$Cls %>% enframe('feature', 'dtw_clust')
dtw_5_kmeans_df <- dtw_5_kmeans$Cls %>% enframe('feature', 'dtw_clust_5')

#### Plot clusters #####
all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
        {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=28:50) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr') %>% 
    group_by(feature, modality) %>% mutate(expr01=scale01(expr), pt_bin_mod=paste0(modality, pt_bin)) 

plot_df <- all_gene_smooths_df %>% 
    inner_join(dtw_kmeans_df) %>%
    inner_join(dtw_5_kmeans_df) %>%
    group_by(feature, modality)

plot_df$feature <- factor(plot_df$feature, levels=h_clust(plot_df, feature, pt_bin_mod, expr01))


ggplot(plot_df, aes(pt_bin, expr, color=modality)) +
    geom_smooth(alpha=0.5) +
    scale_color_manual(values=modality_colors) +
    facet_wrap(~dtw_clust_5) 
ggsave('plots/paper/sfig4/sfig4_ctx_neurogen_cluster_smooth.pdf', width=5, height=3)

ggplot(plot_df, aes(pt_bin, feature, fill=expr)) +
    geom_tile() +
    facet_grid(dtw_clust_5~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100)) +
    article_text()
ggsave('plots/paper/sfig4/sfig4_ctx_neurogen_cluster_heatmap.pdf', width=15, height=30)



#### GO enrichment ####
library(clusterProfiler)
library(org.Hs.eg.db)

feats_use <- rna@assays$RNA@meta.features %>% 
    as_tibble(rownames='feature') %>% 
    top_n(10000, vst.variance)

all_features <- bitr(feats_use$feature, fromType = 'SYMBOL', toType = c('ENSEMBL', 'ENTREZID'), OrgDb = org.Hs.eg.db) %>% 
    as_tibble()

cluster_genes <- plot_df %>% group_by(dtw_clust_5) %>% group_split()

cluster_GO_enrich <- map_dfr(cluster_genes, function(x){
    gene_ids <- filter(all_features, SYMBOL%in%x$feature)
    
    ggo <- enrichGO(
        gene = gene_ids$ENTREZID, 
        universe = all_features$ENTREZID,
        OrgDb = org.Hs.eg.db,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        pAdjustMethod = 'fdr',
        ont = 'BP',
        readable = T
    )
    return(as_tibble(ggo))
  
}, .id='cluster')

cluster_GO_enrich %>% write_tsv('data_/trajectories/ctx_neurogen_clusters_GO_enrich.tsv')






#### Cluster DIEN neurogen peaks ####
H3K27ac_dien_subs[['cRNA_bin']] <- CreateAssayObject((H3K27ac_dien_subs[['cRNA']]@data > 0)*1)
H3K27me3_dien_subs[['cRNA_bin']] <- CreateAssayObject((H3K27me3_dien_subs[['cRNA']]@data > 0)*1)
H3K4me3_dien_subs[['cRNA_bin']] <- CreateAssayObject((H3K4me3_dien_subs[['cRNA']]@data > 0)*1)

H3K27ac_dien_subs <- Pando::aggregate_assay(H3K27ac_dien_subs, assay='cRNA_bin', group_name='pt_bins')
H3K27me3_dien_subs <- Pando::aggregate_assay(H3K27me3_dien_subs, assay='cRNA_bin', group_name='pt_bins')
H3K4me3_dien_subs <- Pando::aggregate_assay(H3K4me3_dien_subs, assay='cRNA_bin', group_name='pt_bins')

mark_feats <- intersect(rownames(H3K27ac_dien_subs[['cRNA']]), rownames(H3K27me3_dien_subs[['cRNA']])) %>% intersect(rownames(H3K4me3_dien_subs[['cRNA']]))

H3K27ac_cluster_detect <- H3K27ac_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]
H3K27me3_cluster_detect <- H3K27me3_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]
H3K4me3_cluster_detect <- H3K4me3_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]

H3K27ac_maxdetect_feats <- colMaxs(H3K27ac_cluster_detect)
H3K27me3_maxdetect_feats <- colMaxs(H3K27me3_cluster_detect)
H3K4me3_maxdetect_feats <- colMaxs(H3K4me3_cluster_detect)

mark_detect_feats <- mark_feats[
    (H3K27ac_maxdetect_feats > 0.03) & (H3K27me3_maxdetect_feats > 0.03) & (H3K4me3_maxdetect_feats > 0.03)
]

rna_dien_subs <- FindVariableFeatures(rna_dien_subs, nfeatures=6000)
H3K27ac_dien_subs <- FindVariableFeatures(H3K27ac_dien_subs, assay='cRNA', nfeatures=6000)
H3K27me3_dien_subs <- FindVariableFeatures(H3K27me3_dien_subs, assay='cRNA', nfeatures=6000)
H3K4me3_dien_subs <- FindVariableFeatures(H3K4me3_dien_subs, assay='cRNA', nfeatures=6000)

all_var_feats <- bind_rows(
    'RNA' = as_tibble(rna_dien_subs[['RNA']]@meta.features, rownames='feature'),
    'H3K27ac' = as_tibble(H3K27ac_dien_subs[['cRNA']]@meta.features, rownames='feature'),
    'H3K27me3' = as_tibble(H3K27me3_dien_subs[['cRNA']]@meta.features, rownames='feature'),
    'H3K4me3' = as_tibble(H3K4me3_dien_subs[['cRNA']]@meta.features, rownames='feature'),
    .id = 'modality'
)

rna_var_feats <- all_var_feats %>% 
    filter(modality=='RNA') %>% 
    top_n(10000, vst.variance) %>% 
    pull(feature) %>% 
    intersect(mark_detect_feats) -> feats_use



#### Smooth with GAMs ####
H3K27ac_clusters <- H3K27ac_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_clusters <- H3K27me3_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_clusters <- H3K4me3_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
rna_clusters <- rna_dien_subs@assays$RNA@misc$summary$pt_bins[as.character(1:20), ]

x <- 1:20
H3K27ac_gams <- H3K27ac_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
rna_smooths <- rna_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27ac_smooths_scale <- H3K27ac_smooths %>% scale()
H3K27me3_smooths_scale <- H3K27me3_smooths %>% scale()
H3K4me3_smooths_scale <- H3K4me3_smooths %>% scale()
rna_smooths_scale <- rna_smooths %>% scale()

all_smooths_list <- map(set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths_scale[,n],
            H3K27me3=H3K27me3_smooths_scale[,n],
            H3K4me3=H3K4me3_smooths_scale[,n],
            rna=rna_smooths_scale[,n]
        )
    )
})

all_smooths_raw_list <- map(set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths[,n],
            H3K27me3=H3K27me3_smooths[,n],
            H3K4me3=H3K4me3_smooths[,n],
            rna=rna_smooths[,n]
        )
    )
})

#### Distances with DTW and L2 (euclidean) ####
library(dtw)
library(doParallel)
registerDoParallel(36)

dtw_smooth_dist <- map_par(names(all_smooths_list), function(i){
    map_dfr(names(all_smooths_list), function(j){
        ij_align <- dtw(all_smooths_list[[i]], all_smooths_list[[j]])
        ij_dist <- ij_align$distance
        return(tibble(
            i=i, j=j, dist=ij_dist
        ))
    })
}, parallel=T)

dtw_smooth_dist_df <- bind_rows(dtw_smooth_dist)
dtw_smooth_dist_mat <- dtw_smooth_dist_df %>%
    pivot_wider(names_from=j, values_from=dist) %>%
    column_to_rownames('i') %>% as.matrix()

# pheatmap::pheatmap(dtw_smooth_dist_mat)


#### K-means clustering ####
library(FCPS)

dtw_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 10
)

dtw_5_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 5
)

dtw_kmeans_df <- dtw_kmeans$Cls %>% enframe('feature', 'dtw_clust')
dtw_5_kmeans_df <- dtw_5_kmeans$Cls %>% enframe('feature', 'dtw_clust_5')

#### Plot clusters #####
all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
        {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=1:20) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr') %>% 
    group_by(feature, modality) %>% mutate(expr01=scale01(expr), pt_bin_mod=paste0(modality, pt_bin)) 

plot_df <- all_gene_smooths_df %>% 
    inner_join(dtw_kmeans_df) %>%
    inner_join(dtw_5_kmeans_df) %>%
    group_by(feature, modality)

plot_df$feature <- factor(plot_df$feature, levels=h_clust(plot_df, feature, pt_bin_mod, expr01))


ggplot(plot_df, aes(pt_bin, expr01, color=modality)) +
    geom_smooth(alpha=0.5) +
    scale_color_manual(values=modality_colors) +
    facet_wrap(~dtw_clust_5) 
ggsave('plots/paper/sfig4/sfig4_dien_neurogen_cluster_smooth.pdf', width=5, height=3)

ggplot(plot_df, aes(pt_bin, feature, fill=expr)) +
    geom_tile() +
    facet_grid(dtw_clust_5~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100)) +
    article_text()
ggsave('plots/paper/sfig4/sfig4_dien_neurogen_cluster_heatmap.pdf', width=15, height=30)




#### GO enrichment ####
library(clusterProfiler)
library(org.Hs.eg.db)

feats_use <- rna@assays$RNA@meta.features %>% 
    as_tibble(rownames='feature') %>% 
    top_n(10000, vst.variance)

all_features <- bitr(feats_use$feature, fromType = 'SYMBOL', toType = c('ENSEMBL', 'ENTREZID'), OrgDb = org.Hs.eg.db) %>% 
    as_tibble()

cluster_genes <- plot_df %>% group_by(dtw_clust_5) %>% group_split()

cluster_GO_enrich <- map_dfr(cluster_genes, function(x){
    gene_ids <- filter(all_features, SYMBOL%in%x$feature)
    
    ggo <- enrichGO(
        gene = gene_ids$ENTREZID, 
        universe = all_features$ENTREZID,
        OrgDb = org.Hs.eg.db,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        pAdjustMethod = 'fdr',
        ont = 'BP',
        readable = T
    )
    return(as_tibble(ggo))
    
}, .id='cluster')

cluster_GO_enrich %>% write_tsv('data_/trajectories/dien_neurogen_clusters_GO_enrich.tsv')





#### Cluster NT neurogen peaks ####
H3K27ac_nt_subs[['cRNA_bin']] <- CreateAssayObject((H3K27ac_nt_subs[['cRNA']]@data > 0)*1)
H3K27me3_nt_subs[['cRNA_bin']] <- CreateAssayObject((H3K27me3_nt_subs[['cRNA']]@data > 0)*1)
H3K4me3_nt_subs[['cRNA_bin']] <- CreateAssayObject((H3K4me3_nt_subs[['cRNA']]@data > 0)*1)

H3K27ac_nt_subs <- Pando::aggregate_assay(H3K27ac_nt_subs, assay='cRNA_bin', group_name='pt_bins')
H3K27me3_nt_subs <- Pando::aggregate_assay(H3K27me3_nt_subs, assay='cRNA_bin', group_name='pt_bins')
H3K4me3_nt_subs <- Pando::aggregate_assay(H3K4me3_nt_subs, assay='cRNA_bin', group_name='pt_bins')

mark_feats <- intersect(rownames(H3K27ac_nt_subs[['cRNA']]), rownames(H3K27me3_nt_subs[['cRNA']])) %>% intersect(rownames(H3K4me3_nt_subs[['cRNA']]))

H3K27ac_cluster_detect <- H3K27ac_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]
H3K27me3_cluster_detect <- H3K27me3_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]
H3K4me3_cluster_detect <- H3K4me3_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]

H3K27ac_maxdetect_feats <- colMaxs(H3K27ac_cluster_detect)
H3K27me3_maxdetect_feats <- colMaxs(H3K27me3_cluster_detect)
H3K4me3_maxdetect_feats <- colMaxs(H3K4me3_cluster_detect)

mark_detect_feats <- mark_feats[
    (H3K27ac_maxdetect_feats > 0.02) & (H3K27me3_maxdetect_feats > 0.02) & (H3K4me3_maxdetect_feats > 0.02)
]

rna_nt_subs <- FindVariableFeatures(rna_nt_subs, nfeatures=6000)
H3K27ac_nt_subs <- FindVariableFeatures(H3K27ac_nt_subs, assay='cRNA', nfeatures=6000)
H3K27me3_nt_subs <- FindVariableFeatures(H3K27me3_nt_subs, assay='cRNA', nfeatures=6000)
H3K4me3_nt_subs <- FindVariableFeatures(H3K4me3_nt_subs, assay='cRNA', nfeatures=6000)

all_var_feats <- bind_rows(
    'RNA' = as_tibble(rna_nt_subs[['RNA']]@meta.features, rownames='feature'),
    'H3K27ac' = as_tibble(H3K27ac_nt_subs[['cRNA']]@meta.features, rownames='feature'),
    'H3K27me3' = as_tibble(H3K27me3_nt_subs[['cRNA']]@meta.features, rownames='feature'),
    'H3K4me3' = as_tibble(H3K4me3_nt_subs[['cRNA']]@meta.features, rownames='feature'),
    .id = 'modality'
)

rna_var_feats <- all_var_feats %>% 
    filter(modality=='RNA') %>% 
    top_n(10000, vst.variance) %>% 
    pull(feature) %>% 
    intersect(mark_detect_feats) -> feats_use



#### Smooth with GAMs ####
H3K27ac_clusters <- H3K27ac_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_clusters <- H3K27me3_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_clusters <- H3K4me3_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
rna_clusters <- rna_nt_subs@assays$RNA@misc$summary$pt_bins[as.character(1:20), ]

x <- 1:20
H3K27ac_gams <- H3K27ac_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = scale01(y) ~ s(x, bs = 'cs'))
    })
rna_smooths <- rna_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27ac_smooths_scale <- H3K27ac_smooths %>% scale()
H3K27me3_smooths_scale <- H3K27me3_smooths %>% scale()
H3K4me3_smooths_scale <- H3K4me3_smooths %>% scale()
rna_smooths_scale <- rna_smooths %>% scale()

all_smooths_list <- map(set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths_scale[,n],
            H3K27me3=H3K27me3_smooths_scale[,n],
            H3K4me3=H3K4me3_smooths_scale[,n],
            rna=rna_smooths_scale[,n]
        )
    )
})

all_smooths_raw_list <- map(set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths[,n],
            H3K27me3=H3K27me3_smooths[,n],
            H3K4me3=H3K4me3_smooths[,n],
            rna=rna_smooths[,n]
        )
    )
})

#### Distances with DTW and L2 (euclidean) ####
library(dtw)
library(doParallel)
registerDoParallel(36)

dtw_smooth_dist <- map_par(names(all_smooths_list), function(i){
    map_dfr(names(all_smooths_list), function(j){
        ij_align <- dtw(all_smooths_list[[i]], all_smooths_list[[j]])
        ij_dist <- ij_align$distance
        return(tibble(
            i=i, j=j, dist=ij_dist
        ))
    })
}, parallel=T)

dtw_smooth_dist_df <- bind_rows(dtw_smooth_dist)
dtw_smooth_dist_mat <- dtw_smooth_dist_df %>%
    pivot_wider(names_from=j, values_from=dist) %>%
    column_to_rownames('i') %>% as.matrix()

# pheatmap::pheatmap(dtw_smooth_dist_mat)


#### K-means clustering ####
library(FCPS)

dtw_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 10
)

dtw_5_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 5
)

dtw_kmeans_df <- dtw_kmeans$Cls %>% enframe('feature', 'dtw_clust')
dtw_5_kmeans_df <- dtw_5_kmeans$Cls %>% enframe('feature', 'dtw_clust_5')

#### Plot clusters #####
all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
        {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=1:20) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr') %>% 
    group_by(feature, modality) %>% mutate(expr01=scale01(expr), pt_bin_mod=paste0(modality, pt_bin)) 

plot_df <- all_gene_smooths_df %>% 
    inner_join(dtw_kmeans_df) %>%
    inner_join(dtw_5_kmeans_df) %>%
    group_by(feature, modality)

plot_df$feature <- factor(plot_df$feature, levels=h_clust(plot_df, feature, pt_bin_mod, expr01))


ggplot(plot_df, aes(pt_bin, expr01, color=modality)) +
    geom_smooth(alpha=0.5) +
    scale_color_manual(values=modality_colors) +
    facet_wrap(~dtw_clust_5) 
ggsave('plots/paper/sfig4/sfig4_nt_neurogen_cluster_smooth.pdf', width=5, height=3)

ggplot(plot_df, aes(pt_bin, feature, fill=expr)) +
    geom_tile() +
    facet_grid(dtw_clust_5~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100)) +
    article_text()
ggsave('plots/paper/sfig4/sfig4_nt_neurogen_cluster_heatmap.pdf', width=15, height=30)



#### GO enrichment ####
library(clusterProfiler)
library(org.Hs.eg.db)

feats_use <- rna@assays$RNA@meta.features %>% 
    as_tibble(rownames='feature') %>% 
    top_n(10000, vst.variance)

all_features <- bitr(feats_use$feature, fromType = 'SYMBOL', toType = c('ENSEMBL', 'ENTREZID'), OrgDb = org.Hs.eg.db) %>% 
    as_tibble()

cluster_genes <- plot_df %>% group_by(dtw_clust_5) %>% group_split()

cluster_GO_enrich <- map_dfr(cluster_genes, function(x){
    gene_ids <- filter(all_features, SYMBOL%in%x$feature)
    
    ggo <- enrichGO(
        gene = gene_ids$ENTREZID, 
        universe = all_features$ENTREZID,
        OrgDb = org.Hs.eg.db,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        pAdjustMethod = 'fdr',
        ont = 'BP',
        readable = T
    )
    return(as_tibble(ggo))
    
}, .id='cluster')

cluster_GO_enrich %>% write_tsv('data_/trajectories/nt_neurogen_clusters_GO_enrich.tsv')












#### DE/DA analysis for neuron vs NPC scatter plots ####
#### Read data ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')

cluster_srt <- read_rds('data/RNA/all_RNA_marks_combined_clusters_srt.rds')
cluster_graph <- read_rds('data/RNA/all_RNA_cluster_graph.rds')


# DE results 
H3K27ac_nt_neuron_da <- read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_nt_NvsNPC.tsv')
H3K27ac_ctx_neuron_da <- read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_ctx_NvsNPC.tsv')
H3K27ac_dien_neuron_da <- read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_dien_NvsNPC.tsv')
H3K27ac_all_neuron_da <- read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_all_NvsNPC.tsv')

H3K27ac_neuron_da <- bind_rows('nt'=H3K27ac_nt_neuron_da, 'ctx'=H3K27ac_ctx_neuron_da, 'dien'=H3K27ac_dien_neuron_da, 'all'=H3K27ac_all_neuron_da, .id='lineage')

H3K27me3_nt_neuron_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_nt_NvsNPC.tsv')
H3K27me3_ctx_neuron_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_ctx_NvsNPC.tsv')
H3K27me3_dien_neuron_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_dien_NvsNPC.tsv')
H3K27me3_all_neuron_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_all_NvsNPC.tsv')

H3K27me3_neuron_da <- bind_rows('nt'=H3K27me3_nt_neuron_da, 'ctx'=H3K27me3_ctx_neuron_da, 'dien'=H3K27me3_dien_neuron_da, 'all'=H3K27me3_all_neuron_da, .id='lineage')

H3K4me3_nt_neuron_da <- read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_nt_NvsNPC.tsv')
H3K4me3_ctx_neuron_da <- read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_ctx_NvsNPC.tsv')
H3K4me3_dien_neuron_da <- read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_dien_NvsNPC.tsv')
H3K4me3_all_neuron_da <- read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_all_NvsNPC.tsv')

H3K4me3_neuron_da <- bind_rows('nt'=H3K4me3_nt_neuron_da, 'ctx'=H3K4me3_ctx_neuron_da, 'dien'=H3K4me3_dien_neuron_da, 'all'=H3K4me3_all_neuron_da, .id='lineage')

neuron_da <- bind_rows('H3K27ac'=H3K27ac_neuron_da, 'H3K27me3'=H3K27me3_neuron_da, 'H3K4me3'=H3K4me3_neuron_da, .id='mark')

peak_genes <- ClosestFeature(marks$H3K27me3, neuron_da$feature) %>% as_tibble() %>% distinct(query_region, gene_name)

neuron_da <- neuron_da %>% 
    inner_join(peak_genes, by=c('feature'='query_region'))

neuron_da %>% write_tsv('data/results/diff_expression/all_marks_DA_NvsNPC.tsv')


lineage_da <- read_tsv('data/results/diff_expression/all_marks_DA_lineage_coarse.tsv')

peak_genes <- ClosestFeature(marks$H3K27me3, lineage_da$feature) %>% as_tibble() %>% distinct(query_region, gene_name)

lineage_da <- lineage_da %>% 
    inner_join(peak_genes, by=c('feature'='query_region'))

lin_pairs_da <- bind_rows(
    H3K27ac=read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_lineages_pairwise.tsv'),
    H3K27me3=read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_lineages_pairwise.tsv'),
    H3K4me3=read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_lineages_pairwise.tsv'),
    .id='mark'
)

peak_genes <- ClosestFeature(marks$H3K27me3, lin_pairs_da$feature) %>% as_tibble() %>% distinct(query_region, gene_name)

lin_pairs_da <- lin_pairs_da %>% 
    inner_join(peak_genes, by=c('feature'='query_region'))

lin_pairs_da %>% write_tsv('data/results/diff_expression/all_marks_DA_peaks_lineages_pairwise.tsv')

#### Do DE for RNA ####

rna_nt <- read_rds('data/trajectories/nt/RNA_nt_dpt_srt.rds')
rna_dien <- read_rds('data/trajectories/dien/RNA_dien_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_ctx_dpt_srt.rds')
rna_ac <- read_rds('data/trajectories/astrocytes/RNA_astrocytes_npcs_srt.rds')


# NPC vs neuron
npc_neuron_nt_de <- de(rna_nt, 'state') %>% filter(group=='neuron')
npc_neuron_dien_de <- de(rna_dien, 'state') %>% filter(group=='neuron')
npc_neuron_ctx_de <- de(rna_ctx, 'state') %>% filter(group=='neuron')

rna_nt$lineage <- 'nt'
subs_npc <- sample(colnames(rna_nt)[rna_nt$state=='npc'], 2500)
subs_neuron <- sample(colnames(rna_nt)[rna_nt$state=='neuron'], 500)
rna_nt_subs <- subset(rna_nt, cells=c(subs_npc, subs_neuron))

rna_dien$lineage <- 'dien'
subs_cells <- sample(colnames(rna_dien), 3000)
rna_dien_subs <- subset(rna_dien, cells=subs_cells)

rna_ctx$lineage <- 'ctx'
subs_npc <- sample(colnames(rna_ctx)[rna_ctx$state=='npc'], 2500)
subs_neuron <- sample(colnames(rna_ctx)[rna_ctx$state=='neuron'], 500)
rna_ctx_subs <- subset(rna_ctx, cells=c(subs_npc, subs_neuron))


rna_lins <- merge(rna_nt_subs, list(rna_dien_subs, rna_ctx_subs))

npc_neuron_all <- de(rna_lins, 'state') %>% 
    filter(group=='neuron') %>% 
    select(feature, 'neuron_fc'=fc, 'neuron_pval'=padj)

lineage_combs <- list(c('nt', 'dien'), c('nt', 'ctx'), c('dien', 'ctx'))
lin_pairs_de <- map_dfr(lineage_combs, function(x){
    srt_use <- subset(rna_lins, lineage%in%c(x))
    srt_use$test_var <- srt_use$lineage==x[2]
    de(srt_use, 'test_var') %>% 
        mutate(group1=x[2], group2=x[1])
})

lin_pairs_de %>% write_tsv('data/results/diff_expression/RNA_DE_lineages_pairwise.tsv')


#### CTX vs DIEN ####
sig_colors <- c('#57B2B5', '#BD71CE', '#1660A1', 'darkgrey')
names(sig_colors) <- c('lin', 'neuron', 'both', 'none')

ctx_dien_lin_de <- lin_pairs_de %>% 
    filter(group1=='ctx', group2=='dien', group==T) %>% 
    select(feature, 'lin_fc'=fc, 'lin_pval'=padj)

ctx_dien_de <- inner_join(npc_neuron_all, ctx_dien_lin_de) %>% 
    filter(abs(lin_fc)>0 | abs(neuron_fc)>0) %>% 
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_fc)>0.25 | abs(lin_fc)>0.25),
        sig_group=case_when(
            !(lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'neuron',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & !(neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'lin',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'both',
            T ~ 'none'
        )
    )

label_data <- filter(ctx_dien_de, abs(lin_fc)>1 | abs(neuron_fc)>1)

p1 <- ggplot(ctx_dien_de, aes(lin_fc, neuron_fc, label=feature, alpha=sig, color=sig_group, size=abs(neuron_fc))) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_point() +
    scale_x_continuous(limits=c(-2,2)) +
    scale_y_continuous(limits=c(-2,2)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    geom_text(data=label_data, color='black', size=2) +
    scale_size_continuous(range=c(0.1,0.5)) +
    scale_color_manual(values=sig_colors) +
    article_text() +
    no_legend() +
    labs(x='<- DIEN vs CTX ->', y='<- NPC vs NEURON ->')


ctx_dien_de_groups <- ctx_dien_de %>% 
    distinct(feature, sig_group) %>% 
    filter(sig_group!='none')

# -> compare with chromatin
ctx_dien_lin_da <- lin_pairs_da %>% 
    filter(group1=='ctx', group2=='dien') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'lin_coef'=coef, 'lin_pval'=padj, 'lin_dr'=log_dr, gene_name)

npc_neuron_da <- neuron_da %>% 
    filter(lineage=='all') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'neuron_coef'=coef, 'neuron_pval'=padj, 'neuron_dr'=log_dr, gene_name)

ctx_dien_da <- inner_join(npc_neuron_da, ctx_dien_lin_da) %>% 
    inner_join(ctx_dien_de_groups, by=c('gene_name'='feature')) %>%
    mutate(lin_coef_clip=clip_abs(lin_coef, 5), neuron_coef_clip=clip_abs(neuron_coef, 5)) %>% 
    filter(abs(lin_coef)>0 | abs(neuron_coef)>0) %>%
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_coef)>0.25 | abs(lin_coef)>0.25)
    )

jitterpos <- position_jitter(seed=111)

right_plot_df <- ctx_dien_da %>% mutate(sig_group=ifelse(sig_group=='lin', 'none', sig_group))
p_right <- ggplot(right_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(right_plot_df, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    rotate_x_text(90) +
    theme_rangeframe() + scale_axis_rangeframe() +
    no_legend() +
    no_label()

top_plot_df <- ctx_dien_da %>% mutate(sig_group=ifelse(sig_group=='neuron', 'none', sig_group))
p_top <- ggplot(top_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(ctx_dien_da, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    coord_flip() +
    no_legend() +
    no_label()

layout <- '
AAAAA#
BBBBBC
BBBBBC
BBBBBC
BBBBBC
BBBBBC
'

p_top + p1 + p_right + plot_layout(design = layout)

ggsave('plots/paper/sfig4/sfig4_ctx_dien_de_da_scatter.pdf', width=10, height=10, units='cm')




#### CTX vs NT ####
ctx_nt_lin_de <- lin_pairs_de %>% 
    filter(group1=='ctx', group2=='nt', group==T) %>% 
    select(feature, 'lin_fc'=fc, 'lin_pval'=padj)

ctx_nt_de <- inner_join(npc_neuron_all, ctx_nt_lin_de) %>% 
    filter(abs(lin_fc)>0 | abs(neuron_fc)>0) %>% 
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_fc)>0.25 | abs(lin_fc)>0.25),
        sig_group=case_when(
            !(lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'neuron',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & !(neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'lin',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'both',
            T ~ 'none'
        )
    )

label_data <- filter(ctx_nt_de, abs(lin_fc)>1 | abs(neuron_fc)>1)

p1 <- ggplot(ctx_nt_de, aes(lin_fc, neuron_fc, label=feature, alpha=sig, color=sig_group, size=abs(neuron_fc))) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_point() +
    scale_x_continuous(limits=c(-2,2)) +
    scale_y_continuous(limits=c(-2,2)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    geom_text(data=label_data, color='black', size=2) +
    scale_size_continuous(range=c(0.1,0.5)) +
    scale_color_manual(values=sig_colors) +
    article_text() +
    no_legend() +
    labs(x='<- DIEN vs NT ->', y='<- NPC vs NEURON ->')


ctx_nt_de_groups <- ctx_nt_de %>% 
    distinct(feature, sig_group) %>% 
    filter(sig_group!='none')

# -> compare with chromatin
ctx_nt_lin_da <- lin_pairs_da %>% 
    filter(group1=='ctx', group2=='nt') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'lin_coef'=coef, 'lin_pval'=padj, 'lin_dr'=log_dr, gene_name)

npc_neuron_da <- neuron_da %>% 
    filter(lineage=='all') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'neuron_coef'=coef, 'neuron_pval'=padj, 'neuron_dr'=log_dr, gene_name)

ctx_nt_da <- inner_join(npc_neuron_da, ctx_nt_lin_da) %>% 
    inner_join(ctx_nt_de_groups, by=c('gene_name'='feature')) %>%
    mutate(lin_coef_clip=clip_abs(lin_coef, 5), neuron_coef_clip=clip_abs(neuron_coef, 5)) %>% 
    filter(abs(lin_coef)>0 | abs(neuron_coef)>0) %>%
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_coef)>0.25 | abs(lin_coef)>0.25)
    )

jitterpos <- position_jitter(seed=111)

right_plot_df <- ctx_nt_da %>% mutate(sig_group=ifelse(sig_group=='lin', 'none', sig_group))
p_right <- ggplot(right_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(right_plot_df, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    rotate_x_text(90) +
    theme_rangeframe() + scale_axis_rangeframe() +
    no_legend() +
    no_label()

top_plot_df <- ctx_nt_da %>% mutate(sig_group=ifelse(sig_group=='neuron', 'none', sig_group))
p_top <- ggplot(top_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(ctx_nt_da, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    coord_flip() +
    no_legend() +
    no_label()

layout <- '
AAAAA#
BBBBBC
BBBBBC
BBBBBC
BBBBBC
BBBBBC
'

p_top + p1 + p_right + plot_layout(design = layout)

ggsave('plots/paper/sfig4/sfig4_ctx_nt_de_da_scatter.pdf', width=10, height=10, units='cm')




#### CTX vs NT ####
dien_nt_lin_de <- lin_pairs_de %>% 
    filter(group1=='dien', group2=='nt', group==T) %>% 
    select(feature, 'lin_fc'=fc, 'lin_pval'=padj)

dien_nt_de <- inner_join(npc_neuron_all, dien_nt_lin_de) %>% 
    filter(abs(lin_fc)>0 | abs(neuron_fc)>0) %>% 
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_fc)>0.25 | abs(lin_fc)>0.25),
        sig_group=case_when(
            !(lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'neuron',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & !(neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'lin',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'both',
            T ~ 'none'
        )
    )

label_data <- filter(dien_nt_de, abs(lin_fc)>1 | abs(neuron_fc)>1)

p1 <- ggplot(dien_nt_de, aes(lin_fc, neuron_fc, label=feature, alpha=sig, color=sig_group, size=abs(neuron_fc))) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_point() +
    scale_x_continuous(limits=c(-2,2)) +
    scale_y_continuous(limits=c(-2,2)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    geom_text(data=label_data, color='black', size=2) +
    scale_size_continuous(range=c(0.1,0.5)) +
    scale_color_manual(values=sig_colors) +
    article_text() +
    no_legend() +
    labs(x='<- DIEN vs NT ->', y='<- NPC vs NEURON ->')


dien_nt_de_groups <- dien_nt_de %>% 
    distinct(feature, sig_group) %>% 
    filter(sig_group!='none')

# -> compare with chromatin
dien_nt_lin_da <- lin_pairs_da %>% 
    filter(group1=='dien', group2=='nt') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'lin_coef'=coef, 'lin_pval'=padj, 'lin_dr'=log_dr, gene_name)

npc_neuron_da <- neuron_da %>% 
    filter(lineage=='all') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'neuron_coef'=coef, 'neuron_pval'=padj, 'neuron_dr'=log_dr, gene_name)

dien_nt_da <- inner_join(npc_neuron_da, dien_nt_lin_da) %>% 
    inner_join(dien_nt_de_groups, by=c('gene_name'='feature')) %>%
    mutate(lin_coef_clip=clip_abs(lin_coef, 5), neuron_coef_clip=clip_abs(neuron_coef, 5)) %>% 
    filter(abs(lin_coef)>0 | abs(neuron_coef)>0) %>%
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_coef)>0.25 | abs(lin_coef)>0.25)
    )

jitterpos <- position_jitter(seed=111)

right_plot_df <- dien_nt_da %>% mutate(sig_group=ifelse(sig_group=='lin', 'none', sig_group))
p_right <- ggplot(right_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(right_plot_df, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    rotate_x_text(90) +
    theme_rangeframe() + scale_axis_rangeframe() +
    no_legend() +
    no_label()

top_plot_df <- dien_nt_da %>% mutate(sig_group=ifelse(sig_group=='neuron', 'none', sig_group))
p_top <- ggplot(top_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(dien_nt_da, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    coord_flip() +
    no_legend() +
    no_label()

layout <- '
AAAAA#
BBBBBC
BBBBBC
BBBBBC
BBBBBC
BBBBBC
'

p_top + p1 + p_right + plot_layout(design = layout)

ggsave('plots/paper/sfig4/sfig4_dien_nt_de_da_scatter.pdf', width=10, height=10, units='cm')







