source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(destiny)
library(dtw)

setwd('~/projects/cutntag/')

rank_cols <- function(x){
    xr <- presto::rank_matrix(t(x))$X_ranked
    colnames(xr) <- colnames(x)
    rownames(xr) <- rownames(x)
    return(xr)
}


#### Read data ###
marks <- read_rds('data/all_marks_list_v3.4lines.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.3lines.rds')
rna_pt_meta <- read_tsv('data/RNA/cellrank/RNA_full_cellrank_probs.tsv') %>% 
    dplyr::rename('cell'=1) %>% 
    select(cell, velocity_pseudotime, pseudotime_ranks) %>% 
    column_to_rownames('cell')

rna <- AddMetaData(rna, rna_pt_meta)

#### Prep data ####
marks <- map(marks, function(x){
    x$lineage2 <- case_when(
        x$celltype_jf=='astrocytes' ~ 'astrocytes',
        T ~ x$lineage
    )
    return(x)
})

rna$lineage2 <- case_when(
    rna$celltype_jf=='astrocytes' ~ 'astrocytes',
    T ~ rna$lineage
)

marks_lin <- marks %>% map(SplitObject, split.by='lineage2')
rna_lin <- rna %>% SplitObject(split.by='lineage2')

marks_lin <- purrr::transpose(marks_lin)

#### Run diffmap for H3K27ac cortex ####
marks_ctx <- marks_lin$ctx

H3K27ac_ctx <- marks_ctx$H3K27ac

H3K27ac_ctx <- H3K27ac_ctx %>% 
    FindTopFeatures(min.cutoff='q90') %>% 
    RunSVD() 

H3K27ac_ctx <- H3K27ac_ctx %>% 
    RunUMAP(dims=2:30, reduction='lsi')

DepthCor(H3K27ac_ctx)

H3K27ac_pca <- H3K27ac_ctx[['lsi']]@cell.embeddings[,2:20]
H3K27ac_diffmap <- DiffusionMap(H3K27ac_pca, verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5]))

H3K27ac_dpt <- DPT(H3K27ac_diffmap)

plot(H3K27ac_dpt, col_by='DPT3')

H3K27ac_dpt_df <- H3K27ac_dpt[,1:5]
colnames(H3K27ac_dpt_df) <- paste0('DPT', 1:5)
rownames(H3K27ac_dpt_df) <- colnames(H3K27ac_ctx)
H3K27ac_dpt_df <- as.data.frame(rank_cols(H3K27ac_dpt_df))

H3K27ac_ctx$dpt_tip <- H3K27ac_dpt@tips[,'Tips1']
H3K27ac_ctx[['diffmap']] <- CreateDimReducObject(
    H3K27ac_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27ac_ctx <- AddMetaData(H3K27ac_ctx, H3K27ac_dc_df)
H3K27ac_ctx <- AddMetaData(H3K27ac_ctx, H3K27ac_dpt_df)

H3K27ac_ctx %>% 
    feature_plot(features=c('STMN2', 'VIM', 'DC1', 'DC2', 'DC3', 'DPT1', 'DPT2', 'DPT3'), order=T)

H3K27ac_ctx %>% write_rds('data/trajectories/ctx/H3K27ac_ctx_dpt_srt.rds')
H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_ctx_dpt_srt.rds')


#### Run diffmap for H3K4me3 cortex ####
H3K4me3_ctx <- marks_ctx$H3K4me3[, !str_detect(marks_ctx$H3K4me3$orig.ident, '_ret')]

H3K4me3_ctx <- H3K4me3_ctx %>% 
    FindTopFeatures(min.cutoff='q90') %>% 
    RunSVD() 

H3K4me3_ctx <- H3K4me3_ctx %>% 
    RunUMAP(dims=2:30, reduction='lsi')

DepthCor(H3K4me3_ctx)

H3K4me3_ctx %>% dim_plot(group.by=c('orig.ident', 'state', 'clusters'), label=T, order=T)
H3K4me3_ctx %>% feature_plot(features=c('NEUROD6', 'EMX1', 'SOX2'), order=T)

H3K4me3_pca <- H3K4me3_ctx[['lsi']]@cell.embeddings[,2:20]
H3K4me3_diffmap <- DiffusionMap(H3K4me3_pca, verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5]))

H3K4me3_dpt <- DPT(H3K4me3_diffmap)

plot(H3K4me3_dpt, col_by='DPT1')

H3K4me3_dpt_df <- H3K4me3_dpt[,1:5]
colnames(H3K4me3_dpt_df) <- paste0('DPT', 1:5)
rownames(H3K4me3_dpt_df) <- colnames(H3K4me3_ctx)
H3K4me3_dpt_df <- as.data.frame(rank_cols(H3K4me3_dpt_df))

H3K4me3_ctx$dpt_tip <- H3K4me3_dpt@tips[,'Tips1']
H3K4me3_ctx[['diffmap']] <- CreateDimReducObject(
    H3K4me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K4me3_ctx <- AddMetaData(H3K4me3_ctx, H3K4me3_dc_df)
H3K4me3_ctx <- AddMetaData(H3K4me3_ctx, H3K4me3_dpt_df)


H3K4me3_ctx %>% 
    feature_plot(features=c('STMN2', 'VIM', 'DC1', 'DC2', 'DC3', 'DPT1', 'DPT2', 'DPT3'), order=T)

H3K4me3_ctx %>% write_rds('data/trajectories/ctx/H3K4me3_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_ctx_dpt_srt.rds')


#### Run diffmap for H3K27me3 cortex ####
H3K27me3_ctx <- marks_ctx$H3K27me3

H3K27me3_ctx <- H3K27me3_ctx %>% 
    FindTopFeatures(min.cutoff='q90') %>% 
    RunSVD() 

H3K27me3_ctx <- H3K27me3_ctx %>% 
    RunUMAP(dims=2:30, reduction='lsi')

DepthCor(H3K27me3_ctx)

H3K27me3_ctx %>% dim_plot(group.by=c('orig.ident', 'state', 'clusters'), label=T, order=T)
H3K27me3_ctx %>% feature_plot(features=c('STMN2', 'VIM', 'SOX2'), order=T, pt.size=2)

H3K27me3_pca <- H3K27me3_ctx[['lsi']]@cell.embeddings[,2:20]
H3K27me3_diffmap <- DiffusionMap(H3K27me3_pca, verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5]))

H3K27me3_dpt <- DPT(H3K27me3_diffmap)

plot(H3K27me3_dpt, col_by='DPT1')

H3K27me3_dpt_df <- H3K27me3_dpt[,1:5]
colnames(H3K27me3_dpt_df) <- paste0('DPT', 1:5)
rownames(H3K27me3_dpt_df) <- colnames(H3K27me3_ctx)
H3K27me3_dpt_df <- as.data.frame(rank_cols(H3K27me3_dpt_df))

H3K27me3_ctx$dpt_tip <- H3K27me3_dpt@tips[,'Tips1']
H3K27me3_ctx[['diffmap']] <- CreateDimReducObject(
    H3K27me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27me3_ctx <- AddMetaData(H3K27me3_ctx, H3K27me3_dc_df)
H3K27me3_ctx <- AddMetaData(H3K27me3_ctx, H3K27me3_dpt_df)


H3K27me3_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'DPT1', 'DPT2', 'DPT3'), order=T)

H3K27me3_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'DPT1', 'DPT2', 'DPT3'), reduction='diffmap', dims=c(1,2), order=T)

H3K27me3_ctx %>% write_rds('data/trajectories/ctx/H3K27me3_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_ctx_dpt_srt.rds')


#### Run diffmap for RNA cortex ####
rna_ctx <- rna_lin$ctx

rna_ctx <- rna_ctx %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

rna_ctx <- rna_ctx %>% 
    RunUMAP(dims=1:20) 

rna_ctx <- rna_ctx %>% 
    RunUMAP(dims=1:ncol(rna_ctx[['css']]), reduction.name='cssumap', reduction='css')

rna_ctx %>% dim_plot(group.by=c('orig.ident', 'state', 'clusters'), label=T, order=T, reduction='cssumap')
rna_ctx %>% feature_plot(features=c('NEUROD6', 'EMX1', 'SOX2', 'pseudotime_ranks'), order=T, reduction='cssumap')

rna_pca <- rna_ctx[['css']]@cell.embeddings
rna_diffmap <- DiffusionMap(rna_pca, verbose=T)
rna_dpt <- DPT(rna_diffmap)
rna_ctx[['diffmap']] <- CreateDimReducObject(
    rna_diffmap@eigenvectors*100, key='DC_', assay='RNA'
)
rna_dc_df <- as.data.frame(rank_cols(rna_diffmap@eigenvectors[,1:5]))
rna_dpt_df <- rna_dpt[,1:5]
colnames(rna_dpt_df) <- paste0('DPT', 1:5)
rownames(rna_dpt_df) <- colnames(rna_ctx)
rna_dpt_df <- as.data.frame(rank_cols(rna_dpt_df))

rna_ctx <- AddMetaData(rna_ctx, rna_dc_df)
rna_ctx <- AddMetaData(rna_ctx, rna_dpt_df)


rna_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'DPT1', 'DPT2', 'DPT3', 'pseudotime_ranks'), order=T, reduction='cssumap')

rna_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'DPT1', 'DPT2', 'DPT3', 'pseudotime_ranks'), reduction='diffmap', dims=c(1,2), order=T)

rna_ctx %>% write_rds('data/trajectories/ctx/RNA_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_ctx_dpt_srt.rds')


#### Align pseudotimes ####
H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_ctx_dpt_srt.rds')


H3K27ac_ctx$ctx_pt <- rank(-H3K27ac_ctx$DC2) / max(rank(-H3K27ac_ctx$DC2))
p1 <- H3K27ac_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

H3K4me3_ctx$ctx_pt <- rank(-H3K4me3_ctx$DC2) / max(rank(-H3K4me3_ctx$DC2))
p2 <- H3K4me3_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

H3K27me3_ctx$ctx_pt <- rank(H3K27me3_ctx$DC2) / max(rank(H3K27me3_ctx$DC2))
p3 <- H3K27me3_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

rna_ctx$ctx_pt <- rank(-rna_ctx$DC1) / max(rank(-rna_ctx$DC1))
p4 <- rna_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T, reduction='cssumap')

p1 / p2 / p3 / p4
ggsave('plots/trajectories/all_diff_pt_umap.png', width=8, height=16)



#### Plot gene expr/activity along pt to test ####
genes_plot <- c('VIM', 'SOX2', 'NEUROD6', 'NEUROD2', 'EOMES', 'EMX1')

H3K27ac_ctx %>% 
    feature_plot(features=genes_plot, order=T)

H3K27me3_ctx %>% 
    feature_plot(features=genes_plot, order=T)

H3K4me3_ctx %>% 
    feature_plot(features=genes_plot, order=T)


H3K27ac_expr <- H3K27ac_ctx[['cRNA']]@data[genes_plot, ] %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='cell', values_to='expr') %>% 
    inner_join(as_tibble(H3K27ac_ctx@meta.data, rownames='cell')) 

H3K27me3_expr <- H3K27me3_ctx[['cRNA']]@data[genes_plot, ] %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='cell', values_to='expr') %>% 
    inner_join(as_tibble(H3K27me3_ctx@meta.data, rownames='cell')) 

H3K4me3_expr <- H3K4me3_ctx[['cRNA']]@data[genes_plot, ] %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='cell', values_to='expr') %>% 
    inner_join(as_tibble(H3K4me3_ctx@meta.data, rownames='cell')) 

rna_expr <- rna_ctx[['RNA']]@data[genes_plot, ] %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='cell', values_to='expr') %>% 
    inner_join(as_tibble(rna_ctx@meta.data, rownames='cell')) 


plot_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr=expr/max(expr))

ggplot(plot_df, aes(ctx_pt, expr, color=modality)) +
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 3)) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(modality~gene, scales='free')


#### Make pseudotime bins ####

H3K27ac_ctx$pt_bins <- as.numeric(cut(H3K27ac_ctx$ctx_pt, 50, labels=1:50))
H3K27me3_ctx$pt_bins <- as.numeric(cut(H3K27me3_ctx$ctx_pt, 50, labels=1:50))
H3K4me3_ctx$pt_bins <- as.numeric(cut(H3K4me3_ctx$ctx_pt, 50, labels=1:50))
rna_ctx$pt_bins <- as.numeric(cut(rna_ctx$ctx_pt, 50, labels=1:50))


H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='cRNA', group_name='pt_bins')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='cRNA', group_name='pt_bins')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='cRNA', group_name='pt_bins')
rna_ctx <- Pando::aggregate_assay(rna_ctx, assay='RNA', group_name='pt_bins')



#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_ctx@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_clusters <- H3K27me3_ctx@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_clusters <- H3K4me3_ctx@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
rna_clusters <- rna_ctx@assays$RNA@misc$summary$pt_bins[as.character(1:50), ]

genes_plot <- c('STMN2', 'BCL11A', 'NEUROD6', 'NEUROD2')

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

ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'gp')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(~gene, scales='free')


#### GAMs for on pt_bin gene expression -> represent each timepoint by all modalities 
# -> multivariate DTW distance between GAM smooths
# -> K-means clustering -> identify lagging/parallel/antiparallel clusters of genes

mark_feats <- intersect(rownames(H3K27ac_ctx[['cRNA']]), rownames(H3K27me3_ctx[['cRNA']])) %>% intersect(rownames(H3K4me3_ctx[['cRNA']]))

rna_ctx <- FindVariableFeatures(rna_ctx, nfeatures=6000)
H3K27ac_ctx <- FindVariableFeatures(H3K27ac_ctx, assay='cRNA', nfeatures=6000)
H3K27me3_ctx <- FindVariableFeatures(H3K27me3_ctx, assay='cRNA', nfeatures=6000)
H3K4me3_ctx <- FindVariableFeatures(H3K4me3_ctx, assay='cRNA', nfeatures=6000)

all_var_feats <- bind_rows(
    'RNA' = as_tibble(rna_ctx[['RNA']]@meta.features, rownames='feature'),
    'H3K27ac' = as_tibble(H3K27ac_ctx[['cRNA']]@meta.features, rownames='feature'),
    'H3K27me3' = as_tibble(H3K27me3_ctx[['cRNA']]@meta.features, rownames='feature'),
    'H3K4me3' = as_tibble(H3K4me3_ctx[['cRNA']]@meta.features, rownames='feature'),
    .id = 'modality'
)

# feats_use <- VariableFeatures(rna_ctx) %>% 
#     intersect(mark_feats)

feats_use <- VariableFeatures(rna_ctx) %>% 
    intersect(VariableFeatures(H3K27ac_ctx, assay='cRNA')) %>% 
    intersect(VariableFeatures(H3K27me3_ctx, assay='cRNA')) %>% 
    intersect(VariableFeatures(H3K4me3_ctx, assay='cRNA')) 

#### Smooth with GAMs ####
x <- 1:50
H3K27ac_gams <- H3K27ac_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'gp'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'gp'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'gp'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'gp'))
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
# registerDoParallel(36)
# dtw_smooth_dist <- map_par(names(all_smooths_list), function(i){
#     map_dfr(names(all_smooths_list), function(j){
#         ij_align <- dtw(all_smooths_list[[i]], all_smooths_list[[j]])
#         ij_dist <- ij_align$distance
#         return(tibble(
#             i=i, j=j, dist=ij_dist
#         ))
#     })
# }, parallel=T)
# 
# dtw_smooth_dist_df <- bind_rows(dtw_smooth_dist)
# dtw_smooth_dist_mat <- dtw_smooth_dist_df %>% 
#     pivot_wider(names_from=j, values_from=dist) %>% 
#     column_to_rownames('i') %>% as.matrix()
# 
# pheatmap::pheatmap(dtw_smooth_dist_mat)

l2_smooth_dist <- map_par(names(all_smooths_list), function(i){
    map_dfr(names(all_smooths_list), function(j){
        ij_dist <- mean(diag(dist(all_smooths_list[[i]], all_smooths_list[[j]])))
        return(tibble(
            i=i, j=j, dist=ij_dist
        ))
    })
}, parallel=T)

l2_smooth_dist_df <- bind_rows(l2_smooth_dist)
l2_smooth_dist_mat <- l2_smooth_dist_df %>% 
    pivot_wider(names_from=j, values_from=dist) %>% 
    column_to_rownames('i') %>% as.matrix()

pheatmap::pheatmap(l2_smooth_dist_mat)


#### K-means clustering ####
library(FCPS)

# dtw_kmeans <- FCPS::kmeansClustering(
#     DataOrDistances = dtw_smooth_dist_mat,
#     ClusterNo = 10
# )

l2_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = l2_smooth_dist_mat,
    ClusterNo = 10
)


#### Plot clusters #####
var_feats_strict <- all_var_feats %>% 
    group_by(modality) %>% 
    top_n(4000, vst.variance.standardized) %>% 
    group_split() %>% map(~pull(.x, 'feature')) %>% 
    purrr::reduce(intersect)

all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
    {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
    as_tibble() %>% 
    mutate(pt_bin=1:50) %>% 
    return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr')

# dtw_kmeans_df <- dtw_kmeans$Cls %>% enframe('feature', 'dtw_clust')
l2_kmeans_df <- l2_kmeans$Cls %>% enframe('feature', 'l2_clust')
hclust_df <- l2_smooth_dist_mat %>% as.dist() %>% hclust() %>% cutree(k=10) %>% enframe('feature', 'h_clust')

plot_df <- all_gene_smooths_df %>% 
    # inner_join(dtw_kmeans_df) %>% 
    inner_join(l2_kmeans_df) %>% 
    inner_join(hclust_df) %>% 
    group_by(feature, modality) %>% 
    mutate(expr01=scale01(expr))
    
p1 <- ggplot(plot_df, aes(pt_bin, expr, color=modality, group=interaction(feature, modality))) +
    geom_line(alpha=0.5, size=0.5) +
    facet_wrap(~l2_clust) +
    no_legend() +
    ggtitle('euclidean')
p2 <- ggplot(plot_df, aes(pt_bin, expr, color=modality, group=interaction(feature, modality))) +
    geom_line(alpha=0.5, size=0.5) +
    facet_wrap(~dtw_clust) +
    ggtitle('DTW')
p1 | p2
ggsave('plots/trajectories/ctx_gam_smooth_clust_line.png', width=20, height=10)


ggplot(plot_df, aes(pt_bin, expr, color=modality, group=interaction(feature, modality))) +
    geom_line(alpha=0.5, size=0.5) +
    facet_grid(modality~l2_clust, scales='free') 
ggsave('plots/trajectories/ctx_gam_smooth_l2_clust_split_line.png', width=20, height=10)


p1 <- ggplot(plot_df, aes(pt_bin, expr, color=modality)) +
    geom_smooth(alpha=0.5, size=0.5) +
    facet_wrap(~l2_clust) +
    no_legend() +
    ggtitle('euclidean')
p2 <- ggplot(plot_df, aes(pt_bin, expr, color=modality)) +
    geom_smooth(alpha=0.5, size=0.5) +
    facet_wrap(~dtw_clust) +
    ggtitle('DTW')
p1 | p2
ggsave('plots/trajectories/ctx_gam_smooth_clust_smooth.png', width=20, height=10)


ggplot(filter(plot_df, feature%in%var_feats_strict), aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(h_clust~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100))
ggsave('plots/trajectories/ctx_gam_smooth_hclust_var_heatmap.png', width=20, height=30)







#### Match bin means by max corr ####
# RNA + H3K27ac
var_feats <- VariableFeatures(rna_ctx) %>% intersect(VariableFeatures(H3K27ac_ctx, assay='cRNA'))

rna_H3K27ac_cor <- Pando::sparse_cor(t(rna_clusters[,var_feats]), t(H3K27ac_clusters[,var_feats]))
pheatmap::pheatmap(t(rna_H3K27ac_cor), cluster_rows=F, cluster_cols=F)

rna_H3K27ac_align <- dtw(
    x = as.matrix(1-rna_H3K27ac_cor),
    keep=T
)
plot(rna_H3K27ac_align)


# RNA + H3K4me3
var_feats <- VariableFeatures(rna_ctx) %>% intersect(VariableFeatures(H3K4me3_ctx, assay='cRNA'))

rna_H3K4me3_cor <- Pando::sparse_cor(t(rna_clusters[,var_feats]), t(H3K4me3_clusters[,var_feats]))
pheatmap::pheatmap(t(rna_H3K4me3_cor), cluster_rows=F, cluster_cols=F)

rna_H3K4me3_align <- dtw(
    x = as.matrix(1-rna_H3K4me3_cor)
)
plot(rna_H3K4me3_align)


# RNA + H3K27me3
var_feats <- VariableFeatures(rna_ctx) %>% intersect(VariableFeatures(H3K27me3_ctx, assay='cRNA'))
var_feats <- VariableFeatures(rna_ctx) %>% intersect(rownames(H3K27me3_ctx[['cRNA']]))
var_feats <- VariableFeatures(H3K27me3_ctx[['cRNA']]) %>% intersect(rownames(rna_ctx[['RNA']]))

rna_H3K27me3_cor <- Pando::sparse_cor(t(rna_clusters[,var_feats]), t(H3K27me3_clusters[,var_feats]))
pheatmap::pheatmap(t(rna_H3K27me3_cor), cluster_rows=F, cluster_cols=F)

rna_H3K27me3_align <- dtw(
    x = as.matrix(rna_H3K27me3_cor)
)
plot(rna_H3K27me3_align)






