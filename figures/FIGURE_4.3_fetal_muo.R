library(tidyverse)
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/de.R')

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

select <- dplyr::select
Matrix <- Matrix::Matrix

setwd('~/projects/cutntag/')

rank_cols <- function(x){
    xr <- presto::rank_matrix(t(x))$X_ranked
    colnames(xr) <- colnames(x)
    rownames(xr) <- rownames(x)
    return(xr)
}




atlas_celltype_colors <- c(
    'neural_progenitor_cell' = '#c38b79',
    'telencephalic_npc' = '#b77f92',
    'dorsal_npc' = '#D19C2C',
    'dorsal_org' = '#e76f51',
    'ventral_npc' = '#ce93d8',
    'ventral_org' = '#9575cd',
    'dien_mesen_npc' = '#80deea',
    'diencephalic_npc' = '#29b6f6',
    'mesencephalic_npc' = '#26a69a',
    'retinal_npc' = '#ffccbc',
    'rhombencephalic_npc' = '#a5d6a7',
    'intermediate_progenitor' = '#ffe082',
    'neuron' = '#264653',
    'cortical_excitatory_neuron' = '#b44a6e',
    'deeper_layer_cortical_excitatory_neuron' = '#ec407a',
    'upper_layer_cortical_excitatory_neuron' =  '#8B1721',
    'non_cortical_excitatory_neuron' = '#f9a825',
    'diencephalic_excitatory_neuron' = '#0277bd',
    'mesencephalic_excitatory_neuron' = '#00897b',
    'rhombencephalic_excitatory_neuron' = '#558b2f',
    'inhibitory_neuron' = '#039be5',
    'ventral_inhibitory_neuron' = '#76448a',
    'mge_inhibitory_neuron' = '#8e44ad' ,
    'lge_inhibitory_neuron' = '#5c6bc0',
    'cge_inhibitory_neuron' = '#2B254B',
    'diencephalic_inhibitory_neuron' = '#2196f3',
    'mesencephalic_inhibitory_neuron' = '#00acc1',
    'rhombencephalic_inhibitory_neuron' = '#4caf50',
    'chrorid_plexus_epithelium' = '#808b96',
    'astrocyte' = '#f06292',
    'oligodendrocyte_lineage' = '#00bf99',
    'oligodendrocyte_precursor_cell' = '#75e4e8',
    'mature_oligodendrocyte' = '#00604d' ,
    'microglia' = '#b2ff00',
    'vascular_endothelial_cell' = '#F20000',
    'mesenchymal_cell' = '#795548',
    'neural_crest' = '#dce775',
    'pns_neurons' = '#9e9d24',
    'immune_cell' = '#a791aa'
)




fetal_colors <- c("mesencephalic_npc"="#8BAED3",
                "dorsal_npc"="#F2A1BD", 
                "cortical_excitatory_neuron"="#D25B78", 
                'intermediate_progenitor' = "#E2C1BE",
                "astrocytes"="#A9DFBF",
                "mesenchymal_cell"="#949398", 
                "choroid_plexus"="#A1C8C9", 
                "nect"="#8A5796",
                "mesen_ex"="#5961A1", 
                "rhom_ex"="#00569B", 
                "diencephalic_npc"="#049DB1", 
                "dien_ex"="#85A0AB",
                "dorsal_org"='#e76f51',
                "oligodendrocyte_precursor_cell"="#46996C", 
                "non_nect"="#D2927D", 
                "psc"="#C6AAAF", 
                "RGC"="#ECC270",
                "RPC"="#F5DF4D", 
                'immune_cell' = '#a791aa',
                'vascular_endothelial_cell' = '#F20000',
                'ventral_inhibitory_neuron' = '#b20689',
                'neuron' = '#264653',
                'neural_progenitor_cell' = '#f3dea1'
)



#### Read Data ####

annotation <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')

data_paths <- list(
    'MUO_w19' = '/nas/groups/treutlein/DATA/sequencing/20230215_P2296_FIDES_HSIU-CHUAN_ATAC/processed/fetal_rep1/outs/',
    'MUO_w10' = '/local1/DATA/sequencing/20230310_P2324_FIDES_ATAC/processed/fetal_w10_rep1/outs/',
    'MUO_w11' = '/local1/DATA/sequencing/20230310_P2324_FIDES_ATAC/processed/fetal_w11_rep1/outs/'
)

count_list <- map(data_paths, function(x){
    Read10X_h5(paste0(x, 'filtered_feature_bc_matrix.h5'))
})

srt_list <- map(names(count_list), function(n){
    print(n)
    
    x <- count_list[[n]]
    path <- data_paths[[n]]
    
    fragpath <- paste0(path, 'atac_fragments.tsv.gz')
    
    srt <- CreateSeuratObject(
        counts = x$`Gene Expression`,
        assay = "RNA"
    )
    
    srt[["ATAC"]] <- CreateChromatinAssay(
        counts = x$Peaks,
        sep = c(":", "-"),
        fragments = fragpath,
        annotation = annotation
    )
    
    return(srt)
})


rna_srt_list <- srt_list %>% map(~DietSeurat(.x, assays='RNA'))
names(rna_srt_list) <- names(count_list)
rna_srt <- merge(rna_srt_list[[1]], rna_srt_list[2:length(rna_srt_list)], add.cell.ids=names(data_paths))

muo_srt <- merge(srt_list[[1]], srt_list[2:length(srt_list)], add.cell.ids=names(srt_list))
muo_srt <- RenameCells(muo_srt, new.names=colnames(rna_srt))

muo_srt@active.assay <- 'RNA'
muo_srt <- muo_srt %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

ElbowPlot(muo_srt)

muo_srt <- muo_srt %>% 
    RunUMAP(dims=1:15)

muo_srt$log_rna_counts <- log10(muo_srt$nCount_RNA)
muo_srt$log_atac_counts <- log10(muo_srt$nCount_ATAC)

muo_srt %>% write_rds('data/fetal/MUO/MUO_d10_d11_d19_merged_v0_srt.rds')

muo_srt$orig.ident <- colnames(muo_srt) %>% str_replace('(.+)_[ATCG]+-\\d', '\\1')
muo_srt$age <- colnames(muo_srt) %>% str_replace('MUO_(.+)_[ATCG]+-\\d', '\\1')

muo_srt <- muo_srt %>% FindNeighbors()
muo_srt <- muo_srt %>% FindClusters()

muo_srt %>% dim_plot(group.by=c('age', 'seurat_clusters'))
muo_srt %>% feature_plot(features=c('NEUROD6', 'SOX2', 'EOMES', 'OTX2', 'HOXB2', 'RSPO3', 'TCF7L2', 'STMN2', 'GRIA2'), order=T)


#### Some QC ####
meta <- muo_srt@meta.data %>% 
    as_tibble(rownames='cell')

p1 <- ggplot(meta, aes(log_rna_counts)) +
    geom_histogram(color='black', fill='grey') +
    facet_grid(orig.ident~.) +
    ggtitle('RNA counts')

p2 <- ggplot(meta, aes(nFeature_RNA)) +
    geom_histogram(color='black', fill='grey') +
    scale_x_continuous(breaks=seq(0,20000,1000)) +
    facet_grid(orig.ident~.) +
    ggtitle('RNA features')

p3 <- ggplot(meta, aes(log_atac_counts)) +
    geom_histogram(color='black', fill='grey') +
    facet_grid(orig.ident~.) +
    ggtitle('ATAC counts')

p4 <- ggplot(meta, aes(nFeature_ATAC)) +
    geom_histogram(color='black', fill='grey') +
    scale_x_continuous(breaks=seq(0,200000,20000)) +
    facet_grid(orig.ident~.) +
    ggtitle('ATAC features')

(p1 / p2) | (p3 / p4)

ggsave('plots/fetal/init_qc_hist.png', width=12, height=8)




#### Apply filters ####
muo_srt_fltr <- muo_srt %>% subset(
    nFeature_RNA > 1000 &
        nFeature_RNA < 6000 &
        nCount_RNA >1500 &
        nFeature_ATAC < 50000 &
        nFeature_ATAC > 1000
)

meta <- muo_srt_fltr@meta.data %>% 
    as_tibble(rownames='cell')

p1 <- ggplot(meta, aes(log_rna_counts)) +
    geom_histogram(color='black', fill='grey') +
    facet_grid(orig.ident~.) +
    ggtitle('RNA counts')

p2 <- ggplot(meta, aes(nFeature_RNA)) +
    geom_histogram(color='black', fill='grey') +
    scale_x_continuous(breaks=seq(0,20000,1000)) +
    facet_grid(orig.ident~.) +
    ggtitle('RNA counts')

p3 <- ggplot(meta, aes(log_atac_counts)) +
    geom_histogram(color='black', fill='grey') +
    facet_grid(orig.ident~.) +
    ggtitle('RNA counts')

p4 <- ggplot(meta, aes(nFeature_ATAC)) +
    geom_histogram(color='black', fill='grey') +
    scale_x_continuous(breaks=seq(0,200000,20000)) +
    facet_grid(orig.ident~.) +
    ggtitle('RNA counts')

(p1 / p2) | (p3 / p4)

ggsave('plots/fetal/filtered_qc_hist.png', width=12, height=8)


muo_srt_fltr <- muo_srt_fltr %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

muo_srt_fltr <- muo_srt_fltr %>% RunUMAP(dims=1:20)
muo_srt_fltr <- muo_srt_fltr %>% FindNeighbors()
muo_srt_fltr <- muo_srt_fltr %>% FindClusters()

p1 <- muo_srt_fltr %>% dim_plot(group.by=c('age', 'seurat_clusters'), label=T)
p2 <- muo_srt_fltr %>% feature_plot(features=c('NEUROD6', 'SOX2', 'EOMES', 'OTX2', 'HOXB2', 'RSPO3', 'TCF7L2', 'STMN2', 'GRIA2', 'FOXG1'), order=T)

(p1 / p2) + plot_layout(heights=c(1,3))
ggsave('plots/fetal/filtered_marker_umap.png', width=8, height=8, bg='white')

muo_srt_fltr %>% write_rds('data/fetal/MUO/MUO_d10_d11_d19_merged_v1filtered_srt.rds')
muo_srt_fltr <- read_rds('data/fetal/MUO/MUO_d10_d11_d19_merged_v1filtered_srt.rds')

gene_activities <- GeneActivity(
    muo_srt_fltr, 
    assay = 'ATAC',
    features = rownames(muo_srt_fltr[['RNA']])
)

muo_srt_fltr[['gene_activity']] <- CreateAssayObject(gene_activities)

muo_srt_fltr@active.assay <- 'gene_activity'
muo_srt_fltr <- muo_srt_fltr %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()


goi <- c('NEUROD6', 'SOX2', 'EOMES', 'OTX2', 'HOXB2', 'RSPO3', 'TCF7L2', 'STMN2', 'GRIA2', 'FOXG1', 'VIM', 'NEUROD2')
muo_srt_fltr@active.assay <- 'gene_activity'
p1 <- muo_srt_fltr %>% feature_plot(features=goi, order=T, reduction='umap') &
    scale_color_gradientn(colors=gyylgnbu())
muo_srt_fltr@active.assay <- 'RNA'
p2 <- muo_srt_fltr %>% feature_plot(features=goi, order=T, reduction='umap')

p1 / p2
ggsave('plots/fetal/rna_ga_marker_umap.png', width=8, height=16, bg='white')

muo_srt_fltr %>% write_rds('data/fetal/MUO/MUO_d10_d11_d19_merged_v1.1ga_srt.rds')






#### Extract ctx trajectory ####
muo_srt_fltr <- read_rds('data_/fetal/MUO/MUO_d10_d11_d19_merged_v1.1ga_srt.rds')

library(destiny)

muo_srt_fltr@active.assay <- 'RNA'
muo_srt_fltr %>% dim_plot(reduction='umap', label=T)
muo_srt_fltr %>% feature_plot(features=c('NEUROD6', 'STMN2', 'NEUROD2', 'MKI67'), order=T, reduction='umap')

muo_ctx <- muo_srt_fltr %>% subset(seurat_clusters%in%c(10,12) & age=='w19')

muo_ctx@active.assay <- 'RNA'
muo_ctx <- muo_ctx %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

muo_ctx <- muo_ctx %>% 
    RunUMAP(dims=1:20)

muo_ctx %>% dim_plot(group.by=c('orig.ident', 'state', 'seurat_clusters'), label=T, order=T)

muo_ctx@active.assay <- 'gene_activity'
muo_ctx %>% feature_plot(features=c('NEUROD6', 'STMN2', 'NEUROD2', 'MKI67', 'GRIA2', 'EOMES'), order=T)
muo_ctx@active.assay <- 'RNA'
muo_ctx %>% feature_plot(features=c('NEUROD6', 'STMN2', 'NEUROD2', 'MKI67', 'GRIA2', 'EOMES'), order=T)


muo_pca <- muo_ctx[['pca']]@cell.embeddings[,1:10]
muo_diffmap <- DiffusionMap(muo_pca, verbose=T)
muo_dpt <- DPT(muo_diffmap)
muo_ctx[['diffmap']] <- CreateDimReducObject(
    muo_diffmap@eigenvectors*100, key='DC_', assay='RNA'
)
muo_dc_df <- as.data.frame(rank_cols(muo_diffmap@eigenvectors[,1:5]))
muo_dpt_df <- muo_dpt[,1:5]
colnames(muo_dpt_df) <- paste0('DPT', 1:5)
rownames(muo_dpt_df) <- colnames(muo_ctx)
muo_dpt_df <- as.data.frame(rank_cols(muo_dpt_df))

muo_ctx <- AddMetaData(muo_ctx, muo_dpt_df)
muo_ctx <- AddMetaData(muo_ctx, muo_dc_df)

muo_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'DPT1', 'DPT2', 'DPT3'), order=T, reduction='umap')

muo_ctx$ctx_pt <- rank(muo_ctx$DC1) / max(rank(muo_ctx$DC1))
muo_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T, reduction='umap')



muo_ctx %>% 
    feature_plot(features=c('nCount_ATAC', 'log_atac_counts'), order=T, reduction='umap')

muo_ctx %>% dim_plot(group.by=c('orig.ident'), label=T, order=T)




muo_ctx %>% dim_plot(group.by=c('celltype_jf'), label=T, order=T) +
    scale_color_manual(values=fetal_colors) + no_legend()
ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_celltype_umap.png', width=7, height=3)
ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_celltype_umap.pdf', width=7, height=3)


muo_ctx %>% FeaturePlot(features=c('ctx_pt'), order=T, reduction='umap') +
    scale_color_gradientn(colors=rev(pals::magma(100)))
ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_pseudotime_umap.png', width=7, height=3)
ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_pseudotime_umap.pdf', width=7, height=3)


muo_ctx@active.assay <- 'RNA'
p1 <- muo_ctx %>% feature_plot(features=c('NEUROD6', 'NEUROD2', 'GRIA2', 'STMN2', 'SOX2', 'VIM', 'HOPX'), order=T, reduction='umap') &
    scale_color_gradientn(colors=gyorrd())

gypu <- function(bias=1){
    return(colorRampPalette(c('#e5e7e9', brewer.pal(n=9, name='BuPu')), bias=bias)(100))
}
muo_ctx@active.assay <- 'gene_activity'
p2 <- muo_ctx %>% feature_plot(features=c('NEUROD6', 'NEUROD2', 'GRIA2', 'STMN2', 'SOX2', 'VIM', 'HOPX'), order=T, reduction='umap') &
    scale_color_gradientn(colors=gypu(1.3))

p1 / p2
ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_ga_rna_feature_umap.png', width=10, height=12)


muo_ctx %>% write_rds('data/fetal/MUO/MUO_w19.rds')



#### Make pseudotime bins ####
muo_ctx$pt_bins <- as.numeric(cut(muo_ctx$ctx_pt, 20, labels=1:20))

muo_ctx <- Pando::aggregate_assay(muo_ctx, assay='gene_activity', group_name='pt_bins')
muo_ctx <- Pando::aggregate_assay(muo_ctx, assay='RNA', group_name='pt_bins')

muo_ctx[['gene_activity_bin']] <- CreateAssayObject((muo_ctx[['gene_activity']]@data > 0)*1)
muo_ctx[['RNA_bin']] <- CreateAssayObject((muo_ctx[['RNA']]@data > 0)*1)

muo_ctx <- Pando::aggregate_assay(muo_ctx, assay='gene_activity_bin', group_name='pt_bins')
muo_ctx <- Pando::aggregate_assay(muo_ctx, assay='RNA_bin', group_name='pt_bins')


#### Plot expression over pt bins ####
atac_clusters <- muo_ctx@assays$gene_activity_bin@misc$summary$pt_bins[as.character(1:20), ]
rna_clusters <- muo_ctx@assays$RNA_bin@misc$summary$pt_bins[as.character(1:20), ]

genes_plot <- c('STMN2', 'NEUROD6', 'NEUROD2', 'BCL11B')

atac_expr <- atac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

plot_df <- bind_rows('ATAC'=atac_expr, 'RNA'=rna_expr, .id='modality') %>% 
    # dplyr::filter(gene=='NEUROD2') %>% 
    dplyr::group_by(modality, gene) %>% 
    dplyr::mutate(expr01=scale01(expr))


ggplot(plot_df, aes(as.numeric(pt_bins), expr01, color=modality)) +
    geom_bar(stat='identity', fill='grey', alpha=0.6, position='dodge', linewidth=0.2) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'gp'), linewidth=0.5) +
    scale_color_manual(values=c('ATAC'='#9575cd', 'RNA'='#FDA044')) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_grid(~gene, scales='free') +
    article_text() +
    labs(x='Pseudotime bins', y='Scaled detection rate')

ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_pt_rna_atac_smooth.png', width=16, height=4, units='cm')
ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_pt_rna_atac_smooth.pdf', width=16, height=4, units='cm')




genes_plot <- c('GRIA2', 'NEUROD6', 'NEUROD2', 'BCL11B', 'VIM', 'SOX2', 'STMN2', 'HOPX', 'GLI3', 'EOMES', 'HOPX', 'NES')

atac_expr <- atac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

plot_df <- bind_rows('ATAC'=atac_expr, 'RNA'=rna_expr, .id='modality') %>% 
    # dplyr::filter(gene=='NEUROD2') %>% 
    dplyr::group_by(modality, gene) %>% 
    dplyr::mutate(expr01=scale01(expr))


ggplot(plot_df, aes(as.numeric(pt_bins), expr01, color=modality)) +
    geom_bar(stat='identity', fill='grey', alpha=0.6, position='dodge', linewidth=0.2) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'gp'), linewidth=0.5) +
    scale_color_manual(values=c('ATAC'='#9575cd', 'RNA'='#FDA044')) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_wrap(~gene, scales='free') +
    article_text() +
    labs(x='Pseudotime bins', y='Scaled detection rate')

ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_pt_rna_atac_more_smooth.png', width=16, height=10, units='cm')
ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_pt_rna_atac_more_smooth.pdf', width=16, height=10, units='cm')


muo_ctx@active.assay <- 'RNA'
p1 <- feature_plot(muo_ctx, features=genes_plot, order=T) &
    scale_color_gradientn(colors=gyorrd())
muo_ctx@active.assay <- 'gene_activity'
p2 <- feature_plot(muo_ctx, features=genes_plot, order=T) &
    scale_color_gradientn(colors=gypu(1.3))
p1 / p2

ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_pt_rna_atac_more_umap.png', bg='white', width=15, height=15)
ggsave('plots/paper/fig3/fig3_fetal_w19_ctx_pt_rna_atac_more_umap.pdf', bg='white', width=15, height=15)





##### Get pt smooths for features ####
gene_clusters <- read_tsv('data/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')

atac_cluster_detect <- muo_ctx@assays$gene_activity_bin@misc$summary$pt_bins[as.character(1:20), ]
rna_cluster_detect <- muo_ctx@assays$RNA_bin@misc$summary$pt_bins[as.character(1:20), ]

feats_use <- gene_clusters %>% pull(feature) %>% unique() %>% intersect(colnames(atac_cluster_detect))


#### Smooth with GAMs ####
x <- 1:20
atac_gams <- atac_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
atac_smooths <- atac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
rna_smooths <- rna_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

atac_smooths_scale <- atac_smooths %>% scale()
rna_smooths_scale <- rna_smooths %>% scale()

all_smooths_list <- map(purrr::set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            atac=atac_smooths_scale[,n],
            rna=rna_smooths_scale[,n]
        )
    )
})

all_smooths_raw_list <- map(purrr::set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            atac=atac_smooths[,n],
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

gene_order <- dtw_smooth_dist_mat %>% as.dist() %>% hclust() %>% {.$labels[.$order]}



all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
        {colnames(.) <- c('ATAC', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=1:20) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(ATAC:RNA, names_to='modality', values_to='expr')


plot_df <- all_gene_smooths_df %>% 
    # inner_join(cluster_hclust_clusts) %>%
    dplyr::group_by(feature, modality) %>% 
    dplyr::mutate(
        expr01=scale01(expr), 
        pt_bin_mod=paste0(modality, pt_bin),
        feature=factor(feature, levels=purrr::reduce(gene_order, c)),
        # dtw_clust=factor(dtw_clust, levels=c(7,1,6,10,4,9,8,2,3,5))
    )

ggplot(plot_df, aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100)) +
    theme(
        panel.border = element_blank()
    )




#### Annotate regions and check switching/bimodal peaks ####
cluster_de <- de(muo_srt_fltr, groups = 'seurat_clusters')

cluster_de %>% filter(padj<1e-4, fc>0.5) %>% group_by(group) %>% top_n(10, fc) %>% View()

cluster_markers <- cluster_de %>% 
    filter(padj<1e-4, fc>0.5) %>% 
    group_by(group) %>% top_n(3, fc) %>% 
    pull(feature) %>% unique()

muo_srt_fltr %>% feature_plot(features=c('DLX5'), order=T)
muo_srt_fltr %>% feature_plot(features=cluster_markers, order=T)
ggsave('plots/fetal/filtered_cluster_marker_umap.png', width=20, height=40, bg='white')


# 5 - vascular
# 7 - RG?
# 8 - mesenchyme/fibroblast
# 13 - immune

muo_srt_fltr$celltype_jf <- case_when(
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(5) ~ 'vascular_endothelial_cell',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(3,16,11) ~ 'neural_progenitor_cell',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(8) ~ 'mesenchymal_cell',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(13) ~ 'immune_cell',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(17) ~ 'mesenchymal_cell',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(2) ~ 'mesencephalic_npc',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(1,15) ~ 'diencephalic_npc',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(6,19,4) ~ 'neuron',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(0) ~ 'rhombencephalic_npc',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(10) ~ 'dorsal_npc',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(18) ~ 'dorsal_org',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(12) ~ 'cortical_excitatory_neuron',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(20) ~ 'ventral_inhibitory_neuron',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(9) ~ 'oligodendrocyte_precursor_cell',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(14) ~ 'intermediate_progenitor',
    muo_srt_fltr$RNA_snn_res.0.8 %in% c(14) ~ 'intermediate_progenitor',
)


p1 <- muo_srt_fltr %>% dim_plot(group.by=c('celltype_jf')) +
    scale_color_manual(values=fetal_colors)
p2 <- muo_srt_fltr %>% dim_plot(group.by=c('age'))

p1 | p2




#### Extract w19 ####
muo_srt_fltr <- read_rds('data_/fetal/MUO/MUO_d10_d11_d19_merged_v1.1ga_srt.rds')
muo_w19 <- muo_srt_fltr %>% subset(age=='w19')

muo_w19@active.assay <- 'RNA'
muo_w19 <- muo_w19 %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

muo_w19 <- muo_w19 %>% 
    FindNeighbors(dims=1:10) %>% 
    FindClusters()

muo_w19 <- muo_w19 %>% RunUMAP(dims=1:10)

muo_w19 %>% dim_plot(group.by=c('celltype_jf', 'seurat_clusters'), label=T) 

muo_w19$celltype_jf <- case_when(
    muo_w19$seurat_clusters %in% c(5,3) ~ 'cortical_excitatory_neuron',
    muo_w19$seurat_clusters %in% c(4) ~ 'intermediate_progenitor',
    muo_w19$seurat_clusters %in% c(1,8,2) ~ 'dorsal_npc',
    muo_w19$seurat_clusters %in% c(0) ~ 'dorsal_org',
    muo_w19$seurat_clusters %in% c(10) ~ 'ventral_inhibitory_neuron',
    muo_w19$seurat_clusters %in% c(9,13,14) ~ 'vascular_endothelial_cell',
    muo_w19$seurat_clusters %in% c(7,15) ~ 'mesenchymal_cell',
    muo_w19$seurat_clusters %in% c(6) ~ 'immune_cell',
    muo_w19$seurat_clusters %in% c(11) ~ 'oligodendrocyte_precursor_cell',
    muo_w19$seurat_clusters %in% c(12) ~ 'astrocytes',
)


muo_w19 %>% dim_plot(group.by=c('celltype_jf')) +
    scale_color_manual(values=fetal_colors)
ggsave('plots/paper/fig3/fetal_w19_celltype_umap.png', bg='white', width=8, height=5)
ggsave('plots/paper/fig3/fetal_w19_celltype_umap.pdf', bg='white', width=8, height=5)


feature_plot(muo_w19, features=c("FOXG1", "NEUROD2", "STMN2", "EMX1", "NFIB", "AQP4", "BCL11B", "EMX1", 'SATB2'), order=T) &
    scale_color_gradientn(colors=gyorrd())
ggsave('plots/paper/fig3/fig3_w19_marker_umap.png', bg='white', width=15, height=15)

feature_plot(muo_w19, features=c('FOXG1', 'SOX2', 'HOPX', 'NEUROD2', 'STMN2', 'EMX1', 'NFIB', 'LMX1A', 'VSX2', 'OLIG1', 'DLX2', 'DCN', 'CLDN5', 'CD53', 'MKI67', 'EOMES'), order=T) &
    scale_color_gradientn(colors=gyorrd())
ggsave('plots/paper/fig3/fig3_w19_marker_umap.png', bg='white', width=15, height=15)



plot_df <- muo_w19@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(plot_df, aes('1', fill=celltype_jf)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=fetal_colors) +
    article_text() +
    no_x_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs('Fraction of cells')
ggsave('plots/paper/fig3/fig3_fetal_w19_celltype_bar.png', bg='white', width=6, height=5, units='cm')
ggsave('plots/paper/fig3/fig3_fetal_w19_celltype_bar.pdf', bg='white', width=6, height=5, units='cm')


muo_w19 %>% write_rds('data_/fetal/MUO/MUO_w19.rds')
muo_w19 <- read_rds('data_/fetal/MUO/MUO_w19.rds')



#### Subcluster neurons ####
muo_neurons <- muo_srt_fltr %>% 
    subset(celltype_jf %in% c('cortical_excitatory_neuron', 'ventral_inhibitory_neuron', 'neuron'))

muo_neurons <- muo_neurons %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

muo_neurons <- muo_neurons %>% FindNeighbors(dims=1:10)
muo_neurons <- muo_neurons %>% FindClusters()
muo_neurons <- muo_neurons %>% RunUMAP(dims=1:10)


muo_neurons$celltype_neuron_jf <- dplyr::case_when(
    muo_neurons$RNA_snn_res.0.8 %in% c(2,0) ~ 'diencephalic_inhibitory_neuron',
    muo_neurons$RNA_snn_res.0.8 %in% c(12) ~ 'diencephalic_excitatory_neuron',
    muo_neurons$RNA_snn_res.0.8 %in% c(8) ~ 'mesencephalic_excitatory_neuron',
    muo_neurons$RNA_snn_res.0.8 %in% c(1) ~ 'mesencephalic_inhibitory_neuron',
)

muo_neurons$celltype_neuron_jf <- dplyr::case_when(
    is.na(muo_neurons$celltype_neuron_jf) ~ muo_neurons$celltype_jf,
    T ~ muo_neurons$celltype_neuron_jf
)

p1 <- muo_neurons %>% dim_plot(group.by=c('celltype_neuron_jf')) +
    scale_color_manual(values=atlas_celltype_colors)
p2 <- muo_neurons %>% dim_plot(group.by=c('seurat_clusters'), label=T)
p_annot <- p1 | p2


goi <- c('NEUROD6', 'SOX2', 'EOMES', 'OTX2', 'HOXB2', 'RSPO3', 'TCF7L2', 'STMN2', 'GRIA2', 'FOXG1', 'NEUROD2', 'GAD2', 'SLC17A6', 'LHX5', 'LHX9')
# muo_neurons@active.assay <- 'gene_activity'
# p1 <- muo_neurons %>% feature_plot(features=goi, order=T, reduction='umap') &
#     scale_color_gradientn(colors=gyylgnbu())
muo_neurons@active.assay <- 'RNA'
p2 <- muo_neurons %>% feature_plot(features=goi, order=T, reduction='umap')
p2

(p_annot / p2) + plot_layout(heights=c(1,2))
ggsave('plots/fetal/neurons_annot_umap.png', bg='white', width=12, height=16)

muo_neurons %>% dim_plot(group.by=c('age'), label=T)



#### Link peaks ####
muo_srt_fltr <- FindVariableFeatures(muo_srt_fltr, nfeatures=8000)
var_feats <- VariableFeatures(muo_srt_fltr)

muo_srt_fltr <- RegionStats(muo_srt_fltr, assay = 'ATAC', genome = BSgenome.Hsapiens.UCSC.hg38)

muo_srt_fltr <- LinkPeaks(
    muo_srt_fltr,
    peak.assay = 'ATAC', 
    expression.assay = 'RNA', 
    distance = 100000,
    genes.use = var_feats
)

muo_srt_fltr %>% write_rds('data/fetal/MUO/MUO_d10_d11_d19_merged_v1.2links_srt.rds')



#### Apply filters on ATAC only ####
muo_srt_atac_fltr <- muo_srt %>% subset(
    nFeature_ATAC < 50000 &
        nFeature_ATAC > 1000
)

muo_srt_atac_fltr@active.assay <- 'ATAC'
muo_srt_atac_fltr <- muo_srt_atac_fltr %>% 
    FindTopFeatures() %>% 
    RunTFIDF() %>% 
    RunSVD()
    
ElbowPlot(muo_srt_atac_fltr)
DepthCor(muo_srt_atac_fltr)
    
muo_srt_atac_fltr <- muo_srt_atac_fltr %>% 
    RunUMAP(dims=2:20)


p1 <- muo_srt_atac_fltr %>% dim_plot(group.by=c('age', 'seurat_clusters'), split.by='orig.ident')
muo_srt_atac_fltr@active.assay <- 'RNA'
p2 <- muo_srt_atac_fltr %>% feature_plot(features=goi, order=T, reduction='umap')
muo_srt_atac_fltr@active.assay <- 'gene_activity'
p3 <- muo_srt_atac_fltr %>% feature_plot(features=goi, order=T, reduction='umap') &
    scale_color_gradientn(colors=gyylgnbu())
p3
(p1 / p2) + plot_layout(heights=c(1,3))
ggsave('plots/fetal/atac_filtered_marker_umap.png', width=8, height=8, bg='white')


gene_activities <- GeneActivity(
    muo_srt_atac_fltr, 
    assay = 'ATAC',
    features = rownames(muo_srt_fltr[['RNA']])
)

muo_srt_atac_fltr[['gene_activity']] <- CreateAssayObject(gene_activities)

muo_srt_fltr %>% write_rds('data/fetal/MUO/MUO_d10_d11_d19_merged_v1.1.2atacfilter_srt.rds')



#### Intersection of fetal peaks with organoid peaks ####
muo_w19 <- read_rds('data_/fetal/MUO/MUO_w19.rds')
all_peaks_srt <- read_rds('data05/regions/all_regions_v2_max01_clust50_srt.rds')

fetal_detect <- rowMeans(muo_w19[['ATAC']]@data > 0)
fetal_peaks <- names(fetal_detect)[fetal_detect>0.05]

fet_atac_ct_olap <- findOverlaps(StringToGRanges(all_peaks_srt$repr_region), StringToGRanges(fetal_peaks))

atac_ct_tbl <- tibble(
    region = fetal_peaks,
    in_organoids = fetal_peaks %in% fetal_peaks[subjectHits(fet_atac_ct_olap)]
)

p2 <- ggplot(atac_ct_tbl, aes('fetalATAC_in_orgCT', fill=in_organoids)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_y_continuous(limits=c(0,50000)) +
    labs(x='') + no_legend() + 
    scale_axis_rangeframe() + theme_rangeframe()

ct_atac_tbl <- tibble(
    region = all_peaks_srt$repr_region,
    in_organoids = all_peaks_srt$repr_region %in%  all_peaks_srt$repr_region[queryHits(fet_atac_ct_olap)]
)

p1 <- ggplot(ct_atac_tbl, aes('orgCT_in_fetalATAC', fill=in_organoids)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_y_continuous(limits=c(0,50000)) +
    labs(x='') + no_legend() + 
    scale_axis_rangeframe() + theme_rangeframe()

p1 | p2
ggsave('plots/paper/CT_compare_fetal_ATAC_intersection_bar.pdf', width=3, height=2)




H3K4me3_bw1 <- import.bw('/nas/groups/treutlein/DATA/sequencing/20230310_P2324_FIDES_ATAC/processed_bulk/mapped/bamCoverage/230117_15_FZ_bCT_fb-w19_H3K4me3_S1.filtered.seq_depth_norm.bw')
H3K4me3_bw2 <- import.bw('/nas/groups/treutlein/DATA/sequencing/20230525_P2426_FIDES/mapped/bamCoverage/23011709_25_FZ_bCT_fb-w19_H3K4me3_S1.filtered.seq_depth_norm.bw')

H3K4me3_bw1_fltr <- H3K4me3_bw1[H3K4me3_bw1$score>20] %>% {intersect(.,.)} %>% {.[width(.)>200]}
H3K4me3_bw2_fltr <- H3K4me3_bw2[H3K4me3_bw2$score>20] %>% {intersect(.,.)} %>% {.[width(.)>200]}

H3K4me3_fetal_peaks <- intersect(H3K4me3_bw1_fltr, H3K4me3_bw2_fltr)
H3K4me3_fetal_peaks$score <- 1
H3K4me3_fetal_peaks %>% export.bw(con = 'data_/CT/H3K4me3_fetal_peaks.bw')
H3K4me3_fetal_peaks %>% export.bed(con = 'data_/CT/H3K4me3_fetal_peaks.bed')

H3K4me3_fetal_peaks_chr <- H3K4me3_fetal_peaks %>% GRangesToString()

fet_H3K4me3_ct_olap <- findOverlaps(StringToGRanges(all_peaks_srt$H3K4me3[!is.na(all_peaks_srt$H3K4me3)]), H3K4me3_fetal_peaks)

H3K4me3_ct_tbl <- tibble(
    region = H3K4me3_fetal_peaks_chr,
    in_organoids = H3K4me3_fetal_peaks_chr %in% H3K4me3_fetal_peaks_chr[subjectHits(fet_H3K4me3_ct_olap)]
)

p1 <- ggplot(H3K4me3_ct_tbl, aes('fetalH3K4me3_in_orgH3K4me3', fill=in_organoids)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_y_continuous(limits=c(0,40000)) +
    labs(x='') + no_legend() + 
    scale_axis_rangeframe() + theme_rangeframe() +
    ggtitle('H3K4me3')

ct_H3K4me3_tbl <- tibble(
    region = all_peaks_srt$repr_region,
    in_organoids = all_peaks_srt$repr_region %in%  all_peaks_srt$repr_region[queryHits(fet_H3K4me3_ct_olap)]
)

p2 <- ggplot(ct_H3K4me3_tbl, aes('orgH3K4me3_in_fetalH3K4me3', fill=in_organoids)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_y_continuous(limits=c(0,40000)) +
    labs(x='') + no_legend() + 
    scale_axis_rangeframe() + theme_rangeframe()

pH3K4me3 <- p1 | p2 
pH3K4me3



H3K27me3_bw1 <- import.bw('/nas/groups/treutlein/DATA/sequencing/20230310_P2324_FIDES_ATAC/processed_bulk/mapped/bamCoverage/230117_13_FZ_bCT_fb-w19_H3K27me3_S1.filtered.seq_depth_norm.bw')
H3K27me3_bw2 <- import.bw('/nas/groups/treutlein/DATA/sequencing/20230525_P2426_FIDES/mapped/bamCoverage/23011709_23_FZ_bCT_fb-w19_H3K27me3_S1.filtered.seq_depth_norm.bw')

H3K27me3_bw1_fltr <- H3K27me3_bw1[H3K27me3_bw1$score>20] %>% {intersect(.,.)} %>% {.[width(.)>200]}
H3K27me3_bw2_fltr <- H3K27me3_bw2[H3K27me3_bw2$score>20] %>% {intersect(.,.)} %>% {.[width(.)>200]}

H3K27me3_fetal_peaks <- intersect(H3K27me3_bw1_fltr, H3K27me3_bw2_fltr)
H3K27me3_fetal_peaks$score <- 1
H3K27me3_fetal_peaks %>% export.bw(con = 'data_/CT/H3K27me3_fetal_peaks.bw')
H3K27me3_fetal_peaks %>% export.bed(con = 'data_/CT/H3K27me3_fetal_peaks.bed')

H3K27me3_fetal_peaks_chr <- H3K27me3_fetal_peaks %>% GRangesToString()

fet_H3K27me3_ct_olap <- findOverlaps(StringToGRanges(all_peaks_srt$H3K27me3[!is.na(all_peaks_srt$H3K27me3)]), H3K27me3_fetal_peaks)

H3K27me3_ct_tbl <- tibble(
    region = H3K27me3_fetal_peaks_chr,
    in_organoids = H3K27me3_fetal_peaks_chr %in% H3K27me3_fetal_peaks_chr[subjectHits(fet_H3K27me3_ct_olap)]
)

p1 <- ggplot(H3K27me3_ct_tbl, aes('fetalH3K27me3_in_orgH3K27me3', fill=in_organoids)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_y_continuous(limits=c(0,40000)) +
    labs(x='') + no_legend() + 
    scale_axis_rangeframe() + theme_rangeframe() +
    ggtitle('H3K27me3')

ct_H3K27me3_tbl <- tibble(
    region = all_peaks_srt$repr_region,
    in_organoids = all_peaks_srt$repr_region %in%  all_peaks_srt$repr_region[queryHits(fet_H3K27me3_ct_olap)]
)

p2 <- ggplot(ct_H3K27me3_tbl, aes('orgH3K27me3_in_fetalH3K4me3', fill=in_organoids)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_y_continuous(limits=c(0,40000)) +
    labs(x='') + no_legend() + 
    scale_axis_rangeframe() + theme_rangeframe()

pH3K27me3 <- p1 | p2 





H3K27ac_bw1 <- import.bw('/nas/groups/treutlein/DATA/sequencing/20230310_P2324_FIDES_ATAC/processed_bulk/mapped/bamCoverage/230117_14_FZ_bCT_fb-w19_H3K27ac_S1.filtered.seq_depth_norm.bw')
H3K27ac_bw2 <- import.bw('/nas/groups/treutlein/DATA/sequencing/20230525_P2426_FIDES/mapped/bamCoverage/23011709_24_FZ_bCT_fb-w19_H3K27ac_S1.filtered.seq_depth_norm.bw')

H3K27ac_bw1_fltr <- H3K27ac_bw1[H3K27ac_bw1$score>20] %>% {intersect(.,.)} %>% {.[width(.)>200]}
H3K27ac_bw2_fltr <- H3K27ac_bw2[H3K27ac_bw2$score>20] %>% {intersect(.,.)} %>% {.[width(.)>200]}

H3K27ac_fetal_peaks <- intersect(H3K27ac_bw1_fltr, H3K27ac_bw2_fltr)
H3K27ac_fetal_peaks$score <- 1
H3K27ac_fetal_peaks %>% export.bw(con = 'data_/CT/H3K27ac_fetal_peaks.bw')
H3K27ac_fetal_peaks %>% export.bed(con = 'data_/CT/H3K27ac_fetal_peaks.bed')

H3K27ac_fetal_peaks_chr <- H3K27ac_fetal_peaks %>% GRangesToString()

fet_H3K27ac_ct_olap <- findOverlaps(StringToGRanges(all_peaks_srt$H3K27ac[!is.na(all_peaks_srt$H3K27ac)]), H3K27ac_fetal_peaks)

H3K27ac_ct_tbl <- tibble(
    region = H3K27ac_fetal_peaks_chr,
    in_organoids = H3K27ac_fetal_peaks_chr %in% H3K27ac_fetal_peaks_chr[subjectHits(fet_H3K27ac_ct_olap)]
)

p1 <- ggplot(H3K27ac_ct_tbl, aes('fetalH3K27ac_in_orgH3K27ac', fill=in_organoids)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_y_continuous(limits=c(0,40000)) +
    labs(x='') + no_legend() + 
    scale_axis_rangeframe() + theme_rangeframe() +
    ggtitle('H3K27ac')

ct_H3K27ac_tbl <- tibble(
    region = all_peaks_srt$repr_region,
    in_organoids = all_peaks_srt$repr_region %in%  all_peaks_srt$repr_region[queryHits(fet_H3K27ac_ct_olap)]
)

p2 <- ggplot(ct_H3K27ac_tbl, aes('orgH3K27ac_in_fetalH3K4me3', fill=in_organoids)) +
    geom_bar(color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_y_continuous(limits=c(0,40000)) +
    labs(x='') + no_legend() + 
    scale_axis_rangeframe() + theme_rangeframe()

pH3K27ac <- p1 | p2


pH3K4me3 / pH3K27me3 / pH3K27ac & article_text()
ggsave('plots/paper/CT_compare_fetal_intersection_bar.pdf', width=2, height=5)



out_tbl <- bind_rows('H3K27ac'=H3K27ac_ct_tbl, 'H3K27me3'=H3K27me3_ct_tbl, 'H3K4me3'=H3K4me3_ct_tbl, .id='mark') %>% 
    group_by(mark) %>% 
    summarize(count=n(), frac_intersect=sum(in_organoids)/n(), n_intersect=sum(in_organoids))

out_tbl %>% write_rds('data_/tables/ed5g_fetal_intersect.tsv')




















