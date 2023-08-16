source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(destiny)
library(dtw)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')

rank_cols <- function(x){
    xr <- presto::rank_matrix(t(x))$X_ranked
    colnames(xr) <- colnames(x)
    rownames(xr) <- rownames(x)
    return(xr)
}

#### Read data ###
marks  <- read_rds('data/all_marks_list_v3.4lines.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.3lines.rds')
rna_pt_meta <- read_tsv('data/RNA/cellrank/RNA_full_cellrank_probs.tsv') %>% 
    dplyr::rename('cell'=1) %>% 
    select(cell, velocity_pseudotime, pseudotime_ranks) %>% 
    column_to_rownames('cell')

rna <- AddMetaData(rna, rna_pt_meta)

gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')

#### Prep data: EB -> N.epi -> CTX ####
marks <- map(marks, function(x){
    x$ctx_traject <- (
        (x$lineage == 'ctx' |
            (x$celltype_jf %in% c('nect', 'psc') & x$age %in% c('15d', 'EB'))) &
            x$stage != 'retina'
    )
    return(x)
})

rna$ctx_traject <- (
    rna$lineage == 'ctx' |
        (rna$celltype_jf %in% c('nect', 'psc') & rna$age %in% c('15d', 'EB')) &
        rna$stage != 'retina'
)

#### Get pseudotimes with diffmaps
#### H3K27ac ####
H3K27ac_ctx <- marks$H3K27ac %>% subset(ctx_traject==TRUE)

H3K27ac_ctx <- H3K27ac_ctx %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

H3K27ac_ctx <- H3K27ac_ctx %>% 
    RunUMAP(dims=2:10, reduction='lsi')

H3K27ac_ctx <- H3K27ac_ctx %>% 
    FindNeighbors(reduction='lsi') %>% 
    FindClusters()

DepthCor(H3K27ac_ctx)
ElbowPlot(H3K27ac_ctx, reduction='lsi')
dim_plot(H3K27ac_ctx, label=T, group.by=c('clusters'))

clusters_use <- H3K27ac_ctx$clusters %>% setdiff(
    c('mid_36', 'EB_13', 'EB_17', 'late_74', 'late_83', 'late_30', 'late_43', 'late_76', 
      'late_36', 'late_77', 'late_25', 'late_33', 'late_57', 'late_88', 'late_20', 'late_4'
      ))
H3K27ac_ctx$ctx_traject <- H3K27ac_ctx$clusters %in% clusters_use
H3K27ac_ctx <- subset(H3K27ac_ctx, clusters%in%clusters_use)

H3K27ac_ctx <- H3K27ac_ctx %>% 
    RunUMAP(dims=2:10, reduction='lsi')

dim_plot(H3K27ac_ctx, label=T, group.by=c('clusters', 'ctx_traject'))

H3K27ac_pca <- H3K27ac_ctx[['lsi']]@cell.embeddings[,2:10]
H3K27ac_diffmap <- DiffusionMap(H3K27ac_pca, verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5]))

H3K27ac_ctx[['diffmap']] <- CreateDimReducObject(
    H3K27ac_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27ac_ctx <- AddMetaData(H3K27ac_ctx, H3K27ac_dc_df)

H3K27ac_ctx %>% 
    feature_plot(features=c('STMN2', 'VIM', 'DC1', 'DC2', 'DC3'), order=T)

H3K27ac_ctx %>% 
    dim_plot(group.by=c('stage', 'celltype_jf', 'seurat_clusters'), order=T)

H3K27ac_ctx %>% write_rds('data/trajectories/ctx/H3K27ac_eb_ctx_dpt_srt.rds')
# H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_dpt_srt.rds')


#### H3K4me3 ####
H3K4me3_ctx <- marks$H3K4me3 %>% subset(ctx_traject==TRUE)

H3K4me3_ctx <- H3K4me3_ctx %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

H3K4me3_ctx <- H3K4me3_ctx %>% 
    RunUMAP(dims=2:30, reduction='lsi')

H3K4me3_ctx <- H3K4me3_ctx %>% 
    FindNeighbors(reduction='lsi') %>% 
    FindClusters()

DepthCor(H3K4me3_ctx)
ElbowPlot(H3K4me3_ctx, reduction='lsi')
dim_plot(H3K4me3_ctx)
    
H3K4me3_pca <- H3K4me3_ctx[['lsi']]@cell.embeddings[,2:10]
H3K4me3_diffmap <- DiffusionMap(H3K4me3_pca, verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5]))

H3K4me3_ctx[['diffmap']] <- CreateDimReducObject(
    H3K4me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K4me3_ctx <- AddMetaData(H3K4me3_ctx, H3K4me3_dc_df)

H3K4me3_ctx %>% 
    feature_plot(features=c('STMN2', 'VIM', 'DC1', 'DC2', 'DC3'), order=F)

H3K4me3_ctx %>% 
    dim_plot(group.by=c('stage', 'celltype_jf', 'seurat_clusters'), order=T)

H3K4me3_ctx %>% write_rds('data/trajectories/ctx/H3K4me3_eb_ctx_dpt_srt.rds')
# H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_dpt_srt.rds')



#### H3K27me3 ####
H3K27me3_ctx <- marks$H3K27me3 %>% subset(ctx_traject==TRUE)

H3K27me3_ctx <- H3K27me3_ctx %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

H3K27me3_ctx <- H3K27me3_ctx %>% 
    RunUMAP(dims=2:30, reduction='lsi')

H3K27me3_ctx <- H3K27me3_ctx %>% 
    FindNeighbors(reduction='lsi', dims=2:30) %>% 
    FindClusters()

DepthCor(H3K27me3_ctx)
ElbowPlot(H3K27me3_ctx, reduction='lsi')
dim_plot(H3K27me3_ctx, label=T, group.by=c('clusters'))

# clusters_use <- H3K27me3_ctx$clusters %>% setdiff(
#     c('mo8_16', 'late_43', 'late_99', 'late_21', 
#       'late_27', 'late_75', 'late_50', 'late_64', 'late_73', 'late_63', 'late_12'
# ))
# H3K27me3_ctx$ctx_traject <- H3K27me3_ctx$clusters %in% clusters_use
# H3K27me3_ctx <- subset(H3K27me3_ctx, clusters%in%clusters_use)
# 
# dim_plot(H3K27me3_ctx, label=T, group.by=c('clusters', 'ctx_traject'))
#     
# H3K27me3_ctx <- H3K27me3_ctx %>% 
#     RunUMAP(dims=2:30, reduction='lsi')

H3K27me3_pca <- H3K27me3_ctx[['lsi']]@cell.embeddings[,2:10]
H3K27me3_diffmap <- DiffusionMap(H3K27me3_pca, verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5]))

H3K27me3_ctx[['diffmap']] <- CreateDimReducObject(
    H3K27me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27me3_ctx <- AddMetaData(H3K27me3_ctx, H3K27me3_dc_df)

H3K27me3_ctx %>% 
    feature_plot(features=c('STMN2', 'VIM', 'DC1', 'DC2', 'DC3'), order=T)

H3K27me3_ctx %>% 
    feature_plot(features=c('STMN2', 'VIM', 'EMX1', 'NEUROD6', 'GLI3'), order=T, pt.size=1)

H3K27me3_ctx %>% 
    dim_plot(group.by=c('stage', 'celltype_jf', 'seurat_clusters'), order=T)

H3K27me3_ctx %>% write_rds('data/trajectories/ctx/H3K27me3_eb_ctx_dpt_srt.rds')
# H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_dpt_srt.rds')


#### Re-compute PCA for ctx trajectory ####
rna_ctx <- rna %>% subset(ctx_traject==TRUE)

rna_ctx <- rna_ctx %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

#### Run diffmap for RNA cortex ####
rna_ctx %>% dim_plot(group.by=c('orig.ident', 'state', 'clusters'), label=T, order=T, reduction='cssumap')
rna_ctx %>% feature_plot(features=c('NEUROD6', 'EMX1', 'SOX2', 'pseudotime_ranks'), order=T, reduction='cssumap')

rna_pca <- rna_ctx[['pca']]@cell.embeddings
# rna_pca <- rna_ctx[['css']]@cell.embeddings
rna_diffmap <- DiffusionMap(rna_pca, verbose=T)
rna_ctx[['diffmap']] <- CreateDimReducObject(
    rna_diffmap@eigenvectors*100, key='DC_', assay='RNA'
)
rna_dc_df <- as.data.frame(rank_cols(rna_diffmap@eigenvectors[,1:5]))
rna_ctx <- AddMetaData(rna_ctx, rna_dc_df)

rna_ctx %>%
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'pseudotime_ranks'), reduction='cssumap')

rna_ctx %>%
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'pseudotime_ranks'), reduction='diffmap')

rna_ctx %>% write_rds('data/trajectories/ctx/RNA_eb_ctx_dpt_srt.rds')
# rna_ctx <- read_rds('data/trajectories/ctx/RNA_eb_ctx_dpt_srt.rds')


#### Align pseudotimes ####
H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_eb_ctx_dpt_srt.rds')
ctx_pt_meta <- read_tsv('data/RNA/cellrank/RNA_ctx_vpt.tsv') %>% 
    dplyr::rename('cell'=1) %>% 
    column_to_rownames('cell')

rna_ctx <- AddMetaData(rna_ctx, ctx_pt_meta)

H3K27ac_ctx$ctx_pt <- rank(-H3K27ac_ctx$DC1) / max(rank(-H3K27ac_ctx$DC1))
p1 <- H3K27ac_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

H3K4me3_ctx$ctx_pt <- rank(H3K4me3_ctx$DC1) / max(rank(H3K4me3_ctx$DC1))
p2 <- H3K4me3_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

H3K27me3_ctx$ctx_pt <- rank(H3K27me3_ctx$DC1) / max(rank(H3K27me3_ctx$DC1))
p3 <- H3K27me3_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

rna_ctx$ctx_pt <- rank(rna_ctx$ctx_vpt_ranks) / max(rank(rna_ctx$ctx_vpt_ranks))
p4 <- rna_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'pseudotime_ranks', 'ctx_pt'), order=T, reduction='cssumap')

p1 / p2 / p3 / p4
ggsave('plots/trajectories/eb_ctx_diff_pt_umap.png', width=8, height=20)

H3K27ac_ctx %>% write_rds('data/trajectories/ctx/H3K27ac_eb_ctx_dpt_srt.rds')
H3K4me3_ctx %>% write_rds('data/trajectories/ctx/H3K4me3_eb_ctx_dpt_srt.rds')
H3K27me3_ctx %>% write_rds('data/trajectories/ctx/H3K27me3_eb_ctx_dpt_srt.rds')
rna_ctx %>% write_rds('data/trajectories/ctx/RNA_eb_ctx_dpt_srt.rds')


#### Subsample modalities to minimal # cells for each timepoint ####
H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_eb_ctx_dpt_srt.rds')

H3K27ac_ctx$age <- ifelse(H3K27ac_ctx$age=='128d', '4mo', H3K27ac_ctx$age)
H3K27me3_ctx$age <- ifelse(H3K27me3_ctx$age=='128d', '4mo', H3K27me3_ctx$age)
H3K4me3_ctx$age <- ifelse(H3K4me3_ctx$age=='128d', '4mo', H3K4me3_ctx$age)

ncells_table <- bind_rows(
    table(H3K27ac_ctx$age),
    table(H3K27me3_ctx$age),
    table(H3K4me3_ctx$age),
    table(rna_ctx$age)
)

ncells_min <- colMins(as.matrix(ncells_table))
ncells_min[ncells_min<100] <- 100
names(ncells_min) <- colnames(ncells_table)


subs_list <- map(c('H3K27ac'=H3K27ac_ctx, 'H3K4me3'=H3K4me3_ctx, 'H3K27me3'=H3K27me3_ctx, 'RNA'=rna_ctx), function(srt){
    meta <- srt@meta.data %>% as_tibble(rownames='cell') %>% 
        group_by(age) %>% group_split() %>% 
        map(~sample_n(.x, size=min(ncells_min[.x$age[1]], nrow(.x)))) %>% 
        bind_rows()
    return(subset(srt, cells=meta$cell))
})

H3K27ac_ctx_subs <- subs_list$H3K27ac
H3K27me3_ctx_subs <- subs_list$H3K27me3
H3K4me3_ctx_subs <- subs_list$H3K4me3
rna_ctx_subs <- subs_list$RNA

# Re-rank pseudotime
H3K27ac_ctx_subs$ctx_pt <- rank(H3K27ac_ctx_subs$ctx_pt) / max(rank(H3K27ac_ctx_subs$ctx_pt))
H3K27me3_ctx_subs$ctx_pt <- rank(H3K27me3_ctx_subs$ctx_pt) / max(rank(H3K27me3_ctx_subs$ctx_pt))
H3K4me3_ctx_subs$ctx_pt <- rank(H3K4me3_ctx_subs$ctx_pt) / max(rank(H3K4me3_ctx_subs$ctx_pt))
rna_ctx_subs$ctx_pt <- rank(rna_ctx_subs$ctx_pt) / max(rank(rna_ctx_subs$ctx_pt))


#### Make pseudotime bins ####
H3K27ac_ctx_subs$pt_bins <- as.numeric(cut(H3K27ac_ctx_subs$ctx_pt, 50, labels=1:50))
H3K27me3_ctx_subs$pt_bins <- as.numeric(cut(H3K27me3_ctx_subs$ctx_pt, 50, labels=1:50))
H3K4me3_ctx_subs$pt_bins <- as.numeric(cut(H3K4me3_ctx_subs$ctx_pt, 50, labels=1:50))
rna_ctx_subs$pt_bins <- as.numeric(cut(rna_ctx_subs$ctx_pt, 50, labels=1:50))

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='cRNA', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='cRNA', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='cRNA', group_name='pt_bins')
rna_ctx_subs <- Pando::aggregate_assay(rna_ctx_subs, assay='RNA', group_name='pt_bins')

H3K27ac_meta <- H3K27ac_ctx_subs@meta.data %>% as_tibble(rownames='cell') 
H3K4me3_meta <- H3K4me3_ctx_subs@meta.data %>% as_tibble(rownames='cell') 
H3K27me3_meta <- H3K27me3_ctx_subs@meta.data %>% as_tibble(rownames='cell') 
rna_meta <- rna_ctx_subs@meta.data %>% as_tibble(rownames='cell') 
    
p1 <- ggplot(H3K27ac_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27ac')

p2 <- ggplot(H3K4me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K4me3')

p3 <- ggplot(H3K27me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27me3')
    
p4 <- ggplot(rna_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('rna')

p1 / p2 / p3 / p4 & theme_void() & scale_fill_manual(values=colours_timescale)
ggsave('plots/trajectories/eb_ctx_subs_age_dist_bar.png', width=8, height=5)

#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_clusters <- H3K27me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_clusters <- H3K4me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
rna_clusters <- rna_ctx_subs@assays$RNA@misc$summary$pt_bins[as.character(1:50), ]

genes_plot <- c('STMN2', 'BCL11A', 'NEUROD6', 'NEUROD2', 'POU5F1', 'SOX2', 'VIM', 'GLI3', 'GRIA2')

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
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(~gene, scales='free')

p2 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(modality~gene, scales='free')

p1 / p2 + plot_layout(heights=c(1,3))

ggsave('plots/trajectories/eb_ctx_subs_gene_expr_line.png', width=10, height=6)


H3K27ac_ctx_subs %>% write_rds('data/trajectories/ctx/H3K27ac_eb_ctx_subs_dpt_srt.rds')
H3K4me3_ctx_subs %>% write_rds('data/trajectories/ctx/H3K4me3_eb_ctx_subs_dpt_srt.rds')
H3K27me3_ctx_subs %>% write_rds('data/trajectories/ctx/H3K27me3_eb_ctx_subs_dpt_srt.rds')
rna_ctx_subs %>% write_rds('data/trajectories/ctx/RNA_eb_ctx_subs_dpt_srt.rds')



#### Subset neurogenesis ####
H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_eb_ctx_dpt_srt.rds')

H3K27ac_ctx_neuro <- H3K27ac_ctx %>% subset(stage=='late')
H3K4me3_ctx_neuro <- H3K4me3_ctx %>% subset(stage=='late')
H3K27me3_ctx_neuro <- H3K27me3_ctx %>% subset(stage=='late')
rna_ctx_neuro <- rna_ctx %>% subset(stage=='late')


#### Annotate Neurons ####
rna_ctx_neuro$ctx_pt <- rank(rna_ctx_neuro$ctx_pt) / max(rank(rna_ctx_neuro$ctx_pt))
rna_ctx_neuro$pt_bins <- as.numeric(cut(rna_ctx_neuro$ctx_pt, 20, labels=1:20))
meta <- rna_ctx_neuro@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K27ac_ctx_neuro$ctx_pt <- rank(H3K27ac_ctx_neuro$ctx_pt) / max(rank(H3K27ac_ctx_neuro$ctx_pt))
H3K27ac_ctx_neuro$pt_bins <- as.numeric(cut(H3K27ac_ctx_neuro$ctx_pt, 20, labels=1:20))
meta <- H3K27ac_ctx_neuro@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K27me3_ctx_neuro$ctx_pt <- rank(H3K27me3_ctx_neuro$ctx_pt) / max(rank(H3K27me3_ctx_neuro$ctx_pt))
H3K27me3_ctx_neuro$pt_bins <- as.numeric(cut(H3K27me3_ctx_neuro$ctx_pt, 20, labels=1:20))
meta <- H3K27me3_ctx_neuro@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K4me3_ctx_neuro$ctx_pt <- rank(H3K4me3_ctx_neuro$ctx_pt) / max(rank(H3K4me3_ctx_neuro$ctx_pt))
H3K4me3_ctx_neuro$pt_bins <- as.numeric(cut(H3K4me3_ctx_neuro$ctx_pt, 20, labels=1:20))
meta <- H3K4me3_ctx_neuro@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()



rna_ctx_neuro$neuron <- rna_ctx_neuro$pt_bins > 15
H3K27ac_ctx_neuro$neuron <- H3K27ac_ctx_neuro$pt_bins > 15
H3K27me3_ctx_neuro$neuron <- H3K27me3_ctx_neuro$pt_bins > 15
H3K4me3_ctx_neuro$neuron <- H3K4me3_ctx_neuro$pt_bins > 15

H3K27ac_ctx_neuro %>% write_rds('data/trajectories/ctx/H3K27ac_ctx_neuro_dpt_srt.rds')
H3K4me3_ctx_neuro %>% write_rds('data/trajectories/ctx/H3K4me3_ctx_neuro_dpt_srt.rds')
H3K27me3_ctx_neuro %>% write_rds('data/trajectories/ctx/H3K27me3_ctx_neuro_dpt_srt.rds')
rna_ctx_neuro %>% write_rds('data/trajectories/ctx/RNA_ctx_neuro_dpt_srt.rds')



#### GAMs for on pt_bin gene expression -> represent each timepoint by all modalities ####
# -> multivariate DTW distance between GAM smooths
# -> K-means clustering -> identify lagging/parallel/antiparallel clusters of genes

H3K27ac_ctx_subs <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_subs_dpt_srt.rds')
H3K4me3_ctx_subs <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_subs_dpt_srt.rds')
H3K27me3_ctx_subs <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_subs_dpt_srt.rds')
rna_ctx_subs <- read_rds('data/trajectories/ctx/RNA_eb_ctx_subs_dpt_srt.rds')

mark_feats <- intersect(rownames(H3K27ac_ctx_subs[['cRNA']]), rownames(H3K27me3_ctx_subs[['cRNA']])) %>% intersect(rownames(H3K4me3_ctx_subs[['cRNA']]))

H3K27ac_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K27ac_ctx_subs[['cRNA']]@data > 0)*1)
H3K27me3_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K27me3_ctx_subs[['cRNA']]@data > 0)*1)
H3K4me3_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K4me3_ctx_subs[['cRNA']]@data > 0)*1)

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='cRNA_bin', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='cRNA_bin', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='cRNA_bin', group_name='pt_bins')

H3K27ac_cluster_detect <- H3K27ac_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), mark_feats]
H3K27me3_cluster_detect <- H3K27me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), mark_feats]
H3K4me3_cluster_detect <- H3K4me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), mark_feats]

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
    top_n(6000, vst.variance) %>% 
    pull(feature) %>% 
    intersect(mark_detect_feats) -> feats_use


#### Smooth with GAMs ####
x <- 1:50
H3K27ac_gams <- H3K27ac_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
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

# l2_smooth_dist <- map_par(names(all_smooths_list), function(i){
#     map_dfr(names(all_smooths_list), function(j){
#         ij_dist <- mean(diag(dist(all_smooths_list[[i]], all_smooths_list[[j]])))
#         return(tibble(
#             i=i, j=j, dist=ij_dist
#         ))
#     })
# }, parallel=T)
# 
# l2_smooth_dist_df <- bind_rows(l2_smooth_dist)
# l2_smooth_dist_mat <- l2_smooth_dist_df %>% 
#     pivot_wider(names_from=j, values_from=dist) %>% 
#     column_to_rownames('i') %>% as.matrix()
# 
# pheatmap::pheatmap(l2_smooth_dist_mat)


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

# l2_kmeans <- FCPS::kmeansClustering(
#     DataOrDistances = l2_smooth_dist_mat,
#     ClusterNo = 10
# )

# hclust_df <- dtw_smooth_dist_mat %>% as.dist() %>% hclust() %>% cutree(k=10) %>% enframe('feature', 'h_clust')
dtw_kmeans_df <- dtw_kmeans$Cls %>% enframe('feature', 'dtw_clust')
dtw_5_kmeans_df <- dtw_5_kmeans$Cls %>% enframe('feature', 'dtw_clust_5')
# l2_kmeans_df <- l2_kmeans$Cls %>% enframe('feature', 'l2_clust')

#### Plot clusters #####
all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
        {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=1:50) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr')

plot_df <- all_gene_smooths_df %>% 
    inner_join(dtw_kmeans_df) %>%
    inner_join(dtw_5_kmeans_df) %>%
    # inner_join(l2_kmeans_df) %>% 
    # inner_join(hclust_df) %>% 
    group_by(feature, modality) %>% 
    mutate(expr01=scale01(expr), pt_bin_mod=paste0(modality, pt_bin)) 

plot_df$feature <- factor(plot_df$feature, levels=h_clust(plot_df, feature, pt_bin_mod, expr01))


ggplot(plot_df, aes(pt_bin, expr01, color=modality, group=interaction(feature, modality))) +
    geom_line(alpha=0.5, size=0.5) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_wrap(~dtw_clust) 
ggsave('plots/trajectories/ctx_subs_gam_smooth_10clust_line.png', width=15, height=6)

ggplot(plot_df, aes(pt_bin, expr01, color=modality, group=interaction(feature, modality))) +
    geom_line(alpha=0.5, size=0.5) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(dtw_clust~modality, scales='free') 
ggsave('plots/trajectories/ctx_subs_gam_smooth_dtw_10clust_split_line.png', width=8, height=10)

ggplot(plot_df, aes(pt_bin, expr01, color=modality)) +
    geom_smooth(alpha=0.5, size=1) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_wrap(~dtw_clust) 
ggsave('plots/trajectories/ctx_subs_gam_smooth_10clust_smooth.png', width=8, height=4)

ggplot(plot_df, aes(pt_bin, expr01, color=modality, group=interaction(feature, modality))) +
    geom_line(alpha=0.5, size=0.5) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_wrap(~dtw_clust_5) 
ggsave('plots/trajectories/ctx_subs_gam_smooth_5clust_line.png', width=12, height=4)

ggplot(plot_df, aes(pt_bin, expr01, color=modality, group=interaction(feature, modality))) +
    geom_line(alpha=0.5, size=0.5) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(dtw_clust_5~modality, scales='free') 
ggsave('plots/trajectories/ctx_subs_gam_smooth_dtw_5clust_split_line.png', width=8, height=5)

ggplot(plot_df, aes(pt_bin, expr01, color=modality)) +
    geom_smooth(alpha=0.5, size=1) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_wrap(~dtw_clust_5) 
ggsave('plots/trajectories/ctx_subs_gam_smooth_5clust_smooth.png', width=6, height=3)


ggplot(plot_df, aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(dtw_clust~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100))
ggsave('plots/trajectories/ctx_gam_smooth_10clust_heatmap.png', width=20, height=35)


ggplot(plot_df, aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(dtw_clust_5~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100))
ggsave('plots/trajectories/ctx_gam_smooth_5clust_heatmap.png', width=20, height=35)



#### Separate enhancer and promoter ####
#### Add activities ####
enhancer_act <- read_rds('data/gene_activity/CT_distal_activity.rds')
prom_act <- read_rds('data/gene_activity/CT_promoter_activity.rds')
gb_act <- read_rds('data/gene_activity/CT_gene_body_activity.rds')


H3K27ac_ctx_subs[['eRNA']] <- CreateAssayObject(enhancer_act$H3K27ac[, colnames(H3K27ac_ctx_subs)])
H3K27me3_ctx_subs[['eRNA']] <- CreateAssayObject(enhancer_act$H3K27me3[, colnames(H3K27me3_ctx_subs)])
H3K4me3_ctx_subs[['eRNA']] <- CreateAssayObject(enhancer_act$H3K4me3[, colnames(H3K4me3_ctx_subs)])

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='eRNA', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='eRNA', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='eRNA', group_name='pt_bins')

H3K27ac_ctx_subs[['pRNA']] <- CreateAssayObject(Pando::aggregate_matrix(prom_act$H3K27ac, rownames(prom_act$H3K27ac))[, colnames(H3K27ac_ctx_subs)])
H3K27me3_ctx_subs[['pRNA']] <- CreateAssayObject(Pando::aggregate_matrix(prom_act$H3K27me3, rownames(prom_act$H3K27me3))[, colnames(H3K27me3_ctx_subs)])
H3K4me3_ctx_subs[['pRNA']] <- CreateAssayObject(Pando::aggregate_matrix(prom_act$H3K4me3, rownames(prom_act$H3K4me3))[, colnames(H3K4me3_ctx_subs)])

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='pRNA', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='pRNA', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='pRNA', group_name='pt_bins')

H3K27ac_ctx_subs[['gbRNA']] <- CreateAssayObject(gb_act$H3K27ac[, colnames(H3K27ac_ctx_subs)])
H3K27me3_ctx_subs[['gbRNA']] <- CreateAssayObject(gb_act$H3K27me3[, colnames(H3K27me3_ctx_subs)])
H3K4me3_ctx_subs[['gbRNA']] <- CreateAssayObject(gb_act$H3K4me3[, colnames(H3K4me3_ctx_subs)])

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='gbRNA', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='gbRNA', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='gbRNA', group_name='pt_bins')


H3K27ac_ctx_subs  %>% write_rds('data/trajectories/ctx/H3K27ac_eb_ctx_subs_dpt_srt.rds')
H3K4me3_ctx_subs  %>% write_rds('data/trajectories/ctx/H3K4me3_eb_ctx_subs_dpt_srt.rds')
H3K27me3_ctx_subs  %>% write_rds('data/trajectories/ctx/H3K27me3_eb_ctx_subs_dpt_srt.rds')



#### Get pt bin summaries ####
H3K27ac_enh_clusters <- H3K27ac_ctx_subs@assays$eRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_enh_clusters <- H3K27me3_ctx_subs@assays$eRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_enh_clusters <- H3K4me3_ctx_subs@assays$eRNA@misc$summary$pt_bins[as.character(1:50), ]

H3K27ac_prom_clusters <- H3K27ac_ctx_subs@assays$pRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_prom_clusters <- H3K27me3_ctx_subs@assays$pRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_prom_clusters <- H3K4me3_ctx_subs@assays$pRNA@misc$summary$pt_bins[as.character(1:50), ]

H3K27ac_gb_clusters <- H3K27ac_ctx_subs@assays$gbRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_gb_clusters <- H3K27me3_ctx_subs@assays$gbRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_gb_clusters <- H3K4me3_ctx_subs@assays$gbRNA@misc$summary$pt_bins[as.character(1:50), ]


genes_plot <- c('STMN2', 'BCL11A', 'NEUROD6', 'NEUROD2', 'POU5F1', 'SOX2', 'VIM', 'GLI3', 'GRIA2')

marks_gene_act_list <- list(
    'H3K27ac' = list(
        'enhancer'=H3K27ac_enh_clusters,
        'promoter'=H3K27ac_prom_clusters,
        'gene_body'=H3K27ac_gb_clusters
    ),
    'H3K27me3' = list(
        'enhancer'=H3K27me3_enh_clusters,
        'promoter'=H3K27me3_prom_clusters,
        'gene_body'=H3K27me3_gb_clusters
    ),
    'H3K4me3' = list(
        'enhancer'=H3K4me3_enh_clusters,
        'promoter'=H3K4me3_prom_clusters,
        'gene_body'=H3K4me3_gb_clusters
    )
)

marks_gene_act_expr <- map_dfr(marks_gene_act_list, function(m){
    map_dfr(m, function(n){
        n[, genes_plot] %>% t() %>% 
            as_tibble(rownames='gene') %>% 
            pivot_longer(!gene, names_to='pt_bins', values_to='expr') %>% 
            return()
    }, .id='window')
}, .id='modality')

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

rna_expr$modality <- 'RNA'

plot_df <- bind_rows(marks_gene_act_expr, rna_expr) %>% 
    group_by(modality, gene, window) %>% 
    mutate(expr=scale01(expr))

p1 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(~gene, scales='free')

p2 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(modality+window~gene, scales='free')

p1 / p2 + plot_layout(heights=c(1,3))



#### Match bin means by max corr ####
rna_ctx_subs <- FindVariableFeatures(rna_ctx_subs, nfeatures=2000)
H3K27ac_ctx_subs <- FindVariableFeatures(H3K27ac_ctx_subs, assay='cRNA', nfeatures=2000)
H3K27me3_ctx_subs <- FindVariableFeatures(H3K27me3_ctx_subs, assay='cRNA', nfeatures=2000)
H3K4me3_ctx_subs <- FindVariableFeatures(H3K4me3_ctx_subs, assay='cRNA', nfeatures=2000)

# RNA + H3K27ac
# var_feats <- VariableFeatures(rna_ctx_subs) %>% intersect(VariableFeatures(H3K27ac_ctx_subs, assay='cRNA'))
var_feats <- unique(filter(plot_df, dtw_clust==5)$feature)
var_feats <- feats_use

rna_H3K27ac_cor <- Pando::sparse_cor(t(rna_clusters[,var_feats]), t(H3K27ac_clusters[,var_feats]))
pheatmap::pheatmap(t(rna_H3K27ac_cor), cluster_rows=F, cluster_cols=F)

rna_H3K27ac_align <- dtw(
    x = as.matrix(1-rna_H3K27ac_cor),
    keep=T
)

rna_H3K27ac_align_df <- tibble(
    RNA_idx = rna_H3K27ac_align$index1,
    H3K27ac_idx = rna_H3K27ac_align$index2
)

ggplot(rna_H3K27ac_align_df, aes(RNA_idx, H3K27ac_idx)) +
    geom_abline() +
    geom_line()


# RNA + H3K4me3
# var_feats <- VariableFeatures(rna_ctx_subs) %>% intersect(VariableFeatures(H3K4me3_ctx_subs, assay='cRNA'))

rna_H3K4me3_cor <- Pando::sparse_cor(t(rna_clusters[,var_feats]), t(H3K4me3_clusters[,var_feats]))
pheatmap::pheatmap(t(rna_H3K4me3_cor), cluster_rows=F, cluster_cols=F)

rna_H3K4me3_align <- dtw(
    x = as.matrix(1-rna_H3K4me3_cor)
)
plot(rna_H3K4me3_align)

rna_H3K4me3_align_df <- tibble(
    RNA_idx = rna_H3K4me3_align$index1,
    H3K4me3_idx = rna_H3K4me3_align$index2
)


ggplot(rna_H3K4me3_align_df, aes(RNA_idx, H3K4me3_idx)) +
    geom_abline() +
    geom_line()


# RNA + H3K27me3
# var_feats <- VariableFeatures(rna_ctx_subs) %>% intersect(VariableFeatures(H3K27me3_ctx_subs, assay='cRNA'))

rna_H3K27me3_cor <- Pando::sparse_cor(t(rna_clusters[,var_feats]), t(H3K27me3_clusters[,var_feats]))
pheatmap::pheatmap(t(rna_H3K27me3_cor), cluster_rows=F, cluster_cols=F)

rna_H3K27me3_align <- dtw(
    x = as.matrix(rna_H3K27me3_cor)
)
plot(rna_H3K27me3_align)


rna_H3K27me3_align_df <- tibble(
    RNA_idx = rna_H3K27me3_align$index1,
    H3K27me3_idx = rna_H3K27me3_align$index2
)

ggplot(rna_H3K27me3_align_df, aes(RNA_idx, H3K27me3_idx)) +
    geom_abline() +
    geom_line()




#### Pt patterns on peak level ####
#### Feature selection ####

H3K27ac_ctx_subs <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_subs_dpt_srt.rds')
H3K4me3_ctx_subs <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_subs_dpt_srt.rds')
H3K27me3_ctx_subs <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_subs_dpt_srt.rds')
rna_ctx_subs <- read_rds('data/trajectories/ctx/RNA_eb_ctx_subs_dpt_srt.rds')

peak_intersects <- read_tsv('data/intersect/all_marks_intersect_matches.tsv')

H3K27ac_ctx_subs[['peaks_bin']] <- CreateAssayObject((H3K27ac_ctx_subs[['peaks']]@data > 0)*1)
H3K27me3_ctx_subs[['peaks_bin']] <- CreateAssayObject((H3K27me3_ctx_subs[['peaks']]@data > 0)*1)
H3K4me3_ctx_subs[['peaks_bin']] <- CreateAssayObject((H3K4me3_ctx_subs[['peaks']]@data > 0)*1)

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='peaks_bin', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='peaks_bin', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='peaks_bin', group_name='pt_bins')

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='peaks', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='peaks', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='peaks', group_name='pt_bins')

H3K27ac_clusters <- H3K27ac_ctx_subs@assays$peaks@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_clusters <- H3K27me3_ctx_subs@assays$peaks@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_clusters <- H3K4me3_ctx_subs@assays$peaks@misc$summary$pt_bins[as.character(1:50), ]
rna_clusters <- rna_ctx_subs@assays$RNA@misc$summary$pt_bins[as.character(1:50), ]

H3K27ac_cluster_detect <- H3K27ac_ctx_subs@assays$peaks_bin@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_cluster_detect <- H3K27me3_ctx_subs@assays$peaks_bin@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_cluster_detect <- H3K4me3_ctx_subs@assays$peaks_bin@misc$summary$pt_bins[as.character(1:50), ]

H3K27ac_maxdetect_feats <- sparseMatrixStats::colMaxs(H3K27ac_cluster_detect)
H3K27me3_maxdetect_feats <- sparseMatrixStats::colMaxs(H3K27me3_cluster_detect)
H3K4me3_maxdetect_feats <- sparseMatrixStats::colMaxs(H3K4me3_cluster_detect)

names(H3K27ac_maxdetect_feats) <- colnames(H3K27ac_cluster_detect)
names(H3K27me3_maxdetect_feats) <- colnames(H3K27me3_cluster_detect)
names(H3K4me3_maxdetect_feats) <- colnames(H3K4me3_cluster_detect)

H3K27ac_maxdetect_df <- H3K27ac_maxdetect_feats %>% 
    enframe('H3K27ac', 'detection')

gene_annot_use <- gene_annot[gene_annot$gene_name%in%rownames(rna)]
isect_genes <- ClosestFeature(H3K27ac_ctx_subs, unique(peak_intersects$isect), annotation=gene_annot_use) %>% 
    as_tibble() %>% select('isect'=query_region, 'RNA'=gene_name) 

isects_annot <- peak_intersects %>% 
    inner_join(isect_genes) %>% 
    inner_join(H3K27ac_maxdetect_df) %>% 
    distinct()

rna_var_feats <- all_var_feats %>% 
    filter(modality=='RNA') %>% 
    top_n(4000, vst.variance) %>% 
    pull(feature) 

isects_use <- isects_annot %>% 
    filter(RNA%in%rna_var_feats, detection>0.02)


#### Plot individual peaks / genes ####
genes_plot <- c('STMN2', 'BCL11A', 'NEUROD6', 'NEUROD2', 'POU5F1', 'SOX2', 'VIM', 'GLI3', 'GRIA2')

isects_plot <- isects_annot %>% 
    filter(RNA%in%genes_plot, detection>0.02)

H3K27ac_expr_plot <- H3K27ac_clusters[, unique(isects_plot$H3K27ac)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(isects_plot, by=c('feature'='H3K27ac'), keep=T)

H3K27me3_expr_plot <- H3K27me3_clusters[, unique(isects_plot$H3K27me3)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(isects_plot, by=c('feature'='H3K27me3'), keep=T)

H3K4me3_expr_plot <- H3K4me3_clusters[, unique(isects_plot$H3K4me3)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(isects_plot, by=c('feature'='H3K4me3'), keep=T)

rna_expr_plot <- rna_clusters[, unique(isects_plot$RNA)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(isects_plot, by=c('feature'='RNA'), keep=T)

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr_plot, 'H3K27me3'=H3K27me3_expr_plot, 'H3K4me3'=H3K4me3_expr_plot, 'RNA'=rna_expr_plot, .id='modality') %>% 
    group_by(isect, RNA, modality, pt_bins) %>% 
    summarize(expr=mean(expr)) %>% 
    group_by(modality, isect) %>% 
    mutate(expr=scale01(expr))

ggplot(filter(plot_df, RNA=='GLI3'), aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=mod_colors) +
    facet_grid(RNA~isect, scales='free')



#### Smooth with GAMs ####
x <- 1:50
H3K27ac_gams <- H3K27ac_clusters[, unique(isects_use$H3K27ac)] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_clusters[, unique(isects_use$H3K27me3)] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_clusters[, unique(isects_use$H3K4me3)] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_clusters[, unique(isects_use$RNA)] %>%
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
rna_smooths <- rna_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27ac_smooths_scale <- H3K27ac_smooths %>% scale() %>% {.[is.na(.)]<-0;.}
H3K27me3_smooths_scale <- H3K27me3_smooths %>% scale() %>% {.[is.na(.)]<-0;.}
H3K4me3_smooths_scale <- H3K4me3_smooths %>% scale() %>% {.[is.na(.)]<-0;.}
rna_smooths_scale <- rna_smooths %>% scale() %>% {.[is.na(.)]<-0;.}

all_smooths_list <- map(set_names(unique(isects_use$isect)), function(i){
    idf <- filter(isects_use, isect==i)
    do.call(
        cbind, 
        list(
            H3K27ac=rowMeans(H3K27ac_smooths_scale[,idf$H3K27ac,drop=F]),
            H3K27me3=rowMeans(H3K27me3_smooths_scale[,idf$H3K27me3,drop=F]),
            H3K4me3=rowMeans(H3K4me3_smooths_scale[,idf$H3K4me3,drop=F]),
            rna=rowMeans(rna_smooths_scale[,idf$RNA,drop=F])
        )
    )
})

all_smooths_list %>% write_rds('data/trajectories/ctx/H3K27ac_detect_peaks_all_smooths.rds')


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

dtw_smooth_dist_mat <- bind_rows(dtw_smooth_dist) %>%
    pivot_wider(names_from=j, values_from=dist) %>%
    column_to_rownames('i') %>% as.matrix()

dtw_smooth_dist_mat %>% write_rds('data/trajectories/ctx/H3K27ac_detect_peaks_dtw_dist.rds')


#### K-means clustering ####
library(FCPS)

dtw_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 10
)

dtw_kmeans_df <- dtw_kmeans$Cls %>% enframe('feature', 'dtw_clust')


#### Plot clusters #####
all_gene_smooths_df <- map_dfr(all_smooths_list, function(x){
    x %>% 
        {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=1:50) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr')

plot_df <- all_gene_smooths_df %>% 
    inner_join(dtw_kmeans_df) %>%
    inner_join(isects_annot, by=c('feature'='isect')) %>%
    group_by(feature, modality) %>% 
    mutate(
        expr01=scale01(expr), 
        expr01=ifelse(is.na(expr01), 0, expr01),
        pt_bin_mod=paste0(modality, pt_bin)
    ) 

ggplot(filter(plot_df, dtw_clust==5), aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(dtw_clust~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100))


ggplot(filter(plot_df, RNA=='NEUROD6'), aes(pt_bin, expr01, color=modality)) +
    geom_smooth(alpha=0.5, size=1) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_wrap(~dtw_clust) 

ggplot(plot_df, aes(pt_bin, expr01, color=modality)) +
    geom_smooth(alpha=0.5, size=1) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_wrap(~dtw_clust) 





