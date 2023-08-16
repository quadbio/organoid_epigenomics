source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(destiny)
library(presto)
library(dtw)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')

mark_colors <- c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')

rank_cols <- function(x){
    xr <- presto::rank_matrix(t(x))$X_ranked
    colnames(xr) <- colnames(x)
    rownames(xr) <- rownames(x)
    return(xr)
}

find_peaks_near_genes <- function(
        peaks, genes, distance = 200000, sep = c('-', '-'), only_tss = FALSE
){
    if (only_tss){
        genes <- resize(x = genes, width = 1, fix = 'start')
        genes_extended <- suppressWarnings(
            expr = Extend(
                genes, upstream = distance, downstream = distance
            )
        )
    } else {
        genes_extended <- suppressWarnings(
            expr = Extend(
                genes, upstream = distance, downstream = distance
            )
        )
    }
    overlaps <- findOverlaps(
        query = peaks,
        subject = genes_extended,
        type = 'any',
        select = 'all'
    )
    hit_matrix <- sparseMatrix(
        i = queryHits(overlaps),
        j = subjectHits(overlaps),
        x = 1,
        dims = c(length(peaks), length(genes_extended))
    )
    rownames(hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
    colnames(hit_matrix) <- genes_extended$gene_name
    return(hit_matrix)
}


#### Read data ####
marks  <- read_rds('data/all_marks_list_v3.4lines.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.3lines.rds')
rna_pt_meta <- read_tsv('data/RNA/cellrank/RNA_full_cellrank_probs.tsv') %>% 
    dplyr::rename('cell'=1) %>% 
    select(cell, velocity_pseudotime, pseudotime_ranks) %>% 
    column_to_rownames('cell')

rna <- AddMetaData(rna, rna_pt_meta)

gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')


#### Prep data: CTX neurogenesis d35 ####
marks <- map(marks, function(x){
    x$ctx_traject <- (
        (x$lineage == 'ctx') &
            (x$age == '128d')
    )
    return(x)
})

rna$ctx_traject <- (
    (rna$lineage == 'ctx') &
        (rna$age == '4mo')
)

#### Get pseudotimes with diffmaps
#### H3K27ac ####
H3K27ac_ctx <- marks$H3K27ac %>% subset(ctx_traject==TRUE)

H3K27ac_ctx <- subset(H3K27ac_ctx, log_peak_counts>3)

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

H3K27ac_pca <- H3K27ac_ctx[['lsi']]@cell.embeddings[,2:10]
H3K27ac_diffmap <- DiffusionMap(H3K27ac_pca, verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5]))

H3K27ac_ctx[['diffmap']] <- CreateDimReducObject(
    H3K27ac_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27ac_ctx <- AddMetaData(H3K27ac_ctx, H3K27ac_dc_df)

H3K27ac_ctx %>% 
    FeaturePlot(features=c('STMN2', 'VIM', 'log_peak_counts', 'DC1', 'DC2', 'DC3'), order=T)

H3K27ac_ctx %>% 
    dim_plot(group.by=c('stage', 'celltype_jf', 'seurat_clusters'), order=T)

H3K27ac_ctx %>% write_rds('data/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
# H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')




#### H3K4me3 ####
H3K4me3_ctx <- marks$H3K4me3 %>% subset(ctx_traject==TRUE)

H3K4me3_ctx$log_peak_counts %>% hist()
H3K4me3_ctx <- subset(H3K4me3_ctx, log_peak_counts>2.5 & log_peak_counts<4)

H3K4me3_ctx <- H3K4me3_ctx %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

H3K4me3_ctx <- H3K4me3_ctx %>% 
    RunUMAP(dims=2:10, reduction='lsi')

H3K4me3_ctx <- H3K4me3_ctx %>% 
    FindNeighbors(reduction='lsi') %>% 
    FindClusters()

DepthCor(H3K4me3_ctx)
ElbowPlot(H3K4me3_ctx, reduction='lsi')
dim_plot(H3K4me3_ctx, label=T, group.by=c('clusters'))

H3K4me3_pca <- H3K4me3_ctx[['lsi']]@cell.embeddings[,2:10]
H3K4me3_diffmap <- DiffusionMap(H3K4me3_pca, verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5]))

H3K4me3_ctx[['diffmap']] <- CreateDimReducObject(
    H3K4me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K4me3_ctx <- AddMetaData(H3K4me3_ctx, H3K4me3_dc_df)

H3K4me3_ctx@active.assay <- 'cRNA'
H3K4me3_ctx %>% 
    FeaturePlot(features=c('STMN2', 'VIM', 'log_peak_counts', 'DC1', 'DC2', 'DC3'), order=T)

H3K4me3_ctx %>% 
    dim_plot(group.by=c('stage', 'celltype_jf', 'seurat_clusters'), order=T)

H3K4me3_ctx %>% write_rds('data/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
# H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')



#### H3K27me3 ####
H3K27me3_ctx <- marks$H3K27me3 %>% subset(ctx_traject==TRUE)

H3K27me3_ctx$log_peak_counts %>% hist()
H3K27me3_ctx <- subset(H3K27me3_ctx, log_peak_counts<3.5)

H3K27me3_ctx <- H3K27me3_ctx %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

H3K27me3_ctx <- H3K27me3_ctx %>% 
    RunUMAP(dims=2:10, reduction='lsi')

H3K27me3_ctx <- H3K27me3_ctx %>% 
    FindNeighbors(reduction='lsi') %>% 
    FindClusters()

DepthCor(H3K27me3_ctx)
ElbowPlot(H3K27me3_ctx, reduction='lsi')
dim_plot(H3K27me3_ctx, label=T, group.by=c('clusters'))

H3K27me3_pca <- H3K27me3_ctx[['lsi']]@cell.embeddings[,2:10]
H3K27me3_diffmap <- DiffusionMap(H3K27me3_pca, verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5]))

H3K27me3_ctx[['diffmap']] <- CreateDimReducObject(
    H3K27me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27me3_ctx <- AddMetaData(H3K27me3_ctx, H3K27me3_dc_df)

H3K27me3_ctx %>% 
    FeaturePlot(features=c('crna_STMN2', 'crna_VIM', 'log_peak_counts', 'DC1', 'DC2', 'DC3'), order=T)

H3K27me3_ctx %>% 
    dim_plot(group.by=c('stage', 'celltype_jf', 'seurat_clusters'), order=T)

H3K27me3_ctx %>% write_rds('data/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
# H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')


#### Re-compute PCA for ctx trajectory ####
rna_ctx <- rna %>% subset(ctx_traject==TRUE)

rna_ctx <- rna_ctx %>% 
    FindVariableFeatures() %>% 
    {VariableFeatures(.) <- setdiff(VariableFeatures(.), cc_genes_all); .} %>% 
    ScaleData(vars.to.regress=c('S.Score', 'G2M.Score')) %>% 
    RunPCA()

rna_ctx <- rna_ctx %>% 
    RunUMAP(dims=1:20)


#### Run diffmap for RNA cortex ####
rna_ctx %>% dim_plot(group.by=c('orig.ident', 'state', 'clusters'), label=T, order=T, reduction='umap')
rna_ctx %>% feature_plot(features=c('NEUROD6', 'EMX1', 'SOX2', 'MKI67', 'pseudotime_ranks'), order=T, reduction='umap')

rna_pca <- rna_ctx[['pca']]@cell.embeddings[,1:20]
rna_diffmap <- DiffusionMap(rna_pca, verbose=T)
rna_ctx[['diffmap']] <- CreateDimReducObject(
    rna_diffmap@eigenvectors*100, key='DC_', assay='RNA'
)
rna_dc_df <- as.data.frame(rank_cols(rna_diffmap@eigenvectors[,1:5]))
rna_ctx <- AddMetaData(rna_ctx, rna_dc_df)

rna_ctx %>%
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'pseudotime_ranks'), reduction='umap')

rna_ctx %>%
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'pseudotime_ranks'), reduction='diffmap')

rna_ctx %>% write_rds('data/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')
# rna_ctx <- read_rds('data/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')



#### Align pseudotimes ####
H3K27ac_ctx <- read_rds('data_/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data_/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data_/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data_/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')

H3K27ac_ctx$ctx_pt <- rank(H3K27ac_ctx$DC1) / max(rank(H3K27ac_ctx$DC1))
p1 <- H3K27ac_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

H3K4me3_ctx$ctx_pt <- rank(-H3K4me3_ctx$DC1) / max(rank(-H3K4me3_ctx$DC1))
p2 <- H3K4me3_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

H3K27me3_ctx$ctx_pt <- rank(H3K27me3_ctx$DC1) / max(rank(H3K27me3_ctx$DC1))
p3 <- H3K27me3_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)

rna_ctx$ctx_pt <- rank(-rna_ctx$DC1) / max(rank(-rna_ctx$DC1))
p4 <- rna_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'pseudotime_ranks', 'ctx_pt'), order=T, reduction='umap')

p1 / p2 / p3 / p4
ggsave('plots/trajectories/2m_ctx_diff_pt_umap_.png', width=8, height=20)

H3K27ac_ctx %>% write_rds('data/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx %>% write_rds('data/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx %>% write_rds('data/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx %>% write_rds('data/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')



#### Make pseudotime bins ####
H3K27ac_ctx$pt_bins2 <- as.numeric(cut(H3K27ac_ctx$ctx_pt, 10, labels=1:10))
H3K27me3_ctx$pt_bins2 <- as.numeric(cut(H3K27me3_ctx$ctx_pt, 10, labels=1:10))
H3K4me3_ctx$pt_bins2 <- as.numeric(cut(H3K4me3_ctx$ctx_pt, 10, labels=1:10))
rna_ctx$pt_bins2 <- as.numeric(cut(rna_ctx$ctx_pt, 10, labels=1:10))

H3K27ac_ctx$pt_bins <- as.numeric(cut(H3K27ac_ctx$ctx_pt, 20, labels=1:20))
H3K27me3_ctx$pt_bins <- as.numeric(cut(H3K27me3_ctx$ctx_pt, 20, labels=1:20))
H3K4me3_ctx$pt_bins <- as.numeric(cut(H3K4me3_ctx$ctx_pt, 20, labels=1:20))
rna_ctx$pt_bins <- as.numeric(cut(rna_ctx$ctx_pt, 20, labels=1:20))

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='cRNA', group_name='pt_bins2')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='cRNA', group_name='pt_bins2')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='cRNA', group_name='pt_bins2')
rna_ctx <- Pando::aggregate_assay(rna_ctx, assay='RNA', group_name='pt_bins2')

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='cRNA', group_name='pt_bins')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='cRNA', group_name='pt_bins')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='cRNA', group_name='pt_bins')
rna_ctx <- Pando::aggregate_assay(rna_ctx, assay='RNA', group_name='pt_bins')

H3K27ac_ctx[['RNA_bin']] <- CreateAssayObject((H3K27ac_ctx[['cRNA']]@data > 0)*1)
H3K27me3_ctx[['RNA_bin']] <- CreateAssayObject((H3K27me3_ctx[['cRNA']]@data > 0)*1)
H3K4me3_ctx[['RNA_bin']] <- CreateAssayObject((H3K4me3_ctx[['cRNA']]@data > 0)*1)
rna_ctx[['RNA_bin']] <- CreateAssayObject((rna_ctx[['RNA']]@data > 0)*1)

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='RNA_bin', group_name='pt_bins2')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='RNA_bin', group_name='pt_bins2')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='RNA_bin', group_name='pt_bins2')
rna_ctx <- Pando::aggregate_assay(rna_ctx, assay='RNA_bin', group_name='pt_bins2')

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='RNA_bin', group_name='pt_bins')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='RNA_bin', group_name='pt_bins')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='RNA_bin', group_name='pt_bins')
rna_ctx <- Pando::aggregate_assay(rna_ctx, assay='RNA_bin', group_name='pt_bins')


#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K27me3_clusters <- H3K27me3_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K4me3_clusters <- H3K4me3_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
rna_clusters <- rna_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]

genes_plot <- c('STMN2', 'BCL11A', 'NEUROD6', 'NEUROD2', 'POU5F1', 'SOX2', 'VIM', 'GLI3', 'GRIA2')

H3K27ac_expr <- H3K27ac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins2', values_to='expr')

H3K27me3_expr <- H3K27me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins2', values_to='expr')

H3K4me3_expr <- H3K4me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins2', values_to='expr')

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins2', values_to='expr')

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr01=scale01(expr))

p1 <- ggplot(plot_df, aes(as.numeric(pt_bins2), expr01, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(~gene, scales='free')

p2 <- ggplot(plot_df, aes(as.numeric(pt_bins2), expr01, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(modality~gene, scales='free')

p1 / p2 + plot_layout(heights=c(1,3))
ggsave('plots/trajectories/2m_ctx_gene_expr_line.png', width=10, height=6)


H3K27ac_ctx %>% write_rds('data/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx %>% write_rds('data/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx %>% write_rds('data/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx %>% write_rds('data/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')



#### Plot only for NEUROD2 for figure ####
ggplot(filter(plot_df, gene=='NEUROD2'), aes(as.numeric(pt_bins2), expr, color=modality, fill=modality)) +
    geom_bar(stat='identity', position='dodge', color='black', width=0.9, linewidth=0.1) +
    # geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    scale_color_manual(values=mark_colors) +
    scale_fill_manual(values=mark_colors) +
    facet_grid(modality~gene, scales='free') +
    article_text() + no_legend() +
    labs(x='Pseudotime bins', y='Detection rate')
ggsave('plots/paper/fig3/fig3_NEUROD2_gene_act_detection_pt_bar.png', width=5, height=5, units='cm', bg='white')
ggsave('plots/paper/fig3/fig3_NEUROD2_gene_act_detection_pt_bar.pdf', width=5, height=5, units='cm', bg='white')



#### Get cluster from the previous clustering ####
H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')

gene_clusters <- read_tsv('data_/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')
neuron_genes <- gene_clusters %>% filter(dtw_clust==27) %>% pull(feature) %>% unique()

all_smooths <- read_rds('data_/trajectories/ctx/ctx_pt_detection_smooths.rds')

    


#### Get peak detection around these genes over pt ####
peak_intersects <- read_tsv('data_/intersect/all_marks_intersect_matches.tsv')
gene_annot_use <- gene_annot[gene_annot$gene_name%in%rownames(rna)]

H3K27ac_ctx[['peaks_bin']] <- CreateAssayObject((H3K27ac_ctx[['peaks']]@data > 0)*1)
H3K27me3_ctx[['peaks_bin']] <- CreateAssayObject((H3K27me3_ctx[['peaks']]@data > 0)*1)
H3K4me3_ctx[['peaks_bin']] <- CreateAssayObject((H3K4me3_ctx[['peaks']]@data > 0)*1)

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='peaks_bin', group_name='pt_bins2')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='peaks_bin', group_name='pt_bins2')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='peaks_bin', group_name='pt_bins2')

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='peaks', group_name='pt_bins2')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='peaks', group_name='pt_bins2')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='peaks', group_name='pt_bins2')

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='peaks_bin', group_name='pt_bins')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='peaks_bin', group_name='pt_bins')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='peaks_bin', group_name='pt_bins')

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='peaks', group_name='pt_bins')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='peaks', group_name='pt_bins')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='peaks', group_name='pt_bins')

H3K27ac_peaks <- H3K27ac_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K27me3_peaks <- H3K27me3_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K4me3_peaks <- H3K4me3_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]

H3K27ac_peaks_detect <- H3K27ac_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K27me3_peaks_detect <- H3K27me3_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K4me3_peaks_detect <- H3K4me3_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]

H3K27ac_peaks_expr <- H3K27ac_ctx@assays$peaks@misc$summary$pt_bins2[as.character(1:10), ]
H3K27me3_peaks_expr <- H3K27me3_ctx@assays$peaks@misc$summary$pt_bins2[as.character(1:10), ]
H3K4me3_peaks_expr <- H3K4me3_ctx@assays$peaks@misc$summary$pt_bins2[as.character(1:10), ]

H3K27ac_maxdetect_feats <- colnames(H3K27ac_peaks_detect)[colMaxs(H3K27ac_peaks_detect)>0.005]
H3K27me3_maxdetect_feats <- colnames(H3K27me3_peaks_detect)[colMaxs(H3K27me3_peaks_detect)>0.005]
H3K4me3_maxdetect_feats <- colnames(H3K4me3_peaks_detect)[colMaxs(H3K4me3_peaks_detect)>0.005]

rna_gene_expr <- rna_ctx@assays$RNA@misc$summary$pt_bins2[as.character(1:10), ]

peak_intersects_detect <- peak_intersects %>% 
    filter(
        H3K27me3 %in% H3K27me3_maxdetect_feats,
        H3K27ac %in% H3K27ac_maxdetect_feats,
        H3K4me3 %in% H3K4me3_maxdetect_feats
    )

isect2genes <- find_peaks_near_genes(
    peaks = StringToGRanges(unique(peak_intersects_detect$isect)), 
    genes = gene_annot_use,
    distance = 100000
)

isect_genes <- aggregate_matrix(t(isect2genes), groups=colnames(isect2genes), fun='sum') %>% 
    {.[rowSums(.)>0,]} %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='isect') %>% 
    filter(value!=0) %>% select(-value)

isects_annot <- peak_intersects_detect %>% 
    inner_join(isect_genes)

isects_annot %>% write_tsv('data_/intersect/all_marks_intersect_matches_annot.tsv')

isects_plot <- isects_annot %>% 
    filter(gene%in%neuron_genes) 

H3K27ac_expr_plot <- H3K27ac_peaks_detect[, unique(isects_plot$H3K27ac)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(isects_plot, by=c('feature'='H3K27ac'), keep=T)

H3K27me3_expr_plot <- H3K27me3_peaks_detect[, unique(isects_plot$H3K27me3)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(isects_plot, by=c('feature'='H3K27me3'), keep=T)

H3K4me3_expr_plot <- H3K4me3_peaks_detect[, unique(isects_plot$H3K4me3)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(isects_plot, by=c('feature'='H3K4me3'), keep=T)

rna_expr_plot <- rna_gene_expr[, unique(isects_plot$gene)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(isects_plot, by=c('feature'='gene'), keep=T)

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr_plot, 'H3K27me3'=H3K27me3_expr_plot, 'H3K4me3'=H3K4me3_expr_plot, 'RNA'=rna_expr_plot, .id='modality') %>% 
    # filter(H3K27ac=='chr17-39610343-39613968' & H3K4me3=='chr17-39611133-39613031') %>% 
    # group_by(isect, gene, modality, pt_bins) %>% 
    # filter(modality!='RNA') %>% 
    # summarize(expr=mean(expr)) %>% 
    group_by(modality, isect) %>% 
    mutate(expr01=scale01(expr)) %>% 
    group_by(gene) %>% 
    mutate(peak_num=as.numeric(factor(isect)))

plot_df %>% write_tsv('data_/trajectories/ctx/ctx_traject_primed_gene_peaks.tsv')

# ggplot(plot_df, aes(as.numeric(pt_bins), expr01, color=modality, group=modality)) +
#     geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
#     geom_point(size=0.2) +
#     scale_color_manual(values=mark_colors) +
#     facet_wrap(gene~isect, scales='free')

# ggplot(filter(plot_df, gene%in%c('NEUROD2', 'MEIS2', 'STMN2', 'BCL11A', 'GRIA2', 'POU3F3')), aes(as.numeric(pt_bins), expr01, color=modality, group=modality)) +
#     geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
#     geom_point(size=0.8) +
#     scale_color_manual(values=mark_colors) +
#     facet_grid(gene+isect~modality, scales='free')



ggplot(filter(plot_df, gene=='NEUROD2'), aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_bar(stat='identity', position='dodge', fill='lightgray') +
    # geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    scale_color_manual(values=mark_colors) +
    facet_grid(modality~gene+isect, scales='free')

region <- FindRegion(marks$H3K27ac, 'chr17-39610343-39613968') %>% 
    Extend(upstream = 1500, downstream = 1500)

H3K27ac_ctx$pt_groups <- H3K27ac_ctx$ctx_pt %>% cut(3, labels=1:3)
CoveragePlot(H3K27ac_ctx, region, group.by='pt_groups', ranges = StringToGRanges('chr17-39610343-39613968')) &
    scale_fill_manual(values=as.character(celltype_colors[c('ctx_npc', 'ctx_ip', 'ctx_ex')]))

H3K4me3_ctx$pt_groups <- H3K4me3_ctx$ctx_pt %>% cut(3, labels=1:3)
CoveragePlot(H3K4me3_ctx, region, group.by='pt_groups', ranges = StringToGRanges('chr17-39610343-39613968')) &
    scale_fill_manual(values=as.character(celltype_colors[c('ctx_npc', 'ctx_ip', 'ctx_ex')]))

H3K27me3_ctx$pt_groups <- H3K27me3_ctx$ctx_pt %>% cut(3, labels=1:3)
CoveragePlot(H3K27me3_ctx, region, group.by='pt_groups', ranges = StringToGRanges('chr17-39610343-39613968')) &
    scale_fill_manual(values=as.character(celltype_colors[c('ctx_npc', 'ctx_ip', 'ctx_ex')]))



#### Plots for figure ####
ggplot(filter(plot_df, gene=='NEUROD2'), aes(as.numeric(pt_bins), expr, color=modality, fill=modality)) +
    geom_bar(stat='identity', position='dodge', color='black', width=0.9, linewidth=0.1, alpha=0.5) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    scale_color_manual(values=mark_colors) +
    scale_fill_manual(values=mark_colors) +
    facet_grid(modality~gene+H3K27ac, scales='free') +
    article_text() + no_legend() +
    labs(x='Pseudotime bins', y='Detection rate')
ggsave('plots/paper/fig3/fig3_NEUROD2_tss_peak_detection_pt_bar.png', width=5, height=6, units='cm', bg='white')
ggsave('plots/paper/fig3/fig3_NEUROD2_tss_peak_detection_pt_bar.pdf', width=5, height=6, units='cm', bg='white')



H3K27ac_ctx$pt_groups <- H3K27ac_ctx$ctx_pt %>% cut(3, labels=1:3)
p1 <- CoveragePlot(H3K27ac_ctx, region, group.by='pt_groups', peaks=F, annotation=F) &
    scale_fill_manual(values=c('#abebc6', '#28b463', '#1d8348'))

H3K4me3_ctx$pt_groups <- H3K4me3_ctx$ctx_pt %>% cut(3, labels=1:3)
p2 <- CoveragePlot(H3K4me3_ctx, region, group.by='pt_groups', peaks=F, annotation=F) &
    scale_fill_manual(values=c('#ce93d8', '#9c27b0', '#6a1b9a'))

H3K27me3_ctx$pt_groups <- H3K27me3_ctx$ctx_pt %>% cut(3, labels=1:3)
p3 <- CoveragePlot(H3K27me3_ctx, region, group.by='pt_groups', ranges = StringToGRanges('chr17-39610343-39613968'), peaks=F) &
    scale_fill_manual(values=c('#aed6f1', '#3498db', '#21618c'))

(p1 / p2 / p3) + plot_layout(heights=c(1,1,2)) & article_text()
ggsave('plots/paper/fig3/fig3_NEUROD2_tss_peak_tracks.pdf', bg='white', width=8, height=12, units='cm')
ggsave('plots/paper/fig3/fig3_NEUROD2_tss_peak_tracks.png', bg='white', width=8, height=12, units='cm')




peak_region <- 'chr17-39580669-39635490'
region <- FindRegion(marks$H3K27ac, peak_region) %>% 
    Extend(upstream = 0, downstream = 0)
H3K27ac_ctx$pt_groups <- H3K27ac_ctx$ctx_pt %>% cut(3, labels=1:3)
p1 <- CoveragePlot(H3K27ac_ctx, region, group.by='pt_groups', peaks=F, window=500, annotation=F) &
    scale_fill_manual(values=c('#abebc6', '#28b463', '#1d8348'))

H3K4me3_ctx$pt_groups <- H3K4me3_ctx$ctx_pt %>% cut(3, labels=1:3)
p2 <- CoveragePlot(H3K4me3_ctx, region, group.by='pt_groups', peaks=F, window=500, annotation=F) &
    scale_fill_manual(values=c('#ce93d8', '#9c27b0', '#6a1b9a'))

H3K27me3_ctx$pt_groups <- H3K27me3_ctx$ctx_pt %>% cut(3, labels=1:3)
p3 <- CoveragePlot(H3K27me3_ctx, region, group.by='pt_groups', ranges = StringToGRanges('chr17-39610343-39613968'), peaks=F, window=500) &
    scale_fill_manual(values=c('#aed6f1', '#3498db', '#21618c'))

(p1 / p2 / p3) + plot_layout(heights=c(1,1,2)) & article_text()
ggsave('plots/paper/fig3/fig3_NEUROD2_full_peak_tracks.pdf', bg='white', width=8, height=12, units='cm')
ggsave('plots/paper/fig3/fig3_NEUROD2_full_peak_tracks.png', bg='white', width=8, height=12, units='cm')




gene_expr <- rna_ctx@assays$RNA@data['NEUROD2',] %>% 
    enframe('cell', 'expr')
rna_ctx$pt_groups <- rna_ctx$ctx_pt %>% cut(3, labels=1:3)
rna_meta <- rna_ctx@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(gene_expr)

ggplot(rna_meta, aes(pt_groups, expr, fill=pt_groups)) +
    geom_bar(stat='summary') +
    geom_errorbar(stat='summary', width=0.05, linewidth=0.2) +
    scale_fill_manual(values=c('#fff59d', '#fdd835', '#f57f17')) +
    article_text() +
    no_legend() +
    labs(y='Expression', x='Pseudotime bins')

ggsave('plots/paper/fig3/fig3_NEUROD2_gene_expr_bar.pdf', bg='white', width=4, height=4, units='cm')
ggsave('plots/paper/fig3/fig3_NEUROD2_gene_expr_bar.png', bg='white', width=4, height=4, units='cm')





#### Only for H3K27ac + H3K27me3 look for more distal elements ####
#### -> not great, don't show
K27_isect_matches <- read_tsv('data_/intersect/H3K27_switches_marks_intersect_matches.tsv')

switch_intersects_detect <- K27_isect_matches %>% 
    filter(
        H3K27me3 %in% H3K27me3_maxdetect_feats,
        H3K27ac %in% H3K27ac_maxdetect_feats,
    )

switch2genes <- find_peaks_near_genes(
    peaks = StringToGRanges(unique(switch_intersects_detect$isect)), 
    genes = gene_annot_use,
    distance = 500000
)

switch_genes <- aggregate_matrix(t(switch2genes), groups=colnames(switch2genes), fun='sum') %>% 
    {.[rowSums(.)>0,]} %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='isect') %>% 
    filter(value!=0) %>% select(-value)

switch_annot <- switch_intersects_detect %>% 
    inner_join(switch_genes)

switch_plot <- switch_annot %>% 
    filter(gene%in%c(neuron_genes, 'NEUROD6'))

H3K27ac_expr_plot <- H3K27ac_peaks_detect[, unique(switch_plot$H3K27ac)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(switch_plot, by=c('feature'='H3K27ac'), keep=T)

H3K27me3_expr_plot <- H3K27me3_peaks_detect[, unique(switch_plot$H3K27me3)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(switch_plot, by=c('feature'='H3K27me3'), keep=T)

rna_expr_plot <- rna_gene_expr[, unique(switch_plot$gene)] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bins', values_to='expr') %>% 
    inner_join(switch_plot, by=c('feature'='gene'), keep=T)

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr_plot, 'H3K27me3'=H3K27me3_expr_plot, 'RNA'=rna_expr_plot, .id='modality') %>% 
    filter(gene%in%c('NEUROD2', 'NEUROD6', 'GRIA2', 'MEF2C', 'STMN2', 'MIB2')) %>%
    group_by(isect, gene, modality, pt_bins) %>% 
    # filter(modality!='RNA') %>%
    # summarize(expr=mean(expr)) %>% 
    group_by(modality, isect) %>% 
    mutate(expr01=scale01(expr)) %>% 
    group_by(gene) %>% 
    mutate(peak_num=as.numeric(factor(isect)))

ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_bar(stat='identity', position='dodge', fill='lightgray') +
    # geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    scale_color_manual(values=mark_colors) +
    facet_grid(modality~gene+isect+H3K27ac, scales='free')

peak_region <- 'chr2-60568028-60570874'
region <- FindRegion(marks$H3K27ac, peak_region) %>% 
    Extend(upstream = 50000, downstream = 50000)

p1 <- CoveragePlot(H3K27ac_ctx, region, group.by='pt_groups', ranges = StringToGRanges(peak_region)) &
    scale_fill_manual(values=as.character(celltype_colors[c('ctx_npc', 'ctx_ip', 'ctx_ex')]))
p2 <- CoveragePlot(H3K27me3_ctx, region, group.by='pt_groups', ranges = StringToGRanges(peak_region)) &
    scale_fill_manual(values=as.character(celltype_colors[c('ctx_npc', 'ctx_ip', 'ctx_ex')]))
p1 / p2

peak_region <- 'chr17-39580669-39635490'
region <- FindRegion(marks$H3K27ac, peak_region) %>% 
    Extend(upstream = 0, downstream = 0)

p1 <- CoveragePlot(H3K27ac_ctx, region, group.by='pt_groups', ranges = StringToGRanges(peak_region)) &
    scale_fill_manual(values=as.character(celltype_colors[c('ctx_npc', 'ctx_ip', 'ctx_ex')]))
p2 <- CoveragePlot(H3K27me3_ctx, region, group.by='pt_groups', ranges = StringToGRanges(peak_region)) &
    scale_fill_manual(values=as.character(celltype_colors[c('ctx_npc', 'ctx_ip', 'ctx_ex')]))
p1 / p2



####

#### Quantify pt lag ####
#### -> at what pt does expression significantly diverge from the 'initial' condition ####

H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')

rna_ctx <- FindVariableFeatures(rna_ctx, nfeatures=1000)
mark_feats <- intersect(rownames(H3K27ac_ctx[['cRNA']]), rownames(H3K27me3_ctx[['cRNA']])) %>% intersect(rownames(H3K4me3_ctx[['cRNA']]))

H3K27ac_cluster_detect <- H3K27ac_ctx@assays$cRNA@misc$summary$pt_bins2[as.character(1:10), mark_feats]
H3K27me3_cluster_detect <- H3K27me3_ctx@assays$cRNA@misc$summary$pt_bins2[as.character(1:10), mark_feats]
H3K4me3_cluster_detect <- H3K4me3_ctx@assays$cRNA@misc$summary$pt_bins2[as.character(1:10), mark_feats]

H3K27ac_maxdetect_feats <- colMaxs(H3K27ac_cluster_detect)
H3K27me3_maxdetect_feats <- colMaxs(H3K27me3_cluster_detect)
H3K4me3_maxdetect_feats <- colMaxs(H3K4me3_cluster_detect)

mark_detect_feats <- mark_feats[
    (H3K27ac_maxdetect_feats > 0.05) | (H3K27me3_maxdetect_feats > 0.05) | (H3K4me3_maxdetect_feats > 0.05)
]


# Select Features
# -> variable in RNA and detected in at least one mark
# Define initial bin(s)
# Iterate over bins and test significance 
# For each gene and modality pick first bin with significant divergence 

gene_clusters <- read_tsv('data_/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')
genes_test <- gene_clusters %>% filter(dtw_clust==27) %>% pull(feature) %>% unique()

H3K27ac_ctx$nCount_RNA <- H3K27ac_ctx$nCount_cRNA
H3K27me3_ctx$nCount_RNA <- H3K27me3_ctx$nCount_cRNA
H3K4me3_ctx$nCount_RNA <- H3K4me3_ctx$nCount_cRNA

H3K27ac_ctx[['RNA']] <- H3K27ac_ctx[['cRNA']]
H3K27me3_ctx[['RNA']] <- H3K27me3_ctx[['cRNA']]
H3K4me3_ctx[['RNA']] <- H3K4me3_ctx[['cRNA']]

Idents(rna_ctx) <- 'pt_bins2'
rna_lr_pt_test <- map_dfr(set_names(as.character(2:10)), function(x){
    features_use <- intersect(genes_test, rownames(rna_ctx[['RNA']]))
    FindMarkers(rna_ctx, ident.1='1', ident.2=x, features=features_use, logfc.threshold=0, test.use='LR', latent.vars='nCount_RNA', assay='RNA') %>% 
        as_tibble(rownames='feature') %>% 
        return()
}, .id='pt_bin')

Idents(H3K27ac_ctx) <- 'pt_bins2'
H3K27ac_lr_pt_test <- map_dfr(set_names(as.character(2:10)), function(x){
    features_use <- intersect(genes_test, rownames(H3K27ac_ctx[['RNA']]))
    FindMarkers(H3K27ac_ctx, ident.1='1', ident.2=x, features=features_use, logfc.threshold=0, test.use='LR', latent.vars='nCount_RNA', assay='RNA') %>% 
        as_tibble(rownames='feature') %>% 
        return()
}, .id='pt_bin')

Idents(H3K27me3_ctx) <- 'pt_bins2'
H3K27me3_lr_pt_test <- map_dfr(set_names(as.character(2:10)), function(x){
    features_use <- intersect(genes_test, rownames(H3K27me3_ctx[['RNA']]))
    FindMarkers(H3K27me3_ctx, ident.1='1', ident.2=x, features=features_use, logfc.threshold=0, test.use='LR', latent.vars='nCount_RNA', assay='RNA') %>% 
        as_tibble(rownames='feature') %>% 
        return()
}, .id='pt_bin')

Idents(H3K4me3_ctx) <- 'pt_bins2'
H3K4me3_lr_pt_test <- map_dfr(set_names(as.character(2:10)), function(x){
    features_use <- intersect(genes_test, rownames(H3K4me3_ctx[['RNA']]))
    FindMarkers(H3K4me3_ctx, ident.1='1', ident.2=x, features=features_use, logfc.threshold=0, test.use='LR', latent.vars='nCount_RNA', assay='RNA') %>% 
        as_tibble(rownames='feature') %>% 
        return()
}, .id='pt_bin')


test_pt_bins <- function(srt_obj, features=NULL, assay='RNA_bin', family='binomial'){
    srt_obj@active.assay <- assay
    srt_obj <- DietSeurat(srt_obj, assays=assay)
    
    test_result <- Pando::map_par(set_names(as.character(2:10)), function(x){
        srt_obj <- subset(srt_obj, pt_bins2%in%c('1', x))
        srt_obj$test_var <- srt_obj$pt_bins2 == x
        
        detection_rates <- Pando::aggregate_matrix(t(srt_obj[[assay]]@counts), srt_obj$test_var)
        detrate_df <- tibble(
            feature = colnames(detection_rates),
            detect_self = as.numeric(detection_rates['TRUE', ]),
            detect_other = as.numeric(detection_rates['FALSE', ])
        ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
        
        featurs_use <- intersect(features, rownames(srt_obj[[assay]]))
        
        test_df <- lr_de(
            object = srt_obj,
            test_var = 'test_var',
            covariates = c('nCount_RNA'),
            family = family,
            slot = 'counts',
            assay = assay,
            features_use = featurs_use,
            parallel = T
        )
        test_results <- inner_join(test_df, detrate_df) %>% 
            return()
    }, parallel=F) %>% bind_rows(.id='pt_bin')
    return(test_result)
}

library(doParallel)
registerDoParallel(20)

rna_binom_pt_test <- test_pt_bins(rna_ctx, features=genes_test, assay='RNA_bin')
H3K27ac_binom_pt_test <- test_pt_bins(H3K27ac_ctx, features=genes_test, assay='RNA_bin')
H3K27me3_binom_pt_test <- test_pt_bins(H3K27me3_ctx, features=genes_test, assay='RNA_bin')
H3K4me3_binom_pt_test <- test_pt_bins(H3K4me3_ctx, features=genes_test, assay='RNA_bin')

rna_lr_pt_test <- rna_lr_pt_test %>% 
    mutate(pval=p_val, log_dr=avg_log2FC)

plot_df <- bind_rows('RNA'=rna_binom_pt_test, 'H3K27ac'=H3K27ac_binom_pt_test, 'H3K27me3'=H3K27me3_binom_pt_test, 'H3K4me3'=H3K4me3_binom_pt_test, .id='modality') %>% 
    group_by(modality) %>% mutate(padj=p.adjust(pval, 'fdr'))

plot_df %>% write_tsv('data/trajectories/ctx/ctx_pt_lr_binom_de.tsv')


all_lr_pt_de <- read_tsv('data/trajectories/ctx/ctx_pt_lr_binom_de.tsv')

plot_df <- all_lr_pt_de %>% 
    mutate(pclip=clip_abs(-log10(padj), 5))

ggplot(filter(plot_df, feature=='GRIA2'), aes(as.numeric(pt_bin), log_dr, color=modality, size=pclip)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept=c(1,-1), color='darkgrey') +
    geom_point() +
    geom_path(size=0.1) +
    facet_wrap(~feature) +
    scale_size_continuous(range=c(1,3)) +
    scale_y_continuous(limits=c(-4,4)) +
    scale_color_manual(values=modality_colors)


# Summary stat?
# -> 1st bin where there is a significance divergence in both directions

plot_df <- all_lr_pt_de %>% 
    filter(pval<0.01, log_dr>log2(1.5), modality!='H3K27me3') %>% 
    group_by(feature) %>% 
    filter(length(unique(modality))>1) %>%
    group_by(modality, feature) %>% 
    summarize(min_bin=min(pt_bin), min_dr=log_dr[which.min(pt_bin)]) %>% 
    mutate(signed_bin=sign(min_dr)*min_bin) 

ggplot(plot_df, aes(min_bin, modality, fill=modality, color=modality)) +
    geom_boxplot(linewidth=0.2, color='black', alpha=0.5) +
    geom_quasirandom(size=0.3) +
    scale_fill_manual(values=modality_colors) + 
    scale_color_manual(values=modality_colors) + 
    scale_x_continuous(limits=c(0,10)) + 
    article_text() + no_legend() +
    labs(x='First divergent bin', y='Modality')

ggsave('plots/paper/fig3/fig3_ctx_traject_priming_quant_boxplot.pdf', width=5, height=2.4, units='cm')
ggsave('plots/paper/fig3/fig3_ctx_traject_priming_quant_boxplot.png', width=5, height=2.4, units='cm')





#### Cluster heatmap ####
gene_clusters <- read_tsv('data_/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')

H3K27ac_ctx <- read_rds('data_/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data_/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data_/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data_/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')

mark_feats <- intersect(rownames(H3K27ac_ctx[['cRNA']]), rownames(H3K27me3_ctx[['cRNA']])) %>% intersect(rownames(H3K4me3_ctx[['cRNA']]))

H3K27ac_cluster_detect <- H3K27ac_ctx@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]
H3K27me3_cluster_detect <- H3K27me3_ctx@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]
H3K4me3_cluster_detect <- H3K4me3_ctx@assays$cRNA@misc$summary$pt_bins[as.character(1:20), mark_feats]
rna_cluster_detect <- rna_ctx@assays$RNA@misc$summary$pt_bins[as.character(1:20),]

H3K27ac_maxdetect_feats <- colMaxs(H3K27ac_cluster_detect)
H3K27me3_maxdetect_feats <- colMaxs(H3K27me3_cluster_detect)
H3K4me3_maxdetect_feats <- colMaxs(H3K4me3_cluster_detect)
rna_maxdetect_feats <- colMaxs(rna_cluster_detect)

mark_detect_feats <- mark_feats[
    (H3K27ac_maxdetect_feats > 0.02) & (H3K27me3_maxdetect_feats > 0.02) & (H3K4me3_maxdetect_feats > 0.02)
]

feats_use <- gene_clusters %>% 
    pull(feature) %>% intersect(mark_detect_feats)


#### Smooth with GAMs ####
x <- 1:20
H3K27ac_gams <- H3K27ac_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
rna_smooths <- rna_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27ac_smooths_scale <- H3K27ac_smooths %>% scale()
H3K27me3_smooths_scale <- H3K27me3_smooths %>% scale()
H3K4me3_smooths_scale <- H3K4me3_smooths %>% scale()
rna_smooths_scale <- rna_smooths %>% scale()

all_smooths_list <- map(purrr::set_names(colnames(rna_smooths)), function(n){
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

all_smooths_raw_list <- map(purrr::set_names(colnames(rna_smooths)), function(n){
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

set.seed(111)
dtw_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 10
)

dtw_kmeans_df <- dtw_kmeans$Cls %>% enframe('feature', 'dtw_clust')


#### Hclust per cluster ####
cluster_hclust_order <- map(set_names(sort(unique(dtw_kmeans_df$dtw_clust))), function(n){
    cl_genes <- filter(dtw_kmeans_df, dtw_clust==n)$feature
    cl_expr <- dtw_smooth_dist_mat[cl_genes, cl_genes]
    if (length(cl_genes)>2){
        cl_expr %>% as.dist() %>% hclust() %>% {.$labels[.$order]} %>% return()
    } else {
        return(cl_genes)
    }
})


#### Plot clusters #####
all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
        {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=1:20) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr')


plot_df <- all_gene_smooths_df %>% 
    inner_join(dtw_kmeans_df) %>%
    # inner_join(cluster_hclust_clusts) %>%
    dplyr::group_by(feature, modality) %>% 
    dplyr::mutate(
        expr01=scale01(expr), 
        pt_bin_mod=paste0(modality, pt_bin),
        feature=factor(feature, levels=purrr::reduce(cluster_hclust_order, c)),
        # dtw_clust=factor(dtw_clust, levels=c(7,1,6,10,4,9,8,2,3,5))
    )


plot_df %>% write_tsv('data/trajectories/ctx/all_mod_4m_ctx_clusters.tsv')


ggplot(plot_df, aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(dtw_clust~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100)) +
    theme(
        panel.border = element_blank()
    )
ggsave('plots/paper/fig3/fig3_ctx_2m_gam_smooth_10clust_labels_heatmap.pdf', width=15, height=20, limitsize=FALSE)
ggsave('plots/paper/fig3/fig3_ctx_2m_gam_smooth_10clust_labels_heatmap.png', width=15, height=20, limitsize=FALSE)



ggplot(plot_df, aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(dtw_clust~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100)) +
    article_text() +
    no_y_text() +
    theme(
        panel.border = element_blank(),
        panel.spacing = unit(0.1, 'cm')
    ) + no_legend() +
    labs(x='Pseudotime bin', y='Feature')
ggsave('plots/paper/fig3/fig3_ctx_2m_gam_smooth_10clust_heatmap.pdf', width=6, height=8, unit='cm')
ggsave('plots/paper/fig3/fig3_ctx_2m_gam_smooth_10clust_heatmap.pdf', width=6, height=8, unit='cm')



#### Filter only for neuronal expression and cluster with H3K27ac ####
H3K27ac_neuronal_idx <- H3K27ac_cluster_detect[18:20, ] %>% colMaxs() %>% {which(.>0.05)}
H3K27ac_neuronal_feats <- colnames(H3K27ac_cluster_detect)[H3K27ac_neuronal_idx]

H3K4me3_neuronal_idx <- H3K4me3_cluster_detect[18:20, ] %>% colMaxs() %>% {which(.>0.05)}
H3K4me3_neuronal_feats <- colnames(H3K4me3_cluster_detect)[H3K4me3_neuronal_idx]

rna_cluster_detect_scaled <- rna_cluster_detect %>% {.[,colMaxs(.)>0]} %>% apply(2,scale01) %>% Matrix::Matrix(sparse=T)
rna_neuronal_idx <- rna_cluster_detect_scaled[18:20, ] %>% colMaxs() %>% {which(.>0.9)}

feats_use <- colnames(rna_cluster_detect_scaled)[rna_neuronal_idx] %>% intersect(H3K27ac_neuronal_feats) %>% intersect(H3K4me3_neuronal_feats)


#### Smooth with GAMs ####
x <- 1:20
H3K27ac_gams <- H3K27ac_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_cluster_detect[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
rna_smooths <- rna_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27ac_smooths_scale <- H3K27ac_smooths %>% scale()
H3K27me3_smooths_scale <- H3K27me3_smooths %>% scale()
H3K4me3_smooths_scale <- H3K4me3_smooths %>% scale()
rna_smooths_scale <- rna_smooths %>% scale()

all_smooths_list <- map(purrr::set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths_scale[,n],
            # H3K27me3=H3K27me3_smooths_scale[,n],
            # H3K4me3=H3K4me3_smooths_scale[,n],
            rna=rna_smooths_scale[,n]
        )
    )
})

all_smooths_raw_list <- map(purrr::set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths[,n],
            # H3K27me3=H3K27me3_smooths[,n],
            # H3K4me3=H3K4me3_smooths[,n],
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

gene_order <- dtw_smooth_dist_mat %>% as.dist() %>% hclust() %>% {.$labels[.$order]} %>% return()


#### Plot clusters #####
all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
        # {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
        {colnames(.) <- c('H3K27ac', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=1:20) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr')


plot_df <- all_gene_smooths_df %>% 
    # inner_join(dtw_kmeans_df) %>%
    # inner_join(cluster_hclust_clusts) %>%
    dplyr::group_by(feature, modality) %>% 
    dplyr::mutate(
        expr01=scale01(expr), 
        pt_bin_mod=paste0(modality, pt_bin),
        feature=factor(feature, levels=gene_order),
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



#### Quantify lag for those genes ####
genes_test <- colnames(rna_cluster_detect_scaled)[rna_neuronal_idx] %>% intersect(H3K27ac_neuronal_feats)

test_pt_bins <- function(srt_obj, features=NULL, assay='RNA_bin', family='binomial'){
    srt_obj@active.assay <- assay
    srt_obj <- DietSeurat(srt_obj, assays=assay)
    
    test_result <- Pando::map_par(set_names(as.character(2:10)), function(x){
        srt_obj <- subset(srt_obj, pt_bins2%in%c('1', x))
        srt_obj$test_var <- srt_obj$pt_bins2 == x
        
        detection_rates <- Pando::aggregate_matrix(t(srt_obj[[assay]]@counts), srt_obj$test_var)
        detrate_df <- tibble(
            feature = colnames(detection_rates),
            detect_self = as.numeric(detection_rates['TRUE', ]),
            detect_other = as.numeric(detection_rates['FALSE', ])
        ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
        
        featurs_use <- intersect(features, rownames(srt_obj[[assay]]))
        
        test_df <- lr_de(
            object = srt_obj,
            test_var = 'test_var',
            covariates = c('nCount_RNA'),
            family = family,
            slot = 'counts',
            assay = assay,
            features_use = featurs_use,
            parallel = T
        )
        test_results <- inner_join(test_df, detrate_df) %>% 
            return()
    }, parallel=F) %>% bind_rows(.id='pt_bin')
    return(test_result)
}

library(doParallel)
registerDoParallel(20)

rna_binom_pt_test <- test_pt_bins(rna_ctx, features=genes_test, assay='RNA_bin')
H3K27ac_binom_pt_test <- test_pt_bins(H3K27ac_ctx, features=genes_test, assay='RNA_bin')
H3K27me3_binom_pt_test <- test_pt_bins(H3K27me3_ctx, features=genes_test, assay='RNA_bin')
H3K4me3_binom_pt_test <- test_pt_bins(H3K4me3_ctx, features=genes_test, assay='RNA_bin')


all_lr_pt_de <- bind_rows('RNA'=rna_binom_pt_test, 'H3K27ac'=H3K27ac_binom_pt_test, 'H3K27me3'=H3K27me3_binom_pt_test, 'H3K4me3'=H3K4me3_binom_pt_test, .id='modality') %>% 
    group_by(modality) %>% mutate(padj=p.adjust(pval, 'fdr'))

all_lr_pt_de %>% write_tsv('data_/trajectories/ctx/ctx_all_neuonal_pt_lr_binom_de.tsv')
all_lr_pt_de <- read_tsv('data_/trajectories/ctx/ctx_all_neuonal_pt_lr_binom_de.tsv')

all_lr_pt_diff <- all_lr_pt_de %>% 
    filter(pval<0.01, log_dr>log2(1.5), modality!='H3K27me3') %>% 
    group_by(feature) %>% 
    filter(length(unique(modality))>1) %>%
    group_by(modality, feature) %>% 
    dplyr::summarize(min_bin=min(as.numeric(pt_bin)), min_dr=log_dr[which.min(pt_bin)]) %>% 
    dplyr::mutate(signed_bin=sign(min_dr)*min_bin) 

ggplot(all_lr_pt_diff, aes(min_bin, modality, fill=modality, color=modality)) +
    geom_boxplot(linewidth=0.2, color='black', alpha=0.5, outlier.shape=NA) +
    geom_quasirandom(size=0.1) +
    scale_fill_manual(values=modality_colors) + 
    scale_color_manual(values=modality_colors) + 
    scale_x_continuous(limits=c(0,10)) + 
    article_text() + no_legend() +
    labs(x='First divergent bin', y='Modality')

ggsave('plots/paper/fig3/rev2_lag_quantification_neuron_genes_boxplot.pdf', width=6, height=3, units='cm')

 H3K27ac_lag <- all_lr_pt_diff %>% 
    filter(modality!='H3K4me3') %>% 
    select(modality, feature, min_bin) %>% 
    pivot_wider(names_from='modality', values_from='min_bin') %>% 
    mutate(bin_diff=RNA-H3K27ac) %>% filter(!is.na(bin_diff)) %>% select(feature, bin_diff)

H3K4me3_lag <- all_lr_pt_diff %>% 
    filter(modality!='H3K27ac') %>% 
    select(modality, feature, min_bin) %>% 
    pivot_wider(names_from='modality', values_from='min_bin') %>% 
    mutate(bin_diff=RNA-H3K4me3) %>% filter(!is.na(bin_diff)) %>% select(feature, bin_diff)

plot_df <- bind_rows('H3K27ac'=H3K27ac_lag, 'H3K4me3'=H3K4me3_lag, .id = 'mark')

ggplot(plot_df, aes(bin_diff, mark, fill=mark)) + 
    geom_density_ridges() +
    geom_vline(xintercept=c(-2,2)) +
    scale_fill_manual(values=modality_colors) +
    article_text() + no_legend() +
    labs(x='Pseudotime bin lag relative to RNA', y='Mark') 

ggsave('plots/paper/fig3/rev2_lag_quantification_neuron_genes_dridges.pdf', width=6, height=3, units='cm')


#### GO enrichment ####
go_enrich <- function(features, background, p_adjust_method='fdr', ontology='BP', ...){
    require(clusterProfiler)
    require(org.Hs.eg.db)
    
    all_features <- bitr(background, fromType='SYMBOL', toType=c('ENTREZID'), OrgDb=org.Hs.eg.db) %>% 
        as_tibble()
    
    gene_ids <- filter(all_features, SYMBOL%in%features)
    
    ggo <- enrichGO(
        gene = unique(gene_ids$ENTREZID), 
        universe = unique(all_features$ENTREZID),
        OrgDb = org.Hs.eg.db,
        readable = T,
        pAdjustMethod = p_adjust_method,
        ont = ontology,
        ...
    )
    res <- as_tibble(ggo@result) %>% 
        rename('count'=Count) %>% 
        mutate(
            fg_in = as.numeric(str_replace(GeneRatio, '(\\d+)/(\\d+)', '\\1')),
            fg_out = as.numeric(str_replace(GeneRatio, '(\\d+)/(\\d+)', '\\2')),
            bg_in = as.numeric(str_replace(BgRatio, '(\\d+)/(\\d+)', '\\1')),
            bg_out = as.numeric(str_replace(BgRatio, '(\\d+)/(\\d+)', '\\2')),
            odds_ratio = (fg_in/fg_out) / (bg_in/bg_out),
            log2_odds_ratio = log2(odds_ratio)
        )
    
    return(res)
}


feats_use <- colnames(rna_cluster_detect_scaled)[rna_neuronal_idx] %>% intersect(H3K27ac_neuronal_feats)
neuron_feats <- colnames(rna_cluster_detect_scaled)[rna_neuronal_idx]
all_feats <- colnames(rna_cluster_detect_scaled)

primed_group <- plot_df %>% 
    filter(bin_diff>=2) %>% pull(feature)
write(primed_group, 'data_/trajectories/ctx/enrichment/primed_group_global.txt')

mid_group <- plot_df %>% 
    filter((bin_diff>(-2)) & (bin_diff<2)) %>% pull(feature) %>% setdiff(primed_group)
write(mid_group, 'data_/trajectories/ctx/enrichment/mid_group_global.txt')

rev_group <- plot_df %>% 
    filter(bin_diff<=(-2)) %>% pull(feature) %>% setdiff(primed_group) %>% setdiff(mid_group)
write(rev_group, 'data_/trajectories/ctx/enrichment/neg_group_global.txt')


# primed_enrich <- go_enrich(features=primed_group, background=feats_use)
# mid_enrich <- go_enrich(features=mid_group, background=feats_use)
# rev_enrich <- go_enrich(features=rev_group, background=feats_use)
#     
# primed_enrich_nrn <- go_enrich(features=primed_group, background=neuron_feats)
# mid_enrich_nrn <- go_enrich(features=mid_group, background=neuron_feats)
# rev_enrich_nrn <- go_enrich(features=rev_group, background=neuron_feats)
    
primed_enrich_glob <- go_enrich(features=primed_group, background=all_feats)
mid_enrich_glob <- go_enrich(features=mid_group, background=all_feats)
rev_enrich_glob <- go_enrich(features=rev_group, background=all_feats)
    
primed_enrich_glob %>% write_tsv('data_/trajectories/ctx/enrichment/primed_global_bg_goenrich.tsv')
mid_enrich_glob %>% write_tsv('data_/trajectories/ctx/enrichment/mid_global_bg_goenrich.tsv')
rev_enrich_glob %>% write_tsv('data_/trajectories/ctx/enrichment/neg_global_bg_goenrich.tsv')



# For the neuron genes
gene_clusters <- read_tsv('data_/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')
neuron_genes <- gene_clusters %>% filter(dtw_clust==27) %>% pull(feature) %>% unique()
neuron_goenrich <- go_enrich(features=neuron_genes, background=gene_clusters$feature)
neuron_goenrich %>% write_tsv('data_/trajectories/ctx/enrichment/neuron_cluster_goenrich.tsv')





#### TF motif enrichment ####
isects_annot <- read_tsv('data_/intersect/all_marks_intersect_matches_annot.tsv')
data(motif2tf)

H3K27ac_regions <- isects_annot %>% 
    filter(gene%in%primed_group) %>% 
    pull(H3K27ac) %>% unique()

H3K27ac_motif_enrich <- FindMotifs(
    object = marks$H3K27ac,
    features = H3K27ac_regions
) %>% as_tibble() %>% inner_join(motif2tf)

H3K27ac_motif_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/primed_gene_H3K27ac_peak_motif_enrich.tsv')

H3K27ac_regions <- isects_annot %>% 
    filter(gene%in%mid_group) %>% 
    pull(H3K27ac) %>% unique()

H3K27ac_motif_enrich <- FindMotifs(
    object = marks$H3K27ac,
    features = H3K27ac_regions
) %>% as_tibble() %>% inner_join(motif2tf)

H3K27ac_motif_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/mid_gene_H3K27ac_peak_motif_enrich.tsv')

H3K4me3_regions <- isects_annot %>% 
    filter(gene%in%primed_group) %>% 
    pull(H3K4me3) %>% unique()

H3K4me3_motif_enrich <- FindMotifs(
    object = marks$H3K4me3,
    features = H3K4me3_regions
) %>% as_tibble() %>% inner_join(motif2tf)

H3K4me3_motif_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/primed_gene_H3K4me3_peak_motif_enrich.tsv')

H3K4me3_regions <- isects_annot %>% 
    filter(gene%in%mid_group) %>% 
    pull(H3K4me3) %>% unique()

H3K4me3_motif_enrich <- FindMotifs(
    object = marks$H3K4me3,
    features = H3K4me3_regions
) %>% as_tibble() %>% inner_join(motif2tf)

H3K4me3_motif_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/mid_gene_H3K4me3_peak_motif_enrich.tsv')


#### Enrichment in merged regions ####
tf_motif_enrichment <- function(regions, motifs, background_regions, parallel=T){
    data(motif2tf)
    group_vec <- background_regions %in% regions
    map_par(colnames(motifs), function(mot){
        ctab <- table(motifs[background_regions, mot], group_vec)
        ftest <- try(fisher.test(ctab))
        if (any(class(ftest)=='try-error')){
            return(tibble())
        }
        return(tibble(
            motif = mot,
            pval = ftest$p.value,
            odds_ratio = ftest$estimate
        ))
    }, parallel=parallel) %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate(padj=p.adjust(pval, 'fdr')) %>%
        dplyr::mutate(logodds=log2(odds_ratio)) %>%
        dplyr::inner_join(motif2tf, multiple = "all") %>%
        return()
}

all_bg <- peak_meta$repr_region %>% unique() %>% StringToGRanges()
all_peaks_srt <- read_rds('data05/regions/all_regions_v2_max01_clust50_srt.rds')
all_motif_mat <- read_rds('data_/all_peaks/all_motifs.rds')

gene_ranges_use <- CollapseToLongestTranscript(gene_annot)
genes_primed <- gene_ranges_use[gene_ranges_use$gene_name%in%primed_group]
genes_mid <- gene_ranges_use[gene_ranges_use$gene_name%in%mid_group]

primed2peaks <- find_peaks_near_genes(all_bg, genes_primed, distance=4000)
primed_peaks <- rownames(primed2peaks)[rowMaxs(primed2peaks)>0] 
mid2peaks <- find_peaks_near_genes(all_bg, genes_mid, distance=4000)
mid_peaks <- rownames(mid2peaks)[rowMaxs(mid2peaks)>0] 

primed_tf_enrich <- tf_motif_enrichment(primed_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
primed_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/primed_tf_enrich.tsv')
primed_tf_enrich <- read_tsv('data_/trajectories/ctx/enrichment/primed_tf_enrich.tsv')

mid_tf_enrich <- tf_motif_enrichment(mid_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
mid_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/mid_tf_enrich.tsv')


primed2peaks <- find_peaks_near_genes(all_bg, genes_primed, distance=4000, only_tss=T)
primed_peaks <- rownames(primed2peaks)[rowMaxs(primed2peaks)>0] 
mid2peaks <- find_peaks_near_genes(all_bg, genes_mid, distance=4000, only_tss=T)
mid_peaks <- rownames(mid2peaks)[rowMaxs(mid2peaks)>0] 

primed_tf_enrich <- tf_motif_enrichment(primed_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
primed_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/primed_tf_promoter_enrich.tsv')

mid_tf_enrich <- tf_motif_enrichment(mid_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
mid_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/mid_tf_promoter_enrich.tsv')


# For the neuron genes
gene_clusters <- read_tsv('data_/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')
neuron_genes <- gene_clusters %>% filter(dtw_clust==27) %>% pull(feature) %>% unique()

neuron_gene_ranges <- gene_ranges_use[gene_ranges_use$gene_name%in%neuron_genes]
neuron2peaks <- find_peaks_near_genes(all_bg, neuron_gene_ranges, distance=4000, only_tss=T)
neuron_peaks <- rownames(neuron2peaks)[rowMaxs(neuron2peaks)>0] 

neuron_tf_enrich <- tf_motif_enrichment(neuron_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
neuron_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/neuron_cluster_tf_promoter_enrich.tsv')


neuron2peaks <- find_peaks_near_genes(all_bg, neuron_gene_ranges, distance=4000)
neuron_peaks <- rownames(neuron2peaks)[rowMaxs(neuron2peaks)>0] 

neuron_tf_enrich <- tf_motif_enrichment(neuron_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
neuron_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/neuron_cluster_tf_enrich.tsv')




#### Check NEUROD2 family tf expression ####
bhlh_tfs <- motif2tf %>% filter(origin=='JASPAR2020', family=='Tal-related factors') %>% pull(tf) %>% unique()
bhlh_mots <- motif2tf %>% filter(origin=='JASPAR2020', family=='Tal-related factors') %>% pull(motif) %>% unique()
feature_plot(rna, features=bhlh_tfs, order=T, reduction='cssumap')
ggsave('plots/paper/bhlh_fam_tf_expression_feature_umap.png', width=10, height=15)

MotifPlot(marks$H3K27me3, motifs = bhlh_mots)
ggsave('plots/paper/bhlh_fam_tf_motif_pwm.png', width=10, height=7)
ggsave('plots/paper/bhlh_fam_tf_motif_pwm.pdf', width=10, height=7)

#####


#### Other plots for the supp figure ####

# nCounts over pseudotime
# Jitter genes + barplots for exemplary genes

H3K27ac_meta <- as_tibble(H3K27ac_ctx@meta.data, rownames='cell')
H3K4me3_meta <- as_tibble(H3K4me3_ctx@meta.data, rownames='cell')
H3K27me3_meta <- as_tibble(H3K27me3_ctx@meta.data, rownames='cell')
rna_meta <- as_tibble(rna_ctx@meta.data, rownames='cell')

all_meta <- bind_rows(
    'H3K27ac'=H3K27ac_meta,
    'H3K4me3'=H3K4me3_meta,
    'H3K27me3'=H3K27me3_meta,
    'rna'=rna_meta,
    .id='modality'
)

ggplot(all_meta, aes(ctx_pt, log10(nCount_RNA))) +
    geom_point(size=0.01, shape=16) +
    facet_grid(modality~.) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='Pseudotime', y='Log10(counts)')
ggsave('plots/paper/fig3/fig3_ctx_2m_ncount_pt_scatter.pdf', width=5, height=4, unit='cm')
ggsave('plots/paper/fig3/fig3_ctx_2m_ncount_pt_scatter.png', width=5, height=4, unit='cm')


ggplot(all_meta, aes(pt_bins, log10(nCount_RNA))) +
    geom_bar(stat='summary', color='black', linewidth=0.1, fill='grey') +
    # geom_errorbar(stat='summary') +
    facet_grid(modality~.) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='Pseudotime', y='Log10(counts)')
ggsave('plots/paper/fig3/fig3_ctx_2m_ncount_pt_bar.pdf', width=5, height=4, unit='cm')
ggsave('plots/paper/fig3/fig3_ctx_2m_ncount_pt_bar.png', width=5, height=4, unit='cm')



genes_plot <- c('GRIA2', 'STMN2', 'NEUROD2', 'SOX2', 'VIM')

H3K27ac_expr_plot <- H3K27ac_ctx[['RNA']]@data[genes_plot, ] %>%
    as_tibble(rownames='feature') %>%
    pivot_longer(!feature, names_to='cell', values_to='expr')

H3K27me3_expr_plot <- H3K27me3_ctx[['RNA']]@data[genes_plot, ] %>%
    as_tibble(rownames='feature') %>%
    pivot_longer(!feature, names_to='cell', values_to='expr')

H3K4me3_expr_plot <- H3K4me3_ctx[['RNA']]@data[genes_plot, ] %>%
    as_tibble(rownames='feature') %>%
    pivot_longer(!feature, names_to='cell', values_to='expr')

rna_expr_plot <- rna_ctx[['RNA']]@data[genes_plot, ] %>%
    as_tibble(rownames='feature') %>%
    pivot_longer(!feature, names_to='cell', values_to='expr')

expr_meta <- bind_rows('H3K27ac'=H3K27ac_expr_plot, 'H3K27me3'=H3K27me3_expr_plot, 'H3K4me3'=H3K4me3_expr_plot, 'rna'=rna_expr_plot, .id='modality') %>% 
    mutate(feature=factor(feature, levels=genes_plot)) %>% 
    inner_join(all_meta)

plot_df_list <- expr_meta %>% 
    group_by(modality) %>% 
    group_split()




gyylorrd2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='YlOrRd')[3:9]), bias=bias)(100))
}

gyylgn2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='YlGn')[3:9]), bias=bias)(100))
}

gybu2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='Blues')[3:9]), bias=bias)(100))
}

gyrdpu2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='RdPu')[3:9]), bias=bias)(100))
}

# pals::brewer.rdpu(100)[0:70] pals::brewer.blues(100)[0:70] pals::brewer.ylgn(100)[0:80]

modality_scales <- list(
    H3K27ac = gyylgn2(1.5),
    rna = gyylorrd2(1.5),
    H3K27me3 = gybu2(1.2),
    H3K4me3 = gyrdpu2(1.5)
)

plots <- map(plot_df_list, function(plot_df){
    mod <- plot_df$modality[1]
    color_scale <- 
    ggplot(arrange(plot_df, expr), aes(ctx_pt, feature, color=expr)) +
        geom_jitter(size=0.1, height=0.3) +
        facet_grid(modality~.) +
        # scale_color_gradientn(colors=grad(pals::brewer.ylgnbu, 1.5)) +
        scale_color_gradientn(colors=modality_scales[[mod]]) +
        article_text() +
        labs(x='Pseudotime', y='Feature', color='Log-normalized counts')
})
wrap_plots(plots, ncol=1)

ggsave('plots/paper/fig3/fig3_ctx_2m_pt_expr_jitter.png', width=10, height=10, units='cm', bg='white')
ggsave('plots/paper/fig3/fig3_ctx_2m_pt_expr_jitter.pdf', width=10, height=10, units='cm', bg='white')



H3K27ac_expr_plot <- H3K27ac_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), genes_plot] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bin', values_to='expr')

H3K27me3_expr_plot <- H3K27me3_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), genes_plot] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bin', values_to='expr')

H3K4me3_expr_plot <- H3K4me3_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), genes_plot] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bin', values_to='expr')

rna_expr_plot <- rna_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), genes_plot] %>% t() %>% 
    as_tibble(rownames='feature') %>% 
    pivot_longer(!feature, names_to='pt_bin', values_to='expr')

expr_meta <- bind_rows('H3K27ac'=H3K27ac_expr_plot, 'H3K27me3'=H3K27me3_expr_plot, 'H3K4me3'=H3K4me3_expr_plot, 'rna'=rna_expr_plot, .id='modality') %>% 
    mutate(feature=factor(feature, levels=genes_plot)) 


plot_df <- expr_meta 

ggplot(plot_df, aes(as.numeric(pt_bin), expr, color=modality, fill=modality)) +
    geom_bar(stat='identity', position='dodge', color='black', width=0.9, linewidth=0.1, alpha=0.5) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    scale_color_manual(values=mark_colors) +
    scale_fill_manual(values=mark_colors) +
    facet_grid(modality~feature, scales='free') +
    article_text() + no_legend() +
    labs(x='Pseudotime bins', y='Detection rate')

ggsave('plots/paper/fig3/fig3_ctx_2m_pt_detection_bar.png', width=12, height=9, units='cm', bg='white')
ggsave('plots/paper/fig3/fig3_ctx_2m_pt_detection_bar.pdf', width=12, height=9, units='cm', bg='white')




#### Check for primed genes DE over timepoints ####
age_de <- read_tsv('data_/results/diff_expression/RNA_DE_age.tsv')
links_grn <- read_rds('data_/GRN/RNA_ATAC_pseudocells_pando_links_srt.rds')
nd2_coef <- read_tsv('data_/GRN/neuro_grn_coefs.tsv')
nd2_coef_status <- read_tsv('data_/GRN/neuro_grn_pt_status.tsv')

nd2_only_coef <- coef(links_grn) %>% 
    filter(padj<1e-4, corr>0.2, estimate>0, ((target=='NEUROD2' & target%in%nd2_coef_status$target) | tf=='NEUROD2')) 

out_genes <- nd2_only_coef %>% filter(tf=='NEUROD2') %>% pull(target) %>% intersect(colnames(H3K27ac_clusters))
in_genes <- nd2_only_coef %>% filter(target=='NEUROD2') %>% pull(tf) %>% intersect(colnames(H3K27ac_clusters))

out_genes_de <- age_de %>% filter(gene%in%out_genes)
in_genes_de <- age_de %>% filter(gene%in%in_genes)
nd2_de <- age_de %>% filter(gene=='NEUROD2')

grn_gene_de <- bind_rows('in'=in_genes_de, 'out'=out_genes_de, 'NEUROD2'=nd2_de, .id='status')

names(age_colors2) <- c('EB', '15d', '35d', '60d', '63d', '65d', '4mo', '8mo')
plot_df <- grn_gene_de %>% 
    mutate(age=factor(cluster, levels=names(age_colors2)[1:8])) %>% 
    filter(!is.na(age), avg_log2FC>0.1) %>% 
    group_by(gene) %>% 
    filter(as.numeric(age)==min(as.numeric(age)))

ggplot(arrange(plot_df, status=='in'), aes(age, avg_log2FC, fill=status)) +
    geom_point(position='dodge', color='black', shape=21, size=3) 
    
ggplot(plot_df, aes(age, fill=status)) +
    geom_bar(alpha=0.5, color='black') +
    labs(x='Time point', y='# genes first upregulated in this timepoint')
    

















