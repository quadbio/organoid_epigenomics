source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(destiny)

rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist
Matrix <- Matrix::Matrix

setwd('~/projects/cutntag/')


rank_cols <- function(x){
    xr <- presto::rank_matrix(t(x))$X_ranked
    colnames(xr) <- colnames(x)
    rownames(xr) <- rownames(x)
    return(xr)
}


#### Read data ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')

rna_pt_meta <- read_tsv('data/RNA/cellrank/RNA_full_cellrank_probs.tsv') %>% 
    dplyr::rename('cell'=1) %>% 
    select(cell, velocity_pseudotime, pseudotime_ranks) %>% 
    column_to_rownames('cell')

rna <- AddMetaData(rna, rna_pt_meta)

cluster_srt <- read_rds('data/RNA/all_RNA_marks_combined_clusters_srt.rds')
cluster_graph <- read_rds('data/RNA/all_RNA_cluster_graph.rds')


#### Prep data: dienscphalic trajectory ####
marks <- map(marks, function(x){
    x$dien_traject <- (
        (x$lineage == 'dien') &
            (x$stage != 'retina') &
            (x$age != '8mo')
    )
    return(x)
})

rna$dien_traject <- (
    (rna$lineage == 'dien') &
        (rna$stage != 'retina') &
        (rna$age != '8mo')
)

dim_plot(rna, group.by='dien_traject', reduction='cssumap')



#### H3K27ac ####
H3K27ac_dien <- marks$H3K27ac %>% subset(dien_traject==TRUE)

H3K27ac_dien@active.assay <- 'peaks'
H3K27ac_dien <- H3K27ac_dien %>%
    RunTFIDF() %>% 
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

lsi <- CreateAssayObject(t(H3K27ac_dien[['lsi']]@cell.embeddings))
H3K27ac_dien[['lsi']] <- NULL
H3K27ac_dien[['lsi']] <- lsi
H3K27ac_dien <- H3K27ac_dien %>% ScaleData(assay='lsi', vars.to.regress=c('log_peak_counts', 'age'))
lsi <- t(H3K27ac_dien[['lsi']]@scale.data)
H3K27ac_dien[['lsi']] <- NULL
H3K27ac_dien[['lsi']] <- CreateDimReducObject(lsi, key='LSI_', assay='peaks')

H3K27ac_dien <- H3K27ac_dien %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K27ac_dien)
ElbowPlot(H3K27ac_dien, reduction='lsi')

dim_plot(H3K27ac_dien, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27ac_dien@active.assay <- 'cRNA'
feature_plot(H3K27ac_dien, label=T, features=c('STMN2', 'VIM', 'NES', 'GRIA2', 'log_peak_counts'), order=T, pt.size=2)

H3K27ac_pca <- H3K27ac_dien[['lsi']]@cell.embeddings[,2:10]
H3K27ac_diffmap <- DiffusionMap(H3K27ac_pca, verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5]))

H3K27ac_dien[['diffmap']] <- CreateDimReducObject(
    H3K27ac_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27ac_dien <- AddMetaData(H3K27ac_dien, H3K27ac_dc_df)

H3K27ac_dien %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'NES', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5'), order=T)

mexp <- as.numeric(H3K27ac_dien[['cRNA']]['STMN2'])
plot(H3K27ac_dien$DC1, mexp, col=factor(H3K27ac_dien$age))

mexp <- as.numeric(H3K27ac_dien[['cRNA']]['VIM'])
plot(H3K27ac_dien$DC1, mexp, col=factor(H3K27ac_dien$age))

plot(H3K27ac_dien$DC1, H3K27ac_dien$log_peak_counts, col=factor(H3K27ac_dien$age))

H3K27ac_dien$dien_pt <- rank(-H3K27ac_dien$DC1) / max(rank(H3K27ac_dien$DC1))
H3K27ac_dien %>% write_rds('data/trajectories/dien/H3K27ac_dien_dpt_lsi_regress_srt.rds')



#### H3K4me3 ####
H3K4me3_dien <- marks$H3K4me3 %>% subset(dien_traject==TRUE)

H3K4me3_dien@active.assay <- 'peaks'
H3K4me3_dien <- H3K4me3_dien %>%
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

lsi <- CreateAssayObject(t(H3K4me3_dien[['lsi']]@cell.embeddings))
H3K4me3_dien[['lsi']] <- NULL
H3K4me3_dien[['lsi']] <- lsi
H3K4me3_dien <- H3K4me3_dien %>% ScaleData(assay='lsi', vars.to.regress=c('log_peak_counts', 'age'))
lsi <- t(H3K4me3_dien[['lsi']]@scale.data)
H3K4me3_dien[['lsi']] <- NULL
H3K4me3_dien[['lsi']] <- CreateDimReducObject(lsi, key='LSI_', assay='peaks')

H3K4me3_dien <- H3K4me3_dien %>% 
    RunUMAP(dims=c(2:10), reduction='lsi')

DepthCor(H3K4me3_dien)
ElbowPlot(H3K4me3_dien, reduction='lsi')

dim_plot(H3K4me3_dien, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K4me3_dien@active.assay <- 'cRNA'
feature_plot(H3K4me3_dien, features=c('STMN2', 'VIM', 'NES', 'GRIA2', 'RSPO3', 'log_peak_counts'), order=T, pt.size=2)

H3K4me3_pca <- H3K4me3_dien[['lsi']]@cell.embeddings[,2:10]
H3K4me3_diffmap <- DiffusionMap(H3K4me3_pca, verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5]))

H3K4me3_dien[['diffmap']] <- CreateDimReducObject(
    H3K4me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K4me3_dien <- AddMetaData(H3K4me3_dien, H3K4me3_dc_df)

H3K4me3_dien %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'log_peak_counts'), order=T)

mexp <- as.numeric(H3K4me3_dien[['cRNA']]['STMN2'])
plot(H3K4me3_dien$DC1, mexp, col=factor(H3K4me3_dien$age))

mexp <- as.numeric(H3K4me3_dien[['cRNA']]['VIM'])
plot(H3K4me3_dien$DC1, mexp, col=factor(H3K4me3_dien$age))

plot(H3K4me3_dien$DC1, H3K4me3_dien$log_peak_counts, col=factor(H3K4me3_dien$age))

H3K4me3_dien$dien_pt <- rank(H3K4me3_dien$DC1) / max(rank(H3K4me3_dien$DC1))
H3K4me3_dien %>% write_rds('data/trajectories/dien/H3K4me3_dien_dpt_lsi_regress_srt.rds')


#### H3K27me3 ####
H3K27me3_dien <- marks$H3K27me3 %>% subset(dien_traject==TRUE)

H3K27me3_dien@active.assay <- 'peaks'
H3K27me3_dien <- H3K27me3_dien %>%
    RunTFIDF(assay='peaks') %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD(assay='peaks')

lsi <- t(H3K27me3_dien[['lsi']]@cell.embeddings)
H3K27me3_dien[['lsi']] <- NULL
H3K27me3_dien[['lsi']] <- CreateAssayObject(lsi)
H3K27me3_dien <- H3K27me3_dien %>% ScaleData(assay='lsi', vars.to.regress=c('log_peak_counts', 'age'))
lsi <- H3K27me3_dien[['lsi']]@scale.data
H3K27me3_dien[['lsi']] <- NULL
H3K27me3_dien[['lsi']] <- CreateDimReducObject(t(lsi), key='LSI_', assay='peaks')

H3K27me3_dien <- H3K27me3_dien %>% 
    RunUMAP(dims=2:15, reduction='lsi')

DepthCor(H3K27me3_dien, reduction='lsi')
ElbowPlot(H3K27me3_dien, reduction='lsi')

dim_plot(H3K27me3_dien, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27me3_dien@active.assay <- 'cRNA'
feature_plot(H3K27me3_dien, features=c('STMN2', 'NES', 'PAX6', 'SOX2', 'VIM', 'log_peak_counts'), order=T, pt.size=2)

H3K27me3_pca <- H3K27me3_dien[['lsi']]@cell.embeddings[,2:15]
H3K27me3_diffmap <- DiffusionMap(H3K27me3_pca, verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5]))

H3K27me3_dien[['diffmap']] <- CreateDimReducObject(
    H3K27me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27me3_dien <- AddMetaData(H3K27me3_dien, H3K27me3_dc_df)

H3K27me3_dien %>% 
    feature_plot(features=c('STMN2', 'VIM', 'RSPO3', 'SOX2', 'DC1', 'DC2', 'DC3', 'DC4', 'log_peak_counts'), order=T)

mexp <- as.numeric(H3K27me3_dien[['cRNA']]['STMN2'])
plot(H3K27me3_dien$DC2, mexp, col=factor(H3K27me3_dien$age))

mexp <- as.numeric(H3K27me3_dien[['cRNA']]['VIM'])
plot(H3K27me3_dien$DC2, mexp, col=factor(H3K27me3_dien$age))

plot(H3K27me3_dien$DC2, H3K27me3_dien$log_peak_counts, col=factor(H3K27me3_dien$age))

H3K27me3_dien$dien_pt <- rank(H3K27me3_dien$DC2) / max(rank(H3K27me3_dien$DC2))
H3K27me3_dien %>% write_rds('data/trajectories/dien/H3K27me3_dien_dpt_lsi_regress_srt.rds')


#### Re-compute PCA for dien trajectory ####
rna_dien <- rna %>% subset(dien_traject==TRUE)

rna_dien <- rna_dien %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

ElbowPlot(rna_dien)

rna_dien <- rna_dien %>% 
    RunUMAP(reduction='css', dims=1:ncol(rna_dien[['css']]))

rna_dien %>% dim_plot(label=T, group.by=c('clusters', 'celltype_jf', 'age'))
rna_dien %>% feature_plot(features=c('STMN2', 'VIM', 'RSPO3', 'SOX2', 'pseudotime_ranks'), order=T)

rna %>% feature_plot(features=c('STMN2', 'VIM', 'RSPO3', 'SOX2', 'pseudotime_ranks'), order=T, reduction='cssumap')

rna_dien$dien_pt <- rank(rna_dien$pseudotime_ranks) / max(rank(rna_dien$pseudotime_ranks))
rna_dien %>% write_rds('data/trajectories/dien/RNA_dien_dpt_srt.rds')


dien_de <- de(rna_dien, groups = 'state')
vulcano_plot(dien_de)




#### Subsample modalities to minimal # cells for each timepoint ####
H3K27ac_dien <- read_rds('data/trajectories/dien/H3K27ac_dien_dpt_lsi_regress_srt.rds')
H3K4me3_dien <- read_rds('data/trajectories/dien/H3K4me3_dien_dpt_lsi_regress_srt.rds')
H3K27me3_dien <- read_rds('data/trajectories/dien/H3K27me3_dien_dpt_lsi_regress_srt.rds')
rna_dien <- read_rds('data/trajectories/dien/RNA_dien_dpt_srt.rds')

H3K27ac_dien$age <- ifelse(H3K27ac_dien$age=='128d', '4mo', H3K27ac_dien$age)
H3K27me3_dien$age <- ifelse(H3K27me3_dien$age=='128d', '4mo', H3K27me3_dien$age)
H3K4me3_dien$age <- ifelse(H3K4me3_dien$age=='128d', '4mo', H3K4me3_dien$age)

ncells_table <- bind_rows(
    table(H3K27ac_dien$age),
    table(H3K27me3_dien$age),
    table(H3K4me3_dien$age),
    table(rna_dien$age)
)

ncells_min <- colMins(as.matrix(ncells_table))
ncells_min[ncells_min<100] <- 100
names(ncells_min) <- colnames(ncells_table)

subs_list <- map(c('H3K27ac'=H3K27ac_dien, 'H3K4me3'=H3K4me3_dien, 'H3K27me3'=H3K27me3_dien, 'RNA'=rna_dien), function(srt){
    meta <- srt@meta.data %>% as_tibble(rownames='cell') %>% 
        group_by(age) %>% group_split() %>% 
        map(~sample_n(.x, size=min(ncells_min[.x$age[1]], nrow(.x)))) %>% 
        bind_rows()
    return(subset(srt, cells=meta$cell))
})

H3K27ac_dien_subs <- subs_list$H3K27ac
H3K27me3_dien_subs <- subs_list$H3K27me3
H3K4me3_dien_subs <- subs_list$H3K4me3
rna_dien_subs <- subs_list$RNA

# Re-rank pseudotime
H3K27ac_dien_subs$dien_pt <- rank(H3K27ac_dien_subs$dien_pt) / max(rank(H3K27ac_dien_subs$dien_pt))
H3K27me3_dien_subs$dien_pt <- rank(H3K27me3_dien_subs$dien_pt) / max(rank(H3K27me3_dien_subs$dien_pt))
H3K4me3_dien_subs$dien_pt <- rank(H3K4me3_dien_subs$dien_pt) / max(rank(H3K4me3_dien_subs$dien_pt))
rna_dien_subs$dien_pt <- rank(rna_dien_subs$dien_pt) / max(rank(rna_dien_subs$dien_pt))


#### Make pseudotime bins ####
H3K27ac_dien_subs$pt_bins <- as.numeric(cut(H3K27ac_dien_subs$dien_pt, 20, labels=1:20))
H3K27me3_dien_subs$pt_bins <- as.numeric(cut(H3K27me3_dien_subs$dien_pt, 20, labels=1:20))
H3K4me3_dien_subs$pt_bins <- as.numeric(cut(H3K4me3_dien_subs$dien_pt, 20, labels=1:20))
rna_dien_subs$pt_bins <- as.numeric(cut(rna_dien_subs$dien_pt, 20, labels=1:20))

H3K27ac_dien_subs <- Pando::aggregate_assay(H3K27ac_dien_subs, assay='cRNA', group_name='pt_bins')
H3K27me3_dien_subs <- Pando::aggregate_assay(H3K27me3_dien_subs, assay='cRNA', group_name='pt_bins')
H3K4me3_dien_subs <- Pando::aggregate_assay(H3K4me3_dien_subs, assay='cRNA', group_name='pt_bins')
rna_dien_subs <- Pando::aggregate_assay(rna_dien_subs, assay='RNA', group_name='pt_bins')

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
ggsave('plots/trajectories/dien_subs_age_dist_bar.png', width=8, height=5)



#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_clusters <- H3K27me3_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_clusters <- H3K4me3_dien_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
rna_clusters <- rna_dien_subs@assays$RNA@misc$summary$pt_bins[as.character(1:20), ]

genes_plot <- c('RSPO2', 'NEUROD2', 'WLS', 'SOX2', 'VIM', 'STMN2', 'GRIA2', 'DCX')

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

ggsave('plots/trajectories/dien_subs_gene_expr_line.png', width=10, height=6)


H3K27ac_dien_subs %>% write_rds('data/trajectories/dien/H3K27ac_dien_subs_dpt_lsi_regress_srt.rds')
H3K4me3_dien_subs %>% write_rds('data/trajectories/dien/H3K4me3_dien_subs_dpt_lsi_regress_srt.rds')
H3K27me3_dien_subs %>% write_rds('data/trajectories/dien/H3K27me3_dien_subs_dpt_lsi_regress_srt.rds')
rna_dien_subs %>% write_rds('data/trajectories/dien/RNA_dien_subs_dpt_srt.rds')



#### Annotate neurons ####

H3K27ac_dien <- read_rds('data/trajectories/dien/H3K27ac_dien_dpt_lsi_regress_srt.rds')
H3K4me3_dien <- read_rds('data/trajectories/dien/H3K4me3_dien_dpt_lsi_regress_srt.rds')
H3K27me3_dien <- read_rds('data/trajectories/dien/H3K27me3_dien_dpt_lsi_regress_srt.rds')
rna_dien <- read_rds('data/trajectories/dien/RNA_dien_dpt_srt.rds')

rna_dien$pt_bins <- as.numeric(cut(rna_dien$dien_pt, 20, labels=1:20))
meta <- rna_dien@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K27ac_dien$pt_bins <- as.numeric(cut(H3K27ac_dien$dien_pt, 20, labels=1:20))
meta <- H3K27ac_dien@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K27me3_dien$pt_bins <- as.numeric(cut(H3K27me3_dien$dien_pt, 20, labels=1:20))
meta <- H3K27me3_dien@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K4me3_dien$pt_bins <- as.numeric(cut(H3K4me3_dien$dien_pt, 20, labels=1:20))
meta <- H3K4me3_dien@meta.data %>% 
    as_tibble(rownames='cell')

ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()


rna_dien$neuron <- rna_dien$pt_bins > 17
H3K27ac_dien$neuron <- H3K27ac_dien$pt_bins > 17
H3K27me3_dien$neuron <- H3K27me3_dien$pt_bins > 17
H3K4me3_dien$neuron <- H3K4me3_dien$pt_bins > 17

H3K27ac_dien %>% write_rds('data/trajectories/dien/H3K27ac_dien_dpt_lsi_regress_srt.rds')
H3K4me3_dien %>% write_rds('data/trajectories/dien/H3K4me3_dien_dpt_lsi_regress_srt.rds')
H3K27me3_dien %>% write_rds('data/trajectories/dien/H3K27me3_dien_dpt_lsi_regress_srt.rds')
rna_dien %>% write_rds('data/trajectories/dien/RNA_dien_dpt_srt.rds')










