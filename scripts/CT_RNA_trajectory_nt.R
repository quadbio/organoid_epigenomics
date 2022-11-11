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


#### Prep data: non-telencephalon trajectory ####
marks <- map(marks, function(x){
    x$nt_traject <- (
        (x$lineage == 'nt') &
            (x$stage != 'retina') &
            (x$age != '8mo')
    )
    return(x)
})

rna$nt_traject <- (
    (rna$lineage == 'nt') &
        (rna$stage != 'retina') &
        (rna$age != '8mo')
)

dim_plot(rna, group.by='nt_traject', reduction='cssumap')



#### H3K27ac ####
H3K27ac_nt <- marks$H3K27ac %>% subset(nt_traject==TRUE)

H3K27ac_nt@active.assay <- 'peaks'
H3K27ac_nt <- H3K27ac_nt %>%
    RunTFIDF() %>% 
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

lsi <- CreateAssayObject(t(H3K27ac_nt[['lsi']]@cell.embeddings))
H3K27ac_nt[['lsi']] <- NULL
H3K27ac_nt[['lsi']] <- lsi
H3K27ac_nt <- H3K27ac_nt %>% ScaleData(assay='lsi', vars.to.regress=c('log_peak_counts', 'nCount_peaks', 'age'))
lsi <- t(H3K27ac_nt[['lsi']]@scale.data)
H3K27ac_nt[['lsi']] <- NULL
H3K27ac_nt[['lsi']] <- CreateDimReducObject(lsi, key='LSI_', assay='peaks')

H3K27ac_nt <- H3K27ac_nt %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K27ac_nt)
ElbowPlot(H3K27ac_nt, reduction='lsi')

dim_plot(H3K27ac_nt, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27ac_nt@active.assay <- 'cRNA'
feature_plot(H3K27ac_nt, label=T, features=c('STMN2', 'VIM', 'NES', 'GRIA2', 'log_peak_counts'), order=T, pt.size=2)

H3K27ac_pca <- H3K27ac_nt[['lsi']]@cell.embeddings[,2:10]
H3K27ac_diffmap <- DiffusionMap(H3K27ac_pca, verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5]))

H3K27ac_nt[['diffmap']] <- CreateDimReducObject(
    H3K27ac_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27ac_nt <- AddMetaData(H3K27ac_nt, H3K27ac_dc_df)

H3K27ac_nt %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'NES', 'DC1', 'DC2', 'DC3', 'DC4', 'log_peak_counts'), order=T)

mexp <- as.numeric(H3K27ac_nt[['cRNA']]['STMN2'])
plot(H3K27ac_nt$DC1, mexp, col=factor(H3K27ac_nt$age))

mexp <- as.numeric(H3K27ac_nt[['cRNA']]['VIM'])
plot(H3K27ac_nt$DC1, mexp, col=factor(H3K27ac_nt$age))

plot(H3K27ac_nt$DC1, H3K27ac_nt$log_peak_counts, col=factor(H3K27ac_nt$age))

H3K27ac_nt$nt_pt <- rank(-H3K27ac_nt$DC1) / max(rank(H3K27ac_nt$DC1))
H3K27ac_nt %>% write_rds('data/trajectories/nt/H3K27ac_nt_dpt_lsi_regress_srt.rds')



#### H3K4me3 ####
H3K4me3_nt <- marks$H3K4me3 %>% subset(nt_traject==TRUE)

H3K4me3_nt@active.assay <- 'peaks'
H3K4me3_nt <- H3K4me3_nt %>%
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

lsi <- CreateAssayObject(t(H3K4me3_nt[['lsi']]@cell.embeddings))
H3K4me3_nt[['lsi']] <- NULL
H3K4me3_nt[['lsi']] <- lsi
H3K4me3_nt <- H3K4me3_nt %>% ScaleData(assay='lsi', vars.to.regress=c('log_peak_counts', 'nCount_peaks', 'age'))
lsi <- t(H3K4me3_nt[['lsi']]@scale.data)
H3K4me3_nt[['lsi']] <- NULL
H3K4me3_nt[['lsi']] <- CreateDimReducObject(lsi, key='LSI_', assay='peaks')

H3K4me3_nt <- H3K4me3_nt %>% 
    RunUMAP(dims=c(2:10), reduction='lsi')

DepthCor(H3K4me3_nt)
ElbowPlot(H3K4me3_nt, reduction='lsi')

dim_plot(H3K4me3_nt, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K4me3_nt@active.assay <- 'cRNA'
feature_plot(H3K4me3_nt, features=c('STMN2', 'VIM', 'NES', 'GRIA2', 'RSPO3', 'log_peak_counts'), order=T, pt.size=2)

H3K4me3_pca <- H3K4me3_nt[['lsi']]@cell.embeddings[,2:10]
H3K4me3_diffmap <- DiffusionMap(H3K4me3_pca, verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5]))

H3K4me3_nt[['diffmap']] <- CreateDimReducObject(
    H3K4me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K4me3_nt <- AddMetaData(H3K4me3_nt, H3K4me3_dc_df)

H3K4me3_nt %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'log_peak_counts'), order=T)

mexp <- as.numeric(H3K4me3_nt[['cRNA']]['STMN2'])
plot(H3K4me3_nt$DC1, mexp, col=factor(H3K4me3_nt$age))

mexp <- as.numeric(H3K4me3_nt[['cRNA']]['VIM'])
plot(H3K4me3_nt$DC1, mexp, col=factor(H3K4me3_nt$age))

plot(H3K4me3_nt$DC1, H3K4me3_nt$log_peak_counts, col=factor(H3K4me3_nt$age))

H3K4me3_nt$nt_pt <- rank(-H3K4me3_nt$DC1) / max(rank(H3K4me3_nt$DC1))
H3K4me3_nt %>% write_rds('data/trajectories/nt/H3K4me3_nt_dpt_lsi_regress_srt.rds')


#### H3K27me3 ####
H3K27me3_nt <- marks$H3K27me3 %>% subset(nt_traject==TRUE)

H3K27me3_nt@active.assay <- 'peaks'
H3K27me3_nt <- H3K27me3_nt %>%
    RunTFIDF(assay='peaks') %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD(assay='peaks')

lsi <- t(H3K27me3_nt[['lsi']]@cell.embeddings)
H3K27me3_nt[['lsi']] <- NULL
H3K27me3_nt[['lsi']] <- CreateAssayObject(lsi)
H3K27me3_nt <- H3K27me3_nt %>% ScaleData(assay='lsi', vars.to.regress=c('log_peak_counts', 'nCount_peaks', 'age'))
lsi <- H3K27me3_nt[['lsi']]@scale.data
H3K27me3_nt[['lsi']] <- NULL
H3K27me3_nt[['lsi']] <- CreateDimReducObject(t(lsi), key='LSI_', assay='peaks')

H3K27me3_nt <- H3K27me3_nt %>% 
    RunUMAP(dims=2:15, reduction='lsi')

DepthCor(H3K27me3_nt, reduction='lsi')
ElbowPlot(H3K27me3_nt, reduction='lsi')

dim_plot(H3K27me3_nt, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27me3_nt@active.assay <- 'cRNA'
feature_plot(H3K27me3_nt, features=c('STMN2', 'TCF7L2', 'NEUROD2', 'SOX2', 'WLS', 'log_peak_counts'), order=T, pt.size=2)

H3K27me3_pca <- H3K27me3_nt[['lsi']]@cell.embeddings[,2:15]
H3K27me3_diffmap <- DiffusionMap(H3K27me3_pca, verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5]))

H3K27me3_nt[['diffmap']] <- CreateDimReducObject(
    H3K27me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27me3_nt <- AddMetaData(H3K27me3_nt, H3K27me3_dc_df)

H3K27me3_nt %>% 
    feature_plot(features=c('STMN2', 'VIM', 'RSPO3', 'SOX2', 'DC1', 'DC2', 'DC3', 'DC4', 'log_peak_counts'), order=T)

mexp <- as.numeric(H3K27me3_nt[['cRNA']]['WLS'])
plot(H3K27me3_nt$DC1, mexp, col=factor(H3K27me3_nt$age))

mexp <- as.numeric(H3K27me3_nt[['cRNA']]['GRIA2'])
plot(H3K27me3_nt$DC1, mexp, col=factor(H3K27me3_nt$age))

plot(H3K27me3_nt$DC1, H3K27me3_nt$log_peak_counts, col=factor(H3K27me3_nt$age))

H3K27me3_nt$nt_pt <- rank(H3K27me3_nt$DC1) / max(rank(H3K27me3_nt$DC1))
H3K27me3_nt %>% write_rds('data/trajectories/nt/H3K27me3_nt_dpt_lsi_regress_srt.rds')


#### Re-compute PCA for nt trajectory ####
rna_nt <- rna %>% subset(nt_traject==TRUE)

rna_nt <- rna_nt %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

ElbowPlot(rna_nt)

rna_nt <- rna_nt %>% 
    RunUMAP(reduction='css', dims=1:ncol(rna_nt[['css']]))

rna_nt %>% dim_plot(label=T, group.by=c('clusters', 'celltype_jf', 'age'))
rna_nt %>% feature_plot(features=c('STMN2', 'VIM', 'RSPO3', 'SOX2', 'pseudotime_ranks'), order=T)

rna %>% feature_plot(features=c('STMN2', 'VIM', 'RSPO3', 'SOX2', 'pseudotime_ranks'), order=T, reduction='cssumap')

rna_nt$nt_pt <- rank(rna_nt$pseudotime_ranks) / max(rank(rna_nt$pseudotime_ranks))
rna_nt %>% write_rds('data/trajectories/nt/RNA_nt_dpt_srt.rds')


nt_de <- de(rna_nt, groups = 'state')

nt_de_ct <- nt_de %>% 
    filter(feature%in%rownames(H3K27ac_nt), group=='neuron')
vulcano_plot(nt_de_ct, top=30)


#### Subsample modalities to minimal # cells for each timepoint ####
H3K27ac_nt <- read_rds('data/trajectories/nt/H3K27ac_nt_dpt_lsi_regress_srt.rds')
H3K4me3_nt <- read_rds('data/trajectories/nt/H3K4me3_nt_dpt_lsi_regress_srt.rds')
H3K27me3_nt <- read_rds('data/trajectories/nt/H3K27me3_nt_dpt_lsi_regress_srt.rds')
rna_nt <- read_rds('data/trajectories/nt/RNA_nt_dpt_srt.rds')

H3K27ac_nt$age <- ifelse(H3K27ac_nt$age=='128d', '4mo', H3K27ac_nt$age)
H3K27me3_nt$age <- ifelse(H3K27me3_nt$age=='128d', '4mo', H3K27me3_nt$age)
H3K4me3_nt$age <- ifelse(H3K4me3_nt$age=='128d', '4mo', H3K4me3_nt$age)

ncells_table <- bind_rows(
    table(H3K27ac_nt$age),
    table(H3K27me3_nt$age),
    table(H3K4me3_nt$age),
    table(rna_nt$age)
)

ncells_min <- colMins(as.matrix(ncells_table))
ncells_min[ncells_min<100] <- 100
names(ncells_min) <- colnames(ncells_table)

subs_list <- map(c('H3K27ac'=H3K27ac_nt, 'H3K4me3'=H3K4me3_nt, 'H3K27me3'=H3K27me3_nt, 'RNA'=rna_nt), function(srt){
    meta <- srt@meta.data %>% as_tibble(rownames='cell') %>% 
        group_by(age) %>% group_split() %>% 
        map(~sample_n(.x, size=min(ncells_min[.x$age[1]], nrow(.x)))) %>% 
        bind_rows()
    return(subset(srt, cells=meta$cell))
})

H3K27ac_nt_subs <- subs_list$H3K27ac
H3K27me3_nt_subs <- subs_list$H3K27me3
H3K4me3_nt_subs <- subs_list$H3K4me3
rna_nt_subs <- subs_list$RNA

# Re-rank pseudotime
H3K27ac_nt_subs$nt_pt <- rank(H3K27ac_nt_subs$nt_pt) / max(rank(H3K27ac_nt_subs$nt_pt))
H3K27me3_nt_subs$nt_pt <- rank(H3K27me3_nt_subs$nt_pt) / max(rank(H3K27me3_nt_subs$nt_pt))
H3K4me3_nt_subs$nt_pt <- rank(H3K4me3_nt_subs$nt_pt) / max(rank(H3K4me3_nt_subs$nt_pt))
rna_nt_subs$nt_pt <- rank(rna_nt_subs$nt_pt) / max(rank(rna_nt_subs$nt_pt))


#### Make pseudotime bins ####
H3K27ac_nt_subs$pt_bins <- as.numeric(cut(H3K27ac_nt_subs$nt_pt, 20, labels=1:20))
H3K27me3_nt_subs$pt_bins <- as.numeric(cut(H3K27me3_nt_subs$nt_pt, 20, labels=1:20))
H3K4me3_nt_subs$pt_bins <- as.numeric(cut(H3K4me3_nt_subs$nt_pt, 20, labels=1:20))
rna_nt_subs$pt_bins <- as.numeric(cut(rna_nt_subs$nt_pt, 20, labels=1:20))

H3K27ac_nt_subs <- Pando::aggregate_assay(H3K27ac_nt_subs, assay='cRNA', group_name='pt_bins')
H3K27me3_nt_subs <- Pando::aggregate_assay(H3K27me3_nt_subs, assay='cRNA', group_name='pt_bins')
H3K4me3_nt_subs <- Pando::aggregate_assay(H3K4me3_nt_subs, assay='cRNA', group_name='pt_bins')
rna_nt_subs <- Pando::aggregate_assay(rna_nt_subs, assay='RNA', group_name='pt_bins')





#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_clusters <- H3K27me3_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_clusters <- H3K4me3_nt_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
rna_clusters <- rna_nt_subs@assays$RNA@misc$summary$pt_bins[as.character(1:20), ]

genes_plot <- c('OTX2', 'LHX5', 'WLS', 'SOX2', 'VIM', 'STMN2', 'GRIA2', 'DCX')

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

ggsave('plots/trajectories/nt_subs_gene_expr_line.png', width=10, height=6)


H3K27ac_nt_subs %>% write_rds('data/trajectories/nt/H3K27ac_nt_subs_dpt_lsi_regress_srt.rds')
H3K4me3_nt_subs %>% write_rds('data/trajectories/nt/H3K4me3_nt_subs_dpt_lsi_regress_srt.rds')
H3K27me3_nt_subs %>% write_rds('data/trajectories/nt/H3K27me3_nt_subs_dpt_lsi_regress_srt.rds')
rna_nt_subs %>% write_rds('data/trajectories/nt/RNA_nt_subs_dpt_srt.rds')



#### Annotate Neurons ####
rna_nt$pt_bins <- as.numeric(cut(rna_nt$nt_pt, 20, labels=1:20))
meta <- rna_nt@meta.data %>% 
    as_tibble(rownames='cell')
    
ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K27ac_nt$pt_bins <- as.numeric(cut(H3K27ac_nt$nt_pt, 20, labels=1:20))
meta <- H3K27ac_nt@meta.data %>% 
    as_tibble(rownames='cell')
    
ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K27me3_nt$pt_bins <- as.numeric(cut(H3K27me3_nt$nt_pt, 20, labels=1:20))
meta <- H3K27me3_nt@meta.data %>% 
    as_tibble(rownames='cell')
    
ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()

H3K4me3_nt$pt_bins <- as.numeric(cut(H3K4me3_nt$nt_pt, 20, labels=1:20))
meta <- H3K4me3_nt@meta.data %>% 
    as_tibble(rownames='cell')
    
ggplot(meta, aes(pt_bins, fill=celltype_jf)) +
    geom_bar()


rna_nt$neuron <- rna_nt$pt_bins > 15
H3K27ac_nt$neuron <- H3K27ac_nt$pt_bins > 15
H3K27me3_nt$neuron <- H3K27me3_nt$pt_bins > 15
H3K4me3_nt$neuron <- H3K4me3_nt$pt_bins > 15

H3K27ac_nt %>% write_rds('data/trajectories/nt/H3K27ac_nt_dpt_lsi_regress_srt.rds')
H3K4me3_nt %>% write_rds('data/trajectories/nt/H3K4me3_nt_dpt_lsi_regress_srt.rds')
H3K27me3_nt %>% write_rds('data/trajectories/nt/H3K27me3_nt_dpt_lsi_regress_srt.rds')
rna_nt %>% write_rds('data/trajectories/nt/RNA_nt_dpt_srt.rds')



















