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


cluster_srt@active.assay <- 'H3K27me3_RNA'
feature_plot(cluster_srt, label=T, features=c('STMN2', 'VIM', 'PAX6', 'VSX2', 'LHX5', 'RSPO3'), order=T, pt.size=2)


#### Prep data: RPC -> RGC ####
marks <- map(marks, function(x){
    x$ret_traject <- (
        x$celltype_jf %in% c('RPC', 'RGC') &
            x$stage == 'retina'
    )
    return(x)
})

rna$ret_traject <- (
    rna$celltype_jf %in% c('RPC', 'RGC') &
        rna$stage == 'retina'
)

#### Get pseudotimes with diffmaps ####
#### H3K27ac ####
H3K27ac_ret <- marks$H3K27ac %>% subset(ret_traject==TRUE)

clusters_use <- H3K27ac_ret$clusters %>% setdiff(
    c('retina_25', 'retina_9'))
H3K27ac_ret$ret_traject <- H3K27ac_ret$clusters %in% clusters_use
H3K27ac_ret <- subset(H3K27ac_ret, clusters%in%clusters_use)

H3K27ac_ret <- H3K27ac_ret %>%
    RunTFIDF() %>% 
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

# H3K27ac_ret[['lsi']] <- H3K27ac_ret[['clsi']]

H3K27ac_ret <- H3K27ac_ret %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K27ac_ret)
ElbowPlot(H3K27ac_ret, reduction='lsi')

dim_plot(H3K27ac_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27ac_ret@active.assay <- 'cRNA'
feature_plot(H3K27ac_ret, label=T, features=c('STMN2', 'PTH2', 'PAX6', 'VSX2', 'GRIA2', 'VIM'), order=T, pt.size=2)

H3K27ac_pca <- H3K27ac_ret[['lsi']]@cell.embeddings[,2:10]
H3K27ac_diffmap <- DiffusionMap(H3K27ac_pca, verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5]))

H3K27ac_ret[['diffmap']] <- CreateDimReducObject(
    H3K27ac_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27ac_ret <- AddMetaData(H3K27ac_ret, H3K27ac_dc_df)

H3K27ac_ret %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5'), order=T)

H3K27ac_ret %>% 
    dim_plot(group.by=c('age', 'celltype_jf', 'clusters'), order=T)

mexp <- as.numeric(H3K27ac_ret[['cRNA']]['STMN2'])
mexp <- as.numeric(H3K27ac_ret[['cRNA']]['EBF3'])
plot(H3K27ac_ret$DC1, mexp, col=factor(H3K27ac_ret$age))
plot(H3K27ac_ret$DC1, H3K27ac_ret$log_peak_counts, col=factor(H3K27ac_ret$age))

cor(mexp, H3K27ac_dc_df)

H3K27ac_ret %>% write_rds('data/trajectories/retina/H3K27ac_retina_dpt_lsi_srt.rds')



#### H3K27ac 6w ####
H3K27ac_ret <- marks$H3K27ac %>% subset(ret_traject==TRUE & age=='ret-6w')

H3K27ac_ret <- H3K27ac_ret %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

# H3K27ac_ret[['lsi']] <- H3K27ac_ret[['clsi']]

clusters_use <- H3K27ac_ret$clusters %>% setdiff(
    c('retina_25', 'retina_9'))
H3K27ac_ret$ret_traject <- H3K27ac_ret$clusters %in% clusters_use
H3K27ac_ret <- subset(H3K27ac_ret, clusters%in%clusters_use)

DepthCor(H3K27ac_ret)
ElbowPlot(H3K27ac_ret, reduction='lsi')

H3K27ac_ret <- H3K27ac_ret %>% 
    RunUMAP(dims=2:10, reduction='lsi')

dim_plot(H3K27ac_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27ac_ret@active.assay <- 'cRNA'
feature_plot(H3K27ac_ret, label=T, features=c('STMN2', 'VIM', 'PAX6', 'VSX2', 'LHX5', 'RSPO3'), order=T, pt.size=2)

H3K27ac_pca <- H3K27ac_ret[['lsi']]@cell.embeddings[,2:10]
H3K27ac_diffmap <- DiffusionMap(H3K27ac_pca, verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5]))

H3K27ac_ret[['diffmap']] <- CreateDimReducObject(
    H3K27ac_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27ac_ret <- AddMetaData(H3K27ac_ret, H3K27ac_dc_df)

H3K27ac_ret %>% 
    feature_plot(features=c('STMN2', 'VIM', 'DC1', 'DC2', 'DC3'), order=T)

H3K27ac_ret %>% 
    dim_plot(group.by=c('age', 'celltype_jf', 'clusters'), order=T)

H3K27ac_ret %>% write_rds('data/trajectories/retina/H3K27ac_retina_6w_dpt_srt.rds')


#### H3K27ac 12w ####
H3K27ac_ret <- marks$H3K27ac %>% subset(ret_traject==TRUE & age=='ret-12w')

clusters_use <- H3K27ac_ret$clusters %>% setdiff(
    c('retina_25', 'retina_9', 'retuna_29', 'retina_27', 'retina_39', 'retina_33'))
H3K27ac_ret$ret_traject <- H3K27ac_ret$clusters %in% clusters_use
H3K27ac_ret <- subset(H3K27ac_ret, clusters%in%clusters_use)

H3K27ac_ret <- H3K27ac_ret %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

# H3K27ac_ret[['lsi']] <- H3K27ac_ret[['clsi']]

DepthCor(H3K27ac_ret)
ElbowPlot(H3K27ac_ret, reduction='lsi')

H3K27ac_ret <- H3K27ac_ret %>% 
    RunUMAP(dims=2:10, reduction='lsi')

dim_plot(H3K27ac_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27ac_ret@active.assay <- 'cRNA'
feature_plot(H3K27ac_ret, label=T, features=c('STMN2', 'VIM', 'PAX6', 'VSX2', 'LHX5', 'RSPO3'), order=T, pt.size=2)

H3K27ac_pca <- H3K27ac_ret[['lsi']]@cell.embeddings[,2:10]
H3K27ac_diffmap <- DiffusionMap(H3K27ac_pca, verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5]))

H3K27ac_ret[['diffmap']] <- CreateDimReducObject(
    H3K27ac_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27ac_ret <- AddMetaData(H3K27ac_ret, H3K27ac_dc_df)

H3K27ac_ret %>% 
    feature_plot(features=c('STMN2', 'VIM', 'DC1', 'DC2', 'DC3', 'log_peak_counts'), order=T)

H3K27ac_ret %>% 
    dim_plot(group.by=c('age', 'celltype_jf', 'clusters'), order=T)

H3K27ac_ret %>% write_rds('data/trajectories/retina/H3K27ac_retina_12w_dpt_srt.rds')





#### H3K4me3 ####
H3K4me3_ret <- marks$H3K4me3 %>% subset(ret_traject==TRUE)

H3K4me3_ret@active.assay <- 'peaks'
H3K4me3_ret <- H3K4me3_ret %>%
    RunTFIDF() %>% 
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

# lsi <- CreateAssayObject(t(H3K4me3_ret[['clsi']]@cell.embeddings))
lsi <- CreateAssayObject(t(H3K4me3_ret[['lsi']]@cell.embeddings))
H3K4me3_ret[['lsi']] <- NULL
H3K4me3_ret[['lsi']] <- lsi
H3K4me3_ret <- H3K4me3_ret %>% ScaleData(assay='lsi', vars.to.regress=c('nCount_peaks', 'log_peak_counts'))
lsi <- t(H3K4me3_ret[['lsi']]@scale.data)
H3K4me3_ret[['lsi']] <- NULL
H3K4me3_ret[['lsi']] <- CreateDimReducObject(lsi, key='LSI_', assay='peaks')

H3K4me3_ret <- H3K4me3_ret %>% 
    RunUMAP(dims=c(2:10), reduction='lsi')

DepthCor(H3K4me3_ret)
ElbowPlot(H3K4me3_ret, reduction='lsi')

dim_plot(H3K4me3_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K4me3_ret@active.assay <- 'cRNA'
feature_plot(H3K4me3_ret, features=c('STMN2', 'PTH2', 'PAX6', 'VSX2', 'GRIA2', 'VIM', 'log_peak_counts'), order=T, pt.size=2)

H3K4me3_pca <- H3K4me3_ret[['lsi']]@cell.embeddings[,2:6]
H3K4me3_diffmap <- DiffusionMap(H3K4me3_pca, verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5]))

H3K4me3_ret[['diffmap']] <- CreateDimReducObject(
    H3K4me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K4me3_ret <- AddMetaData(H3K4me3_ret, H3K4me3_dc_df)

H3K4me3_ret %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5'), order=T)

mexp <- as.numeric(H3K4me3_ret[['cRNA']]['STMN2'])
mexp <- as.numeric(H3K4me3_ret[['cRNA']]['VIM'])
plot(H3K4me3_ret$DC1, as.numeric(H3K4me3_ret[['cRNA']]['STMN2']), col=factor(H3K4me3_ret$age))
plot(H3K4me3_ret$DC1, as.numeric(H3K4me3_ret[['cRNA']]['VSX2']), col=factor(H3K4me3_ret$age))

plot(H3K4me3_ret$DC1, H3K4me3_ret$log_peak_counts, col=factor(H3K4me3_ret$age))

cor(mexp, H3K4me3_dc_df, method = 'spearman')

H3K4me3_ret %>% write_rds('data/trajectories/retina/H3K4me3_retina_dpt_regress_srt.rds')



#### H3K4me3 6w ####
H3K4me3_ret <- marks$H3K4me3 %>% subset(ret_traject==TRUE & age=='ret-6w')

H3K4me3_ret <- H3K4me3_ret %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

# H3K4me3_ret[['lsi']] <- H3K4me3_ret[['clsi']]

H3K4me3_ret <- H3K4me3_ret %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K4me3_ret)
ElbowPlot(H3K4me3_ret, reduction='lsi')

dim_plot(H3K4me3_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K4me3_ret@active.assay <- 'cRNA'
feature_plot(H3K4me3_ret, label=T, features=c('STMN2', 'PTH2', 'PAX6', 'VSX2', 'GRIA2', 'VIM'), order=T, pt.size=2)

H3K4me3_pca <- H3K4me3_ret[['lsi']]@cell.embeddings[,2:10]
H3K4me3_diffmap <- DiffusionMap(H3K4me3_pca, verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5]))

H3K4me3_ret[['diffmap']] <- CreateDimReducObject(
    H3K4me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K4me3_ret <- AddMetaData(H3K4me3_ret, H3K4me3_dc_df)

H3K4me3_ret %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5'), order=T)

H3K4me3_ret %>% 
    dim_plot(group.by=c('age', 'celltype_jf', 'clusters'), order=T)

H3K4me3_ret %>% write_rds('data/trajectories/retina/H3K4me3_retina_6w_dpt_srt.rds')



#### H3K4me3 12w ####
H3K4me3_ret <- marks$H3K4me3 %>% subset(ret_traject==TRUE & age=='ret-12w')

H3K4me3_ret <- H3K4me3_ret %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

# H3K4me3_ret[['lsi']] <- H3K4me3_ret[['clsi']]

H3K4me3_ret <- H3K4me3_ret %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K4me3_ret)
ElbowPlot(H3K4me3_ret, reduction='lsi')

dim_plot(H3K4me3_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K4me3_ret@active.assay <- 'cRNA'
feature_plot(H3K4me3_ret, label=T, features=c('STMN2', 'PTH2', 'PAX6', 'VSX2', 'GRIA2', 'VIM'), order=T, pt.size=2)

H3K4me3_pca <- H3K4me3_ret[['lsi']]@cell.embeddings[,2:10]
H3K4me3_diffmap <- DiffusionMap(H3K4me3_pca, verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5]))

H3K4me3_ret[['diffmap']] <- CreateDimReducObject(
    H3K4me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K4me3_ret <- AddMetaData(H3K4me3_ret, H3K4me3_dc_df)

H3K4me3_ret %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5'), order=T)

H3K4me3_ret %>% 
    dim_plot(group.by=c('age', 'celltype_jf', 'clusters'), order=T)

H3K4me3_ret %>% write_rds('data/trajectories/retina/H3K4me3_retina_12w_dpt_srt.rds')




#### H3K27me3 ####
H3K27me3_ret <- marks$H3K27me3 %>% subset(ret_traject==TRUE)

# H3K27me3_ret <- H3K27me3_ret %>%
#     RunTFIDF() %>%
#     FindTopFeatures(min.cutoff='q80') %>%
#     RunSVD()

H3K27me3_ret[['lsi']] <- H3K27me3_ret[['clsi']]

lsi <- CreateAssayObject(t(H3K27me3_ret[['lsi']]@cell.embeddings))
H3K27me3_ret[['lsi']] <- NULL
H3K27me3_ret[['lsi']] <- lsi
H3K27me3_ret <- H3K27me3_ret %>% ScaleData(assay='lsi', vars.to.regress=c('nCount_peaks', 'log_peak_counts'))
lsi <- t(H3K27me3_ret[['lsi']]@scale.data)
H3K27me3_ret[['lsi']] <- NULL
H3K27me3_ret[['lsi']] <- CreateDimReducObject(lsi, key='LSI_', assay='peaks')

H3K27me3_ret <- H3K27me3_ret %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K27me3_ret, reduction='lsi')
ElbowPlot(H3K27me3_ret, reduction='lsi')

dim_plot(H3K27me3_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27me3_ret@active.assay <- 'cRNA'
feature_plot(H3K27me3_ret, features=c('STMN2', 'PTH2', 'PAX6', 'VSX2', 'GRIA2', 'VIM'), order=T, pt.size=2)

H3K27me3_pca <- H3K27me3_ret[['lsi']]@cell.embeddings[,2:10]
H3K27me3_diffmap <- DiffusionMap(H3K27me3_pca, verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5]))

H3K27me3_ret[['diffmap']] <- CreateDimReducObject(
    H3K27me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27me3_ret <- AddMetaData(H3K27me3_ret, H3K27me3_dc_df)

H3K27me3_ret %>% 
    feature_plot(features=c('STMN2', 'VSX2', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5'), order=T)

mexp <- as.numeric(H3K27me3_ret[['cRNA']]['STMN2'])
mexp <- as.numeric(H3K27me3_ret[['cRNA']]['VSX2'])
plot(H3K27me3_ret$DC3, mexp, col=factor(H3K27me3_ret$age))
plot(H3K27me3_ret$DC3, H3K27me3_ret$log_peak_counts, col=factor(H3K27me3_ret$age))

cor(mexp, H3K27me3_dc_df)

H3K27me3_ret %>% write_rds('data/trajectories/retina/H3K27me3_retina_dpt_lsi_srt.rds')


#### H3K27me3 6w ####
H3K27me3_ret <- marks$H3K27me3 %>% subset(ret_traject==TRUE & age=='ret-6w')

H3K27me3_ret <- H3K27me3_ret %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

# H3K27me3_ret[['lsi']] <- H3K27me3_ret[['clsi']]

H3K27me3_ret <- H3K27me3_ret %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K27me3_ret)
ElbowPlot(H3K27me3_ret, reduction='lsi')

dim_plot(H3K27me3_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27me3_ret@active.assay <- 'cRNA'
feature_plot(H3K27me3_ret, label=T, features=c('STMN2', 'PTH2', 'PAX6', 'VSX2', 'GRIA2', 'VIM'), order=T, pt.size=2)

H3K27me3_pca <- H3K27me3_ret[['lsi']]@cell.embeddings[,2:10]
H3K27me3_diffmap <- DiffusionMap(H3K27me3_pca, verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5]))

H3K27me3_ret[['diffmap']] <- CreateDimReducObject(
    H3K27me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27me3_ret <- AddMetaData(H3K27me3_ret, H3K27me3_dc_df)

H3K27me3_ret %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5'), order=T)

H3K27me3_ret %>% 
    dim_plot(group.by=c('age', 'celltype_jf', 'clusters'), order=T)

H3K27me3_ret %>% write_rds('data/trajectories/retina/H3K27me3_retina_6w_dpt_srt.rds')



#### H3K27me3 12w ####
H3K27me3_ret <- marks$H3K27me3 %>% subset(ret_traject==TRUE & age=='ret-12w')

H3K27me3_ret <- H3K27me3_ret %>%
    FindTopFeatures(min.cutoff='q80') %>%
    RunSVD()

# H3K27me3_ret[['lsi']] <- H3K27me3_ret[['clsi']]

H3K27me3_ret <- H3K27me3_ret %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K27me3_ret)
ElbowPlot(H3K27me3_ret, reduction='lsi')

dim_plot(H3K27me3_ret, label=T, group.by=c('clusters', 'celltype_jf', 'age'))

H3K27me3_ret@active.assay <- 'cRNA'
feature_plot(H3K27me3_ret, label=T, features=c('STMN2', 'PTH2', 'PAX6', 'VSX2', 'GRIA2', 'VIM'), order=T, pt.size=2)

H3K27me3_pca <- H3K27me3_ret[['lsi']]@cell.embeddings[,2:10]
H3K27me3_diffmap <- DiffusionMap(H3K27me3_pca, verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5]))

H3K27me3_ret[['diffmap']] <- CreateDimReducObject(
    H3K27me3_diffmap@eigenvectors*100, key='DC_', assay='peaks'
)
H3K27me3_ret <- AddMetaData(H3K27me3_ret, H3K27me3_dc_df)

H3K27me3_ret %>% 
    feature_plot(features=c('STMN2', 'VIM', 'GRIA2', 'EBF3', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5'), order=T)

H3K27me3_ret %>% 
    dim_plot(group.by=c('age', 'celltype_jf', 'clusters'), order=T)

H3K27me3_ret %>% write_rds('data/trajectories/retina/H3K27me3_retina_12w_dpt_srt.rds')



#### Re-compute PCA for ret trajectory ####
rna_ret <- rna %>% subset(ret_traject==TRUE)

rna_ret <- rna_ret %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

ElbowPlot(rna_ret)

rna_ret <- rna_ret %>% 
    RunUMAP(reduction='css', dims=1:ncol(rna_ret[['css']]))

#### Run diffmap for RNA cortex ####
rna_ret %>% dim_plot(group.by=c('clusters', 'celltype_jf', 'age'), label=T, order=T, reduction='umap')
rna_ret %>% feature_plot(features=c('STMN2', 'VSX2', 'VIM', 'pseudotime_ranks'), order=T, reduction='umap')

# rna_pca <- rna_ret[['css']]@cell.embeddings
rna_diffmap <- DiffusionMap(rna_pca, verbose=T)
rna_ret[['diffmap']] <- CreateDimReducObject(
    rna_diffmap@eigenvectors*100, key='DC_', assay='RNA'
)
rna_dc_df <- as.data.frame(rank_cols(rna_diffmap@eigenvectors[,1:5]))
rna_ret <- AddMetaData(rna_ret, rna_dc_df)

rna_ret %>%
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'pseudotime_ranks'), reduction='cssumap')

rna_ret %>%
    feature_plot(features=c('STMN2', 'SOX2', 'DC1', 'DC2', 'DC3', 'pseudotime_ranks'), reduction='diffmap')

rna_ret %>% write_rds('data/trajectories/retina/RNA_retina_dpt_srt.rds')







#### Align pseudotimes ####
H3K27ac_ret <- read_rds('data/trajectories/retina/H3K27ac_retina_dpt_srt.rds')
H3K4me3_ret <- read_rds('data/trajectories/retina/H3K4me3_retina_dpt_srt.rds')
H3K27me3_ret <- read_rds('data/trajectories/retina/H3K27me3_retina_dpt_srt.rds')

rna_ret <- read_rds('data/trajectories/retina/RNA_retina_dpt_srt.rds')

H3K27ac_ret$ret_pt <- rank(-H3K27ac_ret$DC3) / max(rank(-H3K27ac_ret$DC3))
H3K4me3_ret$ret_pt <- rank(H3K4me3_ret$DC3) / max(rank(H3K4me3_ret$DC3))
H3K27me3_ret$ret_pt <- rank(H3K27me3_ret$DC3) / max(rank(H3K27me3_ret$DC3))
rna_ret$ret_pt <- rank(rna_ret$velocity_pseudotime) / max(rank(rna_ret$velocity_pseudotime))

plot(H3K27me3_ret$DC1, H3K27me3_ret$log_peak_counts)

p1 <- H3K27ac_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ret_pt'), order=T)
p2 <- H3K4me3_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ret_pt'), order=T)
p3 <- H3K27me3_6w_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ret_pt'), order=T)
p4 <- rna_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'pseudotime_ranks', 'ret_pt'), order=T, reduction='umap')

p1 / p2 / p3 / p4
ggsave('plots/trajectories/retina_diff_pt_umap.png', width=8, height=20)

#### Check out expression patterns for different axes ####
#### Make pseudotime bins ####

H3K27ac_ret$pt_bins <- as.numeric(cut(H3K27ac_ret$ret_pt, 20, labels=1:20))
H3K27me3_ret$pt_bins <- as.numeric(cut(H3K27me3_ret$ret_pt, 20, labels=1:20))
H3K4me3_ret$pt_bins <- as.numeric(cut(H3K4me3_ret$ret_pt, 20, labels=1:20))
rna_ret$pt_bins <- as.numeric(cut(rna_ret$ret_pt, 20, labels=1:20))

H3K27ac_ret <- Pando::aggregate_assay(H3K27ac_ret, assay='cRNA', group_name='pt_bins')
H3K27me3_ret <- Pando::aggregate_assay(H3K27me3_ret, assay='cRNA', group_name='pt_bins')
H3K4me3_ret <- Pando::aggregate_assay(H3K4me3_ret, assay='cRNA', group_name='pt_bins')
rna_ret <- Pando::aggregate_assay(rna_ret, assay='RNA', group_name='pt_bins')

H3K27ac_clusters <- H3K27ac_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_clusters <- H3K27me3_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_clusters <- H3K4me3_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
rna_clusters <- rna_ret@assays$RNA@misc$summary$pt_bins[as.character(1:20), ]

genes_plot <- c('LHX2', 'STMN2', 'GRIA2', 'SFRP2', 'HES1', 'VIM', 'SFRP1', 'SOX4', 'MLLT11', 'CCND1', 'FOS')

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

ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(~gene, scales='free')

H3K27ac_ret %>% write_rds('data/trajectories/retina/H3K27ac_retina_dpt_srt.rds')
H3K4me3_ret %>% write_rds('data/trajectories/retina/H3K4me3_retina_dpt_srt.rds')
H3K27me3_ret %>% write_rds('data/trajectories/retina/H3K27me3_retina_dpt_srt.rds')
rna_ret %>% write_rds('data/trajectories/retina/RNA_retina_dpt_srt.rds')




#### Now we try this for H3K27 only ####
H3K27me3_all_ret <- read_rds('data/trajectories/retina/H3K27me3_retina_dpt_srt.rds')

H3K27me3_6w_ret$ret_pt <- rank(H3K27me3_6w_ret$DC1) / max(rank(H3K27me3_6w_ret$DC1))
p1 <- H3K27me3_6w_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ret_pt'), order=T)

H3K27me3_12w_ret$ret_pt <- rank(-H3K27me3_12w_ret$DC1) / max(rank(-H3K27me3_12w_ret$DC1))
p2 <- H3K27me3_12w_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ret_pt'), order=T)

H3K27me3_all_ret$ret_pt <- rank(H3K27me3_all_ret$DC2) / max(rank(H3K27me3_all_ret$DC2))
p3 <- H3K27me3_all_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ret_pt'), order=T)


H3K27me3_6w_ret$pt_bins <- as.numeric(cut(H3K27me3_6w_ret$ret_pt, 20, labels=1:20))
H3K27me3_12w_ret$pt_bins <- as.numeric(cut(H3K27me3_12w_ret$ret_pt, 20, labels=1:20))
H3K27me3_all_ret$pt_bins <- as.numeric(cut(H3K27me3_all_ret$ret_pt, 20, labels=1:20))

H3K27me3_6w_ret <- Pando::aggregate_assay(H3K27me3_6w_ret, assay='cRNA', group_name='pt_bins')
H3K27me3_12w_ret <- Pando::aggregate_assay(H3K27me3_12w_ret, assay='cRNA', group_name='pt_bins')
H3K27me3_all_ret <- Pando::aggregate_assay(H3K27me3_all_ret, assay='cRNA', group_name='pt_bins')

H3K27me3_6w_clusters <- H3K27me3_6w_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_12w_clusters <- H3K27me3_12w_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K27me3_all_clusters <- H3K27me3_all_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]

w6_expr <- H3K27me3_6w_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

w12_expr <- H3K27me3_12w_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

all_expr <- H3K27me3_all_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

plot_df <- bind_rows('w6'=w6_expr, 'w12'=w12_expr, 'all'=all_expr, .id='age') %>% 
    group_by(age, gene) %>% 
    mutate(expr=scale01(expr))

ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=age, group=age)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    facet_grid(~gene, scales='free')


#### Now we try this for H3K4 only ####
H3K4me3_all_ret <- read_rds('data/trajectories/retina/H3K4me3_retina_dpt_srt.rds')

H3K4me3_6w_ret$ret_pt <- rank(H3K4me3_6w_ret$DC1) / max(rank(H3K4me3_6w_ret$DC1))
p1 <- H3K4me3_6w_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ret_pt'), order=T)

H3K4me3_12w_ret$ret_pt <- rank(H3K4me3_12w_ret$DC1) / max(rank(H3K4me3_12w_ret$DC1))
p2 <- H3K4me3_12w_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ret_pt'), order=T)

H3K4me3_all_ret$ret_pt <- rank(-H3K4me3_all_ret$DC1) / max(rank(-H3K4me3_all_ret$DC1))
p3 <- H3K4me3_all_ret %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'DC3'), order=T)

p1 / p2 / p3


H3K4me3_6w_ret$pt_bins <- as.numeric(cut(H3K4me3_6w_ret$ret_pt, 20, labels=1:20))
H3K4me3_12w_ret$pt_bins <- as.numeric(cut(H3K4me3_12w_ret$ret_pt, 20, labels=1:20))
H3K4me3_all_ret$pt_bins <- as.numeric(cut(H3K4me3_all_ret$ret_pt, 20, labels=1:20))

H3K4me3_6w_ret <- Pando::aggregate_assay(H3K4me3_6w_ret, assay='cRNA', group_name='pt_bins')
H3K4me3_12w_ret <- Pando::aggregate_assay(H3K4me3_12w_ret, assay='cRNA', group_name='pt_bins')
H3K4me3_all_ret <- Pando::aggregate_assay(H3K4me3_all_ret, assay='cRNA', group_name='pt_bins')

H3K4me3_6w_clusters <- H3K4me3_6w_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_12w_clusters <- H3K4me3_12w_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]
H3K4me3_all_clusters <- H3K4me3_all_ret@assays$cRNA@misc$summary$pt_bins[as.character(1:20), ]

genes_plot <- c('LHX2', 'EEF1D', 'HMX1', 'SFRP2', 'HES1', 'VIM', 'SFRP1', 'SOX4', 'MLLT11', 'CCND1', 'FOS')

w6_expr <- H3K4me3_6w_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

w12_expr <- H3K4me3_12w_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

all_expr <- H3K4me3_all_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

plot_df <- bind_rows('w6'=w6_expr, 'w12'=w12_expr, 'all'=all_expr, .id='age') %>% 
    group_by(age, gene) %>% 
    mutate(expr=scale01(expr))

ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=age, group=age)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    facet_grid(~gene, scales='free')




#### Subsample modalities to minimal # cells for each timepoint ####
H3K27ac_ret <- read_rds('data/trajectories/retina/H3K27ac_retina_dpt_srt.rds')
H3K4me3_ret <- read_rds('data/trajectories/retina/H3K4me3_retina_dpt_srt.rds')
H3K27me3_ret <- read_rds('data/trajectories/retina/H3K27me3_retina_dpt_srt.rds')
rna_ret <- read_rds('data/trajectories/retina/RNA_retina_dpt_srt.rds')

H3K27ac_ret$age <- H3K27ac_ret$age %>% str_replace('-', '_')
H3K27me3_ret$age <- H3K27me3_ret$age %>% str_replace('-', '_')
H3K4me3_ret$age <- H3K4me3_ret$age %>% str_replace('-', '_')

ncells_table <- bind_rows(
    table(H3K27ac_ret$age),
    table(H3K27me3_ret$age),
    table(H3K4me3_ret$age),
    table(rna_ret$age)
)

ncells_min <- colMins(as.matrix(ncells_table))
ncells_min[ncells_min<100] <- 100
names(ncells_min) <- colnames(ncells_table)

subs_list <- map(c('H3K27ac'=H3K27ac_ret, 'H3K4me3'=H3K4me3_ret, 'H3K27me3'=H3K27me3_ret, 'RNA'=rna_ret), function(srt){
    meta <- srt@meta.data %>% as_tibble(rownames='cell') %>% 
        group_by(age) %>% group_split() %>% 
        map(~sample_n(.x, size=min(ncells_min[.x$age[1]], nrow(.x)))) %>% 
        bind_rows()
    return(subset(srt, cells=meta$cell))
})

H3K27ac_ret_subs <- subs_list$H3K27ac
H3K27me3_ret_subs <- subs_list$H3K27me3
H3K4me3_ret_subs <- subs_list$H3K4me3
rna_ret_subs <- subs_list$RNA

# Re-rank pseudotime
H3K27ac_ret_subs$ret_pt <- rank(H3K27ac_ret_subs$ret_pt) / max(rank(H3K27ac_ret_subs$ret_pt))
H3K27me3_ret_subs$ret_pt <- rank(H3K27me3_ret_subs$ret_pt) / max(rank(H3K27me3_ret_subs$ret_pt))
H3K4me3_ret_subs$ret_pt <- rank(H3K4me3_ret_subs$ret_pt) / max(rank(H3K4me3_ret_subs$ret_pt))
rna_ret_subs$ret_pt <- rank(rna_ret_subs$ret_pt) / max(rank(rna_ret_subs$ret_pt))


#### Make pseudotime bins ####
H3K27ac_ret_subs$pt_bins <- as.numeric(cut(H3K27ac_ret_subs$ret_pt, 50, labels=1:50))
H3K27me3_ret_subs$pt_bins <- as.numeric(cut(H3K27me3_ret_subs$ret_pt, 50, labels=1:50))
H3K4me3_ret_subs$pt_bins <- as.numeric(cut(H3K4me3_ret_subs$ret_pt, 50, labels=1:50))
rna_ret_subs$pt_bins <- as.numeric(cut(rna_ret_subs$ret_pt, 50, labels=1:50))

H3K27ac_ret_subs <- Pando::aggregate_assay(H3K27ac_ret_subs, assay='cRNA', group_name='pt_bins')
H3K27me3_ret_subs <- Pando::aggregate_assay(H3K27me3_ret_subs, assay='cRNA', group_name='pt_bins')
H3K4me3_ret_subs <- Pando::aggregate_assay(H3K4me3_ret_subs, assay='cRNA', group_name='pt_bins')
rna_ret_subs <- Pando::aggregate_assay(rna_ret_subs, assay='RNA', group_name='pt_bins')

H3K27ac_meta <- H3K27ac_ret_subs@meta.data %>% as_tibble(rownames='cell') 
H3K4me3_meta <- H3K4me3_ret_subs@meta.data %>% as_tibble(rownames='cell') 
H3K27me3_meta <- H3K27me3_ret_subs@meta.data %>% as_tibble(rownames='cell') 
rna_meta <- rna_ret_subs@meta.data %>% as_tibble(rownames='cell') 

p1 <- ggplot(H3K27ac_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27ac')

p2 <- ggplot(H3K4me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K4me3')

p3 <- ggplot(H3K27me3_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('H3K27me3')

p4 <- ggplot(rna_meta, aes(pt_bins, fill=age)) +
    geom_bar(position='fill') + ggtitle('rna')

p1 / p2 / p3 / p4 & theme_void() & scale_fill_manual(values=colours_timescale)
ggsave('plots/trajectories/retina_subs_age_dist_bar.png', width=8, height=5, bg='white')

#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_ret_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_clusters <- H3K27me3_ret_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_clusters <- H3K4me3_ret_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
rna_clusters <- rna_ret_subs@assays$RNA@misc$summary$pt_bins[as.character(1:50), ]

genes_plot <- c('STMN2', 'TFAP2B', 'VSX2', 'EBF3', 'PRDM8', 'SOX2', 'VIM', 'PTH2', 'GRIA2')

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

ggsave('plots/trajectories/retina_subs_gene_expr_line.png', width=10, height=6, bg='white')


H3K27ac_ret_subs %>% write_rds('data/trajectories/retina/H3K27ac_retina_subs_dpt_srt.rds')
H3K4me3_ret_subs %>% write_rds('data/trajectories/retina/H3K4me3_retina_subs_dpt_srt.rds')
H3K27me3_ret_subs %>% write_rds('data/trajectories/retina/H3K27me3_retina_subs_dpt_srt.rds')
rna_ret_subs %>% write_rds('data/trajectories/retina/RNA_retina_subs_dpt_srt.rds')



#### Summarize to highres clusters and order those ####
library(tidygraph)
library(ggraph)

cluster_graph <- read_rds('data/RNA/all_RNA_cluster_graph.rds')
cluster_pt_meta <- cluster_graph %N>% as_tibble() %>% select('clusters'='name', pseudotime_ranks, velocity_pseudotime)

cluster_meta <- cluster_srt@meta.data %>% as_tibble(rownames='clusters') %>% inner_join(cluster_pt_meta)

H3K27ac_ret <- read_rds('data/trajectories/retina/H3K27ac_retina_dpt_srt.rds') %>% 
    NormalizeData(assay='cRNA')
H3K4me3_ret <- read_rds('data/trajectories/retina/H3K4me3_retina_dpt_srt.rds') %>% 
    NormalizeData(assay='cRNA')
H3K27me3_ret <- read_rds('data/trajectories/retina/H3K27me3_retina_dpt_srt.rds') %>% 
    NormalizeData(assay='cRNA')
rna_ret <- read_rds('data/trajectories/retina/RNA_retina_dpt_srt.rds') 

H3K27ac_ret <- Pando::aggregate_assay(H3K27ac_ret, assay='cRNA', group_name='clusters')
H3K27me3_ret <- Pando::aggregate_assay(H3K27me3_ret, assay='cRNA', group_name='clusters')
H3K4me3_ret <- Pando::aggregate_assay(H3K4me3_ret, assay='cRNA', group_name='clusters')
rna_ret <- Pando::aggregate_assay(rna_ret, assay='RNA', group_name='clusters')

H3K27ac_clusters <- H3K27ac_ret@assays$cRNA@misc$summary$clusters
H3K27me3_clusters <- H3K27me3_ret@assays$cRNA@misc$summary$clusters
H3K4me3_clusters <- H3K4me3_ret@assays$cRNA@misc$summary$clusters
rna_clusters <- rna_ret@assays$RNA@misc$summary$clusters

genes_plot <- c('STMN2', 'TFAP2B', 'VSX2', 'EBF3', 'PRDM8', 'SOX2', 'VIM', 'PTH2', 'GRIA2')

H3K27ac_expr <- H3K27ac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='H3K27ac_clusters', values_to='expr') %>% 
    inner_join(cluster_meta)

H3K27me3_expr <- H3K27me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='H3K27me3_clusters', values_to='expr') %>% 
    inner_join(cluster_meta)

H3K4me3_expr <- H3K4me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='H3K4me3_clusters', values_to='expr') %>% 
    inner_join(cluster_meta)

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='clusters', values_to='expr') %>% 
    inner_join(cluster_meta)

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr=scale01(expr))
    

p1 <- ggplot(plot_df, aes(rank(pseudotime_ranks), expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(~gene, scales='free')

p2 <- ggplot(plot_df, aes(rank(pseudotime_ranks), expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(modality~gene, scales='free')

p1 / p2 + plot_layout(heights=c(1,3))



#### Try diffmaps on clusters ####
# Get neuron-specific genes 
age_de <- de(rna_ret, groups = 'state')
age_de_marks <- age_de %>% 
    filter(feature%in%colnames(H3K27ac_clusters), group=='neuron')

vulcano_plot(age_de_marks, top_only = F)

traject_markers <- age_de_marks %>% 
    filter(fc!=0) %>% 
    group_by(sign(fc)) %>% 
    top_n(10, abs(fc)) %>% pull(feature) %>% unique()


# Compute diffmap for marks
# H3K27ac
H3K27ac_diffmap <- DiffusionMap(as.matrix(H3K27ac_clusters[, traject_markers]), verbose=T)
H3K27ac_dc_df <- as.data.frame(rank_cols(H3K27ac_diffmap@eigenvectors[,1:5])) %>% 
    as_tibble(rownames='H3K27ac_clusters')

H3K27ac_dc_df$pseudotime <- rank(-H3K27ac_dc_df$DC1) / max(H3K27ac_dc_df$DC1)

plot(H3K27ac_dc_df$pseudotime, H3K27ac_clusters[,'GRIA2'])
plot(H3K27ac_dc_df$pseudotime, H3K27ac_clusters[,'VIM'])

H3K27ac_pt_meta <- H3K27ac_dc_df %>% 
    inner_join(cluster_meta)

ggplot(H3K27ac_pt_meta, aes(pseudotime, pseudotime_ranks, color=DC1)) +
    geom_point(size=5)


# H3K27me3
H3K27me3_diffmap <- DiffusionMap(as.matrix(H3K27me3_clusters[, traject_markers]), verbose=T)
H3K27me3_dc_df <- as.data.frame(rank_cols(H3K27me3_diffmap@eigenvectors[,1:5])) %>% 
    as_tibble(rownames='H3K27me3_clusters')

H3K27me3_dc_df$pseudotime <- rank(H3K27me3_dc_df$DC1) / max(H3K27me3_dc_df$DC1)

plot(H3K27me3_dc_df$pseudotime, H3K27me3_clusters[,'STMN2'])

H3K27me3_pt_meta <- H3K27me3_dc_df %>% 
    inner_join(cluster_meta)

ggplot(H3K27me3_pt_meta, aes(pseudotime, pseudotime_ranks, color=DC1)) +
    geom_point(size=5)


# H3K4me3
H3K4me3_diffmap <- DiffusionMap(as.matrix(H3K4me3_clusters[, traject_markers]), verbose=T)
H3K4me3_dc_df <- as.data.frame(rank_cols(H3K4me3_diffmap@eigenvectors[,1:5])) %>% 
    as_tibble(rownames='H3K4me3_clusters')

H3K4me3_dc_df$pseudotime <- rank(-H3K4me3_dc_df$DC1) / max(H3K4me3_dc_df$DC1)

plot(H3K4me3_dc_df$pseudotime, H3K4me3_clusters[,'STMN2'])

H3K4me3_pt_meta <- H3K4me3_dc_df %>% 
    inner_join(cluster_meta)

ggplot(H3K4me3_pt_meta, aes(pseudotime, pseudotime_ranks, color=DC1)) +
    geom_point(size=5)


### Smooth
H3K27ac_expr <- H3K27ac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='H3K27ac_clusters', values_to='expr') %>% 
    inner_join(H3K27ac_pt_meta)

H3K27me3_expr <- H3K27me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='H3K27me3_clusters', values_to='expr') %>% 
    inner_join(H3K27me3_pt_meta)

H3K4me3_expr <- H3K4me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='H3K4me3_clusters', values_to='expr') %>% 
    inner_join(H3K4me3_pt_meta)

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='clusters', values_to='expr') %>% 
    inner_join(cluster_meta) %>% mutate(pseudotime=rank(pseudotime_ranks)/max(rank(pseudotime_ranks)))

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(
        expr=scale01(expr)
    )



p1 <- ggplot(plot_df, aes(pseudotime, expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(~gene, scales='free')

p2 <- ggplot(plot_df, aes(pseudotime, expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')) +
    facet_grid(modality~gene, scales='free')

p1 / p2 + plot_layout(heights=c(1,3))


















