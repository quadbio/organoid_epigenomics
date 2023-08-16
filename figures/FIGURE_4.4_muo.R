source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')


#### Read data ####
marks <- read_rds('data05/all_marks_list_v3.4lines.rds')
rna <- read_rds('data05/RNA/RNA_all_srt_v2.3lines.rds')

annotation <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')

data_paths <- list(
    'multiome_org_WHB_120d_rep1'='/links/groups/treutlein/DATA/sequencing/20230609_P2454_FIDES/processed/multiome_org_WHB_120d_rep1/outs/', 
    'multiome_org_WHB_120d_rep2'='/links/groups/treutlein/DATA/sequencing/20230609_P2454_FIDES/processed/multiome_org_WHB_120d_rep2/outs/'
)

count_list <- map(data_paths, function(x){
    Read10X_h5(paste0(x, 'filtered_feature_bc_matrix.h5'))
})

srt_list <- map(set_names(names(count_list)), function(n){
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

muo_srt <- merge(srt_list[[1]], srt_list[2:length(srt_list)], add.cell.ids=names(srt_list))

muo_srt@active.assay <- 'RNA'
muo_srt <- muo_srt %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

ElbowPlot(muo_srt)

muo_srt <- muo_srt %>% 
    RunUMAP(dims=1:20)

muo_srt$log_rna_counts <- log10(muo_srt$nCount_RNA)
muo_srt$log_atac_counts <- log10(muo_srt$nCount_ATAC)

muo_srt %>% write_rds('data05/MUO/MUO_120d_merged.rds')

muo_srt$orig.ident <- colnames(muo_srt) %>% str_replace('(.+)_[ATCG]+-\\d', '\\1')

muo_srt <- muo_srt %>% FindNeighbors()
muo_srt <- muo_srt %>% FindClusters()

muo_srt %>% dim_plot(group.by=c('seurat_clusters'))
muo_srt %>% feature_plot(features=c('NEUROD6', 'SOX2', 'EOMES', 'OTX2', 'HOXB2', 'RSPO3', 'TCF7L2', 'STMN2', 'GRIA2'), order=T)

muo_srt$log_rna_counts %>% hist()
muo_srt$log_atac_counts %>% hist()

muo_srt_fltr <- subset(
    muo_srt,
    nCount_RNA > 100 &
        nCount_ATAC > 2000
)

muo_srt@active.assay <- 'RNA'
muo_srt_fltr <- muo_srt_fltr %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    RunUMAP(dims=1:20)

muo_srt_fltr@active.assay <- 'ATAC'
muo_srt_fltr <- muo_srt_fltr %>% 
    FindTopFeatures() %>% 
    RunTFIDF() %>% 
    RunSVD()

ElbowPlot(muo_srt_fltr, reduction = 'lsi')

muo_srt_fltr <- muo_srt_fltr %>% 
    RunUMAP(dims=2:20, reduction='lsi', reduction.name='lsiumap')

muo_srt_fltr %>% dim_plot(group.by=c('seurat_clusters'), reduction='lsiumap')
muo_srt_fltr %>% feature_plot(features=c('NEUROD6', 'SOX2', 'EOMES', 'NFIX', 'EMX1', 'RSPO3', 'GLI3', 'MKI67', 'GRIA2'), order=T, reduction='lsiumap')


#### Update them frags ####
frags <- Fragments(muo_srt_fltr)
new_frags1 <- UpdatePath(frags[[1]], '/local1/scratch/fleckj/projects/cutntag/data05/multiome_org_WHB_120d_rep1/atac_fragments.tsv.gz')
new_frags2 <- UpdatePath(frags[[2]], '/local1/scratch/fleckj/projects/cutntag/data05/multiome_org_WHB_120d_rep2/atac_fragments.tsv.gz')
frags[[1]] <- NULL
frags[[1]] <- new_frags1
frags[[2]] <- NULL
frags[[2]] <- new_frags2
muo_srt_fltr <- SetAssayData(muo_srt_fltr, slot='fragments', new.data=frags)

gene_activities <- GeneActivity(
    muo_srt_fltr, 
    assay = 'ATAC',
    features = rownames(muo_srt_fltr[['RNA']])
)

muo_srt_fltr[['gene_activity']] <- CreateAssayObject(gene_activities)

muo_srt_fltr@active.assay <- 'gene_activity'
muo_srt_fltr <- muo_srt_fltr %>% NormalizeData()

muo_srt_fltr@active.assay <- 'gene_activity'
p1 <- muo_srt_fltr %>% feature_plot(features=c('NEUROD6', 'SOX2', 'EOMES', 'NFIX', 'GLI3', 'NFIA', 'HOPX', 'STMN2', 'EMX1'), order=T, reduction='lsiumap')
muo_srt_fltr@active.assay <- 'RNA'
p2 <- muo_srt_fltr %>% feature_plot(features=c('NEUROD6', 'SOX2', 'EOMES', 'NFIX', 'GLI3', 'NFIA', 'HOPX', 'STMN2', 'EMX1'), order=T, reduction='lsiumap')
p1 / p2

muo_srt_fltr@active.assay <- 'gene_activity'
p1 <- muo_srt_fltr %>% feature_plot(features=c('NEUROD6', 'SOX2', 'EOMES', 'NFIX', 'GLI3', 'NFIA', 'HOPX', 'STMN2', 'EMX1', 'AQP4'), order=T, reduction='umap')
muo_srt_fltr@active.assay <- 'RNA'
p2 <- muo_srt_fltr %>% feature_plot(features=c('NEUROD6', 'SOX2', 'EOMES', 'NFIX', 'GLI3', 'NFIA', 'HOPX', 'STMN2', 'EMX1', 'AQP4'), order=T, reduction='umap')
p1 / p2


muo_srt_fltr <- SCTransform(muo_srt_fltr)
muo_srt_fltr %>% write_rds('data05/MUO/MUO_120d_v1preproc.rds')



#### Integrate with other data ####
muo_srt_fltr@active.assay <- 'RNA'
muo_srt_rna <- DietSeurat(muo_srt_fltr, assays=c('RNA'))

anchors <- FindIntegrationAnchors(list(muo_srt_rna, rna))
anchors %>% write_rds('data05/MUO/MOU_RNA_integration_anchors.rds')

muo_rna_integrated <- IntegrateData(anchors)

muo_rna_integrated@active.assay <- 'integrated'
muo_rna_integrated <- muo_rna_integrated %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

muo_rna_integrated <- muo_rna_integrated %>% 
    RunUMAP(dims=1:20, reduction.name='ccaumap')

muo_rna_integrated %>% feature_plot(features=c('NEUROD6', 'EMX1', 'GRIA2', 'NFIA', 'GLI3', 'NFIX', 'AQP4'), order=T, reduction='ccaumap')
muo_rna_integrated %>% dim_plot(group.by=c('orig.ident'), reduction='ccaumap')
muo_rna_integrated %>% dim_plot(group.by=c('celltype_jf'), reduction='ccaumap')

muo_rna_integrated %>% write_rds('data05/MUO/MUO_120d_v1_RNA_v3.4_integrated.rds')

tanchors <- FindTransferAnchors(reference=rna, query=muo_srt_rna, reduction='pcaproject', k.filter=100, dims=1:20)
tanchors %>% write_rds('data05/MUO/MOU_RNA_transfer_anchors.rds')

celltype_transfer <- TransferData(tanchors, refdata='celltype_jf', reference=rna, k.weight=10)
muo_srt_fltr$celltype_transfer <- celltype_transfer$predicted.id

muo_srt_fltr %>% dim_plot(group.by='celltype_transfer', reduction='umap')



#### Try subset ctx trajectory ####
muo_srt_fltr <- muo_srt_fltr %>% subset(nFeature_ATAC>2000 & nFeature_RNA>100)
muo_srt_fltr <- muo_srt_fltr %>% subset(nFeature_ATAC<40000)

muo_srt_fltr@active.assay <- 'RNA'
muo_srt_fltr <- muo_srt_fltr %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    RunUMAP(dims=1:20)

muo_srt_fltr <- muo_srt_fltr %>% FindNeighbors(dims=1:20, reduction='pca')
muo_srt_fltr <- muo_srt_fltr %>% FindClusters()

muo_srt_fltr %>% dim_plot(group.by=c('seurat_clusters'), reduction='umap', label=T)
muo_srt_fltr %>% feature_plot(features=c('TCF7L2', 'OTX2', 'NEUROD6', 'EMX1', 'GRIA2', 'GLI3', 'AQP4', 'MKI67', 'FOXG1'), reduction='umap', order=T)

muo_srt_fltr %>% dim_plot(group.by=c('seurat_clusters'), reduction='lsiumap', label=T)
muo_srt_fltr %>% feature_plot(features=c('nFeature_RNA'), reduction='lsiumap')

muo_srt_ctx_pre <- muo_srt_fltr %>% subset(seurat_clusters%in%c(12,9,3,11,5,1,4,2))

muo_srt_ctx_pre <- muo_srt_ctx_pre %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    FindNeighbors(dims=1:20, reduction='pca') %>% 
    RunUMAP(dims=1:20, reduction='pca')

muo_srt_ctx_pre <- muo_srt_ctx_pre %>% FindClusters()

muo_srt_ctx_pre %>% dim_plot(group.by=c('seurat_clusters'), reduction='umap')

muo_srt_ctx_pre@active.assay <- 'RNA'
p1 <- muo_srt_ctx_pre %>% feature_plot(features=c('SATB2', 'MKI67', 'STMN2', 'EOMES', 'EMX1', 'SOX2', 'BCL11B', 'NEUROD2', 'GLI3'), reduction='umap', order=T)

muo_srt_ctx_pre@active.assay <- 'gene_activity'
p2 <- muo_srt_ctx_pre %>% feature_plot(features=c('SATB2', 'MKI67', 'STMN2', 'EOMES', 'EMX1', 'SOX2', 'BCL11B', 'NEUROD2', 'GLI3'), reduction='umap', order=T)

p1 / p2

muo_srt_ctx_pre@active.assay <- 'ATAC'
muo_srt_ctx_pre <- muo_srt_ctx_pre %>% 
    FindTopFeatures() %>% 
    RunTFIDF() %>% 
    RunSVD() %>% 
    FindNeighbors(dims=2:20, reduction='lsi') %>% 
    RunUMAP(dims=2:20, reduction='lsi')

muo_srt_ctx_pre %>% feature_plot(features=c('nFeature_ATAC'), reduction='umap', order=T)


#### Infer pseudotime ####
pca_dims <- 1:10

muo_srt_ctx_pre@active.assay <- 'RNA'
muo_srt_ctx <- muo_srt_ctx_pre %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

muo_srt_ctx <- muo_srt_ctx %>% 
    RunUMAP(dims=pca_dims, reduction='pca')

muo_srt_ctx <- muo_srt_ctx %>% 
    FindNeighbors(reduction='pca') %>% 
    FindClusters()

ElbowPlot(muo_srt_ctx, reduction='pca')
dim_plot(muo_srt_ctx, label=T, group.by=c('seurat_clusters'))
feature_plot(muo_srt_ctx, features=c('NEUROD6', 'MKI67', 'STMN2', 'EOMES', 'EMX1', 'SOX2'), order=T)

rna_muo_pca <- muo_srt_ctx[['pca']]@cell.embeddings[,pca_dims]
rna_muo_diffmap <- DiffusionMap(rna_muo_pca, verbose=T)
rna_muo_dc_df <- as.data.frame(rank_cols(rna_muo_diffmap@eigenvectors[,1:5]))

muo_srt_ctx[['diffmap']] <- CreateDimReducObject(
    rna_muo_diffmap@eigenvectors*100, key='DC_', assay='RNA'
)
muo_srt_ctx <- AddMetaData(muo_srt_ctx, rna_muo_dc_df)

muo_srt_ctx %>% 
    FeaturePlot(features=c('STMN2', 'VIM', 'DC1', 'DC2', 'DC3'), order=T)

muo_srt_ctx %>% 
    dim_plot(group.by=c('stage', 'celltype_jf', 'seurat_clusters'), order=T)


muo_srt_ctx %>% write_rds('data_/trajectories/ctx/MOU_RNA_4m_ctx_dpt_srt.rds')

muo_srt_ctx %>% 
    feature_plot(features=c('STMN2', 'SOX2', 'VIM', 'DC1', 'DC2', 'ctx_pt'), order=T)



#### Bin pt ####
muo_srt_ctx <- read_rds('data_/trajectories/ctx/MOU_RNA_4m_ctx_dpt_srt.rds')

muo_srt_ctx$ctx_pt <- rank(-muo_srt_ctx$DC1) / max(rank(muo_srt_ctx$DC1))
muo_srt_ctx$pt_bins2 <- as.numeric(cut(muo_srt_ctx$ctx_pt, 10, labels=1:10))
muo_srt_ctx$pt_bins <- as.numeric(cut(muo_srt_ctx$ctx_pt, 20, labels=1:20))
muo_srt_ctx <- Pando::aggregate_assay(muo_srt_ctx, assay='RNA', group_name='pt_bins2')
muo_srt_ctx <- Pando::aggregate_assay(muo_srt_ctx, assay='RNA', group_name='pt_bins')
muo_srt_ctx[['RNA_bin']] <- CreateAssayObject((muo_srt_ctx[['RNA']]@data > 0)*1)
muo_srt_ctx <- Pando::aggregate_assay(muo_srt_ctx, assay='RNA_bin', group_name='pt_bins2')
muo_srt_ctx <- Pando::aggregate_assay(muo_srt_ctx, assay='RNA_bin', group_name='pt_bins')

muo_srt_ctx <- Pando::aggregate_assay(muo_srt_ctx, assay='gene_activity', group_name='pt_bins2')
muo_srt_ctx <- Pando::aggregate_assay(muo_srt_ctx, assay='gene_activity', group_name='pt_bins')
muo_srt_ctx[['gene_activity_bin']] <- CreateAssayObject((muo_srt_ctx[['gene_activity']]@data > 0)*1)
muo_srt_ctx <- Pando::aggregate_assay(muo_srt_ctx, assay='gene_activity_bin', group_name='pt_bins2')
muo_srt_ctx <- Pando::aggregate_assay(muo_srt_ctx, assay='gene_activity_bin', group_name='pt_bins')

muo_ga_clusters <- muo_srt_ctx@assays$gene_activity_bin@misc$summary$pt_bins2[as.character(1:10), ]
muo_rna_clusters <- muo_srt_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]



#### Load other pseudotimes ####
H3K27ac_ctx <- read_rds('data_/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data_/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data_/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data_/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')

H3K27ac_clusters <- H3K27ac_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K27me3_clusters <- H3K27me3_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K4me3_clusters <- H3K4me3_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
rna_clusters <- rna_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]

genes_plot <- c('NEUROD6', 'NEUROD2', 'GRIA2', 'STMN2', 'VIM', 'SOX2', 'GLI3', 'EOMES', 'NES')

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

muo_rna_expr <- muo_rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins2', values_to='expr')

muo_ga_expr <- muo_ga_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins2', values_to='expr')

all_mod_pt <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, 'MUO_RNA'=muo_rna_expr, 'MUO_ATAC'=muo_ga_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr01=scale01(expr))

modality_colors <- c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044', 'MUO_ATAC'='#9575cd', 'MUO_RNA'='#FDDC44')


plot_df <- all_mod_pt
p1 <- ggplot(plot_df, aes(as.numeric(pt_bins2), expr01, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se = F) +
    geom_point(size=0.2) +
    scale_color_manual(values=modality_cols) +
    facet_grid(~gene, scales='free')

p2 <- ggplot(plot_df, aes(as.numeric(pt_bins2), expr01, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=modality_cols) +
    facet_grid(modality~gene, scales='free')

p1 / p2 + plot_layout(heights=c(1,3))
ggsave('plots/paper/fig3/MUO_4m_ctx_pt_all_mods_smooth.png', width=30, height=16, units='cm')
ggsave('plots/paper/fig3/MUO_4m_ctx_pt_all_mods_smooth.pdf', width=30, height=16, units='cm')


plot_df <- all_mod_pt %>% 
    filter(modality%in%c('RNA', 'MUO_RNA'))
ggplot(plot_df, aes(as.numeric(pt_bins2), expr01, color=modality, group=modality)) +
    geom_bar(stat='identity', fill='grey', alpha=0.6, position='dodge', linewidth=0.2) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    scale_color_manual(values=modality_cols) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_grid(~gene, scales='free') +
    article_text() +
    labs(x='Pseudotime bins', y='Scaled detection rate')




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
MUO_RNA_binom_pt_test <- test_pt_bins(muo_srt_ctx, features=genes_test, assay='RNA_bin')
MUO_ATAC_binom_pt_test <- test_pt_bins(muo_srt_ctx, features=genes_test, assay='gene_activity_bin')


all_lr_pt_de <- bind_rows('RNA'=rna_binom_pt_test, 'H3K27ac'=H3K27ac_binom_pt_test, 'H3K27me3'=H3K27me3_binom_pt_test, 'H3K4me3'=H3K4me3_binom_pt_test, 'MUO_RNA'=MUO_RNA_binom_pt_test, 'MUO_ATAC'=MUO_ATAC_binom_pt_test, .id='modality') %>% 
    group_by(modality) %>% mutate(padj=p.adjust(pval, 'fdr'))

# all_lr_pt_de %>% write_tsv('data_/trajectories/ctx/ctx_all_neuonal_pt_lr_binom_de.tsv')
all_lr_pt_de <- read_tsv('data_/trajectories/ctx/ctx_all_neuonal_pt_lr_binom_de.tsv')

all_lr_pt_diff <- all_lr_pt_de %>% 
    filter(padj<0.05, log_dr>0.25, modality!='H3K27me3') %>% 
    group_by(feature) %>% 
    # filter(length(unique(modality))>1) %>%
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

ggsave('plots/paper/fig3/MUO_lag_quantification_neuron_genes_boxplot.pdf', width=6, height=3, units='cm')



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

ATAC_lag <- all_lr_pt_diff %>% 
    filter(modality!='MUO_ATAC') %>% 
    select(modality, feature, min_bin) %>% 
    pivot_wider(names_from='modality', values_from='min_bin') %>% 
    mutate(bin_diff=RNA-H3K4me3) %>% filter(!is.na(bin_diff)) %>% select(feature, bin_diff)

plot_df <- bind_rows('H3K27ac'=H3K27ac_lag, 'H3K4me3'=H3K4me3_lag, 'MUO_ATAC'=ATAC_lag, .id = 'mark')

ggplot(plot_df, aes(bin_diff, mark, fill=mark)) + 
    geom_density_ridges() +
    geom_vline(xintercept=c(-2,2)) +
    scale_fill_manual(values=modality_colors) +
    article_text() + no_legend() +
    labs(x='Pseudotime bin lag relative to RNA', y='Mark') 

ggsave('plots/paper/fig3/MUO_lag_quantification_neuron_genes_dridges.pdf', width=6, height=3, units='cm')




#### Quantify lag for initial neuron genes ####
gene_clusters <- read_tsv('data_/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')
genes_test <- gene_clusters %>% filter(dtw_clust==27) %>% pull(feature) %>% unique()

library(doParallel)
registerDoParallel(20)

rna_binom_pt_test <- test_pt_bins(rna_ctx, features=genes_test, assay='RNA_bin')
H3K27ac_binom_pt_test <- test_pt_bins(H3K27ac_ctx, features=genes_test, assay='RNA_bin')
H3K27me3_binom_pt_test <- test_pt_bins(H3K27me3_ctx, features=genes_test, assay='RNA_bin')
H3K4me3_binom_pt_test <- test_pt_bins(H3K4me3_ctx, features=genes_test, assay='RNA_bin')
MUO_RNA_binom_pt_test <- test_pt_bins(muo_srt_ctx, features=genes_test, assay='RNA_bin')
MUO_ATAC_binom_pt_test <- test_pt_bins(muo_srt_ctx, features=genes_test, assay='gene_activity_bin')

all_lr_pt_de <- bind_rows('RNA'=rna_binom_pt_test, 'H3K27ac'=H3K27ac_binom_pt_test, 'H3K27me3'=H3K27me3_binom_pt_test, 'H3K4me3'=H3K4me3_binom_pt_test, 'MUO_RNA'=MUO_RNA_binom_pt_test, 'MUO_ATAC'=MUO_ATAC_binom_pt_test, .id='modality') %>% 
    group_by(modality) %>% mutate(padj=p.adjust(pval, 'fdr'))

all_lr_pt_de %>% write_tsv('data_/trajectories/ctx/ctx_neuronal_cluster_pt_lr_binom_de.tsv')

all_lr_pt_diff <- all_lr_pt_de %>% 
    filter(padj<0.05, log_dr>0.2, modality!='H3K27me3') %>% 
    group_by(feature) %>% 
    # filter(length(unique(modality))>1) %>%
    group_by(modality, feature) %>% 
    dplyr::summarize(min_bin=min(as.numeric(pt_bin)), min_dr=log_dr[which.min(pt_bin)]) %>% 
    dplyr::mutate(signed_bin=sign(min_dr)*min_bin) 

ggplot(all_lr_pt_diff, aes(min_bin, modality, fill=modality, color=modality)) +
    geom_boxplot(linewidth=0.2, color='black', alpha=0.5, outlier.shape=NA) +
    geom_quasirandom(size=0.5) +
    scale_fill_manual(values=modality_colors) + 
    scale_color_manual(values=modality_colors) + 
    scale_x_continuous(limits=c(0,10)) + 
    article_text() + no_legend() +
    labs(x='First divergent bin', y='Modality')

ggsave('plots/paper/fig3/MUO_lag_quantification_neuron_cluster_genes_boxplot.pdf', width=6, height=3, units='cm')




#### Jitter with also ATAC ####

genes_plot <- c('GRIA2', 'STMN2', 'NEUROD2', 'SOX2', 'VIM')

purples2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', rev(pals::magma(9))[3:9]), bias=bias)(100))
}


ga_expr_plot <- muo_srt_ctx[['gene_activity']]@data[genes_plot, ] %>%
    as_tibble(rownames='feature') %>%
    pivot_longer(!feature, names_to='cell', values_to='expr') %>% 
    inner_join(as_tibble(muo_srt_ctx@meta.data, rownames='cell')) 

ggplot(arrange(ga_expr_plot, expr), aes(ctx_pt, feature, color=expr)) +
    geom_jitter(size=0.1, height=0.3) +
    # scale_color_gradientn(colors=grad(pals::brewer.ylgnbu, 1.5)) +
    scale_color_gradientn(colors=purples2()) +
    article_text() +
    labs(x='Pseudotime', y='Feature', color='Log-normalized counts')

ggsave('plots/paper/fig3/fig3_ctx_2m_pt_expr_jitter_atac_magma.png', width=10, height=5, units='cm', bg='white')
ggsave('plots/paper/fig3/fig3_ctx_2m_pt_expr_jitter_atac_magma.pdf', width=10, height=5, units='cm', bg='white')

















