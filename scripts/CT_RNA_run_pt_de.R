source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(destiny)
library(dtw)

library(doParallel)
registerDoParallel(20)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')

mark_colors <- c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')


#### Read stuff ####
H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')

rna_ctx <- FindVariableFeatures(rna_ctx, nfeatures=2000)
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


genes_test <- intersect(mark_detect_feats, VariableFeatures(rna_ctx))


H3K27ac_ctx$nCount_RNA <- H3K27ac_ctx$nCount_cRNA
H3K27me3_ctx$nCount_RNA <- H3K27me3_ctx$nCount_cRNA
H3K4me3_ctx$nCount_RNA <- H3K4me3_ctx$nCount_cRNA

H3K27ac_ctx[['RNA']] <- H3K27ac_ctx[['cRNA']]
H3K27me3_ctx[['RNA']] <- H3K27me3_ctx[['cRNA']]
H3K4me3_ctx[['RNA']] <- H3K4me3_ctx[['cRNA']]


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

rna_binom_pt_test <- test_pt_bins(rna_ctx, features=genes_test, assay='RNA_bin')
H3K27ac_binom_pt_test <- test_pt_bins(H3K27ac_ctx, features=genes_test, assay='RNA_bin')
H3K27me3_binom_pt_test <- test_pt_bins(H3K27me3_ctx, features=genes_test, assay='RNA_bin')
H3K4me3_binom_pt_test <- test_pt_bins(H3K4me3_ctx, features=genes_test, assay='RNA_bin')

plot_df <- bind_rows('RNA'=rna_binom_pt_test, 'H3K27ac'=H3K27ac_binom_pt_test, 'H3K27me3'=H3K27me3_binom_pt_test, 'H3K4me3'=H3K4me3_binom_pt_test, .id='modality') %>% 
    group_by(modality) %>% mutate(padj=p.adjust(pval, 'fdr'))

plot_df %>% write_tsv('data/trajectories/ctx/ctx_pt_lr_binom_de.tsv')

