source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(doParallel)
registerDoParallel(36)

setwd('~/projects/cutntag/')


#### Read and prepare stuff ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')

H3K27ac_nt <- read_rds('data/trajectories/nt/H3K27ac_nt_dpt_lsi_regress_srt.rds')
H3K4me3_nt <- read_rds('data/trajectories/nt/H3K4me3_nt_dpt_lsi_regress_srt.rds')
H3K27me3_nt <- read_rds('data/trajectories/nt/H3K27me3_nt_dpt_lsi_regress_srt.rds')

H3K27ac_dien <- read_rds('data/trajectories/dien/H3K27ac_dien_dpt_lsi_regress_srt.rds')
H3K4me3_dien <- read_rds('data/trajectories/dien/H3K4me3_dien_dpt_lsi_regress_srt.rds')
H3K27me3_dien <- read_rds('data/trajectories/dien/H3K27me3_dien_dpt_lsi_regress_srt.rds')

H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_ctx_neuro_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_ctx_neuro_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_ctx_neuro_dpt_srt.rds')


marks_nt <- map(list('H3K27ac'=H3K27ac_nt, 'H3K27me3'=H3K27me3_nt, 'H3K4me3'=H3K4me3_nt), function(srt){
    srt@active.assay <- 'peaks'
    srt <- DietSeurat(srt, assays = c('peaks'))
    srt_peaks_bin <- as(srt@assays$peaks@counts>0, 'dgCMatrix')
    srt[['peaks_bin']] <- CreateAssayObject(srt_peaks_bin)
    return(srt)
})

marks_dien <- map(list('H3K27ac'=H3K27ac_dien, 'H3K27me3'=H3K27me3_dien, 'H3K4me3'=H3K4me3_dien), function(srt){
    srt@active.assay <- 'peaks'
    srt <- DietSeurat(srt, assays = c('peaks'))
    srt_peaks_bin <- as(srt@assays$peaks@counts>0, 'dgCMatrix')
    srt[['peaks_bin']] <- CreateAssayObject(srt_peaks_bin)
    return(srt)
})

marks_ctx <- map(list('H3K27ac'=H3K27ac_ctx, 'H3K27me3'=H3K27me3_ctx, 'H3K4me3'=H3K4me3_ctx), function(srt){
    srt@active.assay <- 'peaks'
    srt <- DietSeurat(srt, assays = c('peaks'))
    srt_peaks_bin <- as(srt@assays$peaks@counts>0, 'dgCMatrix')
    srt[['peaks_bin']] <- CreateAssayObject(srt_peaks_bin)
    return(srt)
})

mark_names <- set_names(names(marks_ctx))

new_neuron_annot <- map(mark_names, function(m){
    meta <- rbind(marks_dien[[m]]@meta.data['neuron'], marks_nt[[m]]@meta.data['neuron'], marks_ctx[[m]]@meta.data['neuron'])
    colnames(meta)[1] <- 'neuron'
    return(meta)
})

marks_test <- map(mark_names, function(m){
    srt <- marks[[m]]
    srt@active.assay <- 'peaks'
    srt <- DietSeurat(srt, assays = c('peaks'))
    srt_peaks_bin <- as(srt@assays$peaks@counts>0, 'dgCMatrix')
    srt[['peaks_bin']] <- CreateAssayObject(srt_peaks_bin)
    srt <- AddMetaData(srt, new_neuron_annot[[m]])
    srt$neuron <- ifelse(is.na(srt$neuron), 'NA', srt$neuron)
    return(srt)
})




# #### NPC vs neuron ####
# map(mark_names, function(m){
#     mark_use <- marks_nt[[m]]
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     mark_use$test_var <- mark_use$neuron
#     
#     detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#     detrate_df <- tibble(
#         feature = colnames(detection_rates),
#         detect_self = as.numeric(detection_rates['TRUE', ]),
#         detect_other = as.numeric(detection_rates['FALSE', ])
#     ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#     
#     features_use <- filter(detrate_df, abs(log_dr)>0.1, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#     
#     test_df <- lr_de(
#         object = mark_use,
#         test_var = 'test_var',
#         covariates = c('nCount_peaks'),
#         family = 'binomial',
#         slot = 'counts',
#         assay = 'peaks_bin',
#         features_use = features_use,
#         parallel = T
#     )
#     test_results <- inner_join(test_df, detrate_df)
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_nt_NvsNPC.tsv'))
# })
# 
# 
# 
# map(mark_names, function(m){
#     mark_use <- marks_ctx[[m]]
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     mark_use$test_var <- mark_use$neuron
#     
#     detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#     detrate_df <- tibble(
#         feature = colnames(detection_rates),
#         detect_self = as.numeric(detection_rates['TRUE', ]),
#         detect_other = as.numeric(detection_rates['FALSE', ])
#     ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#     
#     features_use <- filter(detrate_df, abs(log_dr)>0.1, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#     
#     test_df <- lr_de(
#         object = mark_use,
#         test_var = 'test_var',
#         covariates = c('nCount_peaks'),
#         family = 'binomial',
#         slot = 'counts',
#         assay = 'peaks_bin',
#         features_use = features_use,
#         parallel = T
#     )
#     test_results <- inner_join(test_df, detrate_df)
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_ctx_NvsNPC.tsv'))
# })
# 
# 
# 
# map(mark_names, function(m){
#     mark_use <- marks_dien[[m]]
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     mark_use$test_var <- mark_use$neuron
#     
#     detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#     detrate_df <- tibble(
#         feature = colnames(detection_rates),
#         detect_self = as.numeric(detection_rates['TRUE', ]),
#         detect_other = as.numeric(detection_rates['FALSE', ])
#     ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#     
#     features_use <- filter(detrate_df, abs(log_dr)>0.1, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#     
#     test_df <- lr_de(
#         object = mark_use,
#         test_var = 'test_var',
#         covariates = c('nCount_peaks'),
#         family = 'binomial',
#         slot = 'counts',
#         assay = 'peaks_bin',
#         features_use = features_use,
#         parallel = T
#     )
#     test_results <- inner_join(test_df, detrate_df)
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_dien_NvsNPC.tsv'))
# })

map(mark_names, function(m){
    mark_use <- subset(marks_test[[m]], neuron!='NA')

    peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
    features_detect <- names(peak_counts)[peak_counts>10]

    mark_use$test_var <- mark_use$neuron=='TRUE'

    detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
    detrate_df <- tibble(
        feature = colnames(detection_rates),
        detect_self = as.numeric(detection_rates['TRUE', ]),
        detect_other = as.numeric(detection_rates['FALSE', ])
    ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))

    features_use <- filter(detrate_df, abs(log_dr)>0.1, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature

    test_df <- lr_de(
        object = mark_use,
        test_var = 'test_var',
        covariates = c('nCount_peaks', 'log_peak_counts'),
        family = 'binomial',
        slot = 'counts',
        assay = 'peaks_bin',
        features_use = features_use,
        parallel = T
    )
    test_results <- inner_join(test_df, detrate_df)

    test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_all_NvsNPC.tsv'))
})


#### Lineage 1 vs all ####
# lineage_combs <- list(c('nt', 'dien'), c('nt', 'ctx'), c('dien', 'ctx'))
# map(mark_names, function(m){
#     mark_here <- marks_test[[m]] %>%
#         subset(lineage%in%c('nt', 'dien', 'ctx') & age%in%c('35d', '60d', '128d'))
#     
#     peak_counts <- rowSums(mark_here@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     test_results <- map_dfr(lineage_combs, function(x){
#         mark_use <- mark_here %>% subset(lineage%in%x)
#         
#         mark_use$test_var <- mark_use$lineage == x[2]
#         
#         detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#         detrate_df <- tibble(
#             feature = colnames(detection_rates),
#             detect_self = as.numeric(detection_rates['TRUE', ]),
#             detect_other = as.numeric(detection_rates['FALSE', ])
#         ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#         
#         features_use <- filter(detrate_df, abs(log_dr)>0.1, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#         
#         test_df <- lr_de(
#             object = mark_use,
#             test_var = 'test_var',
#             covariates = c('nCount_peaks'),
#             family = 'binomial',
#             slot = 'counts',
#             assay = 'peaks_bin',
#             features_use = features_use,
#             parallel = T
#         )
#         tres <- inner_join(test_df, detrate_df)
#         tres$group1 <- x[2]
#         tres$group2 <- x[1]
#         return(dplyr::select(tres, group1, group2, everything()))
#     })
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_lineages_pairwise.tsv'))
# })









