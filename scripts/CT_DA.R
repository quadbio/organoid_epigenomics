source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(doParallel)
registerDoParallel(36)

setwd('~/projects/cutntag/')

marks <- read_rds('data/CT/all_marks_list_v3.4lines.rds')

marks_test <- map(marks, function(srt){
    srt@active.assay <- 'peaks'
    srt <- DietSeurat(srt, assays = c('peaks'))
    srt_peaks_bin <- as(srt@assays$peaks@counts>0, 'dgCMatrix')
    srt[['peaks_bin']] <- CreateAssayObject(srt_peaks_bin)
    return(srt)
})

mark_names <- set_names(names(marks_test))


# #### All celltypes ####
# map(mark_names, function(m){
#     mark_use <- marks_test[[m]] 
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     test_vars <- set_names(unique(mark_use$celltype_jf))
#     
#     test_results <- map_dfr(test_vars, function(v){
#         print(v)
#         mark_use$test_var <- mark_use$celltype_jf == v
#         
#         detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#         detrate_df <- tibble(
#             feature = colnames(detection_rates),
#             detect_self = as.numeric(detection_rates['TRUE', ]),
#             detect_other = as.numeric(detection_rates['FALSE', ])
#         ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#         
#         features_use <- filter(detrate_df, abs(log_dr)>0.25, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#         
#         test_df <- lr_de(
#             object = mark_use,
#             test_var = 'test_var',
#             covariates = 'nCount_peaks',
#             family = 'binomial',
#             slot = 'counts',
#             assay = 'peaks_bin',
#             features_use = features_use,
#             parallel = T
#         )
#         return(inner_join(test_df, detrate_df))
#     }, .id='group')
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_celltype.tsv'))
# })
#  
# 
# 
# 
# #### All ages ####
# map(mark_names, function(m){
#     mark_use <- marks_test[[m]] 
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     test_vars <- set_names(unique(mark_use$age))
#     
#     test_results <- map_dfr(test_vars, function(v){
#         print(v)
#         mark_use$test_var <- mark_use$age == v
#         detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#         detrate_df <- tibble(
#             feature = colnames(detection_rates),
#             detect_self = as.numeric(detection_rates['TRUE', ]),
#             detect_other = as.numeric(detection_rates['FALSE', ])
#         ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#         
#         features_use <- filter(detrate_df, abs(log_dr)>0.25, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#         
#         test_df <- lr_de(
#             object = mark_use,
#             test_var = 'test_var',
#             covariates = 'nCount_peaks',
#             family = 'binomial',
#             slot = 'counts',
#             assay = 'peaks_bin',
#             features_use = features_use,
#             parallel = T
#         )
#         return(inner_join(test_df, detrate_df))
#     }, .id='group')
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_age.tsv'))
# })
#     
# 
# 
# #### Lineages all late ####
# map(mark_names, function(m){
#     mark_use <- subset(marks_test[[m]], stage!='EB' & stage!='mid' & lineage!='other')
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     test_vars <- set_names(unique(mark_use$lineage))
#     
#     test_results <- map_dfr(test_vars, function(v){
#         print(v)
#         mark_use$test_var <- mark_use$lineage == v
#         detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#         detrate_df <- tibble(
#             feature = colnames(detection_rates),
#             detect_self = as.numeric(detection_rates['TRUE', ]),
#             detect_other = as.numeric(detection_rates['FALSE', ])
#         ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#         
#         features_use <- filter(detrate_df, abs(log_dr)>0.25, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#         
#         test_df <- lr_de(
#             object = mark_use,
#             test_var = 'test_var',
#             covariates = 'nCount_peaks',
#             family = 'binomial',
#             slot = 'counts',
#             assay = 'peaks_bin',
#             features_use = features_use,
#             parallel = T
#         )
#         return(inner_join(test_df, detrate_df))
#     }, .id='group')
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_lineage_all.tsv'))
# })
#     
# 
# 
# #### Lineages neurons ####
# map(mark_names, function(m){
#     mark_use <- subset(marks_test[[m]], stage!='EB' & stage!='mid' & lineage!='other' & state=='neuron')
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     test_vars <- set_names(unique(mark_use$celltype_jf))
#     
#     test_results <- map_dfr(test_vars, function(v){
#         print(v)
#         mark_use$test_var <- mark_use$celltype_jf == v
#         detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#         detrate_df <- tibble(
#             feature = colnames(detection_rates),
#             detect_self = as.numeric(detection_rates['TRUE', ]),
#             detect_other = as.numeric(detection_rates['FALSE', ])
#         ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#         
#         features_use <- filter(detrate_df, abs(log_dr)>0.25, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#         
#         test_df <- lr_de(
#             object = mark_use,
#             test_var = 'test_var',
#             covariates = 'nCount_peaks',
#             family = 'binomial',
#             slot = 'counts',
#             assay = 'peaks_bin',
#             features_use = features_use,
#             parallel = T
#         )
#         return(inner_join(test_df, detrate_df))
#     }, .id='group')
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_lineage_neurons.tsv'))
# })
#     
# 
# 
# #### Lineages NPCs ####
# map(mark_names, function(m){
#     mark_use <- subset(marks_test[[m]], stage!='EB' & stage!='mid' & lineage!='other' & state=='npc')
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     test_vars <- set_names(unique(mark_use$lineage))
#     
#     test_results <- map_dfr(test_vars, function(v){
#         print(v)
#         mark_use$test_var <- mark_use$lineage == v
#         detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#         detrate_df <- tibble(
#             feature = colnames(detection_rates),
#             detect_self = as.numeric(detection_rates['TRUE', ]),
#             detect_other = as.numeric(detection_rates['FALSE', ])
#         ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#         
#         features_use <- filter(detrate_df, abs(log_dr)>0.25, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#         
#         test_df <- lr_de(
#             object = mark_use,
#             test_var = 'test_var',
#             covariates = 'nCount_peaks',
#             family = 'binomial',
#             slot = 'counts',
#             assay = 'peaks_bin',
#             features_use = features_use,
#             parallel = T
#         )
#         return(inner_join(test_df, detrate_df))
#     }, .id='group')
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_lineage_npcs.tsv'))
# })
# 
# 
# 
# #### Lineages + early ####
# map(mark_names, function(m){
#     mark_use <- marks_test[[m]] %>% subset(lineage%in%c('ctx', 'dien', 'nt', 'retina') | celltype_jf%in%c('psc', 'nect'))
#     mark_use$lineage_early <- case_when(
#         mark_use$lineage != 'other' ~ mark_use$lineage,
#         mark_use$celltype_jf %in% c('psc', 'nect') ~ mark_use$celltype_jf
#     )
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     test_vars <- set_names(unique(mark_use$lineage_early))
#     
#     test_results <- map_dfr(test_vars, function(v){
#         print(v)
#         mark_use$test_var <- mark_use$lineage_early == v
#         detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
#         detrate_df <- tibble(
#             feature = colnames(detection_rates),
#             detect_self = as.numeric(detection_rates['TRUE', ]),
#             detect_other = as.numeric(detection_rates['FALSE', ])
#         ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
#         
#         features_use <- filter(detrate_df, abs(log_dr)>0.25, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
#         
#         test_df <- lr_de(
#             object = mark_use,
#             test_var = 'test_var',
#             covariates = 'nCount_peaks',
#             family = 'binomial',
#             slot = 'counts',
#             assay = 'peaks_bin',
#             features_use = features_use,
#             parallel = T
#         )
#         return(inner_join(test_df, detrate_df))
#     }, .id='group')
#     
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_lineage_early.tsv'))
# })
   

#### Lineages + early + other stuff ####
map(mark_names, function(m){
    mark_use <- marks_test[[m]]
    mark_use$lineage_early <- case_when(
        mark_use$lineage != 'other' ~ mark_use$lineage,
        T ~ mark_use$celltype_jf
    )
    
    peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
    features_detect <- names(peak_counts)[peak_counts>10]
    
    test_vars <- set_names(unique(mark_use$lineage_early))
    
    test_results <- map_dfr(test_vars, function(v){
        print(v)
        mark_use$test_var <- mark_use$lineage_early == v
        detection_rates <- Pando::aggregate_matrix(t(mark_use[['peaks_bin']]@counts), mark_use$test_var)
        detrate_df <- tibble(
            feature = colnames(detection_rates),
            detect_self = as.numeric(detection_rates['TRUE', ]),
            detect_other = as.numeric(detection_rates['FALSE', ])
        ) %>% mutate(log_dr=log2(detect_self) - log2(detect_other))
        
        features_use <- filter(detrate_df, abs(log_dr)>0.25, feature%in%features_detect, detect_self>0.005 | detect_other>0.005)$feature
        
        test_df <- lr_de(
            object = mark_use,
            test_var = 'test_var',
            covariates = 'nCount_peaks',
            family = 'binomial',
            slot = 'counts',
            assay = 'peaks_bin',
            features_use = features_use,
            parallel = T
        )
        return(inner_join(test_df, detrate_df))
    }, .id='group')
    
    test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_lineage_coarse.tsv'))
})
    








