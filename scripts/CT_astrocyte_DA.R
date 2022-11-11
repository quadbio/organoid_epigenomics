source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(Signac)
library(doParallel)
registerDoParallel(36)

setwd('~/projects/cutntag/')


feature_plot(rna, features=c('OPRM1', 'OPRK1', 'OPRD1', 'BTK'), order=T)


H3K27me3_ac_v_n <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_ACvsN.tsv')

peak_genes <-  ClosestFeature(marks$H3K27me3, StringToGRanges(H3K27me3_ac_v_n$feature)) %>% 
    as_tibble()

plot_df <- inner_join(H3K27me3_ac_v_n, peak_genes, by=c('feature'='query_region'))



ggplot(plot_df, aes(log_dr, -log10(pval), label=gene_name)) +
    geom_point() +
    geom_text_repel(data=filter(plot_df, pval<1e-20))

H3K27me3_ac <- read_rds('data/trajectories/astrocytes/H3K27me3_astrocytes_srt.rds')
feature_plot(marks$H3K27me3, features=c('EMX1', 'FOXG1', 'SOX2', 'CELF4', 'KCNK9', 'GALNT9', 'MYCN', 'SFRP1', 'PLCH2'), order=T)
dim_plot(marks$H3K27me3, group.by='celltype_jf')


H3K27me3_ac_v_npc <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_ACvsNPC.tsv')

peak_genes <-  ClosestFeature(marks$H3K27me3, StringToGRanges(H3K27me3_ac_v_npc$feature)) %>% 
    as_tibble()

plot_df <- inner_join(H3K27me3_ac_v_npc, peak_genes, by=c('feature'='query_region'))

ggplot(plot_df, aes(log_dr, -log10(pval), label=gene_name)) +
    geom_point() +
    geom_text_repel(data=filter(plot_df, pval<1e-20))


feature_plot(marks$H3K27me3, features=c('EMX1', 'FOXG1', 'SOX2', 'GRIA2', 'KCNK9', 'GALNT9', 'MYCN', 'SFRP1', 'PLCH2'), order=T)
feature_plot(rna, features=c('EMX1', 'FOXG1', 'SOX2', 'CELF4', 'KCNK9', 'GALNT9', 'MYCN', 'SFRP1', 'PLCH2'), order=T, reduction='cssumap')
dim_plot(marks$H3K27me3, group.by='celltype_jf')


#### Read and prepare stuff ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')

H3K27ac_ac <- read_rds('data/trajectories/astrocytes/H3K27ac_astrocytes_srt.rds')
H3K27me3_ac <- read_rds('data/trajectories/astrocytes/H3K27me3_astrocytes_srt.rds')
H3K4me3_ac <- read_rds('data/trajectories/astrocytes/H3K4me3_astrocytes_srt.rds')

marks_ac <- map(list('H3K27ac'=H3K27ac_ac, 'H3K27me3'=H3K27me3_ac, 'H3K4me3'=H3K4me3_ac), function(srt){
    srt@active.assay <- 'peaks'
    srt <- DietSeurat(srt, assays = c('peaks'))
    srt_peaks_bin <- as(srt@assays$peaks@counts>0, 'dgCMatrix')
    srt[['peaks_bin']] <- CreateAssayObject(srt_peaks_bin)
    return(srt)
})

marks_test <- map(marks, function(srt){
    srt@active.assay <- 'peaks'
    srt <- DietSeurat(srt, assays = c('peaks'))
    srt_peaks_bin <- as(srt@assays$peaks@counts>0, 'dgCMatrix')
    srt[['peaks_bin']] <- CreateAssayObject(srt_peaks_bin)
    return(srt)
})

mark_names <- set_names(names(marks_test))


#### N vs AC ####
map(mark_names, function(m){
    mark_use <- marks_test[[m]] %>%
        subset(celltype_jf%in%c('astrocytes', 'mesen_ex', 'dien_ex', 'ctx_ex', 'rhom_ex') & age%in%c('8mo', '128d', '60d'))

    peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
    features_detect <- names(peak_counts)[peak_counts>10]

    mark_use$test_var <- mark_use$celltype_jf == 'astrocytes'

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
        covariates = c('nCount_peaks'),
        family = 'binomial',
        slot = 'counts',
        assay = 'peaks_bin',
        features_use = features_use,
        parallel = T
    )
    test_results <- inner_join(test_df, detrate_df)

    test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_ACvsN.tsv'))
})


# #### NPC vs AC ####
# map(mark_names, function(m){
#     mark_use <- marks_test[[m]] %>%
#         subset(celltype_jf%in%c('astrocytes', 'nt_npc', 'dien_npc', 'ctx_npc') & age%in%c('8mo', '128d', '60d'))
# 
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
# 
#     mark_use$test_var <- mark_use$celltype_jf == 'astrocytes'
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
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_ACvsNPC.tsv'))
# })
# 
# 
# #### AC telen vs non-telen ####
# map(mark_names, function(m){
#     mark_use <- marks_ac[[m]]
#     
#     peak_counts <- rowSums(mark_use@assays$peaks_bin@counts)
#     features_detect <- names(peak_counts)[peak_counts>10]
#     
#     mark_use$test_var <- mark_use$ac_region == 'telen'
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
#     test_results %>% write_tsv(paste0('data/results/diff_expression/', m, '_DA_peaks_astrocytes_TvsNT.tsv'))
# })









