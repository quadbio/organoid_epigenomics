source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')

library(Pando)

filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist

setwd('~/projects/cutntag/')

#### Read data ####
drugs_all <- read_rds('data/drugs/drugs_all_v1.2lines_srt.rds')

dim_plot(drugs_all, group.by=c('orig.ident', 'DMSO', 'inhib_annotation', 'HTO_classification'))
feature_plot(drugs_all, features=c('nFeature_RNA'))


#### Subset and integrate ####
#### d15 ####
drugs_d15 <- drugs_all %>% subset(orig.ident=='day_15')

drugs_d15 <- drugs_d15 %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

ElbowPlot(drugs_d15)

drugs_d15 <- RunUMAP(
    drugs_d15, 
    reduction='pca', 
    dims=1:20,
    reduction.name='umap'
)

p1 <- dim_plot(drugs_d15, group.by=c('DMSO', 'inhib_annotation'))
p2 <- feature_plot(drugs_d15, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T)
p_nat <- p1 / p2 + plot_layout(heights=c(1,3))



library(simspec)
drugs_d15 <- cluster_sim_spectrum(drugs_d15, label_tag = 'inhib_annotation')
drugs_d15 <- RunUMAP(
    drugs_d15, 
    reduction='css', 
    dims=1:ncol(drugs_d15[['css']]),
    reduction.name='cssumap'
)
p1 <- dim_plot(drugs_d15, group.by=c('DMSO', 'inhib_annotation'), reduction='cssumap')
p2 <- feature_plot(drugs_d15, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='cssumap')
p_css <- p1 / p2 + plot_layout(heights=c(1,3))


library(SeuratWrappers)
library(harmony)
drugs_d15 <- RunHarmony(drugs_d15, group.by.vars = 'orig.ident')
drugs_d15 <- RunUMAP(
    drugs_d15, 
    reduction='harmony', 
    dims=1:ncol(drugs_d15[['harmony']]),
    reduction.name='humap'
)
p1 <- dim_plot(drugs_d15, group.by=c('DMSO', 'inhib_annotation'), reduction='humap')
p2 <- feature_plot(drugs_d15, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='humap')
p_harm <- p1 / p2 + plot_layout(heights=c(1,3))


p_nat | p_css | p_harm
ggsave('plots/drugs/drugs_d15_integration_screen_umap.png', width=16, height=10)

drugs_d15 %>% write_rds('data/drugs/drugs_d15_v1_v1.2lines_srt.rds')


#### d18 ####
drugs_d18 <- drugs_all %>% subset(orig.ident=='day_18')

drugs_d18 <- drugs_d18 %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

ElbowPlot(drugs_d18)

drugs_d18 <- RunUMAP(
    drugs_d18, 
    reduction='pca', 
    dims=1:20,
    reduction.name='umap'
)

p1 <- dim_plot(drugs_d18, group.by=c('DMSO', 'inhib_annotation'))
p2 <- feature_plot(drugs_d18, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T)
p_nat <- p1 / p2 + plot_layout(heights=c(1,3))



library(simspec)
drugs_d18 <- cluster_sim_spectrum(drugs_d18, label_tag = 'inhib_annotation')
drugs_d18 <- RunUMAP(
    drugs_d18, 
    reduction='css', 
    dims=1:ncol(drugs_d18[['css']]),
    reduction.name='cssumap'
)
p1 <- dim_plot(drugs_d18, group.by=c('DMSO', 'inhib_annotation'), reduction='cssumap')
p2 <- feature_plot(drugs_d18, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='cssumap')
p_css <- p1 / p2 + plot_layout(heights=c(1,3))


library(SeuratWrappers)
library(harmony)
drugs_d18 <- RunHarmony(drugs_d18, group.by.vars = 'orig.ident')
drugs_d18 <- RunUMAP(
    drugs_d18, 
    reduction='harmony', 
    dims=1:ncol(drugs_d18[['harmony']]),
    reduction.name='humap'
)
p1 <- dim_plot(drugs_d18, group.by=c('DMSO', 'inhib_annotation'), reduction='humap')
p2 <- feature_plot(drugs_d18, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='humap')
p_harm <- p1 / p2 + plot_layout(heights=c(1,3))


p_nat | p_css | p_harm
ggsave('plots/drugs/drugs_d18_integration_screen_umap.png', width=16, height=10)

drugs_d18 %>% write_rds('data/drugs/drugs_d18_v1_v1.2lines_srt.rds')


#### d15 & d18 ####
drugs_d158 <- drugs_all %>% subset(orig.ident%in%c('day_15', 'day_18'))

drugs_d158 <- drugs_d158 %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

ElbowPlot(drugs_d158)

drugs_d158 <- RunUMAP(
    drugs_d158, 
    reduction='pca', 
    dims=1:20,
    reduction.name='umap'
)

p1 <- dim_plot(drugs_d158, group.by=c('DMSO', 'inhib_annotation'))
p2 <- feature_plot(drugs_d158, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T)
p_nat <- p1 / p2 + plot_layout(heights=c(1,3))



library(simspec)
drugs_d158 <- cluster_sim_spectrum(drugs_d158, label_tag = 'inhib_annotation')
drugs_d158 <- RunUMAP(
    drugs_d158, 
    reduction='css', 
    dims=1:ncol(drugs_d158[['css']]),
    reduction.name='cssumap'
)
p1 <- dim_plot(drugs_d158, group.by=c('DMSO', 'inhib_annotation'), reduction='cssumap')
p2 <- feature_plot(drugs_d158, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='cssumap')
p_css <- p1 / p2 + plot_layout(heights=c(1,3))


library(SeuratWrappers)
library(harmony)
drugs_d158 <- RunHarmony(drugs_d158, group.by.vars = 'orig.ident')
drugs_d158 <- RunUMAP(
    drugs_d158, 
    reduction='harmony', 
    dims=1:ncol(drugs_d158[['harmony']]),
    reduction.name='humap'
)
p1 <- dim_plot(drugs_d158, group.by=c('DMSO', 'inhib_annotation'), reduction='humap')
p2 <- feature_plot(drugs_d158, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='humap')
p_harm <- p1 / p2 + plot_layout(heights=c(1,3))


p_nat | p_css | p_harm
ggsave('plots/drugs/drugs_d15_d18_integration_screen_umap.png', width=16, height=10)


drugs_d158 %>% write_rds('data/drugs/drugs_d15_d18_v1_v1.2lines_srt.rds')


####### -> Continue with d15+d18 ####
drugs_d158 <- read_rds('data/drugs/drugs_d15_d18_v1_v1.2lines_srt.rds')

drugs_d158$inhibitor_target <- case_when(
    str_detect(drugs_d158$inhib_annotation, '485') ~ 'H3K27ac',
    str_detect(drugs_d158$inhib_annotation, '395') ~ 'H3K27me3',
    T ~ 'none'
)

dim_plot(drugs_d158, group.by=c('DMSO', 'inhibitor_target'), reduction='humap')


#### Clustering and enrichment ####

drugs_d158 <- FindNeighbors(drugs_d158, reduction='harmony')
drugs_d158 <- FindClusters(drugs_d158, reduction='harmony', resolution=0.2)

dim_plot(drugs_d158, group.by=c('seurat_clusters', 'inhibitor_target', 'orig.ident'), reduction='humap', label=T)
ggsave('plots/drugs/drugs_d15_d18_clusters_umap.png', width=10, height=5)


d158_enrich_sample <- map_dfr(set_names(c('H3K27ac', 'H3K27me3')), function(x){
    map_dfr(set_names(c('day_15', 'day_18')), function(y){
        sample_idx <- (drugs_d158$inhibitor_target %in% c(x, 'none')) & (drugs_d158$orig.ident == y)
        perturb_vec <- factor(drugs_d158$inhibitor_target[sample_idx], levels=c('none',x))
        cluster_vec <- drugs_d158$seurat_clusters[sample_idx]
        group_vec <- drugs_d158$orig.ident[sample_idx]
        map_dfr(set_names(unique(cluster_vec)), function(c){
            ctest <- mantelhaen_test(perturb_vec, cluster_vec==c, z=group_vec)
            ftest <- fisher.test(perturb_vec, cluster_vec==c)
            return(tibble(
                pval_cmh = ctest$p.value,
                pval_fisher = ftest$p.value,
                common_odds_ratio = ctest$estimate,
                common_logodds = log2(ctest$estimate),
                odds_ratio = ftest$estimate,
                logodds = log2(ftest$estimate)
            ))
        }, .id='group')
    }, .id='sample')
}, .id='inhibitor')


d158_enrich_all <- map_dfr(set_names(c('H3K27ac', 'H3K27me3')), function(x){
    sample_idx <- (drugs_d158$inhibitor_target %in% c(x, 'none')) 
    perturb_vec <- factor(drugs_d158$inhibitor_target[sample_idx], levels=c('none',x))
    cluster_vec <- drugs_d158$seurat_clusters[sample_idx]
    group_vec <- drugs_d158$orig.ident[sample_idx]
    map_dfr(set_names(unique(cluster_vec)), function(c){
        ctest <- mantelhaen_test(perturb_vec, cluster_vec==c, z=group_vec)
        ftest <- fisher.test(perturb_vec, cluster_vec==c)
        return(tibble(
            pval_cmh = ctest$p.value,
            pval_fisher = ftest$p.value,
            common_odds_ratio = ctest$estimate,
            common_logodds = log2(ctest$estimate),
            odds_ratio = ftest$estimate,
            logodds = log2(ftest$estimate)
        ))
    }, .id='group')
}, .id='inhibitor')


d158_cluster_enrich <- bind_rows('sample'=d158_enrich_sample, 'all'=d158_enrich_all, .id='level')

ggplot(d158_enrich_sample, aes(group, common_logodds, fill=sample)) +
    geom_bar(stat='identity', position='dodge') +
    facet_grid(~inhibitor)
ggsave('plots/drugs/drugs_d15_d18_r02_enrich_samples_bar.png', width=5, height=3)


ggplot(d158_enrich_all, aes(group, common_logodds)) +
    geom_bar(stat='identity', position='dodge', width=0.05, size=0) +
    geom_point(aes(fill=-log10(pval_cmh), size=-log10(pval_cmh)), shape=21) +
    facet_grid(~inhibitor)
ggsave('plots/drugs/drugs_d15_d18_r02_enrich_all_bar.png', width=5, height=3)


cluster_de <- de(drugs_d158, groups='seurat_clusters')
vulcano_plot(cluster_de)
ggsave('plots/drugs/drugs_d15_d18_r02_cluster_de_vulcano.png', width=20, height=12)


plot_df <- cluster_de %>% 
    group_by(group) %>% 
    filter(padj<1e-4) %>% 
    summarize(up_ratio=log2(sum(fc>0.25) / sum(fc<(-0.25)))) %>% 
    arrange(desc(up_ratio)) %>% 
    mutate(group=factor(group, levels=unique(.$group)))
    
ggplot(plot_df, aes(up_ratio, group)) +
    geom_bar(position='dodge', stat='identity')

ggsave('plots/drugs/drugs_d15_d18_up_ratio_bar.png', width=3, height=3)


# Write test files
d158_cluster_enrich %>% write_tsv('data/drugs/results/enrichment/drugs_d15_d18_r02_enrich.tsv')
cluster_de %>% write_tsv('data/drugs/results/diff_expression/drugs_d15_d18_r02_cluster_de.tsv')


#### Feature plots of examples 
goi <- cluster_de %>% 
    filter(padj<1e-4, group%in%c(0,5)) %>% 
    group_by(group) %>% 
    top_n(10, fc) %>% 
    pull(feature) %>% unique()

feature_plot(drugs_d158, features=goi, reduction='humap', order=T)
ggsave('plots/drugs/drugs_d15_d18_c05_feature_umap.png', width=8, height=12)




goi <- cluster_de %>% 
    filter(padj<1e-4, group%in%c(0), !feature%in%cc_genes_all) %>% 
    group_by(group) %>% 
    top_n(16, -fc) %>% 
    pull(feature) %>% unique()

feature_plot(drugs_d158, features=goi, reduction='humap', order=T)
ggsave('plots/drugs/drugs_d15_d18_c0_down_feature_umap.png', width=8, height=12)



goi <- cluster_de %>% 
    filter(padj<1e-4, group%in%c(8), !feature%in%cc_genes_all) %>% 
    group_by(group) %>% 
    top_n(16, fc) %>% 
    pull(feature) %>% unique()

feature_plot(drugs_d158, features=goi, reduction='humap', order=T)
ggsave('plots/drugs/drugs_d15_d18_c8_up_feature_umap.png', width=8, height=12)



goi <- cluster_de %>% 
    filter(padj<1e-4, group%in%c(7), !feature%in%cc_genes_all) %>% 
    group_by(group) %>% 
    top_n(16, fc) %>% 
    pull(feature) %>% unique()

feature_plot(drugs_d158, features=goi, reduction='humap', order=T)
ggsave('plots/drugs/drugs_d15_d18_c7_up_feature_umap.png', width=8, height=12)



goi <- cluster_de %>% 
    filter(padj<1e-4, group%in%c(9), !feature%in%cc_genes_all) %>% 
    group_by(group) %>% 
    top_n(16, fc) %>% 
    pull(feature) %>% unique()

feature_plot(drugs_d158, features=goi, reduction='humap', order=T)
ggsave('plots/drugs/drugs_d15_d18_c9_up_feature_umap.png', width=8, height=12)








###### Later stage -> d21 ####

drugs_d21 <- drugs_all %>% subset(orig.ident!='day_15' & orig.ident!='day_18')

drugs_d21 <- drugs_d21 %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()

ElbowPlot(drugs_d21)

drugs_d21 <- RunUMAP(
    drugs_d21, 
    reduction='pca', 
    dims=1:20,
    reduction.name='umap'
)

p1 <- dim_plot(drugs_d21, group.by=c('DMSO', 'inhib_annotation'), reduction='umap')
p2 <- feature_plot(drugs_d21, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T)
p_nat <- p1 / p2 + plot_layout(heights=c(1,3))


library(simspec)
drugs_d21 <- cluster_sim_spectrum(drugs_d21, label_tag = 'inhib_annotation')
drugs_d21 <- RunUMAP(
    drugs_d21, 
    reduction='css', 
    dims=1:ncol(drugs_d21[['css']]),
    reduction.name='cssumap'
)
p1 <- dim_plot(drugs_d21, group.by=c('DMSO', 'inhib_annotation'), reduction='cssumap')
p2 <- feature_plot(drugs_d21, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='cssumap')
p_css <- p1 / p2 + plot_layout(heights=c(1,3))


library(SeuratWrappers)
library(harmony)
drugs_d21 <- RunHarmony(drugs_d21, group.by.vars = 'orig.ident')
drugs_d21 <- RunUMAP(
    drugs_d21, 
    reduction='harmony', 
    dims=1:ncol(drugs_d21[['harmony']]),
    reduction.name='humap'
)
p1 <- dim_plot(drugs_d21, group.by=c('DMSO', 'inhib_annotation'), reduction='humap')
p2 <- feature_plot(drugs_d21, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='humap')
p_harm <- p1 / p2 + plot_layout(heights=c(1,3))


p_nat | p_css | p_harm
ggsave('plots/drugs/drugs_d21_integration_screen_umap.png', width=28, height=18)


drugs_d21 %>% write_rds('data/drugs/drugs_d21_v1_v1.2lines_srt.rds')



####### Reannot and subcluster ####
drugs_d21 <- read_rds('data/drugs/drugs_d21_v1_v1.2lines_srt.rds')

drugs_d21$inhibitor_target <- case_when(
    str_detect(drugs_d21$inhib_annotation, '485') ~ 'H3K27ac',
    str_detect(drugs_d21$inhib_annotation, '395') ~ 'H3K27me3',
    T ~ 'none'
)

drugs_d21$replicate <- case_when(
    str_detect(drugs_d21$orig.ident, 'day21_2') ~ '2',
    T ~ '1'
)

dim_plot(drugs_d21, group.by=c('replicate', 'inhibitor_target'), reduction='humap')


#### Clustering and enrichment ####
drugs_d21 <- FindNeighbors(drugs_d21, reduction='harmony')
drugs_d21 <- FindClusters(drugs_d21, reduction='harmony', resolution=0.4)

dim_plot(drugs_d21, group.by=c('seurat_clusters', 'inhibitor_target', 'replicate'), reduction='humap', label=T)
ggsave('plots/drugs/drugs_d21_clusters_umap.png', width=10, height=5)


d21_enrich_sample <- map_dfr(set_names(c('H3K27ac', 'H3K27me3')), function(x){
    map_dfr(set_names(c('1', '2')), function(y){
        sample_idx <- (drugs_d21$inhibitor_target %in% c(x, 'none')) & (drugs_d21$replicate == y)
        perturb_vec <- factor(drugs_d21$inhibitor_target[sample_idx], levels=c('none',x))
        cluster_vec <- drugs_d21$seurat_clusters[sample_idx]
        group_vec <- drugs_d21$replicate[sample_idx]
        map_dfr(set_names(unique(cluster_vec)), function(c){
            ctest <- mantelhaen_test(perturb_vec, cluster_vec==c, z=group_vec)
            ftest <- fisher.test(perturb_vec, cluster_vec==c)
            return(tibble(
                pval_cmh = ctest$p.value,
                pval_fisher = ftest$p.value,
                common_odds_ratio = ctest$estimate,
                common_logodds = log2(ctest$estimate),
                odds_ratio = ftest$estimate,
                logodds = log2(ftest$estimate)
            ))
        }, .id='group')
    }, .id='sample')
}, .id='inhibitor')


d21_enrich_all <- map_dfr(set_names(c('H3K27ac', 'H3K27me3')), function(x){
    sample_idx <- (drugs_d21$inhibitor_target %in% c(x, 'none')) 
    perturb_vec <- factor(drugs_d21$inhibitor_target[sample_idx], levels=c('none',x))
    cluster_vec <- drugs_d21$seurat_clusters[sample_idx]
    group_vec <- drugs_d21$replicate[sample_idx]
    map_dfr(set_names(unique(cluster_vec)), function(c){
        ctest <- mantelhaen_test(perturb_vec, cluster_vec==c, z=group_vec)
        ftest <- fisher.test(perturb_vec, cluster_vec==c)
        return(tibble(
            pval_cmh = ctest$p.value,
            pval_fisher = ftest$p.value,
            common_odds_ratio = ctest$estimate,
            common_logodds = log2(ctest$estimate),
            odds_ratio = ftest$estimate,
            logodds = log2(ftest$estimate)
        ))
    }, .id='group')
}, .id='inhibitor')


d21_cluster_enrich <- bind_rows('sample'=d21_enrich_sample, 'all'=d21_enrich_all, .id='level')


ggplot(d21_enrich_sample, aes(group, common_logodds, fill=sample)) +
    geom_bar(stat='identity', position='dodge') +
    facet_grid(~inhibitor)
ggsave('plots/drugs/drugs_d21_r02_enrich_samples_bar.png', width=5, height=3)


ggplot(d21_enrich_all, aes(group, common_logodds)) +
    geom_bar(stat='identity', position='dodge', width=0.05, size=0) +
    geom_point(aes(fill=-log10(pval_cmh), size=-log10(pval_cmh)), shape=21) +
    facet_grid(~inhibitor)
ggsave('plots/drugs/drugs_d21_r02_enrich_all_bar.png', width=5, height=3)


cluster_de <- de(drugs_d21, groups='seurat_clusters')
vulcano_plot(cluster_de)
ggsave('plots/drugs/drugs_d21_r04_cluster_de_vulcano.png', width=20, height=12)


plot_df <- cluster_de %>% 
    group_by(group) %>% 
    filter(padj<1e-4) %>% 
    summarize(up_ratio=log2(sum(fc>0.25) / sum(fc<(-0.25)))) %>% 
    arrange(desc(up_ratio)) %>% 
    mutate(group=factor(group, levels=unique(.$group)))

ggplot(plot_df, aes(up_ratio, group)) +
    geom_bar(position='dodge', stat='identity')

ggsave('plots/drugs/drugs_d21_up_ratio_bar.png', width=3, height=3)


# Write test files
d21_cluster_enrich %>% write_tsv('data/drugs/results/enrichment/drugs_d21_r04_enrich.tsv')
cluster_de %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_r04_cluster_de.tsv')












