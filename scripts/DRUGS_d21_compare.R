source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')

library(Pando)
library(harmony)

filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist

setwd('~/projects/cutntag/')

#### Read data ####
drugs_d21 <- read_rds('data/drugs/drugs_d21_v1_v1.2lines_srt.rds')

drugs_d21$inhibitor_target <- case_when(
    str_detect(drugs_d21$inhib_annotation, '485') ~ 'H3K27ac',
    str_detect(drugs_d21$inhib_annotation, '395') ~ 'H3K27me3',
    T ~ 'DMSO'
)

dim_plot(drugs_d21)
feature_plot(drugs_d21, features=c('nCount_RNA'))



#### Subset and re-preproc H3K27me3 inhib ####
drugs_d21_H3K27me3 <- drugs_d21 %>% subset(inhibitor_target%in%c('H3K27me3', 'DMSO'))

drugs_d21_H3K27me3 <- RunHarmony(drugs_d21_H3K27me3, group.by.vars = 'inhibitor_target')
drugs_d21_H3K27me3 <- RunUMAP(
    drugs_d21_H3K27me3, 
    reduction='harmony', 
    dims=1:ncol(drugs_d21_H3K27me3[['harmony']]),
    reduction.name='humap'
)
p1 <- dim_plot(drugs_d21_H3K27me3, group.by=c('DMSO', 'inhib_annotation'), reduction='humap')
p2 <- feature_plot(drugs_d21_H3K27me3, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'APOE', 'SIX3'), order=T, reduction='humap')
p1 / p2 + plot_layout(heights=c(1,3))

drugs_d21_H3K27me3 <- drugs_d21_H3K27me3 %>% 
    FindNeighbors(reduction='harmony', dims=1:ncol(drugs_d21_H3K27me3[['harmony']]))

drugs_d21_H3K27me3 <- drugs_d21_H3K27me3 %>% FindClusters(resolution=0.1)

drugs_d21_H3K27me3 %>% write_rds('data/drugs/drugs_d21_A395_v1_v1.2lines_srt.rds')


#### DE between clusters ####
dim_plot(drugs_d21_H3K27me3, reduction='humap', label=T)

cluster_r04_de <- de(drugs_d21_H3K27me3, groups = 'RNA_snn_res.0.1')

vulcano_plot(cluster_r04_de)

feature_plot(drugs_d21_H3K27me3, features=c('FOXG1', 'WLS', 'POU5F1', 'STMN2', 'SOX10', 'LMO4', 'HOXB9', 'DCN'), order=T, reduction='humap')



#### Subset and re-preproc H3K27ac inhib ####
drugs_d21_H3K27ac <- drugs_d21 %>% subset(inhibitor_target%in%c('H3K27ac', 'DMSO'))

drugs_d21_H3K27ac <- RunHarmony(drugs_d21_H3K27ac, group.by.vars = 'inhibitor_target')
drugs_d21_H3K27ac <- RunUMAP(
    drugs_d21_H3K27ac, 
    reduction='harmony', 
    dims=1:ncol(drugs_d21_H3K27ac[['harmony']]),
    reduction.name='humap'
)
p1 <- dim_plot(drugs_d21_H3K27ac, group.by=c('DMSO', 'inhib_annotation'), reduction='humap')
p2 <- feature_plot(drugs_d21_H3K27ac, features=c('FOXG1', 'WLS', 'POU5F1', 'LIN28A', 'STMN2', 'SIX3'), order=T, reduction='humap')
p1 / p2 + plot_layout(heights=c(1,3))

drugs_d21_H3K27ac <- drugs_d21_H3K27ac %>% 
    FindNeighbors(reduction='harmony', dims=1:ncol(drugs_d21_H3K27ac[['harmony']]))

drugs_d21_H3K27ac <- drugs_d21_H3K27ac %>% FindClusters(resolution=0.1)

drugs_d21_H3K27ac %>% write_rds('data/drugs/drugs_d21_A485_v1_v1.2lines_srt.rds')


#### DE between clusters ####
dim_plot(drugs_d21_H3K27ac, reduction='humap', label=T, group.by='RNA_snn_res.0.1')

cluster_r04_de <- de(drugs_d21_H3K27ac, groups = 'RNA_snn_res.0.1')

vulcano_plot(cluster_r04_de)

feature_plot(drugs_d21_H3K27ac, features=c('FOXG1', 'WLS', 'POU5F1', 'STMN2', 'SOX10', 'LMO4', 'HOXB9', 'DCN', 'DLX6'), order=T, reduction='humap')

feature_plot(drugs_d21_H3K27ac, features=c('FOXG1', 'WLS', 'nFeature_RNA'), order=T, reduction='humap')



#### Test diff enrichment #####
inhib_samples <- unique(drugs_d21_H3K27ac$inhib_annotation[!str_detect(drugs_d21_H3K27ac$inhib_annotation, 'DMSO')])
dmso_samples <- unique(drugs_d21_H3K27ac$inhib_annotation[str_detect(drugs_d21_H3K27ac$inhib_annotation, 'DMSO')])

d21_H3K27ac_enrich <- map_dfr(set_names(inhib_samples), function(x){
    sample_idx <- drugs_d21_H3K27ac$inhib_annotation %in% c(x, dmso_samples)
    cluster_vec <- drugs_d21_H3K27ac$RNA_snn_res.0.1[sample_idx]
    group_vec <- !drugs_d21_H3K27ac$inhib_annotation[sample_idx] %in% dmso_samples
    map_dfr(set_names(unique(cluster_vec)), function(c){
        ftest <- fisher.test(group_vec, cluster_vec==c)
        return(tibble(
            pval = ftest$p.value,
            odds_ratio = ftest$estimate,
            logodds = log2(ftest$estimate)
        ))
    }, .id='group')
}, .id='conc')


p1 <- ggplot(d21_H3K27ac_enrich, aes(group, logodds, fill=conc, alpha=pval<1e-4)) +
    geom_bar(stat='identity', position='dodge')

p2 <- dim_plot(drugs_d21_H3K27ac, reduction='humap', label=T, group.by='RNA_snn_res.0.1')

p1 | p2




inhib_samples <- unique(drugs_d21_H3K27me3$inhib_annotation[!str_detect(drugs_d21_H3K27me3$inhib_annotation, 'DMSO')])
dmso_samples <- unique(drugs_d21_H3K27me3$inhib_annotation[str_detect(drugs_d21_H3K27me3$inhib_annotation, 'DMSO')])

d21_H3K27me3_enrich <- map_dfr(set_names(inhib_samples), function(x){
    sample_idx <- drugs_d21_H3K27me3$inhib_annotation %in% c(x, dmso_samples)
    cluster_vec <- drugs_d21_H3K27me3$RNA_snn_res.0.1[sample_idx]
    group_vec <- !drugs_d21_H3K27me3$inhib_annotation[sample_idx] %in% dmso_samples
    map_dfr(set_names(unique(cluster_vec)), function(c){
        ftest <- fisher.test(group_vec, cluster_vec==c)
        return(tibble(
            pval = ftest$p.value,
            odds_ratio = ftest$estimate,
            logodds = log2(ftest$estimate)
        ))
    }, .id='group')
}, .id='conc')


p1 <- ggplot(d21_H3K27me3_enrich, aes(group, logodds, fill=conc, alpha=pval<1e-4)) +
    geom_bar(stat='identity', position='dodge')

p2 <- dim_plot(drugs_d21_H3K27me3, reduction='humap', label=T, group.by='RNA_snn_res.0.1')
p1 | p2



#### Compare DE with d15 ####

d21_A395_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d21_A395_global_inhib_de.tsv')
d15_A395_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_global_inhib_de.tsv')

d21_A395_global_de %>% filter(p_val_adj<0.01)

comp_de <- d21_A395_global_de %>% 
    full_join(d15_A395_global_de, by=c('feature'))

ggplot(comp_de, aes(avg_log2FC.x, avg_log2FC.y, label=feature)) +
    geom_point(alpha=0.5) +
    # geom_text() +
    scale_size_continuous(range=c(0.1,3)) +
    scale_y_continuous(na.value=0) +
    scale_x_continuous(na.value=0) +
    labs(x='log2FC(day 21)', y='log2FC(day15)')


d21_A485_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d21_A485_global_inhib_de.tsv')
d15_A485_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_global_inhib_de.tsv')

comp_de <- d21_A485_global_de %>% 
    inner_join(d15_A485_global_de, by=c('feature'))

ggplot(comp_de, aes(avg_log2FC.x, avg_log2FC.y, size=-log10(p_val.x))) +
    geom_point(alpha=0.5) +
    scale_size_continuous(range=c(0.1,3)) +
    labs(x='log2FC(day 21)', y='log2FC(day15)')


