source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist

setwd('~/projects/cutntag/')


#### Read data ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.3lines.rds')

neuron_vars <- read_tsv('data/trajectories/ALL_lineage_pt_r2.tsv')
ac_vars <- read_tsv('data/trajectories/ALL_astrocyte_telen_diff_r2.tsv')



#### Compute max branch for all genes ####

rna_clusters_use <- rna$clusters[rna$lineage%in%c('nt', 'dien', 'ctx')] %>% unique()
H3K27ac_clusters_use <- marks$H3K27ac$clusters[marks$H3K27ac$lineage%in%c('nt', 'dien', 'ctx')] %>% unique()
H3K27me3_clusters_use <- marks$H3K27me3$clusters[marks$H3K27me3$lineage%in%c('nt', 'dien', 'ctx')] %>% unique()
H3K4me3_clusters_use <- marks$H3K4me3$clusters[marks$H3K4me3$lineage%in%c('nt', 'dien', 'ctx')] %>% unique()

rna_cluster_expr <- rna@assays$RNA@misc$summary$clusters[rna_clusters_use,]
H3K27ac_cluster_expr <- marks$H3K27ac@assays$peaks@misc$summary$clusters[H3K27ac_clusters_use,]
H3K27me3_cluster_expr <- marks$H3K27me3@assays$peaks@misc$summary$clusters[H3K27me3_clusters_use,]
H3K4me3_cluster_expr <- marks$H3K4me3@assays$peaks@misc$summary$clusters[H3K4me3_clusters_use,]

rna_cluster_max <- rownames(rna_cluster_expr)[apply(rna_cluster_expr, 2, which.max)]
names(rna_cluster_max) <- colnames(rna_cluster_expr)
H3K27ac_cluster_max <- rownames(H3K27ac_cluster_expr)[apply(H3K27ac_cluster_expr, 2, which.max)]
names(H3K27ac_cluster_max) <- colnames(H3K27ac_cluster_expr)
H3K27me3_cluster_max <- rownames(H3K27me3_cluster_expr)[apply(H3K27me3_cluster_expr, 2, which.max)]
names(H3K27me3_cluster_max) <- colnames(H3K27me3_cluster_expr)
H3K4me3_cluster_max <- rownames(H3K4me3_cluster_expr)[apply(H3K4me3_cluster_expr, 2, which.max)]
names(H3K4me3_cluster_max) <- colnames(H3K4me3_cluster_expr)

rna_cluster_lin <- rna@meta.data %>% as_tibble(rownames='cell') %>% distinct(clusters, lineage, celltype_jf)
H3K27ac_cluster_lin <- marks$H3K27ac@meta.data %>% as_tibble(rownames='cell') %>% distinct(clusters, lineage, celltype_jf)
H3K27me3_cluster_lin <- marks$H3K27me3@meta.data %>% as_tibble(rownames='cell') %>% distinct(clusters, lineage, celltype_jf)
H3K4me3_cluster_lin <- marks$H3K4me3@meta.data %>% as_tibble(rownames='cell') %>% distinct(clusters, lineage, celltype_jf)

rna_max_lin <- rna_cluster_max %>% enframe('feature', 'clusters') %>% inner_join(rna_cluster_lin)
H3K27ac_max_lin <- H3K27ac_cluster_max %>% enframe('feature', 'clusters') %>% inner_join(H3K27ac_cluster_lin)
H3K27me3_max_lin <- H3K27me3_cluster_max %>% enframe('feature', 'clusters') %>% inner_join(H3K27me3_cluster_lin)
H3K4me3_max_lin <- H3K4me3_cluster_max %>% enframe('feature', 'clusters') %>% inner_join(H3K4me3_cluster_lin)

max_lin_all <- bind_rows('RNA'=rna_max_lin, 'H3K27ac'=H3K27ac_max_lin, 'H3K27me3'=H3K27me3_max_lin, 'H3K4me3'=H3K4me3_max_lin, .id='modality')


#### Plot gene and peak variance colored by max lineage/celltype ####
lineage_cols <- c('ctx'=pantone_celltype[['ctx_npc']], 'nt'=pantone_celltype[['nt_npc']], 'dien'=pantone_celltype[['dien_npc']])

plot_df <- neuron_vars %>% 
    filter(feature_type=='peak', vars!='lin_pt', !is.na(peak)) %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=vars, values_from=r2) 

ggplot(plot_df, aes(lin, pt)) +
    geom_hex() +
    # geom_text_repel(data=filter(plot_df, top_lin | top_pt), size=2, max.overlaps=99999) +
    facet_grid(modality~.) +
    article_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Variance explained by region', y='Variance explained by pseudotime')

#### Fuck lost the code, have to redo 

plot_df <- neuron_vars %>% 
    filter(feature_type=='peak', vars!='lin_pt', !is.na(peak)) %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=vars, values_from=r2)

plots <- map(unique(plot_df$modality), function(m){
    p_df <- filter(plot_df, modality==m)
    ggplot(p_df, aes(pt, lin, label=gene)) +
        geom_hex(bins=100) +
        scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
        scale_axis_rangeframe() + theme_rangeframe() +
        article_text() +
        scale_x_continuous(limits=c(0,0.5)) +
        scale_y_continuous(limits=c(0,0.5)) +
        facet_grid(modality~.) +
        labs(x='Variance explained by pseudotime', y='Variance explained by region')
})
wrap_plots(plots, ncol=1)
ggsave('plots/paper/fig3/fig3_lin_pt_variance_density.pdf', width=6.3, height=11.5, units='cm')


ggplot(plot_df, aes(pt, lin, label=gene)) +
    geom_hex(bins=100, alpha=0.5) +
    geom_text(size=1) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(limits=c(0,0.5)) +
    scale_y_continuous(limits=c(0,0.5)) +
    facet_grid(modality~.) +
    labs(x='Variance explained by pseudotime', y='Variance explained by region')
ggsave('plots/paper/fig3/fig3_lin_pt_variance_labels.pdf', width=100, height=100, units='cm')


plot_df <- neuron_vars %>% 
    filter(feature_type=='gene', vars!='lin_pt', modality=='RNA') %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=vars, values_from=r2)

ggplot(plot_df, aes(pt, lin, label=gene)) +
    geom_hex(bins=100) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(limits=c(0,0.8)) +
    scale_y_continuous(limits=c(0,1)) +
    facet_grid(modality~.) +
    labs(x='Variance explained by pseudotime', y='Variance explained by region')
ggsave('plots/paper/fig3/fig3_lin_pt_RNA_variance_density.pdf', width=6.3, height=4, units='cm')


ggplot(plot_df, aes(pt, lin, label=gene)) +
    geom_hex(bins=100, alpha=0.5) +
    geom_text(size=1) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    # scale_x_continuous(limits=c(0,0.5)) +
    # scale_y_continuous(limits=c(0,0.5)) +
    facet_grid(modality~.) +
    labs(x='Variance explained by pseudotime', y='Variance explained by region')
ggsave('plots/paper/fig3/fig3_lin_pt_RNA_variance_labels.pdf', width=100, height=100, units='cm')






plot_df <- ac_vars %>% 
    filter(feature_type=='peak', vars!='state_region', !is.na(peak)) %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=vars, values_from=r2)

plots <- map(unique(plot_df$modality), function(m){
    p_df <- filter(plot_df, modality==m)
    ggplot(p_df, aes(state, region, label=gene)) +
        geom_hex(bins=100) +
        scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
        scale_axis_rangeframe() + theme_rangeframe() +
        article_text() +
        scale_x_continuous(limits=c(0,0.4)) +
        scale_y_continuous(limits=c(0,0.4)) +
        facet_grid(modality~.) +
        labs(x='Variance explained by state', y='Variance explained by region')
})
wrap_plots(plots, ncol=1)
ggsave('plots/paper/fig3/fig3_astrocyte_variance_density.pdf', width=6.3, height=11.5, units='cm')


ggplot(plot_df, aes(state, region, label=gene)) +
    geom_hex(bins=100, alpha=0.5) +
    geom_text(size=1) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    # scale_x_continuous(limits=c(0,0.5)) +
    # scale_y_continuous(limits=c(0,0.5)) +
    facet_grid(modality~.) +
    labs(x='Variance explained by pseudotime', y='Variance explained by region')
ggsave('plots/paper/fig3/fig3_astrocyte_variance_labels.pdf', width=100, height=100, units='cm')



plot_df <- ac_vars %>% 
    filter(feature_type=='gene', vars!='state_region', modality=='RNA') %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=vars, values_from=r2)

ggplot(plot_df, aes(state, region, label=gene)) +
    geom_hex(bins=100) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(limits=c(0,0.8)) +
    scale_y_continuous(limits=c(0,1)) +
    labs(x='Variance explained by state', y='Variance explained by region')
ggsave('plots/paper/fig3/fig3_astrocyte_RNA_variance_density.pdf', width=6.3, height=4, units='cm')


ggplot(plot_df, aes(state, region, label=gene)) +
    geom_hex(bins=100, alpha=0.5) +
    geom_text(size=1) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    # scale_x_continuous(limits=c(0,0.5)) +
    # scale_y_continuous(limits=c(0,0.5)) +
    facet_grid(modality~.) +
    labs(x='Variance explained by pseudotime', y='Variance explained by region')
ggsave('plots/paper/fig3/fig3_astrocyte_RNA_variance_labels.pdf', width=100, height=100, units='cm')



#### Boxplots to compare modalities ####
plot_df <- neuron_vars %>% 
    filter(feature_type=='peak', !is.na(peak)) %>% 
    mutate(
        modality=factor(modality, levels=c('H3K27ac', 'H3K27me3', 'H3K4me3')),
        vars=factor(vars, levels=c('lin', 'pt', 'lin_pt'))
    )

ggplot(plot_df, aes(vars, r2, fill=modality)) +
    geom_boxplot(size=0.1, outlier.size = 0.1, outlier.shape = 16, outlier.alpha = 0.5) +
    scale_fill_manual(values=modality_colors) +
    scale_x_discrete(labels=c('Region', 'Differentiation', 'Region+Diff.')) +
    article_text() +
    rotate_x_text(45) +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Variables', y='R2', fill='Modality')

ggsave('plots/paper/fig3/fig3_lin_pt_variance_boxplot.pdf', width=5, height=6, units='cm')

plot_df <- neuron_vars %>% 
    filter(feature_type=='gene', modality=='RNA') %>% 
    mutate(
        vars=factor(vars, levels=c('lin', 'pt', 'lin_pt'))
    )

ggplot(plot_df, aes(vars, r2, fill=modality)) +
    geom_boxplot(size=0.1, outlier.size = 0.1, outlier.shape = 16, outlier.alpha = 0.5) +
    scale_fill_manual(values=modality_colors) +
    scale_x_discrete(labels=c('Region', 'Differentiation', 'Region+Diff.')) +
    article_text() +
    rotate_x_text(45) +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Variables', y='R2', fill='Modality')

ggsave('plots/paper/fig3/fig3_lin_pt_RNA_variance_boxplot.pdf', width=4, height=6, units='cm')



plot_df <- ac_vars %>% 
    filter(feature_type=='peak', !is.na(peak)) %>% 
    mutate(
        modality=factor(modality, levels=c('H3K27ac', 'H3K27me3', 'H3K4me3')),
        vars=factor(vars, levels=c('region', 'state', 'state_region'))
    )

ggplot(plot_df, aes(vars, r2, fill=modality)) +
    geom_boxplot(size=0.1, outlier.size = 0.1, outlier.shape = 16, outlier.alpha = 0.5) +
    scale_fill_manual(values=modality_colors) +
    scale_x_discrete(labels=c('Region', 'Differentiation', 'Region+Diff.')) +
    article_text() +
    rotate_x_text(45) +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Variables', y='R2', fill='Modality')

ggsave('plots/paper/fig3/fig3_astrocyte_variance_boxplot.pdf', width=5, height=6, units='cm')


plot_df <- ac_vars %>% 
    filter(feature_type=='gene', modality=='RNA') %>% 
    mutate(
        vars=factor(vars, levels=c('region', 'state', 'state_region'))
    )

ggplot(plot_df, aes(vars, r2, fill=modality)) +
    geom_boxplot(size=0.1, outlier.size = 0.1, outlier.shape = 16, outlier.alpha = 0.5) +
    scale_fill_manual(values=modality_colors) +
    scale_x_discrete(labels=c('Region', 'Differentiation', 'Region+Diff.')) +
    article_text() +
    rotate_x_text(45) +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Variables', y='R2', fill='Modality')

ggsave('plots/paper/fig3/fig3_astrocyte_RNA_variance_boxplot.pdf', width=4, height=6, units='cm')





#### Quantify genes with var > 0.1 ####
library(ggupset)

rna_neuro %>% write_rds('data/trajectories/RNA_neuro_dpt_srt.rds')
rna_genes <- rna_neuro %>% FindVariableFeatures(nfeatures=10000) %>% VariableFeatures()
       
plot_df <- neuron_vars %>%
    filter(feature_type=='peak' | modality=='RNA', vars=='lin', gene%in%rna_genes) %>%
    group_by(gene, vars, modality) %>% filter(r2==max(r2)) %>% 
    group_by(modality, vars) %>% mutate(q10=r2 > quantile(r2, probs=0.1)) %>% 
    filter(q10) %>%
    # filter(r2>0.1) %>%
    distinct(gene, modality, vars) %>%
    group_by(gene) %>% 
    # summarize(modality_list=list(modality), modality_str=paste0(sort(modality), collapse='_'), 
              # nmod=length(modality), has_rna='RNA'%in%modality) 
    summarize(
        RNA = 'RNA'%in%modality,
        H3K27ac = 'H3K27ac'%in%modality,
        H3K27me3 = 'H3K27me3'%in%modality,
        H3K4me3 = 'H3K4me3'%in%modality,
        RNA_H3K27ac = 'H3K27ac'%in%modality & 'RNA'%in%modality,
        RNA_H3K27me3 = 'H3K27me3'%in%modality & 'RNA'%in%modality,
        RNA_H3K4me3 = 'H3K4me3'%in%modality & 'RNA'%in%modality,
        RNA_exclusive = 'RNA'%in%modality & !('H3K27me3'%in%modality | 'H3K27ac'%in%modality | 'H3K4me3'%in%modality)
    ) %>% 
    summarize(
        RNA = sum(RNA),
        # H3K27ac = sum(H3K27ac),
        # H3K27me3 = sum(H3K27me3),
        # H3K4me3 = sum(H3K4me3),
        RNA_H3K27ac = sum(RNA_H3K27ac),
        RNA_H3K27me3 = sum(RNA_H3K27me3),
        RNA_H3K4me3 = sum(RNA_H3K4me3),
        RNA_exclusive = sum(RNA_exclusive)
    ) %>% pivot_longer(everything(), values_to='count', names_to='mods') %>% 
    mutate(mods=factor(mods, levels=c('RNA', 'RNA_H3K27ac', 'RNA_H3K27me3', 'RNA_H3K4me3', 'RNA_exclusive')))


p1 <- ggplot(plot_df, aes(mods, count)) +
    geom_bar(size=0.1, fill='darkgrey', color='black', stat='identity') +
    article_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Modality', y='# genes', title='Region\n-specific peaks')
p1

plot_df <- neuron_vars %>%
    filter(feature_type=='peak' | modality=='RNA', vars=='pt') %>%
    group_by(gene, vars, modality) %>% filter(r2==max(r2)) %>% 
    group_by(modality, vars) %>% mutate(q10=r2 > quantile(r2, probs=0.1)) %>% 
    filter(q10) %>%
    # filter(r2>0.1) %>%
    distinct(gene, modality, vars) %>%
    group_by(gene) %>% 
    # summarize(modality_list=list(modality), modality_str=paste0(sort(modality), collapse='_'), 
    #           nmod=length(modality), has_rna='RNA'%in%modality) 
    summarize(
        RNA = 'RNA'%in%modality,
        H3K27ac = 'H3K27ac'%in%modality,
        H3K27me3 = 'H3K27me3'%in%modality,
        H3K4me3 = 'H3K4me3'%in%modality,
        RNA_H3K27ac = 'H3K27ac'%in%modality & 'RNA'%in%modality,
        RNA_H3K27me3 = 'H3K27me3'%in%modality & 'RNA'%in%modality,
        RNA_H3K4me3 = 'H3K4me3'%in%modality & 'RNA'%in%modality,
        RNA_exclusive = 'RNA'%in%modality & !('H3K27me3'%in%modality | 'H3K27ac'%in%modality | 'H3K4me3'%in%modality)
    ) %>% 
    summarize(
        RNA = sum(RNA),
        # H3K27ac = sum(H3K27ac),
        # H3K27me3 = sum(H3K27me3),
        # H3K4me3 = sum(H3K4me3),
        RNA_H3K27ac = sum(RNA_H3K27ac),
        RNA_H3K27me3 = sum(RNA_H3K27me3),
        RNA_H3K4me3 = sum(RNA_H3K4me3),
        RNA_exclusive = sum(RNA_exclusive)
    ) %>% pivot_longer(everything(), values_to='count', names_to='mods') %>% 
    mutate(mods=factor(mods, levels=c('RNA', 'RNA_H3K27ac', 'RNA_H3K27me3', 'RNA_H3K4me3', 'RNA_exclusive')))


p2 <- ggplot(plot_df, aes(mods, count)) +
    geom_bar(size=0.1, fill='darkgrey', color='black', stat='identity') +
    article_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Modality', y='# genes', title='Differentiation\n-specific peaks')

(p1 | p2) & rotate_x_text(40)
ggsave('plots/paper/fig3/fig3_lin_pt_q10_specific_genes_bar.pdf', width=6, height=5, units='cm')    
   


  

#### Same for astrocytes ####

rna_ac_npc <- read_rds('data/trajectories/astrocytes/RNA_astrocytes_npcs_srt.rds')
rna_genes <- rna_ac_npc %>% FindVariableFeatures(nfeatures=10000) %>% VariableFeatures()

plot_df <- ac_vars %>%
    filter(feature_type=='peak' | modality=='RNA', vars=='region', gene%in%rna_genes) %>%
    group_by(gene, vars, modality) %>% filter(r2==max(r2)) %>% 
    group_by(modality, vars) %>% mutate(q10=r2 > quantile(r2, probs=0.1)) %>% 
    filter(q10) %>%
    # filter(r2>0.1) %>%
    distinct(gene, modality, vars) %>%
    group_by(gene) %>% 
    # summarize(modality_list=list(modality), modality_str=paste0(sort(modality), collapse='_'), 
    # nmod=length(modality), has_rna='RNA'%in%modality) 
    summarize(
        RNA = 'RNA'%in%modality,
        H3K27ac = 'H3K27ac'%in%modality,
        H3K27me3 = 'H3K27me3'%in%modality,
        H3K4me3 = 'H3K4me3'%in%modality,
        RNA_H3K27ac = 'H3K27ac'%in%modality & 'RNA'%in%modality,
        RNA_H3K27me3 = 'H3K27me3'%in%modality & 'RNA'%in%modality,
        RNA_H3K4me3 = 'H3K4me3'%in%modality & 'RNA'%in%modality,
        RNA_exclusive = 'RNA'%in%modality & !('H3K27me3'%in%modality | 'H3K27ac'%in%modality | 'H3K4me3'%in%modality)
    ) %>% 
    summarize(
        RNA = sum(RNA),
        # H3K27ac = sum(H3K27ac),
        # H3K27me3 = sum(H3K27me3),
        # H3K4me3 = sum(H3K4me3),
        RNA_H3K27ac = sum(RNA_H3K27ac),
        RNA_H3K27me3 = sum(RNA_H3K27me3),
        RNA_H3K4me3 = sum(RNA_H3K4me3),
        RNA_exclusive = sum(RNA_exclusive)
    ) %>% pivot_longer(everything(), values_to='count', names_to='mods') %>% 
    mutate(mods=factor(mods, levels=c('RNA', 'RNA_H3K27ac', 'RNA_H3K27me3', 'RNA_H3K4me3', 'RNA_exclusive')))


p1 <- ggplot(plot_df, aes(mods, count)) +
    geom_bar(size=0.1, fill='darkgrey', color='black', stat='identity') +
    article_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Modality', y='# genes', title='Region\n-specific peaks')
p1

plot_df <- ac_vars %>%
    filter(feature_type=='peak' | modality=='RNA', vars=='state') %>%
    group_by(gene, vars, modality) %>% filter(r2==max(r2)) %>% 
    group_by(modality, vars) %>% mutate(q10=r2 > quantile(r2, probs=0.1)) %>% 
    filter(q10) %>%
    # filter(r2>0.1) %>%
    distinct(gene, modality, vars) %>%
    group_by(gene) %>% 
    # summarize(modality_list=list(modality), modality_str=paste0(sort(modality), collapse='_'), 
    #           nmod=length(modality), has_rna='RNA'%in%modality) 
    summarize(
        RNA = 'RNA'%in%modality,
        H3K27ac = 'H3K27ac'%in%modality,
        H3K27me3 = 'H3K27me3'%in%modality,
        H3K4me3 = 'H3K4me3'%in%modality,
        RNA_H3K27ac = 'H3K27ac'%in%modality & 'RNA'%in%modality,
        RNA_H3K27me3 = 'H3K27me3'%in%modality & 'RNA'%in%modality,
        RNA_H3K4me3 = 'H3K4me3'%in%modality & 'RNA'%in%modality,
        RNA_exclusive = 'RNA'%in%modality & !('H3K27me3'%in%modality | 'H3K27ac'%in%modality | 'H3K4me3'%in%modality)
    ) %>% 
    summarize(
        RNA = sum(RNA),
        # H3K27ac = sum(H3K27ac),
        # H3K27me3 = sum(H3K27me3),
        # H3K4me3 = sum(H3K4me3),
        RNA_H3K27ac = sum(RNA_H3K27ac),
        RNA_H3K27me3 = sum(RNA_H3K27me3),
        RNA_H3K4me3 = sum(RNA_H3K4me3),
        RNA_exclusive = sum(RNA_exclusive)
    ) %>% pivot_longer(everything(), values_to='count', names_to='mods') %>% 
    mutate(mods=factor(mods, levels=c('RNA', 'RNA_H3K27ac', 'RNA_H3K27me3', 'RNA_H3K4me3', 'RNA_exclusive')))


p2 <- ggplot(plot_df, aes(mods, count)) +
    geom_bar(size=0.1, fill='darkgrey', color='black', stat='identity') +
    article_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Modality', y='# genes', title='Differentiation\n-specific peaks')

(p1 | p2) & rotate_x_text(40)
ggsave('plots/paper/fig3/fig3_astrocyte_q10_specific_genes_bar.pdf', width=6, height=5, units='cm')    




#### Compare NPC and Astrocyte regional variance ####

npc_vars <- read_tsv('data/trajectories/ALL_NPC_diff_r2.tsv')
aconly_vars <- read_tsv('data/trajectories/ALL_AConly_diff_r2.tsv')
neurons_vars <- read_tsv('data/trajectories/ALL_Neurons_diff_r2.tsv')

#### NPC vs astrocyte ####
compare_vars <- full_join(npc_vars, aconly_vars, by=c('gene', 'peak', 'vars', 'modality', 'feature_type'))
bind_vars <- bind_rows('NPC'=npc_vars, 'AC'=aconly_vars, .id='state')

plot_df <- compare_vars %>% 
    filter(feature_type=='peak', !is.na(peak)) %>% 
    mutate(
        modality=factor(modality, levels=c('H3K27ac', 'H3K27me3', 'H3K4me3'))
    )

ggplot(plot_df, aes(r2.x, r2.y, label=gene)) +
    geom_text(size=2) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(na.value=0) +
    scale_y_continuous(na.value=0) +
    facet_grid(modality~.) +
    labs(x='Regional variance in NPCs', y='Regional variance in Astrocytes')
ggsave('plots/paper/fig3/fig3_NPCvsAC_chrom_variance_labels.pdf', width=100, height=100, units='cm')

ggplot(plot_df, aes(r2.x, r2.y, label=gene)) +
    geom_point(size=0.1) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(na.value=0) +
    scale_y_continuous(na.value=0) +
    facet_grid(modality~.) +
    labs(x='Regional variance in NPCs', y='Regional variance in Astrocytes')
ggsave('plots/paper/fig3/fig3_NPCvsAC_chrom_variance_scatter.pdf', width=6.3, height=12, units='cm')



plot_df <- compare_vars %>% 
    filter(feature_type=='gene', modality=='RNA')

ggplot(plot_df, aes(r2.x, r2.y, label=gene)) +
    geom_text(size=2) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(na.value=0) +
    scale_y_continuous(na.value=0) +
    labs(x='Regional variance in NPCs', y='Regional variance in Astrocytes')
ggsave('plots/paper/fig3/fig3_NPCvsAC_RNA_variance_labels.pdf', width=100, height=100, units='cm')

ggplot(plot_df, aes(r2.x, r2.y, label=gene)) +
    geom_point(size=0.1) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(na.value=0) +
    scale_y_continuous(na.value=0) +
    labs(x='Regional variance in NPCs', y='Regional variance in Astrocytes')
ggsave('plots/paper/fig3/fig3_NPCvsAC_RNA_variance_scatter.pdf', width=6.3, height=4, units='cm')




#### NPC vs neurons ####
compare_vars <- full_join(npc_vars, neurons_vars, by=c('gene', 'peak', 'vars', 'modality', 'feature_type'))
bind_vars <- bind_rows('NPC'=npc_vars, 'N'=neurons_vars, .id='state')

plot_df <- compare_vars %>% 
    filter(feature_type=='peak', !is.na(peak)) %>% 
    mutate(
        modality=factor(modality, levels=c('H3K27ac', 'H3K27me3', 'H3K4me3'))
    )

ggplot(plot_df, aes(r2.x, r2.y, label=gene)) +
    geom_text(size=2) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(na.value=0) +
    scale_y_continuous(na.value=0) +
    facet_grid(modality~.) +
    labs(x='Regional variance in NPCs', y='Regional variance in Neurons')
ggsave('plots/paper/fig3/fig3_NPCvsN_chrom_variance_labels.pdf', width=100, height=100, units='cm')

ggplot(plot_df, aes(r2.x, r2.y, label=gene)) +
    geom_point(size=0.1) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(na.value=0) +
    scale_y_continuous(na.value=0) +
    facet_grid(modality~.) +
    labs(x='Regional variance in NPCs', y='Regional variance in Neurons')
ggsave('plots/paper/fig3/fig3_NPCvsN_chrom_variance_scatter.pdf', width=6.3, height=12, units='cm')





plot_df <- compare_vars %>% 
    filter(feature_type=='gene', modality=='RNA')

ggplot(plot_df, aes(r2.x, r2.y, label=gene)) +
    geom_text(size=2) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(na.value=0) +
    scale_y_continuous(na.value=0) +
    labs(x='Regional variance in NPCs', y='Regional variance in Neurons')
ggsave('plots/paper/fig3/fig3_NPCvsN_RNA_variance_labels.pdf', width=100, height=100, units='cm')

ggplot(plot_df, aes(r2.x, r2.y, label=gene)) +
    geom_point(size=0.1) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log2') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    scale_x_continuous(na.value=0) +
    scale_y_continuous(na.value=0) +
    labs(x='Regional variance in NPCs', y='Regional variance in Neurons')
ggsave('plots/paper/fig3/fig3_NPCvsN_RNA_variance_scatter.pdf', width=6.3, height=4, units='cm')




