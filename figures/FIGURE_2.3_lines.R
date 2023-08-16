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

group_order <- c('psc', 'non_nect', 'nect', 'ctx', 'dien', 'nt', 'astrocytes', 'retina', 'choroid_plexus')

#### Heatmaps split by line ####
H3K27ac_use <- marks$H3K27ac 
H3K27ac_use <- RunTFIDF(H3K27ac_use, assay='peaks')
H3K27ac_use$lineage_coarse <- case_when(
    H3K27ac_use$lineage != 'other' ~ H3K27ac_use$lineage,
    T ~ H3K27ac_use$celltype_jf
)
H3K27ac_use <- Pando::aggregate_assay(H3K27ac_use, group_name='lineage_coarse', assay='peaks')

line_cells <- H3K27ac_use@meta.data %>% filter(!is.na(line)) %>% rownames()
H3K27ac_use_ <- subset(H3K27ac_use, cells=line_cells)
H3K27ac_split <- SplitObject(H3K27ac_use_, split.by='line')
H3K27ac_split <- map(H3K27ac_split, Pando::aggregate_assay, group_name='lineage_coarse', assay='peaks')

H3K27me3_use <- marks$H3K27me3 
H3K27me3_use <- RunTFIDF(H3K27me3_use, assay='peaks')
H3K27me3_use$lineage_coarse <- case_when(
    H3K27me3_use$lineage != 'other' ~ H3K27me3_use$lineage,
    T ~ H3K27me3_use$celltype_jf
)
H3K27me3_use <- Pando::aggregate_assay(H3K27me3_use, group_name='lineage_coarse', assay='peaks')

line_cells <- H3K27me3_use@meta.data %>% filter(!is.na(line)) %>% rownames()
H3K27me3_use_ <- subset(H3K27me3_use, cells=line_cells)
H3K27me3_split <- SplitObject(H3K27me3_use_, split.by='line')
H3K27me3_split <- map(H3K27me3_split, Pando::aggregate_assay, group_name='lineage_coarse', assay='peaks')

H3K4me3_use <- marks$H3K4me3 
H3K4me3_use <- RunTFIDF(H3K4me3_use, assay='peaks')
H3K4me3_use$lineage_coarse <- case_when(
    H3K4me3_use$lineage != 'other' ~ H3K4me3_use$lineage,
    T ~ H3K4me3_use$celltype_jf
)
H3K4me3_use <- Pando::aggregate_assay(H3K4me3_use, group_name='lineage_coarse', assay='peaks')

line_cells <- H3K4me3_use@meta.data %>% filter(!is.na(line)) %>% rownames()
H3K4me3_use_ <- subset(H3K4me3_use, cells=line_cells)
H3K4me3_split <- SplitObject(H3K4me3_use_, split.by='line')
H3K4me3_split <- map(H3K4me3_split, Pando::aggregate_assay, group_name='lineage_coarse', assay='peaks')


#### Read all DA peaks ####
all_top_peaks <- read_tsv('data_/results/diff_expression/all_marks_top_DA_lineage_coarse.tsv')
all_peaks <- read_tsv('data_/results/diff_expression/all_marks_DA_lineage_coarse.tsv')
H3K27ac_top <- all_top_peaks %>% filter(mark=='H3K27ac')
H3K27me3_top <- all_top_peaks %>% filter(mark=='H3K27me3')
H3K4me3_top <- all_top_peaks %>% filter(mark=='H3K4me3')



#### Peak expression heatmap with all marks ####
H3K27ac_cluster <- H3K27ac_use@assays$peaks@misc$summary$lineage_coarse[, unique(H3K27ac_top$feature)]
H3K27ac_cluster_list <- map(H3K27ac_split, function(srt){
    srt@assays$peaks@misc$summary$lineage_coarse[, unique(H3K27ac_top$feature)]  
})

H3K27ac_peak_order <- H3K27ac_cluster %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

H3K27ac_cluster_expr <- map_dfr(H3K27ac_cluster_list, function(x){
    x %>% as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr') %>% return()
}, .id='line')

H3K27ac_plot_df <- H3K27ac_cluster_expr %>% 
    filter(feature%in%H3K27ac_top$feature, group%in%group_order) %>% 
    inner_join(select(H3K27ac_top, feature, 'top_group'=group)) %>% 
    mutate(
        feature=factor(feature, levels=H3K27ac_peak_order),
        group=factor(group, levels=group_order),
        top_group=factor(top_group, levels=group_order),
        expr_clip=ifelse(expr>1, 1, expr)
    ) 



H3K27me3_cluster <- H3K27me3_use@assays$peaks@misc$summary$lineage_coarse[, unique(H3K27me3_top$feature)]
H3K27me3_cluster_list <- map(H3K27me3_split, function(srt){
    srt@assays$peaks@misc$summary$lineage_coarse[, unique(H3K27me3_top$feature)]  
})

H3K27me3_peak_order <- H3K27me3_cluster %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

H3K27me3_cluster_expr <- map_dfr(H3K27me3_cluster_list, function(x){
    x %>% as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr') %>% return()
}, .id='line')

H3K27me3_plot_df <- H3K27me3_cluster_expr %>% 
    filter(feature%in%H3K27me3_top$feature, group%in%group_order) %>% 
    inner_join(select(H3K27me3_top, feature, 'top_group'=group)) %>% 
    mutate(
        feature=factor(feature, levels=H3K27me3_peak_order),
        group=factor(group, levels=group_order),
        top_group=factor(top_group, levels=group_order),
        expr_clip=ifelse(expr>1, 1, expr)
    ) 




H3K4me3_cluster <- H3K4me3_use@assays$peaks@misc$summary$lineage_coarse[, unique(H3K4me3_top$feature)]
H3K4me3_cluster_list <- map(H3K4me3_split, function(srt){
    srt@assays$peaks@misc$summary$lineage_coarse[, unique(H3K4me3_top$feature)]  
})

H3K4me3_peak_order <- H3K4me3_cluster %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

H3K4me3_cluster_expr <- map_dfr(H3K4me3_cluster_list, function(x){
    x %>% as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr') %>% return()
}, .id='line')

H3K4me3_plot_df <- H3K4me3_cluster_expr %>% 
    filter(feature%in%H3K4me3_top$feature, group%in%group_order) %>% 
    inner_join(select(H3K4me3_top, feature, 'top_group'=group)) %>% 
    mutate(
        feature=factor(feature, levels=H3K4me3_peak_order),
        group=factor(group, levels=group_order),
        top_group=factor(top_group, levels=group_order),
        expr_clip=ifelse(expr>1, 1, expr)
    ) %>% filter(!(group=='retina' & line!='B7'))


p1bu <- ggplot(H3K27ac_plot_df, aes(line, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=bugn(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(top_group~group, space='free', scales='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K27ac') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() 

p1 <- ggplot(H3K27ac_plot_df, aes(line, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=ylgn(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(top_group~group, space='free', scales='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K27ac') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() 


p2 <- ggplot(H3K27me3_plot_df, aes(line, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=blues(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(top_group~group, space='free', scales='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K27me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() 


p3 <- ggplot(H3K4me3_plot_df, aes(line, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=rdpu(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(top_group~group, space='free', scales='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K4me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text()


(p1 | p2 | p3) + plot_layout(guides = 'collect') & 
    theme(
        panel.spacing = unit(0.1, 'cm'),
        plot.margin = unit(rep(0.1,4), 'cm'),
        strip.text = element_blank(),
        panel.spacing.x = unit(0, 'lines')
    )

ggsave('plots/paper/fig2/rev_fig2_top50_DA_peak_line_heatmap.pdf', width=12, height=30)
ggsave('plots/paper/fig2/rev_fig2_top50_DA_peak_line_heatmap.png', width=12, height=30)

(p1bu | p2 | p3) + plot_layout(guides = 'collect') & 
    theme(
        panel.spacing = unit(0.1, 'cm'),
        plot.margin = unit(rep(0.1,4), 'cm'),
        strip.text = element_blank(),
        panel.spacing.x = unit(0, 'lines')
    )

ggsave('plots/paper/fig2/rev_fig2_top50_DA_peak_line_bugn_heatmap.pdf', width=12, height=30)
ggsave('plots/paper/fig2/rev_fig2_top50_DA_peak_line_bugn_heatmap.png', width=12, height=30)




top50_peaks <- bind_rows('H3K27ac'=H3K27ac_plot_df, 'H3K27me3'=H3K27me3_plot_df, 'H3K4me3'=H3K4me3_plot_df, .id='mark')
top50_peaks %>% write_rds('data_/results/diff_expression/all_marks_DA_top50_peaks.tsv')


#### Heatmap with matched genes of all marks ####

# Select based on RNA DE
rna_ct_de <- read_tsv('data_/results/diff_expression/RNA_DE_lineage.tsv')

mark_feats <- intersect(rownames(H3K27ac_use[['cRNA']]), rownames(H3K27me3_use[['cRNA']])) %>% intersect(rownames(H3K4me3_use[['cRNA']]))
groups_use <- c('ctx', 'dien', 'nt', 'retina')

region_markers <- rna_ct_de %>% 
    filter(p_val_adj<1e-5, gene%in%mark_feats, cluster%in%groups_use) %>% 
    group_by(cluster) %>% 
    dplyr::mutate(n_genes=n(), cluster=factor(cluster, levels=group_order)) %>% 
    arrange(dplyr::desc(cluster), dplyr::desc(avg_log2FC)) %>% 
    filter(row_number()<=15)

rna_use <- Pando::aggregate_assay(rna, group_name='lineage_coarse', assay='RNA')
H3K27ac_use <- Pando::aggregate_assay(H3K27ac_use, group_name='lineage_coarse', assay='cRNA')
H3K27me3_use <- Pando::aggregate_assay(H3K27me3_use, group_name='lineage_coarse', assay='cRNA')
H3K4me3_use <- Pando::aggregate_assay(H3K4me3_use, group_name='lineage_coarse', assay='cRNA')

rna_cluster_expr <- rna_use@assays$RNA@misc$summary$lineage_coarse[, unique(region_markers$gene)]
H3K27ac_cluster_expr <- H3K27ac_use@assays$cRNA@misc$summary$lineage_coarse[, unique(region_markers$gene)]
H3K27me3_cluster_expr <- H3K27me3_use@assays$cRNA@misc$summary$lineage_coarse[, unique(region_markers$gene)]
H3K4me3_cluster_expr <- H3K4me3_use@assays$cRNA@misc$summary$lineage_coarse[, unique(region_markers$gene)]

rna_cluster_expr_df <- rna_cluster_expr %>% as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr')
H3K27ac_cluster_expr_df <- H3K27ac_cluster_expr %>% as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr')
H3K27me3_cluster_expr_df <- H3K27me3_cluster_expr %>% as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr')
H3K4me3_cluster_expr_df <- H3K4me3_cluster_expr %>% as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr')


rna_plot_df <- rna_cluster_expr_df %>% 
    filter(group%in%groups_use, group%in%rna_ct_de$cluster) %>% 
    mutate(
        feature=factor(feature, levels=unique(region_markers$gene)),
        group=factor(group, levels=groups_use),
    ) %>% group_by(feature) %>% 
    dplyr::mutate(expr01=scale01(expr))

H3K27ac_plot_df <- H3K27ac_cluster_expr_df %>% 
    filter(group%in%groups_use, group%in%rna_ct_de$cluster) %>% 
    mutate(
        feature=factor(feature, levels=unique(region_markers$gene)),
        group=factor(group, levels=groups_use),
    ) %>% group_by(feature) %>% 
    dplyr::mutate(expr01=scale01(expr))

H3K27me3_plot_df <- H3K27me3_cluster_expr_df %>% 
    filter(group%in%groups_use, group%in%rna_ct_de$cluster) %>% 
    mutate(
        feature=factor(feature, levels=unique(region_markers$gene)),
        group=factor(group, levels=groups_use),
    ) %>% group_by(feature) %>% 
    dplyr::mutate(expr01=scale01(expr))

H3K4me3_plot_df <- H3K4me3_cluster_expr_df %>% 
    filter(group%in%groups_use, group%in%rna_ct_de$cluster) %>% 
    mutate(
        feature=factor(feature, levels=unique(region_markers$gene)),
        group=factor(group, levels=groups_use),
    ) %>% group_by(feature) %>% 
    dplyr::mutate(expr01=scale01(expr))

p1 <- ggplot(rna_plot_df, aes(group, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=ylorrd(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Gene', x='Group', fill='Gene expression', title='RNA') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() +
    theme(
        strip.text = element_blank()
    )

p1l <- ggplot(rna_plot_df, aes(group, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=ylorrd(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Gene', x='Group', fill='Gene expression', title='RNA') +
    article_text() +
    rotate_x_text(40) +
    theme(
        strip.text = element_blank()
    )


p2 <- ggplot(H3K27ac_plot_df, aes(group, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=ylgn(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Gene', x='Group', fill='Gene activity', title='H3K27ac') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() +
    theme(
        strip.text = element_blank(),
        axis.title.y = element_blank()
    )

p3 <- ggplot(H3K27me3_plot_df, aes(group, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=blues(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Gene', x='Group', fill='Gene activity', title='H3K27me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() +
    theme(
        strip.text = element_blank(),
        axis.title.y = element_blank()
    )


p4 <- ggplot(H3K4me3_plot_df, aes(group, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=rdpu(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Gene', x='Group', fill='Gene activity', title='H3K4me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() +
    theme(
        strip.text = element_blank(),
        axis.title.y = element_blank()
    )

(p1 | p2 | p3 | p4) + plot_layout(guides='collect')

ggsave('plots/paper/fig2/rev_fig2_top15_region_DE_unscaled_gene_heatmap.pdf', width=15, height=10, units='cm')
ggsave('plots/paper/fig2/rev_fig2_top15_region_DE_unscaled_gene_heatmap.png', width=15, height=10, units='cm')


(p1l | p2 | p3 | p4) + plot_layout(guides='collect')

ggsave('plots/paper/fig2/rev_fig2_top15_region_DE_labelled_unscaled_gene_heatmap.pdf', width=15, height=12, units='cm')
ggsave('plots/paper/fig2/rev_fig2_top15_region_DE_labelled_unscaled_gene_heatmap.png', width=15, height=12, units='cm')



#### Plot peak class distribution over stages ####
Kme_isect_matches <- read_tsv('data_/intersect/H3Kme_bivalent_marks_intersect_matches.tsv')
K27_isect_matches <- read_tsv('data_/intersect/H3K27_switches_marks_intersect_matches.tsv')
prom_isect_matches <- read_tsv('data_/intersect/H3K27K4_promoters_marks_intersect_matches.tsv')

H3K27ac_ctx <- read_rds('data_/trajectories/ctx/H3K27ac_eb_ctx_subs_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data_/trajectories/ctx/H3K4me3_eb_ctx_subs_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data_/trajectories/ctx/H3K27me3_eb_ctx_subs_dpt_srt.rds')

H3K27ac_ctx[['peaks_bin']] <- CreateAssayObject((H3K27ac_ctx[['peaks']]@data > 0)*1)
H3K4me3_ctx[['peaks_bin']] <- CreateAssayObject((H3K4me3_ctx[['peaks']]@data > 0)*1)
H3K27me3_ctx[['peaks_bin']] <- CreateAssayObject((H3K27me3_ctx[['peaks']]@data > 0)*1)

H3K27ac_ctx$ctx_pt <- rank(H3K27ac_ctx$ctx_pt) / max(rank(H3K27ac_ctx$ctx_pt))
H3K4me3_ctx$ctx_pt <- rank(H3K4me3_ctx$ctx_pt) / max(rank(H3K4me3_ctx$ctx_pt))
H3K27me3_ctx$ctx_pt <- rank(H3K27me3_ctx$ctx_pt) / max(rank(H3K27me3_ctx$ctx_pt))

H3K27ac_ctx$pt_bin2 <- as.numeric(cut(H3K27ac_ctx$ctx_pt, 25, labels=1:25))
H3K27ac_ctx$pt_bin <- as.numeric(cut(H3K27ac_ctx$ctx_pt, 50, labels=1:50))
H3K27ac_ctx$diff_stage <- case_when(
    H3K27ac_ctx$pt_bin <= 4 ~ 'psc',
    H3K27ac_ctx$pt_bin <= 27 ~ 'nepi',
    H3K27ac_ctx$pt_bin <= 44 ~ 'npc',
    T ~ 'neuron'
)
H3K4me3_ctx$pt_bin2 <- as.numeric(cut(H3K4me3_ctx$ctx_pt, 25, labels=1:25))
H3K4me3_ctx$pt_bin <- as.numeric(cut(H3K4me3_ctx$ctx_pt, 50, labels=1:50))
H3K4me3_ctx$diff_stage <- case_when(
    H3K4me3_ctx$pt_bin <= 4 ~ 'psc',
    H3K4me3_ctx$pt_bin <= 27 ~ 'nepi',
    H3K4me3_ctx$pt_bin <= 44 ~ 'npc',
    T ~ 'neuron'
)
H3K27me3_ctx$pt_bin2 <- as.numeric(cut(H3K27me3_ctx$ctx_pt, 25, labels=1:25))
H3K27me3_ctx$pt_bin <- as.numeric(cut(H3K27me3_ctx$ctx_pt, 50, labels=1:50))
H3K27me3_ctx$diff_stage <- case_when(
    H3K27me3_ctx$pt_bin <= 4 ~ 'psc',
    H3K27me3_ctx$pt_bin <= 27 ~ 'nepi',
    H3K27me3_ctx$pt_bin <= 44 ~ 'npc',
    T ~ 'neuron'
)

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='peaks_bin', group_name='clusters')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='peaks_bin', group_name='clusters')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='peaks_bin', group_name='clusters')

H3K27ac_ctx <- Pando::aggregate_assay(H3K27ac_ctx, assay='peaks_bin', group_name='pt_bin2')
H3K4me3_ctx <- Pando::aggregate_assay(H3K4me3_ctx, assay='peaks_bin', group_name='pt_bin2')
H3K27me3_ctx <- Pando::aggregate_assay(H3K27me3_ctx, assay='peaks_bin', group_name='pt_bin2')


p1 <- VlnPlot(H3K27me3_ctx, features='nCount_peaks', group.by='pt_bin2')
p2 <- VlnPlot(H3K4me3_ctx, features='nCount_peaks', group.by='pt_bin2')
p3 <- VlnPlot(H3K27ac_ctx, features='nCount_peaks', group.by='pt_bin2')

p1 / p2 / p3


H3K27ac_ctx_cluster_detect <- H3K27ac_ctx@assays$peaks_bin@misc$summary$clusters %>% {colnames(.)[colMaxs(.)>0.02]}
H3K4me3_ctx_cluster_detect <- H3K4me3_ctx@assays$peaks_bin@misc$summary$clusters %>% {colnames(.)[colMaxs(.)>0.02]}
H3K27me3_ctx_cluster_detect <- H3K27me3_ctx@assays$peaks_bin@misc$summary$clusters %>% {colnames(.)[colMaxs(.)>0.02]}

H3K27ac_ctx_stage_detect <- H3K27ac_ctx@assays$peaks_bin@misc$summary$pt_bin2[, H3K27ac_ctx_cluster_detect] %>% 
    as_tibble(rownames='stage') %>% 
    pivot_longer(!stage, names_to='peak', values_to='detect') 
H3K4me3_ctx_stage_detect <- H3K4me3_ctx@assays$peaks_bin@misc$summary$pt_bin2[, H3K4me3_ctx_cluster_detect] %>% 
    as_tibble(rownames='stage') %>% 
    pivot_longer(!stage, names_to='peak', values_to='detect') 
H3K27me3_ctx_stage_detect <- H3K27me3_ctx@assays$peaks_bin@misc$summary$pt_bin2[, H3K27me3_ctx_cluster_detect] %>% 
    as_tibble(rownames='stage') %>% 
    pivot_longer(!stage, names_to='peak', values_to='detect') 


n_bival_H3K27me3 <- StringToGRanges(Kme_isect_matches$H3K27me3)
n_bival_H3K4me3 <- StringToGRanges(Kme_isect_matches$H3K4me3)
Kme_isect_matches$dist <- IRanges::distance(n_bival_H3K27me3, n_bival_H3K4me3)

n_prom_H3K27ac <- StringToGRanges(prom_isect_matches$H3K27ac)
n_prom_H3K4me3 <- StringToGRanges(prom_isect_matches$H3K4me3)
prom_isect_matches$dist <- IRanges::distance(n_prom_H3K27ac, n_prom_H3K4me3)

all_stage_detect <- filter(prom_isect_matches, dist==0) %>% 
    full_join(filter(Kme_isect_matches, dist==0)) %>% 
    inner_join(H3K27ac_ctx_stage_detect, by=c('H3K27ac'='peak')) %>%
    inner_join(H3K27me3_ctx_stage_detect, by=c('stage', 'H3K27me3'='peak'), suffix=c('_H3K27ac','_H3K27me3')) %>% 
    inner_join(H3K4me3_ctx_stage_detect, by=c('stage', 'H3K4me3'='peak'), suffix=c('','_H3K4me3')) 

H3K27ac_q <- H3K27ac_ctx_stage_detect$detect %>% quantile(0.75)
H3K27me3_q <- H3K27me3_ctx_stage_detect$detect %>% quantile(0.75)
H3K4me3_q <- H3K4me3_ctx_stage_detect$detect %>% quantile(0.75)

# H3K27ac_q <- max(0.01, H3K27ac_q)
# H3K27me3_q <- max(0.01, H3K27me3_q)
# H3K4me3_q <- max(0.01, H3K4me3_q)

all_stage_class <- all_stage_detect %>% 
    mutate(
        detect_H3K27ac=replace_na(detect_H3K27ac, 0),
        detect_H3K27me3=replace_na(detect_H3K27me3, 0),
        detect=replace_na(detect, 0)
    ) %>% 
    # group_by(stage, isect) %>% 
    # summarize(
    #     detect_H3K27ac=max(detect_H3K27ac),
    #     detect_H3K27me3=max(detect_H3K27me3),
    #     detect=max(detect),
    #     max_detect=pmax(detect_H3K27me3, detect)
    # ) %>% 
    mutate(peak_class=case_when(
        (detect_H3K27me3 > H3K27me3_q) & (detect > H3K4me3_q) ~ 'bivalent',
        (detect > H3K4me3_q) & (detect_H3K27ac > H3K27ac_q) ~ 'promoter',
        detect_H3K27me3 > H3K27me3_q ~ 'H3K27me3',
        detect_H3K27ac > H3K27ac_q ~ 'H3K27ac',
        detect > H3K4me3_q ~ 'H3K4me3',
        T ~ 'gone'
    )) %>% 
    filter(peak_class!='gone') %>%
    mutate(peak_class=factor(peak_class, levels=c('H3K27me3', 'bivalent', 'H3K4me3', 'promoter', 'H3K27ac', 'gone'))) %>% 
    # mutate(stage=factor(stage, levels=c('psc', 'nepi', 'npc', 'neuron')))
    # mutate(stage=factor(stage, levels=c('EB', '15d', '35d', '60d', '128d', '4mo', '8mo')))
    # mutate(stage=factor(stage, levels=c('EB', 'mid', 'late', 'mo8')))
    mutate(stage=factor(stage, levels=c(1:25)))


modality_colors['bivalent'] <- '#9C6FF7'
modality_colors['gone'] <- 'black'
modality_colors['promoter'] <- '#E58B4A'


ggplot(all_stage_class, aes(stage, fill=peak_class)) +
    geom_bar(position='fill', alpha=0.8) +
    scale_fill_manual(values=modality_colors) +
    article_text()
ggsave('plots/paper/fig3/fig3_ctx_traject_age_peak_class_bar.pdf', width=6, height=5, units='cm')


ggplot(all_stage_class, aes(stage, fill=peak_class)) +
    geom_bar(alpha=0.8) +
    scale_fill_manual(values=modality_colors) +
    article_text()

ggplot(all_stage_class, aes(stage, max_detect, color=peak_class)) +
    geom_point() +
    scale_fill_manual(values=modality_colors) +
    article_text()

isect_ranges <- all_stage_class$isect %>% unique() %>% StringToGRanges()
isect_genes <- ClosestFeature(marks$H3K27me3, isect_ranges) %>% as_tibble()

all_class_genes <- all_stage_class %>% 
    inner_join(isect_genes, by=c('isect'='query_region'))

all_class_genes$gene_name %>% unique

all_class_genes %>% write_tsv('data_/trajectories/ctx/ctx_traject_ptbin25_0dist_peak_class.tsv')


all_class_genes <- read_tsv('data_/trajectories/ctx/ctx_traject_ptbin25_peak_class.tsv') %>% 
    mutate(peak_class=factor(peak_class, levels=c('H3K27me3', 'bivalent', 'H3K4me3', 'promoter', 'H3K27ac', 'gone')))


ctx_gene_classes <- all_class_genes %>% 
    dplyr::filter(stage>15) %>% 
    dplyr::distinct(stage, peak_class, gene_name) %>% 
    dplyr::group_by(gene_name) %>% 
    dplyr::mutate(
        frac_bivalent=sum(peak_class=='bivalent')/n(),
        frac_H3K27me3=sum(peak_class=='H3K27me3')/n(),
        frac_H3K27ac=sum(peak_class=='H3K27ac')/n(),
        frac_promoter=sum(peak_class=='promoter')/n()
    )

ggplot(all_class_genes, aes(stage, fill=peak_class)) +
    geom_bar(position='fill', alpha=0.8) +
    scale_fill_manual(values=modality_colors) +
    article_text()

b25_bival <- all_class_genes %>% filter(stage%in%c(23:25), peak_class=='bivalent') %>% pull(gene_name) %>% unique
b25_bival %>% write('data_/trajectories/ctx/ctx_traject_ptbin25_bival_neurons.txt')

full_H3K27me3 <- ctx_gene_classes %>% filter(frac_H3K27me3>0.9) %>% pull(gene_name) %>% setdiff(b25_bival) %>% unique()
full_H3K27me3 %>% write('data_/trajectories/ctx/ctx_traject_ptbin25_mostly_H3K27me3.txt')

full_H3K27ac <- ctx_gene_classes %>% filter(frac_H3K27ac>0.9) %>% pull(gene_name) %>% setdiff(b25_bival) %>% unique()
full_H3K27ac %>% write('data_/trajectories/ctx/ctx_traject_ptbin25_mostly_H3K27ac.txt')

dynamic_other <- ctx_gene_classes$gene_name %>% setdiff(b25_bival) %>% setdiff(full_H3K27me3) %>% setdiff(full_H3K27ac) %>% unique()
dynamic_other %>% write('data_/trajectories/ctx/ctx_traject_ptbin25_other_dynamic.txt')


nepi_bival <- all_class_genes %>% filter(stage%in%c(4:14), peak_class=='bivalent') %>% pull(gene_name) %>% unique
npc_bival <- all_class_genes %>% filter(stage%in%c(15:22), peak_class=='bivalent') %>% pull(gene_name) %>% unique

setdiff(b1_bival, b25_bival)

intersect(b25_bival, b1_bival)
intersect(nepi_bival, b1_bival)

setdiff(nepi_bival, b1_bival)




