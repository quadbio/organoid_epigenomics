library(tidyverse)
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(ggraph)
library(tidygraph)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

mutate <- dplyr::mutate
summarize <- dplyr::summarize
group_by <- dplyr::group_by
rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist
Matrix <- Matrix::Matrix

setwd('~/projects/cutntag/')

data(motif2tf)

#### Read data ####
marks <- read_rds('data05/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data05/RNA/RNA_all_srt_v2.2matched.rds')

links_grn <- read_rds('data_/GRN/RNA_ATAC_pseudocells_pando_links_srt.rds')
atac_grn <- read_rds('data_/GRN/RNA_ATAC_pseudocells_pando_atac_srt.rds')

links_grn <- find_modules(links_grn)
atac_grn <- find_modules(atac_grn)


modules <- NetworkModules(links_grn)



##### Show module changes over pt ####
peak_intersects <- read_tsv('data_/intersect/all_marks_intersect_matches.tsv')
switch_regions <- read_tsv('data_/intersect/H3K27_switches_marks_intersect_matches.tsv')
promoter_regions <- read_tsv('data_/intersect/H3K27K4_promoters_marks_intersect_matches.tsv')
bival_regions <- read_tsv('data_/intersect/H3Kme_bivalent_marks_intersect_matches.tsv')

gene_clusters <- read_tsv('data_/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')
neuron_genes <- gene_clusters %>% filter(dtw_clust==27) %>% pull(feature) %>% unique()


timelag_test <- read_tsv('data_/trajectories/ctx/ctx_pt_lr_binom_de.tsv')


min_bin <- timelag_test %>% 
    dplyr::filter(pval<0.01, log_dr>log2(1.5), modality!='H3K27me3') %>% 
    dplyr::group_by(feature) %>% 
    dplyr::filter(length(unique(modality))>1) %>%
    dplyr::group_by(modality, feature) %>% 
    dplyr::summarize(min_bin=min(pt_bin), min_dr=log_dr[which.min(pt_bin)]) %>% 
    dplyr::mutate(signed_bin=sign(min_dr)*min_bin) 

ggplot(filter(min_bin, modality=='RNA'), aes(min_bin, modality, label=feature)) +
    geom_text(position='jitter') +
    scale_fill_manual(values=modality_colors) + 
    scale_color_manual(values=modality_colors) + 
    scale_x_continuous(limits=c(0,10)) + 
    article_text() + no_legend() +
    labs(x='First divergent bin', y='Modality')


neuron_genes <- min_bin %>% arrange(min_bin) %>% pull(feature) %>% unique()

nd2_mod <- modules@meta %>% filter(target%in%neuron_genes)
nd2_regs <- nd2_mod$regions %>% str_split(';') %>% unlist()
nd2_ranges <- nd2_regs %>% StringToGRanges()

nd2_coef <- coef(links_grn) %>% 
    filter(target%in%nd2_mod$target, tf%in%nd2_mod$tf, padj<0.01, corr>0.2, region%in%nd2_regs, target%in%neuron_genes, estimate>0) %>% 
    group_by(tf) %>% 
    filter(length(unique(target))>2)



#### Epigen status of regions ####
#### Get regions ####
H3K27ac_peak_detect <- H3K27ac_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K27me3_peak_detect <- H3K27me3_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K4me3_peak_detect <- H3K4me3_ctx@assays$peaks_bin@misc$summary$pt_bins2[as.character(1:10), ]

H3K27ac_peaks <- H3K27ac_peak_detect %>% {.[,colMaxs(.)>0.05]} %>% colnames()
H3K27ac_ranges <- H3K27ac_peaks %>% StringToGRanges()
H3K27ac_olaps <- findOverlaps(nd2_ranges, H3K27ac_ranges)
H3K27ac_matches <- tibble(
    module_regions = nd2_regs[queryHits(H3K27ac_olaps)],
    H3K27ac_regions = H3K27ac_peaks[subjectHits(H3K27ac_olaps)]
)


H3K27me3_peaks <- H3K27me3_peak_detect %>% {.[,colMaxs(.)>0.05]} %>% colnames()
H3K27me3_ranges <- H3K27me3_peaks %>% StringToGRanges()
H3K27me3_olaps <- findOverlaps(nd2_ranges, H3K27me3_ranges)
H3K27me3_matches <- tibble(
    module_regions = nd2_regs[queryHits(H3K27me3_olaps)],
    H3K27me3_regions = H3K27me3_peaks[subjectHits(H3K27me3_olaps)]
)

H3K4me3_peaks <- H3K4me3_peak_detect %>% {.[,colMaxs(.)>0.05]} %>% colnames()
H3K4me3_ranges <- H3K4me3_peaks %>% StringToGRanges()
H3K4me3_olaps <- findOverlaps(nd2_ranges, H3K4me3_ranges)
H3K4me3_matches <- tibble(
    module_regions = nd2_regs[queryHits(H3K4me3_olaps)],
    H3K4me3_regions = H3K4me3_peaks[subjectHits(H3K4me3_olaps)]
)

all_matches <- full_join(H3K27ac_matches, H3K27me3_matches) %>% 
    full_join(H3K4me3_matches) %>% distinct()


#### Determine detection threshold ####
H3K27ac_match_detect <- H3K27ac_peak_detect[, unique(H3K27ac_matches$H3K27ac_regions)] %>% 
    as_tibble(rownames='pt_bin') %>% 
    pivot_longer(!pt_bin, names_to='H3K27ac_regions', values_to='H3K27ac_detect')

H3K27me3_match_detect <- H3K27me3_peak_detect[, unique(H3K27me3_matches$H3K27me3_regions)] %>% 
    as_tibble(rownames='pt_bin') %>% 
    pivot_longer(!pt_bin, names_to='H3K27me3_regions', values_to='H3K27me3_detect')

H3K4me3_match_detect <- H3K4me3_peak_detect[, unique(H3K4me3_matches$H3K4me3_regions)] %>% 
    as_tibble(rownames='pt_bin') %>% 
    pivot_longer(!pt_bin, names_to='H3K4me3_regions', values_to='H3K4me3_detect')


H3K27ac_thresh <- 0.05
H3K27me3_thresh <- 0.05
H3K4me3_thresh <- 0.05

all_matches_detect <- all_matches %>% 
    left_join(H3K27ac_match_detect) %>% 
    left_join(H3K27me3_match_detect) %>% 
    left_join(H3K4me3_match_detect) 

match_status <- all_matches_detect %>% 
    dplyr::group_by(module_regions, pt_bin) %>% 
    dplyr::summarize(
        H3K27ac=replace_na(H3K27ac_detect>H3K27ac_thresh, FALSE),
        H3K27me3=replace_na(H3K27me3_detect>H3K27me3_thresh, FALSE),
        H3K4me3=replace_na(H3K4me3_detect>H3K4me3_thresh, FALSE)
    ) %>% 
    dplyr::mutate(reg_status=case_when(
        H3K27ac & H3K4me3 ~ 'H3K4me3 + H3K27ac',
        H3K27me3 & H3K4me3 ~ 'H3K4me3 + H3K27me3',
        H3K27me3 ~ 'H3K27me3',
        H3K4me3 ~ 'H3K4me3',
        H3K27ac ~ 'H3K27ac',
        T ~ 'none'
    ))


modality_colors['H3K4me3 + H3K27me3'] <- '#9C6FF7'
modality_colors['H3K27me3'] <- '#027BB8'
modality_colors['none'] <- 'grey'
modality_colors['H3K4me3 + H3K27ac'] <- '#E58B4A'
modality_colors['H3K4me3'] <- '#E56F8E'
modality_colors['H3K27ac'] <- '#7CAF9C'



#### Join with coefs and plot ####
nd2_coef_status <- nd2_coef %>% inner_join(match_status, by=c('region'='module_regions'))

nd2_coef %>% write_tsv('data/GRN/neuro_grn_coefs.tsv')
nd2_coef_status %>% write_tsv('data/GRN/neuro_grn_pt_status.tsv')

nd2_graph <- nd2_coef_status %>%
    filter(estimate>0) %>%
    group_by(region) %>% 
    filter(sum(reg_status=='none')<9) %>% 
    as_tbl_graph() %E>%
    mutate(from_node = .N()$name[from]) %>% 
    filter(!is.na(pt_bin), !is.na(reg_status)) %N>%
    filter(centrality_degree(mode='out')>0 | centrality_degree(mode='in')>0) %>%
    mutate(
        y_coord=case_when(
            name %in% neuron_genes ~ -1,
            !(name %in% neuron_genes) ~ 1,
        ),
        centrality = centrality_eigen()
    ) %>%
    group_by(y_coord) %>%
    mutate(n_genes=n()) %>%
    mutate(x_coord = ifelse(!name%in%neuron_genes, as.numeric(factor(name))-(n_genes/2), as.numeric(factor(name, levels=neuron_genes))),) %>% 
    activate('edges') %>% 
    mutate(pt_bin=factor(as.numeric(pt_bin))) %>% 
    arrange(reg_status!='none')



set.seed(666)
ggraph(nd2_graph, layout='drl') +
    geom_edge_parallel(aes(color=reg_status, alpha = after_stat(index)), width=0.8, sep=unit(0.8, "mm")) +
    scale_edge_alpha('Edge direction', guide = 'edge_direction', range=c(0,1)) +
    geom_node_point(aes(size=centrality), color='black', fill='darkgrey', shape=21, stroke=0.2) +
    geom_node_text(aes(label=name), repel=T, size=5/ggplot2::.pt) +
    scale_edge_color_manual(values=modality_colors) +
    scale_size_continuous(range=c(0.2,2.5)) +
    facet_grid(~factor(as.numeric(pt_bin)))
ggsave('plots/paper/fig3/fig3_neurogen_GRN_all_bins_parallel_graph.png', width=60, height=6, bg='white', units='cm')
ggsave('plots/paper/fig3/fig3_neurogen_GRN_all_bins_parallel_graph.pdf', width=60, height=6, bg='white', units='cm')


set.seed(666)
ggraph(nd2_graph, layout='drl') +
    geom_edge_diagonal(aes(color=reg_status, alpha = after_stat(index)), width=0.8) +
    scale_edge_alpha('Edge direction', guide = 'edge_direction', range=c(0,1)) +
    geom_node_point(aes(size=centrality), color='black', fill='darkgrey', shape=21, stroke=0.2) +
    geom_node_text(aes(label=name), repel=T, size=5/ggplot2::.pt) +
    scale_edge_color_manual(values=modality_colors) +
    scale_size_continuous(range=c(0.2,2.5)) +
    facet_grid(~factor(as.numeric(pt_bin)))
ggsave('plots/paper/fig3/fig3_neurogen_GRN_all_bins_diagonal_graph.png', width=60, height=6, bg='white', units='cm')
ggsave('plots/paper/fig3/fig3_neurogen_GRN_all_bins_diagonal_graph.pdf', width=60, height=6, bg='white', units='cm')


for (i in 1:10){
    set.seed(666)
    ggraph(filter(activate(nd2_graph, 'edges'), pt_bin==i), layout='drl') +
        geom_edge_parallel(aes(color=reg_status, alpha = after_stat(index)), width=0.8, sep=unit(0.8, "mm")) +
        scale_edge_alpha('Edge direction', guide = 'edge_direction', range=c(0,1)) +
        geom_node_point(aes(size=centrality), color='black', fill='darkgrey', shape=21, stroke=0.2) +
        geom_node_text(aes(label=name), repel=T, size=5/ggplot2::.pt) +
        scale_edge_color_manual(values=modality_colors) +
        scale_size_continuous(range=c(0.2,2.5)) +
        facet_grid(~factor(as.numeric(pt_bin))) +
        theme_void() + no_legend()
    ggsave(paste0('plots/paper/fig3/grn_gif/fig3_neurogen_GRN_all_bin_', i, '.png'), width=6, height=6, bg='white', units='cm', dpi=600)
}


set.seed(666)
ggraph(filter(activate(nd2_graph, 'edges'), pt_bin%in%c(1,4,10)), layout='drl') +
    geom_edge_parallel(aes(color=reg_status, alpha = after_stat(index)), width=2) +
    scale_edge_alpha('Edge direction', guide = 'edge_direction', range=c(0,1)) +
    geom_node_point(aes(size=centrality), color='black', fill='darkgrey', shape=21) +
    geom_node_text(aes(label=name), repel=T, size=5/ggplot2::.pt) +
    scale_edge_color_manual(values=modality_colors) +
    facet_grid(~factor(as.numeric(pt_bin)))



set.seed(666)
ggraph(filter(activate(nd2_graph, 'edges'), pt_bin%in%c(1,4,10)), layout='drl') +
    geom_edge_diagonal(aes(color=reg_status, alpha = after_stat(index)), width=2) +
    scale_edge_alpha('Edge direction', guide = 'edge_direction', range=c(0,1)) +
    geom_node_point(aes(size=centrality), color='black', fill='darkgrey', shape=21) +
    geom_node_text(aes(label=name), repel=T) +
    scale_edge_color_manual(values=modality_colors) +
    facet_grid(~factor(as.numeric(pt_bin)))





#### Check expression of module genes ####
H3K27ac_ctx <- read_rds('data_/trajectories/ctx/H3K27ac_2m_ctx_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data_/trajectories/ctx/H3K4me3_2m_ctx_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data_/trajectories/ctx/H3K27me3_2m_ctx_dpt_srt.rds')
rna_ctx <- read_rds('data_/trajectories/ctx/RNA_2m_ctx_dpt_srt.rds')
muo_srt_ctx <- read_rds('data_/trajectories/ctx/MOU_RNA_4m_ctx_dpt_srt.rds')

H3K27ac_clusters <- H3K27ac_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K27me3_clusters <- H3K27me3_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
H3K4me3_clusters <- H3K4me3_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
rna_clusters <- rna_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]
muo_ga_clusters <- muo_srt_ctx@assays$gene_activity_bin@misc$summary$pt_bins2[as.character(1:10), ]
muo_rna_clusters <- muo_srt_ctx@assays$RNA_bin@misc$summary$pt_bins2[as.character(1:10), ]

genes_plot <- c('STMN2', 'BCL11A', 'NEUROD6', 'NEUROD2', 'POU5F1', 'SOX2', 'VIM', 'GLI3', 'GRIA2')

nd2_only_coef <- coef(links_grn) %>% 
    filter(padj<1e-4, corr>0.2, estimate>0, ((target=='NEUROD2' & target%in%nd2_coef_status$target) | tf=='NEUROD2')) 

out_genes <- nd2_only_coef %>% filter(tf=='NEUROD2') %>% pull(target) %>% intersect(colnames(H3K27ac_clusters))
in_genes <- nd2_only_coef %>% filter(target=='NEUROD2') %>% pull(tf) %>% intersect(colnames(H3K27ac_clusters))

genes_plot <- c(in_genes, out_genes, 'NEUROD2')

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

mod_expr_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, 'MUO_RNA'=muo_rna_expr, 'MUO_ATAC'=muo_ga_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr01=scale01(expr))

modality_colors <- c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044', 'MUO_ATAC'='#9575cd', 'MUO_RNA'='#FDDC44')


ggplot(mod_expr_df, aes(as.numeric(pt_bins2), expr01, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=modality_colors) +
    facet_wrap(~gene, scales='free')


expr_df <- rna_expr %>% 
    mutate(modality='RNA') %>% 
    group_by(modality, gene) %>% 
    mutate(expr01=scale01(expr)) 


ggplot(filter(mod_expr_df, gene=='NEUROD2', modality!='MUO_RNA'), aes(as.numeric(pt_bins2), expr01, color=modality)) +
    geom_smooth(mapping=aes(group=modality), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, linewidth=0.6) +
    # geom_point(mapping=aes(group=modality, fill=modality), size=1, shape=21, color='black', stroke=0.1) +
    article_text() +
    scale_color_manual(values=modality_colors) +
    scale_fill_manual(values=modality_colors) +
    scale_x_continuous(breaks=seq(1,10)) +
    labs(x='Pseudotime bins', y='Scaled detection')
ggsave('plots/paper/fig3/fig3_NEUROD2_pt_detection_line.pdf', units='cm', width=7, height=3)




ggplot(filter(mod_expr_df, gene=='NEUROD2', modality!='MUO_RNA'), aes(as.numeric(pt_bins2), expr01, color=modality)) +
    geom_smooth(data=filter(expr_df, gene%in%out_genes), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color='#bbdefb', linewidth=0.5) +
    geom_smooth(data=filter(expr_df, gene%in%in_genes), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color='#e1bee7', linewidth=0.5) +
    # geom_smooth(data=filter(expr_df, gene=='NEUROD2'), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color=modality_colors['RNA'], linewidth=0.6) +
    geom_smooth(data=filter(expr_df, gene=='NEUROD2', modality!='MUO_RNA'), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color='black', linewidth=0.6) +
    # geom_text(data=filter(expr_df, pt_bins2==10), mapping=aes(label=gene), color='black', x=11, size=5/ggplot2::.pt) +
    # scale_x_continuous(limits=c(1,12)) +
    scale_color_manual(values=modality_colors) +
    scale_x_continuous(breaks=seq(1,10)) +
    article_text() +
    labs(x='Pseudotime bins', y='Scaled detection')
ggsave('plots/paper/fig3/fig3_NEUROD2_pt_reg_rna_line.pdf', units='cm', width=5, height=3)





ggplot(filter(mod_expr_df, gene=='NEUROD2', modality!='MUO_RNA'), aes(as.numeric(pt_bins2), expr01, color=modality)) +
    # geom_smooth(data=filter(expr_df, gene%in%out_genes), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color='#bbdefb', linewidth=0.5) +
    geom_smooth(data=filter(expr_df, gene%in%in_genes), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color='#e1bee7', linewidth=0.5) +
    # geom_smooth(data=filter(expr_df, gene=='NEUROD2'), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color=modality_colors['RNA'], linewidth=0.6) +
    geom_smooth(data=filter(mod_expr_df, gene=='NEUROD2', modality!='MUO_RNA'), mapping=aes(group=modality), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, linewidth=0.6) +
    # geom_text(data=filter(expr_df, pt_bins2==10), mapping=aes(label=gene), color='black', x=11, size=5/ggplot2::.pt) +
    scale_color_manual(values=modality_colors) +
    scale_x_continuous(breaks=seq(1,10)) +
    article_text() +
    labs(x='Pseudotime bins', y='Scaled detection')
ggsave('plots/paper/fig3/fig3_NEUROD2_pt_upstream_rna_line.pdf', units='cm', width=7, height=3)



ggplot(filter(mod_expr_df, gene=='NEUROD2', modality!='MUO_RNA'), aes(as.numeric(pt_bins2), expr01, color=modality)) +
    geom_smooth(data=filter(expr_df, gene%in%out_genes), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color='#bbdefb', linewidth=0.5) +
    # geom_smooth(data=filter(expr_df, gene%in%in_genes), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color='#e1bee7', linewidth=0.5) +
    # geom_smooth(data=filter(expr_df, gene=='NEUROD2'), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color=modality_colors['RNA'], linewidth=0.6) +
    geom_smooth(data=filter(mod_expr_df, gene=='NEUROD2', modality!='MUO_RNA'), mapping=aes(group=modality), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, linewidth=0.6) +
    # geom_text(data=filter(expr_df, pt_bins2==10), mapping=aes(label=gene), color='black', x=11, size=5/ggplot2::.pt) +
    scale_color_manual(values=modality_colors) +
    scale_x_continuous(breaks=seq(1,10)) +
    article_text() +
    labs(x='Pseudotime bins', y='Scaled detection')
ggsave('plots/paper/fig3/fig3_NEUROD2_pt_downstream_rna_line.pdf', units='cm', width=7, height=3)




plots <- map(unique(genes_plot), function(g){
    color <- ifelse(g %in% out_genes, '#bbdefb', '#e1bee7')
    p <- ggplot(mod_expr_df, aes(as.numeric(pt_bins2), expr01, color=modality)) +
        geom_smooth(data=filter(mod_expr_df, gene=='NEUROD2'), mapping=aes(group=modality), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, linewidth=0.3) +
        geom_smooth(data=filter(expr_df, gene==g), mapping=aes(group=gene), method=mgcv::gam, formula = y ~ s(x, bs = 'cs'), se=F, color=color, linewidth=1) +
        scale_color_manual(values=modality_colors) +
        scale_x_continuous(breaks=seq(1,10)) +
        article_text() +
        no_legend() +
        labs(x='Pseudotime bins', y='Scaled detection', title=g)
    return(p)
})
wrap_plots(plots)

ggsave('plots/paper/fig3/fig3_NEUROD2_pt_reg_rna_all_genes_line.pdf', units='cm', width=25, height=15)








