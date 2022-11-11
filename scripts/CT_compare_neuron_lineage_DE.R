source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist
Matrix <- Matrix::Matrix

setwd('~/projects/cutntag/')

#### Read data ####
marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')

cluster_srt <- read_rds('data/RNA/all_RNA_marks_combined_clusters_srt.rds')
# cluster_srt %>% write_rds('data/RNA/all_RNA_marks_combined_clusters_srt.rds')
cluster_graph <- read_rds('data/RNA/all_RNA_cluster_graph.rds')

# tree_coords <- cluster_graph %>% as_tibble() %>% select(name, pseudotime_ranks, tree_pos) %>% 
#     column_to_rownames('name') %>% as.matrix()
# colnames(tree_coords) <- c('tree_1', 'tree_2')
# tree_coords[,2] <- 4-tree_coords[,2]
# 
# cluster_srt[['tree']] <- CreateDimReducObject(tree_coords)


# DE results 
H3K27ac_nt_neuron_da <- read_tsv('data/results/diff_expression/')


H3K27ac_nt_neuron_da <- read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_nt_NvsNPC.tsv')
H3K27ac_ctx_neuron_da <- read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_ctx_NvsNPC.tsv')
H3K27ac_dien_neuron_da <- read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_dien_NvsNPC.tsv')
H3K27ac_all_neuron_da <- read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_all_NvsNPC.tsv')

H3K27ac_neuron_da <- bind_rows('nt'=H3K27ac_nt_neuron_da, 'ctx'=H3K27ac_ctx_neuron_da, 'dien'=H3K27ac_dien_neuron_da, 'all'=H3K27ac_all_neuron_da, .id='lineage')

H3K27me3_nt_neuron_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_nt_NvsNPC.tsv')
H3K27me3_ctx_neuron_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_ctx_NvsNPC.tsv')
H3K27me3_dien_neuron_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_dien_NvsNPC.tsv')
H3K27me3_all_neuron_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_all_NvsNPC.tsv')

H3K27me3_neuron_da <- bind_rows('nt'=H3K27me3_nt_neuron_da, 'ctx'=H3K27me3_ctx_neuron_da, 'dien'=H3K27me3_dien_neuron_da, 'all'=H3K27me3_all_neuron_da, .id='lineage')

H3K4me3_nt_neuron_da <- read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_nt_NvsNPC.tsv')
H3K4me3_ctx_neuron_da <- read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_ctx_NvsNPC.tsv')
H3K4me3_dien_neuron_da <- read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_dien_NvsNPC.tsv')
H3K4me3_all_neuron_da <- read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_all_NvsNPC.tsv')

H3K4me3_neuron_da <- bind_rows('nt'=H3K4me3_nt_neuron_da, 'ctx'=H3K4me3_ctx_neuron_da, 'dien'=H3K4me3_dien_neuron_da, 'all'=H3K4me3_all_neuron_da, .id='lineage')

neuron_da <- bind_rows('H3K27ac'=H3K27ac_neuron_da, 'H3K27me3'=H3K27me3_neuron_da, 'H3K4me3'=H3K4me3_neuron_da, .id='mark')

peak_genes <- ClosestFeature(marks$H3K27me3, neuron_da$feature) %>% as_tibble() %>% distinct(query_region, gene_name)

neuron_da <- neuron_da %>% 
    inner_join(peak_genes, by=c('feature'='query_region'))
    
neuron_da %>% write_tsv('data/results/diff_expression/all_marks_DA_NvsNPC.tsv')


lineage_da <- read_tsv('data/results/diff_expression/all_marks_DA_lineage_coarse.tsv')

peak_genes <- ClosestFeature(marks$H3K27me3, lineage_da$feature) %>% as_tibble() %>% distinct(query_region, gene_name)

lineage_da <- lineage_da %>% 
    inner_join(peak_genes, by=c('feature'='query_region'))

lin_pairs_da <- bind_rows(
    H3K27ac=read_tsv('data/results/diff_expression/H3K27ac_DA_peaks_lineages_pairwise.tsv'),
    H3K27me3=read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_lineages_pairwise.tsv'),
    H3K4me3=read_tsv('data/results/diff_expression/H3K4me3_DA_peaks_lineages_pairwise.tsv'),
    .id='mark'
)

peak_genes <- ClosestFeature(marks$H3K27me3, lin_pairs_da$feature) %>% as_tibble() %>% distinct(query_region, gene_name)

lin_pairs_da <- lin_pairs_da %>% 
    inner_join(peak_genes, by=c('feature'='query_region'))

lin_pairs_da %>% write_tsv('data/results/diff_expression/all_marks_DA_peaks_lineages_pairwise.tsv')

#### Do DE for RNA ####

rna_nt <- read_rds('data/trajectories/nt/RNA_nt_dpt_srt.rds')
rna_dien <- read_rds('data/trajectories/dien/RNA_dien_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_ctx_dpt_srt.rds')
rna_ac <- read_rds('data/trajectories/astrocytes/RNA_astrocytes_npcs_srt.rds')


# NPC vs neuron
npc_neuron_nt_de <- de(rna_nt, 'state') %>% filter(group=='neuron')
npc_neuron_dien_de <- de(rna_dien, 'state') %>% filter(group=='neuron')
npc_neuron_ctx_de <- de(rna_ctx, 'state') %>% filter(group=='neuron')

rna_nt$lineage <- 'nt'
subs_npc <- sample(colnames(rna_nt)[rna_nt$state=='npc'], 2500)
subs_neuron <- sample(colnames(rna_nt)[rna_nt$state=='neuron'], 500)
rna_nt_subs <- subset(rna_nt, cells=c(subs_npc, subs_neuron))

rna_dien$lineage <- 'dien'
subs_cells <- sample(colnames(rna_dien), 3000)
rna_dien_subs <- subset(rna_dien, cells=subs_cells)

rna_ctx$lineage <- 'ctx'
subs_npc <- sample(colnames(rna_ctx)[rna_ctx$state=='npc'], 2500)
subs_neuron <- sample(colnames(rna_ctx)[rna_ctx$state=='neuron'], 500)
rna_ctx_subs <- subset(rna_ctx, cells=c(subs_npc, subs_neuron))


rna_lins <- merge(rna_nt_subs, list(rna_dien_subs, rna_ctx_subs))

npc_neuron_all <- de(rna_lins, 'state') %>% 
    filter(group=='neuron') %>% 
    select(feature, 'neuron_fc'=fc, 'neuron_pval'=padj)

lineage_combs <- list(c('nt', 'dien'), c('nt', 'ctx'), c('dien', 'ctx'))
lin_pairs_de <- map_dfr(lineage_combs, function(x){
    srt_use <- subset(rna_lins, lineage%in%c(x))
    srt_use$test_var <- srt_use$lineage==x[2]
    de(srt_use, 'test_var') %>% 
        mutate(group1=x[2], group2=x[1])
})

lin_pairs_de %>% write_tsv('data/results/diff_expression/RNA_DE_lineages_pairwise.tsv')


#### CTX vs DIEN ####
sig_colors <- c('#57B2B5', '#BD71CE', '#1660A1', 'darkgrey')
names(sig_colors) <- c('lin', 'neuron', 'both', 'none')

ctx_dien_lin_de <- lin_pairs_de %>% 
    filter(group1=='ctx', group2=='dien', group==T) %>% 
    select(feature, 'lin_fc'=fc, 'lin_pval'=padj)

ctx_dien_de <- inner_join(npc_neuron_all, ctx_dien_lin_de) %>% 
    filter(abs(lin_fc)>0 | abs(neuron_fc)>0) %>% 
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_fc)>0.25 | abs(lin_fc)>0.25),
        sig_group=case_when(
            !(lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'neuron',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & !(neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'lin',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'both',
            T ~ 'none'
        )
    )


ctx_dien_de %>% group_by(sig_group, sign(neuron_fc), sign(lin_fc)) %>% 
    filter(sig_group!='none') %>% 
    summarize(count=n())



label_data <- filter(ctx_dien_de, abs(lin_fc)>1 | abs(neuron_fc)>1)

p1 <- ggplot(ctx_dien_de, aes(lin_fc, neuron_fc, label=feature, alpha=sig, color=sig_group, size=abs(neuron_fc))) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_point() +
    scale_x_continuous(limits=c(-2,2)) +
    scale_y_continuous(limits=c(-2,2)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    geom_text(data=label_data, color='black', size=2) +
    scale_size_continuous(range=c(0.1,0.5)) +
    scale_color_manual(values=sig_colors) +
    article_text() +
    no_legend() +
    labs(x='<- DIEN vs CTX ->', y='<- NPC vs NEURON ->')


table(ctx_dien_de$sig_group)


ctx_dien_de_groups <- ctx_dien_de %>% 
    distinct(feature, sig_group) %>% 
    filter(sig_group!='none')

# -> compare with chromatin
ctx_dien_lin_da <- lin_pairs_da %>% 
    filter(group1=='ctx', group2=='dien') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'lin_coef'=coef, 'lin_pval'=padj, 'lin_dr'=log_dr, gene_name)

npc_neuron_da <- neuron_da %>% 
    filter(lineage=='all') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'neuron_coef'=coef, 'neuron_pval'=padj, 'neuron_dr'=log_dr, gene_name)

ctx_dien_da <- inner_join(npc_neuron_da, ctx_dien_lin_da) %>% 
    inner_join(ctx_dien_de_groups, by=c('gene_name'='feature')) %>%
    mutate(lin_coef_clip=clip_abs(lin_coef, 5), neuron_coef_clip=clip_abs(neuron_coef, 5)) %>% 
    filter(abs(lin_coef)>0 | abs(neuron_coef)>0) %>%
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_coef)>0.25 | abs(lin_coef)>0.25)
    )

jitterpos <- position_jitter(seed=111)

right_plot_df <- ctx_dien_da %>% mutate(sig_group=ifelse(sig_group=='lin', 'none', sig_group))
p_right <- ggplot(right_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(right_plot_df, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    rotate_x_text(90) +
    theme_rangeframe() + scale_axis_rangeframe() +
    no_legend() +
    no_label()

top_plot_df <- ctx_dien_da %>% mutate(sig_group=ifelse(sig_group=='neuron', 'none', sig_group))
p_top <- ggplot(top_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(ctx_dien_da, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    coord_flip() +
    no_legend() +
    no_label()

layout <- '
AAAAA#
BBBBBC
BBBBBC
BBBBBC
BBBBBC
BBBBBC
'

p_top + p1 + p_right + plot_layout(design = layout)

ggsave('plots/paper/sfig4/sfig4_ctx_dien_de_da_scatter.pdf', width=10, height=10, units='cm')




#### CTX vs NT ####
ctx_nt_lin_de <- lin_pairs_de %>% 
    filter(group1=='ctx', group2=='nt', group==T) %>% 
    select(feature, 'lin_fc'=fc, 'lin_pval'=padj)

ctx_nt_de <- inner_join(npc_neuron_all, ctx_nt_lin_de) %>% 
    filter(abs(lin_fc)>0 | abs(neuron_fc)>0) %>% 
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_fc)>0.25 | abs(lin_fc)>0.25),
        sig_group=case_when(
            !(lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'neuron',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & !(neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'lin',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'both',
            T ~ 'none'
        )
    )

ctx_nt_de %>% group_by(sig_group, sign(neuron_fc), sign(lin_fc)) %>% 
    filter(sig_group!='none') %>% 
    summarize(count=n())

label_data <- filter(ctx_nt_de, abs(lin_fc)>1 | abs(neuron_fc)>1)

p1 <- ggplot(ctx_nt_de, aes(lin_fc, neuron_fc, label=feature, alpha=sig, color=sig_group, size=abs(neuron_fc))) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_point() +
    scale_x_continuous(limits=c(-2,2)) +
    scale_y_continuous(limits=c(-2,2)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    geom_text(data=label_data, color='black', size=2) +
    scale_size_continuous(range=c(0.1,0.5)) +
    scale_color_manual(values=sig_colors) +
    article_text() +
    no_legend() +
    labs(x='<- DIEN vs NT ->', y='<- NPC vs NEURON ->')


ctx_nt_de_groups <- ctx_nt_de %>% 
    distinct(feature, sig_group) %>% 
    filter(sig_group!='none')

# -> compare with chromatin
ctx_nt_lin_da <- lin_pairs_da %>% 
    filter(group1=='ctx', group2=='nt') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'lin_coef'=coef, 'lin_pval'=padj, 'lin_dr'=log_dr, gene_name)

npc_neuron_da <- neuron_da %>% 
    filter(lineage=='all') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'neuron_coef'=coef, 'neuron_pval'=padj, 'neuron_dr'=log_dr, gene_name)

ctx_nt_da <- inner_join(npc_neuron_da, ctx_nt_lin_da) %>% 
    inner_join(ctx_nt_de_groups, by=c('gene_name'='feature')) %>%
    mutate(lin_coef_clip=clip_abs(lin_coef, 5), neuron_coef_clip=clip_abs(neuron_coef, 5)) %>% 
    filter(abs(lin_coef)>0 | abs(neuron_coef)>0) %>%
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_coef)>0.25 | abs(lin_coef)>0.25)
    )

jitterpos <- position_jitter(seed=111)

right_plot_df <- ctx_nt_da %>% mutate(sig_group=ifelse(sig_group=='lin', 'none', sig_group))
p_right <- ggplot(right_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(right_plot_df, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    rotate_x_text(90) +
    theme_rangeframe() + scale_axis_rangeframe() +
    no_legend() +
    no_label()

top_plot_df <- ctx_nt_da %>% mutate(sig_group=ifelse(sig_group=='neuron', 'none', sig_group))
p_top <- ggplot(top_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(ctx_nt_da, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    coord_flip() +
    no_legend() +
    no_label()

layout <- '
AAAAA#
BBBBBC
BBBBBC
BBBBBC
BBBBBC
BBBBBC
'

p_top + p1 + p_right + plot_layout(design = layout)

ggsave('plots/paper/sfig4/sfig4_ctx_nt_de_da_scatter.pdf', width=10, height=10, units='cm')




#### CTX vs NT ####
dien_nt_lin_de <- lin_pairs_de %>% 
    filter(group1=='dien', group2=='nt', group==T) %>% 
    select(feature, 'lin_fc'=fc, 'lin_pval'=padj)

dien_nt_de <- inner_join(npc_neuron_all, dien_nt_lin_de) %>% 
    filter(abs(lin_fc)>0 | abs(neuron_fc)>0) %>% 
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_fc)>0.25 | abs(lin_fc)>0.25),
        sig_group=case_when(
            !(lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'neuron',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & !(neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'lin',
            (lin_pval<1e-4 & abs(lin_fc)>0.25) & (neuron_pval<1e-4 & abs(neuron_fc)>0.25) ~ 'both',
            T ~ 'none'
        )
    )

dien_nt_de %>% group_by(sig_group, sign(neuron_fc), sign(lin_fc)) %>% 
    filter(sig_group!='none') %>% 
    summarize(count=n())

label_data <- filter(dien_nt_de, abs(lin_fc)>1 | abs(neuron_fc)>1)

p1 <- ggplot(dien_nt_de, aes(lin_fc, neuron_fc, label=feature, alpha=sig, color=sig_group, size=abs(neuron_fc))) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(0.25, -0.25), linetype='dashed', color='grey') +
    geom_point() +
    scale_x_continuous(limits=c(-2,2)) +
    scale_y_continuous(limits=c(-2,2)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    geom_text(data=label_data, color='black', size=2) +
    scale_size_continuous(range=c(0.1,0.5)) +
    scale_color_manual(values=sig_colors) +
    article_text() +
    no_legend() +
    labs(x='<- DIEN vs NT ->', y='<- NPC vs NEURON ->')


dien_nt_de_groups <- dien_nt_de %>% 
    distinct(feature, sig_group) %>% 
    filter(sig_group!='none')

# -> compare with chromatin
dien_nt_lin_da <- lin_pairs_da %>% 
    filter(group1=='dien', group2=='nt') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'lin_coef'=coef, 'lin_pval'=padj, 'lin_dr'=log_dr, gene_name)

npc_neuron_da <- neuron_da %>% 
    filter(lineage=='all') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    select(mark, feature, 'neuron_coef'=coef, 'neuron_pval'=padj, 'neuron_dr'=log_dr, gene_name)

dien_nt_da <- inner_join(npc_neuron_da, dien_nt_lin_da) %>% 
    inner_join(dien_nt_de_groups, by=c('gene_name'='feature')) %>%
    mutate(lin_coef_clip=clip_abs(lin_coef, 5), neuron_coef_clip=clip_abs(neuron_coef, 5)) %>% 
    filter(abs(lin_coef)>0 | abs(neuron_coef)>0) %>%
    mutate(
        sig=(neuron_pval<1e-4 | lin_pval<1e-4) & (abs(neuron_coef)>0.25 | abs(lin_coef)>0.25)
    )

jitterpos <- position_jitter(seed=111)

right_plot_df <- dien_nt_da %>% mutate(sig_group=ifelse(sig_group=='lin', 'none', sig_group))
p_right <- ggplot(right_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(right_plot_df, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    rotate_x_text(90) +
    theme_rangeframe() + scale_axis_rangeframe() +
    no_legend() +
    no_label()

top_plot_df <- dien_nt_da %>% mutate(sig_group=ifelse(sig_group=='neuron', 'none', sig_group))
p_top <- ggplot(top_plot_df, aes(mark, neuron_coef_clip, label=gene_name, alpha=sig, color=sig_group, size=abs(neuron_coef_clip))) +
    # geom_vline(xintercept = c(2, -2), linetype='dashed', color='grey') +
    geom_jitter(position=jitterpos, shape=16) +
    # geom_text(data=filter(dien_nt_da, gene_name%in%label_data$feature), color='black', size=2, position=jitterpos) +
    scale_size_continuous(range=c(0.1,0.2)) +
    scale_color_manual(values=sig_colors) +
    theme_rangeframe() + scale_axis_rangeframe() +
    scale_y_continuous(limits=c(-4,4)) +
    article_text() +
    coord_flip() +
    no_legend() +
    no_label()

layout <- '
AAAAA#
BBBBBC
BBBBBC
BBBBBC
BBBBBC
BBBBBC
'

p_top + p1 + p_right + plot_layout(design = layout)

ggsave('plots/paper/sfig4/sfig4_dien_nt_de_da_scatter.pdf', width=10, height=10, units='cm')







#### Compute explained variance for each peak over highres clusters ####
#### Neurogenesis ####
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

rna_nt <- read_rds('data/trajectories/nt/RNA_nt_dpt_srt.rds')
rna_dien <- read_rds('data/trajectories/dien/RNA_dien_dpt_srt.rds')
rna_ctx <- read_rds('data/trajectories/ctx/RNA_ctx_dpt_srt.rds')

H3K27ac_neuro_meta <- bind_rows(
    mutate(H3K27ac_nt@meta.data, neuro_pt=nt_pt), 
    mutate(H3K27ac_dien@meta.data, neuro_pt=dien_pt), 
    mutate(H3K27ac_ctx@meta.data, neuro_pt=ctx_pt)
) 

H3K27me3_neuro_meta <- bind_rows(
    mutate(H3K27me3_nt@meta.data, neuro_pt=nt_pt), 
    mutate(H3K27me3_dien@meta.data, neuro_pt=dien_pt), 
    mutate(H3K27me3_ctx@meta.data, neuro_pt=ctx_pt)
) 

H3K4me3_neuro_meta <- bind_rows(
    mutate(H3K4me3_nt@meta.data, neuro_pt=nt_pt), 
    mutate(H3K4me3_dien@meta.data, neuro_pt=dien_pt), 
    mutate(H3K4me3_ctx@meta.data, neuro_pt=ctx_pt)
) 

rna_neuro_meta <- bind_rows(
    mutate(rna_nt@meta.data, neuro_pt=rank(velocity_pseudotime)/max(rank(velocity_pseudotime))), 
    mutate(rna_dien@meta.data, neuro_pt=rank(velocity_pseudotime)/max(rank(velocity_pseudotime))), 
    mutate(rna_ctx@meta.data, neuro_pt=rank(velocity_pseudotime)/max(rank(velocity_pseudotime)))
) 

H3K27ac_neuro <- marks$H3K27ac %>% 
    subset(cells=rownames(H3K27ac_neuro_meta)) %>% 
    AddMetaData(H3K27ac_neuro_meta)

H3K27me3_neuro <- marks$H3K27me3 %>% 
    subset(cells=rownames(H3K27me3_neuro_meta)) %>% 
    AddMetaData(H3K27me3_neuro_meta)

H3K4me3_neuro <- marks$H3K4me3 %>% 
    subset(cells=rownames(H3K4me3_neuro_meta)) %>% 
    AddMetaData(H3K4me3_neuro_meta)

rna_neuro <- rna %>% 
    subset(cells=rownames(rna_neuro_meta)) %>% 
    AddMetaData(rna_neuro_meta)

H3K27ac_neuro %>% write_rds('data/trajectories/H3K27ac_neuro_dpt_srt.rds')
H3K27me3_neuro %>% write_rds('data/trajectories/H3K27me3_neuro_dpt_srt.rds')
H3K4me3_neuro %>% write_rds('data/trajectories/H3K4me3_neuro_dpt_srt.rds')
rna_neuro %>% write_rds('data/trajectories/RNA_neuro_dpt_srt.rds')


library(doParallel)
registerDoParallel(36)


rna_neuro <- Pando::aggregate_assay(rna_neuro, 'clusters', assay='RNA')
feats_use <- rna_neuro %>% FindVariableFeatures(nfeatures=4000) %>% VariableFeatures()
rna_neuro_cluster_expr <- rna_neuro@assays$RNA@misc$summary$clusters[, feats_use]
colnames(rna_neuro_cluster_expr) <- str_replace_all(colnames(rna_neuro_cluster_expr), '-', '_')
rna_neuro_cluster_meta <- rna_neuro@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    group_by(clusters) %>% 
    summarise(
        neuro_pt = mean(neuro_pt),
        lineage = lineage[1]
    )

model_frame <- data.frame(
    lineage=factor(rna_neuro_cluster_meta$lineage),
    pt=rna_neuro_cluster_meta$neuro_pt
)

module_mat <- rna_neuro_cluster_expr[rna_neuro_cluster_meta$clusters, ]

idx_use <- which(colMaxs(rna_neuro_cluster_expr)>0)
module_fit_list <- Pando::map_par(idx_use, function(i){
    # print(colnames(rna_neuro_cluster_expr)[i])
    model_formula <- reformulate(colnames(model_frame), response=colnames(rna_neuro_cluster_expr)[i])
    model_data <- cbind(rna_neuro_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})

module_lin_fit_list <- Pando::map_par(idx_use, function(i){
    model_formula <- reformulate('lineage', response=colnames(rna_neuro_cluster_expr)[i])
    model_data <- cbind(rna_neuro_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})

module_pt_fit_list <- Pando::map_par(idx_use, function(i){
    model_formula <- reformulate('pt', response=colnames(rna_neuro_cluster_expr)[i])
    model_data <- cbind(rna_neuro_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})


module_fit_r2 <- map_dfr(set_names(module_fit_list, colnames(rna_neuro_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

module_lin_r2 <- map_dfr(set_names(module_lin_fit_list, colnames(rna_neuro_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

module_pt_r2 <- map_dfr(set_names(module_pt_fit_list, colnames(rna_neuro_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

all_r2 <- bind_rows('lin_pt'=module_fit_r2, 'lin'=module_lin_r2, 'pt'=module_pt_r2, .id='vars') %>% 
    mutate(feature=str_replace_all(feature, '_', '-'))

all_r2 %>% write_tsv('data/trajectories/RNA_lineage_pt_r2.tsv')

plot_df <- all_r2 %>% 
    pivot_wider(!deviance, names_from=vars, values_from=r2)

ggplot(plot_df, aes(lin, pt, label=feature)) +
    geom_text()



# For all marks 
peak_var_list <- map(list('H3K27ac'=H3K27ac_neuro, 'H3K27me3'=H3K27me3_neuro, 'H3K4me3'=H3K4me3_neuro), function(mark){
    mark_neuro <- Pando::aggregate_assay(mark, 'clusters', assay='peaks_bin')
    peak_detection <- mark_neuro@assays$peaks_bin@misc$summary$clusters
    feats_use <- colnames(peak_detection)[colMaxs(peak_detection) > 0.05]
    mark_neuro <- Pando::aggregate_assay(mark, 'clusters', assay='peaks')
    mark_neuro_cluster_expr <- mark_neuro@assays$peaks@misc$summary$clusters[, feats_use]
    colnames(mark_neuro_cluster_expr) <- str_replace_all(colnames(mark_neuro_cluster_expr), '-', '_')
    mark_neuro_cluster_meta <- mark_neuro@meta.data %>% 
        as_tibble(rownames='cell') %>% 
        group_by(clusters) %>% 
        summarise(
            neuro_pt = mean(neuro_pt),
            lineage = lineage[1]
        )
    
    model_frame <- data.frame(
        lineage=factor(mark_neuro_cluster_meta$lineage),
        pt=mark_neuro_cluster_meta$neuro_pt
    )
    
    module_mat <- mark_neuro_cluster_expr[mark_neuro_cluster_meta$clusters, ]
    
    idx_use <- which(colMaxs(mark_neuro_cluster_expr)>0)
    module_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate(colnames(model_frame), response=colnames(mark_neuro_cluster_expr)[i])
        model_data <- cbind(mark_neuro_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_lin_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('lineage', response=colnames(mark_neuro_cluster_expr)[i])
        model_data <- cbind(mark_neuro_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_pt_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('pt', response=colnames(mark_neuro_cluster_expr)[i])
        model_data <- cbind(mark_neuro_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_fit_r2 <- map_dfr(set_names(module_fit_list, colnames(mark_neuro_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    module_lin_r2 <- map_dfr(set_names(module_lin_fit_list, colnames(mark_neuro_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    module_pt_r2 <- map_dfr(set_names(module_pt_fit_list, colnames(mark_neuro_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    bind_rows('lin_pt'=module_fit_r2, 'lin'=module_lin_r2, 'pt'=module_pt_r2, .id='vars') %>% 
        mutate(feature=str_replace_all(feature, '_', '-')) %>% 
        return()
})

peak_var_list %>% write_rds('data/trajectories/marks_lineage_pt_peaks_r2.rds')


plot_df <- peak_var_list$H3K4me3 %>% 
    pivot_wider(!deviance, names_from=vars, values_from=r2)

ggplot(plot_df, aes(lin, pt, label=feature)) +
    geom_text()




# For all marks gene activities
gene_var_list <- map(list('H3K27ac'=H3K27ac_neuro, 'H3K27me3'=H3K27me3_neuro, 'H3K4me3'=H3K4me3_neuro), function(mark){
    mark_neuro <- Pando::aggregate_assay(mark, 'clusters', assay='cRNA')
    feats_use <- mark_neuro %>% FindVariableFeatures(nfeatures=4000, assay='cRNA') %>% VariableFeatures(assay='cRNA')
    mark_neuro_cluster_expr <- mark_neuro@assays$cRNA@misc$summary$clusters[, feats_use]
    colnames(mark_neuro_cluster_expr) <- str_replace_all(colnames(mark_neuro_cluster_expr), '-', '_')
    mark_neuro_cluster_meta <- mark_neuro@meta.data %>% 
        as_tibble(rownames='cell') %>% 
        group_by(clusters) %>% 
        summarise(
            neuro_pt = mean(neuro_pt),
            lineage = lineage[1]
        )
    
    model_frame <- data.frame(
        lineage=factor(mark_neuro_cluster_meta$lineage),
        pt=mark_neuro_cluster_meta$neuro_pt
    )
    
    module_mat <- mark_neuro_cluster_expr[mark_neuro_cluster_meta$clusters, ]
    
    idx_use <- which(colMaxs(mark_neuro_cluster_expr)>0)
    module_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate(colnames(model_frame), response=colnames(mark_neuro_cluster_expr)[i])
        model_data <- cbind(mark_neuro_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_lin_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('lineage', response=colnames(mark_neuro_cluster_expr)[i])
        model_data <- cbind(mark_neuro_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_pt_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('pt', response=colnames(mark_neuro_cluster_expr)[i])
        model_data <- cbind(mark_neuro_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_fit_r2 <- map_dfr(set_names(module_fit_list, colnames(mark_neuro_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    module_lin_r2 <- map_dfr(set_names(module_lin_fit_list, colnames(mark_neuro_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    module_pt_r2 <- map_dfr(set_names(module_pt_fit_list, colnames(mark_neuro_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    bind_rows('lin_pt'=module_fit_r2, 'lin'=module_lin_r2, 'pt'=module_pt_r2, .id='vars') %>% 
        mutate(feature=str_replace_all(feature, '_', '-')) %>% 
        return()
})

gene_var_list %>% write_rds('data/trajectories/marks_lineage_pt_genes_r2.rds')


gene_var_list <- read_rds('data/trajectories/marks_lineage_pt_genes_r2.rds')
peak_var_list <- read_rds('data/trajectories/marks_lineage_pt_peaks_r2.rds')


gene_var_df <- bind_rows(gene_var_list, .id='modality') %>% rename('gene'=feature)
peak_var_df <- bind_rows(peak_var_list, .id='modality') %>% rename('peak'=feature)

peak_features <- ClosestFeature(marks$H3K27me3, unique(peak_var_df$peak)) %>% 
    as_tibble() %>% 
    select('gene'=gene_name, 'peak'=query_region)

peak_var_df <- inner_join(peak_var_df, peak_features)

all_r2$modality <- 'RNA'
all_r2 <- rename(all_r2, 'gene'=feature)
all_vars <- bind_rows('gene'=all_r2, 'gene'=gene_var_df, 'peak'=peak_var_df, .id='feature_type') %>% 
    select(gene, peak, vars, r2, deviance, modality, feature_type, everything())

all_vars %>% write_tsv('data/trajectories/ALL_lineage_pt_r2.tsv')



#### Plot some stuff ####
all_vars <- read_tsv('data/trajectories/ALL_lineage_pt_r2.tsv')

plot_df <- all_vars %>% 
    filter(feature_type=='gene', vars!='lin_pt') %>% 
    select(gene, vars, r2, modality) %>% 
    pivot_wider(names_from=modality, values_from=r2)
    
p1 <- ggplot(plot_df, aes(RNA, H3K27ac, label=gene)) +
    geom_point(size=0.4) +
    geom_text_repel(data=filter(plot_df, RNA>0.6), size=2, max.overlaps=99999) +
    facet_grid(~vars)

p2 <- ggplot(plot_df, aes(RNA, H3K4me3, label=gene)) +
    geom_point(size=0.4) +
    geom_text_repel(data=filter(plot_df, RNA>0.6), size=2, max.overlaps=99999) +
    facet_grid(~vars)

p3 <- ggplot(plot_df, aes(RNA, H3K27me3, label=gene)) +
    geom_point(size=0.4) +
    geom_text_repel(data=filter(plot_df, RNA>0.6), size=2, max.overlaps=99999) +
    facet_grid(~vars)

p1 / p2 / p3


peak_plot_df <- all_vars %>% 
    filter((feature_type=='peak' & modality=='H3K27ac'), vars!='lin_pt') %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=modality, values_from=r2)

gene_plot_df <- all_vars %>% 
    filter((feature_type=='gene' & modality=='RNA'), vars!='lin_pt') %>% 
    select(gene, vars, 'RNA'=r2) 

plot_df <- inner_join(peak_plot_df, gene_plot_df)
    
p1 <- ggplot(plot_df, aes(RNA, H3K27ac, label=gene)) +
    geom_point(size=0.4) +
    geom_text_repel(data=filter(plot_df, RNA>0.7), size=2, max.overlaps=99999) +
    facet_grid(~vars)

peak_plot_df <- all_vars %>% 
    filter((feature_type=='peak' & modality=='H3K4me3'), vars!='lin_pt') %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=modality, values_from=r2)

gene_plot_df <- all_vars %>% 
    filter((feature_type=='gene' & modality=='RNA'), vars!='lin_pt') %>% 
    select(gene, vars, 'RNA'=r2) 

plot_df <- inner_join(peak_plot_df, gene_plot_df)

p2 <- ggplot(plot_df, aes(RNA, H3K4me3, label=gene)) +
    geom_point(size=0.4) +
    geom_text_repel(data=filter(plot_df, RNA>0.7), size=2, max.overlaps=99999) +
    facet_grid(~vars)

peak_plot_df <- all_vars %>% 
    filter((feature_type=='peak' & modality=='H3K27me3'), vars!='lin_pt') %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=modality, values_from=r2)

gene_plot_df <- all_vars %>% 
    filter((feature_type=='gene' & modality=='RNA'), vars!='lin_pt') %>% 
    select(gene, vars, 'RNA'=r2) 

plot_df <- inner_join(peak_plot_df, gene_plot_df)

p3 <- ggplot(plot_df, aes(RNA, H3K27me3, label=gene)) +
    geom_point(size=0.4) +
    geom_text_repel(data=filter(plot_df, RNA>0.7), size=2, max.overlaps=99999) +
    facet_grid(~vars)

p1 / p2 / p3



plot_df <- all_vars %>% 
    filter(feature_type=='gene', vars!='lin_pt') %>% 
    select(gene, vars, r2, modality) %>% 
    pivot_wider(names_from=vars, values_from=r2) %>% 
    group_by(modality) %>% arrange(desc(lin)) %>% mutate(top_lin=row_number()<10) %>% 
    arrange(desc(pt)) %>% mutate(top_pt=row_number()<10)

ggplot(plot_df, aes(lin, pt, label=gene)) +
    geom_point(size=0.5) +
    geom_text_repel(data=filter(plot_df, top_lin | top_pt), size=2, max.overlaps=99999) +
    facet_grid(modality~.)



plot_df <- all_vars %>% 
    filter(feature_type=='peak', vars!='lin_pt', !is.na(peak)) %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=vars, values_from=r2) %>% 
    group_by(modality) %>% arrange(desc(lin)) %>% mutate(top_lin=row_number()<10) %>% 
    arrange(desc(pt)) %>% mutate(top_pt=row_number()<10)

ggplot(plot_df, aes(lin, pt, label=gene)) +
    geom_point(size=0.5) +
    geom_text_repel(data=filter(plot_df, top_lin | top_pt), size=2, max.overlaps=99999) +
    facet_grid(modality~.)


plot_df <- all_vars %>% 
    filter(feature_type=='gene')

ggplot(plot_df, aes(vars, r2, fill=modality)) +
    geom_boxplot(outlier.size = 0.2)


plot_df <- all_vars %>% 
    filter(feature_type=='peak', !is.na(peak))

ggplot(plot_df, aes(vars, r2, fill=modality)) +
    geom_boxplot(outlier.size = 0.2)



#### Astrocytes ####
rna_ac <- read_rds('data/trajectories/astrocytes/RNA_astrocytes_srt.rds')
H3K27ac_ac <- read_rds('data/trajectories/astrocytes/H3K27ac_astrocytes_srt.rds')
H3K27me3_ac <- read_rds('data/trajectories/astrocytes/H3K27me3_astrocytes_srt.rds')
H3K4me3_ac <- read_rds('data/trajectories/astrocytes/H3K4me3_astrocytes_srt.rds')

rna_ac_npc <- rna %>% subset(celltype_jf%in%c('ctx_npc', 'nt_npc', 'dien_npc', 'astrocytes'))
H3K27ac_ac_npc <- marks$H3K27ac %>% subset(celltype_jf%in%c('ctx_npc', 'nt_npc', 'dien_npc', 'astrocytes'))
H3K27me3_ac_npc <- marks$H3K27me3 %>% subset(celltype_jf%in%c('ctx_npc', 'nt_npc', 'dien_npc', 'astrocytes'))
H3K4me3_ac_npc <- marks$H3K4me3 %>% subset(celltype_jf%in%c('ctx_npc', 'nt_npc', 'dien_npc', 'astrocytes'))

rna_ac_npc <- AddMetaData(rna_ac_npc, rna_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K27ac_ac_npc <- AddMetaData(H3K27ac_ac_npc, H3K27ac_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K27me3_ac_npc <- AddMetaData(H3K27me3_ac_npc, H3K27me3_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K4me3_ac_npc <- AddMetaData(H3K4me3_ac_npc, H3K4me3_ac@meta.data[c('ac_region', 'seurat_clusters')])

rna_ac_npc$ac_clusters <- paste0(rna_ac_npc$clusters, '_', rna_ac_npc$seurat_clusters)
H3K27ac_ac_npc$ac_clusters <- paste0(H3K27ac_ac_npc$clusters, '_', H3K27ac_ac_npc$seurat_clusters)
H3K27me3_ac_npc$ac_clusters <- paste0(H3K27me3_ac_npc$clusters, '_', H3K27me3_ac_npc$seurat_clusters)
H3K4me3_ac_npc$ac_clusters <- paste0(H3K4me3_ac_npc$clusters, '_', H3K4me3_ac_npc$seurat_clusters)

H3K27ac_ac_npc@active.assay <- 'peaks'
H3K27ac_ac_npc_peaks_bin <- as(H3K27ac_ac_npc@assays$peaks@counts>0, 'dgCMatrix')
H3K27ac_ac_npc[['peaks_bin']] <- CreateAssayObject(H3K27ac_ac_npc_peaks_bin)

H3K27me3_ac_npc@active.assay <- 'peaks'
H3K27me3_ac_npc_peaks_bin <- as(H3K27me3_ac_npc@assays$peaks@counts>0, 'dgCMatrix')
H3K27me3_ac_npc[['peaks_bin']] <- CreateAssayObject(H3K27me3_ac_npc_peaks_bin)

H3K4me3_ac_npc@active.assay <- 'peaks'
H3K4me3_ac_npc_peaks_bin <- as(H3K4me3_ac_npc@assays$peaks@counts>0, 'dgCMatrix')
H3K4me3_ac_npc[['peaks_bin']] <- CreateAssayObject(H3K4me3_ac_npc_peaks_bin)


rna_ac_npc <- Pando::aggregate_assay(rna_ac_npc, 'ac_clusters', assay='RNA')
feats_use <- rna_ac_npc %>% FindVariableFeatures(nfeatures=4000) %>% VariableFeatures()
rna_ac_npc_cluster_expr <- rna_ac_npc@assays$RNA@misc$summary$ac_clusters[, feats_use]
colnames(rna_ac_npc_cluster_expr) <- str_replace(colnames(rna_ac_npc_cluster_expr), '-', '_')
rna_ac_npc_cluster_meta <- rna_ac_npc@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    group_by(ac_clusters) %>% 
    mutate(
        region=case_when(
            ac_region=='telen' ~ 'ctx',
            lineage=='ctx' ~ 'ctx',
            # ac_region=='other' ~ 'other',
            # ac_region%in%c('mesen', 'rhom') ~ 'nt',
            T ~ 'nt'
        ),
        ac_npc=case_when(
            celltype_jf=='astrocytes' ~ 'as',
            T ~ 'npc'
        )
    ) %>% 
    summarise(
        ac_npc = ac_npc[1],
        region = region[1]
    )

model_frame <- data.frame(
    ac_npc=factor(rna_ac_npc_cluster_meta$ac_npc),
    region=factor(rna_ac_npc_cluster_meta$region)
)

module_mat <- rna_ac_npc_cluster_expr[rna_ac_npc_cluster_meta$ac_clusters, ]

idx_use <- which(colMaxs(rna_ac_npc_cluster_expr)>0)
module_fit_list <- Pando::map_par(idx_use, function(i){
    model_formula <- reformulate(colnames(model_frame), response=colnames(rna_ac_npc_cluster_expr)[i])
    model_data <- cbind(rna_ac_npc_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})

module_ac_fit_list <- Pando::map_par(idx_use, function(i){
    model_formula <- reformulate('ac_npc', response=colnames(rna_ac_npc_cluster_expr)[i])
    model_data <- cbind(rna_ac_npc_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})

module_reg_fit_list <- Pando::map_par(idx_use, function(i){
    model_formula <- reformulate('region', response=colnames(rna_ac_npc_cluster_expr)[i])
    model_data <- cbind(rna_ac_npc_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})


module_fit_r2 <- map_dfr(set_names(module_fit_list, colnames(rna_ac_npc_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

module_ac_r2 <- map_dfr(set_names(module_ac_fit_list, colnames(rna_ac_npc_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(rna_ac_npc_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

all_r2 <- bind_rows('state_region'=module_fit_r2, 'state'=module_ac_r2, 'region'=module_reg_r2, .id='vars') %>% 
    mutate(feature=str_replace_all(feature, '_', '-'))

all_r2 %>% write_tsv('data/trajectories/RNA_astrocyte_telen_diff_r2.tsv')

plot_df <- all_r2 %>% 
    pivot_wider(!deviance, names_from=vars, values_from=r2)

ggplot(plot_df, aes(region, state, label=feature)) +
    geom_text()




# For all marks gene activities
gene_var_list <- map(list('H3K27ac'=H3K27ac_ac_npc, 'H3K27me3'=H3K27me3_ac_npc, 'H3K4me3'=H3K4me3_ac_npc), function(mark){
    mark_ac_npc <- Pando::aggregate_assay(mark, 'ac_clusters', assay='cRNA')
    feats_use <- mark_ac_npc %>% FindVariableFeatures(nfeatures=4000, assay='cRNA') %>% VariableFeatures(assay='cRNA')
    mark_ac_npc_cluster_expr <- mark_ac_npc@assays$cRNA@misc$summary$ac_clusters[, feats_use]
    colnames(mark_ac_npc_cluster_expr) <- str_replace_all(colnames(mark_ac_npc_cluster_expr), '-', '_')
    
    mark_ac_npc_cluster_meta <- mark_ac_npc@meta.data %>% 
        as_tibble(rownames='cell') %>% 
        group_by(ac_clusters) %>% 
        mutate(
            region=case_when(
                mode(ac_region)=='telen' ~ 'ctx',
                mode(ac_region)=='nt' ~ 'nt',
                lineage!='ctx' ~ 'nt',
                T ~ lineage
            ),
            ac_npc=case_when(
                celltype_jf=='astrocytes' ~ 'as',
                T ~ 'npc'
            )
        ) %>% 
        summarise(
            ac_npc = ac_npc[1],
            region = region[1]
        )
    
    model_frame <- data.frame(
        ac_npc=factor(mark_ac_npc_cluster_meta$ac_npc),
        region=factor(mark_ac_npc_cluster_meta$region)
    )
    
    module_mat <- mark_ac_npc_cluster_expr[mark_ac_npc_cluster_meta$ac_clusters, ]
    
    idx_use <- which(colMaxs(mark_ac_npc_cluster_expr)>0)
    module_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate(colnames(model_frame), response=colnames(mark_ac_npc_cluster_expr)[i])
        model_data <- cbind(mark_ac_npc_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_ac_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('ac_npc', response=colnames(mark_ac_npc_cluster_expr)[i])
        model_data <- cbind(mark_ac_npc_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_reg_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('region', response=colnames(mark_ac_npc_cluster_expr)[i])
        model_data <- cbind(mark_ac_npc_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    
    module_fit_r2 <- map_dfr(set_names(module_fit_list, colnames(mark_ac_npc_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    module_ac_r2 <- map_dfr(set_names(module_ac_fit_list, colnames(mark_ac_npc_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(mark_ac_npc_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    all_r2 <- bind_rows('state_region'=module_fit_r2, 'state'=module_ac_r2, 'region'=module_reg_r2, .id='vars') %>% 
        mutate(feature=str_replace_all(feature, '_', '-')) %>% 
        return()
    
})

gene_var_list %>% write_rds('data/trajectories/marks_astrocyte_telen_diff_genes_r2.rds')


peak_var_list <- map(list('H3K27ac'=H3K27ac_ac_npc, 'H3K27me3'=H3K27me3_ac_npc, 'H3K4me3'=H3K4me3_ac_npc), function(mark){
    mark_ac_npc <- Pando::aggregate_assay(mark, 'ac_clusters', assay='peaks_bin')
    peak_detection <- mark_ac_npc@assays$peaks_bin@misc$summary$ac_clusters
    feats_use <- colnames(peak_detection)[colMaxs(peak_detection) > 0.05]
    mark_ac_npc <- Pando::aggregate_assay(mark, 'ac_clusters', assay='peaks')
    mark_ac_npc_cluster_expr <- mark_ac_npc@assays$peaks@misc$summary$ac_clusters[, feats_use]
    colnames(mark_ac_npc_cluster_expr) <- str_replace_all(colnames(mark_ac_npc_cluster_expr), '-', '_')
    
    mark_ac_npc_cluster_meta <- mark_ac_npc@meta.data %>% 
        as_tibble(rownames='cell') %>% 
        group_by(ac_clusters) %>% 
        mutate(
            region=case_when(
                ac_region=='telen' ~ 'ctx',
                ac_region=='nt' ~ 'nt',
                lineage!='ctx' ~ 'nt',
                T ~ lineage
            ),
            ac_npc=case_when(
                celltype_jf=='astrocytes' ~ 'as',
                T ~ 'npc'
            )
        ) %>% 
        summarise(
            ac_npc = ac_npc[1],
            region = region[1]
        )
    
    model_frame <- data.frame(
        ac_npc=factor(mark_ac_npc_cluster_meta$ac_npc),
        region=factor(mark_ac_npc_cluster_meta$region)
    )
    
    module_mat <- mark_ac_npc_cluster_expr[mark_ac_npc_cluster_meta$ac_clusters, ]
    
    idx_use <- which(colMaxs(mark_ac_npc_cluster_expr)>0)
    module_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate(colnames(model_frame), response=colnames(mark_ac_npc_cluster_expr)[i])
        model_data <- cbind(mark_ac_npc_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_ac_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('ac_npc', response=colnames(mark_ac_npc_cluster_expr)[i])
        model_data <- cbind(mark_ac_npc_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_reg_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('region', response=colnames(mark_ac_npc_cluster_expr)[i])
        model_data <- cbind(mark_ac_npc_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    
    module_fit_r2 <- map_dfr(set_names(module_fit_list, colnames(mark_ac_npc_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    module_ac_r2 <- map_dfr(set_names(module_ac_fit_list, colnames(mark_ac_npc_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(mark_ac_npc_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    all_r2 <- bind_rows('state_region'=module_fit_r2, 'state'=module_ac_r2, 'region'=module_reg_r2, .id='vars') %>% 
        mutate(feature=str_replace_all(feature, '_', '-')) %>% 
        return()
    
})

peak_var_list %>% write_rds('data/trajectories/marks_astrocyte_telen_diff_peaks_r2.rds')




gene_var_df <- bind_rows(gene_var_list, .id='modality') %>% rename('gene'=feature)
peak_var_df <- bind_rows(peak_var_list, .id='modality') %>% rename('peak'=feature)

peak_features <- ClosestFeature(marks$H3K27me3, unique(peak_var_df$peak)) %>% 
    as_tibble() %>% 
    select('gene'=gene_name, 'peak'=query_region)

peak_var_df <- inner_join(peak_var_df, peak_features)

all_r2$modality <- 'RNA'
all_r2 <- rename(all_r2, 'gene'=feature)
all_vars <- bind_rows('gene'=all_r2, 'gene'=gene_var_df, 'peak'=peak_var_df, .id='feature_type') %>% 
    select(gene, peak, vars, r2, deviance, modality, feature_type, everything())

all_vars %>% write_tsv('data/trajectories/ALL_astrocyte_telen_diff_r2.tsv')



#### Plot some stuff ####
plot_df <- all_vars %>% 
    filter(feature_type=='gene', vars!='state_region') %>% 
    select(gene, vars, r2, modality) %>% 
    pivot_wider(names_from=modality, values_from=r2)

p1 <- ggplot(plot_df, aes(RNA, H3K27ac, label=gene)) +
    geom_point(size=0.4) +
    facet_grid(~vars)

p2 <- ggplot(plot_df, aes(RNA, H3K4me3, label=gene)) +
    geom_point(size=0.4) +
    facet_grid(~vars)

p3 <- ggplot(plot_df, aes(RNA, H3K27me3, label=gene)) +
    geom_point(size=0.4) +
    facet_grid(~vars)

p1 / p2 / p3




plot_df <- all_vars %>% 
    filter(feature_type=='peak', vars!='state_region', !is.na(peak)) %>% 
    select(peak, gene, vars, r2, modality) %>% 
    pivot_wider(names_from=vars, values_from=r2)

ggplot(plot_df, aes(state, region, label=gene)) +
    geom_hex(bins=100) +
    scale_fill_gradientn(colors=grad(pals::inferno, 0.8), trans='log') +
    facet_grid(modality~.)









#### Compute regional variance for NPC and astrocyte separately #####
#### NPC ####
rna_npc <- rna %>% subset(celltype_jf%in%c('ctx_npc', 'nt_npc', 'dien_npc'))
H3K27ac_npc <- marks$H3K27ac %>% subset(celltype_jf%in%c('ctx_npc', 'nt_npc', 'dien_npc'))
H3K27me3_npc <- marks$H3K27me3 %>% subset(celltype_jf%in%c('ctx_npc', 'nt_npc', 'dien_npc'))
H3K4me3_npc <- marks$H3K4me3 %>% subset(celltype_jf%in%c('ctx_npc', 'nt_npc', 'dien_npc'))

rna_npc <- AddMetaData(rna_npc, rna_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K27ac_npc <- AddMetaData(H3K27ac_npc, H3K27ac_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K27me3_npc <- AddMetaData(H3K27me3_npc, H3K27me3_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K4me3_npc <- AddMetaData(H3K4me3_npc, H3K4me3_ac@meta.data[c('ac_region', 'seurat_clusters')])

rna_npc$ac_clusters <- paste0(rna_npc$clusters, '_', rna_npc$seurat_clusters)
H3K27ac_npc$ac_clusters <- paste0(H3K27ac_npc$clusters, '_', H3K27ac_npc$seurat_clusters)
H3K27me3_npc$ac_clusters <- paste0(H3K27me3_npc$clusters, '_', H3K27me3_npc$seurat_clusters)
H3K4me3_npc$ac_clusters <- paste0(H3K4me3_npc$clusters, '_', H3K4me3_npc$seurat_clusters)

H3K27ac_npc@active.assay <- 'peaks'
H3K27ac_npc_peaks_bin <- as(H3K27ac_npc@assays$peaks@counts>0, 'dgCMatrix')
H3K27ac_npc[['peaks_bin']] <- CreateAssayObject(H3K27ac_npc_peaks_bin)

H3K27me3_npc@active.assay <- 'peaks'
H3K27me3_npc_peaks_bin <- as(H3K27me3_npc@assays$peaks@counts>0, 'dgCMatrix')
H3K27me3_npc[['peaks_bin']] <- CreateAssayObject(H3K27me3_npc_peaks_bin)

H3K4me3_npc@active.assay <- 'peaks'
H3K4me3_npc_peaks_bin <- as(H3K4me3_npc@assays$peaks@counts>0, 'dgCMatrix')
H3K4me3_npc[['peaks_bin']] <- CreateAssayObject(H3K4me3_npc_peaks_bin)


rna_npc <- Pando::aggregate_assay(rna_npc, 'ac_clusters', assay='RNA')
feats_use <- rna_npc %>% FindVariableFeatures(nfeatures=4000) %>% VariableFeatures()
rna_npc_cluster_expr <- rna_npc@assays$RNA@misc$summary$ac_clusters[, feats_use]
colnames(rna_npc_cluster_expr) <- str_replace(colnames(rna_npc_cluster_expr), '-', '_')
rna_npc_cluster_meta <- rna_npc@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    group_by(ac_clusters) %>% 
    mutate(
        region=case_when(
            ac_region=='telen' ~ 'ctx',
            lineage=='ctx' ~ 'ctx',
            # ac_region=='other' ~ 'other',
            # ac_region%in%c('mesen', 'rhom') ~ 'nt',
            T ~ 'nt'
        ),
        ac_npc=case_when(
            celltype_jf=='astrocytes' ~ 'as',
            T ~ 'npc'
        )
    ) %>% 
    summarise(
        ac_npc = ac_npc[1],
        region = region[1]
    )

model_frame <- data.frame(
    ac_npc=factor(rna_npc_cluster_meta$ac_npc),
    region=factor(rna_npc_cluster_meta$region)
)

module_mat <- rna_npc_cluster_expr[rna_npc_cluster_meta$ac_clusters, ]

idx_use <- which(colMaxs(rna_npc_cluster_expr)>0)
module_reg_fit_list <- Pando::map_par(idx_use, function(i){
    model_formula <- reformulate('region', response=colnames(rna_npc_cluster_expr)[i])
    model_data <- cbind(rna_npc_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})

module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(rna_npc_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

all_r2 <- bind_rows('region'=module_reg_r2, .id='vars') %>% 
    mutate(feature=str_replace_all(feature, '_', '-'))

all_r2 %>% write_tsv('data/trajectories/RNA_NPC_telen_diff_r2.tsv')



peak_var_list <- map(list('H3K27ac'=H3K27ac_npc, 'H3K27me3'=H3K27me3_npc, 'H3K4me3'=H3K4me3_npc), function(mark){
    mark_npc <- Pando::aggregate_assay(mark, 'ac_clusters', assay='peaks_bin')
    peak_detection <- mark_npc@assays$peaks_bin@misc$summary$ac_clusters
    feats_use <- colnames(peak_detection)[colMaxs(peak_detection) > 0.05]
    mark_npc <- Pando::aggregate_assay(mark, 'ac_clusters', assay='peaks')
    mark_npc_cluster_expr <- mark_npc@assays$peaks@misc$summary$ac_clusters[, feats_use]
    colnames(mark_npc_cluster_expr) <- str_replace_all(colnames(mark_npc_cluster_expr), '-', '_')
    
    mark_npc_cluster_meta <- mark_npc@meta.data %>% 
        as_tibble(rownames='cell') %>% 
        group_by(ac_clusters) %>% 
        mutate(
            region=case_when(
                ac_region=='telen' ~ 'ctx',
                ac_region=='nt' ~ 'nt',
                lineage!='ctx' ~ 'nt',
                T ~ lineage
            ),
            ac_npc=case_when(
                celltype_jf=='astrocytes' ~ 'as',
                T ~ 'npc'
            )
        ) %>% 
        summarise(
            ac_npc = ac_npc[1],
            region = region[1]
        )
    
    model_frame <- data.frame(
        ac_npc=factor(mark_npc_cluster_meta$ac_npc),
        region=factor(mark_npc_cluster_meta$region)
    )
    
    module_mat <- mark_npc_cluster_expr[mark_npc_cluster_meta$ac_clusters, ]
    
    idx_use <- which(colMaxs(mark_npc_cluster_expr)>0)
    module_reg_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('region', response=colnames(mark_npc_cluster_expr)[i])
        model_data <- cbind(mark_npc_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(mark_npc_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    all_r2 <- bind_rows('region'=module_reg_r2, .id='vars') %>% 
        mutate(feature=str_replace_all(feature, '_', '-')) %>% 
        return()
    
})

peak_var_list %>% write_rds('data/trajectories/marks_NPC_diff_peaks_r2.rds')


peak_var_list <- read_rds('data/trajectories/marks_NPC_diff_peaks_r2.rds')
all_r2 <- read_tsv('data/trajectories/RNA_NPC_telen_diff_r2.tsv')

peak_var_df <- bind_rows(peak_var_list, .id='modality') %>% rename('peak'=feature)

peak_features <- ClosestFeature(marks$H3K27me3, unique(peak_var_df$peak)) %>% 
    as_tibble() %>% 
    select('gene'=gene_name, 'peak'=query_region)

peak_var_df <- inner_join(peak_var_df, peak_features)

all_r2$modality <- 'RNA'
all_r2 <- rename(all_r2, 'gene'=feature)
all_vars <- bind_rows('gene'=all_r2, 'peak'=peak_var_df, .id='feature_type') %>% 
    select(gene, peak, vars, r2, deviance, modality, feature_type, everything())

all_vars %>% write_tsv('data/trajectories/ALL_NPC_diff_r2.tsv')





#### Astrocytes ####
rna_ac_only <- rna %>% subset(celltype_jf%in%c('astrocytes'))
H3K27ac_ac_only <- marks$H3K27ac %>% subset(celltype_jf%in%c('astrocytes'))
H3K27me3_ac_only <- marks$H3K27me3 %>% subset(celltype_jf%in%c('astrocytes'))
H3K4me3_ac_only <- marks$H3K4me3 %>% subset(celltype_jf%in%c('astrocytes'))

rna_ac_only <- AddMetaData(rna_ac_only, rna_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K27ac_ac_only <- AddMetaData(H3K27ac_ac_only, H3K27ac_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K27me3_ac_only <- AddMetaData(H3K27me3_ac_only, H3K27me3_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K4me3_ac_only <- AddMetaData(H3K4me3_ac_only, H3K4me3_ac@meta.data[c('ac_region', 'seurat_clusters')])

rna_ac_only$ac_clusters <- paste0(rna_ac_only$clusters, '_', rna_ac_only$seurat_clusters)
H3K27ac_ac_only$ac_clusters <- paste0(H3K27ac_ac_only$clusters, '_', H3K27ac_ac_only$seurat_clusters)
H3K27me3_ac_only$ac_clusters <- paste0(H3K27me3_ac_only$clusters, '_', H3K27me3_ac_only$seurat_clusters)
H3K4me3_ac_only$ac_clusters <- paste0(H3K4me3_ac_only$clusters, '_', H3K4me3_ac_only$seurat_clusters)

H3K27ac_ac_only@active.assay <- 'peaks'
H3K27ac_ac_only_peaks_bin <- as(H3K27ac_ac_only@assays$peaks@counts>0, 'dgCMatrix')
H3K27ac_ac_only[['peaks_bin']] <- CreateAssayObject(H3K27ac_ac_only_peaks_bin)

H3K27me3_ac_only@active.assay <- 'peaks'
H3K27me3_ac_only_peaks_bin <- as(H3K27me3_ac_only@assays$peaks@counts>0, 'dgCMatrix')
H3K27me3_ac_only[['peaks_bin']] <- CreateAssayObject(H3K27me3_ac_only_peaks_bin)

H3K4me3_ac_only@active.assay <- 'peaks'
H3K4me3_ac_only_peaks_bin <- as(H3K4me3_ac_only@assays$peaks@counts>0, 'dgCMatrix')
H3K4me3_ac_only[['peaks_bin']] <- CreateAssayObject(H3K4me3_ac_only_peaks_bin)


rna_ac_only <- Pando::aggregate_assay(rna_ac_only, 'ac_clusters', assay='RNA')
feats_use <- rna_ac_only %>% FindVariableFeatures(nfeatures=4000) %>% VariableFeatures()
rna_ac_only_cluster_expr <- rna_ac_only@assays$RNA@misc$summary$ac_clusters[, feats_use]
colnames(rna_ac_only_cluster_expr) <- str_replace(colnames(rna_ac_only_cluster_expr), '-', '_')
rna_ac_only_cluster_meta <- rna_ac_only@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    group_by(ac_clusters) %>% 
    mutate(
        region=case_when(
            ac_region=='telen' ~ 'ctx',
            lineage=='ctx' ~ 'ctx',
            # ac_region=='other' ~ 'other',
            # ac_region%in%c('mesen', 'rhom') ~ 'nt',
            T ~ 'nt'
        ),
        ac_ac_only=case_when(
            celltype_jf=='astrocytes' ~ 'as',
            T ~ 'npc'
        )
    ) %>% 
    summarise(
        ac_ac_only = ac_ac_only[1],
        region = region[1]
    )

model_frame <- data.frame(
    ac_ac_only=factor(rna_ac_only_cluster_meta$ac_ac_only),
    region=factor(rna_ac_only_cluster_meta$region)
)

module_mat <- rna_ac_only_cluster_expr[rna_ac_only_cluster_meta$ac_clusters, ]

idx_use <- which(colMaxs(rna_ac_only_cluster_expr)>0)
module_reg_fit_list <- Pando::map_par(idx_use, function(i){
    model_formula <- reformulate('region', response=colnames(rna_ac_only_cluster_expr)[i])
    model_data <- cbind(rna_ac_only_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})

module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(rna_ac_only_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

all_r2 <- bind_rows('region'=module_reg_r2, .id='vars') %>% 
    mutate(feature=str_replace_all(feature, '_', '-'))

all_r2 %>% write_tsv('data/trajectories/RNA_AConly_telen_diff_r2.tsv')



peak_var_list <- map(list('H3K27ac'=H3K27ac_ac_only, 'H3K27me3'=H3K27me3_ac_only, 'H3K4me3'=H3K4me3_ac_only), function(mark){
    mark_ac_only <- Pando::aggregate_assay(mark, 'ac_clusters', assay='peaks_bin')
    peak_detection <- mark_ac_only@assays$peaks_bin@misc$summary$ac_clusters
    feats_use <- colnames(peak_detection)[colMaxs(peak_detection) > 0.05]
    mark_ac_only <- Pando::aggregate_assay(mark, 'ac_clusters', assay='peaks')
    mark_ac_only_cluster_expr <- mark_ac_only@assays$peaks@misc$summary$ac_clusters[, feats_use]
    colnames(mark_ac_only_cluster_expr) <- str_replace_all(colnames(mark_ac_only_cluster_expr), '-', '_')
    
    mark_ac_only_cluster_meta <- mark_ac_only@meta.data %>% 
        as_tibble(rownames='cell') %>% 
        group_by(ac_clusters) %>% 
        mutate(
            region=case_when(
                ac_region=='telen' ~ 'ctx',
                ac_region=='nt' ~ 'nt',
                lineage!='ctx' ~ 'nt',
                T ~ lineage
            ),
            ac_ac_only=case_when(
                celltype_jf=='astrocytes' ~ 'as',
                T ~ 'npc'
            )
        ) %>% 
        summarise(
            ac_ac_only = ac_ac_only[1],
            region = region[1]
        )
    
    model_frame <- data.frame(
        ac_ac_only=factor(mark_ac_only_cluster_meta$ac_ac_only),
        region=factor(mark_ac_only_cluster_meta$region)
    )
    
    module_mat <- mark_ac_only_cluster_expr[mark_ac_only_cluster_meta$ac_clusters, ]
    
    idx_use <- which(colMaxs(mark_ac_only_cluster_expr)>0)
    module_reg_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('region', response=colnames(mark_ac_only_cluster_expr)[i])
        model_data <- cbind(mark_ac_only_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(mark_ac_only_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    all_r2 <- bind_rows('region'=module_reg_r2, .id='vars') %>% 
        mutate(feature=str_replace_all(feature, '_', '-')) %>% 
        return()
    
})

peak_var_list %>% write_rds('data/trajectories/marks_AConly_diff_peaks_r2.rds')





peak_var_list <- read_rds('data/trajectories/marks_AConly_diff_peaks_r2.rds')
all_r2 <- read_tsv('data/trajectories/RNA_AConly_telen_diff_r2.tsv')

peak_var_df <- bind_rows(peak_var_list, .id='modality') %>% rename('peak'=feature)

peak_features <- ClosestFeature(marks$H3K27me3, unique(peak_var_df$peak)) %>% 
    as_tibble() %>% 
    select('gene'=gene_name, 'peak'=query_region)

peak_var_df <- inner_join(peak_var_df, peak_features)

all_r2$modality <- 'RNA'
all_r2 <- rename(all_r2, 'gene'=feature)
all_vars <- bind_rows('gene'=all_r2, 'peak'=peak_var_df, .id='feature_type') %>% 
    select(gene, peak, vars, r2, deviance, modality, feature_type, everything())

all_vars %>% write_tsv('data/trajectories/ALL_AConly_diff_r2.tsv')








#### Neurons ####
rna_neurons <- rna %>% subset(celltype_jf%in%c('mesen_ex', 'ctx_ex', 'dien_ex', 'rhom_ex'))
H3K27ac_neurons <- marks$H3K27ac %>% subset(celltype_jf%in%c('mesen_ex', 'ctx_ex', 'dien_ex', 'rhom_ex'))
H3K27me3_neurons <- marks$H3K27me3 %>% subset(celltype_jf%in%c('mesen_ex', 'ctx_ex', 'dien_ex', 'rhom_ex'))
H3K4me3_neurons <- marks$H3K4me3 %>% subset(celltype_jf%in%c('mesen_ex', 'ctx_ex', 'dien_ex', 'rhom_ex'))

rna_neurons <- AddMetaData(rna_neurons, rna_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K27ac_neurons <- AddMetaData(H3K27ac_neurons, H3K27ac_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K27me3_neurons <- AddMetaData(H3K27me3_neurons, H3K27me3_ac@meta.data[c('ac_region', 'seurat_clusters')])
H3K4me3_neurons <- AddMetaData(H3K4me3_neurons, H3K4me3_ac@meta.data[c('ac_region', 'seurat_clusters')])

rna_neurons$ac_clusters <- paste0(rna_neurons$clusters, '_', rna_neurons$seurat_clusters)
H3K27ac_neurons$ac_clusters <- paste0(H3K27ac_neurons$clusters, '_', H3K27ac_neurons$seurat_clusters)
H3K27me3_neurons$ac_clusters <- paste0(H3K27me3_neurons$clusters, '_', H3K27me3_neurons$seurat_clusters)
H3K4me3_neurons$ac_clusters <- paste0(H3K4me3_neurons$clusters, '_', H3K4me3_neurons$seurat_clusters)

H3K27ac_neurons@active.assay <- 'peaks'
H3K27ac_neurons_peaks_bin <- as(H3K27ac_neurons@assays$peaks@counts>0, 'dgCMatrix')
H3K27ac_neurons[['peaks_bin']] <- CreateAssayObject(H3K27ac_neurons_peaks_bin)

H3K27me3_neurons@active.assay <- 'peaks'
H3K27me3_neurons_peaks_bin <- as(H3K27me3_neurons@assays$peaks@counts>0, 'dgCMatrix')
H3K27me3_neurons[['peaks_bin']] <- CreateAssayObject(H3K27me3_neurons_peaks_bin)

H3K4me3_neurons@active.assay <- 'peaks'
H3K4me3_neurons_peaks_bin <- as(H3K4me3_neurons@assays$peaks@counts>0, 'dgCMatrix')
H3K4me3_neurons[['peaks_bin']] <- CreateAssayObject(H3K4me3_neurons_peaks_bin)


rna_neurons <- Pando::aggregate_assay(rna_neurons, 'ac_clusters', assay='RNA')
feats_use <- rna_neurons %>% FindVariableFeatures(nfeatures=4000) %>% VariableFeatures()
rna_neurons_cluster_expr <- rna_neurons@assays$RNA@misc$summary$ac_clusters[, feats_use]
colnames(rna_neurons_cluster_expr) <- str_replace(colnames(rna_neurons_cluster_expr), '-', '_')
rna_neurons_cluster_meta <- rna_neurons@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    group_by(ac_clusters) %>% 
    mutate(
        region=case_when(
            ac_region=='telen' ~ 'ctx',
            lineage=='ctx' ~ 'ctx',
            # ac_region=='other' ~ 'other',
            # ac_region%in%c('mesen', 'rhom') ~ 'nt',
            T ~ 'nt'
        ),
        ac_neurons=case_when(
            celltype_jf=='astrocytes' ~ 'as',
            T ~ 'npc'
        )
    ) %>% 
    summarise(
        ac_neurons = ac_neurons[1],
        region = region[1]
    )

model_frame <- data.frame(
    ac_neurons=factor(rna_neurons_cluster_meta$ac_neurons),
    region=factor(rna_neurons_cluster_meta$region)
)

module_mat <- rna_neurons_cluster_expr[rna_neurons_cluster_meta$ac_clusters, ]

idx_use <- which(colMaxs(rna_neurons_cluster_expr)>0)
module_reg_fit_list <- Pando::map_par(idx_use, function(i){
    model_formula <- reformulate('region', response=colnames(rna_neurons_cluster_expr)[i])
    model_data <- cbind(rna_neurons_cluster_expr[,i,drop=F], model_frame)
    return(glm(formula=model_formula, data=model_data))
})

module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(rna_neurons_cluster_expr)[idx_use]), function(x){
    tibble(
        r2 = 1 - x$deviance/x$null.deviance,
        deviance = x$deviance
    )
}, .id='feature')

all_r2 <- bind_rows('region'=module_reg_r2, .id='vars') %>% 
    mutate(feature=str_replace_all(feature, '_', '-'))

all_r2 %>% write_tsv('data/trajectories/RNA_Neurons_telen_diff_r2.tsv')



peak_var_list <- map(list('H3K27ac'=H3K27ac_neurons, 'H3K27me3'=H3K27me3_neurons, 'H3K4me3'=H3K4me3_neurons), function(mark){
    mark_neurons <- Pando::aggregate_assay(mark, 'ac_clusters', assay='peaks_bin')
    peak_detection <- mark_neurons@assays$peaks_bin@misc$summary$ac_clusters
    feats_use <- colnames(peak_detection)[colMaxs(peak_detection) > 0.05]
    mark_neurons <- Pando::aggregate_assay(mark, 'ac_clusters', assay='peaks')
    mark_neurons_cluster_expr <- mark_neurons@assays$peaks@misc$summary$ac_clusters[, feats_use]
    colnames(mark_neurons_cluster_expr) <- str_replace_all(colnames(mark_neurons_cluster_expr), '-', '_')
    
    mark_neurons_cluster_meta <- mark_neurons@meta.data %>% 
        as_tibble(rownames='cell') %>% 
        group_by(ac_clusters) %>% 
        mutate(
            region=case_when(
                ac_region=='telen' ~ 'ctx',
                ac_region=='nt' ~ 'nt',
                lineage!='ctx' ~ 'nt',
                T ~ lineage
            ),
            ac_neurons=case_when(
                celltype_jf=='astrocytes' ~ 'as',
                T ~ 'npc'
            )
        ) %>% 
        summarise(
            ac_neurons = ac_neurons[1],
            region = region[1]
        )
    
    model_frame <- data.frame(
        ac_neurons=factor(mark_neurons_cluster_meta$ac_neurons),
        region=factor(mark_neurons_cluster_meta$region)
    )
    
    module_mat <- mark_neurons_cluster_expr[mark_neurons_cluster_meta$ac_clusters, ]
    
    idx_use <- which(colMaxs(mark_neurons_cluster_expr)>0)
    module_reg_fit_list <- Pando::map_par(idx_use, function(i){
        model_formula <- reformulate('region', response=colnames(mark_neurons_cluster_expr)[i])
        model_data <- cbind(mark_neurons_cluster_expr[,i,drop=F], model_frame)
        return(glm(formula=model_formula, data=model_data))
    })
    
    module_reg_r2 <- map_dfr(set_names(module_reg_fit_list, colnames(mark_neurons_cluster_expr)[idx_use]), function(x){
        tibble(
            r2 = 1 - x$deviance/x$null.deviance,
            deviance = x$deviance
        )
    }, .id='feature')
    
    all_r2 <- bind_rows('region'=module_reg_r2, .id='vars') %>% 
        mutate(feature=str_replace_all(feature, '_', '-')) %>% 
        return()
    
})

peak_var_list %>% write_rds('data/trajectories/marks_Neurons_diff_peaks_r2.rds')





peak_var_list <- read_rds('data/trajectories/marks_Neurons_diff_peaks_r2.rds')
all_r2 <- read_tsv('data/trajectories/RNA_Neurons_telen_diff_r2.tsv')

peak_var_df <- bind_rows(peak_var_list, .id='modality') %>% rename('peak'=feature)

peak_features <- ClosestFeature(marks$H3K27me3, unique(peak_var_df$peak)) %>% 
    as_tibble() %>% 
    select('gene'=gene_name, 'peak'=query_region)

peak_var_df <- inner_join(peak_var_df, peak_features)

all_r2$modality <- 'RNA'
all_r2 <- rename(all_r2, 'gene'=feature)
all_vars <- bind_rows('gene'=all_r2, 'peak'=peak_var_df, .id='feature_type') %>% 
    select(gene, peak, vars, r2, deviance, modality, feature_type, everything())

all_vars %>% write_tsv('data/trajectories/ALL_Neurons_diff_r2.tsv')










#### AC vs NPC vs Neurons ####

H3K27me3_ACvN_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_ACvsN.tsv')
H3K27me3_ACvNPC_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_ACvsNPC.tsv')
H3K27me3_NvNPC_da <- read_tsv('data/results/diff_expression/H3K27me3_DA_peaks_all_NvsNPC.tsv')


H3K27me3_ACvN_dr <- H3K27me3_ACvN_da %>% 
    filter(detect_self>0.02 | detect_other>0.02, pval<1e-4)

H3K27me3_ACvNPC_dr <- H3K27me3_ACvNPC_da %>% 
    filter(detect_self>0.02 | detect_other>0.02, pval<1e-4)

H3K27me3_NvNPC_dr <- H3K27me3_NvNPC_da %>% 
    filter(detect_self>0.02 | detect_other>0.02, pval<1e-4)

peaks_use <- c(H3K27me3_ACvN_dr$feature, H3K27me3_ACvNPC_dr$feature, H3K27me3_NvNPC_dr$feature) %>% unique()
genes_repr <- ClosestFeature(marks$H3K27me3, StringToGRanges(peaks_use))


rna_use <- rna %>% subset(age!='ret_12w' & age!='ret_6w')

rna_de <- de(rna_use, 'state') %>% filter(group!='other')

rna_markers <- rna_de %>% 
    filter(padj<1e-4, feature%in%rownames(marks$H3K27me3[['cRNA']])) %>% 
    group_by(group) %>% 
    top_n(200, fc) %>% 
    distinct(group, feature)

genes_use <- rna_markers %>% 
    pull(feature) %>% unique()

cluster_meta <- rna@meta.data %>% as_tibble() %>% distinct(clusters, state)

gene_expr <- rna[['RNA']]@misc$summary$clusters[, genes_use] %>% 
    as_tibble(rownames='clusters') %>% 
    pivot_longer(!clusters, names_to='feature', values_to='expr')

plot_df <- gene_expr %>% 
    inner_join(cluster_meta) %>% 
    inner_join(rna_markers) %>% 
    filter(state!='other', !str_detect(clusters, 'retina')) %>% 
    group_by(feature) %>% 
    mutate(
        expr=scale01(expr)
    )
    
p1 <- ggplot(plot_df, aes(clusters, feature, fill=expr)) +
    geom_tile() +
    facet_grid(group~state, space='free', scales='free') +
    scale_fill_gradientn(colors=ylorrd(0.8))




gene_expr <- marks$H3K27me3[['cRNA']]@misc$summary$clusters[, genes_use] %>% 
    as_tibble(rownames='clusters') %>% 
    pivot_longer(!clusters, names_to='feature', values_to='expr')

plot_df <- gene_expr %>% 
    inner_join(cluster_meta) %>% 
    inner_join(rna_markers) %>% 
    filter(state!='other', !str_detect(clusters, 'retina')) %>% 
    group_by(feature) %>% 
    mutate(
        expr=scale(expr)
    ) %>% 
    group_by(feature, state) %>% 
    mutate(mean_expr=mean(expr))
    
ggplot(plot_df, aes(feature, clusters, fill=expr)) +
    geom_tile() +
    facet_grid(state~group, space='free', scales='free') +
    scale_fill_gradientn(colors=blues())

p2 <- ggplot(plot_df, aes(expr, feature, fill=mean_expr)) +
    geom_boxplot(outlier.size=0.2) +
    geom_vline(xintercept = 0) +
    facet_grid(group~state, space='free', scales='free') +
    scale_fill_gradientn(colors=blues()) +
    labs()

p1 | p2


peak_genes <- ClosestFeature(marks$H3K27me3, StringToGRanges(peaks_use)) %>% 
    as_tibble()


plot_df <- H3K27me3_ACvN_da %>% 
    left_join(peak_genes, by=c('feature'='query_region')) %>% 
    left_join(rna_markers, by=c('gene_name'='feature'))

p1 <- ggplot(arrange(plot_df, !is.na(group)), aes(log_dr, -log10(pval), color=group, label=gene_name)) +
    geom_point() +
    geom_text_repel(data=filter(plot_df, !is.na(group)), max.overlaps=9999) +
    ggtitle('Neuron vs AC') +
    scale_y_continuous(limits=c(0,100))
    

plot_df <- H3K27me3_ACvNPC_da %>% 
    left_join(peak_genes, by=c('feature'='query_region')) %>% 
    left_join(rna_markers, by=c('gene_name'='feature'))

p2 <- ggplot(arrange(plot_df, !is.na(group)), aes(log_dr, -log10(pval), color=group, label=gene_name)) +
    geom_point() +
    geom_text_repel(data=filter(plot_df, !is.na(group)), max.overlaps=9999) +
    ggtitle('NPC vs AC') +
    scale_y_continuous(limits=c(0,100))

p1 | p2

















