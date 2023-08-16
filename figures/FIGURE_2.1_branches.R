library(tidyverse)
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist

setwd('~/projects/cutntag/')

mark_colors <- c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')


#### Read data ####
marks <- read_rds('data/all_marks_list_v3.3motifs.rds')

lineage_early_da_files <- list.files('data_/results/diff_expression', pattern='*lineage_early.tsv', full.names=T)
names(lineage_early_da_files) <- str_replace(lineage_early_da_files, '.+expression/(.+)_DA_.+', '\\1')
lineage_early_da_list <- lineage_early_da_files %>% map(read_tsv)

lineage_coarse_da_files <- list.files('data_/results/diff_expression', pattern='*lineage_coarse.tsv', full.names=T)
names(lineage_coarse_da_files) <- str_replace(lineage_coarse_da_files, '.+expression/(.+)_DA_.+', '\\1')
lineage_coarse_da_list <- lineage_coarse_da_files %>% map(read_tsv)

lineage_da_files <- list.files('data_/results/diff_expression', pattern='*lineage_all.tsv', full.names=T)
names(lineage_da_files) <- str_replace(lineage_da_files, '.+expression/(.+)_DA_.+', '\\1')
lineage_da_list <- lineage_da_files %>% map(read_tsv)

gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')

group_order <- c('psc', 'non_nect', 'nect', 'ctx', 'dien', 'nt', 'astrocytes', 'retina', 'choroid_plexus')

#### H3K27ac ####
#### Heatmaps with stage & lineage-specific peaks ####
H3K27ac_use <- marks$H3K27ac 
H3K27ac_use <- RunTFIDF(H3K27ac_use, assay='peaks')
H3K27ac_use$lineage_coarse <- case_when(
    H3K27ac_use$lineage != 'other' ~ H3K27ac_use$lineage,
    T ~ H3K27ac_use$celltype_jf
)
H3K27ac_use <- Pando::aggregate_assay(H3K27ac_use, group_name='lineage_coarse', assay='peaks')

H3K27ac_da <- lineage_coarse_da_list$H3K27ac %>% 
    filter(group%in%group_order) %>% 
    mutate(
        padj=p.adjust(pval, method='fdr'),
        log_dr=ifelse(abs(log_dr)>5, sign(log_dr)*5, log_dr),
        signed_p=sign(coef)*-log10(padj),
        p_clip=ifelse(abs(signed_p)>50, sign(signed_p)*50, signed_p)
    ) %>%
    group_by(feature) %>%
    dplyr::filter(any(padj<1e-4 & detect_self>0.05 & coef>0)) 

H3K27ac_de_mat <- H3K27ac_da %>%
    dplyr::select(group, feature, p_clip) %>%
    pivot_wider(names_from = group, values_from = p_clip, values_fill = 0) %>%
    column_to_rownames('feature') %>% as.matrix()

H3K27ac_peak_order <- H3K27ac_de_mat %>% scale() %>% dist() %>% hclust(method='ward.D2') %>% {.$label[.$order]}

H3K27ac_top <- H3K27ac_da %>% 
    group_by(group) %>% 
    top_n(50, coef)

#### H3K27me3 ####
#### Heatmaps with stage & lineage-specific peaks ####
H3K27me3_use <- marks$H3K27me3
H3K27me3_use <- RunTFIDF(H3K27me3_use, assay='peaks')
H3K27me3_use$lineage_coarse <- case_when(
    H3K27me3_use$lineage != 'other' ~ H3K27me3_use$lineage,
    T ~ H3K27me3_use$celltype_jf
)
H3K27me3_use <- Pando::aggregate_assay(H3K27me3_use, group_name='lineage_coarse', assay='peaks')

H3K27me3_da <- lineage_coarse_da_list$H3K27me3 %>% 
    filter(group%in%group_order) %>% 
    mutate(
        padj=p.adjust(pval, method='fdr'),
        log_dr=ifelse(abs(log_dr)>5, sign(log_dr)*5, log_dr),
        signed_p=sign(coef)*-log10(padj),
        p_clip=ifelse(abs(signed_p)>50, sign(signed_p)*50, signed_p)
    ) %>%
    group_by(feature) %>%
    dplyr::filter(any(padj<1e-4 & detect_self>0.05 & coef>0)) 

H3K27me3_de_mat <- H3K27me3_da %>%
    dplyr::select(group, feature, p_clip) %>%
    pivot_wider(names_from = group, values_from = p_clip, values_fill = 0) %>%
    column_to_rownames('feature') %>% as.matrix()

H3K27me3_peak_order <- H3K27me3_de_mat %>% scale() %>% dist() %>% hclust(method='ward.D2') %>% {.$label[.$order]}

H3K27me3_top <- H3K27me3_da %>% 
    group_by(group) %>% 
    top_n(50, coef)


#### H3K4me3 ####
#### Heatmaps with stage & lineage-specific peaks ####
H3K4me3_use <- marks$H3K4me3
H3K4me3_use <- RunTFIDF(H3K4me3_use, assay='peaks')
H3K4me3_use$lineage_coarse <- case_when(
    H3K4me3_use$lineage != 'other' ~ H3K4me3_use$lineage,
    T ~ H3K4me3_use$celltype_jf
)
H3K4me3_use <- Pando::aggregate_assay(H3K4me3_use, group_name='lineage_coarse', assay='peaks')

H3K4me3_da <- lineage_coarse_da_list$H3K4me3 %>% 
    filter(group%in%group_order) %>% 
    mutate(
        padj=p.adjust(pval, method='fdr'),
        log_dr=ifelse(abs(log_dr)>5, sign(log_dr)*5, log_dr),
        signed_p=sign(coef)*-log10(padj),
        p_clip=ifelse(abs(signed_p)>50, sign(signed_p)*50, signed_p)
    ) %>%
    group_by(feature) %>%
    dplyr::filter(any(padj<1e-4 & detect_self>0.05 & coef>0)) 

H3K4me3_de_mat <- H3K4me3_da %>%
    dplyr::select(group, feature, p_clip) %>%
    pivot_wider(names_from = group, values_from = p_clip, values_fill = 0) %>%
    column_to_rownames('feature') %>% as.matrix()

H3K4me3_peak_order <- H3K4me3_de_mat %>% scale() %>% dist() %>% hclust(method='ward.D2') %>% {.$label[.$order]}

H3K4me3_top <- H3K4me3_da %>% 
    group_by(group) %>% 
    top_n(50, coef)


#### Write all DA peaks ####
all_top_peaks <- bind_rows('H3K27me3'=H3K27me3_top, 'H3K27ac'=H3K27ac_top, 'H3K4me3'=H3K4me3_top, .id='mark')
all_top_peaks %>% write_tsv('data/results/diff_expression/all_marks_top_DA_lineage_coarse.tsv')

all_peaks <- bind_rows('H3K27me3'=H3K27me3_da, 'H3K27ac'=H3K27ac_da, 'H3K4me3'=H3K4me3_da, .id='mark')
all_peaks %>% write_tsv('data/results/diff_expression/all_marks_DA_lineage_coarse.tsv')


#### Peak expression heatmap with all marks ####
H3K27ac_cluster <- H3K27ac_use@assays$peaks@misc$summary$lineage_coarse[, unique(H3K27ac_top$feature)]

H3K27ac_peak_order <- H3K27ac_cluster %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

H3K27ac_cluster_expr <- H3K27ac_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr')

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

H3K27me3_peak_order <- H3K27me3_cluster %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

H3K27me3_cluster_expr <- H3K27me3_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr')

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

H3K4me3_peak_order <- H3K4me3_cluster %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

H3K4me3_cluster_expr <- H3K4me3_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='expr')

H3K4me3_plot_df <- H3K4me3_cluster_expr %>% 
    filter(feature%in%H3K4me3_top$feature, group%in%group_order) %>% 
    inner_join(select(H3K4me3_top, feature, 'top_group'=group)) %>% 
    mutate(
        feature=factor(feature, levels=H3K4me3_peak_order),
        group=factor(group, levels=group_order),
        top_group=factor(top_group, levels=group_order),
        expr_clip=ifelse(expr>1, 1, expr)
    ) 

p1 <- ggplot(H3K27ac_plot_df, aes(group, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=ylgn(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(top_group~., space='free', scales='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K27ac') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() +
    theme(
        strip.text = element_blank()
    )

p2 <- ggplot(H3K27me3_plot_df, aes(group, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=blues(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(top_group~., space='free', scales='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K27me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() +
    theme(
        strip.text = element_blank()
    )

p3 <- ggplot(H3K4me3_plot_df, aes(group, feature, fill=expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=rdpu(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(top_group~., space='free', scales='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K4me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() +
    theme(
        strip.text = element_blank()
    )

(p1 | p2 | p3) & theme(
    panel.spacing = unit(0.1, 'cm'),
    plot.margin = unit(rep(0.1,4), 'cm')
)
ggsave('plots/paper/fig2/fig2_top50_DA_peak_ylgn_heatmap.pdf', width=8, height=10)


#### GREAT enrichment for all peaks ####
library(rtracklayer)

rtracklayer::export.bed(
    StringToGRanges(unique(H3K27me3_da$feature)), 
    'data/results/diff_expression/bed/H3K27me3_background.bed'
)

H3K27me3_top_groups <- H3K27me3_da %>% 
    filter(padj<1e-4 & detect_self>0.05 & coef>0) %>% 
    group_by(group) %>% 
    top_n(50, coef) %>%
    group_split() 
    
map(H3K27me3_top_groups, function(x){
    peaks <- x$feature %>% unique()
    region <- x$group[1]
    rtracklayer::export.bed(
        StringToGRanges(peaks), 
        paste0('data/results/diff_expression/bed/H3K27me3_DA_', region, '_top50.bed')
    )
})

rtracklayer::export.bed(
    StringToGRanges(unique(H3K27ac_da$feature)), 
    'data/results/diff_expression/bed/H3K27ac_background.bed'
)


H3K27ac_top_groups <- H3K27ac_da %>% 
    filter(padj<1e-4 & detect_self>0.05 & coef>0) %>% 
    group_by(group) %>% 
    top_n(50, coef) %>%
    group_split() 
    
map(H3K27ac_top_groups, function(x){
    peaks <- x$feature %>% unique()
    region <- x$group[1]
    rtracklayer::export.bed(
        StringToGRanges(peaks), 
        paste0('data/results/diff_expression/bed/H3K27ac_DA_', region, '_top50.bed')
    )
})


rtracklayer::export.bed(
    StringToGRanges(unique(H3K4me3_da$feature)), 
    'data/results/diff_expression/bed/H3K4me3_background.bed'
)


H3K4me3_top_groups <- H3K4me3_da %>% 
    filter(padj<1e-4 & detect_self>0.05 & coef>0) %>% 
    group_by(group) %>% 
    top_n(50, coef) %>%
    group_split() 
    
map(H3K4me3_top_groups, function(x){
    peaks <- x$feature %>% unique()
    region <- x$group[1]
    rtracklayer::export.bed(
        StringToGRanges(peaks), 
        paste0('data/results/diff_expression/bed/H3K4me3_DA_', region, '_top50.bed')
    )
})


#### Read and analyze GREAT enrichment ####
great_dir <- 'data_/results/diff_expression/GREAT/'
great_results <- list.files(great_dir, full.names = T)
great_all_df <- great_results %>% 
    set_names() %>% 
    map_dfr(read_tsv, skip=3, .id='file') %>% 
    {colnames(.)[2] <- 'Ontology';.} %>% 
    filter(!str_detect(Ontology, '#')) %>% 
    mutate(
        mark=str_replace(file, '.+/GREAT//(H3K(27|4)(ac|me3))_.+_GREAT.tsv', '\\1'),
        celltype=str_replace(file, '.+/GREAT//H3K(27|4)(ac|me3)_(.+)_GREAT.tsv', '\\3'),
        celltype=ifelse(celltype=='actrocytes', 'astrocytes', celltype)
    )

great_all_df %>% write_tsv('data/results/diff_expression/GREAT/all_GREAT_enrich.tsv')


great_gobp <- great_all_df %>% filter(Ontology=='GO Biological Process')
great_gomf <- great_all_df %>% filter(Ontology=='GO Molecular Function')

great_gobp_top <- great_gobp %>% 
    filter(HyperFdrQ<0.05, FgRegionsHit>5) %>%
    group_by(mark, celltype) %>% 
    # top_n(10, RegionFoldEnrich) %>% 
    mutate(
        fe_clip=ifelse(RegionFoldEnrich>10, 10, RegionFoldEnrich),
        desc_abbr=str_trunc(Desc, 50, 'right')
    )

desc_order <- great_gobp_top %>% 
    arrange(desc(mark), RegionFoldEnrich) %>% 
    pull(desc_abbr) %>% unique()

ggplot(great_gobp_top, aes(RegionFoldEnrich, factor(desc_abbr, levels=desc_order))) +
    geom_bar(stat='identity') +
    facet_grid(celltype~mark, scales='free', space='free')
ggsave('plots/cutntag/all_marks_DA_GREAT.pdf', width=6, height=40)






###### Switches: H3K27ac vs H3K27me3 ####
# Get intersects 
H3K27me3_ranges <- StringToGRanges(rownames(marks$H3K27me3))
H3K27me3_ranges_extended <- H3K27me3_ranges %>% 
    Extend(upstream=2000, downstream=2000)

H3K27ac_ranges <- StringToGRanges(rownames(marks$H3K27ac))
H3K27ac_ranges_extended <- H3K27ac_ranges %>% 
    Extend(upstream=2000, downstream=2000)

all_intersect <- IRanges::intersect(H3K27me3_ranges_extended, H3K27ac_ranges_extended) 

H3K27me3_olaps <- findOverlaps(H3K27me3_ranges_extended, all_intersect)
H3K27me3_matches <- tibble(
    isect = GRangesToString(all_intersect)[subjectHits(H3K27me3_olaps)],
    H3K27me3 = rownames(marks$H3K27me3)[queryHits(H3K27me3_olaps)]
)

H3K27ac_olaps <- findOverlaps(H3K27ac_ranges_extended, all_intersect)
H3K27ac_matches <- tibble(
    isect = GRangesToString(all_intersect)[subjectHits(H3K27ac_olaps)],
    H3K27ac = rownames(marks$H3K27ac)[queryHits(H3K27ac_olaps)]
)

K27_isect_matches <- inner_join(H3K27me3_matches, H3K27ac_matches) 
K27_isect_matches %>% write_tsv('data/intersect/H3K27_switches_marks_intersect_matches.tsv')



#### Filter DA peaks by intersections ####
K27_isect_matches <- read_tsv('data_/intersect/H3K27_switches_marks_intersect_matches.tsv')
# all_lineage_da <- read_tsv('data/results/diff_expression/all_marks_lineage_DA.tsv')
# all_DA_df <- bind_rows('H3K27ac'=H3K27ac_da, 'H3K27me3'=H3K27me3_da, 'H3K4me3'=H3K4me3_da, .id='mark') %>% 
#     filter(group!='psc', feature%in%all_lineage_da$feature)
all_DA_df <- read_tsv('data_/results/diff_expression/all_marks_lineage_DA.tsv')


# H3K27ac
H3K27ac_da <- all_DA_df %>% filter(mark=='H3K27ac')
H3K27ac_isect_peaks <- intersect(H3K27ac_da$feature, K27_isect_matches$H3K27ac)

H3K27ac_da_isect <- H3K27ac_da %>%
    dplyr::filter(feature%in%H3K27ac_isect_peaks) %>% 
    mutate(mark='H3K27ac')

# H3K27me3
H3K27me3_da <- all_DA_df %>% filter(mark=='H3K27me3')
H3K27me3_isect_peaks <- intersect(H3K27me3_da$feature, K27_isect_matches$H3K27me3)

H3K27me3_da_isect <- H3K27me3_da %>%
    dplyr::filter(feature%in%H3K27me3_isect_peaks) %>% 
    mutate(mark='H3K27me3')

# Combine
K27_da_isect <- K27_isect_matches %>% 
    inner_join(H3K27ac_da_isect, by=c('H3K27ac'='feature')) %>% 
    inner_join(H3K27me3_da_isect, by=c('H3K27me3'='feature', 'group'), suffix=c('_H3K27ac', '_H3K27me3')) %>% 
    group_by(group, isect) %>% 
    filter(any(sign(p_clip_H3K27ac)!=sign(p_clip_H3K27me3)))

K27_da_isect %>% write_tsv('data/intersect/H3K27_switches_marks_DA_intersect.tsv')

K27_da_rbind <- inner_join(K27_isect_matches, H3K27ac_da_isect, by=c('H3K27ac'='feature')) %>% 
    bind_rows(inner_join(K27_isect_matches, H3K27me3_da_isect, by=c('H3K27me3'='feature'))) %>% 
    filter(isect%in%K27_da_isect$isect)

K27_da_mat <- K27_da_rbind %>%
    filter(group!='nect') %>% 
    mutate(mark_group=paste0(mark, '_', group)) %>% 
    dplyr::select(mark_group, isect, p_clip) %>%
    distinct() %>% 
    pivot_wider(names_from = mark_group, values_from = p_clip, values_fill = 0, values_fn=mean) %>%
    column_to_rownames('isect') %>% as.matrix()

isect_peak_order <- K27_da_mat %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

set.seed(111)
K27_kmeans_clust <- kmeans(K27_da_mat, centers=10)$cluster %>% 
    enframe('isect', 'cluster')



K27_da_isect_ranges <- K27_da_isect$isect %>% 
    StringToGRanges()
K27_da_isect_ranges %>% export.bed('data/intersect/H3K27_switches_marks_DA_intersect.bed')


#### Plot enrichment p-value heatmaps switches####
K27ac_max_df <- K27_da_isect %>% 
    group_by(isect) %>% 
    filter(abs(signed_p_H3K27ac)==max(abs(signed_p_H3K27ac))) %>% 
    filter(row_number()==1) %>% 
    arrange(group, p_clip_H3K27ac) %>% 
    select(isect, 'max_group'=group) 

plot_df <- K27_da_rbind %>% 
    inner_join(K27ac_max_df) %>% 
    inner_join(K27_kmeans_clust) %>% 
    # mutate(isect=factor(isect, levels=unique(K27ac_max_df$isect)))
    mutate(
        isect=factor(isect, levels=isect_peak_order),
        cluster=factor(cluster, levels=c(4,10,8,5,6,2,1,7,9,3)),
        max_group=factor(max_group, levels=c('nect', 'ctx', 'dien', 'nt', 'retina')),
        group=factor(group, levels=c('nect', 'ctx', 'dien', 'nt', 'retina'))
    )

p1 <- ggplot(plot_df, aes(group, isect, fill=p_clip)) +
    geom_tile() +
    facet_grid(cluster~mark, scales='free_y', space='free') +
    scale_fill_gradientn(colors=rev(pals::brewer.rdbu(100)), limits=c(-50, 50)) +
    no_y_text() +
    rotate_x_text(40) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Peak', x='Group', fill='Signed log10(FDR)') +
    no_legend()

p2 <- ggplot(plot_df, aes(mark, isect, fill=p_clip)) +
    geom_tile() +
    facet_grid(cluster~group, scales='free_y', space='free') +
    scale_fill_gradientn(colors=rev(pals::brewer.rdbu(100)), limits=c(-50, 50)) +
    no_y_text() +
    rotate_x_text(40) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Peak', x='Group', fill='Signed log10(FDR)')

p1 | p2

ggsave('plots/paper/fig2/fig2_H3K27me3_H3K27ac_switches_DA_peaks_pval_heatmap.pdf', width=6, height=10)


plot_df %>% write_tsv('data_/CT/switching_peaks.tsv')





#### Peak expression heatmaps switches ####
H3K27ac_cluster <- H3K27ac_use@assays$peaks@misc$summary$lineage_coarse[, unique(K27_da_rbind$H3K27ac)]
H3K27ac_cluster_expr <- H3K27ac_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='H3K27ac_expr')

H3K27me3_cluster <- H3K27me3_use@assays$peaks@misc$summary$lineage_coarse[, unique(K27_da_rbind$H3K27me3)]
H3K27me3_cluster_expr <- H3K27me3_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='H3K27me3_expr')

expr_plot_df <- K27_da_rbind %>% 
    select(-(group:p_clip)) %>% 
    distinct() %>% 
    inner_join(H3K27ac_cluster_expr, by=c('H3K27ac'='feature')) %>% 
    inner_join(H3K27me3_cluster_expr, by=c('H3K27me3'='feature', 'group')) %>% 
    filter(group%in%c('nect', 'ctx', 'dien', 'nt', 'retina')) %>% 
    inner_join(K27ac_max_df) %>% 
    inner_join(K27_kmeans_clust) %>% 
    mutate(
        isect=factor(isect, levels=isect_peak_order),
        cluster=factor(cluster, levels=c(4,10,8,5,6,2,1,7,9,3)),
        group=factor(group, levels=group_order)
    ) %>% group_by(isect, mark, group, max_group, cluster) %>% 
    summarize(H3K27ac_expr=mean(H3K27ac_expr), H3K27me3_expr=mean(H3K27me3_expr))

p1 <- ggplot(filter(expr_plot_df, mark=='H3K27ac'), aes(group, isect, fill=H3K27ac_expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=ylgn(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(cluster~., scales='free_y', space='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K27ac') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() 

p1

p2 <- ggplot(filter(expr_plot_df, mark=='H3K27me3'), aes(group, isect, fill=H3K27me3_expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=blues(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(cluster~., scales='free_y', space='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K27me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() 

p2

(p1 | p2) & theme(
    panel.spacing = unit(0.1, 'cm'),
    plot.margin = unit(rep(0.1,4), 'cm')
)
ggsave('plots/paper/fig2/fig2_H3K27me3_H3K27ac_switches_DA_peaks_expr_heatmap.pdf', width=8, height=10)



#### Get regions and coverage for switches ####
regions_use <- unique(K27_da_rbind$isect) %>% 
    StringToGRanges() %>% 
    # resize(fix='center', width = 1) %>% 
    Extend(upstream=1000, downstream=1000)

peak_genes <- Signac::ClosestFeature(marks$H3K27ac, regions_use) %>% 
    as_tibble() %>% mutate(feature=unique(K27_da_rbind$isect)) %>% 
    distinct(feature, gene_name)

marks$H3K27ac$lin_stage <- case_when(
    marks$H3K27ac$lineage != 'other' ~ marks$H3K27ac$lineage,
    marks$H3K27ac$stage != 'mid' ~ 'nepi',
    T ~ 'other'
)

marks$H3K27me3$lin_stage <- case_when(
    marks$H3K27me3$lineage != 'other' ~ marks$H3K27me3$lineage,
    marks$H3K27me3$stage != 'mid' ~ 'nepi',
    T ~ 'other'
)

H3K27ac_lin_nepi <- subset(marks$H3K27ac, lin_stage!='other')
H3K27me3_lin_nepi <- subset(marks$H3K27me3, lin_stage!='other')

names(regions_use) <- unique(K27_da_rbind$isect)

H3K27ac_cov_df <- get_coverage(H3K27ac_lin_nepi, regions=regions_use, group_by='lin_stage')
H3K27me3_cov_df <- get_coverage(H3K27me3_lin_nepi, regions=regions_use, group_by='lin_stage')


K27_plot_df <- bind_rows('H3K27ac'=H3K27ac_cov_df,'H3K27me3'=H3K27me3_cov_df, .id='mark') %>% 
    inner_join(peak_genes, by=c('region'='feature')) %>% 
    inner_join(K27ac_max_df, by=c('region'='isect')) %>% 
    group_by(gene_name) %>% 
    mutate(feature_num=as.numeric(factor(region))) %>% 
    group_by(region, mark) %>%
    mutate(
        gene_peak=paste0(gene_name, '_', feature_num),
        coverage=coverage/max(coverage), 
        region=factor(region, levels=isect_peak_order),
        group=factor(group, levels=c('nepi', 'ctx', 'dien', 'nt', 'retina'))
    )  %>% 
    arrange(region) %>% 
    mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak)))


K27_cov_mean_mat <- K27_plot_df %>% 
    mutate(group_mark = paste0(group, mark)) %>% 
    group_by(group_mark, region) %>% 
    summarize(coverage=mean(coverage)) %>%
    pivot_wider(names_from=group_mark, values_from=coverage) %>% 
    column_to_rownames('region') %>% as.matrix()

K27_cov_mat <- K27_plot_df %>% 
    mutate(group_pos=paste0(mark, group, poscut)) %>% 
    ungroup() %>% 
    select(region, group_pos, coverage) %>% 
    pivot_wider(names_from = group_pos, values_from = coverage, values_fill=0) %>% 
    column_to_rownames('region') %>% as.matrix()

peak_order <- K27_cov_mat %>% dist() %>% hclust() %>% {.$label[.$order]}
# peak_order <- K27_cov_mean_mat %>% dist() %>% hclust() %>% {.$labels[.$order]}

set.seed(111)
peak_clusts <- K27_cov_mean_mat %>% kmeans(centers = 10) %>% {.$cluster} %>% enframe('region', 'kmeans_clust')

H3K27ac_da_isect <- K27_da_rbind %>% filter(mark=='H3K27ac') %>% 
    distinct(isect, group, p_clip) %>% 
    group_by(isect, group) %>% 
    summarise(p_clip=mean(p_clip))

H3K27me3_da_isect <- K27_da_rbind %>% filter(mark=='H3K27me3') %>% 
    distinct(isect, group, p_clip) %>% 
    group_by(isect, group) %>% 
    summarise(p_clip=mean(p_clip))


H3K27ac_plot_df <- K27_plot_df %>% 
    filter(mark=='H3K27ac') %>% 
    left_join(H3K27ac_da_isect, by=c('region'='isect', 'group')) %>% 
    inner_join(K27_kmeans_clust, by=c('region'='isect')) %>% 
    inner_join(peak_clusts) %>% 
    mutate(
        p_clip=ifelse(is.na(p_clip), 0, p_clip),
        group=factor(group, levels=c('nepi', 'ctx', 'dien', 'nt', 'retina')),
        region=factor(region, levels=peak_order),
        cluster=factor(cluster, levels=c(4,10,8,5,6,2,1,7,9,3)),
        kmeans_clust=factor(kmeans_clust, levels=c(6,5,7,4,2,1,8,10,9,3))
    ) %>% 
    filter(!is.na(group))

H3K27me3_plot_df <- K27_plot_df %>% 
    filter(mark=='H3K27me3') %>% 
    left_join(H3K27me3_da_isect, by=c('region'='isect', 'group')) %>% 
    inner_join(K27_kmeans_clust, by=c('region'='isect')) %>% 
    inner_join(peak_clusts) %>% 
    mutate(
        p_clip=ifelse(is.na(p_clip), 0, p_clip),
        group=factor(group, levels=c('nepi', 'ctx', 'dien', 'nt', 'retina')),
        region=factor(region, levels=peak_order),
        cluster=factor(cluster, levels=c(4,10,8,5,6,2,1,7,9,3)),
        kmeans_clust=factor(kmeans_clust, levels=c(6,5,7,4,2,1,8,10,9,3))
    ) %>% 
    filter(!is.na(group))


p1 <- ggplot(H3K27ac_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.05, scale=1.4, fill=modality_colors['H3K27ac']) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    scale_fill_gradientn(colors=pals::brewer.ylgn(100)[0:80], limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) +
    ggtitle('H3K27ac')
p1

p2 <- ggplot(H3K27me3_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.05, scale=1.4, fill=modality_colors['H3K27me3']) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    scale_fill_gradientn(colors=pals::brewer.blues(100)[0:80], limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    no_y_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) +
    ggtitle('H3K27me3')
p2
(p1 | p2) + plot_layout(guides = 'collect')

ggsave('plots/paper/fig2/fig2_H3K27me3_H3K27ac_switches_DA_peaks_tracks_pval_ridges.pdf', width=8, height=12, bg='white')




## Color by mark
p1 <- ggplot(H3K27ac_plot_df, aes(poscut, gene_peak, height=coverage, fill=mark)) +
    geom_ridgeline(color='white', size=0.05, scale=1.4) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    scale_fill_manual(values=mark_colors) +
    article_text() +
    no_x_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) +
    ggtitle('H3K27ac')

p1
p2 <- ggplot(H3K27me3_plot_df, aes(poscut, gene_peak, height=coverage, fill=mark)) +
    geom_ridgeline(color='white', size=0.05, scale=1.4) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    scale_fill_manual(values=mark_colors) +
    article_text() +
    no_x_text() +
    no_y_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) +
    ggtitle('H3K27me3')
p2
(p1 | p2) + plot_layout(guides = 'collect')
ggsave('plots/paper/fig2/fig2_H3K27me3_H3K27ac_switches_DA_peaks_tracks_const_ridges.pdf', width=8, height=12, bg='white')




#### Heatmap for RNA expression ####
rna$lineage_coarse <- case_when(
    rna$lineage != 'other' ~ rna$lineage,
    T ~ rna$celltype_jf
)
rna <- Pando::aggregate_assay(rna, group_name='lineage_coarse', assay='RNA')

gene_names_use <- unique(peak_genes$gene_name) %>% intersect(rownames(rna))

rna_cluster <- rna@assays$RNA@misc$summary$lineage_coarse[, gene_names_use]
rna_cluster_expr <- rna_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='rna_expr')

rna_plot_df <- rna_cluster_expr %>% 
    inner_join(rename(peak_genes, 'isect'='feature'), by=c('feature'='gene_name')) %>% 
    inner_join(K27_kmeans_clust) %>% 
    group_by(feature) %>% 
    mutate(feature_num=as.numeric(factor(isect))) %>% 
    group_by(isect) %>%
    mutate(
        gene_peak=paste0(feature, '_', feature_num),
        isect=factor(isect, levels=isect_peak_order),
        cluster=factor(cluster, levels=c(4,10,8,5,6,2,1,7,9,3)),
        group=factor(group, levels=c('nect', 'ctx', 'dien', 'nt', 'retina'))
    )  %>% 
    filter(!is.na(group)) %>% 
    arrange(isect) %>% 
    mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak))) %>% 
    group_by(gene_peak) %>% 
    mutate(expr01=scale01(rna_expr))


ggplot(rna_plot_df, aes(group, gene_peak, fill=expr01)) +
    geom_tile() +
    scale_fill_gradientn(colors=orrd(1)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(cluster~., scales='free_y', space='free') +
    labs(y='Gene near peak', x='Group', fill='Relative\nexpression', title='RNA') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() 

ggsave('plots/paper/fig2/fig2_H3K27me3_H3K27ac_poised_DA_peaks_RNA_expr_heatmap.pdf', width=3, height=12, bg='white')





###### Bivalent: H3K4me3 vs H3K27me3 ####
# Get intersects 
H3K27me3_ranges <- StringToGRanges(rownames(marks$H3K27me3))
H3K27me3_ranges_extended <- H3K27me3_ranges %>% 
    Extend(upstream=2000, downstream=2000)

H3K4me3_ranges <- StringToGRanges(rownames(marks$H3K4me3))
H3K4me3_ranges_extended <- H3K4me3_ranges %>% 
    Extend(upstream=2000, downstream=2000)

all_intersect <- IRanges::intersect(H3K27me3_ranges_extended, H3K4me3_ranges_extended) 

H3K27me3_olaps <- findOverlaps(H3K27me3_ranges_extended, all_intersect)
H3K27me3_matches <- tibble(
    isect = GRangesToString(all_intersect)[subjectHits(H3K27me3_olaps)],
    H3K27me3 = rownames(marks$H3K27me3)[queryHits(H3K27me3_olaps)]
)

H3K4me3_olaps <- findOverlaps(H3K4me3_ranges_extended, all_intersect)
H3K4me3_matches <- tibble(
    isect = GRangesToString(all_intersect)[subjectHits(H3K4me3_olaps)],
    H3K4me3 = rownames(marks$H3K4me3)[queryHits(H3K4me3_olaps)]
)

Kme_isect_matches <- inner_join(H3K27me3_matches, H3K4me3_matches) 
Kme_isect_matches %>% write_tsv('data/intersect/H3Kme_bivalent_marks_intersect_matches.tsv')


#### Filter DA peaks by intersections ####
Kme_isect_matches <- read_tsv('data_/intersect/H3Kme_bivalent_marks_intersect_matches.tsv')
# all_lineage_da <- read_tsv('data/results/diff_expression/all_marks_lineage_DA.tsv')
# all_DA_df <- bind_rows('H3K4me3'=H3K4me3_da, 'H3K27me3'=H3K27me3_da, 'H3K4me3'=H3K4me3_da, .id='mark') %>% 
#     filter(group!='psc', feature%in%all_lineage_da$feature)
all_DA_df <- read_tsv('data_/results/diff_expression/all_marks_lineage_DA.tsv')


# H3K4me3
H3K4me3_da <- all_DA_df %>% filter(mark=='H3K4me3')
H3K4me3_isect_peaks <- intersect(H3K4me3_da$feature, Kme_isect_matches$H3K4me3)

H3K4me3_da_isect <- H3K4me3_da %>%
    dplyr::filter(feature%in%H3K4me3_isect_peaks) %>% 
    mutate(mark='H3K4me3')

# H3K27me3
H3K27me3_da <- all_DA_df %>% filter(mark=='H3K27me3')
H3K27me3_isect_peaks <- intersect(H3K27me3_da$feature, Kme_isect_matches$H3K27me3)

H3K27me3_da_isect <- H3K27me3_da %>%
    dplyr::filter(feature%in%H3K27me3_isect_peaks) %>% 
    mutate(mark='H3K27me3')

# Combine
Kme_da_isect <- Kme_isect_matches %>% 
    inner_join(H3K4me3_da_isect, by=c('H3K4me3'='feature')) %>% 
    inner_join(H3K27me3_da_isect, by=c('H3K27me3'='feature', 'group'), suffix=c('_H3K4me3', '_H3K27me3')) %>% 
    group_by(group, isect) %>% 
    filter(any(sign(p_clip_H3K4me3)!=sign(p_clip_H3K27me3)))

Kme_da_isect %>% write_tsv('data/intersect/H3Kme_bivalent_marks_DA_intersect.tsv')

Kme_da_rbind <- inner_join(Kme_isect_matches, H3K4me3_da_isect, by=c('H3K4me3'='feature')) %>% 
    bind_rows(inner_join(Kme_isect_matches, H3K27me3_da_isect, by=c('H3K27me3'='feature'))) %>% 
    filter(isect%in%Kme_da_isect$isect)

Kme_da_mat <- Kme_da_rbind %>%
    filter(group!='nect') %>% 
    mutate(mark_group=paste0(mark, '_', group)) %>% 
    dplyr::select(mark_group, isect, p_clip) %>%
    distinct() %>% 
    pivot_wider(names_from = mark_group, values_from = p_clip, values_fill = 0, values_fn=mean) %>%
    column_to_rownames('isect') %>% as.matrix()

isect_peak_order <- Kme_da_mat %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

set.seed(111)
Kme_kmeans_clust <- kmeans(Kme_da_mat, centers=10)$cluster %>% 
    enframe('isect', 'cluster')

Kme_da_isect_ranges <- Kme_da_isect$isect %>% 
    StringToGRanges()
Kme_da_isect_ranges %>% export.bed('data/intersect/H3Kme_bivalent_marks_DA_intersect.bed')


#### Plot enrichment p-value heatmap ####
K4me_max_df <- Kme_da_isect %>% 
    group_by(isect) %>% 
    filter(abs(signed_p_H3K4me3)==max(abs(signed_p_H3K4me3))) %>% 
    filter(row_number()==1) %>% 
    arrange(group, p_clip_H3K4me3) %>% 
    select(isect, 'max_group'=group) 

plot_df <- Kme_da_rbind %>% 
    inner_join(K4me_max_df) %>% 
    inner_join(Kme_kmeans_clust) %>% 
    # mutate(isect=factor(isect, levels=unique(K27ac_max_df$isect)))
    mutate(
        isect=factor(isect, levels=isect_peak_order),
        cluster=factor(cluster, levels=c(3,10,2,8,5,4,1,9,6,7)),
        max_group=factor(max_group, levels=c('nect', 'ctx', 'dien', 'nt', 'retina')),
        group=factor(group, levels=c('nect', 'ctx', 'dien', 'nt', 'retina'))
    )

p1 <- ggplot(plot_df, aes(group, isect, fill=p_clip)) +
    geom_tile() +
    facet_grid(cluster~mark, scales='free_y', space='free') +
    scale_fill_gradientn(colors=rev(pals::brewer.rdbu(100)), limits=c(-50, 50)) +
    no_y_text() +
    rotate_x_text(40) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Peak', x='Group', fill='Signed log10(FDR)') +
    no_legend()

p2 <- ggplot(plot_df, aes(mark, isect, fill=p_clip)) +
    geom_tile() +
    facet_grid(cluster~group, scales='free_y', space='free') +
    scale_fill_gradientn(colors=rev(pals::brewer.rdbu(100)), limits=c(-50, 50)) +
    no_y_text() +
    rotate_x_text(40) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Peak', x='Group', fill='Signed log10(FDR)')

p1 | p2

ggsave('plots/paper/fig2/fig2_H3K27me3_H3K4me3_bivalent_DA_peaks_pval_heatmap.pdf', width=6, height=10)


plot_df %>% write_tsv('data_/CT/bivalent_peaks.tsv')




## Color by expression
H3K27ac_expr_df <- expr_plot_df %>% 
    filter(mark=='H3K27ac') %>% 
    select('region'=isect, 'expr'=H3K27ac_expr, group, H3K27ac) %>% 
    group_by(region, group) %>% 
    summarize(expr=mean(expr)) %>% 
    group_by(region) %>% mutate(expr01=scale01(expr))

H3K27ac_expr_plot_df <- H3K27ac_plot_df %>% 
    mutate(group=fct_recode(group, 'nect'='nepi')) %>% 
    inner_join(H3K27ac_expr_df)


H3K27me3_expr_df <- expr_plot_df %>% 
    filter(mark=='H3K27me3') %>% 
    select('region'=isect, 'expr'=H3K27me3_expr, group, H3K27me3) %>% 
    group_by(region, group) %>% 
    summarize(expr=mean(expr)) %>% 
    group_by(region) %>% mutate(expr01=scale01(expr))

H3K27me3_expr_plot_df <- H3K27me3_plot_df %>% 
    mutate(group=fct_recode(group, 'nect'='nepi')) %>% 
    inner_join(H3K27me3_expr_df)



#### Peak expression heatmap ####
H3K4me3_cluster <- H3K4me3_use@assays$peaks@misc$summary$lineage_coarse[, unique(Kme_da_rbind$H3K4me3)]
H3K4me3_cluster_expr <- H3K4me3_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='H3K4me3_expr')

H3K27me3_cluster <- H3K27me3_use@assays$peaks@misc$summary$lineage_coarse[, unique(Kme_da_rbind$H3K27me3)]
H3K27me3_cluster_expr <- H3K27me3_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='H3K27me3_expr')

expr_plot_df <- Kme_da_rbind %>% 
    select(-(group:p_clip)) %>% 
    distinct() %>% 
    inner_join(H3K4me3_cluster_expr, by=c('H3K4me3'='feature')) %>% 
    inner_join(H3K27me3_cluster_expr, by=c('H3K27me3'='feature', 'group')) %>% 
    filter(group%in%c('nect', 'ctx', 'dien', 'nt', 'retina')) %>% 
    inner_join(K4me_max_df) %>% 
    inner_join(Kme_kmeans_clust) %>% 
    mutate(
        isect=factor(isect, levels=isect_peak_order),
        cluster=factor(cluster, levels=c(3,10,2,8,5,4,1,9,6,7)),
        group=factor(group, levels=group_order)
    ) %>% group_by(isect, mark, group, max_group, cluster) %>% 
    dplyr::summarize(H3K4me3_expr=mean(H3K4me3_expr), H3K27me3_expr=mean(H3K27me3_expr))

p1 <- ggplot(filter(expr_plot_df, mark=='H3K4me3'), aes(group, isect, fill=H3K4me3_expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=rdpu(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(cluster~., scales='free_y', space='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K4me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() 

p1

p2 <- ggplot(filter(expr_plot_df, mark=='H3K27me3'), aes(group, isect, fill=H3K27me3_expr)) +
    geom_tile() +
    scale_fill_gradientn(colors=blues(1.5)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(cluster~., scales='free_y', space='free') +
    labs(y='Peak', x='Group', fill='TF-IDF-normalized\nexpression', title='H3K27me3') +
    article_text() +
    rotate_x_text(40) +
    no_y_text()

p2

(p1 | p2) & theme(
    panel.spacing = unit(0.1, 'cm'),
    plot.margin = unit(rep(0.1,4), 'cm')
)
ggsave('plots/paper/fig2/fig2_H3K27me3_H3K4me3_bivalent_DA_peaks_expr_heatmap.pdf', width=8, height=10)


expr_plot_df %>% write_tsv('data_/CT/bivalent_peaks_expr.tsv')


#### Get regions and coverage for ridge plots ####
regions_use <- unique(Kme_da_rbind$isect) %>% 
    StringToGRanges() %>% 
    # resize(fix='center', width = 1) %>% 
    Extend(upstream=1000, downstream=1000)

peak_genes <- Signac::ClosestFeature(marks$H3K4me3, regions_use) %>% 
    as_tibble() %>% mutate(feature=unique(Kme_da_rbind$isect)) %>% 
    distinct(feature, gene_name)

marks$H3K4me3$lin_stage <- case_when(
    marks$H3K4me3$lineage != 'other' ~ marks$H3K4me3$lineage,
    marks$H3K4me3$stage != 'mid' ~ 'nepi',
    T ~ 'other'
)

marks$H3K27me3$lin_stage <- case_when(
    marks$H3K27me3$lineage != 'other' ~ marks$H3K27me3$lineage,
    marks$H3K27me3$stage != 'mid' ~ 'nepi',
    T ~ 'other'
)

H3K4me3_lin_nepi <- subset(marks$H3K4me3, lin_stage!='other')
H3K27me3_lin_nepi <- subset(marks$H3K27me3, lin_stage!='other')

names(regions_use) <- unique(Kme_da_rbind$isect)

H3K4me3_cov_df <- get_coverage(H3K4me3_lin_nepi, regions=regions_use, group_by='lin_stage')
H3K27me3_cov_df <- get_coverage(H3K27me3_lin_nepi, regions=regions_use, group_by='lin_stage')


Kme_plot_df <- bind_rows('H3K4me3'=H3K4me3_cov_df,'H3K27me3'=H3K27me3_cov_df, .id='mark') %>% 
    inner_join(peak_genes, by=c('region'='feature')) %>% 
    inner_join(K4me_max_df, by=c('region'='isect')) %>% 
    group_by(gene_name) %>% 
    mutate(feature_num=as.numeric(factor(region))) %>% 
    group_by(region, mark) %>%
    mutate(
        gene_peak=paste0(gene_name, '_', feature_num),
        coverage=coverage/max(coverage), 
        region=factor(region, levels=isect_peak_order),
        group=factor(group, levels=c('nepi', 'ctx', 'dien', 'nt', 'retina'))
    )  %>% 
    arrange(region) %>% 
    mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak)))


Kme_cov_mean_mat <- Kme_plot_df %>% 
    mutate(group_mark = paste0(group, mark)) %>% 
    group_by(group_mark, region) %>% 
    summarize(coverage=mean(coverage)) %>%
    pivot_wider(names_from=group_mark, values_from=coverage) %>% 
    column_to_rownames('region') %>% as.matrix()

Kme_cov_mat <- Kme_plot_df %>% 
    mutate(group_pos=paste0(mark, group, poscut)) %>% 
    ungroup() %>% 
    select(region, group_pos, coverage) %>% 
    pivot_wider(names_from = group_pos, values_from = coverage, values_fill=0) %>% 
    column_to_rownames('region') %>% as.matrix()

peak_order <- Kme_cov_mat %>% dist() %>% hclust() %>% {.$label[.$order]}
# peak_order <- Kme_cov_mean_mat %>% dist() %>% hclust() %>% {.$labels[.$order]}

set.seed(111)
peak_clusts <- Kme_cov_mean_mat %>% kmeans(centers = 10) %>% {.$cluster} %>% enframe('region', 'kmeans_clust')

H3K4me3_da_isect <- Kme_da_rbind %>% filter(mark=='H3K4me3') %>% 
    distinct(isect, group, p_clip) %>% 
    group_by(isect, group) %>% 
    summarise(p_clip=mean(p_clip))

H3K27me3_da_isect <- Kme_da_rbind %>% filter(mark=='H3K27me3') %>% 
    distinct(isect, group, p_clip) %>% 
    group_by(isect, group) %>% 
    summarise(p_clip=mean(p_clip))


H3K4me3_plot_df <- Kme_plot_df %>% 
    filter(mark=='H3K4me3') %>% 
    left_join(H3K4me3_da_isect, by=c('region'='isect', 'group')) %>% 
    inner_join(Kme_kmeans_clust, by=c('region'='isect')) %>% 
    inner_join(peak_clusts) %>% 
    mutate(
        p_clip=ifelse(is.na(p_clip), 0, p_clip),
        group=factor(group, levels=c('nepi', 'ctx', 'dien', 'nt', 'retina')),
        region=factor(region, levels=peak_order),
        cluster=factor(cluster, levels=c(3,10,2,8,5,4,1,9,6,7)),
        kmeans_clust=factor(kmeans_clust, levels=c(6,5,7,4,2,1,8,10,9,3))
    ) %>% 
    filter(!is.na(group))

H3K27me3_plot_df <- Kme_plot_df %>% 
    filter(mark=='H3K27me3') %>% 
    left_join(H3K27me3_da_isect, by=c('region'='isect', 'group')) %>% 
    inner_join(Kme_kmeans_clust, by=c('region'='isect')) %>% 
    inner_join(peak_clusts) %>% 
    mutate(
        p_clip=ifelse(is.na(p_clip), 0, p_clip),
        group=factor(group, levels=c('nepi', 'ctx', 'dien', 'nt', 'retina')),
        region=factor(region, levels=peak_order),
        cluster=factor(cluster, levels=c(3,10,2,8,5,4,1,9,6,7)),
        kmeans_clust=factor(kmeans_clust, levels=c(6,5,7,4,2,1,8,10,9,3))
    ) %>% 
    filter(!is.na(group))


p1 <- ggplot(H3K4me3_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.4) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    scale_fill_gradientn(colors=pals::brewer.rdpu(100)[0:70], limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) +
    ggtitle('H3K4me3')
p1

p2 <- ggplot(H3K27me3_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.4) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    scale_fill_gradientn(colors=pals::brewer.blues(100)[0:70], limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    no_y_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) +
    ggtitle('H3K27me3')
p2
(p1 | p2) + plot_layout(guides = 'collect')

ggsave('plots/paper/fig2/fig2_H3K27me3_H3K4me3_bivalent_DA_peaks_tracks_pval_ridges.pdf', width=8, height=18, bg='white')




## Color by expression
H3K4me3_expr_df <- expr_plot_df %>% 
    filter(mark=='H3K4me3') %>% 
    select('region'=isect, 'expr'=H3K4me3_expr, group, H3K4me3) %>% 
    group_by(region, group) %>% 
    summarize(expr=mean(expr)) %>% 
    group_by(region) %>% mutate(expr01=scale01(expr))

H3K4me3_expr_plot_df <- H3K4me3_plot_df %>% 
    mutate(group=fct_recode(group, 'nect'='nepi')) %>% 
    inner_join(H3K4me3_expr_df)


H3K27me3_expr_df <- expr_plot_df %>% 
    filter(mark=='H3K27me3') %>% 
    select('region'=isect, 'expr'=H3K27me3_expr, group, H3K27me3) %>% 
    group_by(region, group) %>% 
    summarize(expr=mean(expr)) %>% 
    group_by(region) %>% mutate(expr01=scale01(expr))

H3K27me3_expr_plot_df <- H3K27me3_plot_df %>% 
    mutate(group=fct_recode(group, 'nect'='nepi')) %>% 
    inner_join(H3K27me3_expr_df)



p1 <- ggplot(H3K4me3_expr_plot_df, aes(poscut, gene_peak, height=coverage, fill=expr01)) +
    geom_ridgeline(color='white', size=0.1, scale=1.4) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    scale_fill_gradientn(colors=pals::brewer.rdpu(100)[0:70]) +
    article_text() +
    no_x_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) +
    ggtitle('H3K4me3')
p1

p2 <- ggplot(H3K27me3_expr_plot_df, aes(poscut, gene_peak, height=coverage, fill=expr01)) +
    geom_ridgeline(color='white', size=0.1, scale=1.4) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    scale_fill_gradientn(colors=pals::brewer.blues(100)[0:70]) +
    article_text() +
    no_x_text() +
    no_y_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) +
    ggtitle('H3K27me3')
p2
(p1 | p2) + plot_layout(guides = 'collect')

ggsave('plots/paper/fig2/fig2_H3K27me3_H3K4me3_bivalent_DA_peaks_tracks_expr01_ridges.pdf', width=8, height=18, bg='white')


## Color by mark
p1 <- ggplot(H3K4me3_plot_df, aes(poscut, gene_peak, height=coverage)) +
    geom_ridgeline(color='white', size=0.1, scale=1.4, fill=mark_colors['H3K4me3']) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    # scale_color_manual(values=mark_colors) +
    article_text() +
    no_x_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) +
    ggtitle('H3K4me3')
p1

p2 <- ggplot(H3K27me3_plot_df, aes(poscut, gene_peak, height=coverage)) +
    geom_ridgeline(color='white', size=0.1, scale=1.4, fill=mark_colors['H3K4me3']) +
    facet_grid(cluster~group, scales='free_y', space='free_y') +
    # scale_color_manual(values=mark_colors) +
    article_text() +
    no_x_text() +
    no_y_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) +
    ggtitle('H3K27me3')
p2
(p1 | p2) + plot_layout(guides = 'collect')
ggsave('plots/paper/fig2/fig2_H3K27me3_H3K4me3_bivalent_DA_peaks_tracks_const_ridges.pdf', width=8, height=18, bg='white')



#### RNA expresion heatmaps ####
gene_names_use <- unique(peak_genes$gene_name) %>% intersect(rownames(rna))

rna_cluster <- rna@assays$RNA@misc$summary$lineage_coarse[, gene_names_use]
rna_cluster_expr <- rna_cluster %>% 
    as_tibble(rownames='group') %>% pivot_longer(!group, names_to='feature', values_to='rna_expr')

rna_plot_df <- rna_cluster_expr %>% 
    inner_join(rename(peak_genes, 'isect'='feature'), by=c('feature'='gene_name')) %>% 
    inner_join(Kme_kmeans_clust) %>% 
    group_by(feature) %>% 
    mutate(feature_num=as.numeric(factor(isect))) %>% 
    group_by(isect) %>%
    mutate(
        gene_peak=paste0(feature, '_', feature_num),
        isect=factor(isect, levels=isect_peak_order),
        cluster=factor(cluster, levels=c(3,10,2,8,5,4,1,9,6,7)),
        group=factor(group, levels=c('nect', 'ctx', 'dien', 'nt', 'retina'))
    )  %>% 
    filter(!is.na(group)) %>% 
    arrange(isect) %>% 
    mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak))) %>% 
    group_by(gene_peak) %>% 
    mutate(expr01=scale01(rna_expr))



ggplot(rna_plot_df, aes(group, gene_peak, fill=expr01)) +
    geom_tile() +
    scale_fill_gradientn(colors=orrd(1)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(cluster~., scales='free_y', space='free') +
    labs(y='Gene near peak', x='Group', fill='Relative\nexpression', title='RNA') +
    article_text() +
    rotate_x_text(40) +
    no_y_text() 

ggsave('plots/paper/fig2/fig2_H3K27me3_H3K4me3_bivalent_DA_peaks_RNA_expr_heatmap.pdf', width=3, height=18, bg='white')



###### Promoters: H3K4me3 vs H3K27ac ####
# Get intersects 
H3K27ac_ranges <- StringToGRanges(rownames(marks$H3K27ac))
H3K27ac_ranges_extended <- H3K27ac_ranges %>% 
    Extend(upstream=2000, downstream=2000)

H3K4me3_ranges <- StringToGRanges(rownames(marks$H3K4me3))
H3K4me3_ranges_extended <- H3K4me3_ranges %>% 
    Extend(upstream=2000, downstream=2000)

all_intersect <- IRanges::intersect(H3K27ac_ranges_extended, H3K4me3_ranges_extended) 

H3K27ac_olaps <- findOverlaps(H3K27ac_ranges_extended, all_intersect)
H3K27ac_matches <- tibble(
    isect = GRangesToString(all_intersect)[subjectHits(H3K27ac_olaps)],
    H3K27ac = rownames(marks$H3K27ac)[queryHits(H3K27ac_olaps)]
)

H3K4me3_olaps <- findOverlaps(H3K4me3_ranges_extended, all_intersect)
H3K4me3_matches <- tibble(
    isect = GRangesToString(all_intersect)[subjectHits(H3K4me3_olaps)],
    H3K4me3 = rownames(marks$H3K4me3)[queryHits(H3K4me3_olaps)]
)

prom_isect_matches <- inner_join(H3K27ac_matches, H3K4me3_matches) 
prom_isect_matches %>% write_tsv('data/intersect/H3K27K4_promoters_marks_intersect_matches.tsv')



#### Intersect bivalent and switches for all lineages ####
Kme_isect_matches <- read_tsv('data/intersect/H3Kme_bivalent_marks_intersect_matches.tsv')
K27_isect_matches <- read_tsv('data/intersect/H3K27_switches_marks_intersect_matches.tsv')
all_DA_df <- read_tsv('data/results/diff_expression/all_marks_lineage_DA.tsv')


#### Check peaks switch K27 into the branches ####
H3K27ac_use[['peaks_bin']] <- CreateAssayObject((H3K27ac_use[['peaks']]@data > 0)*1)
H3K27me3_use[['peaks_bin']] <- CreateAssayObject((H3K27me3_use[['peaks']]@data > 0)*1)

H3K27ac_use <- Pando::aggregate_assay(H3K27ac_use, assay='peaks_bin', group_name='clusters')
H3K27me3_use <- Pando::aggregate_assay(H3K27me3_use, assay='peaks_bin', group_name='clusters')

H3K27ac_cluster <- H3K27ac_use@assays$peaks_bin@misc$summary$clusters
H3K27ac_nepi_clusters <- unique(H3K27ac_use$clusters[H3K27ac_use$lineage_coarse=='nect'])
H3K27ac_nepi_detection <- H3K27ac_cluster[H3K27ac_nepi_clusters, ] %>% colMaxs()
H3K27ac_nepi_peaks <- colnames(H3K27ac_cluster)[H3K27ac_nepi_detection > 0.5]
H3K27ac_nepi_det <- tibble(
    feature=colnames(H3K27ac_cluster),
    detect=H3K27ac_nepi_detection
)

H3K27me3_cluster <- H3K27me3_use@assays$peaks_bin@misc$summary$clusters
H3K27me3_nepi_clusters <- unique(H3K27me3_use$clusters[H3K27me3_use$lineage_coarse=='nect'])
H3K27me3_nepi_detection <- H3K27me3_cluster[H3K27me3_nepi_clusters, ] %>% colMaxs()
H3K27me3_nepi_peaks <- colnames(H3K27me3_cluster)[H3K27me3_nepi_detection > 0.5]
H3K27me3_nepi_det <- tibble(
    feature=colnames(H3K27me3_cluster),
    detect=H3K27me3_nepi_detection
)


# H3K27ac
H3K27ac_da <- all_DA_df %>% filter(mark=='H3K27ac')
H3K27ac_isect_peaks <- intersect(H3K27ac_da$feature, K27_isect_matches$H3K27ac)

H3K27ac_da_isect <- H3K27ac_da %>%
    dplyr::filter(feature%in%H3K27ac_isect_peaks) %>% 
    mutate(mark='H3K27ac')

# H3K27me3
H3K27me3_da <- all_DA_df %>% filter(mark=='H3K27me3')
H3K27me3_isect_peaks <- intersect(H3K27me3_da$feature, K27_isect_matches$H3K27me3)

H3K27me3_da_isect <- H3K27me3_da %>%
    dplyr::filter(feature%in%H3K27me3_isect_peaks) %>% 
    mutate(mark='H3K27me3')


# Combine
K27_da_isect <- K27_isect_matches %>%
    inner_join(H3K27ac_da_isect, by=c('H3K27ac'='feature')) %>%
    inner_join(H3K27me3_da_isect, by=c('H3K27me3'='feature', 'group'), suffix=c('_H3K27ac', '_H3K27me3')) %>%
    group_by(group, isect) %>%
    filter(any(sign(p_clip_H3K27ac)!=sign(p_clip_H3K27me3)))

nepi_status <- K27_isect_matches %>%
    inner_join(H3K27ac_nepi_det, by=c('H3K27ac'='feature')) %>%
    inner_join(H3K27me3_nepi_det, by=c('H3K27me3'='feature'), suffix=c('_H3K27ac', '_H3K27me3')) %>%
    filter(detect_H3K27ac > 0.05 | detect_H3K27me3 > 0.05) %>% 
    mutate(nepi_status=case_when(
        (detect_H3K27ac/detect_H3K27me3) > 2 ~ 'H3K27ac',
        (detect_H3K27me3/detect_H3K27ac) > 2 ~ 'H3K27me3',
        detect_H3K27ac > 0.05 & detect_H3K27me3 > 0.05 ~ 'switching'
    )) %>% 
    filter(!is.na(nepi_status))



#### Split by mark status 
K27_da_isect_bival_counts <- K27_da_isect %>% 
    inner_join(nepi_status) %>% 
    mutate(mark_status = case_when(
        padj_H3K27ac<1e-4 & detect_self_H3K27ac>0.05 & coef_H3K27ac>0 & detect_self_H3K27me3<0.05 ~ 'H3K27ac',
        padj_H3K27me3<1e-4 & detect_self_H3K27me3>0.05 & coef_H3K27me3>0 & detect_self_H3K27ac<0.05 ~ 'H3K27me3',
        # detect_self_H3K27me3>0.05 & detect_self_H3K27ac>0.05 ~ 'switching',
        detect_self_H3K27me3<0.05 & detect_self_H3K27ac<0.05 ~ 'gone',
        T ~ 'other'
    )) %>% group_by(group, mark_status, nepi_status) %>% 
    # summarize(value=n()) %>% 
    summarize(value=length(unique(isect))) %>%
    mutate(
        switch_group=paste0(mark_status, '_', group),
        nepi='nepi'
    ) 

K27_allu <- gather_set_data(K27_da_isect_bival_counts, c(1,3)) %>% 
    mutate(
        # y=factor(y, levels=c('ctx', 'dien', 'nt', 'retina', 'H3K27ac', 'H3K27me3', 'none', 'switching')),
        x=factor(x, levels=c('nepi_status', 'group')),
        mark_status=factor(mark_status, levels=c('H3K27ac', 'H3K27me3', 'switching', 'gone', 'other'))
    ) %>% 
    filter(!mark_status%in%c('other', 'switching')) %>% filter(nepi_status!='switching')

modality_colors['switching'] <- '#29F4E6'
modality_colors['gone'] <- 'black'

ggplot(K27_allu, aes(x, id = id, split = y, value = value)) +
    geom_parallel_sets(aes(fill = mark_status), alpha = 0.3, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.1) +
    geom_parallel_sets_labels(colour = 'white') +
    scale_fill_manual(values=modality_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    no_x_text() +
    no_label() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave('plots/paper/fig2/fig2_nect_switches_alluvial.pdf', width=5, height=4)






#### Check how bivalent peaks resolve into the branches ####
H3K4me3_use[['peaks_bin']] <- CreateAssayObject((H3K4me3_use[['peaks']]@data > 0)*1)
H3K27me3_use[['peaks_bin']] <- CreateAssayObject((H3K27me3_use[['peaks']]@data > 0)*1)

H3K4me3_use <- Pando::aggregate_assay(H3K4me3_use, assay='peaks_bin', group_name='clusters')
H3K27me3_use <- Pando::aggregate_assay(H3K27me3_use, assay='peaks_bin', group_name='clusters')

H3K4me3_cluster <- H3K4me3_use@assays$peaks_bin@misc$summary$clusters[, unique(Kme_da_rbind$H3K4me3)]
H3K4me3_nepi_clusters <- unique(H3K4me3_use$clusters[H3K4me3_use$lineage_coarse=='nect'])
H3K4me3_nepi_detection <- H3K4me3_cluster[H3K4me3_nepi_clusters, ] %>% colMaxs()
H3K4me3_nepi_peaks <- colnames(H3K4me3_cluster)[H3K4me3_nepi_detection > 0.05]

H3K27me3_cluster <- H3K27me3_use@assays$peaks_bin@misc$summary$clusters[, unique(Kme_da_rbind$H3K27me3)]
H3K27me3_nepi_clusters <- unique(H3K27me3_use$clusters[H3K27me3_use$lineage_coarse=='nect'])
H3K27me3_nepi_detection <- H3K27me3_cluster[H3K27me3_nepi_clusters, ] %>% colMaxs()
H3K27me3_nepi_peaks <- colnames(H3K27me3_cluster)[H3K27me3_nepi_detection > 0.05]

H3K4me3_da_isect <- H3K4me3_da %>%
    dplyr::filter(feature%in%H3K4me3_isect_peaks) %>% 
    mutate(mark='H3K4me3')

# H3K27me3
H3K27me3_da <- all_DA_df %>% filter(mark=='H3K27me3')
H3K27me3_isect_peaks <- intersect(H3K27me3_da$feature, Kme_isect_matches$H3K27me3)

H3K27me3_da_isect <- H3K27me3_da %>%
    dplyr::filter(feature%in%H3K27me3_isect_peaks) %>% 
    mutate(mark='H3K27me3')

# Combine
Kme_da_isect <- Kme_isect_matches %>% 
    inner_join(H3K4me3_da_isect, by=c('H3K4me3'='feature')) %>% 
    inner_join(H3K27me3_da_isect, by=c('H3K27me3'='feature', 'group'), suffix=c('_H3K4me3', '_H3K27me3')) %>% 
    group_by(group, isect) %>% 
    filter(any(sign(p_clip_H3K4me3)!=sign(p_clip_H3K27me3)))

nepi_bival <- Kme_isect_matches %>% filter(H3K4me3%in%H3K4me3_nepi_peaks, H3K27me3%in%H3K27me3_nepi_peaks) %>% pull(isect) %>% unique()

Kme_da_isect$bival_in_nepi <- Kme_da_isect$isect %in% nepi_bival


#### Split by mark status 
Kme_da_isect_bival <- Kme_da_isect %>% 
    filter(bival_in_nepi) %>% 
    mutate(mark_status = case_when(
        padj_H3K4me3<1e-4 & detect_self_H3K4me3<0.05 & coef_H3K4me3<0 & detect_self_H3K27me3>0.05 ~ 'H3K27me3',
        padj_H3K27me3<1e-4 & detect_self_H3K27me3<0.05 & coef_H3K27me3<0 & detect_self_H3K4me3>0.05 ~ 'H3K4me3',
        detect_self_H3K27me3>0.05 & detect_self_H3K4me3>0.05 ~ 'bivalent',
        detect_self_H3K27me3<0.05 & detect_self_H3K4me3<0.05 ~ 'gone',
        T ~ 'other'
    )) 

Kme_da_isect_bival_counts <- Kme_da_isect_bival %>% group_by(group, mark_status, bival_in_nepi) %>% 
    # summarize(value=n()) %>% 
    summarize(value=length(unique(isect))) %>%
    mutate(
        switch_group=paste0(mark_status, '_', group),
        nepi='bivalent'
    ) 


Kme_allu <- gather_set_data(Kme_da_isect_bival_counts, c(1,6)) %>% 
    mutate(
        y=factor(y, levels=c('ctx', 'dien', 'nt', 'retina', 'none', 'bivalent')),
        x=factor(x, levels=c('nepi', 'group')),
        mark_status=factor(mark_status, levels=c('H3K4me3', 'H3K27me3', 'bivalent', 'gone', 'other'))
    ) %>% 
    filter(mark_status!='other')

modality_colors['bivalent'] <- '#9C6FF7'
modality_colors['gone'] <- 'black'

ggplot(Kme_allu, aes(x, id = id, split = y, value = value)) +
    geom_parallel_sets(aes(fill = mark_status), alpha = 0.3, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.1) +
    geom_parallel_sets_labels(colour = 'white') +
    scale_fill_manual(values=modality_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    no_x_text() +
    no_label() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave('plots/paper/fig2/fig2_nect_bivalent_alluvial.pdf', width=5, height=4)



#### Plot fraction of peaks with K27ac ####

H3K27ac_lineage_clusters <- unique(H3K27ac_use$clusters[H3K27ac_use$lineage_coarse%in%c('ctx', 'dien', 'nt', 'retina')])
H3K27ac_lineage_detect <- H3K27ac_cluster[H3K27ac_lineage_clusters, unique(K27_da_rbind$H3K27ac)] %>% 
    as_tibble(rownames='clusters') %>% 
    pivot_longer(!clusters, names_to='H3K27ac', values_to='detect') %>% 
    inner_join(distinct(as_tibble(marks$H3K27ac@meta.data), clusters, lineage)) %>% 
    group_by(H3K27ac, lineage) %>% 
    summarize(detect=max(detect))
    
K4me3_K27_matches <- StringToGRanges(Kme_da_isect_bival$H3K4me3) %>% 
    findOverlaps(StringToGRanges(H3K27ac_lineage_detect$H3K27ac))

K27me3_K27_matches <- StringToGRanges(Kme_da_isect_bival$H3K27me3) %>% 
    findOverlaps(StringToGRanges(H3K27ac_lineage_detect$H3K27ac))

Kme_K27_isects <- tibble(
    isect = Kme_da_isect_bival$isect[c(queryHits(K4me3_K27_matches), queryHits(K27me3_K27_matches))],
    H3K27ac = H3K27ac_lineage_detect$H3K27ac[c(subjectHits(K4me3_K27_matches), subjectHits(K27me3_K27_matches))]
) %>% distinct()

Kme_da_isect_K27 <- Kme_da_isect_bival %>% 
    left_join(Kme_K27_isects) %>% 
    left_join(H3K27ac_lineage_detect) %>% 
    mutate(has_K27ac=detect>0.05 & !is.na(detect)) %>% 
    filter(mark_status!='other', mark_status!='gone') %>% 
    distinct(isect, mark_status, has_K27ac)

ggplot(Kme_da_isect_K27, aes(mark_status, fill=has_K27ac)) +
    geom_bar(color='black') +
    facet_grid(~group) +
    scale_fill_manual(values=c('white', 'darkgrey'))

ggsave('plots/paper/fig2/fig2_nect_bivalent_alluvial_K27ac_bar.pdf', width=5, height=4)

ggplot(Kme_da_isect_K27, aes(mark_status, fill=has_K27ac)) +
    geom_bar(color='black', position='fill') +
    facet_grid(~group) +
    scale_fill_manual(values=c('white', 'darkgrey'))

ggsave('plots/paper/fig2/fig2_nect_bivalent_alluvial_K27ac_fill.pdf', width=5, height=4)





#### Plot peak class distribution over stages ####
Kme_isect_matches <- read_tsv('data/intersect/H3Kme_bivalent_marks_intersect_matches.tsv')
K27_isect_matches <- read_tsv('data/intersect/H3K27_switches_marks_intersect_matches.tsv')
prom_isect_matches <- read_tsv('data/intersect/H3K27K4_promoters_marks_intersect_matches.tsv')

H3K27ac_ctx <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_subs_dpt_srt.rds')
H3K4me3_ctx <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_subs_dpt_srt.rds')
H3K27me3_ctx <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_subs_dpt_srt.rds')

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

all_stage_detect <- prom_isect_matches %>% 
    full_join(Kme_isect_matches) %>% 
    inner_join(H3K27ac_ctx_stage_detect, by=c('H3K27ac'='peak')) %>%
    inner_join(H3K27me3_ctx_stage_detect, by=c('stage', 'H3K27me3'='peak'), suffix=c('_H3K27ac','_H3K27me3')) %>% 
    inner_join(H3K4me3_ctx_stage_detect, by=c('stage', 'H3K4me3'='peak'), suffix=c('','_H3K4me3')) 
    
H3K27ac_q <- H3K27ac_ctx_stage_detect$detect %>% quantile(0.75)
H3K27me3_q <- H3K27me3_ctx_stage_detect$detect %>% quantile(0.75)
H3K4me3_q <- H3K4me3_ctx_stage_detect$detect %>% quantile(0.75)

# H3K27ac_q <- 0.02
# H3K27me3_q <- 0.02
# H3K4me3_q <- 0.02

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

all_class_genes %>% write_tsv('data/trajectories/ctx/ctx_traject_ptbin25_peak_class.tsv')


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





#### GO enrichment ####
library(clusterProfiler)
library(org.Hs.eg.db)

feats_use <- all_class_genes$gene_name %>% unique()

all_features <- bitr(feats_use, fromType = 'SYMBOL', toType = c('ENSEMBL', 'ENTREZID'), OrgDb = org.Hs.eg.db) %>% 
    as_tibble()

genes_test <- setdiff(b1_bival, b25_bival)
gene_ids <- filter(all_features, SYMBOL%in%genes_test)

ggo <- enrichGO(
    gene = gene_ids$ENTREZID, 
    universe = all_features$ENTREZID,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = 'fdr',
    ont = 'BP',
    readable = T
)


#### Check for DIEN lineages ####

H3K27ac_nt <- read_rds('data/trajectories/nt/H3K27ac_nt_subs_dpt_lsi_regress_srt.rds')
H3K4me3_nt <- read_rds('data/trajectories/nt/H3K4me3_nt_subs_dpt_lsi_regress_srt.rds')
H3K27me3_nt <- read_rds('data/trajectories/nt/H3K27me3_nt_subs_dpt_lsi_regress_srt.rds')

H3K27ac_nt[['peaks_bin']] <- CreateAssayObject((H3K27ac_nt[['peaks']]@data > 0)*1)
H3K4me3_nt[['peaks_bin']] <- CreateAssayObject((H3K4me3_nt[['peaks']]@data > 0)*1)
H3K27me3_nt[['peaks_bin']] <- CreateAssayObject((H3K27me3_nt[['peaks']]@data > 0)*1)

H3K27ac_nt$nt_pt <- rank(H3K27ac_nt$nt_pt) / max(rank(H3K27ac_nt$nt_pt))
H3K4me3_nt$nt_pt <- rank(H3K4me3_nt$nt_pt) / max(rank(H3K4me3_nt$nt_pt))
H3K27me3_nt$nt_pt <- rank(H3K27me3_nt$nt_pt) / max(rank(H3K27me3_nt$nt_pt))

H3K27ac_nt$pt_bin2 <- as.numeric(cut(H3K27ac_nt$nt_pt, 10, labels=1:10))
H3K27ac_nt$pt_bin <- as.numeric(cut(H3K27ac_nt$nt_pt, 20, labels=1:20))
H3K27ac_nt$diff_stage <- case_when(
    H3K27ac_nt$pt_bin <= 15 ~ 'npc',
    T ~ 'neuron'
)
H3K4me3_nt$pt_bin2 <- as.numeric(cut(H3K4me3_nt$nt_pt, 10, labels=1:10))
H3K4me3_nt$pt_bin <- as.numeric(cut(H3K4me3_nt$nt_pt, 20, labels=1:20))
H3K4me3_nt$diff_stage <- case_when(
    H3K4me3_nt$pt_bin <= 15 ~ 'npc',
    T ~ 'neuron'
)
H3K27me3_nt$pt_bin2 <- as.numeric(cut(H3K27me3_nt$nt_pt, 10, labels=1:10))
H3K27me3_nt$pt_bin <- as.numeric(cut(H3K27me3_nt$nt_pt, 20, labels=1:20))
H3K27me3_nt$diff_stage <- case_when(
    H3K27me3_nt$pt_bin <= 15 ~ 'npc',
    T ~ 'neuron'
)

H3K27ac_nt <- Pando::aggregate_assay(H3K27ac_nt, assay='peaks_bin', group_name='stage')
H3K4me3_nt <- Pando::aggregate_assay(H3K4me3_nt, assay='peaks_bin', group_name='stage')
H3K27me3_nt <- Pando::aggregate_assay(H3K27me3_nt, assay='peaks_bin', group_name='stage')

H3K27ac_nt <- Pando::aggregate_assay(H3K27ac_nt, assay='peaks_bin', group_name='pt_bin2')
H3K4me3_nt <- Pando::aggregate_assay(H3K4me3_nt, assay='peaks_bin', group_name='pt_bin2')
H3K27me3_nt <- Pando::aggregate_assay(H3K27me3_nt, assay='peaks_bin', group_name='pt_bin2')

H3K27ac_nt_stage_detect <- H3K27ac_nt@assays$peaks_bin@misc$summary$pt_bin2 %>% 
    {.[, colMaxs(.)>0.05]} %>%
    as_tibble(rownames='stage') %>% 
    pivot_longer(!stage, names_to='peak', values_to='detect') 
H3K4me3_nt_stage_detect <- H3K4me3_nt@assays$peaks_bin@misc$summary$pt_bin2 %>% 
    {.[, colMaxs(.)>0.05]} %>%
    as_tibble(rownames='stage') %>% 
    pivot_longer(!stage, names_to='peak', values_to='detect') 
H3K27me3_nt_stage_detect <- H3K27me3_nt@assays$peaks_bin@misc$summary$pt_bin2 %>% 
    {.[, colMaxs(.)>0.05]} %>%
    as_tibble(rownames='stage') %>% 
    pivot_longer(!stage, names_to='peak', values_to='detect') 

all_stage_detect <- prom_isect_matches %>% 
    full_join(Kme_isect_matches) %>% 
    inner_join(H3K27ac_nt_stage_detect, by=c('H3K27ac'='peak')) %>% 
    inner_join(H3K27me3_nt_stage_detect, by=c('stage', 'H3K27me3'='peak'), suffix=c('_H3K27ac','_H3K27me3')) %>% 
    inner_join(H3K4me3_nt_stage_detect, by=c('stage', 'H3K4me3'='peak'), suffix=c('','_H3K4me3')) 

H3K27ac_q <- H3K27ac_nt_stage_detect$detect %>% quantile(0.5)
H3K27me3_q <- H3K27me3_nt_stage_detect$detect %>% quantile(0.5)
H3K4me3_q <- H3K4me3_nt_stage_detect$detect %>% quantile(0.5)

H3K27ac_q <- 0.05
H3K27me3_q <- 0.05
H3K4me3_q <- 0.05

all_stage_class <- all_stage_detect %>% 
    mutate(
        detect_H3K27ac=replace_na(detect_H3K27ac, 0),
        detect_H3K27me3=replace_na(detect_H3K27me3, 0),
        detect=replace_na(detect, 0)
    ) %>%
    mutate(peak_class=case_when(
        (detect_H3K27me3 > H3K27me3_q) & (detect > H3K4me3_q) ~ 'bivalent',
        (detect > H3K4me3_q) & (detect_H3K27ac > H3K27ac_q) ~ 'promoter',
        detect_H3K27me3 > H3K27me3_q ~ 'H3K27me3',
        detect_H3K27ac > H3K27ac_q ~ 'H3K27ac',
        detect > H3K4me3_q ~ 'H3K4me3',
        T ~ 'none'
    )) %>% filter(peak_class!='none') %>% 
    mutate(peak_class=factor(peak_class, levels=c('H3K27me3', 'bivalent', 'H3K4me3', 'promoter', 'H3K27ac'))) %>% 
    # mutate(stage=factor(stage, levels=c('psc', 'nepi', 'npc', 'neuron')))
    # mutate(stage=factor(stage, levels=c('EB', '15d', '35d', '60d', '128d', '4mo', '8mo')))
    # mutate(stage=factor(stage, levels=c('EB', 'mid', 'late', 'mo8')))
    mutate(stage=factor(stage, levels=c(1:20)))

modality_colors['promoter'] <- '#E58B4A'

ggplot(all_stage_class, aes(stage, fill=peak_class)) +
    geom_bar(position='fill', alpha=0.8) +
    scale_fill_manual(values=modality_colors) +
    article_text()
ggsave('plots/paper/fig3/fig3_nt_traject_age_peak_class_bar.pdf', width=6, height=5, units='cm')




#### Repression v Expression ####

marks$H3K27me3[['cRNA_bin']] <- CreateAssayObject((marks$H3K27me3[['cRNA']]@data > 0)*1)
marks$H3K27me3 <- aggregate_assay(marks$H3K27me3, assay='cRNA_bin', group_name='clusters')

marks$H3K27ac[['cRNA_bin']] <- CreateAssayObject((marks$H3K27ac[['cRNA']]@data > 0)*1)
marks$H3K27ac <- aggregate_assay(marks$H3K27ac, assay='cRNA_bin', group_name='clusters')

marks$H3K4me3[['cRNA_bin']] <- CreateAssayObject((marks$H3K4me3[['cRNA']]@data > 0)*1)
marks$H3K4me3 <- aggregate_assay(marks$H3K4me3, assay='cRNA_bin', group_name='clusters')

rna[['RNA_bin']] <- CreateAssayObject((rna[['RNA']]@data > 0)*1)
rna <- aggregate_assay(rna, assay='RNA_bin', group_name='clusters')

repressed_genes <- marks$H3K27me3[['cRNA_bin']]@misc$summary$clusters %>% colMaxs()
names(repressed_genes) <- rownames(marks$H3K27me3[['cRNA_bin']])

expressed_genes <- rna[['RNA_bin']]@misc$summary$clusters %>% colMaxs()
names(expressed_genes) <- rownames(rna)

repressed_peaks_df <- enframe(repressed_genes, 'gene', 'max_repr')
expressed_genes_df <- enframe(expressed_genes, 'gene', 'max_expr')
pc_genes <- gene_annot %>% as_tibble() %>% filter(gene_biotype=='protein_coding') %>% pull(gene_name) %>% unique()

repr_groups <- full_join(repressed_peaks_df, expressed_genes_df) %>% 
    mutate(
        max_expr=replace_na(max_expr, 0), 
        max_repr=replace_na(max_repr, 0),
        is_pc=gene%in%pc_genes,
        class=case_when(
            max_repr>0.05 & max_expr>0.05 ~ 'expressed & repressed',
            max_repr<=0.05 & max_expr>0.05 ~ 'expressed & not repressed',
            max_repr>0.05 & max_expr<=0.05 ~ 'not expressed & repressed',
            T ~ 'not expressed & not repressed',
    )) %>% filter(is_pc)

plot_df <- repr_groups
p1 <- ggplot(plot_df, aes(max_repr, max_expr, color=class)) +
    geom_point(size=0.2, alpha=0.5) +
    geom_hline(yintercept=0.05) +
    geom_vline(xintercept=0.05) +
    no_legend() +
    ggtitle('Max cluster expression vs repression')

p2 <- ggplot(plot_df, aes(class, fill=class)) +
    geom_bar() +
    rotate_x_text(40) +
    ggtitle('Counts of gene state across the dataset')

p1 + p2




#### GO enrichment ####
library(clusterProfiler)
library(org.Hs.eg.db)
go_enrich <- function(features, background, p_adjust_method='fdr', ontology='BP', ...){
    require(clusterProfiler)
    require(org.Hs.eg.db)
    
    all_features <- bitr(background, fromType='SYMBOL', toType=c('ENTREZID'), OrgDb=org.Hs.eg.db) %>% 
        as_tibble()
    
    gene_ids <- filter(all_features, SYMBOL%in%features)
    
    ggo <- enrichGO(
        gene = unique(gene_ids$ENTREZID), 
        universe = unique(all_features$ENTREZID),
        OrgDb = org.Hs.eg.db,
        readable = T,
        pAdjustMethod = p_adjust_method,
        ont = ontology,
        ...
    )
    res <- as_tibble(ggo@result) %>% 
        rename('count'=Count) %>% 
        mutate(
            fg_in = as.numeric(str_replace(GeneRatio, '(\\d+)/(\\d+)', '\\1')),
            fg_out = as.numeric(str_replace(GeneRatio, '(\\d+)/(\\d+)', '\\2')),
            bg_in = as.numeric(str_replace(BgRatio, '(\\d+)/(\\d+)', '\\1')),
            bg_out = as.numeric(str_replace(BgRatio, '(\\d+)/(\\d+)', '\\2')),
            odds_ratio = (fg_in/fg_out) / (bg_in/bg_out),
            log2_odds_ratio = log2(odds_ratio)
        )
    
    return(res)
}

gene_groups <- repr_groups %>% 
    group_by(class) %>% group_split() %>% 
    {names(.) <- map_chr(., ~.x$class[1]);.} %>% 
    map(~unique(.x$gene)) 

names(gene_groups) <- str_replace_all(names(gene_groups), ' ', '_') %>% str_replace_all('_&', '') 

gene_groups_go <- map_dfr(gene_groups, function(x){
    go_enrich(x, background=repr_groups$gene)
}, .id='gene_group')

gene_groups_go %>% write_tsv('data_/results/enrichment/gene_groups_go_enrich.tsv')




#### Motif enrich ####
#### Enrichment in merged regions ####
tf_motif_enrichment <- function(regions, motifs, background_regions, parallel=T){
    data(motif2tf)
    group_vec <- background_regions %in% regions
    map_par(colnames(motifs), function(mot){
        ctab <- table(motifs[background_regions, mot], group_vec)
        ftest <- try(fisher.test(ctab))
        if (any(class(ftest)=='try-error')){
            return(tibble())
        }
        return(tibble(
            motif = mot,
            pval = ftest$p.value,
            odds_ratio = ftest$estimate
        ))
    }, parallel=parallel) %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate(padj=p.adjust(pval, 'fdr')) %>%
        dplyr::mutate(logodds=log2(odds_ratio)) %>%
        dplyr::inner_join(motif2tf, multiple = "all") %>%
        return()
}

all_bg <- peak_meta$repr_region %>% unique() %>% StringToGRanges()
all_peaks_srt <- read_rds('data05/regions/all_regions_v2_max01_clust50_srt.rds')
all_motif_mat <- read_rds('data_/all_peaks/all_motifs.rds')

gene_ranges_use <- CollapseToLongestTranscript(gene_annot)

library(doParallel)
registerDoParallel(40)

gene_groups_tf_mot <- map_dfr(gene_groups, function(x){
    gene_ranges_use <- gene_ranges_use[gene_ranges_use$gene_name%in%x]
    gene2peaks <- find_peaks_near_genes(all_bg, gene_ranges_use, distance=4000, only_tss=T)
    peaks_use <- rownames(gene2peaks)[rowMaxs(gene2peaks)>0] 
    tf_enrich <- tf_motif_enrichment(peaks_use, all_motif_mat, GRangesToString(all_bg), parallel=T)
    return(tf_enrich)
}, .id='gene_group')
gene_groups_tf_mot %>% write_tsv('data_/results/enrichment/gene_groups_tf_promoter_enrich.tsv')
gene_groups_tf_mot <- read_tsv('data_/results/enrichment/gene_groups_tf_promoter_enrich.tsv')


gene_groups_tf_mot <- map_dfr(gene_groups, function(x){
    gene_ranges_use <- gene_ranges_use[gene_ranges_use$gene_name%in%x]
    gene2peaks <- find_peaks_near_genes(all_bg, gene_ranges_use, distance=4000)
    peaks_use <- rownames(gene2peaks)[rowMaxs(gene2peaks)>0] 
    tf_enrich <- tf_motif_enrichment(peaks_use, all_motif_mat, GRangesToString(all_bg), parallel=T)
    return(tf_enrich)
}, .id='gene_group')
gene_groups_tf_mot %>% write_tsv('data_/results/enrichment/gene_groups_tf_enrich.tsv')
gene_groups_tf_mot <- read_tsv('data_/results/enrichment/gene_groups_tf_enrich.tsv')



##### Graph abstraction #####
library(ggraph)
cluster_graph <- read_rds('data/noastro_RNA_cluster_graph.rds')

set.seed(9)
ggraph(cluster_graph, layout='fr') +
    geom_edge_link(alpha=0.1, width=0.3) +
    geom_node_point(aes(color=celltype), size=5) +
    scale_color_manual(values=pantone_celltype) +
    theme_void() +
    no_legend()

ggsave('plots/paper/fig2/fig2_cellrank_lineage_graph_fr.pdf', width=7, height=6)


ggraph(cluster_graph, x=-CSSUMAP_1, y=CSSUMAP_2) +
    geom_edge_link(alpha=0.1, width=0.3) +
    geom_node_point(aes(fill=celltype), size=6, shape=21, color='darkgrey') +
    scale_fill_manual(values=pantone_celltype) +
    theme_void() +
    no_legend()

ggsave('plots/paper/fig2/fig2_cellrank_lineage_graph_umap.pdf', width=6.5, height=5)



rna_meta <- rna@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(rna[['cssumap']]@cell.embeddings, rownames='cell'))

ggplot(rna_meta, aes(-CSSUMAP_1, CSSUMAP_2, color=celltype_jf)) +
    geom_point(alpha=0.1, size=0.5) +
    scale_color_manual(values=pantone_celltype) +
    theme_void() +
    no_legend()

ggsave('plots/paper/fig2/fig2_celltype_jf_umap.png', width=7, height=5)



#### Feature plots on graph abstraction ####
cluster_srt <- read_rds('data/noastro_RNA_marks_combined_clusters_srt.rds')

goi <- c('VSX2', 'LHX5', 'FOXG1', 'NFIX', 'WLS', 'NEUROD2')


# New colors with darker grey
gyylorrd2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='YlOrRd')[3:9]), bias=bias)(100))
}

gyylgn2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='YlGn')[3:9]), bias=bias)(100))
}

gybu2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='Blues')[3:9]), bias=bias)(100))
}

gyrdpu2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='RdPu')[3:9]), bias=bias)(100))
}


cluster_srt@active.assay <- 'RNA'
p_rna <- feature_plot(cluster_srt, features=goi, pt.size=2.5) &
    scale_color_gradientn(colors=gyylorrd2())

cluster_srt@active.assay <- 'H3K27ac_RNA'
p_K27ac <- feature_plot(cluster_srt, features=goi, pt.size=2.5) &
    scale_color_gradientn(colors=gyylgn2())

cluster_srt@active.assay <- 'H3K27me3_RNA'
p_K27me3 <- feature_plot(cluster_srt, features=goi, pt.size=2.5) &
    scale_color_gradientn(colors=gybu2())

cluster_srt@active.assay <- 'H3K4me3_RNA'
p_K4me3 <- feature_plot(cluster_srt, features=goi, pt.size=2.5) &
    scale_color_gradientn(colors=gyrdpu2())

p_K27ac | p_K27me3 | p_K4me3 | p_rna 
ggsave('plots/paper/fig2/rev_fig2_feature_cluster_umap.png', width=15, height=6)
ggsave('plots/paper/fig2/rev_fig2_feature_cluster_umap.pdf', width=15, height=6)


cluster_srt@active.assay <- 'RNA'
p_rna <- feature_plot(cluster_srt, features=goi, pt.size=2.5, reduction='fr') &
    scale_color_gradientn(colors=gyylorrd2())

cluster_srt@active.assay <- 'H3K27ac_RNA'
p_K27ac <- feature_plot(cluster_srt, features=goi, pt.size=2.5, reduction='fr') &
    scale_color_gradientn(colors=gyylgn2())

cluster_srt@active.assay <- 'H3K27me3_RNA'
p_K27me3 <- feature_plot(cluster_srt, features=goi, pt.size=2.5, reduction='fr') &
    scale_color_gradientn(colors=gybu2())

cluster_srt@active.assay <- 'H3K4me3_RNA'
p_K4me3 <- feature_plot(cluster_srt, features=goi, pt.size=2.5, reduction='fr') &
    scale_color_gradientn(colors=gyrdpu2())

p_K27ac | p_K27me3 | p_K4me3 | p_rna 
ggsave('plots/paper/fig2/rev_fig2_feature_cluster_fr.png', width=18, height=7)
ggsave('plots/paper/fig2/rev_fig2_feature_cluster_fr.pdf', width=18, height=7)


nn_graph_fr <- read_rds('data/noastro_RNA_cluster_graph.rds')

ggraph(nn_graph_fr, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=celltype), size=4) +
    scale_color_manual(values=pantone_celltype) +
    theme_void()
ggsave('plots/paper/fig2/rev_fig2_celltype_cluster_fr.png', width=7, height=6)
ggsave('plots/paper/fig2/rev_fig2_celltype_cluster_fr.pdf', width=7, height=6)



#### CTX pseudotime clusters ####
H3K27ac_ctx_subs <- read_rds('data/trajectories/ctx/H3K27ac_eb_ctx_subs_dpt_srt.rds')
H3K4me3_ctx_subs <- read_rds('data/trajectories/ctx/H3K4me3_eb_ctx_subs_dpt_srt.rds')
H3K27me3_ctx_subs <- read_rds('data/trajectories/ctx/H3K27me3_eb_ctx_subs_dpt_srt.rds')
rna_ctx_subs <- read_rds('data/trajectories/ctx/RNA_eb_ctx_subs_dpt_srt.rds')

mark_feats <- intersect(rownames(H3K27ac_ctx_subs[['cRNA']]), rownames(H3K27me3_ctx_subs[['cRNA']])) %>% intersect(rownames(H3K4me3_ctx_subs[['cRNA']]))

mod_colors <- c('H3K4me3'='#F87AA7', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')


#### Plot expression over pt bins ####
H3K27ac_clusters <- H3K27ac_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K27me3_clusters <- H3K27me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
H3K4me3_clusters <- H3K4me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), ]
rna_clusters <- rna_ctx_subs@assays$RNA@misc$summary$pt_bins[as.character(1:50), ]

genes_plot <- c(
    'POU5F1',
    'SOX2', 'VIM', 'GLI3', 'GRIA2',
    'STMN2', 'BCL11A', 'NEUROD6', 'NEUROD2', 
    'OTX1','OTX2','FOXB1','ZIC5','ZIC2','LHX5','WLS','MID1',
    'PRICKLE1','WWTR1','FZD5','FOXO1','COX9','SMG6','ANXA2'
) %>% intersect(mark_feats)

genes_plot <- c('INTU','EGR1','B3GAT2','POU3F3','DCLK2','SNCAIP','HOPX','PAG1',
                'TRIM9','LRRC3B','LHX2','NR2F1','POU3F2','ZEB1','MAPRE2','RCN1',
                'EFNB1','PAX6','EMX2','NR2F2','FEZF2','ARX') %>% intersect(mark_feats)

H3K27ac_expr <- H3K27ac_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

H3K27me3_expr <- H3K27me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

H3K4me3_expr <- H3K4me3_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

rna_expr <- rna_clusters[, genes_plot] %>% t() %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='pt_bins', values_to='expr')

plot_df <- bind_rows('H3K27ac'=H3K27ac_expr, 'H3K27me3'=H3K27me3_expr, 'H3K4me3'=H3K4me3_expr, 'RNA'=rna_expr, .id='modality') %>% 
    group_by(modality, gene) %>% 
    mutate(expr=scale01(expr))

p1 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality, group=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    # geom_point(size=0.2) +
    scale_color_manual(values=mod_colors) +
    facet_grid(~gene, scales='free') +
    labs(x='Pseudotime bins', y='Expression', color='Modality')

p2 <- ggplot(plot_df, aes(as.numeric(pt_bins), expr, color=modality)) +
    geom_smooth(method=mgcv::gam, formula = y ~ s(x, bs = 'cs')) +
    geom_point(size=0.2) +
    scale_color_manual(values=mod_colors) +
    facet_grid(modality~gene, scales='free') +
    labs(x='Pseudotime bins', y='Expression', color='Modality')

p1 / p2 + plot_layout(heights=c(1,3))
ggsave('plots/paper/fig3/fig3_ctx_pt_gene_expr_line.pdf', width=25, height=5)



#### GAMs for on pt_bin gene expression -> represent each timepoint by all modalities 
# -> multivariate DTW distance between GAM smooths
# -> K-means clustering -> identify lagging/parallel/antiparallel clusters of genes

#### Feature selection ####
H3K27ac_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K27ac_ctx_subs[['cRNA']]@data > 0)*1)
H3K27me3_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K27me3_ctx_subs[['cRNA']]@data > 0)*1)
H3K4me3_ctx_subs[['cRNA_bin']] <- CreateAssayObject((H3K4me3_ctx_subs[['cRNA']]@data > 0)*1)
rna_ctx_subs[['RNA_bin']] <- CreateAssayObject((rna_ctx_subs[['RNA']]@data > 0)*1)

H3K27ac_ctx_subs <- Pando::aggregate_assay(H3K27ac_ctx_subs, assay='cRNA_bin', group_name='pt_bins')
H3K27me3_ctx_subs <- Pando::aggregate_assay(H3K27me3_ctx_subs, assay='cRNA_bin', group_name='pt_bins')
H3K4me3_ctx_subs <- Pando::aggregate_assay(H3K4me3_ctx_subs, assay='cRNA_bin', group_name='pt_bins')
rna_ctx_subs <- Pando::aggregate_assay(rna_ctx_subs, assay='RNA_bin', group_name='pt_bins')

H3K27ac_cluster_detect <- H3K27ac_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), mark_feats]
H3K27me3_cluster_detect <- H3K27me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), mark_feats]
H3K4me3_cluster_detect <- H3K4me3_ctx_subs@assays$cRNA@misc$summary$pt_bins[as.character(1:50), mark_feats]
rna_cluster_detect <- rna_ctx_subs@assays$RNA@misc$summary$pt_bins[as.character(1:50),]

H3K27ac_maxdetect_feats <- colMaxs(H3K27ac_cluster_detect)
H3K27me3_maxdetect_feats <- colMaxs(H3K27me3_cluster_detect)
H3K4me3_maxdetect_feats <- colMaxs(H3K4me3_cluster_detect)
rna_maxdetect_feats <- colMaxs(rna_cluster_detect)

mark_detect_feats <- mark_feats[
    (H3K27ac_maxdetect_feats > 0.02) & (H3K27me3_maxdetect_feats > 0.02) & (H3K4me3_maxdetect_feats > 0.02)
]

rna_ctx_subs <- FindVariableFeatures(rna_ctx_subs, nfeatures=6000)
H3K27ac_ctx_subs <- FindVariableFeatures(H3K27ac_ctx_subs, assay='cRNA', nfeatures=6000)
H3K27me3_ctx_subs <- FindVariableFeatures(H3K27me3_ctx_subs, assay='cRNA', nfeatures=6000)
H3K4me3_ctx_subs <- FindVariableFeatures(H3K4me3_ctx_subs, assay='cRNA', nfeatures=6000)

all_var_feats <- bind_rows(
    'RNA' = as_tibble(rna_ctx_subs[['RNA']]@meta.features, rownames='feature'),
    'H3K27ac' = as_tibble(H3K27ac_ctx_subs[['cRNA']]@meta.features, rownames='feature'),
    'H3K27me3' = as_tibble(H3K27me3_ctx_subs[['cRNA']]@meta.features, rownames='feature'),
    'H3K4me3' = as_tibble(H3K4me3_ctx_subs[['cRNA']]@meta.features, rownames='feature'),
    .id = 'modality'
)

feats_use <- all_var_feats %>% 
    filter(modality=='RNA') %>% 
    top_n(6000, vst.variance) %>% 
    pull(feature) %>% intersect(mark_detect_feats)


#### Smooth with GAMs ####
x <- 1:50
H3K27ac_gams <- H3K27ac_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27ac_smooths <- H3K27ac_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27me3_gams <- H3K27me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K27me3_smooths <- H3K27me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K4me3_gams <- H3K4me3_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
H3K4me3_smooths <- H3K4me3_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

rna_gams <- rna_clusters[, feats_use] %>% 
    pbapply(2, function(y){
        mgcv::gam(formula = y ~ s(x, bs = 'cs'))
    })
rna_smooths <- rna_gams %>% map_dfr(~.x$fitted.values) %>% as.matrix()

H3K27ac_smooths_scale <- H3K27ac_smooths %>% scale()
H3K27me3_smooths_scale <- H3K27me3_smooths %>% scale()
H3K4me3_smooths_scale <- H3K4me3_smooths %>% scale()
rna_smooths_scale <- rna_smooths %>% scale()

all_smooths_list <- map(purrr::set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths_scale[,n],
            H3K27me3=H3K27me3_smooths_scale[,n],
            H3K4me3=H3K4me3_smooths_scale[,n],
            rna=rna_smooths_scale[,n]
        )
    )
})

all_smooths_raw_list <- map(purrr::set_names(colnames(rna_smooths)), function(n){
    do.call(
        cbind, 
        list(
            H3K27ac=H3K27ac_smooths[,n],
            H3K27me3=H3K27me3_smooths[,n],
            H3K4me3=H3K4me3_smooths[,n],
            rna=rna_smooths[,n]
        )
    )
})




#### Distances with DTW and L2 (euclidean) ####
library(dtw)
library(doParallel)
registerDoParallel(36)

dtw_smooth_dist <- map_par(names(all_smooths_list), function(i){
    map_dfr(names(all_smooths_list), function(j){
        ij_align <- dtw(all_smooths_list[[i]], all_smooths_list[[j]])
        ij_dist <- ij_align$distance
        return(tibble(
            i=i, j=j, dist=ij_dist
        ))
    })
}, parallel=T)

dtw_smooth_dist_df <- bind_rows(dtw_smooth_dist)
dtw_smooth_dist_mat <- dtw_smooth_dist_df %>%
    pivot_wider(names_from=j, values_from=dist) %>%
    column_to_rownames('i') %>% as.matrix()

pheatmap::pheatmap(dtw_smooth_dist_mat)


#### K-means clustering ####
library(FCPS)

set.seed(111)
dtw_kmeans <- FCPS::kmeansClustering(
    DataOrDistances = dtw_smooth_dist_mat,
    ClusterNo = 30
)

dtw_kmeans_df <- dtw_kmeans$Cls %>% enframe('feature', 'dtw_clust')


#### Hclust per cluster ####
cluster_hclusts <- map(set_names(sort(unique(dtw_kmeans_df$dtw_clust))), function(n){
    cl_genes <- filter(dtw_kmeans_df, dtw_clust==n)$feature
    cl_expr <- dtw_smooth_dist_mat[cl_genes, cl_genes]
    cl_expr %>% as.dist() %>% hclust()
})

cluster_hclust_order <- cluster_hclusts %>% map(~.x$labels[.x$order])
cluster_hclust_clusts <- cluster_hclusts %>% map_dfr(~enframe(cutree(.x, h=150), 'feature', 'hclust'), .id='dtw_clust') %>%
    mutate(dtw_clust=as.numeric(dtw_clust))


#### Plot clusters #####
all_gene_smooths_df <- map_dfr(all_smooths_raw_list, function(x){
    x %>% 
        {colnames(.) <- c('H3K27ac', 'H3K27me3', 'H3K4me3', 'RNA');.} %>% 
        as_tibble() %>% 
        mutate(pt_bin=1:50) %>% 
        return()
}, .id='feature') %>% 
    pivot_longer(H3K27ac:RNA, names_to='modality', values_to='expr')


plot_df <- all_gene_smooths_df %>% 
    inner_join(dtw_kmeans_df) %>%
    # inner_join(cluster_hclust_clusts) %>%
    group_by(feature, modality) %>% 
    mutate(
        expr01=scale01(expr), 
        pt_bin_mod=paste0(modality, pt_bin),
        feature=factor(feature, levels=purrr::reduce(cluster_hclust_order, c)),
        # dtw_clust=factor(dtw_clust, levels=c(7,1,6,10,4,9,8,2,3,5))
    )

ggplot(plot_df, aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(dtw_clust~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100)) +
    theme(
        panel.border = element_blank()
    )
ggsave('plots/paper/fig3/fig3_ctx_gam_smooth_30clust_labels_heatmap.pdf', width=15, height=60, limitsize=FALSE)
ggsave('plots/paper/fig3/fig3_ctx_gam_smooth_30clust_labels_heatmap.png', width=15, height=60, limitsize=FALSE)



ggplot(plot_df, aes(pt_bin, feature, fill=expr01)) +
    geom_tile() +
    facet_grid(dtw_clust~modality, scales='free', space='free') +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::magma(100)) +
    article_text() +
    no_y_text() +
    theme(
        panel.border = element_blank(),
        panel.spacing = unit(0.1, 'cm')
    ) + no_legend() +
    labs(x='Pseudotime bin', y='Feature')
ggsave('plots/paper/fig3/fig3_ctx_gam_smooth_30clust_heatmap.pdf', width=6, height=15, unit='cm')
ggsave('plots/paper/fig3/fig3_ctx_gam_smooth_30clust_heatmap.pdf', width=6, height=15, unit='cm')


plot_df %>% ungroup() %>% distinct(feature, dtw_clust) %>% write_tsv('data/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')


#### GO enrichment of clusters ####
library(clusterProfiler)
library(org.Hs.eg.db)
library(doParallel)
registerDoParallel(10)

gene_clusters <- read_tsv('data/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')

# Against all RNA features
feats_use <- all_var_feats %>% 
    filter(modality=='RNA') %>% 
    top_n(6000, vst.variance) %>% 
    pull(feature) 

all_features <- clusterProfiler::bitr(feats_use, fromType = 'SYMBOL', toType = c('ENSEMBL', 'ENTREZID'), OrgDb = org.Hs.eg.db) %>%
    as_tibble()

cluster_ego <- map_par(set_names(unique(gene_clusters$dtw_clust)), function(x){
    clust_genes <- filter(gene_clusters, dtw_clust==x)$feature %>% unique()
    gene_ids <- filter(all_features, SYMBOL%in%clust_genes)
    ego <- enrichGO(
        gene = gene_ids$ENTREZID,
        universe = all_features$ENTREZID,
        OrgDb = org.Hs.eg.db,
        pvalueCutoff = 0.05,
        # qvalueCutoff = 0.05,
        pAdjustMethod = 'fdr',
        ont = 'ALL',
        readable = T
    )
    print(ego)
    return(ego)
})

cluster_ego_rna_bg <- map_dfr(cluster_ego, ~if (!is.null(.x)){as_tibble(.x@result)} , .id='dtw_clust') %>% 
    # filter(ONTOLOGY=='BP', p.adjust<0.01) %>% 
    mutate(ngenes = as.numeric(str_replace(GeneRatio, '\\d+/(\\d+)', '\\1'))) %>% 
    mutate(ngenes_enrich = as.numeric(str_replace(GeneRatio, '(\\d+)/\\d+', '\\1'))) %>% 
    mutate(bggenes = as.numeric(str_replace(BgRatio, '\\d+/(\\d+)', '\\1'))) %>% 
    mutate(bggenes_enrich = as.numeric(str_replace(BgRatio, '(\\d+)/\\d+', '\\1'))) %>% 
    mutate(oddsratio=(ngenes_enrich/ngenes)/(bggenes_enrich/bggenes)) %>% 
    arrange(oddsratio) %>% 
    mutate(Description=factor(Description, levels=unique(.$Description))) %>% 
    arrange(p.adjust) 

cluster_ego_rna_bg %>% write_tsv('data/trajectories/ctx/all_mod_expr_dtw_30clust_go_enrich_rna_bg.tsv')









