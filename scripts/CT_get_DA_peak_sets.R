source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(doParallel)

setwd('~/projects/cutntag/')

select <- dplyr::select
filter <- dplyr::filter

#### Read data ####
marks <- read_rds('data/all_marks_list_v3.4lines.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.3lines.rds')

celltype_da_files <- list.files('data/results/diff_expression', pattern='*celltype.tsv', full.names=T)
names(celltype_da_files) <- str_replace(celltype_da_files, '.+expression/(.+)_DA_.+', '\\1')
celltype_da_list <- celltype_da_files %>% map(read_tsv)

age_da_files <- list.files('data/results/diff_expression', pattern='*age.tsv', full.names=T)
names(age_da_files) <- str_replace(age_da_files, '.+expression/(.+)_DA_.+', '\\1')
age_da_list <- age_da_files %>% map(read_tsv)

lineage_da_files <- list.files('data/results/diff_expression', pattern='*lineage_all.tsv', full.names=T)
names(lineage_da_files) <- str_replace(lineage_da_files, '.+expression/(.+)_DA_.+', '\\1')
lineage_da_list <- lineage_da_files %>% map(read_tsv)

gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')


#### Select DA peaks and write as bed ####
map(names(celltype_da_list), function(x){
    da_df <- celltype_da_list[[x]]
    ct_peaks <- da_df %>%
        mutate(padj=p.adjust(pval, method='fdr')) %>%
        filter(padj<1e-4, detect_self>0.01, log_dr>1) %>%
        pull(feature) %>% unique()

    export.bed(StringToGRanges(ct_peaks), paste0('data/results/diff_expression/', x, '_celltype_DA_peaks.bed'))
})

map(names(age_da_list), function(x){
    da_df <- age_da_list[[x]]
    ct_peaks <- da_df %>%
        mutate(padj=p.adjust(pval, method='fdr')) %>%
        filter(padj<1e-4, detect_self>0.01, log_dr>1) %>%
        pull(feature) %>% unique()

    export.bed(StringToGRanges(ct_peaks), paste0('data/results/diff_expression/', x, '_age_DA_peaks.bed'))
})

ct_peaks <- age_da_list$H3K27ac %>%
    mutate(padj=p.adjust(pval, method='fdr')) %>%
    filter(padj<1e-4, detect_self>0.1, log_dr>1) %>%
    top_n(100, log_dr) %>%
    pull(feature) %>% unique()

export.bed(StringToGRanges(ct_peaks), paste0('data/results/diff_expression/test_DA_peaks.bed'))


#### Export all peaks ####
map(names(marks), function(x){
    srt <- marks[[x]]
    export.bed(StringToGRanges(rownames(srt[['peaks']])), paste0('data/results/', x, '_all_cellranger_peaks.bed'))
    export.bed(StringToGRanges(rownames(srt[['narrow_peaks']])), paste0('data/results/', x, '_all_MACS_peaks.bed'))
})


#### Make heatmaps with DA results ####

plot_peak_heatmap <- function(de_df){
    de_mat <- de_df %>%
        select(group, feature, p_clip) %>%
        pivot_wider(names_from = group, values_from = p_clip, values_fill = 0) %>%
        column_to_rownames('feature') %>% as.matrix()

    peak_order <- de_mat %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}
    group_order <- de_mat %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}


    plot_df <- de_df %>%
        mutate(feature=factor(feature, levels=peak_order), group=factor(group, levels=group_order))


    p <- ggplot(plot_df, aes(group, feature, fill=p_clip)) +
        geom_tile() +
        scale_fill_gradientn(colors=pals::brewer.puor(100), limits=c(-50, 50)) +
        no_y_text() +
        scale_x_discrete(expand=c(0,0)) +
        labs(y='Peak', x='Group', fill='Signed log10(FDR)')
    return(p)
}



#### H3K27ac ####
H3K27ac_use <- marks$H3K27ac %>% subset(lineage%in%c('ctx', 'dien', 'nt', 'retina'))

H3K27ac_df <- lineage_da_list$H3K27ac %>%
    mutate(
        padj=p.adjust(pval, method='fdr'),
        log_dr=ifelse(abs(log_dr)>5, sign(log_dr)*5, log_dr),
        signed_p=sign(coef)*-log10(padj),
        p_clip=ifelse(abs(signed_p)>50, sign(signed_p)*50, signed_p)
    ) %>%
    group_by(feature) %>%
    dplyr::filter(any(padj<1e-4 & detect_self>0.05 & coef>0))

H3K27ac_de_mat <- H3K27ac_df %>%
    dplyr::select(group, feature, p_clip) %>%
    pivot_wider(names_from = group, values_from = p_clip, values_fill = 0) %>%
    column_to_rownames('feature') %>% as.matrix()

H3K27ac_peak_order <- H3K27ac_de_mat %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}
H3K27ac_group_order <- H3K27ac_de_mat %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}
    
H3K27ac_top <- H3K27ac_df %>% 
    group_by(group) %>% 
    top_n(100, coef)

regions_use <- unique(H3K27ac_top$feature) %>% 
    StringToGRanges() %>% 
    # resize(fix='center', width = 1) %>% 
    Extend(upstream=2000, downstream=2000)

peak_genes <- Signac::ClosestFeature(marks$H3K27ac, regions_use) %>% 
    as_tibble() %>% mutate(feature=unique(H3K27ac_top$feature)) %>% 
    distinct(feature, gene_name)

H3K27ac_cov_plots <- map_par(1:length(regions_use), function(i){
    reg <- regions_use[i]
    p <- CoveragePlot(H3K27ac_use, region = reg, annotation = F, peaks = F, group.by = 'lineage', window = 500)
    return(p)
})
names(H3K27ac_cov_plots) <- unique(H3K27ac_top$feature)

H3K27ac_cov_df <- map_dfr(H3K27ac_cov_plots, ~.x$data, .id='feature') %>% 
    inner_join(H3K27ac_df) %>% 
    group_by(group, feature) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, feature, poscut) %>% 
    summarize(coverage=mean(coverage), p_clip=p_clip[1])

H3K27ac_cov_mat <- H3K27ac_cov_df %>% 
    mutate(group_pos=paste0(group, poscut)) %>% 
    ungroup() %>% 
    select(feature, group_pos, coverage) %>% 
    pivot_wider(names_from = group_pos, values_from = coverage, values_fill=0) %>% 
    column_to_rownames('feature') %>% as.matrix()

peak_order <- H3K27ac_cov_mat %>% dist() %>% hclust() %>% {.$label[.$order]}

max_group_df <- H3K27ac_top %>% 
    select(feature, 'max_group'=group)
    
plot_df <- H3K27ac_cov_df %>% 
    inner_join(peak_genes) %>% 
    inner_join(max_group_df) %>% 
    group_by(gene_name) %>% 
    mutate(feature_num=as.numeric(factor(feature))) %>% 
    group_by(feature) %>%
    mutate(
        gene_peak=paste0(gene_name, '_', feature_num),
        coverage=coverage/max(coverage), 
        feature=factor(feature, levels=peak_order),
        group=factor(group, levels=c('ctx', 'dien', 'nt', 'retina')),
        max_group=factor(max_group, levels=c('ctx', 'dien', 'nt', 'retina'))
    )  %>% 
    arrange(feature) %>% mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak)))
    
ggplot(plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.2) +
    facet_grid(max_group~group, scales='free') +
    scale_fill_gradientn(colors=pals::brewer.ylgn(100), limits=c(-50, 50)) +
    no_x_text() +
    theme_void() +
    article_text() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) 

ggsave('plots/cutntag/H3K27ac_lineage_DA_top100_ridges.png', width=5, height=30, bg='white')






#### H3K27me3 ####
H3K27me3_use <- marks$H3K27me3 %>% subset(lineage%in%c('ctx', 'dien', 'nt', 'retina'))

H3K27me3_df <- lineage_da_list$H3K27me3 %>%
    mutate(
        padj=p.adjust(pval, method='fdr'),
        log_dr=ifelse(abs(log_dr)>5, sign(log_dr)*5, log_dr),
        signed_p=sign(coef)*-log10(padj),
        p_clip=ifelse(abs(signed_p)>50, sign(signed_p)*50, signed_p)
    ) %>%
    group_by(feature) %>%
    dplyr::filter(any(padj<1e-4 & detect_self>0.05 & coef>0))

H3K27me3_de_mat <- H3K27me3_df %>%
    dplyr::select(group, feature, p_clip) %>%
    pivot_wider(names_from = group, values_from = p_clip, values_fill = 0) %>%
    column_to_rownames('feature') %>% as.matrix()

H3K27me3_peak_order <- H3K27me3_de_mat %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}
H3K27me3_group_order <- H3K27me3_de_mat %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

H3K27me3_top <- H3K27me3_df %>% 
    group_by(group) %>% 
    top_n(100, coef)

regions_use <- unique(H3K27me3_top$feature) %>% 
    StringToGRanges() %>% 
    # resize(fix='center', width = 1) %>% 
    Extend(upstream=2000, downstream=2000)

peak_genes <- Signac::ClosestFeature(marks$H3K27me3, regions_use) %>% 
    as_tibble() %>% mutate(feature=unique(H3K27me3_top$feature)) %>% 
    distinct(feature, gene_name)

H3K27me3_cov_plots <- map_par(1:length(regions_use), function(i){
    reg <- regions_use[i]
    p <- CoveragePlot(H3K27me3_use, region = reg, annotation = F, peaks = F, group.by = 'lineage', window = 500)
    return(p)
})
names(H3K27me3_cov_plots) <- unique(H3K27me3_top$feature)

H3K27me3_cov_df <- map_dfr(H3K27me3_cov_plots, ~.x$data, .id='feature') %>% 
    inner_join(H3K27me3_df) %>% 
    group_by(group, feature) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, feature, poscut) %>% 
    summarize(coverage=mean(coverage), p_clip=p_clip[1])

H3K27me3_cov_mat <- H3K27me3_cov_df %>% 
    mutate(group_pos=paste0(group, poscut)) %>% 
    ungroup() %>% 
    select(feature, group_pos, coverage) %>% 
    pivot_wider(names_from = group_pos, values_from = coverage, values_fill=0) %>% 
    column_to_rownames('feature') %>% as.matrix()

peak_order <- H3K27me3_cov_mat %>% dist() %>% hclust() %>% {.$label[.$order]}

max_group_df <- H3K27me3_top %>% 
    select(feature, 'max_group'=group)

plot_df <- H3K27me3_cov_df %>% 
    inner_join(peak_genes) %>% 
    inner_join(max_group_df) %>% 
    group_by(gene_name) %>% 
    mutate(feature_num=as.numeric(factor(feature))) %>% 
    group_by(feature) %>%
    mutate(
        gene_peak=paste0(gene_name, '_', feature_num),
        coverage=coverage/max(coverage), 
        feature=factor(feature, levels=peak_order),
        group=factor(group, levels=c('ctx', 'dien', 'nt', 'retina')),
        max_group=factor(max_group, levels=c('ctx', 'dien', 'nt', 'retina'))
    )  %>% 
    arrange(feature) %>% mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak)))

ggplot(plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.2) +
    facet_grid(max_group~group, scales='free') +
    scale_fill_gradientn(colors=pals::brewer.ylgnbu(100), limits=c(-50, 50)) +
    no_x_text() +
    theme_void() +
    article_text() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) 

ggsave('plots/cutntag/H3K27me3_lineage_DA_top100_ridges.png', width=5, height=30, bg='white')






#### H3K4me3 ####
H3K4me3_use <- marks$H3K4me3 %>% subset(lineage%in%c('ctx', 'dien', 'nt', 'retina'))

H3K4me3_df <- lineage_da_list$H3K4me3 %>%
    mutate(
        padj=p.adjust(pval, method='fdr'),
        log_dr=ifelse(abs(log_dr)>5, sign(log_dr)*5, log_dr),
        signed_p=sign(coef)*-log10(padj),
        p_clip=ifelse(abs(signed_p)>50, sign(signed_p)*50, signed_p)
    ) %>%
    group_by(feature) %>%
    dplyr::filter(any(padj<1e-4 & detect_self>0.05 & coef>0))

H3K4me3_de_mat <- H3K4me3_df %>%
    dplyr::select(group, feature, p_clip) %>%
    pivot_wider(names_from = group, values_from = p_clip, values_fill = 0) %>%
    column_to_rownames('feature') %>% as.matrix()

H3K4me3_peak_order <- H3K4me3_de_mat %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}
H3K4me3_group_order <- H3K4me3_de_mat %>% t() %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

H3K4me3_top <- H3K4me3_df %>% 
    group_by(group) %>% 
    top_n(100, coef)

regions_use <- unique(H3K4me3_top$feature) %>% 
    StringToGRanges() %>% 
    # resize(fix='center', width = 1) %>% 
    Extend(upstream=2000, downstream=2000)

peak_genes <- Signac::ClosestFeature(marks$H3K4me3, regions_use) %>% 
    as_tibble() %>% mutate(feature=unique(H3K4me3_top$feature)) %>% 
    distinct(feature, gene_name)

H3K4me3_cov_plots <- map_par(1:length(regions_use), function(i){
    reg <- regions_use[i]
    p <- CoveragePlot(H3K4me3_use, region = reg, annotation = F, peaks = F, group.by = 'lineage', window = 500)
    return(p)
})
names(H3K4me3_cov_plots) <- unique(H3K4me3_top$feature)

H3K4me3_cov_df <- map_dfr(H3K4me3_cov_plots, ~.x$data, .id='feature') %>% 
    inner_join(H3K4me3_df) %>% 
    group_by(group, feature) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, feature, poscut) %>% 
    summarize(coverage=mean(coverage), p_clip=p_clip[1])

H3K4me3_cov_mat <- H3K4me3_cov_df %>% 
    mutate(group_pos=paste0(group, poscut)) %>% 
    ungroup() %>% 
    select(feature, group_pos, coverage) %>% 
    pivot_wider(names_from = group_pos, values_from = coverage, values_fill=0) %>% 
    column_to_rownames('feature') %>% as.matrix()

peak_order <- H3K4me3_cov_mat %>% dist() %>% hclust() %>% {.$label[.$order]}

max_group_df <- H3K4me3_top %>% 
    select(feature, 'max_group'=group)

plot_df <- H3K4me3_cov_df %>% 
    inner_join(peak_genes) %>% 
    inner_join(max_group_df) %>% 
    group_by(gene_name) %>% 
    mutate(feature_num=as.numeric(factor(feature))) %>% 
    group_by(feature) %>%
    mutate(
        gene_peak=paste0(gene_name, '_', feature_num),
        coverage=coverage/max(coverage), 
        feature=factor(feature, levels=peak_order),
        group=factor(group, levels=c('ctx', 'dien', 'nt', 'retina')),
        max_group=factor(max_group, levels=c('ctx', 'dien', 'nt', 'retina'))
    ) %>% 
    arrange(feature) %>% mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak)))

ggplot(plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.2) +
    facet_grid(max_group~group, scales='free') +
    scale_fill_gradientn(colors=pals::brewer.orrd(100), limits=c(-50, 50)) +
    no_x_text() +
    theme_void() +
    article_text() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) 

ggsave('plots/cutntag/H3K4me3_lineage_DA_top100_ridges.png', width=5, height=30, bg='white')


all_DA_df <- bind_rows('H3K4me3'=H3K4me3_df, 'H3K27ac'=H3K27ac_df, 'H3K27me3'=H3K27me3_df, .id='mark')
all_DA_df %>% write_tsv('data/results/diff_expression/all_marks_lineage_DA.tsv')


######## Switches for all three marks ####
# Get intersects 
H3K27me3_ranges <- StringToGRanges(rownames(marks$H3K27me3))
H3K27me3_ranges_extended <- H3K27me3_ranges %>% 
    Extend(upstream=1000, downstream=1000)

H3K27ac_ranges <- StringToGRanges(rownames(marks$H3K27ac))
H3K27ac_ranges_extended <- H3K27ac_ranges %>% 
    Extend(upstream=1000, downstream=1000)

H3K4me3_ranges <- StringToGRanges(rownames(marks$H3K4me3))
H3K4me3_ranges_extended <- H3K4me3_ranges %>% 
    Extend(upstream=1000, downstream=1000)

all_intersect <- IRanges::intersect(H3K27me3_ranges_extended, H3K27ac_ranges_extended) %>% 
    IRanges::intersect(H3K4me3_ranges_extended)

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

H3K4me3_olaps <- findOverlaps(H3K4me3_ranges_extended, all_intersect)
H3K4me3_matches <- tibble(
    isect = GRangesToString(all_intersect)[subjectHits(H3K4me3_olaps)],
    H3K4me3 = rownames(marks$H3K4me3)[queryHits(H3K4me3_olaps)]
)

all_isect_matches <- inner_join(H3K27me3_matches, H3K27ac_matches) %>% 
    inner_join(H3K4me3_matches)

all_isect_matches %>% write_tsv('data/intersect/all_marks_intersect_matches.tsv')



#### Filter DA peaks by intersections ####
all_isect_matches <- read_tsv('data/intersect/all_marks_intersect_matches.tsv')
all_DA_df <- read_tsv('data/results/diff_expression/all_marks_lineage_DA.tsv')

# H3K27ac
H3K27ac_da <- all_DA_df %>% filter(mark=='H3K27ac')
H3K27ac_isect_peaks <- intersect(H3K27ac_da$feature, all_isect_matches$H3K27ac)

H3K27ac_da_isect <- H3K27ac_da %>%
    dplyr::filter(feature%in%H3K27ac_isect_peaks) %>% 
    mutate(mark='H3K27ac')


# H3K27me3
H3K27me3_da <- all_DA_df %>% filter(mark=='H3K27me3')
H3K27me3_isect_peaks <- intersect(H3K27me3_da$feature, all_isect_matches$H3K27me3)

H3K27me3_da_isect <- H3K27me3_da %>%
    dplyr::filter(feature%in%H3K27me3_isect_peaks) %>% 
    mutate(mark='H3K27me3')


# H3K4me3
H3K4me3_da <- all_DA_df %>% filter(mark=='H3K4me3')
H3K4me3_isect_peaks <- intersect(H3K4me3_da$feature, all_isect_matches$H3K4me3)

H3K4me3_da_isect <- H3K4me3_da %>%
    dplyr::filter(feature%in%H3K4me3_isect_peaks) %>% 
    mutate(mark='H3K4me3')

# Combine
all_da_isect <- all_isect_matches %>% 
    inner_join(H3K27ac_da_isect, by=c('H3K27ac'='feature')) %>% 
    inner_join(H3K27me3_da_isect, by=c('H3K27me3'='feature', 'group'), suffix=c('_H3K27ac', '_H3K27me3')) %>% 
    inner_join(H3K4me3_da_isect, by=c('H3K4me3'='feature', 'group'), suffix=c('', '_H3K4me3'))

all_da_rbind <-  inner_join(all_isect_matches, H3K27ac_da_isect, by=c('H3K27ac'='feature')) %>% 
    bind_rows(inner_join(all_isect_matches, H3K27me3_da_isect, by=c('H3K27me3'='feature'))) %>% 
    bind_rows(inner_join(all_isect_matches, H3K4me3_da_isect, by=c('H3K4me3'='feature'))) %>% 
    filter(isect%in%all_da_isect$isect)
    
all_da_mat <- all_da_rbind %>%
    mutate(mark_group=paste0(mark, '_', group)) %>% 
    dplyr::select(mark_group, isect, p_clip) %>%
    distinct() %>% 
    pivot_wider(names_from = mark_group, values_from = p_clip, values_fill = 0, values_fn=mean) %>%
    column_to_rownames('isect') %>% as.matrix()

isect_peak_order <- all_da_mat %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}

#### Plot switches with heatmap ####
plot_df <- all_da_rbind %>% 
    mutate(isect=factor(isect, levels=isect_peak_order))

ggplot(plot_df, aes(mark, isect, fill=p_clip)) +
    geom_tile() +
    facet_grid(~group) +
    scale_fill_gradientn(colors=pals::brewer.puor(100), limits=c(-50, 50)) +
    no_y_text() +
    rotate_x_text(40) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Peak', x='Group', fill='Signed log10(FDR)')



#### Plot switches with ridge plot ####
regions_use <- unique(all_da_rbind$isect) %>% 
    StringToGRanges() %>% 
    # resize(fix='center', width = 1) %>% 
    Extend(upstream=1000, downstream=1000)

peak_genes <- Signac::ClosestFeature(marks$H3K27ac, regions_use) %>% 
    as_tibble() %>% mutate(feature=unique(all_da_rbind$isect)) %>% 
    distinct(feature, gene_name)

H3K27ac_cov_plots <- map_par(1:length(regions_use), function(i){
    reg <- regions_use[i]
    p <- CoveragePlot(H3K27ac_use, region = reg, annotation = F, peaks = F, group.by = 'lineage', window = 500)
    return(p)
})
names(H3K27ac_cov_plots) <- unique(all_da_rbind$isect)

H3K27me3_cov_plots <- map_par(1:length(regions_use), function(i){
    reg <- regions_use[i]
    p <- CoveragePlot(H3K27me3_use, region = reg, annotation = F, peaks = F, group.by = 'lineage', window = 500)
    return(p)
})
names(H3K27me3_cov_plots) <- unique(all_da_rbind$isect)

H3K4me3_cov_plots <- map_par(1:length(regions_use), function(i){
    reg <- regions_use[i]
    p <- CoveragePlot(H3K4me3_use, region = reg, annotation = F, peaks = F, group.by = 'lineage', window = 500)
    return(p)
})
names(H3K4me3_cov_plots) <- unique(all_da_rbind$isect)

H3K27ac_cov_df <- map_dfr(H3K27ac_cov_plots, ~.x$data, .id='isect') 
H3K27me3_cov_df <- map_dfr(H3K27me3_cov_plots, ~.x$data, .id='isect') 
H3K4me3_cov_df <- map_dfr(H3K4me3_cov_plots, ~.x$data, .id='isect') 

all_cov_df <- bind_rows('H3K27ac'=H3K27ac_cov_df,'H3K27me3'=H3K27me3_cov_df,'H3K4me3'=H3K4me3_cov_df, .id='mark') %>% 
    group_by(group, isect, mark) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, isect, poscut, mark) %>% 
    summarize(coverage=mean(coverage))

all_cov_mat <- all_cov_df %>% 
    mutate(group_pos=paste0(mark, group, poscut)) %>% 
    ungroup() %>% 
    select(isect, group_pos, coverage) %>% 
    pivot_wider(names_from = group_pos, values_from = coverage, values_fill=0) %>% 
    column_to_rownames('isect') %>% as.matrix()

peak_order <- all_cov_mat %>% dist() %>% hclust() %>% {.$label[.$order]}


H3K27ac_cov_fmt <- H3K27ac_cov_df %>% 
    group_by(group, isect) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, isect, poscut) %>% 
    summarize(coverage=mean(coverage))

H3K27me3_cov_fmt <- H3K27me3_cov_df %>% 
    group_by(group, isect) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, isect, poscut) %>% 
    summarize(coverage=mean(coverage))

H3K4me3_cov_fmt <- H3K4me3_cov_df %>% 
    group_by(group, isect) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, isect, poscut) %>% 
    summarize(coverage=mean(coverage))


all_plot_df <- bind_rows('H3K27ac'=H3K27ac_cov_fmt,'H3K27me3'=H3K27me3_cov_fmt,'H3K4me3'=H3K4me3_cov_fmt, .id='mark') %>% 
    inner_join(peak_genes, by=c('isect'='feature')) %>% 
    group_by(gene_name) %>% 
    mutate(feature_num=as.numeric(factor(isect))) %>% 
    group_by(isect, mark) %>%
    mutate(
        gene_peak=paste0(gene_name, '_', feature_num),
        coverage=coverage/max(coverage), 
        isect=factor(isect, levels=isect_peak_order),
        group=factor(group, levels=c('ctx', 'dien', 'nt', 'retina'))
    )  %>% 
    arrange(isect) %>% 
    mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak)))



ggplot(all_plot_df, aes(poscut, gene_peak, height=coverage, fill=mark)) +
    geom_ridgeline(color='white', size=0.1, scale=1.2) +
    facet_grid(~group+mark, scales='free') +
    scale_fill_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B')) +
    no_x_text() +
    theme_void() +
    article_text() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) 

ggsave('plots/cutntag/all_marks_DA_isect_ridges.png', width=10, height=20, bg='white')



H3K27ac_da_isect <- all_da_rbind %>% filter(mark=='H3K27ac') %>% 
    distinct(isect, group, p_clip) %>% 
    group_by(isect, group) %>% 
    summarise(p_clip=mean(p_clip))
    
H3K27me3_da_isect <- all_da_rbind %>% filter(mark=='H3K27me3') %>% 
    distinct(isect, group, p_clip) %>% 
    group_by(isect, group) %>% 
    summarise(p_clip=mean(p_clip))
    
H3K4me3_da_isect <- all_da_rbind %>% filter(mark=='H3K4me3') %>% 
    distinct(isect, group, p_clip) %>% 
    group_by(isect, group) %>% 
    summarise(p_clip=mean(p_clip))
    

H3K27ac_plot_df <- all_plot_df %>% 
    filter(mark=='H3K27ac') %>% 
    inner_join(H3K27ac_da_isect)
    
H3K27me3_plot_df <- all_plot_df %>% 
    filter(mark=='H3K27me3') %>% 
    inner_join(H3K27me3_da_isect)
    
H3K4me3_plot_df <- all_plot_df %>% 
    filter(mark=='H3K4me3') %>% 
    inner_join(H3K4me3_da_isect)
    

p1 <- ggplot(H3K27ac_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.2) +
    facet_grid(~group, scales='free') +
    scale_fill_gradientn(colors=pals::brewer.ylgn(100), limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) +
    ggtitle('H3K27ac')

p2 <- ggplot(H3K27me3_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.2) +
    facet_grid(~group, scales='free') +
    scale_fill_gradientn(colors=pals::brewer.ylgnbu(100), limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    no_y_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) +
    ggtitle('H3K27me3')

p3 <- ggplot(H3K4me3_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.2) +
    facet_grid(~group, scales='free') +
    scale_fill_gradientn(colors=pals::brewer.purd(100), limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    no_y_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) +
    ggtitle('H3K4me3')

(p1 | p2 | p3) + plot_layout(guides = 'collect')
ggsave('plots/cutntag/all_marks_DA_isect_split_ridges.png', width=10, height=15, bg='white')




######## H3K27ac - H3K27me3 switches + nepi ####
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
K27_isect_matches %>% write_tsv('data/intersect/H3K27_poised_marks_intersect_matches.tsv')


#### Filter DA peaks by intersections ####
K27_isect_matches <- read_tsv('data/intersect/H3K27_poised_marks_intersect_matches.tsv')
all_DA_df <- read_tsv('data/results/diff_expression/all_marks_lineage_DA.tsv')

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

K27_da_rbind <- inner_join(K27_isect_matches, H3K27ac_da_isect, by=c('H3K27ac'='feature')) %>% 
    bind_rows(inner_join(K27_isect_matches, H3K27me3_da_isect, by=c('H3K27me3'='feature'))) %>% 
    filter(isect%in%K27_da_isect$isect)

K27_da_mat <- K27_da_rbind %>%
    mutate(mark_group=paste0(mark, '_', group)) %>% 
    dplyr::select(mark_group, isect, p_clip) %>%
    distinct() %>% 
    pivot_wider(names_from = mark_group, values_from = p_clip, values_fill = 0, values_fn=mean) %>%
    column_to_rownames('isect') %>% as.matrix()

isect_peak_order <- K27_da_mat %>% scale() %>% dist() %>% hclust() %>% {.$label[.$order]}


#### Plot switches with heatmap ####
K27ac_max_df <- K27_da_isect %>% 
    group_by(isect) %>% 
    filter(abs(signed_p_H3K27ac)==max(abs(signed_p_H3K27ac))) %>% 
    filter(row_number()==1) %>% 
    arrange(group, p_clip_H3K27ac) %>% 
    select(isect, 'max_group'=group) 

plot_df <- K27_da_rbind %>% 
    inner_join(K27ac_max_df) %>% 
    # mutate(isect=factor(isect, levels=unique(K27ac_max_df$isect)))
    mutate(isect=factor(isect, levels=isect_peak_order))

ggplot(plot_df, aes(mark, isect, fill=p_clip)) +
    geom_tile() +
    facet_grid(~group, scales='free_y') +
    scale_fill_gradientn(colors=pals::brewer.puor(100), limits=c(-50, 50)) +
    no_y_text() +
    rotate_x_text(40) +
    scale_x_discrete(expand=c(0,0)) +
    labs(y='Peak', x='Group', fill='Signed log10(FDR)')


#### Plot switches with ridge plot ####
regions_use <- unique(K27_da_rbind$isect) %>% 
    StringToGRanges() %>% 
    # resize(fix='center', width = 1) %>% 
    Extend(upstream=1000, downstream=1000)

peak_genes <- Signac::ClosestFeature(marks$H3K27ac, regions_use) %>% 
    as_tibble() %>% mutate(feature=unique(K27_da_rbind$isect)) %>% 
    distinct(feature, gene_name)

H3K27ac_cov_plots <- map_par(1:length(regions_use), function(i){
    reg <- regions_use[i]
    p <- CoveragePlot(H3K27ac_use, region = reg, annotation = F, peaks = F, group.by = 'lineage', window = 500)
    return(p)
})
names(H3K27ac_cov_plots) <- unique(K27_da_rbind$isect)

H3K27me3_cov_plots <- map_par(1:length(regions_use), function(i){
    reg <- regions_use[i]
    p <- CoveragePlot(H3K27me3_use, region = reg, annotation = F, peaks = F, group.by = 'lineage', window = 500)
    return(p)
})
names(H3K27me3_cov_plots) <- unique(K27_da_rbind$isect)


H3K27ac_cov_df <- map_dfr(H3K27ac_cov_plots, ~.x$data, .id='isect') 
H3K27me3_cov_df <- map_dfr(H3K27me3_cov_plots, ~.x$data, .id='isect') 


K27_cov_df <- bind_rows('H3K27ac'=H3K27ac_cov_df,'H3K27me3'=H3K27me3_cov_df, .id='mark') %>% 
    group_by(group, isect, mark) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, isect, poscut, mark) %>% 
    summarize(coverage=mean(coverage))

K27_cov_mat <- K27_cov_df %>% 
    mutate(group_pos=paste0(mark, group, poscut)) %>% 
    ungroup() %>% 
    select(isect, group_pos, coverage) %>% 
    pivot_wider(names_from = group_pos, values_from = coverage, values_fill=0) %>% 
    column_to_rownames('isect') %>% as.matrix()

peak_order <- K27_cov_mat %>% dist() %>% hclust() %>% {.$label[.$order]}


H3K27ac_cov_fmt <- H3K27ac_cov_df %>% 
    group_by(group, isect) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, isect, poscut) %>% 
    summarize(coverage=mean(coverage))

H3K27me3_cov_fmt <- H3K27me3_cov_df %>% 
    group_by(group, isect) %>% 
    mutate(posnorm=scale01(position)) %>% 
    mutate(poscut=as.numeric(cut(posnorm, 500, labels=1:500))) %>% 
    group_by(group, isect, poscut) %>% 
    summarize(coverage=mean(coverage))


K27_plot_df <- bind_rows('H3K27ac'=H3K27ac_cov_fmt,'H3K27me3'=H3K27me3_cov_fmt, .id='mark') %>% 
    inner_join(peak_genes, by=c('isect'='feature')) %>% 
    inner_join(K27ac_max_df) %>% 
    group_by(gene_name) %>% 
    mutate(feature_num=as.numeric(factor(isect))) %>% 
    group_by(isect, mark) %>%
    mutate(
        gene_peak=paste0(gene_name, '_', feature_num),
        coverage=coverage/max(coverage), 
        isect=factor(isect, levels=isect_peak_order),
        group=factor(group, levels=c('ctx', 'dien', 'nt', 'retina'))
    )  %>% 
    arrange(isect) %>% 
    mutate(gene_peak=factor(gene_peak, levels=unique(.$gene_peak)))


ggplot(K27_plot_df, aes(poscut, gene_peak, height=coverage, fill=mark)) +
    geom_ridgeline(color='white', size=0.1, scale=1.2) +
    facet_grid(~group+mark, scales='free') +
    scale_fill_manual(values=c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B')) +
    no_x_text() +
    theme_void() +
    article_text() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) 

ggsave('plots/cutntag/K27_marks_DA_isect_ridges.png', width=10, height=20, bg='white')



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
    inner_join(H3K27ac_da_isect)

H3K27me3_plot_df <- K27_plot_df %>% 
    filter(mark=='H3K27me3') %>% 
    inner_join(H3K27me3_da_isect)


p1 <- ggplot(H3K27ac_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.3) +
    facet_grid(~group, scales='free') +
    scale_fill_gradientn(colors=pals::brewer.ylgn(100), limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines'),
        axis.text.y = element_text(size=5, hjust=1)
    ) +
    ggtitle('H3K27ac')

p2 <- ggplot(H3K27me3_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.3) +
    facet_grid(~group, scales='free') +
    scale_fill_gradientn(colors=pals::brewer.ylgnbu(100), limits=c(-50, 50)) +
    article_text() +
    no_x_text() +
    no_y_text() +
    theme_void() +
    theme(
        panel.spacing = unit(0.1,'lines')
    ) +
    ggtitle('H3K27me3')

(p1 | p2) + plot_layout(guides = 'collect')
ggsave('plots/cutntag/K27_marks_DA_isect_split_ridges.png', width=8, height=12, bg='white')



#### Add nepi tracks ####
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
    inner_join(peak_clusts) %>% 
    mutate(
        p_clip=ifelse(is.na(p_clip), 0, p_clip),
        group=factor(group, levels=c('nepi', 'ctx', 'dien', 'nt', 'retina')),
        region=factor(region, levels=peak_order),
        kmeans_clust=factor(kmeans_clust, levels=c(6,5,7,4,2,1,8,10,9,3))
    ) %>% 
    filter(!is.na(group))

H3K27me3_plot_df <- K27_plot_df %>% 
    filter(mark=='H3K27me3') %>% 
    left_join(H3K27me3_da_isect, by=c('region'='isect', 'group')) %>% 
    inner_join(peak_clusts) %>% 
    mutate(
        p_clip=ifelse(is.na(p_clip), 0, p_clip),
        group=factor(group, levels=c('nepi', 'ctx', 'dien', 'nt', 'retina')),
        region=factor(region, levels=peak_order),
        kmeans_clust=factor(kmeans_clust, levels=c(6,5,7,4,2,1,8,10,9,3))
    ) %>% 
    filter(!is.na(group))


p1 <- ggplot(H3K27ac_plot_df, aes(poscut, gene_peak, height=coverage, fill=p_clip)) +
    geom_ridgeline(color='white', size=0.1, scale=1.4) +
    facet_grid(kmeans_clust~group, scales='free_y', space='free_y') +
    scale_fill_gradientn(colors=pals::brewer.ylgn(100), limits=c(-50, 50)) +
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
    geom_ridgeline(color='white', size=0.1, scale=1.4) +
    facet_grid(kmeans_clust~group, scales='free_y', space='free_y') +
    scale_fill_gradientn(colors=pals::brewer.ylgnbu(100), limits=c(-50, 50)) +
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

ggsave('plots/cutntag/K27_marks_DA_isect_nepi_split_ridges.png', width=8, height=12, bg='white')



