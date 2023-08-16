library(tidyverse)
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rGREAT)


rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
dist <- stats::dist

data(SCREEN.ccRE.UCSC.hg38)
setwd('~/projects/cutntag/')

#### FUNC ####
mark_colors <- c('H3K4me3'='#CB9ACA', 'H3K27me3'='#3AAFC3', 'H3K27ac'='#5FBE9B', 'RNA'='#FDA044')
ct_order <- c('psc', 'non_nect', 'nect', 'ctx_npc', 'ctx_ip', 'ctx_ex', 'dien_npc', 'dien_ex', 
              'nt_npc', 'mesen_ex', 'rhom_ex', 'RPC', 'RGC', 'astrocytes', 'OPC', 'choroid_plexus')

tf_motif_enrichment <- function(regions, motifs, background_regions, parallel=T){
    data(motif2tf)
    group_vec <- background_regions %in% regions
    map_par(colnames(motifs), function(mot){
        ctab <- table(motifs[background_regions, mot], group_vec)
        ftest <- try(fisher.test(ctab))
        if (any(class(ftest)=='try-error')){
            print('asda')
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


#### Read data ####
marks <- read_rds('data05/all_marks_list_v3.3motifs.rds')
all_DA_df <- read_tsv('data_/results/diff_expression/all_marks_lineage_DA.tsv')
gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')
gene_collapse <- CollapseToLongestTranscript(gene_annot)
prom_ranges <- IRanges::promoters(gene_collapse)

marks$H3K27ac[['peaks_bin']] <- CreateAssayObject((marks$H3K27ac[['peaks']]@data > 0)*1)
marks$H3K27me3[['peaks_bin']] <- CreateAssayObject((marks$H3K27me3[['peaks']]@data > 0)*1)
marks$H3K4me3[['peaks_bin']] <- CreateAssayObject((marks$H3K4me3[['peaks']]@data > 0)*1)

marks$H3K27ac <- Pando::aggregate_assay(marks$H3K27ac, assay='peaks_bin', group_name='clusters')
marks$H3K27me3 <- Pando::aggregate_assay(marks$H3K27me3, assay='peaks_bin', group_name='clusters')
marks$H3K4me3 <- Pando::aggregate_assay(marks$H3K4me3, assay='peaks_bin', group_name='clusters')

marks$H3K27ac <- Pando::aggregate_assay(marks$H3K27ac, assay='peaks', group_name='clusters', slot='counts')
marks$H3K27me3 <- Pando::aggregate_assay(marks$H3K27me3, assay='peaks', group_name='clusters', slot='counts')
marks$H3K4me3 <- Pando::aggregate_assay(marks$H3K4me3, assay='peaks', group_name='clusters', slot='counts')

H3K27ac_cluster_agg <- marks$H3K27ac@assays$peaks_bin@misc$summary$clusters
H3K27ac_cluster_max <- colMaxs(H3K27ac_cluster_agg)
H3K27ac_cluster_num <- colSums2((H3K27ac_cluster_agg>0)*1)

H3K27me3_cluster_agg <- marks$H3K27me3@assays$peaks_bin@misc$summary$clusters
H3K27me3_cluster_max <- colMaxs(H3K27me3_cluster_agg)
H3K27me3_cluster_num <- colSums2((H3K27me3_cluster_agg>0)*1)

H3K4me3_cluster_agg <- marks$H3K4me3@assays$peaks_bin@misc$summary$clusters
H3K4me3_cluster_max <- colMaxs(H3K4me3_cluster_agg)
H3K4me3_cluster_num <- colSums2((H3K4me3_cluster_agg>0)*1)


min_clust <- 50
min_detect <- 0.1

H3K27ac_cluster_peaks <- H3K27ac_cluster_agg[,H3K27ac_cluster_max>min_detect & H3K27ac_cluster_num>=min_clust]
H3K27me3_cluster_peaks <- H3K27me3_cluster_agg[,H3K27me3_cluster_max>min_detect & H3K27me3_cluster_num>=min_clust]
H3K4me3_cluster_peaks <- H3K4me3_cluster_agg[,H3K4me3_cluster_max>min_detect & H3K4me3_cluster_num>=min_clust]

H3K27ac_cluster_counts <- marks$H3K27ac@assays$peaks@misc$summary$clusters[,H3K27ac_cluster_max>min_detect & H3K27ac_cluster_num>=min_clust]
H3K27me3_cluster_counts <- marks$H3K27me3@assays$peaks@misc$summary$clusters[,H3K27me3_cluster_max>min_detect & H3K27me3_cluster_num>=min_clust]
H3K4me3_cluster_counts <- marks$H3K4me3@assays$peaks@misc$summary$clusters[,H3K4me3_cluster_max>min_detect & H3K4me3_cluster_num>=min_clust]




#### Called cells per celltype ####
plot_df <- marks %>% map_dfr(function(x){
    rowSums(aggregate_matrix(t(x@assays$peaks_bin@data), groups=x$celltype_jf) > 0) %>% enframe('celltype', 'count')
}, .id = 'mark') %>% 
    mutate(celltype=factor(celltype, levels=ct_order)) %>% 
    filter(!is.na(celltype))

ggplot(plot_df, aes(celltype, count, fill=mark)) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=modality_colors) +
    article_text() + rotate_x_text(45) +
    labs(x='Cell type', y='# Called peaks', fill='Modality')
    
ggsave('plots/paper/fig_regs/peaks_per_ct.pdf', width=10, height=6, units='cm')

plot_df %>% write_tsv('data_/tables/fig3c_called_peaks.tsv')



#### Match all 3 marks ####
K27_isect_matches <- read_tsv('data_/intersect/H3K27_switches_marks_intersect_matches.tsv')
Kme_isect_matches <- read_tsv('data_/intersect/H3Kme_bivalent_marks_intersect_matches.tsv')
prom_isect_matches <- read_tsv('data_/intersect/H3K27K4_promoters_marks_intersect_matches.tsv')

H3K27ac_peak_ranges <- colnames(H3K27ac_cluster_peaks) %>% StringToGRanges()
H3K27me3_peak_ranges <- colnames(H3K27me3_cluster_peaks) %>% StringToGRanges()
H3K4me3_peak_ranges <- colnames(H3K4me3_cluster_peaks) %>% StringToGRanges()

K27_olaps <- findOverlaps(H3K27ac_peak_ranges, H3K27me3_peak_ranges)
Kme_olaps <- findOverlaps(H3K4me3_peak_ranges, H3K27me3_peak_ranges)
prom_olaps <- findOverlaps(H3K27ac_peak_ranges, H3K4me3_peak_ranges)

K27_olaps_tbl <- tibble(
    H3K27ac = colnames(H3K27ac_cluster_peaks)[queryHits(K27_olaps)],
    H3K27me3 = colnames(H3K27me3_cluster_peaks)[subjectHits(K27_olaps)]
)
K27_olaps_tbl$K27_isect <- GRangesToString(pintersect(StringToGRanges(K27_olaps_tbl$H3K27ac), StringToGRanges(K27_olaps_tbl$H3K27me3)))

Kme_olaps_tbl <- tibble(
    H3K4me3 = colnames(H3K4me3_cluster_peaks)[queryHits(Kme_olaps)],
    H3K27me3 = colnames(H3K27me3_cluster_peaks)[subjectHits(Kme_olaps)]
)
Kme_olaps_tbl$Kme_isect <- GRangesToString(pintersect(StringToGRanges(Kme_olaps_tbl$H3K4me3), StringToGRanges(Kme_olaps_tbl$H3K27me3)))

prom_olaps_tbl <- tibble(
    H3K27ac = colnames(H3K27ac_cluster_peaks)[queryHits(prom_olaps)],
    H3K4me3 = colnames(H3K4me3_cluster_peaks)[subjectHits(prom_olaps)]
)
prom_olaps_tbl$prom_isect <- GRangesToString(pintersect(StringToGRanges(prom_olaps_tbl$H3K27ac), StringToGRanges(prom_olaps_tbl$H3K4me3)))

all_olap_tbl <- K27_olaps_tbl %>% 
    full_join(Kme_olaps_tbl, multiple='all') %>% 
    full_join(prom_olaps_tbl, multiple='all') %>% 
    full_join(tibble(H3K27ac=colnames(H3K27ac_cluster_peaks))) %>% 
    full_join(tibble(H3K27me3=colnames(H3K27me3_cluster_peaks))) %>% 
    full_join(tibble(H3K4me3=colnames(H3K4me3_cluster_peaks))) %>% 
    mutate(region_status=case_when(
        !is.na(K27_isect) & !is.na(Kme_isect) & !is.na(prom_isect) ~ 'all',
        !is.na(K27_isect) ~ 'K27ac_K27me3',
        !is.na(Kme_isect) ~ 'K4me3_K27me3',
        !is.na(prom_isect) ~ 'K27ac_K4me3',
        !is.na(H3K27ac) ~ 'H3K27ac',
        !is.na(H3K27me3) ~ 'H3K27me3',
        !is.na(H3K4me3) ~ 'H3K4me3'
    )) %>% 
    mutate(repr_region=case_when(
        !is.na(K27_isect) ~ K27_isect,
        !is.na(Kme_isect) ~ Kme_isect,
        !is.na(prom_isect) ~ prom_isect,
        !is.na(H3K27ac) ~ H3K27ac,
        !is.na(H3K27me3) ~ H3K27me3,
        !is.na(H3K4me3) ~ H3K4me3
    ))

all_olap_tbl$region_name <- paste0(all_olap_tbl$region_status, '_', 1:nrow(all_olap_tbl))
all_olap_tbl <- select(all_olap_tbl, region_name, everything())

all_olap_tbl %>% write_tsv('data_/intersect/all_olap_matches.tsv')


# Make combined matrix
H3K27ac_common_peaks <- H3K27ac_cluster_peaks
colnames(H3K27ac_common_peaks) <- all_olap_tbl$region_name[match(colnames(H3K27ac_cluster_peaks), all_olap_tbl$H3K27ac)]

H3K27me3_common_peaks <- H3K27me3_cluster_peaks
colnames(H3K27me3_common_peaks) <- all_olap_tbl$region_name[match(colnames(H3K27me3_cluster_peaks), all_olap_tbl$H3K27me3)]

H3K4me3_common_peaks <- H3K4me3_cluster_peaks
colnames(H3K4me3_common_peaks) <- all_olap_tbl$region_name[match(colnames(H3K4me3_cluster_peaks), all_olap_tbl$H3K4me3)]

all_cluster_peaks <- rbind.fill.matrix(H3K27ac_common_peaks, H3K27me3_common_peaks, H3K4me3_common_peaks)
all_cluster_peaks[is.na(all_cluster_peaks)] <- 0
rownames(all_cluster_peaks) <- c(
    paste0('H3K27ac_', rownames(H3K27ac_common_peaks)), 
    paste0('H3K27me3_', rownames(H3K27me3_common_peaks)), 
    paste0('H3K4me3_', rownames(H3K4me3_common_peaks))
)
all_cluster_peaks <- Matrix::Matrix(all_cluster_peaks, sparse=T)


H3K27ac_common_counts <- H3K27ac_cluster_counts
colnames(H3K27ac_common_counts) <- all_olap_tbl$region_name[match(colnames(H3K27ac_cluster_counts), all_olap_tbl$H3K27ac)]

H3K27me3_common_counts <- H3K27me3_cluster_counts
colnames(H3K27me3_common_counts) <- all_olap_tbl$region_name[match(colnames(H3K27me3_cluster_counts), all_olap_tbl$H3K27me3)]

H3K4me3_common_counts <- H3K4me3_cluster_counts
colnames(H3K4me3_common_counts) <- all_olap_tbl$region_name[match(colnames(H3K4me3_cluster_counts), all_olap_tbl$H3K4me3)]

all_cluster_counts <- rbind.fill.matrix(H3K27ac_common_counts, H3K27me3_common_counts, H3K4me3_common_counts)
all_cluster_counts[is.na(all_cluster_counts)] <- 0
rownames(all_cluster_counts) <- c(
    paste0('H3K27ac_', rownames(H3K27ac_common_counts)), 
    paste0('H3K27me3_', rownames(H3K27me3_common_counts)), 
    paste0('H3K4me3_', rownames(H3K4me3_common_counts))
)
all_cluster_counts <- Matrix::Matrix(all_cluster_counts, sparse=T)



# Make seurat obj and preproc
all_peaks_srt <- CreateSeuratObject(scale(all_cluster_peaks))
print(all_peaks_srt)

all_peaks_srt <- all_peaks_srt %>% 
    FindVariableFeatures() %>% 
    ScaleData(vars.to.regress='nFeature_RNA') %>% 
    RunPCA()

print(all_peaks_srt)
ElbowPlot(all_peaks_srt)

all_peaks_srt <- all_peaks_srt %>% 
    RunUMAP(dims=1:20)

all_peaks_srt <- all_peaks_srt %>% FindNeighbors(dims=1:20)
all_peaks_srt <- all_peaks_srt %>% FindClusters()

dim_plot(all_peaks_srt, raster=F, label=T)

peak_annot <- ClosestFeature(marks$H3K27ac, all_olap_tbl$repr_region) %>% 
    as_tibble() %>% 
    inner_join(all_olap_tbl, by=c('query_region'='repr_region'), multiple='all') %>% 
    select(region_name, gene_name, 'repr_region'=query_region, distance) %>% distinct()

all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(peak_annot, 'region_name'))



#### Region specificity ####
all_olap_tbl <- read_tsv('data_/intersect/all_olap_matches.tsv')

H3K27ac_de_meta <- all_DA_df %>% 
    filter(mark=='H3K27ac', padj<1e-4 & detect_self>0.05 & coef>0) %>%
    dplyr::group_by(feature) %>% 
    dplyr::filter(log_dr==max(log_dr)) %>% 
    select('H3K27ac_group'=group, 'H3K27ac'=feature) %>% 
    inner_join(all_olap_tbl, multiple='all') %>% distinct(region_name, H3K27ac_group)

all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(H3K27ac_de_meta, 'region_name'))

H3K27me3_de_meta <- all_DA_df %>% 
    filter(mark=='H3K27me3', padj<1e-4 & detect_self>0.05 & coef>0) %>%
    dplyr::group_by(feature) %>% 
    dplyr::filter(log_dr==max(log_dr)) %>% 
    select('H3K27me3_group'=group, 'H3K27me3'=feature) %>% 
    inner_join(all_olap_tbl, multiple='all') %>% distinct(region_name, H3K27me3_group)

all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(H3K27me3_de_meta, 'region_name'))

H3K4me3_de_meta <- all_DA_df %>% 
    filter(mark=='H3K4me3', padj<1e-4 & detect_self>0.05 & coef>0) %>%
    dplyr::group_by(feature) %>% 
    dplyr::filter(log_dr==max(log_dr)) %>% 
    select('H3K4me3_group'=group, 'H3K4me3'=feature) %>% 
    inner_join(all_olap_tbl, multiple='all') %>% distinct(region_name, H3K4me3_group)

all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(H3K4me3_de_meta, 'region_name'))
all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(all_olap_tbl, 'region_name'))

promoter_peaks <- colnames(all_peaks_srt)[unique(subjectHits(findOverlaps(prom_ranges, StringToGRanges(all_peaks_srt$repr_region))))]

all_peaks_srt$log_clust <- log10(all_peaks_srt$nFeature_RNA)
all_peaks_srt$log_dist <- log10(all_peaks_srt$distance+1)
all_peaks_srt$chrom <- colnames(all_peaks_srt) %>% str_replace('(chr[\\d\\w]+)-.+', '\\1')
all_peaks_srt$max_detect <- colMaxs(all_cluster_peaks)
all_peaks_srt$distal <- all_peaks_srt$distance > 3000
all_peaks_srt$promoter <- colnames(all_peaks_srt) %in% promoter_peaks

all_peaks_srt$ctx_specific_mark <- case_when(
    (all_peaks_srt$H3K27ac_group == 'ctx') ~ 'H3K27ac',
    (all_peaks_srt$H3K27me3_group == 'ctx') ~ 'H3K27me3',
    (all_peaks_srt$H3K4me3_group == 'ctx') ~ 'H3K4me3',
)

all_peaks_srt$dien_specific_mark <- case_when(
    (all_peaks_srt$H3K27ac_group == 'dien') ~ 'H3K27ac',
    (all_peaks_srt$H3K27me3_group == 'dien') ~ 'H3K27me3',
    (all_peaks_srt$H3K4me3_group == 'dien') ~ 'H3K4me3',
)

all_peaks_srt$nt_specific_mark <- case_when(
    (all_peaks_srt$H3K27ac_group == 'nt') ~ 'H3K27ac',
    (all_peaks_srt$H3K27me3_group == 'nt') ~ 'H3K27me3',
    (all_peaks_srt$H3K4me3_group == 'nt') ~ 'H3K4me3',
)

all_peaks_srt$retina_specific_mark <- case_when(
    (all_peaks_srt$H3K27ac_group == 'retina') ~ 'H3K27ac',
    (all_peaks_srt$H3K27me3_group == 'retina') ~ 'H3K27me3',
    (all_peaks_srt$H3K4me3_group == 'retina') ~ 'H3K4me3',
)





#### Age specificity ####
H3K27me3_age_DA <- read_tsv('data_/results/diff_expression/H3K27me3_DA_peaks_age.tsv') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    filter(padj<1e-4 & detect_self>0.05 & coef>0) %>%     
    dplyr::group_by(feature) %>% 
    dplyr::filter(coef==max(coef)) %>% 
    select('H3K27me3_age_group'=group, 'H3K27me3'=feature) %>% 
    inner_join(all_olap_tbl, multiple='all') %>% distinct(region_name, H3K27me3_age_group)

H3K27ac_age_DA <- read_tsv('data_/results/diff_expression/H3K27ac_DA_peaks_age.tsv') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    filter(padj<1e-4 & detect_self>0.05 & coef>0) %>%     
    dplyr::group_by(feature) %>% 
    dplyr::filter(coef==max(coef)) %>% 
    select('H3K27ac_age_group'=group, 'H3K27ac'=feature) %>% 
    inner_join(all_olap_tbl, multiple='all') %>% distinct(region_name, H3K27ac_age_group)

H3K4me3_age_DA <- read_tsv('data_/results/diff_expression/H3K4me3_DA_peaks_age.tsv') %>% 
    mutate(padj=p.adjust(pval, method='fdr')) %>% 
    filter(padj<1e-4 & detect_self>0.05 & coef>0) %>%     
    dplyr::group_by(feature) %>% 
    dplyr::filter(coef==max(coef)) %>% 
    select('H3K4me3_age_group'=group, 'H3K4me3'=feature) %>% 
    inner_join(all_olap_tbl, multiple='all') %>% distinct(region_name, H3K4me3_age_group)

all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(H3K27me3_age_DA, 'region_name'))
all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(H3K27ac_age_DA, 'region_name'))
all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(H3K4me3_age_DA, 'region_name'))


all_peaks_srt$nepi_specific_mark <- case_when(
    (all_peaks_srt$H3K27ac_age_group == '15d') ~ 'H3K27ac',
    (all_peaks_srt$H3K27me3_age_group == '15d') ~ 'H3K27me3',
    (all_peaks_srt$H3K4me3_age_group == '15d') ~ 'H3K4me3',
)

all_peaks_srt$ipsc_specific_mark <- case_when(
    (all_peaks_srt$H3K27ac_age_group == 'EB') ~ 'H3K27ac',
    (all_peaks_srt$H3K27me3_age_group == 'EB') ~ 'H3K27me3',
    (all_peaks_srt$H3K4me3_age_group == 'EB') ~ 'H3K4me3',
)




#### Plot stuff ####
age_names <- names(age_colors2)
age_colors <-c("#6FC4EE", "#2C86BD", "#005E97", "#6065AE", "#6065AE", "#6065AE", "#47439B", "#253773", "#F1BF66", "#DC9C29")
names(age_colors) <- age_names

modality_colors['all'] <- '#546E7A'
modality_colors['H3K27me3'] <- '#027BB8'
modality_colors['H3K4me3'] <- '#E56F8E'
modality_colors['H3K27ac'] <- '#7CAF9C'
modality_colors['K27ac_K4me3'] <- '#E58B4A'
modality_colors['K4me3_K27me3'] <- '#9C6FF7'
modality_colors['K27ac_K27me3'] <- '#0097A7'

p1 <- dim_plot(all_peaks_srt, group.by='region_status', order=T, pt.size=0.5) +
    scale_color_manual(values=modality_colors)

p2 <- dim_plot(all_peaks_srt, group.by=c('promoter', 'distal'), order=T, pt.size=0.1) &
    scale_color_manual(values=c('grey', 'black'))

p1 / p2 + plot_layout(heights=c(2,1))
ggsave('plots/paper/fig_regs/region_status_umap.png', width=8, height=8)
ggsave('plots/paper/fig_regs/region_status_umap.pdf', width=8, height=8)


lineage_colors2['nt'] <- lineage_colors2['mesen_rhom']
p1 <- dim_plot(all_peaks_srt, group.by=c('H3K4me3_group', 'H3K27ac_group', 'H3K27me3_group'), order=T, pt.size=0.4) &
    scale_color_manual(values=lineage_colors2, na.value='grey')

p2 <- dim_plot(all_peaks_srt, group.by=c('H3K4me3_age_group', 'H3K27ac_age_group', 'H3K27me3_age_group'), order=T, pt.size=0.4) &
    scale_color_manual(values=age_colors, na.value='grey')

p1 / p2
ggsave('plots/paper/fig_regs/region_groups_umap.png', width=18, height=10)
ggsave('plots/paper/fig_regs/region_groups_umap.pdf', width=18, height=10)


dim_plot(all_peaks_srt, group.by=c('ctx_specific_mark', 'dien_specific_mark', 'nt_specific_mark', 'retina_specific_mark'), order=T, pt.size=0.5) &
    scale_color_manual(values=modality_colors, na.value='grey')
ggsave('plots/paper/fig_regs/region_region_marks_umap.png', width=12, height=8)
ggsave('plots/paper/fig_regs/region_region_marks_umap.pdf', width=12, height=8)


feature_plot(all_peaks_srt, features=c('max_detect', 'log_dist', 'log_clust'), order=T) &
    theme(legend.position = 'right')
ggsave('plots/paper/fig_regs/region_features_umap.png', width=12, height=8)
ggsave('plots/paper/fig_regs/region_features_umap.pdf', width=12, height=8)


dim_plot(all_peaks_srt, label=T)
ggsave('plots/paper/fig_regs/region_clusters_umap.png', width=8, height=5)
ggsave('plots/paper/fig_regs/region_clusters_umap.pdf', width=8, height=5)

dim_plot(all_peaks_srt)
ggsave('plots/paper/fig_regs/region_clusters_nolabs_umap.png', width=8, height=5)
ggsave('plots/paper/fig_regs/region_clusters_nolabs_umap.pdf', width=8, height=5)


peak_meta <- all_peaks_srt@meta.data %>% 
    as_tibble(rownames='region') %>% 
    mutate(region_dist=factor(case_when(
        promoter ~ 'promoter',
        distal ~'distal',
        T ~ 'gene_body'
    ), levels=c('distal', 'gene_body', 'promoter')))

ggplot(peak_meta, aes(region_status, fill=region_status)) +
    geom_bar() +
    scale_fill_manual(values=modality_colors) +
    article_text() +
    theme_rangeframe() + scale_axis_rangeframe() +
    rotate_x_text(40)
ggsave('plots/paper/fig_regs/region_status_counts_bar.pdf', width=7, height=4, units='cm')


ggplot(peak_meta, aes(region_status, fill=region_dist)) +
    geom_bar(linewidth=0.2, color='black') + 
    scale_fill_manual(values=c('white', 'grey', 'black')) +
    article_text() +
    theme_rangeframe() + scale_axis_rangeframe() +
    rotate_x_text(40)
ggsave('plots/paper/fig_regs/region_promoter_counts_bar.pdf', width=7, height=4, units='cm')


# Barplots with region-specificity per cluster
plot_df <- peak_meta %>% 
    pivot_longer(cols=c('H3K4me3_group', 'H3K27me3_group', 'H3K27ac_group'))
    
ggplot(plot_df, aes(seurat_clusters, fill=value)) +
    geom_bar() +
    facet_grid(~name)

plot_df <- peak_meta %>% 
    pivot_longer(cols=c('H3K4me3_group', 'H3K27me3_group', 'H3K27ac_group')) %>% 
    mutate(name=str_remove(name, '_group')) %>% 
    filter(!is.na(value))
    
ggplot(plot_df, aes(seurat_clusters, fill=value)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=lineage_colors2) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(~name) +
    article_text() +
    labs(x='Louvain clusters', y='Fraction', fill='Branch')

ggsave('plots/paper/fig_regs/region_specificity_cluster_bar.pdf', width=12, height=4, units='cm')
ggsave('plots/paper/fig_regs/region_specificity_cluster_bar.png', width=12, height=4, units='cm')



plot_df <- peak_meta %>% 
    pivot_longer(cols=c('H3K4me3_group', 'H3K27me3_group', 'H3K27ac_group')) %>% 
    mutate(name=str_remove(name, '_group'))
    
ggplot(plot_df, aes(seurat_clusters, fill=!is.na(value))) +
    geom_bar(position='fill') +
    scale_fill_manual(values=c('grey', 'black')) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(~name) +
    article_text() +
    labs(x='Louvain clusters', y='Fraction', fill='Region-specific')

ggsave('plots/paper/fig_regs/region_spec_frac_cluster_bar.pdf', width=12, height=4, units='cm')
ggsave('plots/paper/fig_regs/region_spec_frac_cluster_bar.png', width=12, height=4, units='cm')


ggplot(plot_df, aes(seurat_clusters, fill=!is.na(value))) +
    geom_bar(color='black', linewidth=0.2) +
    scale_fill_manual(values=c('white', 'black')) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(~name) +
    article_text() +
    labs(x='Louvain clusters', y='Count', fill='Region-specific')

ggsave('plots/paper/fig_regs/region_spec_counts_cluster_bar.pdf', width=12, height=4, units='cm')
ggsave('plots/paper/fig_regs/region_spec_counts_cluster_bar.png', width=12, height=4, units='cm')



ggplot(plot_df, aes(name, fill=!is.na(value))) +
    geom_bar(color='black', linewidth=0.2, position='fill') +
    scale_fill_manual(values=c('white', 'black')) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(~seurat_clusters) +
    article_text() + rotate_x_text(40) +
    theme(
        panel.spacing.x = unit(0.1, 'lines')
    ) +
    labs(x='Louvain clusters', y='Fraction', fill='Region-specific')

ggsave('plots/paper/fig_regs/region_spec_frac_cluster_marks_bar.pdf', width=15, height=4, units='cm')
ggsave('plots/paper/fig_regs/region_spec_frac_cluster_marks_bar.png', width=15, height=4, units='cm')



ggplot(plot_df, aes(seurat_clusters, fill=!is.na(value))) +
    geom_bar(color='black', linewidth=0.2) +
    scale_fill_manual(values=c('white', 'black')) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    facet_grid(~name) +
    article_text() +
    labs(x='Louvain clusters', y='Count', fill='Region-specific')

ggsave('plots/paper/fig_regs/region_spec_counts_cluster_bar.pdf', width=12, height=4, units='cm')
ggsave('plots/paper/fig_regs/region_spec_counts_cluster_bar.png', width=12, height=4, units='cm')



plot_df <- peak_meta %>% 
    mutate(lin_specific=!is.na(H3K27ac_group) | !is.na(H3K4me3_group) | !is.na(H3K27me3_group))

ggplot(plot_df, aes(seurat_clusters, fill=lin_specific)) +
    geom_bar(color='black', linewidth=0.2, position='fill') +
    scale_fill_manual(values=c('white', 'black')) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    article_text() +
    labs(x='Louvain clusters', y='Fraction of branch-specific\nregulatory regions', fill='Region-specific')

ggsave('plots/paper/fig_regs/region_spec_frac_all_cluster_bar.pdf', width=8, height=4, units='cm')
ggsave('plots/paper/fig_regs/region_spec_frac_all_cluster_bar.png', width=8, height=4, units='cm')




# Save 
# all_peaks_srt %>% write_rds('data05/regions/all_regions_v2_max01_clust50_srt.rds')
# all_peaks_srt %>% write_rds('data_/all_peaks/all_regions_v2_max01_clust50_srt.rds')
all_peaks_srt <- read_rds('data05/regions/all_regions_v2_max01_clust50_srt.rds')
all_peaks_srt$repr_region %>% StringToGRanges() %>% export.bed('data_/regions/all_regions_v2_max01_clust50.bed')

#### Region groups ####
# - region promoters per region
# - region distal per region
# - only repressed / only active clusters


# * ctx-H3K27ac-distal
# * ctx-H3K27ac-promoter
# * ctx-H3K27me3-distal
# * ctx-H3K27me3-promoter
# * ctx-H3K4me3-distal
# * ctx-H3K4me3-promoter
# -> for each region ...

dim_plot(all_peaks_srt, label=T)

peak_meta <- all_peaks_srt@meta.data %>% as_tibble(rownames='region')

H3K27ac_groups <- peak_meta %>% 
    filter(!is.na(H3K27ac_group)) %>% 
    mutate(dist_group=case_when(
        promoter ~ 'promoter',
        distal ~ 'distal',
        distance==0 ~ 'gene_body'
    )) %>% filter(!is.na(dist_group)) %>% 
    group_by(H3K27ac_group, dist_group) %>% group_split() %>% 
    {names(.) <- map_chr(., function(.x){paste0('H3K27ac_', .x$H3K27ac_group[1], '_', .x$dist_group[1])});.} %>% 
    map(~.x$H3K27ac)

H3K27me3_groups <- peak_meta %>% 
    filter(!is.na(H3K27me3_group)) %>% 
    mutate(dist_group=case_when(
        promoter ~ 'promoter',
        distal ~ 'distal',
        distance==0 ~ 'gene_body'
    )) %>% filter(!is.na(dist_group)) %>% 
    group_by(H3K27me3_group, dist_group) %>% group_split() %>% 
    {names(.) <- map_chr(., function(.x){paste0('H3K27me3_', .x$H3K27me3_group[1], '_', .x$dist_group[1])});.} %>% 
    map(~.x$H3K27me3)

H3K4me3_groups <- peak_meta %>% 
    filter(!is.na(H3K4me3_group)) %>% 
    mutate(dist_group=case_when(
        promoter ~ 'promoter',
        distal ~ 'distal',
        distance==0 ~ 'gene_body'
    )) %>% filter(!is.na(dist_group)) %>% 
    group_by(H3K4me3_group, dist_group) %>% group_split() %>% 
    {names(.) <- map_chr(., function(.x){paste0('H3K4me3_', .x$H3K4me3_group[1], '_', .x$dist_group[1])});.} %>% 
    map(~.x$H3K4me3)

repressed_groups <- peak_meta %>% 
    filter(region_status=='H3K27me3', RNA_snn_res.0.8%in%c(7,11,4)) %>% 
    group_by(RNA_snn_res.0.8) %>% group_split() %>%     
    {names(.) <- map_chr(., function(.x){paste0('H3K27me3_repressed_', .x$RNA_snn_res.0.8[1])});.} %>% 
    map(~.x$H3K27me3)

active_groups <- peak_meta %>% 
    filter(region_status=='H3K27ac', RNA_snn_res.0.8%in%c(10,12,6,5,0,9)) %>% 
    group_by(RNA_snn_res.0.8) %>% group_split() %>%     
    {names(.) <- map_chr(., function(.x){paste0('H3K27ac_active_', .x$RNA_snn_res.0.8[1])});.} %>% 
    map(~.x$H3K27ac)

region_groups <- c(H3K27ac_groups, H3K27me3_groups, H3K4me3_groups, repressed_groups, active_groups)
H3K27ac_region_groups <- c(H3K27ac_groups, active_groups)
H3K27me3_region_groups <- c(H3K27me3_groups, repressed_groups)
H3K4me3_region_groups <- H3K4me3_groups

all_cluster_groups <- peak_meta %>% 
    group_by(RNA_snn_res.0.8) %>% group_split() %>%  
    {names(.) <- map_chr(., function(.x){as.character(.x$RNA_snn_res.0.8)[1]});.} %>% 
    map(~.x$repr_region)

all_status_groups <- peak_meta %>% 
    group_by(region_status) %>% group_split() %>%  
    {names(.) <- map_chr(., function(.x){as.character(.x$region_status)[1]});.} %>% 
    map(~.x$repr_region)



#### TF motif enrichment ####
library(doParallel)
registerDoParallel(40)

H3K27ac_motif_mat <- marks$H3K27ac@assays$peaks@motifs@data
H3K27ac_background <- colnames(H3K27ac_cluster_peaks)

H3K27me3_motif_mat <- marks$H3K27me3@assays$peaks@motifs@data
H3K27me3_background <- colnames(H3K27me3_cluster_peaks)

H3K4me3_motif_mat <- marks$H3K4me3@assays$peaks@motifs@data
H3K4me3_background <- colnames(H3K4me3_cluster_peaks)


H3K27ac_tf_enrich <- map_dfr(set_names(names(H3K27ac_region_groups)), function(x){
    print(x)
    reg <- H3K27ac_region_groups[[x]]
    tf_motif_enrichment(reg, H3K27ac_motif_mat, H3K27ac_background, parallel=T)
}, .id='region_group')

H3K27ac_tf_enrich %>% write_tsv('data_/all_peaks/H3K27ac_peak_groups_tf_enrich.tsv')


H3K27me3_tf_enrich <- map_dfr(set_names(names(H3K27me3_region_groups)), function(x){
    print(x)
    reg <- H3K27me3_region_groups[[x]]
    tf_motif_enrichment(reg, H3K27me3_motif_mat, H3K27me3_background, parallel=T)
}, .id='region_group')

H3K27me3_tf_enrich %>% write_tsv('data_/all_peaks/H3K27me3_peak_groups_tf_enrich.tsv')


H3K4me3_tf_enrich <- map_dfr(set_names(names(H3K4me3_region_groups)), function(x){
    print(x)
    reg <- H3K4me3_region_groups[[x]]
    tf_motif_enrichment(reg, H3K4me3_motif_mat, H3K4me3_background, parallel=T)
}, .id='region_group')

H3K4me3_tf_enrich %>% write_tsv('data_/all_peaks/H3K4me3_peak_groups_tf_enrich.tsv')




H3K27ac_tf_enrich <- map_dfr(set_names(names(H3K27ac_region_groups)), function(x){
    print(x)
    reg <- H3K27ac_region_groups[[x]]
    enrich_df <- as_tibble(FindMotifs(marks$H3K27ac, reg)) %>% inner_join(motif2tf)
    return(enrich_df)
}, .id='region_group')

H3K27ac_tf_enrich %>% write_tsv('data_/all_peaks/H3K27ac_peak_groups_tf_enrich_signac.tsv')


H3K27me3_tf_enrich <- map_dfr(set_names(names(H3K27me3_region_groups)), function(x){
    print(x)
    reg <- H3K27me3_region_groups[[x]]
    enrich_df <- as_tibble(FindMotifs(marks$H3K27me3, reg)) %>% inner_join(motif2tf)
    return(enrich_df)
}, .id='region_group')

H3K27me3_tf_enrich %>% write_tsv('data_/all_peaks/H3K27me3_peak_groups_tf_enrich_signac.tsv')


H3K4me3_tf_enrich <- map_dfr(set_names(names(H3K4me3_region_groups)), function(x){
    print(x)
    reg <- H3K4me3_region_groups[[x]]
    enrich_df <- as_tibble(FindMotifs(marks$H3K4me3, reg)) %>% inner_join(motif2tf)
    return(enrich_df)
}, .id='region_group')

H3K4me3_tf_enrich %>% write_tsv('data_/all_peaks/H3K4me3_peak_groups_tf_enrich_signac.tsv')






H3K27ac_tf_enrich %>% filter(padj<0.01, origin=='JASPAR2020', logodds>1, !is.infinite(logodds)) %>% filter(region_group=='H3K27ac_ctx_promoter') %>% 
    arrange(dplyr::desc(logodds)) %>% pull(tf)
H3K27ac_tf_enrich %>% filter(padj<0.01, origin=='JASPAR2020', logodds>1) %>% filter(region_group=='H3K27ac_active_13')


H3K27ac_tf_enrich <- read_tsv('data_/all_peaks/H3K27ac_peak_groups_tf_enrich.tsv')
H3K27me3_tf_enrich <- read_tsv('data_/all_peaks/H3K27me3_peak_groups_tf_enrich.tsv')
H3K4me3_tf_enrich <- read_tsv('data_/all_peaks/H3K4me3_peak_groups_tf_enrich.tsv')







all_bg <- peak_meta$repr_region %>% unique() %>% StringToGRanges()
all_motif_match <- motifmatchr::matchMotifs(motifs, all_bg, genome=BSgenome.Hsapiens.UCSC.hg38)
all_motif_mat <- all_motif_match@assays@data$motifMatches
rownames(all_motif_mat) <- GRangesToString(all_bg)
colnames(all_motif_mat) <- names(motifs)
all_motif_mat %>% write_rds('data_/all_peaks/all_motifs.rds')


all_motif_mat <- read_rds('data_/all_peaks/all_motifs.rds')
cluster_tf_enrich <- map_dfr(set_names(names(all_cluster_groups)), function(x){
    print(x)
    reg <- all_cluster_groups[[x]]
    tf_motif_enrichment(reg, all_motif_mat, GRangesToString(all_bg), parallel=T)
}, .id='region_group')

cluster_tf_enrich %>% rename('cluster'='region_group') %>% write_tsv('data_/all_peaks/clusters_tf_enrich.tsv')



library(doParallel)
registerDoParallel(40)

# Switching / bival
bival_peaks <- read_tsv('data_/CT/bivalent_peaks.tsv')$H3K4me3 %>% unique()
bival_bg <- H3K4me3_cluster_peaks %>% colnames()
bivalent_tf_enrich <- tf_motif_enrichment(bival_peaks, Motifs(marks$H3K4me3)@data, bival_bg, parallel=T)
bivalent_tf_enrich %>% write_tsv('data_/all_peaks/bivalent_tf_enrich.tsv')


switch_peaks <- read_tsv('data_/CT/switching_peaks.tsv')$H3K27ac %>% unique() 
switch_bg <- H3K27ac_cluster_peaks %>% colnames()
switching_tf_enrich <- tf_motif_enrichment(switch_peaks, Motifs(marks$H3K27ac)@data, switch_bg, parallel=T)
switching_tf_enrich %>% write_tsv('data_/all_peaks/switching_tf_enrich.tsv')


all_bg <- peak_meta$repr_region %>% unique() %>% StringToGRanges()
status_tf_enrich <- map_dfr(set_names(names(all_status_groups)), function(x){
    print(x)
    reg <- all_status_groups[[x]]
    tf_motif_enrichment(reg, all_motif_mat, GRangesToString(all_bg), parallel=T)
}, .id='status')

status_tf_enrich %>% write_tsv('data_/all_peaks/status_tf_enrich.tsv')




#### GREAT enrichment ####
H3K27ac_bg <- StringToGRanges(H3K27ac_background)
H3K27ac_great_enrich <- map_dfr(set_names(names(H3K27ac_region_groups)), function(x){
    print(x)
    reg <- H3K27ac_region_groups[[x]]
    res = great(reg, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", background=H3K27ac_bg)
    res@table %>% as_tibble() %>% return()
}, .id='region_group')

H3K27me3_bg <- StringToGRanges(H3K27me3_background)
H3K27me3_great_enrich <- map_dfr(set_names(names(H3K27me3_region_groups)), function(x){
    print(x)
    reg <- H3K27me3_region_groups[[x]]
    res = great(reg, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", background=H3K27me3_bg)
    res@table %>% as_tibble() %>% return()
}, .id='region_group')

H3K4me3_bg <- StringToGRanges(H3K4me3_background)
H3K4me3_great_enrich <- map_dfr(set_names(names(H3K4me3_region_groups)), function(x){
    print(x)
    reg <- H3K4me3_region_groups[[x]]
    res = great(reg, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", background=H3K4me3_bg)
    res@table %>% as_tibble() %>% return()
}, .id='region_group')


H3K27ac_great_enrich %>% write_tsv('data_/all_peaks/H3K27ac_peak_groups_great_enrich.tsv')
H3K27me3_great_enrich %>% write_tsv('data_/all_peaks/H3K27me3_peak_groups_great_enrich.tsv')
H3K4me3_great_enrich %>% write_tsv('data_/all_peaks/H3K4me3_peak_groups_great_enrich.tsv')



all_bg <- peak_meta$repr_region %>% unique() %>% StringToGRanges()
cluster_great_enrich <- map_dfr(set_names(names(all_cluster_groups)), function(x){
    print(x)
    reg <- StringToGRanges(all_cluster_groups[[x]])
    res <- great(reg, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", background=all_bg)
    res@table %>% as_tibble() %>% return()
}, .id='region_group')

cluster_great_enrich %>% rename('cluster'='region_group') %>% write_tsv('data_/all_peaks/clusters_great_enrich.tsv')

status_great_enrich <- map_dfr(set_names(names(all_status_groups)), function(x){
    print(x)
    reg <- StringToGRanges(all_status_groups[[x]])
    res <- great(reg, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", background=all_bg)
    res@table %>% as_tibble() %>% return()
}, .id='status')

status_great_enrich %>% write_tsv('data_/all_peaks/status_great_enrich.tsv')






#### Compare with known peaks ####
library(rtracklayer)
data("SCREEN.ccRE.UCSC.hg38")

ch = import.chain('~/resources/hg19ToHg38.over.chain')

all_vista_enhancers <- read_tsv('~/resources/vista_enhancers.tsv', col_names=F)
fb_vista_enhancers <- read_tsv('~/resources/fb_vista_enhancers.tsv', col_names=F)
hb_vista_enhancers <- read_tsv('~/resources/hb_vista_enhancers.tsv', col_names=F)
mb_vista_enhancers <- read_tsv('~/resources/mb_vista_enhancers.tsv', col_names=F)
nt_vista_enhancers <- read_tsv('~/resources/nt_vista_enhancers.tsv', col_names=F)

vista_list <- list(
    'all' = all_vista_enhancers,
    'forebrain' = fb_vista_enhancers,
    'hindbrain' = hb_vista_enhancers,
    'midbrain' = mb_vista_enhancers,
    'neuraltube' = nt_vista_enhancers
)

vista_ranges <- vista_list %>% map_par(function(x){
    x$X1 %>% 
        StringToGRanges(sep=c(':','-')) %>% 
        liftOver(ch) %>% 
        purrr::reduce(c) %>% 
        return()
}, parallel=T)

vista_chr <- map(vista_ranges, ~GRangesToString(.x))

vista_meta <- tibble(region=vista_chr$all)
vista_meta$all <- vista_meta$region %in% vista_chr$all
vista_meta$forebrain <- vista_meta$region %in% vista_chr$forebrain
vista_meta$midbrain <- vista_meta$region %in% vista_chr$midbrain
vista_meta$hindbrain <- vista_meta$region %in% vista_chr$hindbrain
vista_meta$neuraltube <- vista_meta$region %in% vista_chr$neuraltube

olaps <- StringToGRanges(vista_meta$region) %>% findOverlaps(all_bg) %>% queryHits()

vista_meta$is_covered <- vista_meta$region %in% vista_meta$region[olaps]

plot_df <- vista_meta %>% 
    pivot_longer(forebrain:neuraltube) %>% 
    filter(value) %>% mutate(name=factor(name, levels=c('neuraltube', 'forebrain', 'midbrain', 'hindbrain')))

ggplot(plot_df, aes(name, fill=is_covered)) +
    geom_bar(position='fill', color='black', size=0.2) +
    scale_fill_manual(values=c('white', '#555555')) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() + rotate_x_text(45) +
    labs(x='Brain region', y='Fraction of enhancers', fill='Is covered')
ggsave('plots/paper/fig_regs/region_vista_fraction_bar.pdf', width=5, height=4, units='cm')

seqlevelsStyle(all_bg)<- 'NCBI'
peek_annot <- ChIPseeker::annotatePeak(all_bg, tssRegion=c(2000, 200), TxDb=EnsDb.Hsapiens.v86)
seqlevels(peek_annot@anno) <- paste0('chr', seqlevels(peek_annot@anno))
peek_annot_df <- peek_annot@anno %>% 
    as_tibble() %>% mutate(repr_region=GRangesToString(peek_annot@anno))

peak_meta_annot <- peak_meta %>% 
    left_join(peek_annot_df) %>% 
    mutate(
        annot_short = case_when(
            str_detect(annotation, 'Intron') ~ 'Intron',
            str_detect(annotation, 'Exon') ~ 'Exon',
            str_detect(annotation, 'Promoter') ~ 'Promoter',
            .default = annotation
        )
    ) %>% 
    mutate(dist_group=case_when(
        promoter ~ 'promoter',
        distal ~ 'distal',
        distance==0 ~ 'gene_body'
    ))

plot_df <- peak_meta_annot
ggplot(plot_df, aes(annot_short, fill=region_status)) +
    geom_bar() +
    scale_fill_manual(values=modality_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    facet_grid(~region_status) +
    article_text() + rotate_x_text(45) +
    labs(x='Annotation', y='Count', fill='Region status')
ggsave('plots/paper/fig_regs/region_status_annot_bar.pdf', width=15, height=6, units='cm')



plot_df <- peak_meta_annot %>% 
    filter(!is.na(H3K27ac))
ggplot(plot_df, aes(annot_short)) +
    geom_bar(fill=modality_colors['H3K27ac']) +
    scale_fill_manual(values=modality_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() + rotate_x_text(45) +
    labs(x='Annotation', y='Count', fill='Region status')
ggsave('plots/paper/fig_regs/H3K27ac_region_annot_bar.pdf', width=7, height=6, units='cm')


plot_df <- peak_meta_annot %>% 
    filter(!is.na(H3K27me3))
ggplot(plot_df, aes(annot)) +
    geom_bar(fill=modality_colors['H3K27me3']) +
    scale_fill_manual(values=modality_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() + rotate_x_text(45) +
    labs(x='Annotation', y='Count', fill='Region status')
ggsave('plots/paper/fig_regs/H3K27me3_region_annot_bar.pdf', width=7, height=6, units='cm')


plot_df <- peak_meta_annot %>% 
    filter(!is.na(H3K4me3))
ggplot(plot_df, aes(annot)) +
    geom_bar(fill=modality_colors['H3K4me3']) +
    scale_fill_manual(values=modality_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() + rotate_x_text(45) +
    labs(x='Annotation', y='Count', fill='Region status')
ggsave('plots/paper/fig_regs/H3K4me3_region_annot_bar.pdf', width=7, height=6, units='cm')

peak_meta_annot %>% mutate(dist_group=case_when(
    promoter ~ 'promoter',
    distal ~ 'distal',
    distance==0 ~ 'gene_body',
)) %>% write_tsv('data_/all_peaks/region_meta_annot.tsv')


#### Check distal/promoter for region-specific peaks ####
peak_meta_annot <- read_tsv('data_/all_peaks/region_meta_annot.tsv')

all_peaks_srt <- AddMetaData(all_peaks_srt, column_to_rownames(peak_meta_annot, 'region'))
all_peaks_srt$annot_short
all_peaks_srt$dist_group

plot_df <- peak_meta_annot %>% 
    filter(!is.na(H3K27ac_group), !is.na(dist_group))

p1 <- ggplot(plot_df, aes(H3K27ac_group, fill=dist_group)) +
    geom_bar(color='black', position='fill') +
    scale_fill_manual(values=c('white', 'grey', 'black')) 

plot_df <- peak_meta_annot %>% 
    filter(!is.na(H3K27me3_group), !is.na(dist_group))    

p2 <- ggplot(plot_df, aes(H3K27me3_group, fill=dist_group)) +
    geom_bar(color='black', position='fill') +
    scale_fill_manual(values=c('white', 'grey', 'black')) 

plot_df <- peak_meta_annot %>% 
    filter(!is.na(H3K4me3_group), !is.na(dist_group))    

p3 <- ggplot(plot_df, aes(H3K4me3_group, fill=dist_group)) +
    geom_bar(color='black', position='fill') +
    scale_fill_manual(values=c('white', 'grey', 'black')) 
    
p1 | p2 | p3
    


#### All other TF motif enrichments now also here ####
#### Enrichment in merged regions ####
library(doParallel)
registerDoParallel(36)

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
genes_primed <- gene_ranges_use[gene_ranges_use$gene_name%in%primed_group]
genes_mid <- gene_ranges_use[gene_ranges_use$gene_name%in%mid_group]

primed2peaks <- find_peaks_near_genes(all_bg, genes_primed, distance=4000)
primed_peaks <- rownames(primed2peaks)[rowMaxs(primed2peaks)>0] 
mid2peaks <- find_peaks_near_genes(all_bg, genes_mid, distance=4000)
mid_peaks <- rownames(mid2peaks)[rowMaxs(mid2peaks)>0] 

primed_tf_enrich <- tf_motif_enrichment(primed_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
primed_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/primed_tf_enrich.tsv')

mid_tf_enrich <- tf_motif_enrichment(mid_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
mid_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/mid_tf_enrich.tsv')

primed2peaks <- find_peaks_near_genes(all_bg, genes_primed, distance=4000, only_tss=T)
primed_peaks <- rownames(primed2peaks)[rowMaxs(primed2peaks)>0] 
mid2peaks <- find_peaks_near_genes(all_bg, genes_mid, distance=4000, only_tss=T)
mid_peaks <- rownames(mid2peaks)[rowMaxs(mid2peaks)>0] 

primed_tf_enrich <- tf_motif_enrichment(primed_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
primed_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/primed_tf_promoter_enrich.tsv')

mid_tf_enrich <- tf_motif_enrichment(mid_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
mid_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/mid_tf_promoter_enrich.tsv')



# For the neuron gene cluster
rna_ct_de <- read_tsv('data_/results/diff_expression/RNA_DE_celltype.tsv')
gene_clusters <- read_tsv('data_/trajectories/ctx/all_mod_pseudotime_genes_expr_dtw_30clust.tsv')
neuron_genes <- gene_clusters %>% filter(dtw_clust==27) %>% pull(feature) %>% unique()

neuron_gene_ranges <- gene_ranges_use[gene_ranges_use$gene_name%in%neuron_genes]
neuron2peaks <- find_peaks_near_genes(all_bg, neuron_gene_ranges, distance=4000, only_tss=T)
neuron_peaks <- rownames(neuron2peaks)[rowMaxs(neuron2peaks)>0] 

neuron_tf_enrich <- tf_motif_enrichment(neuron_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
neuron_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/neuron_cluster_tf_promoter_enrich.tsv')


neuron2peaks <- find_peaks_near_genes(all_bg, neuron_gene_ranges, distance=4000)
neuron_peaks <- rownames(neuron2peaks)[rowMaxs(neuron2peaks)>0] 

neuron_tf_enrich <- tf_motif_enrichment(neuron_peaks, all_motif_mat, GRangesToString(all_bg), parallel=T)
neuron_tf_enrich %>% write_tsv('data_/trajectories/ctx/enrichment/neuron_cluster_tf_enrich.tsv')


neuron_prom_enrich <- read_tsv('data_/trajectories/ctx/enrichment/primed_tf_promoter_enrich.tsv') %>% 
    filter(origin=='JASPAR2020')

neuron_enrich <- read_tsv('data_/trajectories/ctx/enrichment/primed_tf_enrich.tsv') %>% 
    filter(origin=='JASPAR2020')

ctx_de <- rna_ct_de %>% filter(cluster=='ctx_ex')


plot_df <- inner_join(neuron_enrich, ctx_de, by=c('tf'='gene'))

ggplot(plot_df, aes(logodds, avg_log2FC, label=tf, alpha=padj<0.05)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    theme_rangeframe() + scale_axis_rangeframe() +
    geom_text(size=4/ggplot2::.pt) + article_text() +
    labs(x='Log2 fold enrichment\n(motif)', y='Log2 fold change\n(expression in ctx neurons)')

ggsave('plots/paper/fig_regs/primed_motif_enrichment_fc_scatter_text.pdf', width=10, height=8, units='cm')


ggplot(plot_df, aes(logodds, avg_log2FC, label=tf, alpha=padj<0.05)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    theme_rangeframe() + scale_axis_rangeframe() +
    geom_point(size=0.5) + article_text() +
    labs(x='Log2 fold enrichment\n(motif)', y='Log2 fold change\n(expression in ctx neurons)')

ggsave('plots/paper/fig_regs/primed_motif_enrichment_fc_scatter.pdf', width=10, height=8, units='cm')



#### Plot primed regions in UMAP ####
primed2peaks <- find_peaks_near_genes(all_bg, genes_primed, distance=4000, only_tss=F)
primed_peaks <- rownames(primed2peaks)[rowMaxs(primed2peaks)>0] 

all_peaks_srt$is_primed <- all_peaks_srt$repr_region%in%primed_peaks
all_peaks_srt %>% dim_plot(group.by='is_primed', pt.size=0.5, order=T)





#### Bivalent / switching ####
rna_lin_de <- read_tsv('data_/results/diff_expression/RNA_DE_lineage.tsv')
bival_peaks <- read_tsv('data_/CT/bivalent_peaks.tsv')
switch_peaks <- read_tsv('data_/CT/switching_peaks.tsv')

bival_motif_enrich <- map_dfr(set_names(unique(bival_peaks$max_group)), function(x){
    print(x)
    mpeaks <- bival_peaks %>% filter(max_group==x, mark=='H3K4me3') %>% pull(H3K4me3) %>% unique()
    puse <- peak_meta_annot %>% filter(H3K4me3%in%mpeaks) %>% pull(repr_region) %>% unique()
    return(tf_motif_enrichment(puse, all_motif_mat, GRangesToString(all_bg), parallel=T))
}, .id='group')

switch_motif_enrich <- map_dfr(set_names(unique(switch_peaks$max_group)), function(x){
    print(x)
    mpeaks <- switch_peaks %>% filter(max_group==x, mark=='H3K27ac') %>% pull(H3K27ac) %>% unique()
    puse <- peak_meta_annot %>% filter(H3K27ac%in%mpeaks) %>% pull(repr_region)
    return(tf_motif_enrichment(puse, all_motif_mat, GRangesToString(all_bg), parallel=T))
}, .id='group')

bival_df <- bival_motif_enrich %>% 
    inner_join(motif2tf) %>%
    inner_join(rna_lin_de, by=c('group'='cluster', 'tf'='gene')) %>% 
    filter(origin=='JASPAR2020')

bival_df %>% write_tsv('data_/all_peaks/bival_motif_enrich.tsv')

switch_df <- switch_motif_enrich %>% 
    inner_join(motif2tf) %>%
    inner_join(rna_lin_de, by=c('group'='cluster', 'tf'='gene')) %>% 
    filter(origin=='JASPAR2020')

switch_df %>% write_tsv('data_/all_peaks/switch_motif_enrich.tsv')


p1 <- ggplot(bival_df, aes(logodds, avg_log2FC, label=tf, alpha=padj<0.05)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_text(size=4/ggplot2::.pt) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_grid(group~.) + article_text() +
    labs(x='Log2 fold enrichment\n(motif)', y='Log2 fold change\n(expression in branch)', title='Bivalent')

p2 <- ggplot(switch_df, aes(logodds, avg_log2FC, label=tf, alpha=padj<0.05)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_text(size=4/ggplot2::.pt) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_grid(group~.) + article_text() +
    labs(x='Log2 fold enrichment\n(motif)', y='Log2 fold change\n(expression in branch)', title='Switching')

p1 | p2
ggsave('plots/paper/fig_regs/bival_switch_motif_enrichment_fc_scatter_text.pdf', width=15, height=10, units='cm')



p1 <- ggplot(bival_df, aes(logodds, avg_log2FC, label=tf, alpha=padj<0.05)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_point(size=0.5) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_grid(group~.) +
    labs(x='Log2 fold enrichment\n(motif)', y='Log2 fold change\n(expression in branch)', title='Bivalent')

p2 <- ggplot(switch_df, aes(logodds, avg_log2FC, label=tf, alpha=padj<0.05)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_point(size=0.5) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_grid(group~.) +
    labs(x='Log2 fold enrichment\n(motif)', y='Log2 fold change\n(expression in branch)', title='Switching')

p1 | p2
ggsave('plots/paper/fig_regs/bival_switch_motif_enrichment_fc_scatter.pdf', width=25, height=20, units='cm')




#### Top lineage peaks ####
all_top_peaks <- read_tsv('data_/results/diff_expression/all_marks_top_DA_lineage_coarse.tsv')

H3K27ac_motif_mat <- marks$H3K27ac@assays$peaks@motifs@data
H3K27ac_background <- colnames(H3K27ac_cluster_peaks)

H3K27me3_motif_mat <- marks$H3K27me3@assays$peaks@motifs@data
H3K27me3_background <- colnames(H3K27me3_cluster_peaks)

H3K4me3_motif_mat <- marks$H3K4me3@assays$peaks@motifs@data
H3K4me3_background <- colnames(H3K4me3_cluster_peaks)

bg_list <- list(
    H3K27ac = H3K27ac_background,
    H3K27me3 = H3K27me3_background,
    H3K4me3 = H3K4me3_background
)

lineage_motif_enrich <- map_dfr(set_names(unique(all_top_peaks$group)), function(x){
    map_dfr(set_names(unique(all_top_peaks$mark)), function(y){
        puse <- all_top_peaks %>% filter(mark==y, group==x) %>% pull(feature) %>% unique()
        return(FindMotifs(marks[[y]], puse, background=bg_list[[y]]))
    }, .id='mark')
}, .id='group')

lineage_motif_df <- lineage_motif_enrich %>% 
    as_tibble() %>% 
    # group_by(group, mark) %>% 
    mutate(padj=p.adjust(pvalue, method='fdr')) %>% 
    inner_join(motif2tf) %>%
    inner_join(rna_lin_de, by=c('group'='cluster', 'tf'='gene')) %>% 
    filter(origin=='JASPAR2020')

lineage_motif_df %>% write_tsv('data_/all_peaks/mark_lineage_motif_enrich.tsv')

ggplot(lineage_motif_df, aes(log2(fold.enrichment), avg_log2FC, label=tf, alpha=padj<0.01)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_text(size=3/ggplot2::.pt) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_grid(group~mark) + article_text() +
    labs(x='Log2 fold enrichment\n(motif)', y='Log2 fold change\n(expression in branch)')

ggsave('plots/paper/fig_regs/mark_lineage_motif_enrichment_fc_scatter_text.pdf', width=25, height=20, units='cm')

ggplot(lineage_motif_df, aes(log2(fold.enrichment), avg_log2FC, label=tf, alpha=padj<0.01)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_point(size=0.5) +
    theme_rangeframe() + scale_axis_rangeframe() +
    facet_grid(group~mark) + article_text() +
    labs(x='Log2 fold enrichment\n(motif)', y='Log2 fold change\n(expression in branch)')

ggsave('plots/paper/fig_regs/mark_lineage_motif_enrichment_fc_scatter.pdf', width=25, height=20, units='cm')









