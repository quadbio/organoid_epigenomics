library(tidyverse)
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/wrapper.R')
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
marks <- read_rds('data/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')

drugs_d158_H3K27me3 <- read_rds('data/drugs/drugs_d15_d18_A395_v1_v1.2lines_srt.rds')
drugs_d158_H3K27ac <- read_rds('data/drugs/drugs_d15_d18_A485_v1_v1.2lines_srt.rds')

#### Feature plots ####
goi <- c('DLX3', 'TFAP2A', 'NEUROD1', 'DLX6', 'LHX5', 'RSPO3', 'WLS', 'WNT3A', 'FOXG1')
feature_plot(rna, features=goi, reduction='cssumap', order=T)
feature_plot(drugs_d158_H3K27me3, features=goi, order=T, reduction='humap')

dim_plot(drugs_d158_H3K27me3, reduction='humap', label=T)


#### Get detected peaks ####
marks$H3K27ac[['peaks_bin']] <- CreateAssayObject((marks$H3K27ac[['peaks']]@data > 0)*1)
marks$H3K27ac <- Pando::aggregate_assay(marks$H3K27ac, assay='peaks_bin', group_name='clusters')
H3K27ac_peak_detect <- marks$H3K27ac@assays$peaks_bin@misc$summary$clusters
H3K27ac_detected <- colnames(H3K27ac_peak_detect)[colMaxs(H3K27ac_peak_detect)>0.05]

marks$H3K27me3[['peaks_bin']] <- CreateAssayObject((marks$H3K27me3[['peaks']]@data > 0)*1)
marks$H3K27me3 <- Pando::aggregate_assay(marks$H3K27me3, assay='peaks_bin', group_name='clusters')
H3K27me3_peak_detect <- marks$H3K27me3@assays$peaks_bin@misc$summary$clusters
H3K27me3_detected <- colnames(H3K27me3_peak_detect)[colMaxs(H3K27me3_peak_detect)>0.05]

marks$H3K4me3[['peaks_bin']] <- CreateAssayObject((marks$H3K4me3[['peaks']]@data > 0)*1)
marks$H3K4me3 <- Pando::aggregate_assay(marks$H3K4me3, assay='peaks_bin', group_name='clusters')
H3K4me3_peak_detect <- marks$H3K4me3@assays$peaks_bin@misc$summary$clusters
H3K4me3_detected <- colnames(H3K4me3_peak_detect)[colMaxs(H3K4me3_peak_detect)>0.05]


#### Get peaks detected early ####
H3K27ac_early_clusters <- marks$H3K27ac@meta.data %>% 
    filter(celltype_jf%in%c('non_nect', 'psc', 'nect') & age%in%c('EB', '15d', '35d')) %>% 
    pull(clusters) %>% unique()
H3K27ac_early_detect <- marks$H3K27ac@assays$peaks_bin@misc$summary$clusters[H3K27ac_early_clusters, ]
H3K27ac_detected_early <- colnames(H3K27ac_early_detect)[colMaxs(H3K27ac_early_detect)>0.01]

H3K27me3_early_clusters <- marks$H3K27me3@meta.data %>% 
    filter(celltype_jf%in%c('non_nect', 'psc', 'nect') & age%in%c('EB', '15d', '35d')) %>% 
    pull(clusters) %>% unique()
H3K27me3_early_detect <- marks$H3K27me3@assays$peaks_bin@misc$summary$clusters[H3K27me3_early_clusters, ]
H3K27me3_detected_early <- colnames(H3K27me3_early_detect)[colMaxs(H3K27me3_early_detect)>0.01]

H3K4me3_early_clusters <- marks$H3K4me3@meta.data %>% 
    filter(celltype_jf%in%c('non_nect', 'psc', 'nect') & age%in%c('EB', '15d', '35d')) %>% 
    pull(clusters) %>% unique()
H3K4me3_early_detect <- marks$H3K4me3@assays$peaks_bin@misc$summary$clusters[H3K4me3_early_clusters, ]
H3K4me3_detected_early <- colnames(H3K4me3_early_detect)[colMaxs(H3K4me3_early_detect)>0.01]



#### Get bigwigs and summarize over peaks ####
bw_path <- '/local2/USERS/jfleck/projects/cutntag/spike_norm/bamcoverage'
bw_names <- list.files(bw_path)
bw_rpkm_files <- list.files(bw_path, recursive=T, pattern='*_rpkm.bw', full.names=T)

names(bw_rpkm_files) <- bw_names

bw_ranges_list <- map(bw_rpkm_files, rtracklayer::import)

bw_H3K27ac_A485_list <- bw_ranges_list[str_detect(names(bw_ranges_list), 'H3K27ac.*A485.*(18|15)')]
bw_H3K27me3_A395_list <- bw_ranges_list[str_detect(names(bw_ranges_list), 'H3K27me3.*A395.*(18|15)')]

bw_H3K27ac_DMSO_list <- bw_ranges_list[str_detect(names(bw_ranges_list), 'H3K27ac.*DMSO.*(18|15)')]
bw_H3K27me3_DMSO_list <- bw_ranges_list[str_detect(names(bw_ranges_list), 'H3K27me3.*DMSO.*(18|15)')]


H3K27me3_A395_peak_score_mat <- bw_H3K27me3_A395_list %>% map(function(x){
    olaps <- findOverlaps(StringToGRanges(H3K27me3_detected), x)
    score_mat <- as.matrix(x$score[subjectHits(olaps)])
    peaks_scores <- Pando::aggregate_matrix(score_mat, groups = queryHits(olaps)) %>% 
        as.numeric()
    names(peaks_scores) <- H3K27me3_detected
    return(peaks_scores)
}) %>% purrr::reduce(rbind) %>% Matrix(sparse=T)

rownames(H3K27me3_A395_peak_score_mat) <- names(bw_H3K27me3_A395_list)


H3K27ac_A485_peak_score_mat <- bw_H3K27ac_A485_list %>% map(function(x){
    olaps <- findOverlaps(StringToGRanges(H3K27ac_detected), x)
    score_mat <- as.matrix(x$score[subjectHits(olaps)])
    peaks_scores <- Pando::aggregate_matrix(score_mat, groups = queryHits(olaps)) %>% 
        as.numeric()
    names(peaks_scores) <- H3K27ac_detected
    return(peaks_scores)
}) %>% purrr::reduce(rbind) %>% Matrix(sparse=T)

rownames(H3K27ac_A485_peak_score_mat) <- names(bw_H3K27ac_A485_list)


H3K27me3_DMSO_peak_score_mat <- bw_H3K27me3_DMSO_list %>% map(function(x){
    olaps <- findOverlaps(StringToGRanges(H3K27me3_detected), x)
    score_mat <- as.matrix(x$score[subjectHits(olaps)])
    peaks_scores <- Pando::aggregate_matrix(score_mat, groups = queryHits(olaps)) %>% 
        as.numeric()
    names(peaks_scores) <- H3K27me3_detected
    return(peaks_scores)
}) %>% purrr::reduce(rbind) %>% Matrix(sparse=T)

rownames(H3K27me3_DMSO_peak_score_mat) <- names(bw_H3K27me3_DMSO_list)


H3K27ac_DMSO_peak_score_mat <- bw_H3K27ac_DMSO_list %>% map(function(x){
    olaps <- findOverlaps(StringToGRanges(H3K27ac_detected), x)
    score_mat <- as.matrix(x$score[subjectHits(olaps)])
    peaks_scores <- Pando::aggregate_matrix(score_mat, groups = queryHits(olaps)) %>% 
        as.numeric()
    names(peaks_scores) <- H3K27ac_detected
    return(peaks_scores)
}) %>% purrr::reduce(rbind) %>% Matrix(sparse=T)

rownames(H3K27ac_DMSO_peak_score_mat) <- names(bw_H3K27ac_DMSO_list)


#### Calculate mean logFC to control ####
#### A395 H3K27me3 ####

H3K27me3_A395_15d_mat <- H3K27me3_A395_peak_score_mat[str_detect(rownames(H3K27me3_A395_peak_score_mat), '15d'), ]
H3K27me3_A395_18d_mat <- H3K27me3_A395_peak_score_mat[str_detect(rownames(H3K27me3_A395_peak_score_mat), '18d$'), ]
H3K27me3_A395_18d2_mat <- H3K27me3_A395_peak_score_mat[str_detect(rownames(H3K27me3_A395_peak_score_mat), '18d_rep2'), ]

H3K27me3_DMSO_15d_mat <- H3K27me3_DMSO_peak_score_mat[str_detect(rownames(H3K27me3_DMSO_peak_score_mat), '15d'), ]
H3K27me3_DMSO_18d_mat <- H3K27me3_DMSO_peak_score_mat[str_detect(rownames(H3K27me3_DMSO_peak_score_mat), '18d$'), ]
H3K27me3_DMSO_18d2_mat <- H3K27me3_DMSO_peak_score_mat[str_detect(rownames(H3K27me3_DMSO_peak_score_mat), '18d_rep2'), ]

H3K27me3_A395_15d_fc <- (t(log2(H3K27me3_A395_15d_mat+1)) - log2(H3K27me3_DMSO_15d_mat+1)) %>% 
    as_tibble(rownames='peak') %>% pivot_longer(!peak, names_to='sample', values_to='fc')
H3K27me3_A395_18d_fc <- (t(log2(H3K27me3_A395_18d_mat+1)) - log2(H3K27me3_DMSO_18d_mat+1)) %>% 
    as_tibble(rownames='peak') %>% pivot_longer(!peak, names_to='sample', values_to='fc')
H3K27me3_A395_18d2_fc <- (t(log2(H3K27me3_A395_18d2_mat+1)) - log2(H3K27me3_DMSO_18d2_mat+1)) %>% 
    as_tibble(rownames='peak') %>% pivot_longer(!peak, names_to='sample', values_to='fc')

H3K27me3_peak_genes <- Signac::ClosestFeature(marks$H3K27ac, StringToGRanges(H3K27me3_detected)) %>% 
    as_tibble()

H3K27me3_A395_fc <- bind_rows(H3K27me3_A395_18d2_fc, H3K27me3_A395_18d_fc, H3K27me3_A395_15d_fc) %>% 
    filter(peak%in%H3K27me3_detected_early) %>% 
    mutate(
        age=str_replace(sample, '.+A395_(\\d+d).*', '\\1'),
        conc=factor(as.numeric(str_replace(sample, '.+HW_(\\d+)-A395_\\d+d.*', '\\1'))),
        replicate=as.numeric(ifelse(str_detect(sample, 'rep2$'), '2', '1'))
    ) %>% inner_join(H3K27me3_peak_genes, by=c('peak'='query_region'))
    
    
H3K27me3_A395_mean_fc <- H3K27me3_A395_fc %>% 
    group_by(peak) %>% summarize(fc=mean(fc)) %>% arrange(desc(fc)) %>% mutate(peak=factor(peak, levels=unique(.$peak))) %>% 
    inner_join(H3K27me3_peak_genes, by=c('peak'='query_region'))

H3K27me3_A395_fc$peak <- factor(H3K27me3_A395_fc$peak, levels=H3K27me3_A395_mean_fc$peak)

H3K27me3_A395_fc %>% write_tsv('data/drugs/drugs_CT_A395_H3K27me3_peaks_logfc.tsv')



#### Join with DE results ####
H3K27me3_A395_fc <- read_tsv('data/drugs/drugs_CT_A395_H3K27me3_peaks_logfc.tsv')

A395_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_global_inhib_de.tsv')
A395_global_inhib_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_concs_global_inhib_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))
A395_nepi_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_nepi_inhib_de.tsv')
A395_cluster_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_cluster_de.tsv')

top_global_de <- A395_global_de %>% filter(p_val_adj<1e-2, avg_log2FC>0.1) %>% pull(feature) %>% unique()
top_nepi_de <- A395_nepi_de %>% filter(p_val_adj<1e-2, avg_log2FC>0.1) %>% pull(feature) %>% unique()

plot_df <- H3K27me3_A395_fc
p_s <- ggplot(plot_df, aes(peak, fc, color=conc)) +
    geom_point(size=0.1, alpha=0.5) +
    geom_hline(yintercept = 0, size=0.3) +
    geom_line(data=H3K27me3_A395_mean_fc, mapping=aes(x=as.numeric(peak)), color='#34495e', size=1) +
    scale_color_manual(values=conc_colors) +
    scale_x_discrete(expand=c(0.02,0)) +
    scale_y_continuous(breaks=seq(-4,4,2), limits=c(-4,4)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    no_x_text() +
    no_legend() +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    labs(y='log2 fold change', x='Peak') 

plot_df <- H3K27me3_A395_mean_fc
p_hist <- ggplot(plot_df, aes(fc)) +
    geom_histogram(color='black', fill='grey') +
    scale_x_continuous(breaks=seq(-4,4,2), limits=c(-4,4)) +
    theme_void() +
    coord_flip()

(p_s | p_hist) + plot_layout(widths=c(5,1))


conc_colors <- c('1'='#A2D9DE', '3'='#3A90B0', '10'='#0B4771')

plot_df <- inner_join(A395_global_inhib_de, H3K27me3_A395_mean_fc, by=c('feature'='gene_name')) %>% 
    mutate(
        conc=str_replace(inhib_annotation, '.+-(\\d+)', '\\1'),
        age=str_replace(inhib_annotation, '.+_(d\\d+)-\\d+', '\\1')
    ) %>% 
    group_by(conc, peak) %>% 
    summarize(fc=mean(avg_log2FC)) %>% 
    mutate(
        fc_clip=clip_abs(fc, 1),
        peak=factor(peak, levels=H3K27me3_A395_mean_fc$peak),
        conc=factor(conc, levels=rev(c(1,3,10)))
    ) 

ggplot(plot_df, aes(peak, conc, fill=fc_clip)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(bigrad(pals::brewer.rdbu, bias=1.5)[20:180])) + 
    scale_x_discrete(expand=c(0.02,0)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    no_x_text() +
    theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()
    ) 



#### UMAPs ####

dim_plot(drugs_d158_H3K27me3, reduction='humap', group.by='inhib_annotation') +
    scale_color_manual(values=inhibitors_colors_15)

ggsave('plots/paper/fig4/fig4_drugs_d158_H3K27me3_samples_umap.png', width=6,height=5)


cluster_colors <- c('#6c2d76', '#ab4f93', '#b376a3', '#715aa3','#398fa5' , '#bc92be', '#f78f1e' ,'#044972')
dim_plot(drugs_d158_H3K27me3, reduction='humap', group.by='RNA_snn_res.0.2') +
    scale_color_manual(values=cluster_colors)

ggsave('plots/paper/fig4/fig4_drugs_d158_H3K27me3_samples_umap.png', width=6,height=5)




#### Get DEG sets and compute module score ####

cluster_srt <- read_rds('data/all_RNA_marks_combined_clusters_srt.rds')
cluster_meta <- read_tsv('data/all_RNA_cluster_meta.tsv')

A395_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_global_inhib_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))
A395_global_inhib_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_concs_global_inhib_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))
A395_nepi_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_nepi_inhib_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))
A395_cluster_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A395_cluster_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))

nepi_deg <- A395_nepi_de %>% filter(padj<1e-4, avg_log2FC>0.2) %>% pull(feature) %>% unique()
global_deg <- A395_global_de %>% filter(padj<1e-4, avg_log2FC>0.2) %>% pull(feature) %>% unique()

cluster_srt@active.assay <- 'RNA'
cluster_srt <- AddModuleScore(cluster_srt, features=list('nepi_DEG'=nepi_deg, 'global_DEG'=global_deg))
cluster_srt$RNA_nepi_DEG_score <- cluster_srt$Cluster1
cluster_srt$RNA_global_DEG_score <- cluster_srt$Cluster2

cluster_srt@active.assay <- 'H3K27me3_RNA'
cluster_srt <- AddModuleScore(cluster_srt, features=list('nepi_DEG'=nepi_deg, 'global_DEG'=global_deg))
cluster_srt$H3K27me3_nepi_DEG_score <- cluster_srt$Cluster1
cluster_srt$H3K27me3_global_DEG_score <- cluster_srt$Cluster2

cluster_srt@active.assay <- 'H3K27ac_RNA'
cluster_srt <- AddModuleScore(cluster_srt, features=list('nepi_DEG'=nepi_deg, 'global_DEG'=global_deg))
cluster_srt$H3K27ac_nepi_DEG_score <- cluster_srt$Cluster1
cluster_srt$H3K27ac_global_DEG_score <- cluster_srt$Cluster2

cluster_srt@active.assay <- 'H3K4me3_RNA'
cluster_srt <- AddModuleScore(cluster_srt, features=list('nepi_DEG'=nepi_deg, 'global_DEG'=global_deg))
cluster_srt$H3K4me3_nepi_DEG_score <- cluster_srt$Cluster1
cluster_srt$H3K4me3_global_DEG_score <- cluster_srt$Cluster2

cluster_srt$bivalent_nepi_DEG_score <- scale(cluster_srt$H3K4me3_nepi_DEG_score) * scale(cluster_srt$H3K27me3_nepi_DEG_score) * sign(scale(cluster_srt$H3K4me3_nepi_DEG_score) - scale(cluster_srt$H3K27me3_nepi_DEG_score))
cluster_srt$bivalent_global_DEG_score <- scale(cluster_srt$H3K4me3_global_DEG_score) * scale(cluster_srt$H3K27me3_global_DEG_score) * sign(scale(cluster_srt$H3K4me3_global_DEG_score) - scale(cluster_srt$H3K27me3_global_DEG_score))

cluster_srt$switch_nepi_DEG_score <- scale(cluster_srt$H3K27ac_nepi_DEG_score) * scale(cluster_srt$H3K27me3_nepi_DEG_score) * sign(scale(cluster_srt$H3K27ac_nepi_DEG_score) - scale(cluster_srt$H3K27me3_nepi_DEG_score))
cluster_srt$switch_global_DEG_score <- scale(cluster_srt$H3K27ac_global_DEG_score) * scale(cluster_srt$H3K27me3_global_DEG_score) * sign(scale(cluster_srt$H3K27ac_global_DEG_score) - scale(cluster_srt$H3K27me3_global_DEG_score))

cluster_meta <- cluster_srt@meta.data %>% 
    as_tibble(rownames='cluster') %>% 
    inner_join(as_tibble(cluster_srt[['tree']]@cell.embeddings, rownames='cluster'))

cluster_graph_meta <- cluster_meta %>% 
    mutate(
        H3K4me3_scale = scale(H3K4me3_global_DEG_score),
        H3K27me3_scale = scale(H3K27me3_global_DEG_score),
        H3K4me3_score = H3K4me3_scale - min(H3K4me3_scale),
        H3K27me3_score = H3K27me3_scale - min(H3K27me3_scale),
        bivalent_score = clip_abs(H3K4me3_score * H3K27me3_score, 4),
        bivalent_color = clip_abs(scale(H3K4me3_score) - scale(H3K27me3_score), 2)
    )

bival_cols <- colorRampPalette(c('#F773A5', '#810F7C', '#75B3D8'))(100)
print_scale(bival_cols)


ggplot(cluster_graph_meta, aes(x=tree_1, y=-tree_2)) +
    geom_point(aes(fill=H3K27me3_global_DEG_score, size=H3K27me3_global_DEG_score), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.blues)) +
    scale_y_reverse() +
    theme_void()

p1 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=bivalent_color, size=abs(bivalent_score)), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=rev(bival_cols)) +
    scale_y_reverse() +
    theme_void()

p2 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=H3K27me3_global_DEG_score, size=H3K27me3_global_DEG_score), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.blues)) +
    scale_y_reverse() +
    theme_void()

p3 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=H3K4me3_global_DEG_score, size=H3K4me3_global_DEG_score), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.rdpu)) +
    scale_y_reverse() +
    theme_void()

p4 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=RNA_global_DEG_score, size=RNA_global_DEG_score), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.orrd)) +
    scale_y_reverse() +
    theme_void()

p4 / p1 / p2 / p3
ggsave('plots/paper/fig4/rev_fig4_global_deg_divalent_scores_umap.pdf', width=6, height=7)
    


p1 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=RNA_global_DEG_score, size=RNA_global_DEG_score), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.orrd)) +
    scale_y_reverse() +
    theme_void()

p2 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=H3K27me3_global_DEG_score, size=H3K27me3_global_DEG_score), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.blues)) +
    scale_y_reverse() +
    theme_void()

p3 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=H3K4me3_global_DEG_score, size=H3K4me3_global_DEG_score), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.rdpu)) +
    scale_y_reverse() +
    theme_void()

p4 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=H3K27ac_global_DEG_score, size=H3K27ac_global_DEG_score), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.ylgn)) +
    scale_y_reverse() +
    theme_void()

p1 / p2 / p3 / p4 + plot_layout(guides='collect')
ggsave('plots/paper/fig4/rev_fig4_global_deg_scores_umap.pdf', width=6, height=7)




p1 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=bivalent_color, size=abs(bivalent_global_DEG_score)), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4), limits = c(0,2)) +
    scale_fill_gradientn(colors=rev(bival_cols)) +
    scale_y_reverse() +
    theme_void()

p2 <- ggplot(cluster_graph_meta, aes(x=CSSUMAP_1, y=CSSUMAP_2)) +
    geom_point(aes(fill=switch_global_DEG_score, size=abs(switch_global_DEG_score)), shape=21, stroke=0.2) +
    scale_size_continuous(range=c(0.2,4)) +
    scale_fill_gradientn(colors=grad(pals::brewer.bugn)) +
    scale_y_reverse() +
    theme_void()

p1 / p2
ggsave('plots/paper/fig4/rev_fig4_global_deg_comb_scores_umap.pdf', width=6, height=3.5)



#### Barplots with scores ####

plot_df <- cluster_graph_meta %>% 
    mutate(
        is_neuro = celltype_jf%in%c('non_nect', 'astrocytes', 'OPC', 'choroid_plexus'),
        celltype_coarse = case_when(
            celltype_jf%in%c('rhom_ex', 'mesen_ex', 'dien_ex', 'ctx_ex') ~ 'neuron',
            celltype_jf%in%c('nt_npc', 'dien_npc', 'ctx_npc', 'ctx_ip') ~ 'npc',
            T ~ celltype_jf
        ),
        celltype_fac = factor(celltype_coarse,
                              levels=c(
                                  'psc', 'nect', 'npc', 'neuron',
                                  'non_nect', 'choroid_plexus', 'astrocytes', 'OPC'
                                  )),
        # celltype_fac = factor(celltype_jf, 
        #                       levels=c(
        #                           'psc', 'nect', 'ctx_npc','ctx_ip', 'ctx_ex',
        #                           'non_nect', 'choroid_plexus', 'astrocytes', 'OPC'
        #                           ))
    ) %>% 
    filter(!is.na(celltype_fac))
    


ggplot(plot_df, aes(cluster, H3K27me3_global_DEG_score, fill=celltype_jf)) +
    labs(y='Activity score', x='', title='H3K27me3') +
    no_x_text() +
    geom_bar(stat='identity')


p1 <- ggplot(plot_df, aes(celltype_fac, H3K27me3_global_DEG_score, fill=celltype_fac)) +
    labs(y='Activity score', x='', title='H3K27me3') +
    no_x_text()
 
p2 <- ggplot(plot_df, aes(celltype_fac, H3K4me3_global_DEG_score, fill=celltype_fac)) +
    labs(y='Activity score', x='', title='H3K4me3') +
    no_x_text()
 
p3 <- ggplot(plot_df, aes(celltype_fac, RNA_global_DEG_score, fill=celltype_fac)) +
    labs(y='Activity score', x='', title='RNA') +
    no_x_text()
 
p4 <- ggplot(plot_df, aes(celltype_fac, bivalent_color, fill=celltype_fac)) +
    labs(y='Activity score', x='Celltype', title='Bivalency') +
    rotate_x_text(40)
 
(p1 / p2 / p3 / p4) & 
    no_legend() &    
    geom_bar(stat='summary') &
    geom_errorbar(stat='summary', width=0.2, col='darkgrey') &
    geom_beeswarm(size=0.1, shape=16, corral.width = 1.2) &
    scale_fill_manual(values=pantone_celltype) &
    facet_grid(~is_neuro, space='free', scales='free') &
    rotate_x_text(40) &
    article_text() &
    theme(
        strip.text = element_blank()
    ) 
ggsave('plots/paper/fig4/fig4_global_deg_scores_coarse_types_barplot.pdf', width=6, height=15, units='cm')
ggsave('plots/paper/fig4/fig4_global_deg_scores_coarse_types_barplot.png', width=6, height=15, units='cm')
   

ggplot(plot_df, aes(celltype_fac, H3K27me3_global_DEG_score, fill=celltype_jf)) +
    geom_boxplot() +
    scale_fill_manual(values=pantone_celltype)
    

ggplot(plot_df, aes(celltype_fac, H3K27me3_global_DEG_score, fill=celltype_jf)) +
    geom_bar(stat='summary') +
    geom_errorbar(stat='summary', width=0.2) +
    scale_fill_manual(values=pantone_celltype)
    





#### Compare nepi DE with cluster DE ####
nepi_deg <- A395_nepi_de %>% filter(padj<1e-4, avg_log2FC>0.2) %>% pull(feature) %>% unique()
global_deg <- A395_global_de %>% filter(padj<1e-4, avg_log2FC>0.2) %>% pull(feature) %>% unique()
cluster_deg <- A395_cluster_de %>% filter(padj<1e-4, avg_log2FC>0.2) %>% pull(feature) %>% unique()
ectopic_deg <- A395_cluster_de %>% mutate(cluster=cluster+1) %>% filter(padj<1e-4, avg_log2FC>0.2, cluster%in%c(7,5,8,4)) %>% pull(feature) %>% unique()

cluster_colors <- c('#6c2d76', '#ab4f93', '#b376a3', '#715aa3','#398fa5' , '#bc92be', '#f78f1e' ,'#044972')
names(cluster_colors) <- 1:8

ectopic_de <- A395_cluster_de %>% 
    mutate(cluster=cluster+1) %>% 
    filter(padj<1e-4, avg_log2FC>1, cluster%in%c(7,5,8,4)) %>% 
    group_by(feature) %>% 
    filter(avg_log2FC==max(avg_log2FC)) %>% 
    distinct(feature, cluster, avg_log2FC)


plot_df <- A395_nepi_de %>% left_join(ectopic_de, by='feature') %>% 
    mutate(avg_log2FC.y = replace_na(avg_log2FC.y, 0))

ggplot(arrange(plot_df, !is.na(cluster)), aes(avg_log2FC.x, -log10(p_val), color=factor(cluster), label=feature, size=avg_log2FC.y)) +
    geom_text_repel(data=filter(plot_df, !is.na(cluster) & padj<1e-4 & avg_log2FC.x>0.2), size=5/ggplot2::.pt, max.overlaps=9999) +
    geom_point() +
    scale_size_continuous(range=c(0.1, 0.8)) +
    scale_color_manual(values=cluster_colors) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    labs(x='log2(fold change)', y='-log10(p-value)', title='DE in Neuroepithelium', color='Cluster', size='log2FC (cluster DE)')

ggsave('plots/paper/fig4/fig4_nepi_deg_vulcano.pdf', width=6, height=4, units='cm')



drugs_d158_H3K27me3$concentration <- drugs_d158_H3K27me3$inhib_annotation %>% str_replace('.+-(\\d+)', '\\1') 
drugs_d158_H3K27me3$concentration <- ifelse(str_detect(drugs_d158_H3K27me3$concentration, 'DMSO'), 'DMSO', drugs_d158_H3K27me3$concentration)
drugs_d158_H3K27me3$concentration <- factor(drugs_d158_H3K27me3$concentration, levels=c('DMSO', '1', '3', '10'))
feature_plot(drugs_d158_H3K27me3, features=c('LMO4', 'MEIS2', 'NEFM', 'MAP2', 'ANXA2', 'S100A10', 'AHNAK', 'STMN2', 'POU5F1'), reduction='humap', order=T, split.by='concentration', pt.size=0.1) &
    scale_color_gradientn(colors=gyorrd())

ggsave('plots/paper/fig4/fig4_nepi_deg_split_feature_umap.png', width=10, height=14, bg='white')

feature_plot(drugs_d158_H3K27me3, features=c('LMO4', 'MEIS2', 'NEFM', 'MAP2', 'ANXA2', 'S100A10', 'AHNAK', 'STMN2', 'POU5F1'), reduction='humap', order=T, pt.size=0.1) &
    scale_color_gradientn(colors=gyorrd())

ggsave('plots/paper/fig4/fig4_nepi_deg_feature_umap.png', width=8, height=8, bg='white')





#### TF motif enrichment for peaks 
library(doParallel)
registerDoParallel(36)
data(motif2tf)

H3K27me3_A395_fc <- read_tsv('data/drugs/drugs_CT_A395_H3K27me3_peaks_logfc.tsv')

# genes_use <- filter(plot_df, cluster==7 & padj<1e-4 & avg_log2FC.x>0.2)$feature
genes_use <- filter(A395_cluster_de, cluster==6 & padj<1e-4 & avg_log2FC>1)$feature

peaks_in <- H3K27me3_A395_fc %>% 
    filter(gene_name%in%genes_use, fc < -1, distance < 10000) %>% pull(peak) %>% unique()

peaks_in %>% write('data/drugs/results/drugs_cl7_regions_motif_enrichment.txt')

peaks_all <- H3K27me3_A395_fc$peak %>% unique()

group_vec <- peaks_all %in% peaks_in

motif_mat <- marks$H3K27me3@assays$peaks@motifs@data
motif_ftest <- map_par(colnames(motif_mat), function(mot){
    ftest <- table(motif_mat[peaks_all,mot], group_vec) %>% fisher.test()
    out_tbl <- tibble(
        motif = mot,
        pval = ftest$p.value,
        odds_ratio = ftest$estimate
    ) %>% return()
}, parallel=T) %>% 
    bind_rows() %>% 
    mutate(padj=p.adjust(pval, 'fdr')) %>% 
    mutate(logodds=log2(odds_ratio)) 

motif_sig <- motif_ftest %>% 
    filter(padj<0.05) %>% 
    inner_join(motif2tf) %>% 
    filter(origin=='JASPAR2020') %>% 
    arrange(desc(logodds))

feature_plot(drugs_d158_H3K27me3, features=unique(motif_sig$tf), reduction='humap', order=T, pt.size=0.1) &
    scale_color_gradientn(colors=gyorrd())

feature_plot(drugs_d158_H3K27me3, features=genes_use[1:10], reduction='humap', order=T, pt.size=0.1) &
    scale_color_gradientn(colors=gyorrd())



plot_df <- motif_sig %>% 
    inner_join(filter(A395_cluster_de, cluster==6), by=c('tf'='feature')) %>% 
    group_by(tf) %>% summarize(avg_log2FC=mean(avg_log2FC), logodds=mean(logodds)) %>% 
    arrange(logodds) %>% mutate(tf=factor(tf, levels=unique(.$tf))) 

ggplot(plot_df, aes(avg_log2FC, logodds, label=tf)) +
    geom_text()

ggplot(plot_df, aes(logodds, tf, label=tf, fill=avg_log2FC)) +
    geom_bar(stat='identity') +
    scale_fill_gradientn(colors=rev(bigrad(pals::brewer.rdbu, 1.5)), limits=c(-2,2)) +
    article_text() +
    theme_rangeframe() + scale_axis_rangeframe() +
    labs(x='Motif log fold einrichment', y='Transcription factor', fill='log2 fold change\nin cluster 7')

ggsave('plots/paper/fig4/fig4_motif_enrichment_bar.pdf', width=6, height=8, units='cm')


feature_plot(drugs_d158_H3K27me3, features=as.character(plot_df$tf), reduction='humap', order=T, pt.size=0.3) &
    scale_color_gradientn(colors=gyorrd())
ggsave('plots/paper/fig4/fig4_motif_enrichment_tf_feature_umap.png', width=20, height=35)



#### Box- and bar plots to compare conditions ####
conc_colors <- c('1'='#A2D9DE', '3'='#3A90B0', '10'='#0B4771', 'DMSO'='#f78f1e')

genes_plot <- c('AHNAK', 'STMN2', 'ANXA2', 'NEFM', 'S100A10', 'POU5F1')
gene_expr <- drugs_d158_H3K27me3@assays$RNA@data[genes_plot, ] %>% 
    as_tibble(rownames='gene') %>% 
    pivot_longer(!gene, names_to='cell', values_to='expr') 

plot_df <- drugs_d158_H3K27me3@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(gene_expr) %>% 
    mutate(
        age=str_replace(inhib_annotation, '.+A395_(\\d+d).*', '\\1'),
        conc=str_replace(str_replace(inhib_annotation, '.+-(\\d+)', '\\1'), '.*DMSO.*', 'DMSO'),
    ) %>% 
    mutate(
        conc=factor(replace_na(conc, 'DMSO'), levels=c('DMSO', '1', '3', '10'))
    )

ggplot(plot_df, aes(conc, expr)) +
    geom_violin(size=0.2) +
    geom_boxplot(aes(fill=conc), outlier.size=0.05, outlier.shape=16, width=0.5, alpha=0.5, size=0.2) +
    facet_wrap(~gene) +
    scale_fill_manual(values=conc_colors) +
    article_text() +
    labs(x='Concentration', y='Expression', fill='Concentation')

ggsave('plots/paper/fig4/fig4_condition_expr_boxplot.pdf', width=7, height=6, units='cm')


ggplot(plot_df, aes(seurat_clusters, expr, fill=conc)) +
    # geom_violin(size=0.2, width=0.5) +
    geom_boxplot(outlier.size=0.05, outlier.shape=16, width=0.8, alpha=0.5, size=0.2) +
    facet_wrap(~gene) +
    scale_fill_manual(values=conc_colors) +
    article_text() +
    labs(x='Cluster', y='Expression', fill='Concentration')

ggsave('plots/paper/fig4/fig4_condition_expr_cluster_boxplot.pdf', width=20, height=10, units='cm')




#### Cluster proportions per conc and replicate ####

d158_enrich_sample <- map_dfr(set_names(c('day_15', 'day_18')), function(y){
    sample_idx <- (drugs_d158_H3K27me3$inhibitor_target %in% c('H3K27me3', 'DMSO')) & (drugs_d158_H3K27me3$orig.ident == y)
    perturb_vec <- factor(drugs_d158_H3K27me3$inhibitor_target[sample_idx], levels=c('H3K27me3', 'DMSO'))
    cluster_vec <- drugs_d158_H3K27me3$seurat_clusters[sample_idx]
    group_vec <- drugs_d158_H3K27me3$orig.ident[sample_idx]
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



ggplot(d158_enrich_sample, aes(group, -logodds, fill=sample)) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('#90c47d', '#3a90b0')) +
    article_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Cluster', y='Log fold enrichment', fill='Replicate')

ggsave('plots/paper/fig4/fig4_condition_replicate_enrich_bar.pdf', width=6, height=4, units='cm')




inhib_samples <- unique(drugs_d158_H3K27me3$inhib_annotation) %>% 
    {.[!str_detect(., 'DMSO')]} %>% set_names()

d158_enrich_sample <- map_dfr(inhib_samples, function(y){
    sample_idx <- (drugs_d158_H3K27me3$inhibitor_target == 'DMSO') | (drugs_d158_H3K27me3$inhib_annotation == y)
    perturb_vec <- factor(!str_detect(drugs_d158_H3K27me3$inhib_annotation[sample_idx], 'DMSO'))
    cluster_vec <- drugs_d158_H3K27me3$seurat_clusters[sample_idx]
    map_dfr(set_names(unique(cluster_vec)), function(c){
        ftest <- fisher.test(perturb_vec, cluster_vec==c)
        return(tibble(
            pval_fisher = ftest$p.value,
            odds_ratio = ftest$estimate,
            logodds = log2(ftest$estimate)
        ))
    }, .id='group')
}, .id='sample')


ggplot(d158_enrich_sample, aes(group, logodds)) +
    geom_boxplot(alpha=0.5, width=0.5, size=0.2, outlier.shape = NA) +
    geom_point(aes(fill=sample), shape=21, size=0.8, stroke=0.1) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values=inhibitors_colors_15) +
    article_text() +
    scale_axis_rangeframe() + theme_rangeframe() +
    labs(x='Cluster', y='Log fold enrichment', fill='Replicate')

ggsave('plots/paper/fig4/fig4_condition_replicate_enrich_boxplot.pdf', width=6, height=4, units='cm')



#### 

conc_colors <- c('1'='#A2D9DE', '3'='#3A90B0', '10'='#0B4771', 'DMSO'='#f78f1e')

library(entropy)

perturb_idx <- !str_detect(drugs_d158_H3K27me3$inhib_annotation, 'DMSO')
comp_ctrl <- drugs_d158_H3K27me3$seurat_clusters[!perturb_idx]
comp_perturb <- drugs_d158_H3K27me3$seurat_clusters[perturb_idx]


drugs_d158_H3K27me3 %>% dim_plot()

entropy(table(comp_ctrl))
entropy(table(comp_perturb))

inhib_samples <- unique(drugs_d158_H3K27me3$inhib_annotation) %>% set_names()


drugs_d158_H3K27me3 <- FindClusters(drugs_d158_H3K27me3, resolution = 0.1)

d158_entropy_sample <- map_dfr(inhib_samples, function(y){
    sample_idx <- drugs_d158_H3K27me3$inhib_annotation == y
    cluster_vec <- drugs_d158_H3K27me3$seurat_clusters[sample_idx]
    return(tibble(
        entropy = entropy(table(cluster_vec), method='MM')
    ))
}, .id='sample')


plot_df <- d158_entropy_sample %>% 
    mutate(
        perturb=str_replace(sample, '.+_(.+)_d.+', '\\1'),
        conc=str_replace(sample, '.+_.+_d.+-(\\d+)?', '\\1')
    ) %>% 
    mutate(conc=factor(ifelse(perturb=='DMSO', 'DMSO', conc), levels=c('DMSO', '1', '3', '10')))

ggplot(plot_df, aes(conc, entropy, fill=conc)) +
    geom_bar(stat='summary') +
    # geom_point(size=3) +
    stat_summary(geom = 'errorbar') +
    scale_y_continuous(limits=c(0,NA)) +
    scale_fill_manual(values=conc_colors)


drugs_d158_H3K27me3$age=str_replace(drugs_d158_H3K27me3$inhib_annotation, '.+_(d\\d+).*', '\\1')

drugs_d158_H3K27me3$conc=str_replace(str_replace(drugs_d158_H3K27me3$inhib_annotation, '.+-(\\d+)', '\\1'), '.*DMSO.*', 'DMSO')
concs <- unique(drugs_d158_H3K27me3$conc) %>% set_names()

d158_entropy_conc <- map_dfr(concs, function(y){
    sample_idx <- drugs_d158_H3K27me3$conc == y
    cluster_vec <- drugs_d158_H3K27me3$seurat_clusters[sample_idx]
    return(tibble(
        entropy = entropy(table(cluster_vec))
    ))
}, .id='conc')


plot_df <- d158_entropy_conc %>% 
    mutate(conc=factor(conc, levels=names(conc_colors)))

ggplot(plot_df, aes(conc=='DMSO', entropy, color=conc)) +
    # stat_summary(color='grey', size=0.1) +
    geom_point(size=3) +
    scale_y_continuous(limits=c(0,2)) +
    scale_color_manual(values=conc_colors)



#### Cytotrace for entropy calc ####
library(CytoTRACE)
library(entropy)

expr_mat <- drugs_d158_H3K27me3[['RNA']]@counts
expr_bin <- (drugs_d158_H3K27me3[['RNA']]@data > 0)*1

ctrace <- CytoTRACE::CytoTRACE(as.matrix(expr_mat), ncores=16)

ctrace_df <- ctrace$CytoTRACE %>% 
    enframe('cell', 'ctrace') %>% 
    column_to_rownames('cell')

drugs_d158_H3K27me3 <- AddMetaData(drugs_d158_H3K27me3, ctrace_df)

meta <- drugs_d158_H3K27me3@meta.data %>% 
    as_tibble(rownames='cell') 
    
ggplot(meta, aes(conc, ctrace, fill=age)) +
    geom_boxplot() 

feature_plot(drugs_d158_H3K27me3, features='ctrace', reduction='humap')



















#### A485 H3K27ac ####
H3K27ac_detected_early <- colnames(H3K27ac_early_detect)[colMaxs(H3K27ac_early_detect)>0.05]

H3K27ac_A485_15d_mat <- H3K27ac_A485_peak_score_mat[str_detect(rownames(H3K27ac_A485_peak_score_mat), '15d'), ]
H3K27ac_A485_18d_mat <- H3K27ac_A485_peak_score_mat[str_detect(rownames(H3K27ac_A485_peak_score_mat), '18d$'), ]
H3K27ac_A485_18d2_mat <- H3K27ac_A485_peak_score_mat[str_detect(rownames(H3K27ac_A485_peak_score_mat), '18d_rep2'), ]

H3K27ac_DMSO_15d_mat <- H3K27ac_DMSO_peak_score_mat[str_detect(rownames(H3K27ac_DMSO_peak_score_mat), '15d'), ]
H3K27ac_DMSO_18d_mat <- H3K27ac_DMSO_peak_score_mat[str_detect(rownames(H3K27ac_DMSO_peak_score_mat), '18d$'), ]
H3K27ac_DMSO_18d2_mat <- H3K27ac_DMSO_peak_score_mat[str_detect(rownames(H3K27ac_DMSO_peak_score_mat), '18d_rep2'), ]

H3K27ac_A485_15d_fc <- (t(log2(H3K27ac_A485_15d_mat+1)) - log2(H3K27ac_DMSO_15d_mat+1)) %>% 
    as_tibble(rownames='peak') %>% pivot_longer(!peak, names_to='sample', values_to='fc')
H3K27ac_A485_18d_fc <- (t(log2(H3K27ac_A485_18d_mat+1)) - log2(H3K27ac_DMSO_18d_mat+1)) %>% 
    as_tibble(rownames='peak') %>% pivot_longer(!peak, names_to='sample', values_to='fc')
H3K27ac_A485_18d2_fc <- (t(log2(H3K27ac_A485_18d2_mat+1)) - log2(H3K27ac_DMSO_18d2_mat+1)) %>% 
    as_tibble(rownames='peak') %>% pivot_longer(!peak, names_to='sample', values_to='fc')

H3K27ac_peak_genes <- Signac::ClosestFeature(marks$H3K27ac, StringToGRanges(H3K27ac_detected)) %>% 
    as_tibble()

H3K27ac_A485_fc <- bind_rows(H3K27ac_A485_18d2_fc, H3K27ac_A485_18d_fc, H3K27ac_A485_15d_fc) %>% 
    filter(peak%in%H3K27ac_detected_early) %>% 
    mutate(
        age=str_replace(sample, '.+A485_(\\d+d).*', '\\1'),
        conc=factor(as.numeric(str_replace(sample, '.+HW_(\\d+)-A485_\\d+d.*', '\\1'))),
        replicate=as.numeric(ifelse(str_detect(sample, 'rep2$'), '2', '1'))
    ) %>% inner_join(H3K27ac_peak_genes, by=c('peak'='query_region'))


H3K27ac_A485_mean_fc <- H3K27ac_A485_fc %>% 
    group_by(peak) %>% summarize(fc=mean(fc)) %>% arrange(desc(fc)) %>% 
    inner_join(H3K27ac_peak_genes, by=c('peak'='query_region')) %>% 
    mutate(peak=factor(peak, levels=unique(.$peak)))

H3K27ac_A485_fc$peak <- factor(H3K27ac_A485_fc$peak, levels=H3K27ac_A485_mean_fc$peak)

H3K27ac_A485_fc %>% write_tsv('data/drugs/drugs_CT_A485_H3K27ac_peaks_logfc.tsv')


#### Join with DE results ####
H3K27ac_A485_fc <- read_tsv('data/drugs/drugs_CT_A485_H3K27ac_peaks_logfc.tsv')
A485_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_global_inhib_de.tsv')
A485_global_inhib_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_concs_global_inhib_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))
A485_cluster_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_cluster_de.tsv')

top_global_de <- A485_global_de %>% filter(p_val_adj<1e-2, avg_log2FC>0.1) %>% pull(feature) %>% unique()

conc_colors <- c('1'='#90c47d', '3'='#0d5e30')

plot_df <- H3K27ac_A485_fc
p_s <- ggplot(plot_df, aes(peak, fc, color=conc)) +
    geom_point(size=0.01, alpha=0.3, shape=16) +
    geom_hline(yintercept = 0, size=0.1) +
    geom_line(data=H3K27ac_A485_mean_fc, mapping=aes(x=as.numeric(peak)), color='#34495e', size=0.2) +
    scale_color_manual(values=conc_colors) +
    scale_x_discrete(expand=c(0.02,0)) +
    scale_y_continuous(breaks=seq(-4,4,2), limits=c(-4,4)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    no_legend() +
    theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    labs(y='log2 fold change', x='Peak') 

plot_df <- H3K27ac_A485_mean_fc
p_hist <- ggplot(plot_df, aes(fc)) +
    geom_histogram(color='black', fill='grey', size=0.1) +
    scale_x_continuous(breaks=seq(-4,4,2), limits=c(-4,4)) +
    article_text() +
    theme_void() +
    coord_flip()

(p_s | p_hist) + plot_layout(widths=c(5,1))


plot_df <- inner_join(A485_global_inhib_de, H3K27ac_A485_mean_fc, by=c('feature'='gene_name')) %>% 
    mutate(
        conc=str_replace(inhib_annotation, '.+-(\\d+)', '\\1'),
        age=str_replace(inhib_annotation, '.+_(d\\d+)-\\d+', '\\1')
    ) %>% 
    group_by(conc, peak) %>% 
    summarize(fc=mean(avg_log2FC), feature=feature[1]) %>% 
    group_by(conc, feature) %>% 
    mutate(gene_num=paste0(feature, '_', row_number())) %>% 
    mutate(
        fc_clip=clip_abs(fc, 1),
        peak=factor(peak, levels=H3K27ac_A485_mean_fc$peak),
        conc=factor(conc, levels=rev(c(1,3,10)))
    ) %>% 
    arrange(peak) %>% mutate(gene_num=factor(gene_num, levels=unique(.$gene_num)))

p_heat <- ggplot(plot_df, aes(peak, conc, fill=fc_clip)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(bigrad(pals::brewer.rdbu, bias=1.5)[20:180])) + 
    scale_x_discrete(expand=c(0.02,0)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    no_x_text() +
    no_legend() +
    theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()
    ) 

layout <- '
AAAAAAB
AAAAAAB
AAAAAAB
AAAAAAB
AAAAAAB
CCCCCC#
'
    
p_s + p_hist + p_heat + plot_layout(design = layout) & no_margin()
ggsave('plots/paper/sfig4/fig4_A485_H3K27ac_detect05_peak_gene_fc_splot.pdf', width=10, height=6, unit='cm')


ggplot(plot_df, aes(gene_num, conc, fill=fc_clip)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(bigrad(pals::brewer.rdbu, bias=1.5)[20:180])) + 
    scale_x_discrete(expand=c(0.02,0)) +
    theme_rangeframe() + scale_axis_rangeframe() +
    article_text() +
    rotate_x_text(90)

ggsave('plots/paper/sfig4/fig4_A485_H3K27ac_detect05_gene_fc_labels.pdf', width=800, height=3, limitsize = FALSE)




#### UMAP for A485 ####

dim_plot(drugs_d158_H3K27ac, group.by=c('inhib_annotation', 'RNA_snn_res.0.2'), reduction='humap')

feature_plot(drugs_d158_H3K27ac, features=c('PAX6', 'VIM', 'NES', 'TOP2A', 'MKI67', 'KRT8', 'ITM2B'), reduction='humap')

filter(A485_cluster_de, celltype==6) %>% top_n(20, avg_log2FC)


plot_df <- A485_global_de %>% 
    filter(feature%in%H3K27ac_A485_fc$gene_name)

ggplot(plot_df, aes(avg_log2FC, -log10(p_val))) +
    geom_point(size=0.2)




#### DE module score for A485 ####

cluster_graph <- read_rds('data/all_RNA_cluster_graph.rds')
cluster_meta <- read_tsv('data/all_RNA_cluster_meta.tsv')

A485_global_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_global_inhib_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))
A485_global_inhib_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_concs_global_inhib_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))
A485_cluster_de <- read_tsv('data/drugs/results/diff_expression/drugs_d15_d18_A485_cluster_de.tsv') %>% mutate(padj=p.adjust(p_val, 'fdr'))

global_deg <- A485_global_de %>% filter(padj<1e-4, avg_log2FC<(-0.2)) %>% pull(feature) %>% unique()
global_deg <- A485_global_inhib_de %>% filter(padj<1e-4, avg_log2FC<(-0.2), inhib_annotation!='10_A485_d18-3') %>% pull(feature) %>% unique()

cluster_srt@active.assay <- 'RNA'
cluster_srt <- AddModuleScore(cluster_srt, features=list('global_DEG'=global_deg))
cluster_srt$RNA_global_DEG_score <- cluster_srt$Cluster1

cluster_srt@active.assay <- 'H3K27me3_RNA'
cluster_srt <- AddModuleScore(cluster_srt, features=list('global_DEG'=global_deg))
cluster_srt$H3K27me3_global_DEG_score <- cluster_srt$Cluster1

cluster_srt@active.assay <- 'H3K27ac_RNA'
cluster_srt <- AddModuleScore(cluster_srt, features=list('global_DEG'=global_deg))
cluster_srt$H3K27ac_global_DEG_score <- cluster_srt$Cluster1

cluster_srt@active.assay <- 'H3K4me3_RNA'
cluster_srt <- AddModuleScore(cluster_srt, features=list('global_DEG'=global_deg))
cluster_srt$H3K4me3_global_DEG_score <- cluster_srt$Cluster1

cluster_srt$bivalent_global_DEG_score <- scale(cluster_srt$H3K4me3_global_DEG_score) + scale(cluster_srt$H3K27me3_global_DEG_score)
cluster_srt$switch_global_DEG_score <- scale(cluster_srt$H3K27ac_global_DEG_score) + scale(cluster_srt$H3K27me3_global_DEG_score)

cluster_meta <- cluster_srt@meta.data %>% 
    as_tibble(rownames='cluster')

cluster_graph_meta <- cluster_graph %N>% 
    inner_join(cluster_meta)

ggraph(cluster_graph_meta, x=tree_pos, y=pseudotime_ranks) +
    geom_edge_diagonal(color='grey', alpha=0.3) +
    geom_node_point(aes(color=celltype), size=3) +
    scale_color_manual(values=pantone_celltype) +
    theme_void()

p1 <- ggraph(cluster_graph_meta, x=pseudotime_ranks, y=tree_pos) +
    geom_edge_diagonal(color='grey', alpha=0.3) +
    geom_node_point(aes(fill=RNA_global_DEG_score), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=grad(pals::brewer.orrd)) +
    scale_y_reverse() +
    theme_void()

p2 <- ggraph(cluster_graph_meta, x=pseudotime_ranks, y=tree_pos) +
    geom_edge_diagonal(color='grey', alpha=0.3) +
    geom_node_point(aes(fill=H3K27me3_global_DEG_score), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=grad(pals::brewer.blues)) +
    scale_y_reverse() +
    theme_void()

p3 <- ggraph(cluster_graph_meta, x=pseudotime_ranks, y=tree_pos) +
    geom_edge_diagonal(color='grey', alpha=0.3) +
    geom_node_point(aes(fill=H3K4me3_global_DEG_score), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=grad(pals::brewer.rdpu)) +
    scale_y_reverse() +
    theme_void()

p4 <- ggraph(cluster_graph_meta, x=pseudotime_ranks, y=tree_pos) +
    geom_edge_diagonal(color='grey', alpha=0.3) +
    geom_node_point(aes(fill=H3K27ac_global_DEG_score), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=grad(pals::brewer.ylgn)) +
    scale_y_reverse() +
    theme_void()

p1 / p2 / p3 / p4 + plot_layout(guides='collect')
# ggsave('plots/paper/fig4/fig4_global_deg_scores_tree.pdf', width=6, height=5)


p1 <- ggraph(cluster_graph_meta, x=pseudotime_ranks, y=tree_pos) +
    geom_edge_diagonal(color='grey', alpha=0.3) +
    geom_node_point(aes(fill=bivalent_global_DEG_score), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=grad(pals::brewer.bupu)) +
    scale_y_reverse() +
    theme_void()

p2 <- ggraph(cluster_graph_meta, x=pseudotime_ranks, y=tree_pos) +
    geom_edge_diagonal(color='grey', alpha=0.3) +
    geom_node_point(aes(fill=switch_global_DEG_score), size=4, shape=21, stroke=0.2) +
    scale_fill_gradientn(colors=grad(pals::brewer.bugn)) +
    scale_y_reverse() +
    theme_void()

p1 / p2
ggsave('plots/paper/fig4/fig4_global_deg_comb_scores_tree.pdf', width=6, height=2.6)



#### Plot cellrank results ####

cellrank_meta <- read_tsv('data/drugs/drugs_d15_d18_A395_cellrank_probs.tsv')
colnames(cellrank_meta)[1] <- 'cell'
cellrank_meta_df <- column_to_rownames(cellrank_meta, 'cell')

drugs_d158_H3K27me3 <- drugs_d158_H3K27me3[, cellrank_meta$cell]
drugs_d158_H3K27me3 <- AddMetaData(drugs_d158_H3K27me3, cellrank_meta_df)


#### DR in cellrank space ####
library(irlba)

cr_space <- cellrank_meta %>%
    select(cell, to_psc_ranks, to_nne_ranks, to_neural_crest_ranks, to_neurons_ranks) %>%
    column_to_rownames('cell') %>% as.matrix()

cr_pca <- pca(cr_space, n=2, to_df=T) %>%
    inner_join(cellrank_meta)

p1 <- ggplot(cr_pca, aes(PC1, PC2, color=to_psc_ranks)) +
    geom_point(size=0.5, alpha=0.8) +
    scale_color_gradientn(colors=gyylgnbu()) +
    theme_void()
p2 <- ggplot(cr_pca, aes(PC1, PC2, color=to_nne_ranks)) +
    geom_point(size=0.5, alpha=0.8) +
    scale_color_gradientn(colors=gyylgnbu()) +
    theme_void()
p1 | p2




#### Cellrank viz ####
meta <- drugs_d158_H3K27me3@meta.data %>% 
    as_tibble(rownames='cell')

cr_df <- cellrank_meta %>% 
    column_to_rownames('cell')

cr_space <- cellrank_meta %>% 
    column_to_rownames('cell') %>% as.matrix()


lin_order <- c('to_nne_ranks', 'to_neural_crest_ranks', 'to_neurons_ranks')
probs <- cr_space[,lin_order] 
probs <- probs / colMeans(probs)
probs <- probs / rowSums(probs)

angle_vec = seq(0, 2*pi, length.out=4)[1:3]
angle_vec_sin = cos(angle_vec)
angle_vec_cos = sin(angle_vec)

x = rowSums(t(apply(probs, 1, function(x)x*angle_vec_sin)))
y = rowSums(t(apply(probs, 1, function(x)x*angle_vec_cos)))

meta$C1 <- x
meta$C2 <- y

circle_dr <- cbind(x,y)
colnames(circle_dr) <- NULL

drugs_d158_H3K27me3[['circular']] <- CreateDimReducObject(circle_dr, key='CR_')



conc_colors <- c('1'='#A2D9DE', '3'='#3A90B0', '10'='#0B4771', 'DMSO'='#f78f1e')


drugs_d158_H3K27me3$conc <- str_replace(str_replace(drugs_d158_H3K27me3$inhib_annotation, '.+-(\\d+)', '\\1'), '.*DMSO.*', 'DMSO')

drugs_d158_H3K27me3 <- AddMetaData(drugs_d158_H3K27me3, cr_df)

dim_plot(drugs_d158_H3K27me3, reduction='circular', pt.size=1, split.by='conc') +
    scale_color_manual(values=cluster_colors)


dim_plot(drugs_d158_H3K27me3, reduction='circular', pt.size=1, group.by='conc') +
    scale_color_manual(values=conc_colors)


goi <- c('TFAP2A', 'POU5F1', 'DLX6', 'VIM', 'WLS', 'STMN2')
feature_plot(drugs_d158_H3K27me3, features=goi, reduction='circular', order=T) &
    scale_color_gradientn(colors=gyorrd())
ggsave('plots/paper/fig4/fig4_cellrank_features_circular.pdf', width=10, height=15)
ggsave('plots/paper/fig4/fig4_cellrank_features_circular.png', width=10, height=15)


feature_plot(drugs_d158_H3K27me3, features=c('to_nne', 'to_neural_crest', 'to_neurons'), reduction='circular')
feature_plot(drugs_d158_H3K27me3, features=c('to_nne_ranks', 'to_neural_crest_ranks', 'to_neurons_ranks'), reduction='circular')
feature_plot(drugs_d158_H3K27me3, features=c('to_psc', 'to_nne', 'to_neural_crest', 'to_neurons'), reduction='humap')





#### Density plots on circular ####
source('~/scripts/single_cell/cell_density_diff.R')

# knn <- RANN::nn2(Embeddings(drugs_d158_H3K27me3, 'circular')[,1:2], k = 100)
# knngraph <- FindNeighbors(Embeddings(drugs_d158_H3K27me3, 'circular')[,1:2], k.param = 100)

knn <- RANN::nn2(probs, k = 100)
knngraph <- FindNeighbors(probs, k.param = 100)

knn_enrich <- get_knn_cmh_enrichment(knn, drugs_d158_H3K27me3$conc != 'DMSO', drugs_d158_H3K27me3$conc == 'DMSO', stratums = drugs_d158_H3K27me3$orig.ident)
smooth_enrich <- random_walk_with_restart(knn_enrich, knn_mat = as(knngraph$snn, 'Matrix'), num_rounds = 100, alpha = 0.1)

plot_df <- meta %>% 
    inner_join(enframe(smooth_enrich, 'cell', 'enrich'))

p <- ggplot(plot_df, aes(C1, C2, z=enrich, color=enrich)) +
    # stat_summary_hex(bins=50, fun='mean') +
    geom_point(size=0.5) +
    theme_void() 

pb <- ggplot_build(p)
clim <- 18
# p + scale_fill_gradientn(colors=rev(bigrad(pals::brewer.rdbu, bias=1.5)), limits=c(-clim,clim)) 
p + scale_color_gradientn(colors=rev(bigrad(pals::brewer.rdbu, bias=1.5)), limits=c(-clim,clim))  +
    labs(color='Local\nenrichment')

ggsave('plots/paper/fig4/fig4_cellrank_treatment_density_circular.png', width=6, height=4)
ggsave('plots/paper/fig4/fig4_cellrank_treatment_density_circular.pdf', width=6, height=4)
















