source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')


setwd('~/projects/cutntag/')

library(SeuratDisk)


select <- dplyr::select


#### Read stuff ####
marks <- read_rds('data/all_marks_list_v3.1clustered.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.1annotated.rds')
good_genes <- read_tsv('~/resources/gene_sets/protein_coding_genes_goodlist.txt', col_names = F)$X1


#### Summarize clusters ####
rna <- Pando::aggregate_assay(rna, assay='RNA', slot='data', group_name='clusters')
marks <- map(marks, Pando::aggregate_assay, assay='cRNA', slot='data', group_name='clusters')
marks <- map(marks, Pando::aggregate_assay, assay='peaks', slot='data', group_name='clusters')

cluster_meta <- rna@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    distinct(clusters, stage, celltype_jf)



#### Write out matrices for bipartite matching ####
#### Select markers and write summary matrices for each stage ####
stages <- unique(rna$stage)

for (s in stages){
    H3K27ac_ <- subset(marks$H3K27ac, stage==s)
    H3K4me3_ <- subset(marks$H3K4me3, stage==s)
    H3K27me3_ <- subset(marks$H3K27me3, stage==s)
    rna_ <- subset(rna, stage==s)
    
    rna_clusters_use <- unique(cluster_meta$clusters)
    
    genes_use <- intersect(rownames(marks$H3K27ac[['cRNA']]), rownames(marks$H3K4me3[['cRNA']])) %>% intersect(rownames(marks$H3K27me3[['cRNA']]))
    rna_de <- de(rna_, groups='clusters')
    rna_markers <- rna_de %>% filter(padj<1e-5, feature%in%genes_use, !feature%in%cc_genes_all, feature%in%good_genes) %>% group_by(group) %>% 
        top_n(30, fc) %>% pull(feature) %>% unique()
    
    H3K27ac_de <- de(H3K27ac_, groups='clusters', assay='cRNA')
    
    H3K27ac_markers <- H3K27ac_de %>% filter(padj<1e-5, feature%in%rna_markers) %>% group_by(group) %>% 
        top_n(5, fc) %>% pull(feature) %>% unique()
    
    H3K27ac_use <- t(H3K27ac_@assays$cRNA@misc$summary$clusters[, H3K27ac_markers])
    H3K27ac_use <- H3K27ac_use[, str_detect(colnames(H3K27ac_use), s)]
    colnames(H3K27ac_use) <- paste0('H3K27ac_', colnames(H3K27ac_use))
    
    rna_use <- t(rna_@assays$RNA@misc$summary$clusters[rna_clusters_use, H3K27ac_markers])
    rna_use <- rna_use[, str_detect(colnames(rna_use), s)]
    colnames(rna_use) <- paste0('rna_', colnames(rna_use))
    
    H3K27ac_rna <- cbind(H3K27ac_use, rna_use)
    H3K27ac_rna_srt <- CreateSeuratObject(H3K27ac_rna)
    H3K27ac_rna_srt$modality <- c(rep('H3K27ac', ncol(H3K27ac_use)), rep('rna', ncol(rna_use)))
    
    SaveH5Seurat(H3K27ac_rna_srt, filename = paste0('data/RNA/integration/H3K27ac_rna_clusters', s, '.h5seurat'), overwrite=T)
    Convert(source = paste0('data/RNA/integration/H3K27ac_rna_clusters', s, '.h5seurat'), dest='h5ad', overwrite=T)
    
    
    H3K4me3_de <- de(H3K4me3_, groups='clusters', assay='cRNA')
    H3K4me3_markers <- H3K4me3_de %>% filter(padj<1e-5, feature%in%rna_markers) %>% group_by(group) %>% 
        top_n(5, fc) %>% pull(feature) %>% unique()
    
    H3K4me3_use <- t(H3K4me3_@assays$cRNA@misc$summary$clusters[, H3K4me3_markers])
    H3K4me3_use <- H3K4me3_use[, str_detect(colnames(H3K4me3_use), s)]
    colnames(H3K4me3_use) <- paste0('H3K4me3_', colnames(H3K4me3_use))
    
    rna_use <- t(rna_@assays$RNA@misc$summary$clusters[rna_clusters_use, H3K4me3_markers])
    rna_use <- rna_use[, str_detect(colnames(rna_use), s)]
    colnames(rna_use) <- paste0('rna_', colnames(rna_use))
    
    H3K4me3_rna <- cbind(H3K4me3_use, rna_use)
    H3K4me3_rna_srt <- CreateSeuratObject(H3K4me3_rna)
    H3K4me3_rna_srt$modality <- c(rep('H3K4me3', ncol(H3K4me3_use)), rep('rna', ncol(rna_use)))
    
    SaveH5Seurat(H3K4me3_rna_srt, filename = paste0('data/RNA/integration/H3K4me3_rna_clusters', s, '.h5seurat'), overwrite=T)
    Convert(source = paste0('data/RNA/integration/H3K4me3_rna_clusters', s, '.h5seurat'), dest='h5ad', overwrite=T)
    
    
    
    H3K27me3_de <- de(H3K27me3_, groups='clusters', assay='cRNA')
    H3K27me3_markers <- H3K27me3_de %>% filter(padj<1e-5, feature%in%rna_markers) %>% group_by(group) %>% 
        top_n(5, -fc) %>% pull(feature) %>% unique()
    
    H3K27me3_use <- t(H3K27me3_@assays$cRNA@misc$summary$clusters[, H3K27me3_markers])
    H3K27me3_use <- H3K27me3_use[, str_detect(colnames(H3K27me3_use), s)]
    colnames(H3K27me3_use) <- paste0('H3K27me3_', colnames(H3K27me3_use))
    
    rna_use <- t(rna_@assays$RNA@misc$summary$clusters[rna_clusters_use, H3K27me3_markers])
    rna_use <- rna_use[, str_detect(colnames(rna_use), s)]
    colnames(rna_use) <- paste0('rna_', colnames(rna_use))
    
    H3K27me3_rna <- cbind(H3K27me3_use, rna_use)
    H3K27me3_rna_srt <- CreateSeuratObject(H3K27me3_rna)
    H3K27me3_rna_srt$modality <- c(rep('H3K27me3', ncol(H3K27me3_use)), rep('rna', ncol(rna_use)))
    
    SaveH5Seurat(H3K27me3_rna_srt, filename = paste0('data/RNA/integration/H3K27me3_rna_clusters', s, '.h5seurat'), overwrite=T)
    Convert(source = paste0('data/RNA/integration/H3K27me3_rna_clusters', s, '.h5seurat'), dest='h5ad', overwrite=T)
}





#### Get bipartite matches ####
##### H3K27ac ####
H3K27ac_match_files <- list(
    EB='data/RNA/integration/H3K27ac_rna_EB_matches.tsv', 
    mid='data/RNA/integration/H3K27ac_rna_mid_matches.tsv',
    late='data/RNA/integration/H3K27ac_rna_late_matches.tsv',
    mo8='data/RNA/integration/H3K27ac_rna_mo8_matches.tsv',
    retina='data/RNA/integration/H3K27ac_rna_retina_matches.tsv'
)
rna_H3K27ac_matches <- map_dfr(H3K27ac_match_files, read_tsv, .id='stage') %>% 
    mutate(
        rna_cluster = str_replace(rna, 'rna_(.+\\d+)', '\\1'),
        H3K27ac_cluster = str_replace(H3K27ac, 'H3K27ac_(.+\\d+)', '\\1')
    )

umap_coords <- marks$H3K27ac[['umap']]@cell.embeddings %>% 
    as_tibble(rownames='cell')
cluster_coords <- Pando::aggregate_matrix(marks$H3K27ac[['umap']]@cell.embeddings, as.character(marks$H3K27ac$clusters)) %>% 
    {colnames(.) <- c('UMAP1', 'UMAP2');.} %>% 
    as_tibble(rownames='clusters') 

rna_H3K27ac_assign <- rna_H3K27ac_matches %>% 
    inner_join(cluster_coords, by=c('H3K27ac_cluster'='clusters')) %>% 
    full_join(cluster_meta, by=c('rna_cluster'='clusters'))

marks$H3K27ac$celltype_jf <- rna_H3K27ac_assign$celltype_jf[match(marks$H3K27ac$clusters, rna_H3K27ac_assign$H3K27ac_cluster)]

p1 <- dim_plot(marks$H3K27ac, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors, na.value = 'grey') + ggtitle('H3K27ac')
p2 <- dim_plot(marks$H3K27ac, group.by=c('age')) +
    scale_color_manual(values=age_colors2)
p_H3K27ac <- p1 | p2
p_H3K27ac



##### H3K27me3 ####
H3K27me3_match_files <- list(
    EB='data/RNA/integration/H3K27me3_rna_EB_matches.tsv', 
    mid='data/RNA/integration/H3K27me3_rna_mid_matches.tsv',
    late='data/RNA/integration/H3K27me3_rna_late_matches.tsv',
    mo8='data/RNA/integration/H3K27me3_rna_mo8_matches.tsv',
    retina='data/RNA/integration/H3K27me3_rna_retina_matches.tsv'
)
rna_H3K27me3_matches <- map_dfr(H3K27me3_match_files, read_tsv, .id='stage') %>% 
    mutate(
        rna_cluster = str_replace(rna, 'rna_(.+\\d+)', '\\1'),
        H3K27me3_cluster = str_replace(H3K27me3, 'H3K27me3_(.+\\d+)', '\\1')
    )

umap_coords <- marks$H3K27me3[['umap']]@cell.embeddings %>% 
    as_tibble(rownames='cell')
cluster_coords <- Pando::aggregate_matrix(marks$H3K27me3[['umap']]@cell.embeddings, as.character(marks$H3K27me3$clusters)) %>% 
    {colnames(.) <- c('UMAP1', 'UMAP2');.} %>% 
    as_tibble(rownames='clusters') 

rna_H3K27me3_assign <- rna_H3K27me3_matches %>% 
    inner_join(cluster_coords, by=c('H3K27me3_cluster'='clusters')) %>% 
    full_join(cluster_meta, by=c('rna_cluster'='clusters', 'stage'))

marks$H3K27me3$celltype_jf <- rna_H3K27me3_assign$celltype_jf[match(marks$H3K27me3$clusters, rna_H3K27me3_assign$H3K27me3_cluster)]

p1 <- dim_plot(marks$H3K27me3, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors, na.value = 'grey') + ggtitle('H3K27me3')
p2 <- dim_plot(marks$H3K27me3, group.by=c('age')) +
    scale_color_manual(values=age_colors2)
p_H3K27me3 <- p1 | p2
p_H3K27me3



##### H3K4me3 ####
H3K4me3_match_files <- list(
    EB='data/RNA/integration/H3K4me3_rna_EB_matches.tsv', 
    mid='data/RNA/integration/H3K4me3_rna_mid_matches.tsv',
    late='data/RNA/integration/H3K4me3_rna_late_matches.tsv',
    mo8='data/RNA/integration/H3K4me3_rna_mo8_matches.tsv',
    retina='data/RNA/integration/H3K4me3_rna_retina_matches.tsv'
)
rna_H3K4me3_matches <- map_dfr(H3K4me3_match_files, read_tsv, .id='stage') %>% 
    mutate(
        rna_cluster = str_replace(rna, 'rna_(.+\\d+)', '\\1'),
        H3K4me3_cluster = str_replace(H3K4me3, 'H3K4me3_(.+\\d+)', '\\1')
    )

umap_coords <- marks$H3K4me3[['umap']]@cell.embeddings %>% 
    as_tibble(rownames='cell')
cluster_coords <- Pando::aggregate_matrix(marks$H3K4me3[['umap']]@cell.embeddings, as.character(marks$H3K4me3$clusters)) %>% 
    {colnames(.) <- c('UMAP1', 'UMAP2');.} %>% 
    as_tibble(rownames='clusters') 

rna_H3K4me3_assign <- rna_H3K4me3_matches %>% 
    inner_join(cluster_coords, by=c('H3K4me3_cluster'='clusters')) %>% 
    full_join(cluster_meta, by=c('rna_cluster'='clusters', 'stage'))

marks$H3K4me3$celltype_jf <- rna_H3K4me3_assign$celltype_jf[match(marks$H3K4me3$clusters, rna_H3K4me3_assign$H3K4me3_cluster)]

p1 <- dim_plot(marks$H3K4me3, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors, na.value='grey') + ggtitle('H3K4me3')
p2 <- dim_plot(marks$H3K4me3, group.by=c('age')) +
    scale_color_manual(values=age_colors2)
p_H3K4me3 <- p1 | p2
p_H3K4me3

p_H3K27ac / p_H3K4me3 / p_H3K27me3
ggsave('plots/RNA/integration/highres_bipartite_integration_umap.png', width=12, height=18)


#### Write combined ####
all_matches <- bind_rows(
    'H3K27ac'=rename(rna_H3K27ac_matches, 'mark_cluster'=5),
    'H3K27me3'=rename(rna_H3K27me3_matches, 'mark_cluster'=5),
    'H3K4me3'=rename(rna_H3K4me3_matches, 'mark_cluster'=5),
    .id='mark'
)

all_matches$mark_cluster[is.na(all_matches$mark_cluster)] <- 'none'
all_matches %>% write_tsv('data/RNA/integration/cluster_bipartite_matches.tsv')
# all_matches <- read_tsv('data/RNA/integration/cluster_bipartite_matches.tsv')




#### Match clusters to RNA based correlation ####
genes_use <- intersect(rownames(marks$H3K27ac[['cRNA']]), rownames(marks$H3K4me3[['cRNA']])) %>% intersect(rownames(marks$H3K27me3[['cRNA']]))
rna_de <- de(rna, groups='clusters')
rna_markers <- rna_de %>% filter(padj<0.01, feature%in%genes_use, !feature%in%cc_genes_all, feature%in%good_genes) %>% group_by(group) %>% 
    top_n(30, fc) %>% pull(feature) %>% unique()

H3K27ac_de <- de(marks$H3K27ac, groups='clusters', assay='cRNA')
H3K27ac_markers <- H3K27ac_de %>% filter(padj<0.01, feature%in%rna_markers) %>% group_by(group) %>% 
    top_n(5, fc) %>% pull(feature) %>% unique()

H3K27me3_de <- de(marks$H3K27me3, groups='clusters', assay='cRNA')
H3K27me3_markers <- H3K27me3_de %>% filter(padj<0.01, feature%in%rna_markers) %>% group_by(group) %>% 
    top_n(5, -fc) %>% pull(feature) %>% unique()

H3K4me3_de <- de(marks$H3K4me3, groups='clusters', assay='cRNA')
H3K4me3_markers <- H3K4me3_de %>% filter(padj<0.01, feature%in%rna_markers) %>% group_by(group) %>% 
    top_n(5, fc) %>% pull(feature) %>% unique()

rna_H3K27ac_cor <- Pando::sparse_cor(t(rna@assays$RNA@misc$summary$clusters)[H3K27ac_markers, ], t(marks$H3K27ac@assays$cRNA@misc$summary$clusters)[H3K27ac_markers, ])
pheatmap::pheatmap(t(rna_H3K27ac_cor), width=20, height=15)

rna_H3K27me3_cor <- Pando::sparse_cor(t(rna@assays$RNA@misc$summary$clusters)[H3K27me3_markers, ], t(marks$H3K27me3@assays$cRNA@misc$summary$clusters)[H3K27me3_markers, ])
pheatmap::pheatmap(t(-rna_H3K27me3_cor), width=20, height=15)

rna_H3K4me3_cor <- Pando::sparse_cor(t(rna@assays$RNA@misc$summary$clusters)[H3K4me3_markers, ], t(marks$H3K4me3@assays$cRNA@misc$summary$clusters)[H3K4me3_markers, ])
pheatmap::pheatmap(t(rna_H3K4me3_cor), width=20, height=15)

cluster_corr_list <- list(
    H3K27ac = rna_H3K27ac_cor,
    H3K27me3 = rna_H3K27me3_cor,
    H3K4me3 = rna_H3K4me3_cor
) 
cluster_corr_list %>% write_rds('data/RNA/integration/cluster_correlation_list.rds')
cluster_corr_list <- read_rds('data/RNA/integration/cluster_correlation_list.rds')


#### Merge with metadata and plot ####
#### H3K27ac ####
H3K27ac_clust <- rna_H3K27ac_cor %>% t() %>% dist() %>% hclust()
H3K27ac_order <- H3K27ac_clust$labels[H3K27ac_clust$order]
rna_clust <- rna_H3K27ac_cor %>% dist() %>% hclust()
rna_order <- rna_clust$labels[rna_clust$order]

rna_H3K27ac_cor_df <- rna_H3K27ac_cor %>% 
    as_tibble(rownames='rna_clusters') %>% 
    pivot_longer(!rna_clusters, names_to='H3K27ac_clusters', values_to='corr') %>% 
    inner_join(cluster_meta, by=c('rna_clusters'='clusters')) %>% 
    mutate(
        H3K27ac_clusters = factor(H3K27ac_clusters, levels=H3K27ac_order),
        rna_clusters = factor(rna_clusters, levels=unique(.$rna_clusters)),
        H3K27ac_stage = factor(str_replace(H3K27ac_clusters, '(\\w+)_\\d+', '\\1'), levels=c('EB', 'mid', 'late', 'retina', 'mo8')),
        rna_stage = factor(str_replace(rna_clusters, '(\\w+)_\\d+', '\\1'), levels=c('EB', 'mid', 'late', 'retina', 'mo8'))
    ) 

ggplot(rna_H3K27ac_cor_df, aes(rna_clusters, H3K27ac_clusters, fill=corr)) +
    geom_tile() +
    scale_fill_gradientn(limits=c(-0.6, 0.6), colors=rev(pals::brewer.rdylbu(100))) +
    facet_grid(H3K27ac_stage~celltype_jf, scales = 'free', space='free') +
    rotate_x_text(40)

ggsave('plots/RNA/integration/rna_H3K27ac_cor_heatmap.pdf', width=24, height=14)
ggsave('plots/RNA/integration/rna_H3K27ac_cor_heatmap.png', width=24, height=14)


rna_H3K27ac_cor_max <- rna_H3K27ac_cor_df %>% 
    filter(H3K27ac_stage==rna_stage) %>% 
    group_by(H3K27ac_clusters) %>% 
    filter(corr==max(corr))

marks$H3K27ac$celltype_jf <- rna_H3K27ac_cor_max$celltype_jf[match(marks$H3K27ac$clusters, rna_H3K27ac_cor_max$H3K27ac_clusters)]
marks$H3K27ac$max_cor <- rna_H3K27ac_cor_max$corr[match(marks$H3K27ac$clusters, rna_H3K27ac_cor_max$H3K27ac_clusters)]

p1 <- dim_plot(marks$H3K27ac, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors) + ggtitle('H3K27ac')
p2 <- dim_plot(marks$H3K27ac, group.by=c('age')) +
    scale_color_manual(values=age_colors2)
p_H3K27ac <- p1 | p2
p_H3K27ac


#### H3K4me3 #### 
H3K4me3_clust <- rna_H3K4me3_cor %>% t() %>% dist() %>% hclust()
H3K4me3_order <- H3K4me3_clust$labels[H3K4me3_clust$order]
rna_clust <- rna_H3K4me3_cor %>% dist() %>% hclust()
rna_order <- rna_clust$labels[rna_clust$order]

rna_H3K4me3_cor_df <- rna_H3K4me3_cor %>% 
    as_tibble(rownames='rna_clusters') %>% 
    pivot_longer(!rna_clusters, names_to='H3K4me3_clusters', values_to='corr') %>% 
    inner_join(cluster_meta, by=c('rna_clusters'='clusters')) %>% 
    mutate(
        H3K4me3_clusters = factor(H3K4me3_clusters, levels=H3K4me3_order),
        rna_clusters = factor(rna_clusters, levels=unique(.$rna_clusters)),
        H3K4me3_stage = factor(str_replace(H3K4me3_clusters, '(\\w+)_\\d+', '\\1'), levels=c('EB', 'mid', 'late', 'retina', 'mo8')),
        rna_stage = factor(str_replace(rna_clusters, '(\\w+)_\\d+', '\\1'), levels=c('EB', 'mid', 'late', 'retina', 'mo8')),
    ) 

ggplot(rna_H3K4me3_cor_df, aes(rna_clusters, H3K4me3_clusters, fill=corr)) +
    geom_tile() +
    scale_fill_gradientn(limits=c(-0.8, 0.8), colors=rev(pals::brewer.rdylbu(100))) +
    facet_grid(H3K4me3_stage~celltype_jf, scales = 'free', space='free') +
    rotate_x_text(40)

ggsave('plots/RNA/integration/rna_H3K4me3_cor_heatmap.pdf', width=24, height=18)
ggsave('plots/RNA/integration/rna_H3K4me3_cor_heatmap.png', width=24, height=18)



rna_H3K4me3_cor_max <- rna_H3K4me3_cor_df %>% 
    filter(H3K4me3_stage==rna_stage) %>% 
    group_by(H3K4me3_clusters) %>% 
    filter(corr==max(corr))

marks$H3K4me3$celltype_jf <- rna_H3K4me3_cor_max$celltype_jf[match(marks$H3K4me3$clusters, rna_H3K4me3_cor_max$H3K4me3_clusters)]
marks$H3K4me3$max_cor <- rna_H3K4me3_cor_max$corr[match(marks$H3K4me3$clusters, rna_H3K4me3_cor_max$H3K4me3_clusters)]

p1 <- dim_plot(marks$H3K4me3, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors) + ggtitle('H3K4me3')
p2 <- dim_plot(marks$H3K4me3, group.by=c('age')) +
    scale_color_manual(values=age_colors2)
p_H3K4me3 <- p1 | p2
p_H3K4me3



#### H3K27me3 #### 
H3K27me3_clust <- rna_H3K27me3_cor %>% t() %>% dist() %>% hclust()
H3K27me3_order <- H3K27me3_clust$labels[H3K27me3_clust$order]
rna_clust <- rna_H3K27me3_cor %>% dist() %>% hclust()
rna_order <- rna_clust$labels[rna_clust$order]

rna_H3K27me3_cor_df <- rna_H3K27me3_cor %>% 
    as_tibble(rownames='rna_clusters') %>% 
    pivot_longer(!rna_clusters, names_to='H3K27me3_clusters', values_to='corr') %>% 
    inner_join(cluster_meta, by=c('rna_clusters'='clusters')) %>% 
    mutate(
        H3K27me3_clusters = factor(H3K27me3_clusters, levels=H3K27me3_order),
        rna_clusters = factor(rna_clusters, levels=unique(.$rna_clusters)),
        H3K27me3_stage = factor(str_replace(H3K27me3_clusters, '(\\w+)_\\d+', '\\1'), levels=c('EB', 'mid', 'late', 'retina', 'mo8')),
        rna_stage = factor(str_replace(rna_clusters, '(\\w+)_\\d+', '\\1'), levels=c('EB', 'mid', 'late', 'retina', 'mo8')),
    ) 

ggplot(rna_H3K27me3_cor_df, aes(rna_clusters, H3K27me3_clusters, fill=corr)) +
    geom_tile() +
    scale_fill_gradientn(limits=c(-0.5, 0.5), colors=rev(pals::brewer.rdylbu(100))) +
    facet_grid(H3K27me3_stage~celltype_jf, scales = 'free', space='free') +
    rotate_x_text(40)

ggsave('plots/RNA/integration/rna_H3K27me3_cor_heatmap.pdf', width=24, height=14)
ggsave('plots/RNA/integration/rna_H3K27me3_cor_heatmap.png', width=24, height=14)

rna_H3K27me3_cor_max <- rna_H3K27me3_cor_df %>% 
    filter(H3K27me3_stage==rna_stage) %>% 
    group_by(H3K27me3_clusters) %>% 
    filter(corr==min(corr))

marks$H3K27me3$celltype_jf <- rna_H3K27me3_cor_max$celltype_jf[match(marks$H3K27me3$clusters, rna_H3K27me3_cor_max$H3K27me3_clusters)]
marks$H3K27me3$max_cor <- rna_H3K27me3_cor_max$corr[match(marks$H3K27me3$clusters, rna_H3K27me3_cor_max$H3K27me3_clusters)]

p1 <- dim_plot(marks$H3K27me3, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors) + ggtitle('H3K27me3')
p2 <- dim_plot(marks$H3K27me3, group.by=c('age')) +
    scale_color_manual(values=age_colors2)
p_H3K27me3 <- p1 | p2
p_H3K27me3

p_H3K27ac / p_H3K4me3 / p_H3K27me3
ggsave('plots/RNA/integration/highres_maxcorr_integration_umap.png', width=12, height=18)


#### Write correlation df ####

all_corr_df <- bind_rows(
    'H3K27ac'=rna_H3K27ac_cor_df,
    'H3K4me3'=rna_H3K4me3_cor_df,
    'H3K27me3'=rna_H3K27me3_cor_df,
    .id='mark'
)

all_corr_df %>% write_tsv('data/RNA/integration/cluster_corr_df.tsv')
# all_corr_df <- read_tsv('data/RNA/integration/cluster_corr_df.tsv')








#### Format and write stuff ####
#### Highres cluster multimodal object ####
genes_use <- intersect(rownames(marks$H3K27ac[['cRNA']]), rownames(marks$H3K4me3[['cRNA']])) %>% intersect(rownames(marks$H3K27me3[['cRNA']])) %>% intersect(rownames(rna[['RNA']]))
rna_clusters <- rna@assays$RNA@misc$summary$clusters[, genes_use]
H3K27ac_rna_clust <- marks$H3K27ac@assays$cRNA@misc$summary$clusters[, genes_use]
H3K4me3_rna_clust <- marks$H3K4me3@assays$cRNA@misc$summary$clusters[, genes_use]
H3K27me3_rna_clust <- marks$H3K27me3@assays$cRNA@misc$summary$clusters[, genes_use]

H3K27ac_peaks_clust <- marks$H3K27ac@assays$peaks@misc$summary$clusters
H3K4me3_peaks_clust <- marks$H3K4me3@assays$peaks@misc$summary$clusters
H3K27me3_peaks_clust <- marks$H3K27me3@assays$peaks@misc$summary$clusters

H3K27ac_npeaks_clust <- marks$H3K27ac@assays$peaks@misc$summary$clusters
H3K4me3_npeaks_clust <- marks$H3K4me3@assays$peaks@misc$summary$clusters
H3K27me3_npeaks_clust <- marks$H3K27me3@assays$peaks@misc$summary$clusters


#### Select matched clusters ####
H3K27ac_matches_matched <- filter(all_corr_df, mark=='H3K27ac') %>% 
    inner_join(filter(all_matches, mark=='H3K27ac'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K27ac_clusters'='mark_cluster')) %>% 
    group_by(rna_clusters) %>% 
    filter(corr==max(corr))

H3K27ac_matches_missing <- filter(all_corr_df, mark=='H3K27ac') %>% 
    left_join(filter(all_matches, mark=='H3K27ac'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K27ac_clusters'='mark_cluster')) %>% 
    filter(!rna_clusters%in%H3K27ac_matches_matched$rna_clusters) %>% 
    filter(H3K27ac_stage==rna_stage) %>% 
    group_by(rna_clusters) %>% 
    filter(corr==max(corr))

H3K27ac_matches_final <- bind_rows(H3K27ac_matches_matched, H3K27ac_matches_missing) %>% 
    rename('mark_clusters'='H3K27ac_clusters')

H3K27me3_matches_matched <- filter(all_corr_df, mark=='H3K27me3') %>% 
    inner_join(filter(all_matches, mark=='H3K27me3'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K27me3_clusters'='mark_cluster')) %>% 
    group_by(rna_clusters) %>% 
    filter(corr==min(corr))

H3K27me3_matches_missing <- filter(all_corr_df, mark=='H3K27me3') %>% 
    left_join(filter(all_matches, mark=='H3K27me3'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K27me3_clusters'='mark_cluster')) %>% 
    filter(!rna_clusters%in%H3K27me3_matches_matched$rna_clusters) %>% 
    filter(H3K27me3_stage==rna_stage) %>% 
    group_by(rna_clusters) %>% 
    filter(corr==min(corr))

H3K27me3_matches_final <- bind_rows(H3K27me3_matches_matched, H3K27me3_matches_missing) %>% 
    rename('mark_clusters'='H3K27me3_clusters')

H3K4me3_matches_matched <- filter(all_corr_df, mark=='H3K4me3') %>% 
    inner_join(filter(all_matches, mark=='H3K4me3'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K4me3_clusters'='mark_cluster')) %>% 
    group_by(rna_clusters) %>% 
    filter(corr==max(corr))

H3K4me3_matches_missing <- filter(all_corr_df, mark=='H3K4me3') %>% 
    left_join(filter(all_matches, mark=='H3K4me3'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K4me3_clusters'='mark_cluster')) %>% 
    filter(!rna_clusters%in%H3K4me3_matches_matched$rna_clusters) %>% 
    filter(H3K4me3_stage==rna_stage) %>% 
    group_by(rna_clusters) %>% 
    filter(corr==max(corr))

H3K4me3_matches_final <- bind_rows(H3K4me3_matches_matched, H3K4me3_matches_missing) %>% 
    rename('mark_clusters'='H3K4me3_clusters')

all_matches_final <- bind_rows(H3K27ac_matches_final, H3K27me3_matches_final, H3K4me3_matches_final) %>% 
    dplyr::select(!c(H3K27ac, rna, H3K27me3, H3K4me3, H3K4me3_clusters, H3K27me3_clusters, H3K27ac_clusters, corr)) %>% 
    pivot_wider(names_from = mark, values_from = mark_clusters)

all_matches_final %>% write_tsv('data/RNA/integration/cluster_matches_rna.tsv')



#### Get cluster matches for all modalities ####
H3K4me3_matched <- filter(all_corr_df, mark=='H3K4me3') %>% 
    inner_join(filter(all_matches, mark=='H3K4me3'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K4me3_clusters'='mark_cluster')) %>% 
    group_by(H3K4me3_clusters) %>% 
    filter(corr==max(corr))

H3K4me3_missing <- filter(all_corr_df, mark=='H3K4me3') %>% 
    left_join(filter(all_matches, mark=='H3K4me3'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K4me3_clusters'='mark_cluster')) %>% 
    filter(!H3K4me3_clusters%in%H3K4me3_matched$H3K4me3_clusters) %>% 
    filter(H3K4me3_stage==rna_stage) %>%
    group_by(H3K4me3_clusters) %>% 
    filter(corr==max(corr))

H3K4me3_final <- bind_rows(H3K4me3_matched, H3K4me3_missing) %>% 
    rename('mark_clusters'='H3K4me3_clusters') %>% 
    select(!c(H3K27ac, rna, H3K27me3, H3K4me3, H3K27me3_clusters, H3K27ac_clusters)) 

H3K4me3_final %>% write_tsv('data/RNA/integration/cluster_matches_H3K4me3.tsv')



H3K27me3_matched <- filter(all_corr_df, mark=='H3K27me3') %>% 
    inner_join(filter(all_matches, mark=='H3K27me3'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K27me3_clusters'='mark_cluster')) %>% 
    group_by(H3K27me3_clusters) %>% 
    filter(corr==min(corr))

H3K27me3_missing <- filter(all_corr_df, mark=='H3K27me3') %>% 
    left_join(filter(all_matches, mark=='H3K27me3'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K27me3_clusters'='mark_cluster')) %>% 
    filter(!H3K27me3_clusters%in%H3K27me3_matched$H3K27me3_clusters) %>% 
    filter(H3K27me3_stage==rna_stage) %>%
    group_by(H3K27me3_clusters) %>% 
    filter(corr==min(corr))


H3K27me3_final <- bind_rows(H3K27me3_matched, H3K27me3_missing) %>% 
    rename('mark_clusters'='H3K27me3_clusters') %>% 
    select(!c(H3K27ac, rna, H3K27me3, H3K4me3, H3K4me3_clusters, H3K27ac_clusters)) 


H3K27me3_final %>% write_tsv('data/RNA/integration/cluster_matches_H3K27me3.tsv')



H3K27ac_matched <- filter(all_corr_df, mark=='H3K27ac') %>% 
    inner_join(filter(all_matches, mark=='H3K27ac'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K27ac_clusters'='mark_cluster')) %>% 
    group_by(H3K27ac_clusters) %>% 
    filter(corr==max(corr))

H3K27ac_missing <- filter(all_corr_df, mark=='H3K27ac') %>% 
    left_join(filter(all_matches, mark=='H3K27ac'), by=c('mark', 'rna_clusters'='rna_cluster', 'H3K27ac_clusters'='mark_cluster')) %>% 
    filter(!H3K27ac_clusters%in%H3K27ac_matched$H3K27ac_clusters) %>% 
    filter(H3K27ac_stage==rna_stage) %>%
    group_by(H3K27ac_clusters) %>% 
    filter(corr==max(corr))


H3K27ac_final <- bind_rows(H3K27ac_matched, H3K27ac_missing) %>% 
    rename('mark_clusters'='H3K27ac_clusters') %>% 
    select(!c(H3K27ac, rna, H3K27me3, H3K4me3, H3K4me3_clusters, H3K27me3_clusters))

H3K27ac_final %>% write_tsv('data/RNA/integration/cluster_matches_H3K27ac.tsv')


#### Join with objects ####

H3K27ac_meta <- H3K27ac_final[match(marks$H3K27ac$clusters, H3K27ac_final$mark_clusters), ] %>% 
    ungroup() %>% 
    select(corr, rna_clusters, celltype_jf) %>% 
    as.data.frame()

rownames(H3K27ac_meta) <- colnames(marks$H3K27ac)

marks$H3K27ac <- AddMetaData(marks$H3K27ac, metadata=H3K27ac_meta)

dim_plot(marks$H3K27ac, group.by='celltype_jf') +
    scale_color_manual(values=celltype_colors)

H3K27me3_meta <- H3K27me3_final[match(marks$H3K27me3$clusters, H3K27me3_final$mark_clusters), ] %>% 
    ungroup() %>% 
    select(corr, rna_clusters, celltype_jf) %>% 
    as.data.frame()

rownames(H3K27me3_meta) <- colnames(marks$H3K27me3)

marks$H3K27me3 <- AddMetaData(marks$H3K27me3, metadata=H3K27me3_meta)

dim_plot(marks$H3K27me3, group.by='celltype_jf') +
    scale_color_manual(values=celltype_colors)


H3K4me3_meta <- H3K4me3_final[match(marks$H3K4me3$clusters, H3K4me3_final$mark_clusters), ] %>% 
    ungroup() %>% 
    select(corr, rna_clusters, celltype_jf) %>% 
    as.data.frame()

rownames(H3K4me3_meta) <- colnames(marks$H3K4me3)

marks$H3K4me3 <- AddMetaData(marks$H3K4me3, metadata=H3K4me3_meta)

dim_plot(marks$H3K4me3, group.by='celltype_jf')


rna_meta <- all_matches_final[match(rna$clusters, all_matches_final$rna_clusters), ] %>% 
    ungroup() %>% 
    select(H3K27ac, H3K27me3, H3K4me3) %>% 
    as.data.frame()
rownames(rna_meta) <- colnames(rna)

rna <- AddMetaData(rna, metadata=rna_meta)


#### Write objects (and one last plot) ####

p1 <- dim_plot(marks$H3K27ac, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors, na.value='grey') + ggtitle('H3K27ac')
p2 <- dim_plot(marks$H3K27me3, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors, na.value='grey') + ggtitle('H3K27me3')
p3 <- dim_plot(marks$H3K4me3, group.by=c('celltype_jf')) +
    scale_color_manual(values=celltype_colors, na.value='grey') + ggtitle('H3K4me3')
p4 <- dim_plot(rna, group.by=c('celltype_jf'), reduction='cssumap') +
    scale_color_manual(values=celltype_colors, na.value='grey') + ggtitle('rna')

(p1 | p2) / (p3 | p4)
ggsave('plots/RNA/integration/highres_combined_integration_umap.png', width=16, height=12, bg='white')



marks <- map(marks, function(srt){
    srt$lineage <- case_when(
        srt$celltype_jf %in% c('ctx_ip', 'ctx_npc', 'ctx_ex') ~ 'ctx',
        srt$celltype_jf %in% c('dien_ex', 'dien_npc') ~ 'dien',
        srt$celltype_jf %in% c('nt_npc', 'mesen_ex', 'rhom_ex') ~ 'nt',
        srt$celltype_jf %in% c('RGC', 'RPC') ~ 'retina',
        T ~ 'other'
    )
    
    srt$state <- case_when(
        srt$celltype_jf %in% c('ctx_ip', 'ctx_npc', 'RPC', 'dien_npc', 'nt_npc') ~ 'npc',
        srt$celltype_jf %in% c('ctx_ex', 'dien_ex', 'rhom_ex', 'mesen_ex', 'RGC') ~ 'neuron',
        srt$celltype_jf %in% c('astrocytes') ~ 'astrocytes',
        T ~ 'other'
    )
    
    srt$age <- case_when(
        srt$age %in% c('65d', '63d') ~ '60d',
        T ~ srt$age
    )
    return(srt)
})

dim_plot(marks$H3K27ac, group.by=c('lineage', 'state'))

rna$lineage <- case_when(
    rna$celltype_jf %in% c('ctx_ip', 'ctx_npc', 'ctx_ex') ~ 'ctx',
    rna$celltype_jf %in% c('dien_ex', 'dien_npc') ~ 'dien',
    rna$celltype_jf %in% c('nt_npc', 'mesen_ex', 'rhom_ex') ~ 'nt',
    rna$celltype_jf %in% c('RGC', 'RPC') ~ 'retina',
    T ~ 'other'
)

rna$state <- case_when(
    rna$celltype_jf %in% c('ctx_ip', 'ctx_npc', 'RPC', 'dien_npc', 'nt_npc') ~ 'npc',
    rna$celltype_jf %in% c('ctx_ex', 'dien_ex', 'rhom_ex', 'mesen_ex', 'RGC') ~ 'neuron',
    rna$celltype_jf %in% c('astrocytes') ~ 'astrocytes',
    T ~ 'other'
)


rna$age <- case_when(
    rna$age %in% c('65d', '63d', '66d') ~ '60d',
    T ~ rna$age
)

dim_plot(marks$H3K27ac, group.by=c('lineage', 'state'))
dim_plot(rna, group.by=c('lineage', 'state'))


marks %>% write_rds('data/all_marks_list_v3.2annot.rds')
rna %>% write_rds('data/RNA/RNA_all_srt_v2.2matched.rds')


#### Write as h5ad ####
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')
library(SeuratDisk)

rna_ <- DietSeurat(rna, dimreducs = c('pca', 'umap', 'css', 'cssumap'))
SaveH5Seurat(rna_, filename='data/RNA/RNA_all_srt_v2.2matched.h5seurat', overwrite=T)
Convert('data/RNA/RNA_all_srt_v2.2matched.h5seurat', dest='h5ad', overwrite=T)





















