source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

summarize <- dplyr::summarize

setwd('~/projects/cutntag/')

marks <- read_rds('data/all_marks_list_v3.3motifs.rds')
rna_full <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')

meta <- rna_full@meta.data %>%
    as_tibble(rownames='cell') %>%
    inner_join(as_tibble(rna_full[['cssumap']]@cell.embeddings, rownames='cell'))

cluster_meta <- meta %>%
    group_by(clusters) %>%
    dplyr::summarize(CSSUMAP_1=mean(CSSUMAP_1), CSSUMAP_2=mean(CSSUMAP_2), celltype_jf=celltype_jf[1])


ggplot(meta, aes(CSSUMAP_1, CSSUMAP_2, color=celltype_jf, fill=celltype_jf)) +
    geom_point(size=0.1, alpha=0.1) +
    geom_point(data=cluster_meta, size=4, shape=21, color='darkgrey') +
    scale_color_manual(values=pantone_celltype) +
    scale_fill_manual(values=pantone_celltype) +
    theme_void() +
    no_legend()

ggsave('plots/RNA/cluster_celltype_umap.png', width=7, height=5)
ggsave('plots/RNA/cluster_celltype_umap.pdf', width=7, height=5)

ggplot(meta, aes(CSSUMAP_1, CSSUMAP_2, color=celltype_jf, fill=celltype_jf)) +
    geom_point(size=0.1, alpha=0.1) +
    # geom_point(data=cluster_meta, size=4, shape=21, color='darkgrey') +
    scale_color_manual(values=pantone_celltype) +
    scale_fill_manual(values=pantone_celltype) +
    theme_void() +
    no_legend()

ggsave('plots/RNA/celltype_umap.png', width=7, height=5)
ggsave('plots/RNA/celltype_umap.pdf', width=7, height=5)



#### Subset astro and recompute CSS ####

noastro_cells <- colnames(rna_full)[!rna_full$celltype_jf%in%c('astrocytes', 'OPC', 'other', 'choroid_plexus', 'non_nect')]
rna_noastro <- rna_full %>% subset(cells=noastro_cells)

rna_EB <- subset(rna_noastro, orig.ident=='210705_5_FZ_scRNA_HWW4CB_EB')
rna_ret <- subset(rna_noastro, orig.ident%in%c('22011001_2_FZ_scRNA_ret_6w_B7', '22010501_1_FZ_scRNA_ret_12w_B7'))
rna_8mo <- subset(rna_noastro, orig.ident%in%c('211207_9_FZ_scRNA_HWW4CB_8mo'))
rna_other <- subset(rna_noastro, 
                    orig.ident!= '210705_5_FZ_scRNA_HWW4CB_EB' & 
                        orig.ident!= '22011001_2_FZ_scRNA_ret_6w_B7' & 
                        orig.ident!= '22010501_1_FZ_scRNA_ret_12w_B7' &
                        orig.ident!= '211207_9_FZ_scRNA_HWW4CB_8mo'
)

rna_other <- FindVariableFeatures(rna_other)
rna_other <- CellCycleScoring(rna_other, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
VariableFeatures(rna_other) <- VariableFeatures(rna_other) %>% setdiff(cc_genes_all)
rna_other <- ScaleData(rna_other, vars.to.regress = c('S.Score', 'G2M.Score'))
rna_other <- rna_other %>% RunPCA() 
rna_other[['css']] <- NULL

n <- 10
r <- 0.6
d <- 0.4
s <- 1
rna_other <- cluster_sim_spectrum(
    rna_other,
    label_tag='orig.ident',
    use_dr='pca', dims_use = 1:n,
    reduction.name = 'css', reduction.key = 'CSS_'
)

rna_other <- RunUMAP(
    object = rna_other,
    spread = s,
    min.dist = d,
    reduction = 'css', 
    dims = 1:ncol(rna_other[['css']]),
    reduction.name = 'cssumap',
    reduction.key = 'CSSUMAP_'
)

dim_plot(rna_other, group.by=c('orig.ident', 'celltype_jf'), reduction='cssumap')
feature_plot(rna_other, features=all_markers, order=T, reduction='cssumap')

css_model <- rna_other[['css']]@misc$model
EB_project <- css_project(rna_EB, model=css_model)[['css_proj']]@cell.embeddings
ret_project <- css_project(rna_ret, model=css_model)[['css_proj']]@cell.embeddings
m8_project <- css_project(rna_8mo, model=css_model)[['css_proj']]@cell.embeddings

colnames(EB_project) <- colnames(css_model$sim2profiles)
colnames(ret_project) <- colnames(css_model$sim2profiles)
colnames(m8_project) <- colnames(css_model$sim2profiles)

# combined_css <- rbind(EB_project, ret_project, css_model$sim2profiles)
combined_css <- rbind(EB_project, ret_project, m8_project, css_model$sim2profiles)

css_ <- CreateAssayObject(t(combined_css))
rna_noastro[['css_']] <- css_
rna_noastro <- ScaleData(rna_noastro, vars.to.regress = c('S.Score', 'G2M.Score'), assay='css_')
css_cc <- t(rna_noastro[['css_']]@scale.data)
rna_noastro[['css_']] <- NULL
rna_noastro[['css']] <- CreateDimReducObject(as.matrix(css_cc), key = 'CSS_')


d <- 0.2
s <- 0.8
rna_noastro <- RunUMAP(
    object = rna_noastro,
    spread = s,
    min.dist = d,
    reduction = 'css', 
    dims = 1:ncol(rna_noastro[['css']]),
    reduction.name = 'cssumap',
    reduction.key = 'CSSUMAP_'
)

dim_plot(rna_noastro, group.by=c('orig.ident', 'celltype_jf'), reduction='cssumap', pt.size=0.01)

library(SeuratDisk)
rna_noastro %>% write_rds('data/RNA/RNA_all_srt_v3noastro.rds')
rna_noastro <- read_rds('data/RNA/RNA_all_srt_v3noastro.rds')

rna_noastro_ <- DietSeurat(rna_noastro, scale.data=F, dimreducs=c('css', 'cssumap', 'pca'))
VariableFeatures(rna_noastro_) <- NULL
SeuratDisk::SaveH5Seurat(rna_noastro_, 'data/RNA/RNA_all_srt_v3noastro.h5seurat', overwrite=T)
SeuratDisk::Convert('data/RNA/RNA_all_srt_v3noastro.h5seurat', dest='h5ad', overwrite=T)




#### Get cellrank results ####
cellrank_meta <- read_tsv('data/RNA/cellrank/RNA_noastro_cellrank_probs.tsv')
colnames(cellrank_meta)[1] <- 'cell'
cellrank_meta_df <- column_to_rownames(cellrank_meta, 'cell')

rna <- rna_full[, cellrank_meta$cell]
rna <- AddMetaData(rna, cellrank_meta_df)

stage_clusters <- read_tsv('data/all_clusters_stage.tsv')
rna <- AddMetaData(rna, column_to_rownames(stage_clusters, 'cell'))

cr_trans_mat <- readMM('data/RNA/cellrank/RNA_noastro_velo_cr_transition.mtx')
colnames(cr_trans_mat) <- colnames(rna)
rownames(cr_trans_mat) <- colnames(rna)

cr_cluster_mat <- Pando::aggregate_matrix(cr_trans_mat, groups=as.character(rna$clusters))
cr_cluster_mat <- Pando::aggregate_matrix(t(cr_cluster_mat), groups=as.character(rna$clusters))

paga_conn_mat <- readMM('data/RNA/cellrank/RNA_noastro_velo_paga_connectivities.mtx')
colnames(paga_conn_mat) <- unique(rna$clusters)
rownames(paga_conn_mat) <- unique(rna$clusters)


#### DR in cellrank space ####
cr_space <- cellrank_meta %>%
    select(cell, to_ctx_ranks, to_mesen_ranks, to_dien_ranks, to_rhom_ranks, to_retina_ranks, pseudotime_ranks) %>%
    column_to_rownames('cell') %>% as.matrix()

cr_pca <- pca(cr_space, n = 2, to_df = T) %>%
    inner_join(cellrank_meta)

p1 <- ggplot(cr_pca, aes(PC1, PC2, color=to_ctx)) +
    geom_point(size=0.5, alpha=0.8) +
    scale_color_gradientn(colors=gyylgnbu()) +
    theme_void()
p2 <- ggplot(cr_pca, aes(PC1, PC2, color=to_mesen)) +
    geom_point(size=0.5, alpha=0.8) +
    scale_color_gradientn(colors=gyylgnbu()) +
    theme_void()
p1 | p2



#### Reduce to clusters ####
cluster_counts <- as.numeric(table(rna$clusters))
cr_hc_space <- Pando::aggregate_matrix(column_to_rownames(cr_pca, 'cell'), groups = as.character(rna$clusters))
cr_hc_df <- cr_hc_space %>% as_tibble(rownames='clusters')

cluster_counts <- table(rna$clusters, rna$age)[cr_hc_df$clusters, ]

age_clusters <- map(as.character(sort(unique(rna$age))), function(x){
    rownames(cluster_counts)[cluster_counts[, x]>0]
})

ggplot(cr_hc_df, aes(x=.panel_x, y=.panel_y, color=pseudotime_ranks)) +
    geom_point(alpha=0.8, size=3) +
    geom_autodensity() +
    scale_color_viridis() +
    facet_matrix(vars(to_ctx, to_dien, to_rhom, to_mesen, to_retina, pseudotime_ranks), layer.diag=2) +
    theme_article()
ggsave('plots/RNA/cellrank/cellrank_lineage_probs_facet.png', width=12, height=12)




#### Plot stuff ####
p1 <- feature_plot(rna, features=c('pseudotime_ranks'), reduction='cssumap') &
    scale_color_viridis(option = 'magma', direction = -1)
p2 <- dim_plot(rna, group.by=c('celltype_jf'), label=T, reduction='cssumap') +
    scale_color_manual(values=celltype_colors)
names(age_colors2) <- names(age_colors2) %>% str_replace('-', '_')
p3 <- dim_plot(rna, group.by=c('age'), label=T, reduction='cssumap') +
    scale_color_manual(values=age_colors2)
p1 | p2 | p3
ggsave('plots/RNA/rna_pt_celltype_age_umap.png', width=12, height=5)

umap_coords <- rna[['cssumap']]@cell.embeddings %>%
    as_tibble(rownames='cell')
cluster_coords <- aggregate_matrix(rna[['cssumap']]@cell.embeddings, as.character(rna$clusters)) %>%
    as_tibble(rownames='clusters')

meta <- rna@meta.data %>% as_tibble(rownames='cell') %>%
    inner_join(umap_coords)

cluster_meta <- meta %>%
    group_by(clusters) %>%
    summarize(
        CSSUMAP_1 = mean(CSSUMAP_1),
        CSSUMAP_2 = mean(CSSUMAP_2),
        mode_age = mode(age),
        celltype = celltype_jf[1]
    ) %>%
    inner_join(cr_hc_df) %>%
    mutate(lineage=case_when(
        celltype %in% c('mesen_ex', 'rhom_ex', 'nt_npc') ~ 'mesen_rhom',
        celltype %in% c('ctx_npc', 'ctx_ip', 'ctx_ex') ~ 'ctx',
        celltype %in% c('RPC', 'RGC') ~ 'retina',
        celltype %in% c('dien_npc', 'dien_ex') ~ 'dien',
        celltype %in% c('nect') ~ 'nect',
        celltype %in% c('psc') ~ 'EB',
        celltype %in% c('astrocytes') ~ 'astro',
        celltype %in% c('non_nect') ~ 'nn',
        celltype %in% c('choroid_plexus') ~ 'chp',
        celltype %in% c('OPC') ~ 'OPC'
    ))


ggplot(meta, aes(CSSUMAP_1, CSSUMAP_2)) +
    geom_point(color='grey') +
    geom_point(data=cluster_meta, mapping=aes(fill=celltype), size=6, shape=21) +
    scale_fill_manual(values=celltype_colors) +
    theme_void()

ggplot(meta, aes(CSSUMAP_1, CSSUMAP_2)) +
    geom_point(color='grey') +
    geom_point(data=cluster_meta, mapping=aes(fill=lineage), size=6, shape=21) +
    scale_fill_manual(values=lineage_colors2) +
    theme_void()

cluster_meta %>% write_tsv('data/RNA/rna_highres_cluster_meta.tsv')

clusters_use <- cluster_meta %>%
    filter(celltype!='other') %>%
    pull(clusters)


#### Plot paga graph ####
paga_graph <- as_tbl_graph(graph_from_adjacency_matrix(paga_conn_mat[clusters_use,clusters_use], weighted = T)) %>%
    inner_join(cluster_meta, by=c('name'='clusters')) %E>%
    filter(weight==1)

p1 <- ggraph(paga_graph, x=CSSUMAP_1, y=CSSUMAP_2) +
    geom_edge_link() +
    geom_node_point(mapping=aes(fill=celltype), size=6, shape=21) +
    scale_fill_manual(values=celltype_colors)

cellrank_graph <- as_tbl_graph(graph_from_adjacency_matrix(cr_cluster_mat[clusters_use,clusters_use], weighted = T)) %>%
    inner_join(cluster_meta, by=c('name'='clusters')) %E>%
    filter(weight>0.00001)

p2 <- ggraph(cellrank_graph, x=CSSUMAP_1, y=CSSUMAP_2) +
    geom_edge_link() +
    geom_node_point(mapping=aes(fill=celltype), size=6, shape=21) +
    scale_fill_manual(values=celltype_colors)

(p1 | p2) & no_legend()
ggsave('plots/RNA/cellrank/cellrank_paga_cluster_graph_umap.png', width=12, height=5)



#### Get NN graph ####
cr_hc_df <- cr_hc_df %>%
    filter(clusters%in%clusters_use)
cr_hc_df$clusters <- as.character(cr_hc_df$clusters)
pt_order_idx <- order(cr_hc_df$pseudotime_ranks)
pt_order <- as.character(cr_hc_df$clusters[pt_order_idx])
cr_hc_mat <- column_to_rownames(cr_hc_df, 'clusters')[, 5:(ncol(cr_hc_df)-1)]
cr_hc_mat <- cr_hc_mat[pt_order, ]

# Try finding highest transition neighbors instead
k <- 11
# Get 5 highest transitions
cr_cluster_mat_fwd <- as.matrix(cr_cluster_mat[pt_order, pt_order])
cr_cluster_mat_fwd[upper.tri(cr_cluster_mat_fwd)] <- 0
cr_nn_mat <- apply(as.matrix(cr_cluster_mat_fwd), 1, function(x){
    top_idx <- order(x, decreasing = T)[1:k]
    x[-top_idx] <- 0
    return(x)
})

nn_mat <- cr_nn_mat[pt_order, pt_order]
diag(nn_mat) <- 0

# nn_mat <- cr_cluster_mat
nn_mat[nn_mat>0] <- 1

pheatmap::pheatmap(nn_mat, cluster_cols = F, cluster_rows = F)
nn_graph <- as_tbl_graph(graph_from_adjacency_matrix(nn_mat)) %>%
    inner_join(cluster_meta, by=c('name'='clusters'))


#### Plot assigned labels ####
nn_assign_pruned <- nn_graph %E>%
    mutate(from_celltype=.N()$celltype[from], to_celltype=.N()$celltype[to]) %>%
    mutate(from_lineage=.N()$lineage[from], to_lineage=.N()$lineage[to]) %>%
    mutate(prune=case_when(
        (from_lineage == 'ctx') & (to_lineage != 'ctx') ~ T,
        (from_lineage == 'dien') & (to_lineage != 'dien') ~ T,
        (from_lineage == 'mesen_rhom') & (to_lineage != 'mesen_rhom') ~ T,
        (from_lineage == 'retina') & (to_lineage != 'retina') ~ T,
        (from_lineage == 'nn') & (to_lineage != 'nn') ~ T,
        # (to_lineage == 'astro') & (from_lineage == 'nect') ~ T,
        (from_lineage == 'EB') & (!to_lineage %in% c('nn', 'nect', 'EB')) ~ T,
        (from_celltype %in% c('ctx_ex', 'mesen_ex', 'rhom_ex', 'dien_ex', 'RGC')) & from_celltype!=to_celltype ~ T,
        # (from_celltype=='OPC') & (to_celltype %in% c('ctx_ex', 'mesen_ex', 'rhom_ex', 'dien_ex', 'RGC')) ~ T,
        # (from_celltype=='OPC') & (to_celltype %in% c('ctx_ex', 'mesen_ex', 'rhom_ex', 'dien_ex', 'RGC')) ~ T,
        # (from_celltype=='ctx_npc') & (to_celltype == 'ctx_ex') ~ T,
        # (!from_celltype%in%c('dien_npc', 'choroid_plexus')) & (to_celltype == 'choroid_plexus') ~ T,
        # (!to_celltype%in%c('dien_npc', 'choroid_plexus')) & (from_celltype == 'choroid_plexus') ~ T,
        T ~ F
    )) %>%
    filter(!prune)

p1 <- ggraph(nn_assign_pruned, x=CSSUMAP_1, y=CSSUMAP_2) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=celltype), size=4) +
    scale_color_manual(values=pantone_celltype) +
    theme_void() +
    no_legend()

set.seed(55)
p2 <- ggraph(nn_assign_pruned, layout='fr') +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=celltype), size=4) +
    scale_color_manual(values=pantone_celltype) +
    theme_void()
p2

p1 | p2
ggsave('plots/RNA/cellrank/cellrank_lineage_assign_graph.png', width=15, height=7)


nn_assign_pruned %>% write_rds('data/noastro_RNA_cluster_graph.rds')

fr_df <- p2$data[,c('name','x','y')] %>%
    as_tibble()
colnames(fr_df)[2:3] <- c('FR1', 'FR2')
fr_df$FR2[fr_df$FR2<(-6)] <- fr_df$FR2[fr_df$FR2<(-6)] + 3

nn_graph_fr <- nn_assign_pruned %N>%
    inner_join(fr_df)

p2 <- ggraph(nn_graph_fr, x=FR1, y=FR2) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=celltype), size=4) +
    scale_color_manual(values=pantone_celltype) +
    theme_void()
p2

cluster_meta_out <- cluster_meta %>% inner_join(fr_df, by=c('clusters'='name'))
cluster_meta_out %>% write_tsv('data/noastro_RNA_cluster_meta.tsv')



ggplot(cluster_meta_out, aes(FR1, FR2, color=celltype)) +
    geom_point() +
    theme_minimal()


set.seed(55)
ggraph(nn_assign_pruned, layout='fr') +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(color=celltype), size=4) +
    scale_color_manual(values=pantone_celltype) +
    theme_void()

ggsave('plots/RNA/cellrank/cellrank_lineage_assign_graph.png', width=7, height=7)
ggsave('plots/RNA/cellrank/cellrank_lineage_assign_graph.pdf', width=7, height=7)


ggraph(nn_assign_pruned, x=CSSUMAP_1, y=CSSUMAP_2) +
    geom_edge_link(alpha=0.1) +
    geom_node_point(aes(fill=celltype), size=4, shape=21, color='darkgrey') +
    scale_fill_manual(values=pantone_celltype) +
    theme_void() +
    no_legend()

ggsave('plots/RNA/cellrank/cellrank_lineage_assign_umap.png', width=7, height=5)
ggsave('plots/RNA/cellrank/cellrank_lineage_assign_umap.pdf', width=7, height=5)

















