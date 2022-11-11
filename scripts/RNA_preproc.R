source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

setwd('~/projects/cutntag/')

select <- dplyr::select

cc_genes_all <- purrr::reduce(cc.genes.updated.2019, union)

#### Read data ####
counts_path <- '/local2/USERS/jfleck/projects/cutntag/counts/RNA'
counts_list <- list.files(counts_path, full.names = T)
counts_names <- list.files(counts_path)
names(counts_list) <- counts_names

matrix_list <- map(counts_list, function(x){
    path=paste0(x, '/raw_feature_bc_matrix.h5')
    print(path)
    return(Read10X_h5(path))
})

# Some light filtering
matrix_list <- map(matrix_list, function(x){
    x[, colSums(x)>100]
})

srt_list <- map(matrix_list, function(x){
    srt <- CreateSeuratObject(
        counts=x,
        assay='RNA'
    )
    return(srt)
})


#### QC and proper filtering ####
srt_list <- map(srt_list, basic_qc)

meta <- map_dfr(srt_list, ~as_tibble(.x@meta.data, rownames='cell'), .id='sample')

p1 <- ggplot(meta, aes(nCount_RNA)) + scale_x_log10()
p2 <- ggplot(meta, aes(nFeature_RNA)) + scale_x_log10()
p3 <- ggplot(meta, aes(percent_mito))

p1 / p2 / p3 & geom_histogram() & facet_grid(~sample)
ggsave('plots/RNA/all_basic_qc.png', width=20, height=10)


ggplot(meta, aes(nCount_RNA, nFeature_RNA)) +
    geom_point(size=0.5) +
    facet_wrap(~sample)

ggplot(meta, aes(nCount_RNA, percent_mito)) +
    geom_point(size=0.5) +
    scale_x_log10() +
    facet_wrap(~sample)


srt_list %>% write_rds('data/RNA/RNA_all_list_v0raw.rds')


srt_list_fltr <- map(srt_list, 
    ~subset(
        .x, 
        nCount_RNA>2000 & 
            nCount_RNA<1.5e5 & 
            nFeature_RNA>1000 & 
            percent_mito<0.2 
))

meta <- map_dfr(srt_list_fltr, ~as_tibble(.x@meta.data, rownames='cell'), .id='sample')
p1 <- ggplot(meta, aes(nCount_RNA)) 
p2 <- ggplot(meta, aes(nFeature_RNA)) 
p3 <- ggplot(meta, aes(percent_mito))

p1 / p2 / p3 & geom_histogram() & facet_grid(~sample, scales='free')
ggsave('plots/RNA/all_filtered_qc.png', width=20, height=10)


srt_list_fltr %>% write_rds('data/RNA/RNA_all_list_v1filtered.rds')


#### Merge and preproc ####
srt_list_fltr <- map(set_names(names(srt_list_fltr)), function(x) {srt_list_fltr[[x]]$orig.ident <- x; return(srt_list_fltr[[x]])})
rna <- merge(srt_list_fltr[[1]], srt_list_fltr[2:length(srt_list_fltr)], add.cell.ids=names(srt_list_fltr))

rna <- rna %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

ElbowPlot(rna)

rna <- rna %>% 
    FindNeighbors(dims=1:20) %>% 
    FindClusters() %>% 
    RunUMAP(dims=1:20)

p1 <- dim_plot(rna, group.by='orig.ident')
p2 <- dim_plot(rna)
p3 <- feature_plot(rna, features=all_markers, order=T)

(p1 | p2) / p3 + plot_layout(heights=c(1,5))
ggsave('plots/RNA/all_joined_init_umap.png', width=15, height=30, bg='white')

rna %>% write_rds('data/RNA/RNA_all_srt_v1preproc.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v1preproc.rds')


#### Glycolysis filtering ####
rna_data <- GetAssayData(rna, assay='RNA', slot='data')
nowa_fit <- read_rds('~/resources/models/glycolysis/crop3_glyc.rds')

nowa_prediction <- nowakowski_predict(nowa_fit, t(rna_data))
nowa_meta <- nowa_prediction %>%
    group_by(cell) %>%
    filter(p_pred==max(p_pred)) %>%
    dplyr::select('nowakowski_prediction'=cell_type, cell) %>%
    column_to_rownames('cell')

rna <- AddMetaData(rna, nowa_meta)
rna$glyc <- rna$nowakowski_prediction=='Glyc'
dim_plot(rna, group.by='glyc', reduction='umap')

rna <- subset(rna, glyc==FALSE)




#### Integrate ####
library(SeuratWrappers)
library(harmony)

#### MNN ####
rna_split <- SplitObject(rna, 'orig.ident')
rna_mnn <- RunFastMNN(
    rna_split,
    assay='RNA',
    features=2000
)

rna_mnn <- RunUMAP(
    object = rna_mnn,
    reduction = 'mnn', dims = 1:ncol(Reductions(rna_mnn, 'mnn')),
    reduction.name = 'mnnumap',
    reduction.key = 'MNNUMAP_'
)

dim_plot(rna_mnn, group.by=c('group', 'line'), reduction='mnnumap')


#### Harmony ####
rna <- RunHarmony(rna, group.by.vars = 'orig.ident')
rna <- RunUMAP(
    object = rna,
    reduction = 'harmony', dims = 1:30,
    reduction.name = 'humap',
    reduction.key = 'HUMAP_'
)
dim_plot(rna, group.by=c('orig.ident'), reduction='humap')
feature_plot(rna, features=c('FOXG1', 'ZIC2', 'WLS', 'RSPO3', 'SIX6', 'KRT8', 'MKI67'), order=T, reduction='humap')


#### Seurat ####
anchors <- FindIntegrationAnchors(rna_split)
rna_int <- IntegrateData(anchors)
rna_int <- rna_int %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    RunUMAP(dims=1:20)

rna_int <- rna_int %>% 
    FindNeighbors(dims=1:20) %>% 
    FindClusters()

dim_plot(rna_int, group.by=c('group', 'line'), reduction='umap')
feature_plot(rna_int, features=c('FOXG1', 'ZIC2', 'WLS', 'RSPO3', 'SIX6', 'KRT8', 'MKI67'), order=T)


rna_int %>% write_rds('data/more_lines/lines12_RNATAC16_26_integrated_srt.rds')





#### Integrate timepoints/samples ####
rna <- CellCycleScoring(rna, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)

rna_EB <- subset(rna, orig.ident=='210705_5_FZ_scRNA_HWW4CB_EB')
rna_ret <- subset(rna, orig.ident%in%c('22011001_2_FZ_scRNA_ret_6w_B7', '22010501_1_FZ_scRNA_ret_12w_B7'))
rna_8mo <- subset(rna, orig.ident%in%c('211207_9_FZ_scRNA_HWW4CB_8mo'))
rna_other <- subset(rna, 
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

dim_plot(rna_other, group.by='orig.ident', reduction='cssumap')
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
rna[['css_']] <- css_
rna <- ScaleData(rna, vars.to.regress = c('S.Score', 'G2M.Score'), assay='css_')
css_cc <- t(rna[['css_']]@scale.data)
rna[['css_']] <- NULL
rna[['css']] <- CreateDimReducObject(as.matrix(css_cc), key = 'CSS_')


d <- 0.2
s <- 0.8
rna <- RunUMAP(
    object = rna,
    spread = s,
    min.dist = d,
    reduction = 'css', 
    dims = 1:ncol(rna[['css']]),
    reduction.name = 'cssumap',
    reduction.key = 'CSSUMAP_'
)

dim_plot(rna, group.by='orig.ident', reduction='cssumap', pt.size=0.01)

rna <- FindNeighbors(
    rna,
    reduction = 'css', 
    dims = 1:ncol(rna[['css']])
)

rna <- FindClusters(rna, resolution=2)


#### Remove low qual cluster and redo ####
rna <- rna %>% subset(seurat_clusters!=48 & seurat_clusters!=33 & seurat_clusters!=60)
rna <- RunUMAP(
    object = rna,
    spread = s,
    min.dist = d,
    reduction = 'css', 
    dims = 1:ncol(rna[['css']]),
    reduction.name = 'cssumap',
    reduction.key = 'CSSUMAP_'
)

dim_plot(rna, group.by='seurat_clusters', reduction='cssumap', pt.size=0.01)

rna <- FindNeighbors(
    rna,
    reduction = 'css', 
    dims = 1:ncol(rna[['css']])
)

rna <- FindClusters(rna, resolution=2)

p1 <- dim_plot(rna, group.by='orig.ident', reduction='cssumap')
p2 <- dim_plot(rna, group.by='seurat_clusters', reduction='cssumap', label=T)
p3 <- feature_plot(rna, features=all_markers, order=T, reduction='cssumap')

(p1 | p2) / p3 + plot_layout(heights=c(1,5))
ggsave('plots/RNA/all_joined_css_final_umap.png', width=15, height=30, bg='white')


more_markers <- c('APOE', 'S100B', 'OLIG2', 'HOPX', 'ALDH1L1', 'SATB2', 'BCL11A', 'PLP1', 'PCDH15')
feature_plot(rna, features=more_markers, order=T, reduction='cssumap')
ggsave('plots/RNA/all_joined_css_more_markers_umap.png', width=10, height=9, bg='white')

more_markers <- c('STMN2', 'NEUROG1', 'WIF1', 'HES5', 'CNTN2')
feature_plot(rna, features=more_markers, order=T, reduction='cssumap')
ggsave('plots/RNA/all_joined_css_more_markers2_umap.png', width=7, height=6, bg='white')


rna %>% write_rds('data/RNA/RNA_all_srt_v2integrated.rds')


#### Annotate stages ####
rna <- read_rds('data/RNA/RNA_all_srt_v2integrated.rds')

sample <- str_replace(colnames(rna), '(.+)_[ATCG]+-1', '\\1')
rna$orig.ident <- sample
rna$age <- str_replace(sample, '.+_(.+)', '\\1')
rna$age <- ifelse(str_detect(sample, 'rep2_35d'), '65d', rna$age)
rna$age <- ifelse(str_detect(sample, 'ret_6w'), 'ret_6w', rna$age)
rna$age <- ifelse(str_detect(sample, 'ret_12w'), 'ret_12w', rna$age)
rna$mark <- str_replace(sample, '.+(H3[0-9a-zA-Z]+)_.+', '\\1')

rna$stage <-  case_when(
    rna$age == 'EB' ~ 'EB',
    rna$age == '15d' ~ 'mid',
    str_detect(rna$age, 'ret_') ~ 'retina',
    str_detect(rna$age, '8mo') ~ 'mo8',
    T ~ 'late'
)


#### Cluster stages ####
rna_stages <- SplitObject(rna, split.by='stage')

rna_stages <- map(rna_stages, function(srt){
    srt <- FindNeighbors(srt, reduction='css', dims=1:ncol(srt[['css']]))
    return(srt)
})

rna_stages$EB <- FindClusters(rna_stages$EB, resolution=1)
rna_stages$mid <- FindClusters(rna_stages$mid, resolution=2)
rna_stages$late <- FindClusters(rna_stages$late, resolution=5)
rna_stages$retina <- FindClusters(rna_stages$retina, resolution=2)
rna_stages$mo8 <- FindClusters(rna_stages$mo8, resolution=2)

rna_stages$EB$clusters <- rna_stages$EB$RNA_snn_res.1
rna_stages$mid$clusters <- rna_stages$mid$RNA_snn_res.2
rna_stages$late$clusters <- rna_stages$late$RNA_snn_res.5
rna_stages$retina$clusters <- rna_stages$retina$RNA_snn_res.2
rna_stages$mo8$clusters <- rna_stages$mo8$RNA_snn_res.2

rna_clusters <- purrr::reduce(map(rna_stages, function(x) paste0(x$stage, '_', x$clusters)), c)
names(rna_clusters) <- purrr::reduce(map(rna_stages, colnames), c)
rna <- rna %>% AddMetaData(rna_clusters, col.name = 'clusters')

dim_plot(rna, group.by='clusters', reduction='cssumap', label=T, cols=many_more) +
    no_legend()
ggsave('plots/RNA/all_joined_css_highres_clusters_umap.png', width=12, height=12, bg='white')


#### Annotate clusters ####
rna$celltype_jf <- case_when(
    rna$clusters%in%c('retina_17','late_25','late_23','mid_6','mid_2','mid_7','mid_1','mid_8','mid_4','mid_12','mid_3','mid_9','mid_5','mid_14','mid_10') ~ 'nect',
    rna$clusters%in%c('retina_7','retina_6','retina_1','retina_3','retina_10','retina_5','retina_2','retina_9','retina_0','late_50') ~ 'RPC',
    rna$clusters%in%c('retina_14','retina_8','retina_11','retina_4') ~ 'RGC',
    rna$clusters%in%c('late_11', 'late_8') ~ 'rhom_ex',
    rna$clusters%in%c('late_32', 'late_5') ~ 'mesen_ex',
    rna$clusters%in%c('late_48','late_13','late_6','late_30','late_4','late_16','late_45','late_0','late_38','late_39','late_41','late_15','late_35','late_16','late_1','late_19') ~ 'nt_npc',
    rna$clusters%in%c('late_26', 'late_47', 'mo8_18', 'retina_21') ~ 'dien_ex',
    rna$clusters%in%c('late_9', 'late_17', 'mo8_5') ~ 'ctx_ex',
    rna$clusters%in%c('mo8_6', 'mo8_10') ~ 'OPC',
    rna$clusters%in%c('late_44', 'late_46', 'late_9') ~ 'ctx_ip',
    rna$clusters%in%c('mo8_14','late_21','late_40','late_12','late_7','late_22','late_2','late_37') ~ 'ctx_npc',
    rna$clusters%in%c('retina_12','retina_22','retina_13','late_20','late_36','late_43','late_3','late_24','late_31','late_18','mo8_15') ~ 'dien_npc',
    rna$clusters%in%c('late_29','late_49','late_28','late_42','mo8_0','retina_18') ~ 'choroid_plexus',
    rna$clusters%in%c('late_34','late_10','retina_19','mo8_9','mo8_13','mo8_7','mo8_3','mo8_1','mo8_12','mo8_2','mo8_8') ~ 'astrocytes',
    rna$clusters%in%c('mid_11','mid_13','mo8_19','late_51','mid_0') ~ 'non_nect',
    str_detect(rna$clusters, 'EB') ~ 'psc',
    T  ~ 'other'
)

dim_plot(rna, group.by=c('celltype_jf'), reduction='cssumap', label=T) +
    scale_color_manual(values=c(celltype_colors, retina_colors2)) + no_legend()
ggsave('plots/RNA/all_joined_css_celltype_umap.png', width=6, height=6, bg='white') 

rna %>% write_rds('data/RNA/RNA_all_srt_v2.1annotated.rds')










