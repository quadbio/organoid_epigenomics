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
cluster_graph <- read_rds('data/RNA/all_RNA_cluster_graph.rds')


#### Subset astrocytes in RNA ####

# rna_ac <- rna %>% subset(celltype_jf%in%c('astrocytes', 'nt_npc', 'dien_npc', 'ctx_npc') & age%in%c('8mo', '4mo', '60d'))
rna_ac <- rna %>% subset(celltype_jf%in%c('astrocytes') & age%in%c('8mo'))

rna_ac <- rna_ac %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

ElbowPlot(rna_ac)

rna_ac <- rna_ac %>% 
    FindNeighbors(dims=1:20) %>% 
    RunUMAP(dims=1:20)

rna_ac <- rna_ac %>% 
    FindClusters(resolution=0.6)

dim_plot(rna_ac)

p1 <- dim_plot(rna_ac)
p2 <- feature_plot(rna_ac, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'LHX2', 'APOE'), order=T)

p1 / p2

rna_ac$ac_region <- case_when(
    rna_ac$RNA_snn_res.0.6%in%c(6) ~ 'dien',
    rna_ac$RNA_snn_res.0.6%in%c(2,4,10) ~ 'mesen',
    rna_ac$RNA_snn_res.0.6%in%c(0) ~ 'rhom',
    rna_ac$RNA_snn_res.0.6%in%c(7) ~ 'other',
    T ~ 'telen'
)

p1 <- dim_plot(rna_ac, group.by=c('ac_region', 'RNA_snn_res.0.6'))
p2 <- feature_plot(rna_ac, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'WLS', 'GBX2', 'NHLH2', 'APOE', 'S100B'), order=T)

p1 / p2 + plot_layout(heights=c(1,3))
ggsave('plots/RNA/astrocytes_subcluster_umap.png', width=8, height=10, bg='white')


#### Subset astrocytes + progenitors in RNA ####
rna_ac_npc <- rna %>% subset(celltype_jf%in%c('astrocytes', 'nt_npc', 'dien_npc', 'ctx_npc') & age%in%c('8mo', '4mo', '60d'))

rna_ac_npc <- rna_ac_npc %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA() 

ElbowPlot(rna_ac_npc)

rna_ac_npc <- rna_ac_npc %>% 
    FindNeighbors(dims=1:20) %>% 
    RunUMAP(dims=1:20)

rna_ac_npc <- rna_ac_npc %>% 
    FindClusters(resolution=0.6)

dim_plot(rna_ac_npc)

ac_meta <- rna_ac@meta.data['ac_region']

rna_ac_npc <- AddMetaData(rna_ac_npc, ac_meta)

p1 <- dim_plot(rna_ac_npc, group.by=c('ac_region', 'age', 'celltype_jf'), label=T) +
    scale_color_manual(values=pantone_celltype) &
    no_legend()
p2 <- feature_plot(rna_ac_npc, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'WLS', 'GBX2', 'NHLH2', 'APOE', 'S100B'), order=T)

p1 / p2 + plot_layout(heights=c(1,3))
ggsave('plots/RNA/astrocytes_npcs_subcluster_umap.png', width=8, height=10, bg='white')


dim_plot(rna_ac_npc, group.by=c('ac_region', 'clusters', 'celltype_jf'), label=T) +
    scale_color_manual(values=pantone_celltype) &
    no_legend()


rna_ac %>% write_rds('data/trajectories/astrocytes/RNA_astrocytes_srt.rds')
rna_ac_npc %>% write_rds('data/trajectories/astrocytes/RNA_astrocytes_npcs_srt.rds')



#### Subset astrocytes in H3K27ac ####
H3K27ac_ac <- marks$H3K27ac %>% subset(celltype_jf%in%c('astrocytes') & age%in%c('8mo'))

H3K27ac_ac <- H3K27ac_ac %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

ElbowPlot(H3K27ac_ac)

H3K27ac_ac <- H3K27ac_ac %>% 
    FindNeighbors(dims=2:10, reduction='lsi') %>% 
    RunUMAP(dims=2:10, reduction='lsi')

feature_plot(H3K27ac_ac, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'WLS', 'GBX2', 'NHLH2', 'APOE', 'S100B'), order=T)
ggsave('plots/cutntag/H3K27ac_astrocytes_feature_umap.png', width=6, height=6, bg='white')


H3K27ac_ac <- FindClusters(H3K27ac_ac)

H3K27ac_ac$ac_region <- case_when(
    H3K27ac_ac$seurat_clusters %in% c(2,4,7) ~ 'telen',
    T ~ 'nt'
)

dim_plot(H3K27ac_ac, label=T, group.by=c('seurat_clusters', 'ac_region'))

H3K27ac_ac %>% write_rds('data/trajectories/astrocytes/H3K27ac_astrocytes_srt.rds')
H3K27ac_ac <- read_rds('data/trajectories/astrocytes/H3K27ac_astrocytes_srt.rds')



#### Subset astrocytes in H3K4me3 ####
H3K4me3_ac <- marks$H3K4me3 %>% subset(celltype_jf%in%c('astrocytes') & age%in%c('8mo'))

H3K4me3_ac <- H3K4me3_ac %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

ElbowPlot(H3K4me3_ac)

H3K4me3_ac <- H3K4me3_ac %>% 
    FindNeighbors(dims=2:10, reduction='lsi') %>% 
    RunUMAP(dims=2:10, reduction='lsi')

feature_plot(H3K4me3_ac, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'WLS', 'GBX2', 'NHLH2', 'APOE', 'S100B'), order=T)
ggsave('plots/cutntag/H3K4me3_astrocytes_feature_umap.png', width=6, height=6, bg='white')

H3K4me3_ac <- FindClusters(H3K4me3_ac)

H3K4me3_ac$ac_region <- case_when(
    H3K4me3_ac$seurat_clusters %in% c(2) ~ 'telen',
    T ~ 'nt'
)

dim_plot(H3K4me3_ac, label=T, group.by=c('seurat_clusters', 'ac_region'))

H3K4me3_ac %>% write_rds('data/trajectories/astrocytes/H3K4me3_astrocytes_srt.rds')
H3K4me3_ac <- read_rds('data/trajectories/astrocytes/H3K4me3_astrocytes_srt.rds')




#### Subset astrocytes in H3K27me3 ####
H3K27me3_ac <- marks$H3K27me3 %>% subset(celltype_jf%in%c('astrocytes') & age%in%c('8mo'))

H3K27me3_ac <- H3K27me3_ac %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

ElbowPlot(H3K27me3_ac)

H3K27me3_ac <- H3K27me3_ac %>% 
    FindNeighbors(dims=2:10, reduction='lsi') %>% 
    RunUMAP(dims=2:10, reduction='lsi')

feature_plot(H3K27me3_ac, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'WLS', 'GBX2', 'NHLH2', 'APOE', 'S100B'), order=T)
ggsave('plots/cutntag/H3K27me3_astrocytes_feature_umap.png', width=6, height=6, bg='white')

H3K27me3_ac <- FindClusters(H3K27me3_ac)

H3K27me3_ac$ac_region <- case_when(
    H3K27me3_ac$seurat_clusters %in% c(2,3,4) ~ 'telen',
    T ~ 'nt'
)

dim_plot(H3K27me3_ac, label=T, group.by=c('seurat_clusters', 'ac_region'))

H3K27me3_ac %>% write_rds('data/trajectories/astrocytes/H3K27me3_astrocytes_srt.rds')
H3K27me3_ac <- read_rds('data/trajectories/astrocytes/H3K27me3_astrocytes_srt.rds')




cluster_srt@active.assay <- 'H3K27ac_RNA'
p1 <- dim_plot(cluster_srt, group.by=c('age', 'celltype_jf'), pt.size=4)
p2 <- feature_plot(cluster_srt, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'WLS', 'GBX2', 'NHLH2', 'APOE', 'S100B'), order=T, pt.size=4)

p1 / p2 + plot_layout(heights=c(1,3))


cluster_srt@active.assay <- 'H3K4me3_RNA'
p1 <- dim_plot(cluster_srt, group.by=c('age', 'celltype_jf'), pt.size=4)
p2 <- feature_plot(cluster_srt, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'WLS', 'GBX2', 'NHLH2', 'APOE', 'S100B'), order=T, pt.size=4)

p1 / p2 + plot_layout(heights=c(1,3))


cluster_srt@active.assay <- 'H3K27me3_RNA'
p1 <- dim_plot(cluster_srt, group.by=c('age', 'celltype_jf'), pt.size=4)
p2 <- feature_plot(cluster_srt, features=c('FOXG1', 'RSPO3', 'HOXB2', 'OTX2', 'WLS', 'GBX2', 'NHLH2', 'APOE', 'S100B'), order=T, pt.size=4)

p1 / p2 + plot_layout(heights=c(1,3))







