source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

setwd('~/projects/cutntag/')


#### Read Data ####
all_marks <- read_rds('data/all_marks_list_v2preproc.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.1annotated.rds')

#### Add sample/age columns ####
all_marks <- all_marks %>%
    map(function(srt){
        sample <- str_replace(colnames(srt), '(.+)_[ATCG]+-1', '\\1')
        srt$orig.ident <- sample
        srt$age <- str_replace(sample, '.+_(.+)', '\\1')
        srt$age <- ifelse(str_detect(sample, 'rep2_35d'), '65d', srt$age)
        srt$age <- ifelse(str_detect(sample, 'rep2_35d'), '65d', srt$age)
        srt$mark <- str_replace(sample, '.+(H3[0-9a-zA-Z]+)_.+', '\\1')
        srt$sample <- paste0(srt$mark, '_', srt$age)
        srt$log_peak_counts <- log10(srt$nCount_peaks)
        srt$log_peaks <- log10(srt$nFeature_peaks)
        return(srt)
    })

p1_H3K27me3 <- all_marks$H3K27me3 %>% dim_plot(group.by=c('age'))
p1_H3K27ac <- all_marks$H3K27ac %>% dim_plot(group.by=c('age'))
p1_H3K4me3 <- all_marks$H3K4me3 %>% dim_plot(group.by=c('age'))

p2_H3K27me3 <- all_marks$H3K27me3 %>% feature_plot(features='log_peak_counts', order=T)
p2_H3K27ac <- all_marks$H3K27ac %>% feature_plot(features='log_peak_counts', order=T)
p2_H3K4me3 <- all_marks$H3K4me3 %>% feature_plot(features='log_peak_counts', order=T)

p3_H3K27me3 <- all_marks$H3K27me3 %>% feature_plot(features=all_markers, order=T)
p3_H3K27ac <- all_marks$H3K27ac %>% feature_plot(features=all_markers, order=T)
p3_H3K4me3 <- all_marks$H3K4me3 %>% feature_plot(features=all_markers, order=T)


(p1_H3K27me3 | p2_H3K27me3) / p3_H3K27me3 + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K27me3_all_init_umap.png', width=10, height=30)
(p1_H3K27ac | p2_H3K27ac) / p3_H3K27ac + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K27ac_all_init_umap.png', width=10, height=30)
(p1_H3K4me3 | p2_H3K4me3) / p3_H3K4me3 + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K4me3_all_init_umap.png', width=10, height=30)


H3K27me3_meta <- all_marks$H3K27me3@meta.data %>% as_tibble()
H3K27ac_meta <- all_marks$H3K27ac@meta.data %>% as_tibble()
H3K4me3_meta <- all_marks$H3K4me3@meta.data %>% as_tibble()

ggplot(H3K27me3_meta, aes(log_peak_counts)) +
    geom_vline(xintercept = log10(100)) +
    geom_histogram()

ggplot(H3K27ac_meta, aes(log_peak_counts)) +
    geom_vline(xintercept = log10(200)) +
    geom_histogram()

ggplot(H3K4me3_meta, aes(log_peak_counts)) +
    geom_vline(xintercept = log10(200)) +
    geom_histogram()


dim_plot(all_marks$H3K27me3)
dim_plot(all_marks$H3K27ac)
dim_plot(all_marks$H3K4me3)


#### Filter again ####
all_marks$H3K27me3 <- all_marks$H3K27me3 %>% subset(nCount_peaks>100 & seurat_clusters!=9 & seurat_clusters!=5)
all_marks$H3K27ac <- all_marks$H3K27ac %>% subset(nCount_peaks>200 & seurat_clusters!=17)
all_marks$H3K4me3 <- all_marks$H3K4me3 %>% subset(nCount_peaks>200 & seurat_clusters!=17)



#### Basic preproc for narrow_peaks ####
all_marks <- all_marks %>% 
    map(function(srt){
        # Dim reduc
        DefaultAssay(srt) <- 'peaks'
        srt <- RunTFIDF(srt)
        srt <- FindTopFeatures(srt, min.cutoff = 'q50')
        srt <- RunSVD(
            object = srt,
            assay = 'peaks',
            reduction.key = 'cLSI_',
            reduction.name = 'clsi'
        )
        
        srt <- RunUMAP(object = srt, reduction = 'clsi', dims = 2:30)
        srt <- FindNeighbors(object = srt, reduction = 'clsi', dims = 2:30)
        srt <- FindClusters(object = srt, verbose = FALSE, algorithm = 1)
        print(DimPlot(object = srt, label = TRUE) + NoLegend())
        ggsave('plots/cutntag/all_cpeaks_filter_umap.png', width=6, height=6)
        
        return(srt)
    })

all_marks %>% write_rds('data/all_marks_list_v3filtered.rds')


p1_H3K27me3 <- all_marks$H3K27me3 %>% dim_plot(group.by=c('age'))
p1_H3K27ac <- all_marks$H3K27ac %>% dim_plot(group.by=c('age'))
p1_H3K4me3 <- all_marks$H3K4me3 %>% dim_plot(group.by=c('age'))

p2_H3K27me3 <- all_marks$H3K27me3 %>% feature_plot(features='log_peak_counts', order=T)
p2_H3K27ac <- all_marks$H3K27ac %>% feature_plot(features='log_peak_counts', order=T)
p2_H3K4me3 <- all_marks$H3K4me3 %>% feature_plot(features='log_peak_counts', order=T)

p3_H3K27me3 <- all_marks$H3K27me3 %>% feature_plot(features=all_markers, order=T)
p3_H3K27ac <- all_marks$H3K27ac %>% feature_plot(features=all_markers, order=T)
p3_H3K4me3 <- all_marks$H3K4me3 %>% feature_plot(features=all_markers, order=T)


(p1_H3K27me3 | p2_H3K27me3) / p3_H3K27me3 + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K27me3_all_filter_umap.png', width=10, height=30)
(p1_H3K27ac | p2_H3K27ac) / p3_H3K27ac + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K27ac_all_filter_umap.png', width=10, height=30)
(p1_H3K4me3 | p2_H3K4me3) / p3_H3K4me3 + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K4me3_all_filter_umap.png', width=10, height=30)



#### Split into three parts: EB, mid (15d), late ####
all_marks <- read_rds('data/all_marks_list_v3filtered.rds')

all_marks$H3K27ac$stage <-  case_when(
    all_marks$H3K27ac$age == 'EB' ~ 'EB',
    all_marks$H3K27ac$age == '15d' ~ 'mid',
    str_detect(all_marks$H3K27ac$age, 'ret-') ~ 'retina',
    str_detect(all_marks$H3K27ac$age, '8mo') ~ 'mo8',
    T ~ 'late'
)

all_marks$H3K27me3$stage <-  case_when(
    all_marks$H3K27me3$age == 'EB' ~ 'EB',
    all_marks$H3K27me3$age == '15d' ~ 'mid',
    str_detect(all_marks$H3K27me3$age, 'ret-') ~ 'retina',
    str_detect(all_marks$H3K27me3$age, '8mo') ~ 'mo8',
    T ~ 'late'
)

all_marks$H3K4me3$stage <-  case_when(
    all_marks$H3K4me3$age == 'EB' ~ 'EB',
    all_marks$H3K4me3$age == '15d' ~ 'mid',
    str_detect(all_marks$H3K4me3$age, 'ret-') ~ 'retina',
    str_detect(all_marks$H3K4me3$age, '8mo') ~ 'mo8',
    T ~ 'late'
)

all_marks$H3K27me3 %>% dim_plot(group.by=c('stage'))
all_marks$H3K27ac %>% dim_plot(group.by=c('stage'))
all_marks$H3K4me3 %>% dim_plot(group.by=c('stage'))

all_marks$H3K27me3 %>% dim_plot(group.by=c('seurat_clusters'))
all_marks$H3K27ac %>% dim_plot(group.by=c('seurat_clusters'))
all_marks$H3K4me3 %>% dim_plot(group.by=c('seurat_clusters'))


#### Remove low count cluster and preproc again ####
# all_marks$H3K4me3 <- all_marks$H3K4me3 %>% subset(seurat_clusters!=14)
all_marks$H3K4me3 <- all_marks$H3K4me3 %>% subset(seurat_clusters!=16)
all_marks$H3K4me3 <- FindTopFeatures(all_marks$H3K4me3, min.cutoff = 'q80')
all_marks$H3K4me3 <- RunSVD(
    object = all_marks$H3K4me3,
    assay = 'peaks',
    reduction.key = 'cLSI_',
    reduction.name = 'clsi'
)

all_marks$H3K4me3 <- RunUMAP(object = all_marks$H3K4me3, reduction = 'clsi', dims = 2:30)
all_marks$H3K4me3 <- FindNeighbors(object = all_marks$H3K4me3, reduction = 'clsi', dims = 2:30)
all_marks$H3K4me3 <- FindClusters(object = all_marks$H3K4me3, verbose = FALSE, algorithm = 1)
print(DimPlot(object = all_marks$H3K4me3, label = TRUE) + NoLegend())
all_marks$H3K4me3 %>% feature_plot(features='log_peak_counts', order=T)

dim_plot(all_marks$H3K4me3, group.by='stage')

#### Split into stages and cluster ####
H3K27ac_stages <- SplitObject(all_marks$H3K27ac, split.by='stage')
H3K27me3_stages <- SplitObject(all_marks$H3K27me3, split.by='stage')
H3K4me3_stages <- SplitObject(all_marks$H3K4me3, split.by='stage')

H3K27ac_stages <- map(H3K27ac_stages, function(srt){
    srt <- FindNeighbors(srt, reduction='clsi', dims=2:30)
    return(srt)
})

H3K27me3_stages <- map(H3K27me3_stages, function(srt){
    srt <- FindNeighbors(srt, reduction='clsi', dims=2:30)
    return(srt)
})

H3K4me3_stages <- map(H3K4me3_stages, function(srt){
    srt <- FindNeighbors(srt, reduction='clsi', dims=2:30)
    return(srt)
})


H3K27ac_stages$EB <- FindClusters(H3K27ac_stages$EB, resolution=2)
H3K27ac_stages$mid <- FindClusters(H3K27ac_stages$mid, resolution=5)
H3K27ac_stages$late <- FindClusters(H3K27ac_stages$late, resolution=10)
H3K27ac_stages$retina <- FindClusters(H3K27ac_stages$retina, resolution=5)
H3K27ac_stages$mo8 <- FindClusters(H3K27ac_stages$mo8, resolution=5)

H3K27ac_stages$EB$clusters <- H3K27ac_stages$EB$peaks_snn_res.2
H3K27ac_stages$mid$clusters <- H3K27ac_stages$mid$peaks_snn_res.5
H3K27ac_stages$late$clusters <- H3K27ac_stages$late$peaks_snn_res.10
H3K27ac_stages$retina$clusters <- H3K27ac_stages$retina$peaks_snn_res.5
H3K27ac_stages$mo8$clusters <- H3K27ac_stages$mo8$peaks_snn_res.5




H3K27me3_stages$EB <- FindClusters(H3K27me3_stages$EB, resolution=2)
H3K27me3_stages$mid <- FindClusters(H3K27me3_stages$mid, resolution=5)
H3K27me3_stages$late <- FindClusters(H3K27me3_stages$late, resolution=10)
H3K27me3_stages$retina <- FindClusters(H3K27me3_stages$retina, resolution=5)
H3K27me3_stages$mo8 <- FindClusters(H3K27me3_stages$mo8, resolution=5)

H3K27me3_stages$EB$clusters <- H3K27me3_stages$EB$peaks_snn_res.2
H3K27me3_stages$mid$clusters <- H3K27me3_stages$mid$peaks_snn_res.5
H3K27me3_stages$late$clusters <- H3K27me3_stages$late$peaks_snn_res.10
H3K27me3_stages$retina$clusters <- H3K27me3_stages$retina$peaks_snn_res.5
H3K27me3_stages$mo8$clusters <- H3K27me3_stages$mo8$peaks_snn_res.5




H3K4me3_stages$EB <- FindClusters(H3K4me3_stages$EB, resolution=2)
H3K4me3_stages$mid <- FindClusters(H3K4me3_stages$mid, resolution=5)
H3K4me3_stages$late <- FindClusters(H3K4me3_stages$late, resolution=10)
H3K4me3_stages$retina <- FindClusters(H3K4me3_stages$retina, resolution=5)
H3K4me3_stages$mo8 <- FindClusters(H3K4me3_stages$mo8, resolution=5)

H3K4me3_stages$EB$clusters <- H3K4me3_stages$EB$peaks_snn_res.2
H3K4me3_stages$mid$clusters <- H3K4me3_stages$mid$peaks_snn_res.5
H3K4me3_stages$late$clusters <- H3K4me3_stages$late$peaks_snn_res.10
H3K4me3_stages$retina$clusters <- H3K4me3_stages$retina$peaks_snn_res.5
H3K4me3_stages$mo8$clusters <- H3K4me3_stages$mo8$peaks_snn_res.5



#### Add combined clusters to object ####
H3K27ac_clusters <- purrr::reduce(map(H3K27ac_stages, function(x) paste0(x$stage, '_', x$clusters)), c)
names(H3K27ac_clusters) <- purrr::reduce(map(H3K27ac_stages, colnames), c)
all_marks$H3K27ac <- all_marks$H3K27ac %>% AddMetaData(H3K27ac_clusters, col.name = 'clusters')

H3K27me3_clusters <- purrr::reduce(map(H3K27me3_stages, function(x) paste0(x$stage, '_', x$clusters)), c)
names(H3K27me3_clusters) <- purrr::reduce(map(H3K27me3_stages, colnames), c)
all_marks$H3K27me3 <- all_marks$H3K27me3 %>% AddMetaData(H3K27me3_clusters, col.name = 'clusters')

H3K4me3_clusters <- purrr::reduce(map(H3K4me3_stages, function(x) paste0(x$stage, '_', x$clusters)), c)
names(H3K4me3_clusters) <- purrr::reduce(map(H3K4me3_stages, colnames), c)
all_marks$H3K4me3 <- all_marks$H3K4me3 %>% AddMetaData(H3K4me3_clusters, col.name = 'clusters')


#### Plot clusters ####
p1 <- dim_plot(all_marks$H3K27ac, group.by=c('clusters'), reduction='umap', cols=many_more) & 
    no_legend()
p2 <- dim_plot(all_marks$H3K27me3, group.by=c('clusters'), reduction='umap', cols=many_more) & 
    no_legend()
p3 <- dim_plot(all_marks$H3K4me3, group.by=c('clusters'), reduction='umap', cols=many_more) & 
    no_legend()
p1 / p2 / p3 

ggsave('plots/cnt_RNA_highres_stage_clustering.png', width=4, height=8)

all_clusters_stage <- bind_rows(
    enframe(all_marks$H3K27ac$clusters, 'cell', 'cluster'),
    enframe(all_marks$H3K27me3$clusters, 'cell', 'cluster'),
    enframe(all_marks$H3K4me3$clusters, 'cell', 'cluster')
)

all_clusters_stage %>% write_tsv('data/cnt_all_clusters_stage.tsv')
all_marks %>% write_rds('data/all_marks_list_v3.1clustered.rds')




