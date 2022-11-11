source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

setwd('~/projects/cutntag/')


#### Read stuff ####
H3K27ac <- read_rds('data/H3K27ac/H3K27ac_srt_v2macs.rds')
H3K4me3 <- read_rds('data/H3K4me3/H3K4me3_srt_v2macs.rds')
H3K27me3 <- read_rds('data/H3K27me3/H3K27me3_srt_v2macs.rds')

srt_list <- list(
    H3K27me3 = H3K27me3,
    H3K4me3 = H3K4me3,
    H3K27ac = H3K27ac
)

#### Basic preproc for peaks ####
srt_list <- srt_list %>% 
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
        ggsave('plots/cutntag/all_cpeaks_init_umap.png', width=6, height=6)
        
        # compute gene activities
        gene_activities <- GeneActivity(srt, assay = 'peaks')
        # add the gene activity matrix to the Seurat object as a new assay
        srt[['cRNA']] <- CreateAssayObject(counts = gene_activities)
        srt <- NormalizeData(
            object = srt,
            assay = 'cRNA',
            normalization.method = 'LogNormalize',
            scale.factor = median(srt$nCount_cRNA)
        )
        return(srt)
    })

srt_list %>% write_rds('data/all_marks_list_v2preproc.rds')
# srt_list <- read_rds('data/all_marks_list_v2preproc.rds')

#### Basic preproc for narrow_peaks ####
srt_list <- srt_list %>% 
    map(function(srt){
        # Dim reduc
        DefaultAssay(srt) <- 'narrow_peaks'
        srt <- RunTFIDF(srt)
        srt <- FindTopFeatures(srt, min.cutoff = 'q50')
        srt <- RunSVD(
            object = srt,
            assay = 'narrow_peaks',
            reduction.key = 'mLSI_',
            reduction.name = 'mlsi'
        )
        
        srt <- RunUMAP(object = srt, reduction = 'mlsi', dims = 2:30)
        srt <- FindNeighbors(object = srt, reduction = 'mlsi', dims = 2:30)
        srt <- FindClusters(object = srt, verbose = FALSE, algorithm = 1)
        print(DimPlot(object = srt, label = TRUE) + NoLegend())
        ggsave('plots/cutntag/all_mpeaks_init_umap.png', width=6, height=6)
        
        # compute gene activities
        gene_activities <- GeneActivity(srt, assay = 'narrow_peaks')
        # add the gene activity matrix to the Seurat object as a new assay
        srt[['mRNA']] <- CreateAssayObject(counts = gene_activities)
        srt <- NormalizeData(
            object = srt,
            assay = 'mRNA',
            normalization.method = 'LogNormalize',
            scale.factor = median(srt$nCount_mRNA)
        )
        return(srt)
    })

srt_list %>% write_rds('data/all_marks_list_v2preproc.rds')


#### Add sample/age columns ####
srt_list_annot <- srt_list %>%
    map(function(srt){
        sample <- str_replace(colnames(srt), '(.+)_[ATCG]+-1', '\\1')
        srt$orig.ident <- sample
        srt$age <- str_replace(sample, '.+_(.+)', '\\1')
        srt$age <- ifelse(str_detect(sample, 'rep2_35d'), '65d', srt$age)
        srt$age <- ifelse(str_detect(sample, 'rep2_35d'), '65d', srt$age)
        srt$mark <- str_replace(sample, '.+(H3[0-9a-zA-Z]+)_.+', '\\1')
        srt$sample <- paste0(srt$mark, '_', srt$age)
        srt$log_peak_counts <- log(srt$nCount_peaks)
        srt$log_peaks <- log(srt$nFeature_peaks)
        return(srt)
    })

p1_H3K27me3 <- srt_list_annot$H3K27me3 %>% dim_plot(group.by=c('age'))
p1_H3K27ac <- srt_list_annot$H3K27ac %>% dim_plot(group.by=c('age'))
p1_H3K4me3 <- srt_list_annot$H3K4me3 %>% dim_plot(group.by=c('age'))

p2_H3K27me3 <- srt_list_annot$H3K27me3 %>% feature_plot(features='log_peak_counts', order=T)
p2_H3K27ac <- srt_list_annot$H3K27ac %>% feature_plot(features='log_peak_counts', order=T)
p2_H3K4me3 <- srt_list_annot$H3K4me3 %>% feature_plot(features='log_peak_counts', order=T)

p3_H3K27me3 <- srt_list_annot$H3K27me3 %>% feature_plot(features=all_markers, order=T)
p3_H3K27ac <- srt_list_annot$H3K27ac %>% feature_plot(features=all_markers, order=T)
p3_H3K4me3 <- srt_list_annot$H3K4me3 %>% feature_plot(features=all_markers, order=T)


(p1_H3K27me3 | p2_H3K27me3) / p3_H3K27me3 + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K27me3_all_init_umap.png', width=10, height=30)
(p1_H3K27ac | p2_H3K27ac) / p3_H3K27ac + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K27ac_all_init_umap.png', width=10, height=30)
(p1_H3K4me3 | p2_H3K4me3) / p3_H3K4me3 + plot_layout(heights=c(1,5))
ggsave('plots/cutntag/H3K4me3_all_init_umap.png', width=10, height=30)
















