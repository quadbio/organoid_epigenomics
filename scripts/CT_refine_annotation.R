source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(destiny)
library(dtw)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')


#### Read data ###
marks <- read_rds('data/all_marks_list_v3.4lines.rds')


#### H3K27ac ####
H3K27ac_lins <- marks$H3K27ac %>% subset(lineage%in%c('ctx', 'dien', 'nt'))

H3K27ac_lins <- H3K27ac_lins %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

H3K27ac_lins <- H3K27ac_lins %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K27ac_lins)
ElbowPlot(H3K27ac_lins, reduction='lsi')

lineage_colors2['nt'] <- lineage_colors2['mesen_rhom']
p1 <- dim_plot(H3K27ac_lins, label=T, group.by=c('celltype_jf', 'lineage')) +
    scale_color_manual(values=lineage_colors2)
p2 <- feature_plot(H3K27ac_lins, features=c('RSPO3', 'ZIC1', 'WLS', 'FOXG1', 'EMX1', 'LHX5', 'HOXB2', 'OTX2', 'STMN2'), order=T, pt.size=1)
p1 / p2 + plot_layout(heights=c(1,2))
ggsave('plots/cutntag/H3K27ac_pre_refine_annot_umap.png', width=10, height=8)

dim_plot(H3K27ac_lins, label=T, group.by=c('clusters'), cols=many_more) +
    no_legend()
ggsave('plots/cutntag/H3K27ac_pre_refine_clusters_umap.png', width=8, height=8)



#### H3K4me3 ####
H3K4me3_lins <- marks$H3K4me3 %>% subset(lineage%in%c('ctx', 'dien', 'nt'))

H3K4me3_lins <- H3K4me3_lins %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

H3K4me3_lins <- H3K4me3_lins %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K4me3_lins)
ElbowPlot(H3K4me3_lins, reduction='lsi')

p1 <- dim_plot(H3K4me3_lins, label=T, group.by=c('celltype_jf', 'lineage')) +
    scale_color_manual(values=lineage_colors2)
p2 <- feature_plot(H3K4me3_lins, features=c('RSPO3', 'ZIC1', 'WLS', 'FOXG1', 'EMX1', 'LHX5', 'HOXB2', 'OTX2', 'STMN2'), order=T, pt.size=1)
p1 / p2 + plot_layout(heights=c(1,2))
ggsave('plots/cutntag/H3K4me3_pre_refine_annot_umap.png', width=10, height=8)

dim_plot(H3K4me3_lins, label=T, group.by=c('clusters')) +
    no_legend()
ggsave('plots/cutntag/H3K4me3_pre_refine_clusters_umap.png', width=8, height=8)


#### H3K27me3 ####
H3K27me3_lins <- marks$H3K27me3 %>% subset(lineage%in%c('ctx', 'dien', 'nt'))

H3K27me3_lins <- H3K27me3_lins %>% 
    FindTopFeatures(min.cutoff='q80') %>% 
    RunSVD() 

H3K27me3_lins <- H3K27me3_lins %>% 
    RunUMAP(dims=2:10, reduction='lsi')

DepthCor(H3K27me3_lins)
ElbowPlot(H3K27me3_lins, reduction='lsi')

p1 <- dim_plot(H3K27me3_lins, label=T, group.by=c('celltype_jf', 'lineage')) +
    scale_color_manual(values=lineage_colors2)
p2 <- feature_plot(H3K27me3_lins, features=c('RSPO3', 'ZIC1', 'WLS', 'FOXG1', 'EMX1', 'LHX5', 'HOXB2', 'OTX2', 'STMN2'), order=T, pt.size=1)
p1 / p2 + plot_layout(heights=c(1,2))
ggsave('plots/cutntag/H3K27me3_pre_refine_annot_umap.png', width=10, height=8)

dim_plot(H3K27me3_lins, label=T, group.by=c('clusters')) +
    no_legend()
ggsave('plots/cutntag/H3K27me3_pre_refine_clusters_umap.png', width=8, height=8)









