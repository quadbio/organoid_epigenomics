source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')


#### Read data ####
marks <- read_rds('data05/all_marks_list_v3.4lines.rds')
rna <- read_rds('data05/RNA/RNA_all_srt_v2.3lines.rds')

rna %>% dim_plot(reduction='cssumap', group.by=c('line'))
rna %>% feature_plot(reduction='cssumap', features=c('AQP4', 'HOPX', 'GFAP', 'VIM', 'NES', 'SOX2', 'GLI3', 'NR2F1', 'ALDH1L1'), order=T)


#### Unintegrated umap split by timepoint ####
rna_split <- rna %>% SplitObject('age')
rna_split <- map(rna_split, function(srt){
    srt %>% RunUMAP(reduction='pca', dims=1:20) %>% return()
})

plots <- map(names(rna_split), function(x){
    dim_plot(rna_split[[x]], group.by='celltype_jf') + ggtitle(x) %>% return()
})
plots %>% names
wrap_plots(plots) & scale_color_manual(values=pantone_celltype)
ggsave('plots/paper/sfig1/unintegrated_umap_timepoint.png', width=15, height=12, bg='white')
ggsave('plots/paper/sfig1/unintegrated_umap_timepoint.pdf', width=15, height=12, bg='white')

rna_split[[1]]$celltype_jf %>% unique %>% setdiff(names(pantone_celltype))


#### Subset astrocytes and npcs ####
rna_npc <- subset(rna, state%in%c('astrocytes', 'npc') & celltype_jf!='RPC')
rna_npc <- rna_npc %>% RunUMAP(reduction='css', dims=1:ncol(rna_npc[['css']]))    

pantone_celltype <- c('#8BAED3', '#F2A1BD', '#C76BC1', '#D25B78', '#a9dfbf',  '#949398', '#A1C8C9', 
                      '#8A5796', '#5961A1', '#00569B', '#049DB1', '#85A0AB','#E2C1BE', '#46996C', 
                      '#D2927D', '#C6AAAF', '#ECC270', '#F5DF4D')
celltypes <- c('nt_npc', 'ctx_npc', 'oRG', 'ctx_ex', 'astrocytes',  'mesenchyme', 
               'choroid_plexus', 'nect', 'mesen_ex', 'rhom_ex', 'dien_npc', 'dien_ex','ctx_ip', 
               'OPC', 'non_nect', 'psc', 'RGC', 'RPC')
names(pantone_celltype) <- celltypes

p1 <- rna_npc %>% dim_plot(group.by=c('celltype_jf', 'clusters', 'age'), label=T) & no_legend() 
p2 <- rna_npc %>% feature_plot(features=c('AQP4', 'HOPX', 'VIM', 'NES', 'SOX2', 'GLI3', 'NR2F1', 'ALDH1L1', 'TNC', 'BCAN', 'GFAP'), order=T)

p1 / p2 + plot_layout(heights=c(1,3))

rna_npc$celltype_jf[rna_npc$clusters%in%c('late_2', 'late_37')] <- 'oRG'

rna_npc %>% dim_plot(group.by=c('celltype_jf')) & scale_color_manual(values=pantone_celltype)
ggsave('plots/paper/sfig_astrocytes/npc_astro_annot.png', width=5, height=5, bg='white')
ggsave('plots/paper/sfig_astrocytes/npc_astro_annot.pdf', width=5, height=5, bg='white')

# New colors with darker grey
gyylorrd2 <- function(bias=1){
    return(colorRampPalette(c('#cccccc', brewer.pal(n=9, name='YlOrRd')[3:9]), bias=bias)(100))
}

rna_npc %>% feature_plot(features=c('AQP4', 'HOPX', 'VIM', 'NES', 'SOX2', 'GLI3', 'NR2F1', 'ALDH1L1', 'TNC', 'BCAN', 'GFAP'), order=T) &
    scale_color_gradientn(colors=gyylorrd2())

ggsave('plots/paper/sfig_astrocytes/npc_astro_features.png', width=12, height=10, bg='white')
ggsave('plots/paper/sfig_astrocytes/npc_astro_features.pdf', width=12, height=10, bg='white')

rna_npc %>% write_rds('data_/RNA/RNA_npc_astro_srt_v2.3lines.rds')

#### Old astrocyte annotation ####
rna_astro <- read_rds('data_/trajectories/astrocytes/RNA_astrocytes_srt.rds')
rna_npc <- read_rds('data_/RNA/RNA_npc_astro_srt_v2.3lines.rds')

rna_npc <- AddMetaData(rna_npc, rna_astro$ac_region, col.name='ac_region')
rna_npc %>% dim_plot(group.by=c('ac_region')) 

rna_npc %>% write_rds('data_/RNA/RNA_npc_astro_srt_v2.3lines.rds')


#### Subset neurons to check inhub/excit ####
rna_neuron <- subset(rna, state%in%c('neuron'))
rna_neuron <- rna_neuron %>% RunUMAP(reduction='css', dims=1:ncol(rna_neuron[['css']]))    

p1 <- rna_neuron %>% dim_plot(group.by=c('celltype_jf', 'clusters', 'age'), label=T) & no_legend()
p2 <- rna_neuron %>% feature_plot(features=c('SLC17A7', 'SLC17A6', 'SLC32A1', 'SATB2', 'GAD2', 'LHX5', 'LHX1', 'ISL1', 'BCL11B'), order=T)

p1 / p2 + plot_layout(heights=c(1,3))

rna_neuron$neuron_type <- case_when(
    rna_neuron$clusters %in% c('retina_11', 'retina_8', 'late_47', 'late_26', 'late_46', 'retina_21', 'late_5') ~ 'excitatory',
    rna_neuron$celltype_jf %in% c('ctx_ex') ~ 'excitatory',
    T ~ 'inhibitory'
)

p1 <- rna_neuron %>% dim_plot(group.by=c('neuron_type')) & scale_color_manual(values=many[c(2,1)])
p2 <- rna_neuron %>% dim_plot(group.by=c('celltype_jf')) & scale_color_manual(values=pantone_celltype)

p1 | p2

ggsave('plots/paper/sfig_astrocytes/neuron_annot.png', width=8, height=4, bg='white')
ggsave('plots/paper/sfig_astrocytes/neuron_annot.pdf', width=8, height=4, bg='white')

rna_neuron %>% write_rds('data_/RNA/RNA_neuron_srt_v2.3lines.rds')

#### QC plots also split by line ####
rna_meta <- rna@meta.data %>% 
    as_tibble(rownames='cell') 

rna_summary <- rna_meta %>% 
    dplyr::group_by(line) %>% 
    dplyr::summarize(
        median_feats=median(nFeature_RNA),
        median_counts=median(nCount_RNA)
    )

ggplot(rna_meta, aes(nFeature_RNA)) +
    geom_histogram(color='black', fill='grey') +
    geom_vline(data=rna_summary, mapping=aes(xintercept=median_feats),col='red',size=1)+
    facet_grid(line~.)

ggplot(rna_meta, aes(nCount_RNA)) +
    geom_histogram(color='black', fill='grey') +
    geom_vline(data=rna_summary, mapping=aes(xintercept=median_counts),col='red',size=1)+
    facet_grid(line~.)



#### C&T QC plots also split by line ####
ct_meta <- map_dfr(marks, function(x){
    x@meta.data %>% as_tibble(rownames='cell') 
}, .id='mark')

ct_summary <- ct_meta %>% 
    dplyr::group_by(line, mark) %>% 
    dplyr::summarize(
        median_feats=median(nFeature_peaks),
        median_counts=median(nCount_peaks)
    )

ggplot(ct_meta, aes(nCount_peaks)) +
    geom_histogram(color='black', fill='grey', size=0.2) +
    geom_vline(data=ct_summary, mapping=aes(xintercept=median_counts),col='red',size=0.3) +
    scale_x_continuous(limits=c(0,6000)) +
    facet_grid(line~mark, scales='free') +
    article_text() +
    labs(x='# fragments in peaks')

ggsave('plots/paper/sfig1/QC_per_line_hist.pdf', width=10, height=8, bg='white', units='cm')















