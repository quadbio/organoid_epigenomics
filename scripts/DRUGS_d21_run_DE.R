source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/perturbator/enrichment.R')

library(Pando)
library(SeuratWrappers)
library(harmony)
library(future)

plan('multiprocess', workers=36)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')

#### Read stuff ####

drugs_d21_H3K27me3 <- read_rds('data/drugs/drugs_d21_A395_v1_v1.2lines_srt.rds')
drugs_d21_H3K27ac <- read_rds('data/drugs/drugs_d21_A485_v1_v1.2lines_srt.rds')


#### Do DE for H3K27me3 inhib ####
Idents(drugs_d21_H3K27me3) <- 'inhibitor_target'
de_cluster_inhib <- map_dfr(set_names(unique(drugs_d21_H3K27me3$RNA_snn_res.0.1)), function(clust){
    print(clust)
    srt <- subset(drugs_d21_H3K27me3, RNA_snn_res.0.1==clust)
    de_df <- try(FindMarkers(
        srt,
        ident.1='H3K27me3', ident.2='DMSO',
        min.pct=0.05, logfc.threshold=0.1,
        test.use='LR',
        latent.vars=c('nCount_RNA')
    ))
    if (any(class(de_df)=='try-error')){
        return(tibble())
    }
    return(as_tibble(de_df, rownames='feature'))
}, .id='celltype')

de_cluster_inhib %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A395_cluster_inhib_de.tsv')

de_df <- FindMarkers(
    drugs_d21_H3K27me3,
    ident.1='H3K27me3', ident.2='DMSO',
    min.pct=0.05, logfc.threshold=0.1,
    test.use='LR',
    latent.vars=c('nCount_RNA', 'RNA_snn_res.0.1')
)

as_tibble(de_df, rownames='feature') %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A395_global_inhib_de.tsv')


Idents(drugs_d21_H3K27me3) <- 'RNA_snn_res.0.1'
de_cluster <- map_dfr(set_names(unique(drugs_d21_H3K27me3$RNA_snn_res.0.1)), function(clust){
    print(clust)
    de_df <- try(FindMarkers(
        drugs_d21_H3K27me3,
        ident.1=clust,
        min.pct=0.05, logfc.threshold=0.1,
        test.use='LR',
        latent.vars=c('nCount_RNA')
    ))
    if (any(class(de_df)=='try-error')){
        return(tibble())
    }
    return(as_tibble(de_df, rownames='feature'))
}, .id='celltype')

de_cluster %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A395_cluster_de.tsv')



#### Test concentrations separately ####
Idents(drugs_d21_H3K27me3) <- 'inhib_annotation'

all_inhib_conds <- drugs_d21_H3K27me3$inhib_annotation %>% unique %>%
    {.[!str_detect(., 'DMSO')]} %>% set_names()
all_dmso_conds <- drugs_d21_H3K27me3$inhib_annotation %>% unique %>% {.[str_detect(., 'DMSO')]}

conc_de_df <- map_dfr(all_inhib_conds, function(conc){
    de_df <- FindMarkers(
        drugs_d21_H3K27me3,
        ident.1=conc, ident.2=all_dmso_conds,
        min.pct=0.05, logfc.threshold=0.1,
        test.use='LR',
        latent.vars=c('nCount_RNA', 'RNA_snn_res.0.1')
    )
    return(as_tibble(de_df, rownames='feature'))
}, .id='inhib_annotation')

conc_de_df %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A395_concs_global_inhib_de.tsv')




#### Do DE for H3K27ac inhib ####
Idents(drugs_d21_H3K27ac) <- 'inhibitor_target'
de_cluster_inhib <- map_dfr(set_names(unique(drugs_d21_H3K27ac$RNA_snn_res.0.1)), function(clust){
    print(clust)
    srt <- subset(drugs_d21_H3K27ac, RNA_snn_res.0.1==clust)
    de_df <- try(FindMarkers(
        srt,
        ident.1='H3K27ac', ident.2='DMSO',
        min.pct=0.05, logfc.threshold=0.1,
        test.use='LR',
        latent.vars=c('nCount_RNA')
    ))
    if (any(class(de_df)=='try-error')){
        return(tibble())
    }
    return(as_tibble(de_df, rownames='feature'))
}, .id='celltype')

de_cluster_inhib %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A485_cluster_inhib_de.tsv')


de_df <- FindMarkers(
    drugs_d21_H3K27ac,
    ident.1='H3K27ac', ident.2='DMSO',
    min.pct=0.05, logfc.threshold=0.1,
    test.use='LR',
    latent.vars=c('nCount_RNA', 'RNA_snn_res.0.1')
)

as_tibble(de_df, rownames='feature') %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A485_global_inhib_de.tsv')


de_df <- FindMarkers(
    drugs_d21_H3K27ac,
    ident.1='H3K27ac', ident.2='DMSO',
    min.pct=0.05, logfc.threshold=0.1,
    test.use='LR',
    latent.vars=c('nCount_RNA')
)

as_tibble(de_df, rownames='feature') %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A485_global2_inhib_de.tsv')


Idents(drugs_d21_H3K27ac) <- 'RNA_snn_res.0.1'
de_cluster <- map_dfr(set_names(unique(drugs_d21_H3K27ac$RNA_snn_res.0.1)), function(clust){
    print(clust)
    de_df <- try(FindMarkers(
        drugs_d21_H3K27ac,
        ident.1=clust,
        min.pct=0.05, logfc.threshold=0.1,
        test.use='LR',
        latent.vars=c('nCount_RNA')
    ))
    if (any(class(de_df)=='try-error')){
        return(tibble())
    }
    return(as_tibble(de_df, rownames='feature'))
}, .id='celltype')

de_cluster %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A485_cluster_de.tsv')



#### Test concentrations separately ####
Idents(drugs_d21_H3K27ac) <- 'inhib_annotation'

all_inhib_conds <- drugs_d21_H3K27ac$inhib_annotation %>% unique %>%
    {.[!str_detect(., 'DMSO')]} %>% set_names()
all_dmso_conds <- drugs_d21_H3K27ac$inhib_annotation %>% unique %>% {.[str_detect(., 'DMSO')]}

conc_de_df <- map_dfr(all_inhib_conds, function(conc){
    de_df <- FindMarkers(
        drugs_d21_H3K27ac,
        ident.1=conc, ident.2=all_dmso_conds,
        min.pct=0.05, logfc.threshold=0.1,
        test.use='LR',
        latent.vars=c('nCount_RNA', 'RNA_snn_res.0.1')
    )
    return(as_tibble(de_df, rownames='feature'))
}, .id='inhib_annotation')

conc_de_df %>% write_tsv('data/drugs/results/diff_expression/drugs_d21_A485_concs_global_inhib_de.tsv')









