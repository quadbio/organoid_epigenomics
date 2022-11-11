source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(destiny)

setwd('~/projects/cutntag/')



#### Read data ###
marks <- read_rds('data/all_marks_list_v3.2annot.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')

marks <- map(marks, function(x){
    x$lineage2 <- case_when(
        x$celltype_jf=='astrocytes' ~ 'astrocytes',
        T ~ x$lineage
    )
    return(x)
})

rna$lineage2 <- case_when(
    rna$celltype_jf=='astrocytes' ~ 'astrocytes',
    T ~ rna$lineage
)

marks_lin <- marks %>% map(SplitObject, split.by='lineage2')
rna_lin <- rna %>% SplitObject(split.by='lineage2')

marks_lin <- purrr::transpose(marks_lin)
marks_ctx <- marks_lin$ctx

H3K27ac_ctx <- marks_ctx$H3K27ac

H3K27ac_ctx <- H3K27ac_ctx %>% 
    FindTopFeatures(min.cutoff='q5') %>% 
    RunSVD() 

H3K27ac_ctx <- H3K27ac_ctx %>% 
    RunUMAP(dims=2:30, reduction='lsi')

pca <- H3K27ac_ctx[['lsi']]@cell.embeddings[,2:30]
diffmap <- DiffusionMap(pca, verbose=T)
H3K27ac_ctx[['diffmap']] <- CreateDimReducObject(
    diffmap@eigenvectors, key='DC_', assay='peaks'
)

H3K27ac_ctx %>% dim_plot(group.by='sample')
H3K27ac_ctx %>% dim_plot(group.by='sample', reduction='lsi', dims=c(1,4))
H3K27ac_ctx %>% dim_plot(group.by='state', reduction='diffmap', dims=c(1,1))




rna_ctx <- rna_lin$ctx

rna_ctx <- rna_ctx %>% 
    FindTopFeatures(min.cutoff='q5') %>% 
    RunSVD() 

rna_ctx <- rna_ctx %>% 
    RunUMAP(dims=2:30, reduction='lsi')

pca <- rna_ctx[['lsi']]@cell.embeddings[,2:30]
diffmap <- DiffusionMap(pca, verbose=T)
rna_ctx[['diffmap']] <- CreateDimReducObject(
    diffmap@eigenvectors, key='DC_', assay='RNA'
)

rna_ctx %>% dim_plot(group.by='state')
rna_ctx %>% dim_plot(group.by='state', reduction='lsi', dims=c(1,4))
rna_ctx %>% dim_plot(group.by='state', reduction='diffmap', dims=c(1,1))


rna_lins_use <- map(rna_lins_use, function(x){
    pca <- x[['lsi']]@cell.embeddings[,2:30]
    diffmap <- DiffusionMap(pca, verbose=T)
    x[['diffmap']] <- CreateDimReducObject(
        diffmap@eigenvectors, key='DC_'
    )
    return(x)
})

    




