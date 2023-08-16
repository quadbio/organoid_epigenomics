source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')

options(future.globals.maxSize=999999999999999)
library(future)
plan('multicore', workers = 32)

#### Read data ####
rna <- read_rds('data05/RNA/RNA_all_srt_v2.3lines.rds')

# # Select based on RNA DE
# Idents(rna) <- 'celltype_jf'
# celltype_de_rna <- FindAllMarkers(
#     rna, logfc.threshold = 0.25, test.use = 'LR', latent.vars = 'nCount_RNA'
# )
# 
# celltype_de_rna %>% write_tsv('data_/results/diff_expression/RNA_DE_celltype.tsv')


# rna$lineage_coarse <- case_when(
#     rna$lineage != 'other' ~ rna$lineage,
#     T ~ rna$celltype_jf
# )
# 
# Idents(rna) <- 'lineage_coarse'
# lineage_de_rna <- FindAllMarkers(
#     rna, logfc.threshold = 0.25, test.use = 'LR', latent.vars = 'nCount_RNA'
# )
# 
# lineage_de_rna %>% write_tsv('data_/results/diff_expression/RNA_DE_lineage.tsv')


rna_ <- subset(rna, age!='ret_12w' & age!='ret_6w' & lineage=='ctx')

Idents(rna_) <- 'age'
lineage_de_rna <- FindAllMarkers(
    rna_, logfc.threshold=0, test.use='LR', latent.vars='nCount_RNA'
)

lineage_de_rna %>% write_tsv('data_/results/diff_expression/RNA_DE_ctx_age.tsv')
