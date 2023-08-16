library(tidyverse)
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(ChIPseeker)
library(hdf5r)
library(anndata)


rename <- dplyr::rename
filter <- dplyr::filter
select <- dplyr::select
summarize <- dplyr::summarize
dist <- stats::dist

setwd('~/projects/cutntag/')

#### FUNC ####
braun_ad <- hdf5r::h5file('data05/braun_2022_fetal_brain_v2_common_hv2k.h5ad')
braun_mat <- sparseMatrix(
    i = braun_ad[['X/indices']][]+1,
    x = braun_ad[['X/data']][],
    p = braun_ad[['X/indptr']][]
)
rownames(braun_mat) <- braun_ad[['obs/CellID']][]
colnames(braun_mat) <- braun_ad[['var/_index']][]

braun_srt <- CreateSeuratObject(t(braun_mat))

CellClass <- braun_ad[['obs/CellClass/codes']][] + 1
CellClass_idx <- braun_ad[['obs/CellClass/categories']][]
braun_srt$CellClass <- CellClass_idx[CellClass]

Region <- braun_ad[['obs/Region/codes']][] + 1
Region_idx <- braun_ad[['obs/Region/categories']][]
braun_srt$Region <- Region_idx[Region]

Subregion <- braun_ad[['obs/Subregion/codes']][] + 1
Subregion_idx <- braun_ad[['obs/Subregion/categories']][]
braun_srt$Subregion <- Subregion_idx[Subregion]

Donor <- braun_ad[['obs/Donor/codes']][] + 1
Donor_idx <- braun_ad[['obs/Donor/categories']][]
braun_srt$Donor <- Donor_idx[Donor]

Age <- braun_ad[['obs/Age/']][]
braun_srt$Age <- Age

braun_srt %>% write_rds('data05/braun_2022_fetal_brain_v2_common_hv2k_srt.rds')

braun_ctx_neuronal <- braun_srt %>% subset(
    CellClass%in%c('Neuron', 'Neuronal IPC', 'Radial glia') & Subregion == "Cortex"
)

layer_genes <- c('SOX2', 'EOMES', 'TBR1', 'BCL11B', 'SATB2', 'CUX1', 'FOXP2', 'NEUROD6', 'NEUROD2', 'FEZF2') %>% 
    intersect(rownames(braun_ctx_neuronal))
layer_expr <- braun_ctx_neuronal@assays$RNA@data[layer_genes,]
layer_expr <- layer_expr[,colSums(layer_expr)>0]

layer_expr_df <- layer_expr %>% 
    t() %>% as_tibble(rownames='cell') %>% 
    pivot_longer(!cell) %>% 
    inner_join(as_tibble(braun_ctx_neuronal@meta.data, rownames='cell'))

layer_frac_df <- layer_expr_df %>% 
    group_by(name, Age, Donor) %>% 
    dplyr::summarize(frac_pos=sum(value>0)/n()) %>% 
    filter(name%in%c('BCL11B', 'SOX2', 'SATB2', 'NEUROD2'))

ggplot(layer_frac_df, aes(factor(as.integer(Age)), frac_pos)) +
    geom_point() +
    stat_summary() +
    facet_grid(~name)



muo_w19 <- read_rds('data_/fetal/MUO/MUO_w19.rds')
muo_w19_use <- muo_w19 %>% subset(celltype_jf %in% c('cortical_excitatory_neuron', 'dorsal_npc'))

fet_layer_expr <- muo_w19_use@assays$RNA@data[layer_genes,]
fet_layer_expr <- fet_layer_expr[,colSums(fet_layer_expr)>0]

fet_layer_expr_df <- fet_layer_expr %>% 
    t() %>% as_tibble(rownames='cell') %>% 
    pivot_longer(!cell) %>% 
    inner_join(as_tibble(muo_w19_use@meta.data, rownames='cell'))

fet_layer_frac_df <- fet_layer_expr_df %>% 
    group_by(name) %>% 
    dplyr::summarize(frac_pos=sum(value>0)/n()) %>% 
    mutate(Age=19) %>% 
    filter(name%in%c('BCL11B', 'SOX2', 'SATB2', 'NEUROD2'))

plot_df <- bind_rows(layer_frac_df, fet_layer_frac_df)

ggplot(plot_df, aes(factor(as.integer(Age)), frac_pos)) +
    geom_point() +
    stat_summary() +
    facet_grid(~name)




rna <- read_rds('data_/RNA/RNA_all_srt_v2.3lines.rds')
org_layer_expr <- rna@assays$RNA@data[layer_genes,]
org_layer_expr <- org_layer_expr[,colSums(org_layer_expr)>0]

org_layer_expr_df <- org_layer_expr %>% 
    t() %>% as_tibble(rownames='cell') %>% 
    pivot_longer(!cell) %>% 
    inner_join(as_tibble(rna@meta.data, rownames='cell'))

org_layer_frac_df <- org_layer_expr_df %>% 
    filter(!age%in%c('ret_12w', 'ret_6w'), celltype_jf%in%c('ctx_npc', 'ctx_ex')) %>% 
    group_by(name, age, orig.ident) %>% 
    dplyr::summarize(frac_pos=sum(value>0)/n()) %>% 
    mutate(
        Age=as.numeric(as.character(fct_recode(age, 
            '60'='60d', '240'='8mo', '120'='4mo', '30'='35d'
    )))/30*4) %>% 
    filter(name%in%c('BCL11B', 'SOX2', 'SATB2', 'NEUROD2'))

plot_df <- bind_rows('Braun'=layer_frac_df, 'Organoid'=org_layer_frac_df, .id='origin')

ggplot(plot_df, aes(factor(as.integer(Age)), frac_pos, fill=origin)) +
    stat_summary(color='grey') +
    geom_point(shape=21, size=3) +
    facet_grid(~name+origin, scales='free', space='free') +
    theme(
        panel.grid.major.y = element_line(color='grey', linewidth=0.2)
    ) +
    labs(x='Age [weeks]', y='Fraction expressing cells')
ggsave('plots/paper/revplot_org_vs_braun.pdf', width=15, height=4)











