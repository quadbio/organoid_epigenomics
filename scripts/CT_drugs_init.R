source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

setwd('~/projects/cutntag/')


#### Read stuff ####

drugs_init <- read_rds('/links/groups/treutlein/USERS/fides/CutNTag_analysis/hashing/hashing/202205_drugs_all_combined/20220414_annotation_Jonas_all_data_no_unmap.rds')
# drugs_init %>% write_rds('data/drugs/drugs_all_v1fides_srt.rds')

old_names <- drugs_init %>% colnames
barcodes <- old_names %>% str_remove('_\\d') 

bc_ids1 <- old_names %>% str_extract('_\\d') 
bc_ids2 <- barcodes %>% str_remove('[ATCG]+') 

table(drugs_init$orig.ident, bc_ids1)
table(drugs_init$orig.ident, bc_ids2)

sample_names <- case_when(
    bc_ids2=='_4' & bc_ids1=='_1' ~ '210412_1_FZ_scRNA_WH_drug_15d',
    bc_ids2=='_4' & bc_ids1=='_2' ~ '210412_2_FZ_scRNA_WH_drug_15d',
    bc_ids1=='_1' ~ '21033005_1_FZ_scRNA_DMSO_HWW4B_21d',
    bc_ids1=='_2' ~ '21033014_3_FZ_scRNA_A485_HWW4B_21d',
    bc_ids1=='_3' ~ '21033008_2_FZ_scRNA_A395_HWW4B_21d',
    bc_ids1=='_5' ~ '220321_1_FZ_scRNA_drugs_d18_WH',
    bc_ids1=='_6' ~ '220321_2_FZ_scRNA_drugs_d21_WHWB4'
)

barcodes_new <- barcodes %>% str_replace('([ATCG]+)(-1|_4)?', '\\1') 
cellnames_new <- paste0(sample_names, '_', barcodes_new)

drugs_init$sample <- sample_names
drugs_init$barcode <- barcodes_new
drugs_init <- RenameCells(drugs_init, new.names=cellnames_new)

drugs_init %>% write_rds('data/drugs/drugs_all_v1.1ids_srt.rds')



