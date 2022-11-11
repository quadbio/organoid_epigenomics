source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

setwd('~/projects/cutntag/')


#### Read data ####
marks <- read_rds('data/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')
drugs <- read_rds('data/drugs/drugs_all_v1.1ids_srt.rds')


#### Demux drugs ####
drugs_demux_path <- '/local2/USERS/jfleck/projects/cutntag/demux/results/drugs/'
drugs_demux_files <- list.files(drugs_demux_path, pattern='*.best', full.names = T)
names(drugs_demux_files) <- list.files(drugs_demux_path, pattern='*.best') %>% str_remove('.best')

drugs_demux_all <- drugs_demux_files %>% 
    map_dfr(read_tsv, .id='sample')

drugs_demux <- drugs_demux_all %>% 
    filter(BARCODE!='.') %>% 
    mutate(BARCODE=str_replace(BARCODE, '([ATCG]+).+', '\\1')) %>% 
    mutate(cell=paste0(sample, '_', BARCODE))

drugs_demux %>% write_tsv('data/drugs/drugs_demux_all.tsv')

drugs_demux_meta <- drugs_demux %>% 
    filter(cell%in%colnames(drugs)) %>% 
    dplyr::select(cell, sample, 'line'=SNG.1ST, 'p_singlet'=PRB.SNG1, 'best'=BEST) %>% 
    mutate(status = str_replace(best, '(\\w+)-.+', '\\1')) %>% 
    dplyr::filter(!is.na(p_singlet)) 

drugs_demux_meta %>% ggplot(aes(status)) + geom_bar() + facet_wrap(~sample)
drugs_demux_meta %>% ggplot(aes(line)) + geom_bar() + facet_wrap(~sample)

intersect(drugs_demux_meta$sample, drugs$sample)

drugs <- AddMetaData(drugs, column_to_rownames(drugs_demux_meta, 'cell'))

drugs %>% dim_plot(group.by='line', reduction='cpumap')
drugs %>% write_rds('data/drugs/drugs_all_v1.2lines_srt.rds')


#### Demux marks ####
marks_demux_path <- '/local2/USERS/jfleck/projects/cutntag/demux/results/CT'
marks_demux_files <- list.files(marks_demux_path, pattern='*.best', full.names = T, recursive = T)
names(marks_demux_files) <- list.files(marks_demux_path)

marks_demux_all <- marks_demux_files %>% 
    map_dfr(read_tsv, .id='sample')

marks_demux <- marks_demux_all %>% 
    filter(BARCODE!='.') %>% 
    mutate(BARCODE=str_replace(BARCODE, '([ATCG]+).+', '\\1')) %>% 
    mutate(cell=paste0(sample, '_', BARCODE, '-1'))

marks_demux %>% write_tsv('data/CT_marks_demux_all.tsv')

marks_cells <- map(marks, colnames) %>% reduce(c)

marks_demux_meta <- marks_demux %>% 
    filter(cell%in%marks_cells) %>% 
    dplyr::select(cell, sample, 'line'=SNG.1ST, 'p_singlet'=PRB.SNG1, 'best'=BEST) %>% 
    mutate(status = str_replace(best, '(\\w+)-.+', '\\1')) %>% 
    dplyr::filter(!is.na(p_singlet)) %>% 
    mutate(mark=str_replace(sample, '.+(H3K27ac|H3K27me3|H3K4me3).+', '\\1'))

marks_demux_meta %>% ggplot(aes(status)) + geom_bar() + facet_wrap(mark~sample)
marks_demux_meta %>% ggplot(aes(line)) + geom_bar() + facet_wrap(mark~sample)


marks_new <- marks %>% map(function(x){
    x <- AddMetaData(x, column_to_rownames(marks_demux_meta, 'cell'))
    x$line[str_detect(x$orig.ident, '_ret')] <- 'B7'
    return(x)
})

marks_new$H3K27me3 %>% dim_plot(group.by='line')
marks_new %>% write_rds('data/all_marks_list_v3.4lines.rds')


plots <- map(marks_new, dim_plot, group.by='line')
wrap_plots(plots)


marks_meta <- map_dfr(marks_new, ~.x@meta.data) %>% 
    as_tibble(rownames='cell') %>% 
    filter(!is.na(mark), lineage!='retina')

ggplot(marks_meta, aes(lineage, fill=line)) +
    geom_bar(position='fill') +
    facet_grid(~mark)


#### Demux RNA ####
rna_demux_path <- '/local2/USERS/jfleck/projects/cutntag/demux/results/RNA/'
rna_demux_files <- list.files(rna_demux_path, pattern='*.best', full.names = T, recursive = T)
names(rna_demux_files) <- list.files(rna_demux_path)

rna_demux_all <- rna_demux_files %>% 
    map_dfr(read_tsv, .id='sample')

rna_demux <- rna_demux_all %>% 
    filter(BARCODE!='.') %>% 
    mutate(BARCODE=str_replace(BARCODE, '([ATCG]+).+', '\\1')) %>% 
    mutate(cell=paste0(sample, '_', BARCODE, '-1'))

rna_demux %>% write_tsv('data/RNA/RNA_demux_all.tsv')

rna_demux_meta <- rna_demux %>% 
    filter(cell%in%colnames(rna)) %>% 
    dplyr::select(cell, sample, 'line'=SNG.1ST, 'p_singlet'=PRB.SNG1, 'best'=BEST) %>% 
    mutate(status = str_replace(best, '(\\w+)-.+', '\\1')) %>% 
    dplyr::filter(!is.na(p_singlet)) 

rna_demux_meta %>% ggplot(aes(status)) + geom_bar() + facet_wrap(~sample)
rna_demux_meta %>% ggplot(aes(line)) + geom_bar() + facet_wrap(~sample)

rna <- AddMetaData(rna, column_to_rownames(rna_demux_meta, 'cell'))
rna$line[str_detect(rna$orig.ident, '_ret')] <- 'B7'
rna %>% dim_plot(group.by='line', reduction='cssumap')

rna %>% write_rds('data/RNA/RNA_all_srt_v2.3lines.rds')


marks_meta <- rna@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    filter(!is.na(mark), lineage!='retina')

ggplot(marks_meta, aes(lineage, fill=line)) +
    geom_bar(position='fill') 

