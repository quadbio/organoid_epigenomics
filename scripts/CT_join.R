source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

setwd('~/projects/cutntag/')

counts_path <- '/local2/USERS/jfleck/projects/cutntag/counts/CT'
counts_list <- list.files(counts_path, full.names = T)
counts_names <- list.files(counts_path)
names(counts_list) <- counts_names


frags_path <- '/local2/USERS/jfleck/projects/cutntag/fragments'
frags_list <- list.files(frags_path, full.names = T)
frags_names <- list.files(frags_path)
frags_files <- paste0(frags_list, '/fragments.tsv.gz')
names(frags_files) <- frags_names


matrix_list <- map(counts_list, function(x){
    path=paste0(x, '/raw_peak_bc_matrix.h5')
    print(path)
    return(Read10X_h5(path))
})

srt_list <- map(set_names(names(matrix_list)), function(x){
    counts <- matrix_list[[x]]
    rownames(counts) <- str_replace(rownames(counts), ':', '-')
    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        genome = 'hg38',
        fragments = frags_files[x]
    )
    srt <- CreateSeuratObject(
        counts=chrom_assay,
        assay='peaks'
    )
    return(srt)
})

srt_list_annot <- map(set_names(names(counts_list)), function(x){
    path=paste0(counts_list[x], '/singlecell.csv')
    meta <- read_csv(path) %>% column_to_rownames('barcode')
    srt <- AddMetaData(srt_list[[x]], meta)
    srt$pct_reads_in_peaks <- srt$peak_region_fragments / srt$passed_filters * 100
    srt$blacklist_ratio <- srt$blacklist_region_fragments / srt$peak_region_fragments
    return(srt)
})

chr_mark <- str_replace(names(srt_list_annot), '.+(H3[0-9a-zA-Z]+)_.+', '\\1')

H3K27me3_srt_list <- srt_list_annot[chr_mark == 'H3K27me3']
H3K4me3_srt_list <- srt_list_annot[chr_mark == 'H3K4me3']
H3K27ac_srt_list <- srt_list_annot[chr_mark == 'H3K27ac']
H3K9ac_srt_list <- srt_list_annot[chr_mark == 'H3K9ac']

H3K27me3_srt_list %>% write_rds('data/H3K27me3/H3K27me3_srt_list_unfiltered.rds')
H3K4me3_srt_list %>% write_rds('data/H3K4me3/H3K4me3_srt_list_unfiltered.rds')
H3K27ac_srt_list %>% write_rds('data/H3K27ac/H3K27ac_srt_list_unfiltered.rds')
H3K9ac_srt_list %>% write_rds('data/H3K9ac/H3K9ac_srt_list_unfiltered.rds')


#### Preproc everything ####
gene_ranges <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')

do_preproc <- function(sampl, sampl_list){
    srt <- sampl_list[[sampl]]
    print(srt)
    srt <- subset(srt, nCount_peaks>0)
    plot_dir <- paste0('plots/cutntag/individual_samples/', sampl, '/')
    data_dir <- paste0('data/individual_samples/', sampl, '/')
    dir.create(plot_dir)
    dir.create(data_dir)

    Annotation(srt) <- gene_ranges
    srt <- NucleosomeSignal(object = srt)

    meta <- as_tibble(srt@meta.data, rownames='cell')

    write_tsv(meta, paste0(data_dir, sampl, '_meta.tsv'))

    p1 <- ggplot(meta, aes(pct_reads_in_peaks)) +
        geom_histogram() +
        geom_vline(xintercept = 30)
    p2 <- ggplot(meta, aes(peak_region_fragments)) +
        geom_histogram() +
        geom_vline(xintercept = 5000) +
        geom_vline(xintercept = 80000)
    p3 <- ggplot(meta, aes(blacklist_ratio)) +
        geom_histogram() +
        geom_vline(xintercept = 1e-3)
    p4 <- ggplot(meta, aes(nucleosome_signal)) +
        geom_histogram() +
        geom_vline(xintercept = 10)
    p5 <- ggplot(meta, aes(nCount_peaks)) +
        scale_x_log10() +
        geom_histogram()
    p6 <- ggplot(meta, aes(nFeature_peaks)) +
        scale_x_log10() +
        geom_histogram()

    pout <- (p1 + p2)/(p3 + p4)/(p5 + p6) & theme_article()
    ggsave(paste0(plot_dir, 'basic_metrics.pdf'), plot=pout, width=10, height=15)


    srt$nucleosome_group <- ifelse(srt$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

    pout <- FragmentHistogram(object = srt, group.by = 'nucleosome_group') & theme_article()
    ggsave(paste0(plot_dir, 'fragment_hist.pdf'), plot=pout)

    write_rds(srt, paste0(data_dir, sampl, '.rds'))

    return(srt)
}

H3K27me3_list <- map(names(H3K27me3_srt_list), do_preproc, sampl_list=H3K27me3_srt_list)
names(H3K27me3_list) <- names(H3K27me3_srt_list)
H3K4me3_list <- map(names(H3K4me3_srt_list), do_preproc, sampl_list=H3K4me3_srt_list)
names(H3K4me3_list) <- names(H3K4me3_srt_list)
H3K27ac_list <- map(names(H3K27ac_srt_list), do_preproc, sampl_list=H3K27ac_srt_list)
names(H3K27ac_list) <- names(H3K27ac_srt_list)


H3K27me3_list %>% write_rds('data/H3K27me3/H3K27me3_srt_list_preproc.rds')
H3K4me3_list %>% write_rds('data/H3K4me3/H3K4me3_srt_list_preproc.rds')
H3K27ac_list %>% write_rds('data/H3K27ac/H3K27ac_srt_list_preproc.rds')



#### Check filters again ####
H3K27me3_meta <- map_dfr(H3K27me3_list, function(x){
    as_tibble(x@meta.data, rownames='cell')
}, .id='sample')

H3K4me3_meta <- map_dfr(H3K4me3_list, function(x){
    as_tibble(x@meta.data, rownames='cell')
}, .id='sample')

H3K27ac_meta <- map_dfr(H3K27ac_list, function(x){
    as_tibble(x@meta.data, rownames='cell')
}, .id='sample')


 all_meta <- bind_rows('H3K27me3'=H3K27me3_meta, 'H3K4me3'=H3K4me3_meta, 'H3K27ac'=H3K27ac_meta, .id='mark') %>%
    mutate(age=str_replace(sample, '.+_(.+)', '\\1')) %>%
    mutate(age=ifelse(str_detect(sample, 'rep2_35d'), '65d', age)) %>%
    filter(sample!='210616_3_FZ_scCT_H3K27ac_wetting_HWB_128d')

p1 <- ggplot(all_meta, aes(nFeature_peaks, fill=mark)) +
    geom_histogram() +
    facet_grid(age~mark) +
    scale_x_log10(breaks=c(1,100,1000,10000)) +
    theme_bw()

p2 <- ggplot(all_meta, aes(nCount_peaks, fill=mark)) +
    geom_histogram() +
    facet_grid(age~mark) +
    scale_x_log10(breaks=c(1,100,1000,10000)) +
    theme_bw()

p1 + p2 & geom_vline(xintercept = 20, color='black') & geom_vline(xintercept = c(50, 100, 200), color='darkgrey')
ggsave('plots/cutntag/basic_qc.png', width=15, height=8)



#### Apply filters ####

# H3K27ac_list <- H3K27ac_list[-1]
H3K27ac_list_fltr <- map(H3K27ac_list, function(x){
    subset(x, nCount_peaks>50)
})

H3K4me3_list_fltr <- map(H3K4me3_list, function(x){
    subset(x, nCount_peaks>50)
})

H3K27me3_list_fltr <- map(H3K27me3_list, function(x){
    subset(x, nCount_peaks>50)
})

H3K27me3_list_fltr$`210617_7_FZ_scCT_H3K27me3_HWW4CB_rep2_35d` <- H3K27me3_list_fltr$`210617_7_FZ_scCT_H3K27me3_HWW4CB_rep2_35d` %>%
    subset(nCount_peaks>200)



#### Replot ####
H3K27me3_meta <- map_dfr(H3K27me3_list_fltr, function(x){
    as_tibble(x@meta.data, rownames='cell')
}, .id='sample')

H3K4me3_meta <- map_dfr(H3K4me3_list_fltr, function(x){
    as_tibble(x@meta.data, rownames='cell')
}, .id='sample')

H3K27ac_meta <- map_dfr(H3K27ac_list_fltr, function(x){
    as_tibble(x@meta.data, rownames='cell')
}, .id='sample')

all_meta <- bind_rows('H3K27me3'=H3K27me3_meta, 'H3K4me3'=H3K4me3_meta, 'H3K27ac'=H3K27ac_meta, .id='mark') %>%
    mutate(age=str_replace(sample, '.+_(.+)', '\\1')) %>%
    mutate(age=ifelse(str_detect(sample, 'rep2_35d'), '65d', age))

p1 <- ggplot(all_meta, aes(nFeature_peaks, fill=mark)) +
    geom_histogram() +
    facet_grid(age~mark) +
    scale_x_log10(breaks=c(1,100,1000,10000)) +
    theme_bw()

p2 <- ggplot(all_meta, aes(nCount_peaks, fill=mark)) +
    geom_histogram() +
    facet_grid(age~mark) +
    scale_x_log10(breaks=c(1,100,1000,10000)) +
    theme_bw()

p1 + p2 & geom_vline(xintercept = 20, color='black') & geom_vline(xintercept = c(50, 100, 200), color='darkgrey')
ggsave('plots/basic_qc_filtered.png', width=15, height=6)


H3K27me3_list_fltr %>% write_rds('data/H3K27me3/H3K27me3_srt_list_filtered.rds')
H3K4me3_list_fltr %>% write_rds('data/H3K4me3/H3K4me3_srt_list_filtered.rds')
H3K27ac_list_fltr %>% write_rds('data/H3K27ac/H3K27ac_srt_list_filtered.rds')



# #### Merge #### 
H3K27me3_list_fltr <- read_rds('data/H3K27me3/H3K27me3_srt_list_filtered.rds')
H3K4me3_list_fltr <- read_rds('data/H3K4me3/H3K4me3_srt_list_filtered.rds')
H3K27ac_list_fltr <- read_rds('data/H3K27ac/H3K27ac_srt_list_filtered.rds')


H3K27ac_fltr <- merge(H3K27ac_list_fltr[[1]], H3K27ac_list_fltr[2:length(H3K27ac_list_fltr)], add.cell.ids=names(H3K27ac_list_fltr))
H3K27ac_fltr %>% write_rds('data/H3K27ac/H3K27ac_srt_v1filtered.rds')

H3K4me3_fltr <- merge(H3K4me3_list_fltr[[1]], H3K4me3_list_fltr[2:length(H3K4me3_list_fltr)], add.cell.ids=names(H3K4me3_list_fltr))
H3K4me3_fltr %>% write_rds('data/H3K4me3/H3K4me3_srt_v1filtered.rds')

H3K27me3_fltr <- merge(H3K27me3_list_fltr[[1]], H3K27me3_list_fltr[2:length(H3K27me3_list_fltr)], add.cell.ids=names(H3K27me3_list_fltr))
H3K27me3_fltr %>% write_rds('data/H3K27me3/H3K27me3_srt_v1filtered.rds')





