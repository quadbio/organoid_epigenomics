source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/utils/colors.R')
source('~/scripts/utils/utils.R')

library(GenomeInfoDb)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd('~/projects/cutntag/')

print('Reading input')
H3K27me3_srt <- read_rds('data/H3K27me3/H3K27me3_srt_v1filtered.rds')
H3K4me3_srt <- read_rds('data/H3K4me3/H3K4me3_srt_v1filtered.rds')
H3K27ac_srt <- read_rds('data/H3K27ac/H3K27ac_srt_v1filtered.rds')
gene_annot <- read_rds('~/resources/EnsDb.Hsapiens.v86_gene_annot_UCSC.hg38.rds')


#### H3K27me3 ####
print('H3K27me3')
H3K27me3_srt$orig.ident <- str_replace(colnames(H3K27me3_srt), '(.+)_[ATCG]+.*', '\\1')
H3K27me3_split <- SplitObject(H3K27me3_srt, split.by='orig.ident')

print('Adding fragments')
frag_path <- '/local2/USERS/jfleck/projects/cutntag/fragments/'
H3K27me3_split <- map(H3K27me3_split, function(srt){
    cell_names <- str_replace(colnames(srt), '.+_([ATCG]+.*)', '\\1')
    sample_name <- srt$orig.ident[[1]]
    names(cell_names) <- colnames(srt)
    path <- paste0(frag_path, sample_name, '/fragments.tsv.gz')
    frags <- CreateFragmentObject(
        path = path,
        cells = cell_names
    )
    Fragments(srt) <- NULL
    Fragments(srt) <- frags
    return(srt)
})

macs_path <- '/local2/USERS/jfleck/miniconda3/envs/signac/bin/macs2'

print('Running MACS2: narrow peaks')
H3K27me3_peaks <- map(H3K27me3_split, function(srt){
    print(srt$orig.ident[[1]])
    print(srt)
    # call peaks using MACS2
    peaks <- CallPeaks(srt, macs2.path = macs_path, outdir='/local1/TMP/jfleck/')
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
    peaks <- subsetByOverlaps(peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    # quantify counts in each peak
    return(peaks)
})

H3K27me3_peaks %>% write_rds('data/H3K27me3/H3K27me3_narrow_peaks_macs2.rds')
# H3K27me3_peaks <- read_rds('data/ATAC/ATAC_all_merged/ATAC_all_peaks_macs2.rds')

print('Running MACS2: broad peaks')
H3K27me3_broad_peaks <- map(H3K27me3_split, function(srt){
    print(srt$orig.ident[[1]])
    print(srt)
    # call peaks using MACS2
    peaks <- CallPeaks(srt, macs2.path = macs_path, outdir='/local1/TMP/jfleck/', broad=T)
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
    peaks <- subsetByOverlaps(peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    # quantify counts in each peak
    return(peaks)
})

H3K27me3_broad_peaks %>% write_rds('data/H3K27me3/H3K27me3_broad_peaks_macs2.rds')


print('Merging peaks')
# Create a unified set of peaks to quantify in each dataset
H3K27me3_peaks <- reduce(purrr::reduce(H3K27me3_peaks, c))
# Filter out bad peaks based on length
peakwidths <- width(H3K27me3_peaks)
H3K27me3_peaks <- H3K27me3_peaks[peakwidths  < 10000 & peakwidths > 10]

cells_use <- map(Fragments(H3K27me3_srt), function(x){
    names(x@cells)
}) %>% purrr::reduce(c)



print('Counting reads in peaks')
peak_counts <- FeatureMatrix(
    fragments = Fragments(H3K27me3_srt),
    features = H3K27me3_peaks,
    cells = cells_use,
    process_n = 5000
)
chassay <- CreateChromatinAssay(
    peak_counts,
    min.cells=-1,
    min.features=-1,
    # fragments = Fragments(H3K27me3_srt),
    annotation = gene_annot,
)
H3K27me3_srt[['narrow_peaks']] <- chassay

print('Writing output')
H3K27me3_srt %>% write_rds('data/H3K27me3/H3K27me3_srt_v2macs.rds')




#### H3K4me3 ####
print('H3K4me3')
H3K4me3_srt$orig.ident <- str_replace(colnames(H3K4me3_srt), '(.+)_[ATCG]+.*', '\\1')
H3K4me3_split <- SplitObject(H3K4me3_srt, split.by='orig.ident')

print('Adding fragments')
frag_path <- '/local2/USERS/jfleck/projects/cutntag/fragments/'
H3K4me3_split <- map(H3K4me3_split, function(srt){
    cell_names <- str_replace(colnames(srt), '.+_([ATCG]+.*)', '\\1')
    sample_name <- srt$orig.ident[[1]]
    names(cell_names) <- colnames(srt)
    path <- paste0(frag_path, sample_name, '/fragments.tsv.gz')
    frags <- CreateFragmentObject(
        path = path,
        cells = cell_names
    )
    Fragments(srt) <- NULL
    Fragments(srt) <- frags
    return(srt)
})

macs_path <- '/local2/USERS/jfleck/miniconda3/envs/signac/bin/macs2'

print('Running MACS2: narrow peaks')
H3K4me3_peaks <- map(H3K4me3_split, function(srt){
    print(srt$orig.ident[[1]])
    print(srt)
    # call peaks using MACS2
    peaks <- CallPeaks(srt, macs2.path = macs_path, outdir='/local1/TMP/jfleck/')
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
    peaks <- subsetByOverlaps(peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    # quantify counts in each peak
    return(peaks)
})

H3K4me3_peaks %>% write_rds('data/H3K4me3/H3K4me3_narrow_peaks_macs2.rds')
# H3K4me3_peaks <- read_rds('data/ATAC/ATAC_all_merged/ATAC_all_peaks_macs2.rds')

print('Running MACS2: broad peaks')
H3K4me3_broad_peaks <- map(H3K4me3_split, function(srt){
    print(srt$orig.ident[[1]])
    print(srt)
    # call peaks using MACS2
    peaks <- CallPeaks(srt, macs2.path = macs_path, outdir='/local1/TMP/jfleck/', broad=T)
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
    peaks <- subsetByOverlaps(peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    # quantify counts in each peak
    return(peaks)
})

H3K4me3_broad_peaks %>% write_rds('data/H3K4me3/H3K4me3_broad_peaks_macs2.rds')



print('Merging peaks')
# Create a unified set of peaks to quantify in each dataset
H3K4me3_peaks <- reduce(purrr::reduce(H3K4me3_peaks, c))
# Filter out bad peaks based on length
peakwidths <- width(H3K4me3_peaks)
H3K4me3_peaks <- H3K4me3_peaks[peakwidths  < 10000 & peakwidths > 10]

cells_use <- map(Fragments(H3K4me3_srt), function(x){
    names(x@cells)
}) %>% purrr::reduce(c)



print('Counting reads in peaks')
peak_counts <- FeatureMatrix(
    fragments = Fragments(H3K4me3_srt),
    features = H3K4me3_peaks,
    cells = cells_use,
    process_n = 5000
)
chassay <- CreateChromatinAssay(
    peak_counts,
    min.cells=-1,
    min.features=-1,
    fragments = Fragments(H3K4me3_srt),
    annotation = gene_annot,
)
H3K4me3_srt[['narrow_peaks']] <- chassay


print('Writing output')
H3K4me3_srt %>% write_rds('data/H3K4me3/H3K4me3_srt_v2macs.rds')



#### H3K27ac ####
print('H3K27ac')
H3K27ac_srt$orig.ident <- str_replace(colnames(H3K27ac_srt), '(.+)_[ATCG]+.*', '\\1')
H3K27ac_split <- SplitObject(H3K27ac_srt, split.by='orig.ident')

print('Adding fragments')
frag_path <- '/local2/USERS/jfleck/projects/cutntag/fragments/'
H3K27ac_split <- map(H3K27ac_split, function(srt){
    cell_names <- str_replace(colnames(srt), '.+_([ATCG]+.*)', '\\1')
    sample_name <- srt$orig.ident[[1]]
    names(cell_names) <- colnames(srt)
    path <- paste0(frag_path, sample_name, '/fragments.tsv.gz')
    frags <- CreateFragmentObject(
        path = path,
        cells = cell_names
    )
    Fragments(srt) <- NULL
    Fragments(srt) <- frags
    return(srt)
})

macs_path <- '/local2/USERS/jfleck/miniconda3/envs/signac/bin/macs2'

print('Running MACS2: narrow peaks')
H3K27ac_peaks <- map(H3K27ac_split, function(srt){
    print(srt$orig.ident[[1]])
    print(srt)
    # call peaks using MACS2
    peaks <- CallPeaks(srt, macs2.path = macs_path, outdir='/local1/TMP/jfleck/')
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
    peaks <- subsetByOverlaps(peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    # quantify counts in each peak
    return(peaks)
})

H3K27ac_peaks %>% write_rds('data/H3K27ac/H3K27ac_narrow_peaks_macs2.rds')
# H3K27ac_peaks <- read_rds('data/ATAC/ATAC_all_merged/ATAC_all_peaks_macs2.rds')


print('Running MACS2: broad peaks')
H3K27ac_broad_peaks <- map(H3K27ac_split, function(srt){
    print(srt$orig.ident[[1]])
    print(srt)
    # call peaks using MACS2
    peaks <- CallPeaks(srt, macs2.path = macs_path, outdir='/local1/TMP/jfleck/', broad=T)
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
    peaks <- subsetByOverlaps(peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    # quantify counts in each peak
    return(peaks)
})

H3K27ac_broad_peaks %>% write_rds('data/H3K27ac/H3K27ac_broad_peaks_macs2.rds')



print('Merging peaks')
# Create a unified set of peaks to quantify in each dataset
H3K27ac_peaks <- reduce(purrr::reduce(H3K27ac_peaks, c))
# Filter out bad peaks based on length
peakwidths <- width(H3K27ac_peaks)
H3K27ac_peaks <- H3K27ac_peaks[peakwidths  < 10000 & peakwidths > 10]

cells_use <- map(Fragments(H3K27ac_srt), function(x){
    names(x@cells)
}) %>% purrr::reduce(c)



print('Counting reads in peaks')
peak_counts <- FeatureMatrix(
    fragments = Fragments(H3K27ac_srt),
    features = H3K27ac_peaks,
    cells = cells_use,
    process_n = 5000
)
chassay <- CreateChromatinAssay(
    peak_counts,
    min.cells=-1,
    min.features=-1,
    fragments = Fragments(H3K27ac_srt),
    annotation = gene_annot,
)
H3K27ac_srt[['narrow_peaks']] <- chassay


print('Writing output')
H3K27ac_srt %>% write_rds('data/H3K27ac/H3K27ac_srt_v2macs.rds')



