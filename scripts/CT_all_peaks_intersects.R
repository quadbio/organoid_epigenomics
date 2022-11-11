source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')


setwd('~/projects/cutntag/')

marks <- read_rds('data/CT/all_marks_list_v3.3motifs.rds')

detect_peaks <- map(marks, function(m){
    m[['peaks_bin']] <- CreateAssayObject((m[['peaks']]@data > 0)*1)
    m <- Pando::aggregate_assay(m, assay='peaks_bin', group_name='clusters')
    return(m@assays$peaks_bin@misc$summary$clusters %>% {.[, colMaxs(.)>0.05]} %>% colnames())
})


H3K27ac_peak_ranges <- StringToGRanges(rownames(marks$H3K27ac[['peaks']][detect_peaks$H3K27ac]))
H3K27me3_peak_ranges <- StringToGRanges(rownames(marks$H3K27me3[['peaks']][detect_peaks$H3K27me3]))
H3K4me3_peak_ranges <- StringToGRanges(rownames(marks$H3K4me3[['peaks']][detect_peaks$H3K4me3]))

AB_intersect <- IRanges::findOverlaps(H3K27ac_peak_ranges, H3K27me3_peak_ranges)
H3K27ac_AB <- H3K27ac_peak_ranges[unique(queryHits(AB_intersect))]
H3K27me3_AB <- H3K27me3_peak_ranges[unique(subjectHits(AB_intersect))]

AC_intersect <- IRanges::findOverlaps(H3K27ac_peak_ranges, H3K4me3_peak_ranges)
H3K27ac_AC <- H3K27ac_peak_ranges[unique(queryHits(AC_intersect))]
H3K4me3_AC <- H3K4me3_peak_ranges[unique(subjectHits(AC_intersect))]

BC_intersect <- IRanges::findOverlaps(H3K27me3_peak_ranges, H3K4me3_peak_ranges)
H3K27me3_BC <- H3K27me3_peak_ranges[unique(queryHits(BC_intersect))]
H3K4me3_BC <- H3K4me3_peak_ranges[unique(subjectHits(BC_intersect))]

ABC_intersect <- IRanges::intersect(H3K27ac_peak_ranges, H3K27me3_peak_ranges) %>% IRanges::intersect(H3K4me3_peak_ranges)
H3K27ac_ABC <- H3K27ac_peak_ranges[unique(queryHits(findOverlaps(H3K27ac_peak_ranges, ABC_intersect)))]
H3K27me3_ABC <- H3K27me3_peak_ranges[unique(queryHits(findOverlaps(H3K27me3_peak_ranges, ABC_intersect)))]
H3K4me3_ABC <- H3K4me3_peak_ranges[unique(queryHits(findOverlaps(H3K4me3_peak_ranges, ABC_intersect)))]


H3K27ac_ex <- H3K27ac_peak_ranges[-union(queryHits(AB_intersect), queryHits(AC_intersect))]
H3K27me3_ex <- H3K27me3_peak_ranges[-union(subjectHits(AB_intersect), queryHits(BC_intersect))]
H3K4me3_ex <- H3K4me3_peak_ranges[-union(subjectHits(AC_intersect), subjectHits(BC_intersect))]

intersect_counts <- tibble(
    A = length(H3K27ac_ex),
    B = length(H3K27me3_ex),
    C = length(H3K4me3_ex),
    AB = length(H3K27ac_AB),
    # BA = length(H3K27me3_AB),
    AC = length(H3K27ac_AC),
    # CA = length(H3K4me3_AC),
    BC = length(H3K27me3_BC),
    # CB = length(H3K4me3_BC),
    ABC = length(H3K27ac_ABC),
    # BCA = length(H3K27me3_ABC),
    # CAB = length(H3K4me3_ABC)
) %>% pivot_longer(everything()) %>% 
    mutate(name=factor(name, levels=unique(.$name)))

ggplot(intersect_counts, aes(name, value/1000)) +
    geom_bar(stat='identity', size=0.2, color='black', fill='grey') +
    scale_x_discrete(labels=c('K27ac', 'K27me3', 'K4me3', 'K27ac + K27me3', 'K27ac + K4me3', 'K27me3 + K4me3', 'K27ac + K27me3 + K4me3')) +
    article_text() +
    rotate_x_text(40) +
    labs(x='Intersection', y='Count in thousends')
ggsave('plots/paper/fig2/fig2_peak_intersects_counts_detect05_bar.pdf', width=5, height=5, units='cm')


intersect_counts <- tibble(
    H3K27ac = length(H3K27ac_peak_ranges),
    H3K27me3 = length(H3K27me3_peak_ranges),
    H3K4me3 = length(H3K4me3_peak_ranges)
) %>% pivot_longer(everything()) %>% 
    mutate(name=factor(name, levels=unique(.$name)))

ggplot(intersect_counts, aes(name, value/1000)) +
    geom_bar(stat='identity', size=0.2, color='black', fill='grey') +
    scale_x_discrete(labels=c('K27ac', 'K27me3', 'K4me3', 'K27ac + K27me3', 'K27ac + K4me3', 'K27me3 + K4me3', 'K27ac + K27me3 + K4me3')) +
    article_text() +
    rotate_x_text(40) +
    labs(x='Modality', y='Count in thousends')

ggsave('plots/paper/fig2/fig2_peak_counts_detect05_bar.pdf', width=3, height=5, units='cm')



intersect_counts <- tibble(
    A = length(H3K27ac_ex),
    B = length(H3K27me3_ex),
    C = length(H3K4me3_ex),
    AB = length(H3K27ac_AB),
    BA = length(H3K27me3_AB),
    AC = length(H3K27ac_AC),
    CA = length(H3K4me3_AC),
    BC = length(H3K27me3_BC),
    CB = length(H3K4me3_BC),
    ABC = length(H3K27ac_ABC),
    BCA = length(H3K27me3_ABC),
    CAB = length(H3K4me3_ABC)
) %>% pivot_longer(everything()) %>% 
    mutate(name=factor(name, levels=unique(.$name)))

ggplot(intersect_counts, aes(name, value/1000)) +
    geom_bar(stat='identity') +
    rotate_x_text(40) +
    labs(x='Intersection', y='Count in thousends')


peaks_sets <- list(
    H3K27ac_ex = -union(queryHits(AB_intersect), queryHits(AC_intersect)),
    H3K27me3_ex = -union(subjectHits(AB_intersect), queryHits(BC_intersect)),
    H3K4me3_ex = -union(subjectHits(AC_intersect), subjectHits(BC_intersect)),
    H3K27ac_AB = queryHits(AB_intersect),
    H3K27ac_AC = queryHits(AC_intersect),
    H3K27me3_AB = subjectHits(AB_intersect),
    H3K27me3_BC = queryHits(BC_intersect),
    H3K4me3_AC = subjectHits(AC_intersect),
    H3K4me3_BC = queryHits(BC_intersect),
    H3K27ac_AC = queryHits(AC_intersect),
    H3K27ac_AC = queryHits(AC_intersect),
    H3K27ac_ABC = queryHits(findOverlaps(H3K27ac_peak_ranges, ABC_intersect)),
    H3K27me3_ABC = queryHits(findOverlaps(H3K27me3_peak_ranges, ABC_intersect)),
    H3K4me3_ABC = queryHits(findOverlaps(H3K4me3_peak_ranges, ABC_intersect))
)

peaks_sets %>% write_rds('data/intersect/all_peak_intersects.rds')



