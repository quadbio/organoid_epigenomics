source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

library(ArchR)
library(parallel)
addArchRThreads(threads = 16) 
addArchRGenome('hg38')

setwd('~/projects/cutntag/')


#### Read stuff ####
marks <- read_rds('data/all_marks_list_v3.2annot.rds')

H3K27me3_proj <- loadArchRProject('data/arrow/H3K27me3/')
H3K27ac_proj <- loadArchRProject('data/arrow/H3K27ac/')
H3K4me3_proj <- loadArchRProject('data/arrow/H3K4me3/')


#### Transfer and export bigwigs ####
arch_meta <- marks$H3K27me3@meta.data[H3K27me3_proj$srtNames, ]
H3K27me3_proj$celltype_jf <- arch_meta$celltype_jf

arch_meta <- marks$H3K27ac@meta.data[H3K27ac_proj$srtNames, ]
H3K27ac_proj$celltype_jf <- arch_meta$celltype_jf

arch_meta <- marks$H3K4me3@meta.data[H3K4me3_proj$srtNames, ]
H3K4me3_proj$celltype_jf <- arch_meta$celltype_jf

getGroupBW(H3K27me3_proj, groupBy = 'celltype_jf')
getGroupBW(H3K27ac_proj, groupBy = 'celltype_jf')
getGroupBW(H3K4me3_proj, groupBy = 'celltype_jf')



arch_meta <- marks$H3K27me3@meta.data[H3K27me3_proj$srtNames, ]
H3K27me3_proj$lineage <- arch_meta$lineage

arch_meta <- marks$H3K27ac@meta.data[H3K27ac_proj$srtNames, ]
H3K27ac_proj$lineage <- arch_meta$lineage

arch_meta <- marks$H3K4me3@meta.data[H3K4me3_proj$srtNames, ]
H3K4me3_proj$lineage <- arch_meta$lineage

getGroupBW(H3K27me3_proj, groupBy = 'lineage')
getGroupBW(H3K27ac_proj, groupBy = 'lineage')
getGroupBW(H3K4me3_proj, groupBy = 'lineage')


arch_meta <- marks$H3K27me3@meta.data[H3K27me3_proj$srtNames, ]
H3K27me3_proj$state <- arch_meta$state

arch_meta <- marks$H3K27ac@meta.data[H3K27ac_proj$srtNames, ]
H3K27ac_proj$state <- arch_meta$state

arch_meta <- marks$H3K4me3@meta.data[H3K4me3_proj$srtNames, ]
H3K4me3_proj$state <- arch_meta$state

getGroupBW(H3K27me3_proj, groupBy = 'state')
getGroupBW(H3K27ac_proj, groupBy = 'state')
getGroupBW(H3K4me3_proj, groupBy = 'state')


arch_meta <- marks$H3K27me3@meta.data[H3K27me3_proj$srtNames, ]
H3K27me3_proj$age <- arch_meta$age

arch_meta <- marks$H3K27ac@meta.data[H3K27ac_proj$srtNames, ]
H3K27ac_proj$age <- arch_meta$age

arch_meta <- marks$H3K4me3@meta.data[H3K4me3_proj$srtNames, ]
H3K4me3_proj$age <- arch_meta$age

getGroupBW(H3K27me3_proj, groupBy = 'age')
getGroupBW(H3K27ac_proj, groupBy = 'age')
getGroupBW(H3K4me3_proj, groupBy = 'age')


arch_meta <- marks$H3K27me3@meta.data[H3K27me3_proj$srtNames, ]
H3K27me3_proj$stage <- arch_meta$stage

arch_meta <- marks$H3K27ac@meta.data[H3K27ac_proj$srtNames, ]
H3K27ac_proj$stage <- arch_meta$stage

arch_meta <- marks$H3K4me3@meta.data[H3K4me3_proj$srtNames, ]
H3K4me3_proj$stage <- arch_meta$stage

getGroupBW(H3K27me3_proj, groupBy = 'stage')
getGroupBW(H3K27ac_proj, groupBy = 'stage')
getGroupBW(H3K4me3_proj, groupBy = 'stage')






