source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

library(ChIPseeker)
library(ArchR)
library(parallel)
addArchRThreads(threads = 16) 
addArchRGenome('hg38')



setwd('~/projects/cutntag/')
marks <- read_rds('data/all_marks_list_v3.2annot.rds')

setwd('~/projects/cutntag/data/arrow')
map(names(marks), function(n){
    print(n)
    mark <- marks[[n]]
    mark_frags <- Fragments(mark)
    names(mark_frags) <- map_chr(mark_frags, ~ str_split(.x@path, '/')[[1]][8])
    mark_paths <- map_chr(mark_frags, ~.x@path)
    
    mark_arrow <- createArrowFiles(
        inputFiles = mark_paths,
        minFrags = 0,
        maxFrags = 1e8,
        minTSS = 0
    )
    
    mark_proj_ <- ArchRProject(
        ArrowFiles = mark_arrow, 
        outputDirectory = n,
        copyArrows = TRUE 
    )
    
    mark_proj_$srtNames <- mark_proj_$cellNames %>% str_replace('#', '_')
    mark_proj <- mark_proj_[mark_proj_$srtNames%in%colnames(mark), ]
    arch_meta <- mark@meta.data[mark_proj$srtNames, ]
    
    mark_proj$lineage <- arch_meta$lineage
    mark_proj$celltype_jf <- arch_meta$celltype_jf
    mark_proj$age <- arch_meta$age
    mark_proj$state <- arch_meta$state
    mark_proj$stage <- arch_meta$stage
    mark_proj$clusters <- arch_meta$clusters
    
    saveArchRProject(mark_proj)
})



